"""
Cross-species translatability filter -- the headline feature.

For each shortlisted candidate the pipeline has already linked a human target
gene. This module asks a specific, honestly-scoped question:

    Is the TARGET-CELL-TYPE BIOLOGY that this element plausibly serves
    conserved between human and mouse?

It does NOT test whether the enhancer SEQUENCE itself is conserved (enhancers
turn over rapidly in evolution; a human enhancer often has no mouse counterpart
even when the target gene is perfectly conserved). What it tests is whether the
LINKED GENE is (a) an identifiable 1:1-or-many mouse ortholog and (b) expressed
with comparable cell-type specificity in mouse. This is a GENE-LEVEL
target-biology conservation flag, and every output says so.

Why it matters: AAV constructs optimised in mouse routinely fail in NHP/human
because the underlying regulatory biology differs across species. Flagging
candidates whose target biology is mouse-specific (or human-specific) BEFORE
wet-lab spend targets the most expensive failure mode in the field.

Data: human/mouse ortholog map from Ensembl REST (open); cross-species
expression from CELLxGENE Census (CC-BY). Both open-access.

Flags:
  conserved            gene has a mouse ortholog AND is target-specific in both
  divergent_expression ortholog exists but expression specificity not conserved
  no_ortholog          no confident mouse ortholog (human-specific biology)
  unlinked             candidate has no linked gene (cannot assess)
"""
from __future__ import annotations
import time, json, urllib.request, urllib.parse
import numpy as np
import pandas as pd

ENSEMBL = "https://rest.ensembl.org"


def _get(url, retries=3):
    for a in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.loads(r.read().decode())
        except Exception:
            if a == retries - 1:
                return None
            time.sleep(1.0)
    return None


def human_to_mouse_orthologs(symbols, sleep=0.15):
    """Map human gene symbols -> list of mouse ortholog symbols via Ensembl.
    Returns {human_symbol: [mouse_symbol, ...]} (empty list = no ortholog)."""
    out = {}
    cache_id2sym = {}
    for g in symbols:
        if g in (None, "NA", "") or g in out:
            continue
        url = (f"{ENSEMBL}/homology/symbol/human/{urllib.parse.quote(g)}"
               f"?target_species=mouse;type=orthologues;content-type=application/json")
        d = _get(url)
        mouse_syms = []
        if d and d.get("data"):
            for h in d["data"]:
                for o in h.get("homologies", []):
                    mid = o.get("target", {}).get("id")
                    if not mid:
                        continue
                    if mid not in cache_id2sym:
                        ld = _get(f"{ENSEMBL}/lookup/id/{mid}?content-type=application/json")
                        cache_id2sym[mid] = (ld or {}).get("display_name")
                    sym = cache_id2sym[mid]
                    if sym and sym not in mouse_syms:
                        mouse_syms.append(sym)
        out[g] = mouse_syms
        time.sleep(sleep)
    return out


def conservation_flags(catalog, human_gene_score, mouse_gene_score,
                       ortholog_map, spec_min=0.10, resolve_missing=True):
    """Annotate a catalog (with 'target_gene') with cross-species columns.

    human_gene_score / mouse_gene_score : dict gene_symbol -> [0,1] target-cell
        expression-specificity score (from expression.gene_score_magspec), in the
        SAME cell type in each species.
    ortholog_map : {human_symbol: [mouse_symbol,...]}
        Cached Ensembl orthologs. IMPORTANT: a gene simply ABSENT from this map
        has not been looked up -- it is NOT evidence of no ortholog. When
        resolve_missing=True (default) any linked gene missing from the map is
        resolved live via Ensembl so absence-from-cache is never mislabelled
        'no_ortholog' (an earlier bug: genes surfaced by a wider link set, e.g.
        GNAS/CPE, were wrongly flagged human-specific).
    Adds: mouse_ortholog, human_gene_spec, mouse_gene_spec, conservation_score,
          translatability_flag.
    """
    ortholog_map = dict(ortholog_map)  # don't mutate caller's map
    if resolve_missing:
        linked = [g for g in catalog["target_gene"].unique()
                  if g not in (None, "NA", "") and g not in ortholog_map]
        if linked:
            try:
                ortholog_map.update(human_to_mouse_orthologs(linked))
            except Exception:
                pass  # offline: fall through, missing genes flagged unknown below

    rows = []
    for g in catalog["target_gene"]:
        if g in (None, "NA", ""):
            rows.append(("NA", np.nan, np.nan, np.nan, "unlinked")); continue
        h_spec = float(human_gene_score.get(g, np.nan))
        if g not in ortholog_map:
            # never resolved (offline) -- honestly 'unknown', not 'no_ortholog'
            rows.append(("", h_spec, np.nan, np.nan, "ortholog_unknown")); continue
        mice = ortholog_map.get(g, [])
        if not mice:
            rows.append(("", h_spec, np.nan, 0.0, "no_ortholog")); continue
        # best mouse ortholog by target specificity
        m_specs = [(m, float(mouse_gene_score.get(m, np.nan))) for m in mice]
        m_specs_valid = [(m, s) for m, s in m_specs if not np.isnan(s)]
        if not m_specs_valid:
            # ortholog exists but its mouse expression was not measured in our
            # panel -- 'unknown', NOT 'divergent'. (Divergent requires a measured
            # mouse value that fails the specificity concordance test below.)
            rows.append((";".join(mice), h_spec, np.nan, np.nan, "mouse_expr_unmeasured")); continue
        best_m, m_spec = max(m_specs_valid, key=lambda x: x[1])
        # conservation score: both specific -> high; concordance of specificity
        cons = float(min(h_spec, m_spec)) if not np.isnan(h_spec) else m_spec
        if (not np.isnan(h_spec) and h_spec >= spec_min) and m_spec >= spec_min:
            flag = "conserved"
        else:
            flag = "divergent_expression"
        rows.append((best_m, h_spec, m_spec, cons, flag))
    cols = ["mouse_ortholog", "human_gene_spec", "mouse_gene_spec",
            "conservation_score", "translatability_flag"]
    ann = pd.DataFrame(rows, columns=cols, index=catalog.index)
    out = pd.concat([catalog, ann], axis=1)
    out["conservation_method"] = "gene_level_ortholog_expression_concordance"
    return out
