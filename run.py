#!/usr/bin/env python
"""
Single entry point for the reusable AAV regulatory-element pipeline.

    python run.py "Pancreatic_Beta_Cell"
    python run.py "Small_Intestinal_Enterocyte"

The SAME command runs any cell type present in the CATLAS atlas -- the only
argument that changes between tissues is the cell-type name. There is no
per-tissue code path: cell-type knowledge lives entirely in editable config
(tf_config.py identity TFs, controls.py marker/housekeeping sets). Adding a new
tissue = a config entry + its JASPAR PFMs, no pipeline edits.

Outputs (per cell type, written to catalogs/):
    <ct>_catalog.csv   ranked, size-capped, cross-species-annotated parts
    <ct>_funnel.csv    per-stage candidate counts + fold reductions
    <ct>_controls.json positive + negative control results

Requires the cached data backbone in data/ (see build_backbone.py) and the
human/mouse expression panels + ortholog map (built once via Census/Ensembl).
All inputs open-access: CATLAS (Zhang 2021), JASPAR (CC-BY), GENCODE, hg38,
CELLxGENE Census (CC-BY), Ensembl.
"""
from __future__ import annotations
import os, sys, json, argparse
import numpy as np
import pandas as pd

import pipeline as P
import expression as EX
import cross_species as CS
import motif as MO
import tf_config as TC
import controls as CTRL

DATA = os.path.join(os.path.dirname(__file__), "data")
CATA = os.path.join(os.path.dirname(__file__), "catalogs")

# Map a CATLAS cell-type name to the Census scRNA label used for expression.
# (Data, not logic -- extend for new tissues.)
CENSUS_LABEL = {
    "beta":       "type B pancreatic cell",
    "enterocyte": "enterocyte of epithelium of small intestine",
}


def census_label_for(cell_type):
    q = cell_type.lower()
    for k, v in CENSUS_LABEL.items():
        if k in q:
            return v
    return None


def make_motif_scorer(cell_type):
    tfs = TC.tfs_for(cell_type)
    pfm_dir = os.path.join(DATA, "motifs")
    return lambda seqs: MO.motif_scores(seqs, tfs, pfm_dir, TC.TF_PFM)


def run(cell_type, size_cap=800, top_n=25, tau_min=0.90, fetch_seqs=True,
        link_max_dist=1_000_000, impute_unlinked=True, require_gene_link=False):
    bb = P.Backbone(data_dir=DATA)
    linker = P.NearestTSSLinker(os.path.join(DATA, "gene_tss.tsv"))

    # Expression evidence (human + mouse panels are prebuilt & cached)
    hpanel = pd.read_parquet(os.path.join(DATA, "human_expr_panel.parquet"))
    mpanel = pd.read_parquet(os.path.join(DATA, "mouse_expr_panel.parquet"))
    clabel = census_label_for(cell_type)
    human_gs = EX.gene_score_magspec(hpanel, clabel).to_dict() if clabel and clabel in hpanel.columns else None

    # Run the funnel (accessibility -> specificity -> candidate -> size -> motif)
    res = P.run_pipeline(cell_type, bb, linker=linker, gene_score=human_gs,
                         motif_scorer=make_motif_scorer(cell_type),
                         size_cap=size_cap, top_n=top_n, tau_min=tau_min,
                         fetch_seqs=fetch_seqs, verbose=True,
                         link_max_dist=link_max_dist,
                         impute_unlinked=impute_unlinked,
                         require_gene_link=require_gene_link)

    # Cross-species translatability filter (gene-level target-biology conservation)
    catalog = res["catalog"]
    if os.path.exists(os.path.join(DATA, "ortholog_map.json")) and human_gs is not None:
        ortho = json.load(open(os.path.join(DATA, "ortholog_map.json")))
        # mouse gene score for the matching mouse cell-type column
        m_col = None
        for c in mpanel.columns:
            if clabel and (clabel.split()[0].lower() in c.lower() or "enterocyte" in c.lower() and "enterocyte" in clabel):
                m_col = c; break
        if m_col is None and len(mpanel.columns):
            m_col = mpanel.columns[0]
        mouse_gs = EX.gene_score_magspec(mpanel, m_col).to_dict()
        catalog = CS.conservation_flags(catalog, human_gs, mouse_gs, ortho)

    # Controls
    pos = CTRL.positive_control(catalog, cell_type)
    neg = CTRL.negative_control(bb, cell_type)
    neg_e = CTRL.negative_control_expression(human_gs) if human_gs else {"pass": None}
    controls = {"positive": pos,
                "negative_chromatin": {k: v for k, v in neg.items() if k != "table"},
                "negative_expression": {k: v for k, v in neg_e.items() if k != "table"}}

    # Write outputs
    os.makedirs(CATA, exist_ok=True)
    tag = cell_type.replace(" ", "_")
    catalog.to_csv(os.path.join(CATA, f"{tag}_catalog.csv"), index=False)
    res["funnel"].to_csv(os.path.join(CATA, f"{tag}_funnel.csv"), index=False)
    json.dump(controls, open(os.path.join(CATA, f"{tag}_controls.json"), "w"),
              indent=2, default=str)

    print(f"\n[{cell_type}] positive control: {pos['markers_found']} (pass={pos['pass']})")
    print(f"[{cell_type}] negative control (chromatin promoter Tau median): "
          f"{neg['median_promoter_tau']:.3f} (pass={neg['pass']})")
    return {"catalog": catalog, "funnel": res["funnel"], "controls": controls}


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Reusable AAV regulatory-element pipeline")
    ap.add_argument("cell_type", help="CATLAS cell-type name, e.g. 'Pancreatic_Beta_Cell'")
    ap.add_argument("--size-cap", type=int, default=800)
    ap.add_argument("--top-n", type=int, default=25)
    ap.add_argument("--tau-min", type=float, default=0.90)
    ap.add_argument("--no-seqs", action="store_true", help="skip hg38 sequence fetch")
    ap.add_argument("--link-max-dist", type=int, default=1_000_000,
                    help="max element-to-TSS distance (bp) counted as a gene link "
                         "(default 1Mb; use 50000 for a defensible enhancer-range gate)")
    ap.add_argument("--strict-links", action="store_true",
                    help="honest mode: 50kb link gate, no 0.5 imputation, "
                         "gene-linked candidates only in the ranked catalog")
    a = ap.parse_args()
    if a.strict_links:
        run(a.cell_type, size_cap=a.size_cap, top_n=a.top_n, tau_min=a.tau_min,
            fetch_seqs=not a.no_seqs, link_max_dist=50_000,
            impute_unlinked=False, require_gene_link=True)
    else:
        run(a.cell_type, size_cap=a.size_cap, top_n=a.top_n, tau_min=a.tau_min,
            fetch_seqs=not a.no_seqs, link_max_dist=a.link_max_dist)
