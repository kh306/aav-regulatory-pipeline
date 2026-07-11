"""
Reusable pipeline for tissue-specific AAV regulatory-element candidate design.

ONE entry point -- run_pipeline(cell_type, ...) -- applied UNCHANGED across cell
types. There is deliberately no per-tissue branching anywhere in this module:
`cell_type` is only ever an argument used to select a column of a shared
accessibility matrix. Running a new tissue = calling run_pipeline with a
different string.

Funnel stages:
  0. union cCRE catalogue (CATLAS, ~1.14M elements x 222 human cell types)
  1. ACCESSIBILITY  : keep cCREs open in the target cell type ("H" strategy)
  2. SPECIFICITY    : Tau tissue-specificity index across all 222 cell types
  3. CANDIDATE BUILD: merge neighbouring co-accessible cCREs into blocks
  4. SIZE           : AAV packaging cap on candidate block length
  (motif + cross-species layers live in separate modules and plug in on top)

Data backbone (built once by build_backbone.py, cached to disk):
  catlas_ccre_coords.parquet        chrom,start,end,ccre_id,width
  catlas_accessibility_matrix.npz   sparse binary  cCRE x cell_type
  catlas_celltypes.csv              column index -> cell-type name
  gene_tss.tsv                      protein-coding gene TSS (GENCODE v44)

All inputs are open-access: CATLAS (Zhang et al. 2021, Cell), GENCODE, hg38.
"""
from __future__ import annotations
import os, time, urllib.parse, urllib.request, json
import numpy as np
import pandas as pd
import scipy.sparse as sp

# -----------------------------------------------------------------------------
# Backbone loading (cached artefacts; no per-tissue logic)
# -----------------------------------------------------------------------------
class Backbone:
    def __init__(self, data_dir="data"):
        self.dir = data_dir
        self.coords = pd.read_parquet(os.path.join(data_dir, "catlas_ccre_coords.parquet"))
        self.M = sp.load_npz(os.path.join(data_dir, "catlas_accessibility_matrix.npz")).tocsc()
        ct = pd.read_csv(os.path.join(data_dir, "catlas_celltypes.csv"), index_col=0)
        self.celltypes = list(ct["cell_type"])
        self.ct_index = {c: i for i, c in enumerate(self.celltypes)}
        self.breadth = np.asarray(self.M.sum(axis=1)).ravel().astype(np.int32)
        self.n_ct = len(self.celltypes)
        # Tau on binary data reduces to (N - breadth) / (N - 1); precompute once.
        self.tau_all = (self.n_ct - self.breadth) / (self.n_ct - 1)

    def resolve(self, cell_type):
        """Map a requested cell_type to matrix column(s). Substring, case-insensitive.
        Returns (list_of_column_indices, list_of_matched_names)."""
        q = cell_type.lower().replace(" ", "_")
        hits = [i for i, c in enumerate(self.celltypes) if q in c.lower()]
        if not hits:
            raise ValueError(f"cell_type '{cell_type}' matched no CATLAS column. "
                             f"Examples: {self.celltypes[:5]}")
        return hits, [self.celltypes[i] for i in hits]


# -----------------------------------------------------------------------------
# Gene linking (pluggable). Default = nearest-TSS within a window.
# ABC-score linking can be swapped in here later without touching run_pipeline.
# -----------------------------------------------------------------------------
class NearestTSSLinker:
    method_name = "nearest_TSS_gencode_v44"
    def __init__(self, tss_path, max_dist=1_000_000):
        t = pd.read_csv(tss_path, sep="\t")
        self.max_dist = max_dist
        self.by_chrom = {c: d.sort_values("tss").reset_index(drop=True)
                         for c, d in t.groupby("chrom")}
    def link(self, chrom, start, end):
        d = self.by_chrom.get(chrom)
        if d is None:
            return ("NA", np.nan)
        mid = (start + end) // 2
        i = int(np.searchsorted(d["tss"].values, mid))
        best, bestd = "NA", np.inf
        for k in (i - 1, i):
            if 0 <= k < len(d):
                dist = abs(int(d["tss"].values[k]) - mid)
                if dist < bestd:
                    bestd, best = dist, d["gene_name"].values[k]
        if bestd > self.max_dist:
            return ("NA", float(bestd))
        return (best, float(bestd))


# -----------------------------------------------------------------------------
# hg38 sequence fetch (UCSC REST API) -- only for the final shortlist.
# -----------------------------------------------------------------------------
def fetch_hg38(chrom, start, end, retries=3):
    url = ("https://api.genome.ucsc.edu/getData/sequence?"
           + urllib.parse.urlencode({"genome": "hg38", "chrom": chrom,
                                     "start": int(start), "end": int(end)}))
    for a in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                return json.loads(r.read().decode())["dna"].upper()
        except Exception:
            if a == retries - 1:
                return ""
            time.sleep(1.5)
    return ""


# -----------------------------------------------------------------------------
# THE reusable entry point
# -----------------------------------------------------------------------------
def run_pipeline(cell_type, backbone, linker=None, gene_score=None,
                 motif_scorer=None, *,
                 size_cap=800, top_n=25, merge_gap=500,
                 tau_min=0.90, fetch_seqs=True, verbose=True):
    """Run the accessibility->specificity->candidate->size funnel for one cell type.

    Parameters
    ----------
    cell_type : str        target cell-type name (substring of a CATLAS column)
    backbone  : Backbone   shared, pre-loaded data backbone
    linker    : gene linker with .link(chrom,start,end) and .method_name
    size_cap  : int        max candidate block length in bp (AAV packaging cost)
    top_n     : int        size of the ranked shortlist
    merge_gap : int        merge target-accessible cCREs within this gap (bp)
    tau_min   : float      minimum Tau specificity to enter candidate building

    Returns
    -------
    dict with 'catalog' (DataFrame), 'funnel' (per-stage counts),
    'candidates_all' (all size-capped candidates), 'matched_columns'.
    """
    bb = backbone
    cols, names = bb.resolve(cell_type)
    if verbose:
        print(f"[{cell_type}] matched CATLAS columns: {names}")

    # Stage 0: universe
    n0 = bb.coords.shape[0]

    # Stage 1: accessibility -- open in ANY matched target column
    target_mask = np.asarray(bb.M[:, cols].sum(axis=1)).ravel() > 0
    acc_idx = np.where(target_mask)[0]
    n1 = acc_idx.size

    # Stage 2: specificity (Tau across all 222 cell types), restricted to accessible
    tau = bb.tau_all[acc_idx]
    keep = tau >= tau_min
    spec_idx = acc_idx[keep]
    spec_tau = tau[keep]
    n2 = spec_idx.size

    # Stage 3: candidate build -- merge neighbouring specific cCREs into blocks
    sub = bb.coords.iloc[spec_idx].copy()
    sub["tau"] = spec_tau
    sub = sub.sort_values(["chrom", "start"]).reset_index(drop=True)
    blocks, cur = [], None
    for r in sub.itertuples():
        if cur is None:
            cur = dict(chrom=r.chrom, start=int(r.start), end=int(r.end),
                       taus=[r.tau], ids=[r.ccre_id])
        elif r.chrom == cur["chrom"] and int(r.start) - cur["end"] <= merge_gap:
            cur["end"] = int(r.end); cur["taus"].append(r.tau); cur["ids"].append(r.ccre_id)
        else:
            blocks.append(cur)
            cur = dict(chrom=r.chrom, start=int(r.start), end=int(r.end),
                       taus=[r.tau], ids=[r.ccre_id])
    if cur:
        blocks.append(cur)
    cand = pd.DataFrame([{
        "chrom": b["chrom"], "start": b["start"], "end": b["end"],
        "length_bp": b["end"] - b["start"],
        "n_ccre": len(b["ids"]),
        "tau_specificity": float(np.mean(b["taus"])),
        "ccre_ids": ",".join(b["ids"]),
    } for b in blocks])
    n3 = cand.shape[0]

    # Stage 4: AAV size cap
    cand = cand[cand["length_bp"] <= size_cap].reset_index(drop=True)
    n4 = cand.shape[0]

    # Gene linking for ALL size-capped candidates (pluggable linker; fast binary search)
    if linker is not None and len(cand):
        genes, dists = zip(*[linker.link(r.chrom, r.start, r.end)
                             for r in cand.itertuples()])
        cand["target_gene"] = list(genes)
        cand["gene_dist_bp"] = list(dists)
        cand["gene_link_method"] = linker.method_name
    else:
        cand["target_gene"] = "NA"; cand["gene_dist_bp"] = np.nan
        cand["gene_link_method"] = "none"

    # Gene-expression evidence: how target-SPECIFIC is the linked gene's expression?
    # gene_score maps gene_symbol -> [0,1] specificity in the target cell type
    # (computed from a scRNA panel outside this module; passed in so logic stays
    # cell-type agnostic). Genes absent from the panel get a neutral 0.5.
    if gene_score is not None and len(cand):
        gs = cand["target_gene"].map(lambda g: gene_score.get(g, np.nan))
        cand["gene_expr_specificity"] = gs.fillna(0.5)
    else:
        cand["gene_expr_specificity"] = np.nan

    # Composite score. Chromatin specificity (Tau) is the backbone; when gene
    # expression evidence is available it re-weights among the many high-Tau ties,
    # surfacing candidates whose target gene is expressed-specifically in the tissue.
    cand["specificity_per_bp"] = cand["tau_specificity"] / cand["length_bp"]
    if gene_score is not None:
        cand["composite_score"] = cand["tau_specificity"] * (0.5 + 0.5 * cand["gene_expr_specificity"])
        sort_keys = ["composite_score", "tau_specificity", "specificity_per_bp"]
    else:
        cand["composite_score"] = cand["tau_specificity"]
        sort_keys = ["tau_specificity", "specificity_per_bp"]
    cand = cand.sort_values(sort_keys, ascending=False).reset_index(drop=True)

    # Take a slightly wider slice so the optional motif layer can re-rank within it
    pool = cand.head(max(top_n * 2, top_n)).copy()
    if fetch_seqs and len(pool):
        pool["sequence"] = [fetch_hg38(r.chrom, r.start, r.end)
                            for r in pool.itertuples()]

    # Optional motif layer: reward candidates carrying target identity-TF motifs.
    # Kept as a plug-in (motif_scorer is a callable seqs->DataFrame); if absent,
    # the composite score is unchanged and motif columns are NaN.
    if motif_scorer is not None and fetch_seqs and len(pool) and "sequence" in pool:
        mdf = motif_scorer(pool["sequence"].tolist())
        mdf.index = pool.index
        pool = pd.concat([pool, mdf], axis=1)
        mt = pool["motif_total"].astype(float)
        motif_norm = mt / mt.max() if mt.max() > 0 else mt * 0
        pool["motif_score"] = motif_norm
        # blend motif signal into the composite (small weight; chromatin+gene lead)
        pool["composite_score"] = pool["composite_score"] * (0.8 + 0.2 * motif_norm)
        pool = pool.sort_values("composite_score", ascending=False)
    else:
        pool["motif_score"] = np.nan

    shortlist = pool.head(top_n).copy().reset_index(drop=True)
    shortlist["composite_rank"] = np.arange(1, len(shortlist) + 1)
    shortlist.insert(0, "cell_type", cell_type)

    funnel = pd.DataFrame({
        "stage": ["0_union_cCRE", "1_accessible_in_target",
                  f"2_specific_tau>={tau_min}", "3_candidate_blocks",
                  f"4_size<={size_cap}bp", "5_shortlist"],
        "n": [n0, n1, n2, n3, n4, len(shortlist)],
    })
    funnel["fold_reduction_vs_universe"] = n0 / funnel["n"].clip(lower=1)
    if verbose:
        print(funnel.to_string(index=False))
    return {"catalog": shortlist, "funnel": funnel, "matched_columns": names,
            "candidates_all": cand}
