"""
Cross-species expression evidence for the AAV regulatory-element pipeline.

Two jobs, both cell-type-agnostic (cell type is only ever an argument):
  1. gene_score_magspec(panel, target_label)
        -> per-gene [0,1] score combining magnitude and target-concentration
           of expression; used to re-rank high-Tau chromatin candidates so the
           ones tied to genes actually expressed-specifically in the target
           cell type rise to the top.
  2. build_expr_panel(census, organism, labels)
        -> gene x cell-type mean-expression table streamed from CELLxGENE
           Census (open-access, CC-BY). Cached to parquet; slow to build
           (per-label S3 streaming), fast to reuse.

The cross-species translatability filter (cross_species.py) reuses job 1 on
the SAME target cell type in mouse to test whether the target-gene expression
signal is conserved or mouse-specific.
"""
from __future__ import annotations
import numpy as np
import pandas as pd


def gene_score_magspec(panel: pd.DataFrame, target_col: str) -> pd.Series:
    """Per-gene target-cell-type expression score in [0,1].

    score = log1p(mean_expr_in_target) * (mean_expr_in_target / sum across panel)
    i.e. a gene must be BOTH highly expressed AND concentrated in the target to
    score high. Normalised to max=1 for interpretability.
    """
    X = panel.values.astype(float)
    ti = list(panel.columns).index(target_col)
    tgt = X[:, ti]
    rowsum = X.sum(axis=1).astype(float)
    rowsum[rowsum == 0] = np.nan
    conc = tgt / rowsum
    mag = np.log1p(tgt)
    s = pd.Series(mag * conc, index=panel.index).fillna(0.0)
    if s.max() > 0:
        s = s / s.max()
    return s


def build_expr_panel(census, organism, labels, cap=4000, seed=0):
    """Stream mean expression per cell-type label from CELLxGENE Census.

    Returns a genes x labels DataFrame (missing genes filled 0). Streaming is
    slow (~1-3 min/label over S3); build once and cache to parquet.
    """
    import cellxgene_census
    cols = {}
    for lb in labels:
        a = cellxgene_census.get_anndata(
            census, organism=organism,
            obs_value_filter=f'cell_type == "{lb}" and is_primary_data == True',
            var_column_names=["feature_name"], obs_column_names=["cell_type"])
        if a.n_obs == 0:
            continue
        if a.n_obs > cap:
            idx = np.random.RandomState(seed).choice(a.n_obs, cap, replace=False)
            a = a[idx].copy()
        m = np.asarray(a.X.mean(axis=0)).ravel()
        cols[lb] = pd.Series(m, index=a.var["feature_name"].values).groupby(level=0).mean()
    return pd.DataFrame(cols).fillna(0.0)
