"""
Validation controls for the AAV regulatory-element pipeline.

A funnel that finds the right things is only half the evidence; it must also
REJECT the wrong things. This module runs both checks for any cell type:

  POSITIVE control -- do known identity/marker genes for the cell type appear
      among the gene-linked shortlist candidates? (e.g. INS/IAPP for beta,
      OLFM4/APOA4/RBP2 for enterocyte). Passing = the funnel surfaces real
      cell-type biology.

  NEGATIVE control -- does a ubiquitously-expressed housekeeping locus score
      LOW tissue-specificity and get filtered OUT before the shortlist? We take
      cCREs at classic housekeeping genes (ACTB, GAPDH, B2M, ...) and confirm
      their Tau specificity in the target is low and they do not reach the
      shortlist. Passing = the funnel discards non-specific elements.

Marker/housekeeping sets are data (dicts below), editable per cell type. This
keeps the module cell-type-agnostic: it validates whatever cell type it is given.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

# Known positive-control marker genes per cell type (substring-keyed).
POSITIVE_MARKERS = {
    "beta":       ["INS", "IAPP", "PDX1", "NKX6-1", "MAFA", "PCSK1", "CHGA", "SCGN"],
    "enterocyte": ["OLFM4", "APOA4", "RBP2", "FABP6", "PIGR", "ANPEP", "VIL1", "CDX2"],
}

# Housekeeping / ubiquitous loci for the negative control (species: human hg38).
HOUSEKEEPING = ["ACTB", "GAPDH", "B2M", "TBP", "PGK1", "RPL13A"]
HOUSEKEEPING_TSS = {  # hg38 TSS (approx), used to locate their cCREs
    "ACTB":   ("chr7", 5530601),
    "GAPDH":  ("chr12", 6534517),
    "B2M":    ("chr15", 44711477),
    "TBP":    ("chr6", 170554333),
    "PGK1":   ("chrX", 78104248),
    "RPL13A": ("chr19", 49487603),
}


def markers_for(cell_type):
    q = cell_type.lower()
    for k, v in POSITIVE_MARKERS.items():
        if k in q:
            return v
    return []


def positive_control(catalog, cell_type):
    """Which marker genes appear among gene-linked shortlist candidates?"""
    markers = markers_for(cell_type)
    linked = set(str(g) for g in catalog.get("target_gene", []) if g not in (None, "NA", ""))
    found = [m for m in markers if m in linked]
    return {
        "markers_expected": markers,
        "markers_found": found,
        "n_found": len(found),
        "n_expected": len(markers),
        "pass": len(found) >= 1,
        "top_marker_rank": int(catalog.loc[catalog["target_gene"].isin(found), "composite_rank"].min())
                           if found else None,
    }


def negative_control(backbone, cell_type, tau_specific_threshold=0.90):
    """Does the PROMOTER cCRE of each housekeeping gene score LOW target
    specificity (low Tau / broad accessibility) -- i.e. is it correctly NOT
    treated as a tissue-specific candidate?

    We take the cCRE nearest each housekeeping TSS (its promoter). Housekeeping
    promoters are open in most cell types, so their Tau should be well below the
    specificity threshold used to build candidates. Distal enhancers near the
    same gene can be cell-type-restricted, so we test the promoter specifically,
    not the max over a window.
    """
    bb = backbone
    cols, _ = bb.resolve(cell_type)
    c = bb.coords
    rows = []
    for g, (chrom, tss) in HOUSEKEEPING_TSS.items():
        near = c[(c["chrom"] == chrom) & (c["start"] > tss - 3000) & (c["end"] < tss + 3000)]
        if not len(near):
            rows.append((g, np.nan, np.nan, np.nan, np.nan)); continue
        mids = ((near["start"] + near["end"]) // 2).values
        prom_row = near.iloc[int(np.argmin(np.abs(mids - tss)))]
        i = prom_row.name
        acc = bool(np.asarray(bb.M[i, cols].sum()).ravel()[0] > 0)
        rows.append((g, acc, int(bb.breadth[i]), float(bb.tau_all[i]), int(prom_row["start"])))
    df = pd.DataFrame(rows, columns=["gene", "promoter_accessible_in_target",
                                     "promoter_breadth", "promoter_tau", "promoter_start"])
    below = (df["promoter_tau"] < tau_specific_threshold)
    return {
        "table": df,
        "median_promoter_tau": float(df["promoter_tau"].median()),
        "median_promoter_breadth": float(df["promoter_breadth"].median()),
        "frac_below_specific_threshold": float(below.mean()),
        "threshold": tau_specific_threshold,
        "pass": bool(below.mean() >= 0.5),
    }


def negative_control_expression(gene_score, housekeeping=None, spec_max=0.30):
    """Housekeeping genes should score LOW target expression-specificity.
    Complements the chromatin-level negative control with an expression-level one.
    """
    hk = housekeeping or HOUSEKEEPING
    vals = {g: float(gene_score.get(g, np.nan)) for g in hk}
    s = pd.Series(vals)
    below = (s < spec_max)
    return {
        "table": s.rename("expr_specificity").reset_index().rename(columns={"index": "gene"}),
        "median_expr_specificity": float(s.median()),
        "frac_below": float(below.mean()),
        "threshold": spec_max,
        "pass": bool(below.mean() >= 0.5),
    }
