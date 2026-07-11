"""
Motif-enrichment layer ("T" / TFBS-enriched strategy) for the AAV pipeline.

For a target cell type, score each candidate element by how many of the target
cell type's identity-TF motifs it contains (PWM hits above a log-odds
threshold), then test whether the shortlist is enriched for those motifs
relative to a background sample of cCREs (hypergeometric test).

Cell-type knowledge lives entirely in tf_config.py (a data dict). This module
is generic: give it a TF list + PFMs and any sequences, it scans them.

JASPAR PFMs are open-access (CC-BY). Scanning uses Biopython's PSSM search.
"""
from __future__ import annotations
import os
import numpy as np
import pandas as pd
from Bio import motifs as biomotifs
from Bio.Seq import Seq


def load_pwm(pfm_path, pseudocount=0.5):
    """Load a JASPAR PFM and return a Biopython PSSM (log-odds, uniform bg)."""
    with open(pfm_path) as fh:
        m = biomotifs.read(fh, "jaspar")
    m.pseudocounts = pseudocount
    m.background = {"A": .25, "C": .25, "G": .25, "T": .25}
    return m.pssm


def scan_seq(pssm, seq, thresh=None):
    """Count PWM hits (both strands) above `thresh` log-odds in one sequence.
    Default threshold = PSSM.max * 0.8 (a standard 'strong site' heuristic)."""
    if not seq or len(seq) < pssm.length:
        return 0
    s = seq.upper().replace("N", "A")
    if thresh is None:
        thresh = pssm.max * 0.8
    n = 0
    for _pos, score in pssm.search(Seq(s), threshold=thresh, both=True):
        n += 1
    return n


def motif_scores(seqs, tfs, pfm_dir, tf_pfm):
    """Return a DataFrame: one row per sequence, one column of hit-counts per TF,
    plus 'motif_total' (sum) and 'motif_tf_breadth' (# distinct TFs with >=1 hit).
    `seqs` is a list of DNA strings (candidate elements)."""
    pssms = {}
    for tf in tfs:
        fn = tf_pfm.get(tf)
        if fn and os.path.exists(os.path.join(pfm_dir, fn)):
            pssms[tf] = load_pwm(os.path.join(pfm_dir, fn))
    out = {tf: [] for tf in pssms}
    for s in seqs:
        for tf, pssm in pssms.items():
            out[tf].append(scan_seq(pssm, s))
    df = pd.DataFrame(out, index=range(len(seqs)))
    if df.shape[1]:
        df["motif_total"] = df.sum(axis=1)
        df["motif_tf_breadth"] = (df[list(pssms)] > 0).sum(axis=1)
    else:
        df["motif_total"] = 0
        df["motif_tf_breadth"] = 0
    return df


def hypergeom_enrichment(shortlist_hits, background_hits):
    """Hypergeometric test: are motif-positive elements over-represented in the
    shortlist vs a background sample? Returns (odds-ratio-like, p-value)."""
    from scipy.stats import hypergeom
    K = int((np.asarray(background_hits) > 0).sum())          # bg positives
    N = len(background_hits)                                   # bg total
    k = int((np.asarray(shortlist_hits) > 0).sum())           # shortlist positives
    n = len(shortlist_hits)                                    # shortlist total
    if N == 0 or n == 0:
        return (np.nan, np.nan)
    # P(X >= k) under hypergeometric
    p = hypergeom.sf(k - 1, N, K, n)
    frac_short = k / n if n else 0
    frac_bg = K / N if N else 0
    enr = (frac_short / frac_bg) if frac_bg > 0 else np.nan
    return (enr, p)
