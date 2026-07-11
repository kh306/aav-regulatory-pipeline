"""
Editable identity-TF configuration for the motif-enrichment layer.

Maps a target cell type to the transcription factors whose binding motifs we
reward when found in a candidate element. This is the ONLY place cell-type
knowledge is encoded, and it is data (a dict), not branching logic -- the
pipeline reads the entry for whatever `cell_type` it is given. Adding a new
tissue = adding a dict entry + its JASPAR PFM files, no code change.

Matching is by substring against the CATLAS column name (case-insensitive),
so 'Pancreatic_Beta_Cell_1/2' both map via the 'beta' key.

JASPAR matrix IDs (latest versions, CORE vertebrates, open-access CC-BY):
  PDX1 MA0132.3, NKX6-1 MA0674.2, PAX6 MA0069.1, NEUROD1 MA1109.2,
  CDX2 MA0465.3, HNF4A MA1494.2, KLF5 MA0599.1
"""
IDENTITY_TFS = {
    "beta":       ["PDX1", "NKX6-1", "PAX6", "NEUROD1"],   # pancreatic beta cell
    "enterocyte": ["CDX2", "HNF4A", "KLF5"],               # intestinal enterocyte
}

# TF -> JASPAR PFM filename (under data/motifs/)
TF_PFM = {
    "PDX1":    "PDX1_MA0132.3.jaspar",
    "NKX6-1":  "NKX6-1_MA0674.2.jaspar",
    "PAX6":    "PAX6_MA0069.1.jaspar",
    "NEUROD1": "NEUROD1_MA1109.2.jaspar",
    "CDX2":    "CDX2_MA0465.3.jaspar",
    "HNF4A":   "HNF4A_MA1494.2.jaspar",
    "KLF5":    "KLF5_MA0599.1.jaspar",
}


def tfs_for(cell_type):
    """Return the identity-TF list for a cell type by substring match, or []."""
    q = cell_type.lower()
    for key, tfs in IDENTITY_TFS.items():
        if key in q:
            return tfs
    return []
