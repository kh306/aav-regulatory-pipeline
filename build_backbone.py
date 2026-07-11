#!/usr/bin/env python
"""
Rebuild the cached data backbone from public sources.

Produces (in data/):
  catlas_ccre_coords.parquet        chrom,start,end,ccre_id,width
  catlas_accessibility_matrix.npz   sparse binary cCRE x cell_type
  catlas_celltypes.csv              column index -> cell-type name
  gene_tss.tsv                      protein-coding TSS (GENCODE v44)

The 222 per-cell-type cCRE BED files (~800 MB) are downloaded to a scratch dir
and unioned into a single ~1.14M x 222 binary accessibility matrix. The
expression panels (human/mouse) and ortholog map are built separately (Census /
Ensembl) — see expression.py and cross_species.py; they are cached in data/ too.

All sources open-access: CATLAS (Zhang et al. 2021, Cell); GENCODE v44.
"""
from __future__ import annotations
import os, glob, subprocess, urllib.parse
import numpy as np, pandas as pd, scipy.sparse as sp

DATA = os.path.join(os.path.dirname(__file__), "data")
SCRATCH = os.path.join(DATA, "cCRE_by_cell_type")
CATLAS = ("https://decoder-genetics.wustl.edu/catlasv1/catlas_downloads/"
          "humantissues/cCRE_by_cell_type")
GENCODE = ("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
           "release_44/gencode.v44.basic.annotation.gtf.gz")


def download_beds():
    os.makedirs(SCRATCH, exist_ok=True)
    # fetch the directory listing, then each .bed
    idx = subprocess.run(["curl", "-sSLk", CATLAS + "/"],
                         capture_output=True, text=True).stdout
    import re
    beds = re.findall(r'href="([^"?/][^"]*\.bed)"', idx)
    for i, f in enumerate(sorted(set(beds)), 1):
        out = os.path.join(SCRATCH, f)
        if os.path.exists(out) and os.path.getsize(out) > 0:
            continue
        url = CATLAS + "/" + urllib.parse.quote(f)
        subprocess.run(["curl", "-sSLk", "--max-time", "120", "-o", out, url])
        if i % 50 == 0:
            print(f"  downloaded {i} BEDs")
    return sorted(glob.glob(os.path.join(SCRATCH, "*.bed")))


def build_matrix(files):
    ccre_index, coords, celltypes, per_ct = {}, [], [], []
    for fp in files:
        celltypes.append(os.path.basename(fp)[:-4])
        rows = []
        with open(fp) as fh:
            for line in fh:
                p = line.split("\t")
                if len(p) < 4:
                    continue
                key = (p[0], int(p[1]), int(p[2]))
                idx = ccre_index.get(key)
                if idx is None:
                    idx = len(coords); ccre_index[key] = idx
                    coords.append((p[0], int(p[1]), int(p[2]), p[3]))
                rows.append(idx)
        per_ct.append(np.array(rows, dtype=np.int32))
    n_ccre, n_ct = len(coords), len(celltypes)
    ra = np.concatenate(per_ct)
    ca = np.concatenate([np.full(a.shape[0], j, np.int32) for j, a in enumerate(per_ct)])
    M = sp.csr_matrix((np.ones(ra.shape[0], np.uint8), (ra, ca)), shape=(n_ccre, n_ct))
    M.data[:] = 1
    cdf = pd.DataFrame(coords, columns=["chrom", "start", "end", "ccre_id"])
    cdf["width"] = cdf["end"] - cdf["start"]
    cdf.to_parquet(os.path.join(DATA, "catlas_ccre_coords.parquet"), index=False)
    sp.save_npz(os.path.join(DATA, "catlas_accessibility_matrix.npz"), M)
    pd.Series(celltypes, name="cell_type").to_csv(os.path.join(DATA, "catlas_celltypes.csv"))
    print(f"  matrix: {n_ccre} cCREs x {n_ct} cell types, nnz={M.nnz}")


def build_tss():
    gz = os.path.join(DATA, "gencode.gtf.gz")
    if not os.path.exists(gz):
        subprocess.run(["curl", "-sSL", "--max-time", "180", "-o", gz, GENCODE])
    import gzip, re
    rows = []
    with gzip.open(gz, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            gt = re.search(r'gene_type "([^"]+)"', f[8])
            gn = re.search(r'gene_name "([^"]+)"', f[8])
            if not gt or gt.group(1) != "protein_coding":
                continue
            tss = int(f[3]) if f[6] == "+" else int(f[4])
            rows.append((f[0], tss, f[6], gn.group(1) if gn else "NA", "protein_coding"))
    pd.DataFrame(rows, columns=["chrom", "tss", "strand", "gene_name", "gene_type"]) \
        .to_csv(os.path.join(DATA, "gene_tss.tsv"), sep="\t", index=False)
    print(f"  gene_tss.tsv: {len(rows)} protein-coding genes")


if __name__ == "__main__":
    os.makedirs(DATA, exist_ok=True)
    print("1/3 downloading CATLAS cCRE BEDs ...")
    files = download_beds()
    print(f"2/3 building union accessibility matrix from {len(files)} cell types ...")
    build_matrix(files)
    print("3/3 building GENCODE TSS table ...")
    build_tss()
    print("done. backbone cached in data/")
