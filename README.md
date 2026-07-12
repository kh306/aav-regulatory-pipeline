# A reusable, cross-species-validated pipeline for tissue-specific AAV regulatory-element design

**One pipeline. Any cell type. One argument changes.**

This project builds a single, reusable computational pipeline that mines a
human single-cell chromatin-accessibility atlas for **tissue-specific regulatory
elements** (candidate promoters/enhancers) suitable for cell-type-restricted
AAV gene-therapy cassettes — and then **flags each candidate for cross-species
translatability** before any wet-lab spend. The contribution is the *reusable,
translation-aware method*, not any single element: the same code, run unchanged,
produces a ranked "parts catalog" for pancreatic **beta cells** and intestinal
**enterocytes**, and would extend to a new cell type by adding a config entry.

---

## Why this matters (real-world framing)

- **Tissue-restricted expression is a safety and dose lever.** Confining a
  therapeutic transgene to the target cell type lowers off-target toxicity and
  immune exposure and reduces the effective dose — e.g. beta-cell-restricted
  expression for diabetes gene therapy that spares exocrine pancreas and liver.
- **Cross-species failure is the field's most expensive failure mode.** AAV
  constructs optimised in mouse routinely fail in NHP/human because the
  underlying regulatory biology differs across species. Flagging likely
  non-translating candidates *computationally, up front* targets that failure
  mode directly and reduces wasted animal cohorts (a 3Rs benefit).
- **A "parts catalog" is the synthetic-biology registry model applied to gene
  therapy** — reusable, characterised regulatory parts instead of bespoke
  one-offs.

---

## The pipeline (7 stages, one entry point)

```
python run.py "Pancreatic_Beta_Cell"
python run.py "Small_Intestinal_Enterocyte"
```

| # | Stage | What it does |
|---|-------|--------------|
| 0 | **Union cCRE catalogue** | 1,143,424 candidate cis-regulatory elements × 222 human cell types (CATLAS) |
| 1 | **Accessibility** | keep elements open in the target cell type |
| 2 | **Specificity** | Tau tissue-specificity index across all 222 cell types (τ=1 → open only in target) |
| 3 | **Candidate build** | merge neighbouring co-accessible elements into enhancer blocks |
| 4 | **AAV size cap** | enforce a packaging-cost budget (default ≤800 bp) |
| 5 | **Ranking** | composite of chromatin specificity × target-gene expression evidence (+ optional motif enrichment) |
| 6 | **Cross-species translatability** | link → mouse ortholog → mouse expression → conservation flag |

**Reusability is structural, not claimed.** `run_pipeline(cell_type, ...)` has
no per-tissue branching anywhere — `cell_type` is only ever a key that selects a
column of one shared accessibility matrix. Cell-type knowledge lives entirely in
editable *data* (`tf_config.py` identity TFs; `controls.py` marker/housekeeping
sets). The beta and enterocyte runs above differ by **zero lines of code** — only
the argument string (`catalogs/reusability_proof.json`).

---

## Results

**Funnel (same thresholds, both tissues).** The universe of ~1.14M elements
narrows to a size-capped specific set of ~65k (beta, 18×) / ~32k (enterocyte,
36×), from which a ranked shortlist is taken. See `figures/reusability_funnel.png`.

**Positive controls — the funnel finds the right biology.**
- Beta: top-ranked candidates sit at the **INS (insulin)** locus, followed by
  **IAPP** — the two beta-cell hormones.
- Enterocyte: top parts link to **OLFM4, REG1A, APOA4, PIGR, FABP6** — canonical
  enterocyte/brush-border markers.

**Negative controls — the funnel rejects the wrong things.** *(added when
incorporating reviewer feedback.)* Housekeeping promoters (ACTB, GAPDH, B2M,
TBP, PGK1, RPL13A) are accessible in the target but score τ≈0.02 (open in
~218/222 cell types) and expression-specificity ≈0.04 — correctly filtered out
long before the shortlist.

**Motif enrichment.** The beta shortlist is **4.7× enriched** for beta
identity-TF motifs (PDX1/NKX6-1/PAX6/NEUROD1) versus a matched null of
ubiquitous elements (hypergeometric p = 1.6×10⁻¹³; 72% carry a PDX1 site vs
12% of background).

**Cross-species translatability — the headline feature.** Each shortlisted,
gene-linked candidate is flagged `conserved`, `divergent_expression`,
`no_ortholog`, or `unlinked`. Beta INS→*Ins2* and IAPP→*Iapp* flag **conserved**;
in the enterocyte catalog **PHGR1 flags `no_ortholog`** — a human-specific target
a mouse-only workflow would have wrongly trusted. That single flag is the value
proposition: it removes likely-non-translating parts before wet-lab validation
(`figures/reusability_funnel.png`, panel c).

---

## Honest scope & limitations

This is a **candidate-prioritisation and de-risking** tool, not a promoter
guarantee. The pipeline was self-audited against known AI-overclaim failure
modes; see **`AUDIT_RESPONSE.md`** for the review, the bugs it found, and the
`--strict-links` mode that fixes them (recommended for reporting). Key
limitations, stated plainly:

1. **Accessible ≠ active ≠ a working promoter.** Open chromatin includes
   insulators, silencers, and poised-but-inactive elements. Some shortlisted
   elements will be accessible but not activating.
2. **A cCRE is not an autonomous promoter.** These are enhancer-like elements;
   dropped into a cassette on a minimal promoter, an enhancer may not fire or may
   lose specificity (enhancer–promoter compatibility is real).
3. **Specificity is measured only against CATLAS's 222 cell types.** An element
   that looks target-specific could be active in an unsampled cell type →
   residual off-target risk.
4. **Gene linking is nearest-TSS, not ABC.** The CATLAS ABC (activity-by-contact)
   enhancer→gene scores were not downloadable from the public mirror at build
   time, so we use nearest-TSS linking (GENCODE v44) and record
   `gene_link_method` in every catalog row. The linker is pluggable — ABC can be
   dropped in without touching the pipeline.
5. **Cross-species conservation is GENE-LEVEL, not enhancer-sequence-level.**
   Enhancers turn over rapidly in evolution; a human enhancer often has no mouse
   counterpart even when its target gene is perfectly conserved. Our flag asks
   *"is the target biology this element serves conserved between human and
   mouse?"* — not *"is this enhancer sequence conserved?"*. Every output carries
   a `conservation_method` column saying so. Enhancer-sequence conservation
   (liftOver/phastCons) is future work.

**Future work:** ABC-based gene linking; enhancer-level sequence conservation;
NHP/human cross-species expression (Census also carries macaque, marmoset,
chimp); and DL activity prediction (e.g. Borzoi) to address limitation #1.

---

## Data sources (all open-access)

| Data | Source | License |
|------|--------|---------|
| Single-cell chromatin accessibility (cCRE × cell type) | **CATLAS** human tissues (Zhang et al. 2021, *Cell*) | open-access |
| TF binding motifs (PWMs) | **JASPAR** 2024 CORE vertebrates | CC-BY 4.0 |
| Gene models / TSS | **GENCODE** v44 (hg38) | open |
| Reference genome sequence | **hg38** via UCSC REST | open |
| Cross-species single-cell expression | **CELLxGENE Census** (stable 2025-11-08) | CC-BY 4.0 |
| Human↔mouse orthologs | **Ensembl** REST homology | open |

Pipeline code: **MIT** (see `LICENSE`).

---

## Repository layout

```
aav_pipeline/
├── run.py               # single entry point:  python run.py "<cell_type>"
├── pipeline.py          # accessibility → specificity → candidate → size funnel
├── expression.py        # gene expression-specificity scoring (Census)
├── cross_species.py     # ortholog mapping + gene-level conservation flags
├── motif.py             # JASPAR PWM scanning + enrichment
├── tf_config.py         # editable per-cell-type identity-TF sets  (DATA, not code)
├── controls.py          # positive + negative validation controls
├── build_backbone.py    # rebuild the cached data backbone from public sources
├── data/                # cached backbone + motif PWMs + expression panels
├── catalogs/            # per-cell-type + unified parts catalogs, funnels, controls
├── figures/             # reusability funnel + per-tissue specificity figures
├── requirements.txt
├── LICENSE              # MIT
└── .gitignore
```

## Reproducing

```bash
pip install -r requirements.txt
python build_backbone.py          # downloads CATLAS cCRE BEDs + GENCODE, builds matrix
python run.py "Pancreatic_Beta_Cell"
python run.py "Small_Intestinal_Enterocyte"
```

**Performance note (CELLxGENE Census).** Whole-matrix reads combining a
`var_value_filter` with scattered `obs_coords` are extremely slow over the S3
proxy. The working pattern (used here) is a plain `obs_value_filter` for the
cell type, fetching all genes, then subsampling cells and restricting to the
genes of interest *in memory*. Per-cell-type panels are cached to parquet so this
cost is paid once.
