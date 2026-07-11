# Submission — Built with Claude: Life Sciences

**Project:** A reusable, cross-species-validated pipeline for tissue-specific AAV regulatory-element design
**Track:** Research

---

## One-paragraph abstract (submission form)

AAV gene therapies need regulatory elements (promoters/enhancers) that restrict
transgene expression to a target cell type — for safety, dose, and to spare
off-target tissue. Designing these is currently a bespoke, per-tissue effort, and
elements optimised in mouse routinely fail in human because the underlying
regulatory biology differs across species. We built a **single reusable
pipeline** that mines a human single-cell chromatin-accessibility atlas (CATLAS,
1.14M candidate elements × 222 cell types) for tissue-specific regulatory
elements, ranks them by chromatin specificity, target-gene expression, motif
content, and AAV packaging size, and — the differentiating feature — **flags each
candidate for cross-species translatability** by asking whether the target
biology it serves is conserved between human and mouse (CELLxGENE Census +
Ensembl orthologs). The same code, run unchanged (only a cell-type argument
changes), produced validated "parts catalogs" for pancreatic beta cells and
intestinal enterocytes. Positive controls recover the correct markers (INS/IAPP;
OLFM4/APOA4/PIGR), negative controls reject housekeeping loci, motif enrichment
is 4.7× (p=1.6×10⁻¹³), and the cross-species filter caught a human-specific
enterocyte target (PHGR1) that a mouse-only workflow would have wrongly trusted.
The contribution is the reusable, translation-aware **method** — a
synthetic-biology "parts registry" approach to gene-therapy regulatory design.

---

## What's novel

- **Reusability as the deliverable.** Not one promoter — a pipeline with *zero
  per-tissue code paths*, proven by running two tissues off the same code with
  only an argument change (`catalogs/reusability_proof.json`).
- **Cross-species translatability as a first-class filter.** Most element-design
  work stops at "specific in human." We add an explicit, honestly-scoped
  (gene-level) conservation flag that de-risks the most expensive failure mode in
  the field before any wet-lab spend.

## What it produces (per cell type, one command)

- A ranked, size-constrained, cross-species-annotated **parts catalog** (CSV).
- **Specificity + funnel figures**, and a cross-tissue reusability/cost funnel.
- Positive + negative **control** reports.

## Honest limits (stated in README)

Accessible ≠ active promoter; a cCRE is not an autonomous promoter; specificity
is measured only against 222 cell types; gene linking is nearest-TSS (ABC scores
were unavailable on the public mirror); conservation is **gene-level, not
enhancer-sequence-level**. Each is a labelled column or a Future Work item, not a
hidden assumption.

## Data & license

All inputs open-access (CATLAS, JASPAR CC-BY, GENCODE, hg38, CELLxGENE Census
CC-BY, Ensembl); pipeline code MIT.

---

## 60-second demo script (for judges / live walkthrough)

```bash
pip install -r requirements.txt

# 1. The whole pipeline for beta cells — one command
python run.py "Pancreatic_Beta_Cell"
#    → catalogs/Pancreatic_Beta_Cell_catalog.csv  (top parts at the INS locus)
#    → prints: positive control [INS, IAPP] pass; negative control pass

# 2. The SAME command, different tissue — zero code changes
python run.py "Small_Intestinal_Enterocyte"
#    → top parts at OLFM4/REG1A/APOA4/PIGR; PHGR1 flagged no_ortholog (high-risk)

# 3. The headline: open the enterocyte catalog and sort by translatability_flag
#    conserved  = translation-confident   |   no_ortholog = human-specific, de-prioritise
```

The reusability claim is checkable in one line: `diff` the two invocations —
only the argument string differs.
