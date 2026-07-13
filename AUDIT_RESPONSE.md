# Self-audit and corrections

This pipeline was reviewed against the failure modes documented in Gulluoglu
et al., *Practical Use of Advanced AI Frameworks on Real-Life Scientific
Problems: Three Case Studies* (bioRxiv preprint,
doi:10.64898/2026.06.23.734132; PMC13345091; CC BY) — a study finding that
agentic AI systems tend to **overclaim relative to what their output supports**
and produce **plausible-looking steps that aren't justified**. Only the
abstract was retrievable at review time, so we used its thesis (not its detailed
methods) as the audit lens, re-verifying every load-bearing claim in our own
pipeline against source data. Findings and fixes below.

## What we verified holds (not hallucinated)

- **PHGR1 has no mouse ortholog.** Queried Ensembl directly: PHGR1
  (ENSG00000233041, "proline histidine and glycine rich 1") is a real human
  protein-coding gene, and Ensembl returns **zero** mouse orthologues, while the
  control gene OLFM4 correctly returns one. The headline cross-species finding is
  real.
- **INS is the top beta candidate**, 69 bp from the INS TSS. IAPP candidates are
  13-46 kb out. Credible.

## Problems found, and the fixes

**1. Long-range nearest-TSS "links" were presented as facts, not guesses.** In
the original (permissive) run the top enterocyte candidate was **845 kb** from
OLFM4, and another was 879 kb — nearest-TSS assigns *some* gene to every element
regardless of distance, and the conservation flag then inherited that
assignment. The issue is *not* that enhancers can't act over long ranges — they
demonstrably can (the ZRS enhancer drives SHH from ~1 Mb away, inside an intron
of LMBR1). The issue is that nearest-TSS measures **linear distance, not
physical contact**, so it cannot distinguish a genuine long-range enhancer from
a coincidental nearest neighbour in a gene desert. At 845 kb, "OLFM4" is an
unsupported guess dressed as a fact.
→ **Fix:** added a confident-link **distance gate** (`link_max_dist`, default
1 Mb reproduces old behaviour; `--strict-links` sets 50 kb — the range where
"nearest gene" is defensible on proximity alone). Beyond the gate the gene name
is dropped and the element is `unlinked`; `gene_dist_bp` is always retained for
transparency. This is a **precision/recall trade-off**, not a biological law: it
deliberately drops real long-range enhancers (ZRS-like) to avoid overclaiming.
The correct long-term fix is contact-based linking (Hi-C / ABC scores) in the
pipeline's existing pluggable-linker slot — unavailable here only because the
public ABC-score files on the CATLAS mirror were empty.

**2. The beta catalog was mostly gene-desert noise.** Only 7/25 beta candidates
were gene-linked; the other 18 sat a mean of **1.43 Mb** from any TSS and were
lifted into the shortlist by an imputed constant (see #3).
→ **Fix:** `require_gene_link` (on in `--strict-links`) restricts the ranked
catalog to gene-linked candidates. Beta went from **7/25 -> 25/25** gene-linked,
surfacing real islet genes the padding had masked (GNAS, ERO1B, CPE, GAD2).

**3. Unlinked candidates got an imputed `gene_expr_specificity = 0.5`.** This
"half credit for free" is not a measured quantity and materially reordered the
catalog, floating high-Tau gene deserts above real marker enhancers.
→ **Fix:** `impute_unlinked=False` (on in `--strict-links`) leaves them NaN;
they score on chromatin specificity alone and rank below any gene-linked
candidate of equal Tau.

**4. A false-negative bug in the conservation flag (found while fixing #2).**
`conservation_flags` treated *absent from the cached ortholog map* as
`no_ortholog`. Once the distance gate surfaced new genes (GNAS, CPE, ...), they
were wrongly flagged human-specific — the very overclaim we were auditing for,
reproduced in our own fix.
→ **Fix:** the flagger now resolves any missing gene **live via Ensembl** before
judging, and distinguishes three honest states: `no_ortholog` (Ensembl confirms
none), `mouse_expr_unmeasured` (ortholog exists but not in our mouse panel), and
`ortholog_unknown` (offline, never resolved). GNAS/CPE/GAD2/ERO1B correctly moved
to `mouse_expr_unmeasured`; PHGR1 stayed `no_ortholog` (verified).

## Result: permissive vs strict

| metric | beta permissive | beta strict | enterocyte permissive | enterocyte strict |
|---|---|---|---|---|
| gene-linked / 25 | 7 | 25 | 25 | 25 |
| max link distance | 46 kb | 48 kb | **879 kb** | **29 kb** |
| median link distance | 18 kb | 20 kb | 12 kb | 5 kb |
| conserved | 7 | 14 | 22 | 20 |
| flagged high-risk | 0 | 0 | 3 | 5 |

The enterocyte max link distance is the headline correction: **879 kb -> 29 kb**,
every candidate now within the range where a nearest-gene assignment is
defensible on proximity alone. Beta's max rose slightly
(46 -> 48 kb) only because the gate admitted more genuine within-50 kb genes once
the desert padding was removed; all 25 are now real gene-linked parts. The
PHGR1 `no_ortholog` differentiator survives in strict mode (ranks 9, 17, 20).
MALRD1 now correctly carries a `divergent_expression` flag (mouse ortholog
exists but its mouse enterocyte specificity is 0.07 vs 0.53 in human).

## Reproducing

```bash
python run.py "Small_Intestinal_Enterocyte"                 # original (permissive)
python run.py "Small_Intestinal_Enterocyte" --strict-links  # audited (recommended)
```

Both modes ship. The permissive default keeps the original results reproducible;
`--strict-links` is the defensible configuration and the one we recommend
reporting. Outputs: `*_catalog_strict.csv`, `permissive_vs_strict_comparison.csv`,
`reusability_funnel_strict.png`, `*_specificity_strict.png`.
