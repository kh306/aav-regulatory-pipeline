# Can we computationally find — and pre-screen for cross-species translatability — tissue-specific AAV "on-switches" for the tissues where AAV actually matters?

**Research track submission.** This is a research project answering a scientific
question. **OrthoGate** — the pipeline in this repository — is the *method* we
built to answer it; the answer is in the results. (The name: it screens
regulatory elements at the *ortholog* gate — does the target biology carry across
species — before any wet-lab commitment.)

### The research question

> **Can tissue-specific candidate AAV regulatory elements for *extrahepatic* cell
> types be identified genome-wide from single-cell chromatin accessibility alone —
> and can their mouse↔human translatability be predicted computationally, before
> any wet-lab work?**

Two sub-questions follow from it: (1) does a single, tissue-agnostic method
recover the *known* identity biology of each cell type it's pointed at (a
validation that it's finding real signal, not noise)? and (2) is cross-species
translatability a filter that actually changes which candidates you'd pursue?

### What we found (short version)

1. **Yes to identification — in 5 of the 7 cell types tested.** One specificity →
   size → motif funnel, run unchanged, recovered the canonical identity locus as a
   top hit in **five of seven** cell types — INS (beta), GCG (alpha), RYR2 (heart),
   ALB (liver), OLFM4/PIGR (intestine) — with housekeeping loci correctly rejected
   and motif enrichment of 4.7× (p = 1.6×10⁻¹³). In the other two (pancreatic
   ductal, fetal metanephric/kidney) the positive control **failed**: the funnel
   still returned specific candidates (SERPINA3/TACSTD2/TM4SF1 for ductal, GPX3 for
   kidney), but not the curated marker genes. We report these failures rather than
   hide them — they show the positive control is a genuine test that can fail, and
   they flag cell types where either the markers or the gene-linking need work.
2. **Yes, translatability matters.** A meaningful share of the top human
   candidates serve genes whose biology is *not* conserved in mouse — most
   starkly **PHGR1**, an intestinal target with **no mouse ortholog at all**. A
   mouse-first program would chase these into a dead end. Predicting this
   computationally, up front, is the headline result.

### Background

Gene therapy works by delivering a healthy gene into a patient's cells, packed
inside a harmless virus called **AAV**. For safety, you usually want that gene
to switch on in **only one type of cell** — not everywhere in the body. The DNA
"on-switch" that controls this is called a **regulatory element** (a promoter or
enhancer). Finding a *cell-type-specific* on-switch, and knowing in advance
whether it will behave the same in a mouse experiment as in a human, is the
problem this project studies. The method we built takes a cell type and returns a
short, ranked, cross-species-annotated list of candidate on-switches — and,
crucially, **the same method works for any cell type without being rewritten**,
which is what let us ask the question across seven cell types instead of one.

### Why AAV, and why *these* tissues (the LNP contrast)

There are two main ways to deliver a gene today: **lipid nanoparticles (LNP)** and
**AAV**. LNP has largely won the **liver** — lipid nanoparticles naturally drain
to and are taken up by liver cells (hepatocytes), which is why the mRNA vaccines
and the first approved RNA-interference drugs are LNP-based and liver-centric. If
your target is the liver, you increasingly don't need AAV at all.

**AAV's real, durable advantage is everything *outside* the liver** — heart,
brain and nervous system, muscle, retina, pancreas. These "extrahepatic" tissues
are exactly the ones LNP struggles to reach, and they're where tissue-specific
AAV regulatory elements matter most. **This tool is therefore aimed squarely at
extrahepatic cell types.** The flagship examples below — pancreatic islet cells
(diabetes) and heart muscle — are LNP-hard-to-reach tissues where a
cell-type-specific AAV switch is a genuine advance, not a reinvention of what LNP
already does. (Liver/hepatocyte is included only as a technical test case, not a
therapeutic pitch.)

### Why anyone should care

- **Safety and dose.** An on-switch that fires in only the target cell means less
  off-target effect, a lower dose, and a smaller immune reaction. (Example: a
  beta-cell-only switch for a diabetes therapy that spares the exocrine pancreas.)
- **It targets where AAV actually wins.** By focusing on extrahepatic cell types
  (islet, heart, and beyond), the tool addresses the delivery problems LNP can't
  — rather than competing with LNP on its home turf, the liver.
- **It saves the most expensive kind of failure.** Most gene-therapy testing
  starts in mice. A switch that works in a mouse often *fails* in humans or
  monkeys — and you usually find out only after months of animal work. Our tool
  **predicts which switches won't carry over to humans, before any lab work
  starts.** (That also means fewer wasted animals.)
- **It's a reusable "parts list."** Instead of designing every switch from
  scratch, you build a catalog of characterized, reusable parts — the same way
  engineers reuse standard components.

---

## How the tool works — the "funnel"

You run it with a single command, and the *only* thing that changes between cell
types is what you name.

### What do I type for the cell type? (three ways)

**You do NOT need to know the atlas's internal names.** The atlas (CATLAS) uses
technical labels — beta cells are `Pancreatic_Beta_Cell`, kidney cells are
`Fetal_Metanephric_Cell` — that you can't be expected to guess. There are three
ways to tell the tool which cell type you want, easiest first:

**1. Plain English — `--auto` (recommended).** Describe the cell type however you
naturally would. Claude reads your phrase and picks the correct name *from the
fixed list of 222 atlas cell types* (it cannot invent one), then the pipeline
runs on it:

```
python run.py --auto "insulin-producing cells"
python run.py --auto "glucagon-producing pancreatic cells"
python run.py --auto "heart muscle"
```

If nothing genuinely matches your phrase, it says so rather than guessing wrong —
just rephrase more specifically. (This is the ONE place the tool calls Claude,
and only to translate your words into an atlas name; see "Which steps use AI"
below.)

**2. Browse the menu — `--list`.** Prints the valid names so you can pick one. A
`*` marks cell types that have full ranking + cross-species ready:

```
python run.py --list            # all 222 names
python run.py --list pancrea    # only names containing "pancrea"
```

**3. Type the name yourself.** The original way still works, and a unique
substring is enough (you don't need the full string):

```
python run.py "Hepatocyte"
python run.py "beta"            # substring — finds Pancreatic_Beta_Cell
```

### Recommended flag: `--strict-links`

Add `--strict-links` for the more conservative, defensible results (this is the
setting behind the numbers reported below and in the `*_strict` output files):

```
python run.py --auto "glucagon-producing pancreatic cells" --strict-links
python run.py "Pancreatic_Beta_Cell" --strict-links
```

> **Honest note on names.** `--auto` resolves the *atlas* name (the step that was
> actually hard to guess). A second naming system — the expression database's
> label, used for the ranking and cross-species steps — is currently handled by a
> built-in lookup covering the cell types that have prebuilt expression panels
> (beta, enterocyte, pancreatic alpha/ductal, hepatocyte, cardiomyocyte, kidney).
> For a brand-new cell type outside that set, the atlas name still resolves and
> the funnel still runs, but ranking falls back to chromatin-specificity only
> until an expression panel is built for it. This is a documented next step, not
> a hidden failure.

Starting from ~1.1 million candidate DNA regions, it narrows down through a
series of filters — like a funnel — keeping only what could actually work:

| Step | What it keeps | Plain-language reason |
|---|---|---|
| **Start** | ~1.1 million DNA regions across 222 human cell types | A free public atlas (CATLAS) that records, for each region, which cell types it's active in. |
| **1. Active** | regions switched on in your target cell type | An off-switch is useless. |
| **2. Specific** | ...and switched *off* in most of the other 221 cell types | We want one-cell-type switches, not everywhere-switches. |
| **3. Assemble** | neighbouring active regions merged into candidate switches | Real switches are stretches of DNA, not single points. |
| **4. Size** | only candidates small enough to fit in AAV | AAV is a tiny container (~4,700 DNA letters total); the switch must fit alongside the therapeutic gene, so we cap it (~800 letters). |
| **5. Rank** | the top 25, scored and ordered | Best candidates first. |
| **6. Translatability** | each one flagged: will it carry over to mouse/human? | The headline feature (below). |

**Why "reusable" is a real claim, not marketing.** Under the hood there is *no*
cell-type-specific code — no "if beta cell do this, if intestine do that." The
cell-type name is just a label that picks one column out of a shared data table.
Everything a given cell type "knows" lives in editable data files, not in the
program logic. We checked this literally: the two runs above differ by **zero
lines of code** — only the text in quotes (documented in
`catalogs/reusability_proof.json`).

---

## What we found (results)

**Summary across all seven cell types.** The same pipeline, run unchanged (only
the cell-type argument differs), produced a ranked shortlist for each. The
positive control asks whether the tissue's known identity gene shows up as the
**#1-ranked candidate**:

| Cell type | Tissue | AAV vs LNP | Positive control | #1-ranked candidate gene | Config source |
|---|---|---|---|---|---|
| Beta | pancreatic islet | extrahepatic | ✅ pass | **INS** (insulin) | hand-curated (anchor) |
| Alpha | pancreatic islet | extrahepatic | ✅ pass | **GCG** (glucagon) | **AI-generated** (held-out demo) |
| Cardiomyocyte | heart | extrahepatic | ✅ pass | **RYR2** (cardiac Ca²⁺ channel) | hand-curated |
| Enterocyte | intestine | extrahepatic | ✅ pass | **OLFM4** (+ APOA4, FABP6, PIGR) | hand-curated (anchor) |
| Hepatocyte | liver | *LNP territory* | ✅ pass | **ALB** (+ SERPINA1) | hand-curated |
| Pancreatic ductal | pancreas | extrahepatic | ❌ fail | SERPINA3 / TACSTD2 / TM4SF1 | hand-curated |
| Metanephric (kidney) | kidney | extrahepatic | ❌ fail | GPX3 | hand-curated |

**Read this honestly:** in **5 of the 7** cell types the pipeline recovered the
canonical identity gene as its top hit. In the two "fail" rows the funnel still
returned specific candidate genes — they simply weren't the curated markers, so
the positive control correctly reports a miss rather than a pass. We keep these
failures visible: they show the control is a real test that *can* fail, and they
mark cell types where the markers or the gene-linking need more work. One row
(alpha) used a config **Claude proposed with no human curation**, and the control
still passed — a held-out demonstration that the deterministic pipeline validates
an AI-proposed hypothesis rather than trusting it.

**The funnel behaves the same for every cell type.** ~1.1 million regions narrow
down to a ranked shortlist using identical settings. See
`figures/reusability_funnel_strict.png`.

**It finds the *right* biology (positive control).** A good test of any method
is whether it rediscovers what we already know:
- **Beta cells:** the #1 candidate sits right at the **insulin (INS)** gene,
  followed by **IAPP** — the two hormones beta cells are famous for. Exactly
  right.
- **Intestine:** the validated positive-control markers recovered are
  **OLFM4, APOA4, FABP6, PIGR** — all classic intestinal markers. (Other real
  intestinal genes such as RBP2 and ANPEP also appear high in the list.)

**It rejects the *wrong* things (negative control).** "Housekeeping" genes
(like ACTB, GAPDH) are switched on in *every* cell type. The tool correctly
throws their switches out — they score as almost completely non-specific and
never reach the shortlist.

**It recognizes the master regulators (motif check).** Beta cells are controlled
by specific proteins (PDX1, NKX6-1, and others) that leave a recognizable
"footprint" on DNA. The beta shortlist is **4.7× richer** in those footprints
than random regions would be — a very strong statistical signal (p = 1.6×10⁻¹³).
The tool is picking up the known control machinery on its own.

**The headline: it predicts cross-species translatability.** For every candidate,
the tool asks *"is the gene this switch controls even present in mice?"* If it
isn't, no amount of mouse testing can ever validate that switch. The standout
example: in the intestine list, one top candidate is driven by **PHGR1 — a gene
that has no mouse counterpart at all** (we confirmed this directly against a gene
database). A mouse-only study would chase it and hit a dead end. Our tool flags
it up front, for free. Every candidate is labelled **conserved** (safe to test in
mouse), **flagged high-risk** (won't carry over), or **unknown**.

---

## Which steps use AI, and which are deterministic

The scientific results are produced by deterministic code. AI (Claude) is used in
exactly one place — translating your plain-English cell-type request into the
atlas's technical name and proposing which genes to *expect* — and never touches
a score, a rank, or a conclusion.

| Step | AI or deterministic? | Notes |
|---|---|---|
| Resolve "glucagon-producing pancreatic cells" → `Pancreatic_Alpha_Cell_1` (`--auto`) | **AI** | Constrained to the 222 real atlas names; can't invent one. One-time, cached. |
| Propose which TFs / marker genes to expect | **AI** | A *hypothesis*, then checked by the deterministic steps below. Validated against real gene symbols; used only for held-out cell types, never for the beta/enterocyte anchors. |
| Accessibility → specificity (Tau) → candidate blocks → size cap | deterministic | Pure data. No AI. |
| Ranking, motif enrichment, ortholog calls, cross-species flags | deterministic | Pure data. No AI. |
| Positive / negative controls | deterministic | Tests whether the data recovered the expected genes. |

**Why this is safe:** the AI only *proposes* (a name, a hypothesis); the
deterministic pipeline *disposes* (finds the switches, scores them, validates
against known biology). If the AI proposed a wrong marker, the positive control
simply wouldn't recover it. The shipped configs are cached to disk, so the
pipeline runs **fully offline** — no API key or internet needed to reproduce the
results. See `autoconfig.py` for the generate-then-validate implementation and
`configs/logs/` for the auditable record of every AI response.

**One honest caveat:** Claude learned this biology from the same literature that
defines the "right answer", so when a control confirms an AI-proposed marker,
that is corroboration, not fully independent proof. This is why the beta and
enterocyte configs are kept 100% human-curated as clean validation anchors, and
AI-proposed configs are used only for additional (held-out) cell types.

---

## What this tool does *not* do (honest limitations)

This is a tool for **narrowing down and prioritizing candidates** — a smart
short-list you take *into* the lab, not a guaranteed working promoter. It was
also stress-tested against its own blind spots (see **`AUDIT_RESPONSE.md`** for
the self-audit and the more conservative `--strict-links` mode, which is the one
we recommend). The honest limits, in plain terms:

1. **"Switched on" is not the same as "actively driving a gene."** The atlas tells
   us a DNA region is *accessible* (open for business), but some open regions are
   brakes or bystanders, not accelerators. A candidate could be open yet not
   actually turn a gene on. Only a lab test confirms it.
2. **These are switch *parts*, not finished switches.** A candidate enhancer still
   has to be paired with a basic promoter in a real construct, and that
   combination doesn't always work as expected. Assembly and testing are still
   needed.
3. **"Specific" is measured against 222 cell types — not every cell in the body.**
   A candidate that looks target-specific here could still be active in some cell
   type the atlas didn't include. Some off-target risk always remains.
4. **We link each switch to its likely gene by nearest neighbour, not by
   measured 3D contact.** Ideally you'd use experimental data showing which
   switch physically touches which gene; those data weren't publicly available
   when we built this, so we use "closest gene" as a stand-in and label it
   honestly in every result. Note this is a proximity guess: real switches can
   act over very long distances (the famous ZRS switch controls its gene from
   ~1 million letters away), which "closest gene" can't confirm. The tool ships
   two settings: the plain run allows links up to ~1 million letters (permissive,
   keeps long-range candidates but includes weak guesses), and the recommended
   `--strict-links` run tightens this to a conservative 50,000-letter limit where
   "closest gene" is a defensible call. The reported/"strict" results use the
   50 kb setting. The design also lets a better, contact-based method drop in
   later.
5. **The cross-species check is at the *gene* level, not the DNA-sequence level.**
   We ask "is the gene this switch serves also present and active in mice?" — not
   "does this exact DNA sequence exist in mice?" Switch *sequences* evolve fast
   and often differ between species even when the gene is identical, so the
   gene-level question is the more meaningful one for translatability. Every
   result says which check was used.
6. **The AI config check verifies that a gene *exists*, not that it's the *right*
   gene.** When you use `--auto`, Claude proposes which transcription factors and
   marker genes to expect, and we filter that list against the real gene database
   — so invented/misspelled symbols are dropped. But a proposed symbol that is a
   *real gene which simply isn't a genuine marker for that cell type* passes the
   check (e.g. a valid symbol that belongs to an unrelated gene). The downstream
   positive control is the backstop — a wrong marker just won't be recovered, so
   it can't inflate the result — but it means the machine-generated expected-gene
   list should be eyeballed, not trusted blindly. This is exactly why the
   beta/enterocyte anchors stay human-curated and only held-out cell types use
   the generated config.

**Planned improvements:** semantic (not just existence) validation of
AI-proposed markers — e.g. cross-checking each against an expression database
before accepting it; contact-based gene linking; comparing the actual DNA
sequences across species; extending the translatability check to monkeys (the
data already include macaque, marmoset, and chimp); and adding an AI model that
predicts whether a switch is actively driving a gene (addressing limit #1).

---

## Where to look (if you don't want to run any code)

You can understand everything this project found without touching the command
line — just open these files:

- **`figures/summary_slide_strict.png`** — the one-page visual summary. Start here.
- **`figures/reusability_funnel_strict.png`** — the funnel and the results, as a
  picture.
- **`catalogs/parts_catalog_strict.csv`** — the actual answer: the ranked list of
  candidate switches for both cell types. Opens in Excel or Google Sheets. Each
  row is one candidate; useful columns are `target_gene` (the gene it likely
  controls), `translatability_flag` (conserved / high-risk / unknown), and
  `tau_specificity` (how cell-type-specific it is, 0–1).
- **`AUDIT_RESPONSE.md`** — how we stress-tested our own results and what we fixed.

*Everything below is for someone who wants to re-run the analysis on a computer.
You don't need it to read the results.*

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
