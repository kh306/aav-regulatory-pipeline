# Putting this on GitHub (first-timer walkthrough)

You don't need the command line or `git` for this. The browser-upload path
below is the simplest way to get a working repo, and it's fine for the
hackathon.

## Path A — browser upload (recommended, no CLI)

1. **Make a GitHub account** at github.com if you don't have one (free).
2. **Create the repo:** click the **+** (top-right) → **New repository**.
   - Repository name: `aav-regulatory-pipeline`
   - Description: *Reusable, cross-species-validated pipeline for tissue-specific AAV regulatory-element design*
   - **Public** (required so the judges can see it).
   - Do **not** tick "Add a README" (we already have one).
   - Click **Create repository**.
3. **Unzip** `aav_regulatory_pipeline.zip` on your computer. You'll get a folder
   `aav_pipeline/`. Open it so you see `README.md`, `run.py`, `data/`, etc.
4. On the empty repo page, click **uploading an existing file** (the link in
   "Get started by …"). 
5. **Drag the *contents* of `aav_pipeline/`** (not the folder itself) into the
   browser — select all files and subfolders inside `aav_pipeline/` and drop
   them in. GitHub preserves the subfolders (`data/`, `catalogs/`, `figures/`).
   - The upload is ~29 MB; it may take a minute.
   - GitHub's browser upload caps individual files at 25 MB. All files here are
     under that — the largest is `data/catlas_ccre_coords.parquet` at ~19 MB,
     followed by `data/catlas_accessibility_matrix.npz` at ~13 MB — so you're
     fine. If any single file is ever rejected for size, just delete it before
     uploading — `build_backbone.py` regenerates the `data/` files.
6. Scroll down, add a commit message ("Initial commit"), click **Commit changes**.
7. Done — your README renders on the repo home page automatically.

## Path B — git command line (optional)

If you'd rather use the CLI (needs a Personal Access Token, *not* your
password — create one at Settings → Developer settings → Personal access
tokens):

```bash
cd aav_pipeline
git init
git add .
git commit -m "Initial commit"
git branch -M main
git remote add origin https://github.com/<your-username>/aav-regulatory-pipeline.git
git push -u origin main      # paste the token when prompted for a password
```

## After it's up — quick checks

- The **README** shows on the home page with the tables and figure links.
- Click into `figures/` — the PNGs preview inline.
- Open `catalogs/parts_catalog.csv` — GitHub renders CSVs as a sortable table;
  sort by `translatability_flag` to show the conserved / high-risk split live.

## What's in the repo vs. rebuildable

The zip ships the **compact cached backbone** (`data/*.parquet`, `*.npz`, motif
PWMs, expression panels) so `python run.py` and `python demo.py` work
immediately after cloning. The **800 MB of raw CATLAS BED files are *not*
included** — they're excluded by `.gitignore` and regenerated with
`python build_backbone.py` if anyone wants to rebuild from scratch.

## Tip for the demo

For a judge who just cloned the repo, `python demo.py` prints the whole story
(reusability + controls + the cross-species headline) in a few seconds without
recomputing anything.
