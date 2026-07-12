# Updating your GitHub repo with the audited version

You already uploaded the first version through the browser. This is how to push
this round of changes (the self-audit: `--strict-links` mode, the four fixes,
the audit write-up, and the strict figures) the same way — no command line, no
git install.

## What changed since your first upload

**New files** (didn't exist in v1):

| File | What it is |
|---|---|
| `AUDIT_RESPONSE.md` | The self-audit: what was reviewed, the 4 fixes, permissive-vs-strict table |
| `catalogs/parts_catalog_strict.csv` | The defensible catalog (all parts gene-linked within 50 kb) |
| `catalogs/permissive_vs_strict_comparison.csv` | Before/after numbers |
| `catalogs/Pancreatic_Beta_Cell_catalog_strict.csv` | Beta strict catalog (+ `_funnel_strict.csv`, `_controls_strict.json`) |
| `catalogs/Small_Intestinal_Enterocyte_catalog_strict.csv` | Enterocyte strict catalog (+ funnel/controls) |
| `figures/reusability_funnel_strict.png` | Strict funnel (shows the 879 kb → 29 kb distance-gate fix) |
| `figures/beta_specificity_strict.png` | Strict beta specificity figure |
| `figures/enterocyte_specificity_strict.png` | Strict enterocyte specificity figure |
| `figures/summary_slide_strict.png` / `.pdf` | The corrected video slide |

**Changed files** (overwrite the v1 copies):

| File | Why it changed |
|---|---|
| `pipeline.py` | Added `link_max_dist`, `impute_unlinked`, `require_gene_link` |
| `run.py` | Added the `--strict-links` flag |
| `cross_species.py` | Fixed the false `no_ortholog` bug (now resolves genes live) |
| `README.md` | Points to `AUDIT_RESPONSE.md`, documents strict mode |

Everything else (the other `.py` modules, `LICENSE`, `requirements.txt`,
`.gitignore`, `data/`, `SUBMISSION.md`, `demo.py`) is unchanged — you don't need
to touch it.

## The easy way: re-upload the whole folder (recommended)

GitHub's web uploader overwrites files that already exist and adds new ones, so
you can just drop the entire updated folder in and it reconciles both new and
changed files in one commit.

1. Download the fresh **`aav_regulatory_pipeline.zip`** from this session and
   unzip it. This is the complete, current repo.
2. Go to your repo on github.com. Click **Add file → Upload files**.
3. Open the unzipped `aav_pipeline` folder on your computer, select **all** its
   contents, and drag them into the browser upload area. (Dragging the folder's
   *contents*, not the folder itself, keeps the same layout as your first
   upload.)
4. Wait for every file to finish uploading (the count stops climbing).
5. In the **Commit changes** box at the bottom, type a message —
   e.g. `Add self-audit + strict-links mode (distance gate, ortholog-flag fix)`.
6. Click **Commit changes**. Done — new files are added, changed files are
   overwritten, and it's one clean commit in your history.

## If you'd rather upload only what changed

Same **Add file → Upload files** screen, but drag in just the 4 changed files
plus the new ones from the tables above. For files inside `catalogs/` or
`figures/`, first navigate *into* that folder on GitHub, then Upload — otherwise
they land at the top level. The whole-folder method above avoids this fiddliness,
which is why it's recommended.

## After it's up — 30-second check

- Open `AUDIT_RESPONSE.md` on GitHub; it should render the permissive-vs-strict
  table.
- Open `README.md`; the line pointing to `AUDIT_RESPONSE.md` and the
  `--strict-links` instructions should be there.
- Click one strict figure (e.g. `figures/reusability_funnel_strict.png`) to
  confirm images display.
- Optional: the repo still runs the zero-setup demo — `python demo.py` reads the
  shipped catalogs and prints the headline results, no downloads needed.

## Reproducing the audited results (for anyone cloning the repo)

```bash
python run.py "Small_Intestinal_Enterocyte" --strict-links
python run.py "Pancreatic_Beta_Cell"        --strict-links
```

Plain `python run.py "<cell type>"` still reproduces the original (permissive)
results, so both versions are in the history and nothing is lost.
