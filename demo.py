#!/usr/bin/env python
"""
Zero-setup demo: prints the pipeline's story from the shipped catalogs.

Unlike `run.py` (which recomputes from the data backbone), this reads the
pre-built catalogs in catalogs/ and narrates the three headline results, so a
judge can see the payoff in <5 seconds without downloading anything:

    python demo.py

Shows: (1) reusability -- same command, two tissues; (2) positive/negative
controls; (3) the cross-species translatability flag, including the
human-specific part a mouse-only workflow would miss.
"""
import os, json
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CAT = os.path.join(HERE, "catalogs")


def load(tag):
    cat = pd.read_csv(os.path.join(CAT, f"{tag}_catalog.csv"))
    ctl = json.load(open(os.path.join(CAT, f"{tag}_controls.json")))
    return cat, ctl


def show(tag, nice):
    cat, ctl = load(tag)
    print(f"\n{'='*64}\n  {nice}   (python run.py \"{tag}\")\n{'='*64}")
    # top gene-linked parts
    linked = cat[cat["target_gene"].astype(str).ne("nan") & cat["target_gene"].astype(str).ne("NA")]
    top = linked.sort_values("composite_rank").head(5)
    print("  Top gene-linked candidate parts:")
    for _, r in top.iterrows():
        print(f"    rank {int(r['composite_rank']):>2}  {r['chrom']}:{int(r['start'])}"
              f"  →  {r['target_gene']:<8}  τ={r['tau_specificity']:.3f}"
              f"  [{r['translatability_flag']}]")
    pos = ctl["positive"]; neg = ctl["negative_chromatin"]
    print(f"  Positive control: found {pos['markers_found']}  → PASS={pos['pass']}")
    print(f"  Negative control: housekeeping promoter median τ="
          f"{neg['median_promoter_tau']:.3f}  → PASS={neg['pass']}")
    flags = cat["translatability_flag"].value_counts().to_dict()
    print(f"  Cross-species flags: {flags}")
    return cat


def main():
    print("\n" + "#"*64)
    print("#  Reusable, cross-species-validated AAV regulatory-element pipeline")
    print("#  ONE pipeline · ANY cell type · only the argument changes")
    print("#"*64)

    show("Pancreatic_Beta_Cell", "PANCREATIC BETA CELL")
    ent = show("Small_Intestinal_Enterocyte", "SMALL-INTESTINAL ENTEROCYTE")

    print(f"\n{'='*64}\n  THE HEADLINE: cross-species translatability\n{'='*64}")
    risk = ent[ent["translatability_flag"].isin(["no_ortholog", "divergent_expression"])]
    if len(risk):
        print("  These enterocyte parts are HIGH-RISK for mouse→human translation")
        print("  (target gene has no confident mouse ortholog / divergent expression):")
        for _, r in risk.iterrows():
            print(f"    {r['chrom']}:{int(r['start'])}  →  {r['target_gene']}"
                  f"  [{r['translatability_flag']}]")
        print("  A mouse-only workflow would have trusted these. The filter flags them.")

    proof = json.load(open(os.path.join(CAT, "reusability_proof.json")))
    print(f"\n  Reusability proof: code lines changed between the two runs above ="
          f" {proof['code_lines_changed_between_runs']}")
    print(f"  Difference between runs: {proof['difference']}\n")


if __name__ == "__main__":
    main()
