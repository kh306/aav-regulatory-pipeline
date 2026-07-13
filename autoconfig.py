"""
autoconfig.py -- LLM-proposed, deterministically-validated per-cell-type config.

This is the ONE place the pipeline touches a language model, and it is an
OFFLINE, ONE-TIME AUTHORING AID -- not a runtime dependency. The flow:

    free-text cell type  --host.llm-->  proposed {catlas_match, identity_tfs,
                                        positive_controls}  (JSON)
         |
         v  validate() against GENCODE symbols + the CATLAS column list
    validated config written to configs/generated/<key>.json
         |
         v  the deterministic pipeline (pipeline.py / motif.py / controls.py)
            reads that cached file -- NO network, NO API key, reproducible.

Why this is defensible, in one line: the LLM only PROPOSES a hypothesis (which
TFs / markers to expect); the deterministic pipeline DISPOSES of it (motif
enrichment and the positive control independently confirm or reject it using
public chromatin/expression data). No LLM output ever sets a score, a rank, a
Tau value, or a conservation call.

IMPORTANT CAVEAT (also stated in the README): the model learned this biology
from the same literature that defines the "right answer", so a passing check is
corroboration, not fully independent proof. The human-curated beta and
enterocyte configs (tf_config.py / controls.py) are therefore kept as the clean
validation anchors; generated configs are used for held-out cell types.

Generated configs are cached to disk and SHIPPED with the repo, so end users run
the pipeline fully offline. Regenerating requires host.llm (Claude Science kernel).
"""
from __future__ import annotations
import json, os, re, hashlib

HERE = os.path.dirname(__file__)
GEN_DIR = os.path.join(HERE, "configs", "generated")
LOG_DIR = os.path.join(HERE, "configs", "logs")


def _catlas_names(data_dir=None):
    import pandas as pd
    data_dir = data_dir or os.path.join(HERE, "data")
    return pd.read_csv(os.path.join(data_dir, "catlas_celltypes.csv")).iloc[:, -1].tolist()


def _gencode_symbols(data_dir=None):
    import pandas as pd
    data_dir = data_dir or os.path.join(HERE, "data")
    return set(pd.read_csv(os.path.join(data_dir, "gene_tss.tsv"), sep="\t")["gene_name"])


def build_prompt(free_text, catlas_names):
    """Prompt Claude for a structured config. The CATLAS list is supplied so the
    model must choose an existing column (name resolution) rather than invent one."""
    return (
        "You are a molecular-biology config generator for a regulatory-genomics "
        "pipeline. Given a free-text cell-type request, return ONLY a JSON object "
        "(no prose) with these keys:\n\n"
        '  "catlas_match": exactly one string from ALLOWED_CATLAS_NAMES (best match), or null.\n'
        '  "identity_tfs": 3-6 lineage-defining transcription factors (HGNC symbols, uppercase).\n'
        '  "positive_controls": 4-8 canonical cell-type-restricted marker genes (HGNC symbols).\n\n'
        "Rules: standard current HGNC symbols only; catlas_match copied verbatim "
        "from ALLOWED_CATLAS_NAMES or null; strictly valid JSON, nothing else.\n\n"
        f'REQUEST: "{free_text}"\n\n'
        f"ALLOWED_CATLAS_NAMES = {json.dumps(catlas_names)}\n"
    )


def _extract_json(text):
    m = re.search(r"```(?:json)?\s*(\{.*?\})\s*```", text, re.S) or re.search(r"(\{.*\})", text, re.S)
    if not m:
        raise ValueError("no JSON object found in model response")
    return json.loads(m.group(1))


def validate(proposed, catlas_names, gencode_symbols):
    """Deterministically filter the model's proposal. Drops any gene symbol not
    in GENCODE (catches hallucinated/aliased symbols) and rejects an out-of-list
    catlas_match. Returns (validated_config_dict, dropped_dict)."""
    cm = proposed.get("catlas_match")
    if cm not in catlas_names:
        cm = None
    tfs_ok = [g for g in proposed.get("identity_tfs", []) if g in gencode_symbols]
    tfs_bad = [g for g in proposed.get("identity_tfs", []) if g not in gencode_symbols]
    pos_ok = [g for g in proposed.get("positive_controls", []) if g in gencode_symbols]
    pos_bad = [g for g in proposed.get("positive_controls", []) if g not in gencode_symbols]
    return (
        {"catlas_match": cm, "identity_tfs": tfs_ok, "positive_controls": pos_ok},
        {"identity_tfs": tfs_bad, "positive_controls": pos_bad},
    )


def generate_config(free_text, key=None, host=None, data_dir=None, cache=True):
    """Generate + validate + cache a config for a free-text cell type.

    Requires `host` (the Claude Science kernel singleton) for the one LLM call.
    Writes configs/generated/<key>.json (validated) and configs/logs/<key>.json
    (full raw response, for audit). Returns the validated config dict.
    """
    if host is None:
        raise RuntimeError("generate_config needs the Claude Science `host` (host.llm). "
                           "Run inside a Claude Science kernel, or hand-write the config.")
    key = key or re.sub(r"[^a-z0-9]+", "_", free_text.lower()).strip("_")[:24]
    names = _catlas_names(data_dir)
    syms = _gencode_symbols(data_dir)
    prompt = build_prompt(free_text, names)
    r = host.llm(prompt, max_tokens=500)
    proposed = _extract_json(r["text"])
    validated, dropped = validate(proposed, names, syms)
    validated.update({
        "request": free_text,
        "dropped_invalid": dropped,
        "provenance": "machine-generated (host.llm) + validated against GENCODE v44 gene symbols",
        "model": r.get("model"),
        "generated": True,
    })
    if cache:
        os.makedirs(GEN_DIR, exist_ok=True)
        os.makedirs(LOG_DIR, exist_ok=True)
        json.dump(validated, open(os.path.join(GEN_DIR, f"{key}.json"), "w"), indent=2)
        json.dump({"key": key, "request": free_text, "model": r.get("model"),
                   "prompt_sha": hashlib.sha1(prompt.encode()).hexdigest()[:12],
                   "response": r["text"]},
                  open(os.path.join(LOG_DIR, f"{key}.json"), "w"), indent=2)
    return validated


def load_generated(key):
    """Load a cached generated config (offline; no LLM)."""
    p = os.path.join(GEN_DIR, f"{key}.json")
    return json.load(open(p)) if os.path.exists(p) else None
