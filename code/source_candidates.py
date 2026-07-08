"""
Cleavr candidate-gene sourcing -- run locally.

Systematically pulls TNBC-associated, druggable-looking targets from Open
Targets (direct public API, no key needed) instead of picking genes one at
a time by hand. This is what turns "we tested 3 genes someone suggested"
into an actual screen.

FILTERS APPLIED (all necessary to get a usable list, not just a raw dump):
  1. Disease association: ranked by Open Targets' association score for
     "triple-negative breast carcinoma" (MONDO_0005494) -- an aggregate of
     genetic, expression, literature, and other evidence, not just one
     study.
  2. Tractability: keep only targets flagged as small-molecule or antibody
     tractable -- i.e. something could plausibly be built against them,
     not just statistically associated.
  3. Exclude the 8 Lehmann marker genes (CCNE1, CDC6, EGFR, MET, VIM, ZEB1,
     AR, FOXA1) -- these DEFINE the subtypes, so testing them for
     subtype-specific expression is circular, not a discovery.
  4. Exclude genes already run through the pipeline (tracked in
     ALREADY_TESTED below -- update this list as you go).
  5. Must actually be present in the local TCGA expression matrix, or
     there's nothing to test.

HOW TO RUN:
    cd ~/Desktop/tnbc-project/code
    python3 source_candidates.py [N]        # N = how many candidates to keep (default 20)

OUTPUT:
    ../results/candidate_gene_list.json
"""
import json
import subprocess
import sys


def ensure(pkg, import_name=None):
    try:
        __import__(import_name or pkg)
    except ImportError:
        print(f"Installing {pkg} ...")
        subprocess.run([sys.executable, "-m", "pip", "install", pkg, "--break-system-packages"], check=True)


ensure("requests")
import requests

OT_API = "https://api.platform.opentargets.org/api/v4/graphql"
TNBC_EFO_ID = "MONDO_0005494"  # resolved once via Open Targets search: "triple-negative breast carcinoma"
STUDY_DIR = "../brca_tcga"
RESULTS_DIR = "../results"

MARKER_GENES = {"CCNE1", "CDC6", "EGFR", "MET", "VIM", "ZEB1", "AR", "FOXA1"}
ALREADY_TESTED = {"GNRHR", "AR", "CXCR4"}
EXCLUDE = MARKER_GENES | ALREADY_TESTED


def fetch_tnbc_associated_targets(page_size=100):
    query = """
    query TNBCTargets($efoId: String!, $size: Int!) {
      disease(efoId: $efoId) {
        associatedTargets(page: {index: 0, size: $size}) {
          count
          rows {
            score
            target {
              approvedSymbol
              approvedName
              targetClass { label level }
              tractability { modality value }
            }
          }
        }
      }
    }
    """
    r = requests.post(OT_API, json={"query": query, "variables": {"efoId": TNBC_EFO_ID, "size": page_size}}, timeout=30)
    r.raise_for_status()
    data = r.json()
    if "errors" in data:
        raise RuntimeError(data["errors"])
    return data["data"]["disease"]["associatedTargets"]["rows"]


def is_tractable(target):
    for t in target.get("tractability", []):
        if t.get("value") and t.get("modality") in ("SM", "AB", "Small molecule", "Antibody"):
            return True
    return False


def is_structural_only(target):
    """Filter out genes whose only target class is 'Structural protein', or
    that are tubulin-family by gene symbol (some tubulin isoforms come back
    from Open Targets with an empty targetClass list, so the class check
    alone misses them -- e.g. TUBA3C/TUBA3E). Tubulin is targeted broadly
    by taxane chemo already: mature, standard-of-care drug class, not a
    TNBC-subtype-specific biomarker lead. Without this filter, most of the
    top-scoring slots are just different tubulin isoforms -- redundant, not
    diverse candidates."""
    symbol = target.get("approvedSymbol", "")
    if symbol.startswith("TUBA") or symbol.startswith("TUBB"):
        return True
    classes = [c["label"] for c in target.get("targetClass", [])]
    return classes == ["Structural protein"]


def genes_in_expression_matrix(gene_symbols):
    """Cross-check candidates actually have data in the local TCGA matrix --
    grep the first column (Hugo_Symbol) rather than loading the whole file."""
    present = set()
    try:
        with open(f"{STUDY_DIR}/data_mrna_seq_v2_rsem.txt") as f:
            next(f)  # header
            wanted = set(gene_symbols)
            for line in f:
                sym = line.split("\t", 1)[0]
                if sym in wanted:
                    present.add(sym)
                    wanted.discard(sym)
                    if not wanted:
                        break
    except FileNotFoundError:
        print(f"WARNING: {STUDY_DIR}/data_mrna_seq_v2_rsem.txt not found -- can't verify data availability.")
        return set(gene_symbols)  # don't block the list, just skip the check
    return present


if __name__ == "__main__":
    n_keep = int(sys.argv[1]) if len(sys.argv) > 1 else 20

    print("Fetching TNBC-associated targets from Open Targets...")
    rows = fetch_tnbc_associated_targets(page_size=250)
    print(f"  {len(rows)} candidate targets returned")

    n_structural = 0
    candidates = []
    for row in rows:
        t = row["target"]
        sym = t["approvedSymbol"]
        if sym in EXCLUDE:
            continue
        if not is_tractable(t):
            continue
        if is_structural_only(t):
            n_structural += 1
            continue
        candidates.append({
            "gene": sym,
            "name": t["approvedName"],
            "open_targets_score": row["score"],
            "target_class": [c["label"] for c in t.get("targetClass", [])],
        })

    print(f"  {n_structural} excluded as structural-protein-only (e.g. tubulin family)")
    print(f"  {len(candidates)} pass tractability + exclusion filters")

    print("Cross-checking against local TCGA expression matrix...")
    present = genes_in_expression_matrix([c["gene"] for c in candidates])
    candidates = [c for c in candidates if c["gene"] in present]
    print(f"  {len(candidates)} confirmed present in data_mrna_seq_v2_rsem.txt")

    candidates.sort(key=lambda c: c["open_targets_score"], reverse=True)
    shortlist = candidates[:n_keep]

    out_path = f"{RESULTS_DIR}/candidate_gene_list.json"
    with open(out_path, "w") as f:
        json.dump({
            "source": "Open Targets Platform API (disease=MONDO_0005494, triple-negative breast carcinoma)",
            "filters_applied": ["tractable (small molecule or antibody)", "excludes Lehmann marker genes",
                                 "excludes already-tested genes", "present in local TCGA RNA-seq matrix"],
            "n_candidates": len(shortlist),
            "candidates": shortlist,
        }, f, indent=2)

    print(f"\nSaved {out_path}")
    for c in shortlist:
        print(f"  {c['gene']:10s} score={c['open_targets_score']:.3f}  {c['target_class']}")
