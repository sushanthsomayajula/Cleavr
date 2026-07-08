"""
Cleavr candidate-gene scoring -- run locally, after biomarker_pipeline.py.

Takes a gene that's already been run through biomarker_pipeline.py (so
../results/<gene>_analysis_results.json exists) and adds two more axes of
evidence, pulled live from public APIs (no API key needed):

  1. STATISTICAL SIGNAL   -- already computed: does expression differ by
                              TNBC subtype, is it survival-associated?
                              (read straight from the existing results JSON)

  2. DRUGGABILITY          -- from ChEMBL: does a known ligand/drug already
                              exist for this target, and how potent is it?
  3. LITERATURE GAP        -- from PubMed: how many TNBC-specific papers
                              already exist on this gene? Low count + real
                              statistical signal = an under-explored lead,
                              which is the actual thing this tool is for.

This deliberately does NOT collapse the three axes into one composite
"score" -- a single number would imply a precision the underlying evidence
doesn't have (see docs/ROADMAP.md on p-hacking discipline). Instead it
prints a scorecard and leaves the judgment call to a human.

HOW TO RUN (after running biomarker_pipeline.py for the same gene(s) first):
    cd ~/Desktop/tnbc-project/code
    python3 candidate_scoring.py GNRHR AR CXCR4

OUTPUT:
    ../results/candidate_scorecard.json  -- one entry per gene
"""
import json
import subprocess
import sys
import time
import urllib.parse


def ensure(pkg, import_name=None):
    try:
        __import__(import_name or pkg)
    except ImportError:
        print(f"Installing {pkg} ...")
        subprocess.run([sys.executable, "-m", "pip", "install", pkg, "--break-system-packages"], check=True)


ensure("requests")
import requests

RESULTS_DIR = "../results"
CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
PUBMED_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

HEADERS = {"User-Agent": "Cleavr-TNBC-biomarker-tool (research use, contact: sushanthsomayajula@gmail.com)"}


def get_stats_signal(gene_symbol):
    """Read the stats this gene already has from biomarker_pipeline.py's output."""
    path = f"{RESULTS_DIR}/{gene_symbol.lower()}_analysis_results.json"
    try:
        with open(path) as f:
            r = json.load(f)
    except FileNotFoundError:
        print(f"  No {path} found -- run biomarker_pipeline.py for {gene_symbol} first.")
        return None
    return {
        "subtype_kruskal_p": r["kruskal_wallis"]["p"],
        "subtype_significant": r["kruskal_wallis"]["p"] < 0.05,
        "survival_logrank_p": r["survival_logrank_p"],
        "survival_significant": r["survival_logrank_p"] < 0.05,
        "n_samples": r["n_tnbc_with_expression"],
    }


def get_druggability(gene_symbol):
    """Query ChEMBL directly: is there a known target entry, known potent
    ligands, and/or an approved drug for this gene?"""
    try:
        # target_synonym__iexact, not icontains -- a substring match on a short
        # symbol like "AR" matches unrelated targets (e.g. "carbonic anhydrase"
        # contains "ar"). Caught this by spot-checking AR's result before
        # trusting it: icontains returned Carbonic Anhydrase 2, not the
        # androgen receptor. Exact match on the gene symbol only.
        r = requests.get(
            f"{CHEMBL_BASE}/target.json",
            params={"target_synonym__iexact": gene_symbol, "organism": "Homo sapiens",
                    "target_type": "SINGLE PROTEIN", "limit": 5},
            headers=HEADERS, timeout=20,
        )
        r.raise_for_status()
        targets = r.json().get("targets", [])
    except requests.RequestException as e:
        return {"error": str(e)}

    if not targets:
        return {"chembl_target_found": False, "n_potent_ligands": 0, "has_approved_drug": False}

    target_id = targets[0]["target_chembl_id"]
    target_name = targets[0].get("pref_name")

    # Potent ligands: pChEMBL >= 6 (roughly 1 uM or better)
    try:
        r = requests.get(
            f"{CHEMBL_BASE}/activity.json",
            params={"target_chembl_id": target_id, "pchembl_value__gte": 6, "limit": 1},
            headers=HEADERS, timeout=20,
        )
        r.raise_for_status()
        n_potent = r.json().get("page_meta", {}).get("total_count", 0)
    except requests.RequestException:
        n_potent = None

    # Approved-drug mechanism check
    try:
        r = requests.get(
            f"{CHEMBL_BASE}/mechanism.json",
            params={"target_chembl_id": target_id, "limit": 10},
            headers=HEADERS, timeout=20,
        )
        r.raise_for_status()
        mechanisms = r.json().get("mechanisms", [])
    except requests.RequestException:
        mechanisms = []

    drug_names = [m.get("molecule_chembl_id") for m in mechanisms]
    action_types = sorted(set(m.get("action_type") for m in mechanisms if m.get("action_type")))

    return {
        "chembl_target_found": True,
        "chembl_target_id": target_id,
        "chembl_target_name": target_name,
        "n_potent_ligands": n_potent,
        "n_known_mechanisms": len(mechanisms),
        "known_action_types": action_types,
        "example_molecule_chembl_ids": drug_names[:5],
    }


def get_literature_gap(gene_symbol):
    """Query PubMed E-utilities directly: how many papers connect this gene
    to TNBC specifically? Low count = under-explored in this exact context,
    even if the gene is well-studied generally."""
    query = f'{gene_symbol}[Title/Abstract] AND "triple negative breast"[Title/Abstract]'
    try:
        r = requests.get(
            f"{PUBMED_BASE}/esearch.fcgi",
            params={"db": "pubmed", "term": query, "retmode": "json"},
            headers=HEADERS, timeout=20,
        )
        r.raise_for_status()
        count = int(r.json()["esearchresult"]["count"])
    except (requests.RequestException, KeyError, ValueError) as e:
        return {"error": str(e)}

    if count <= 5:
        category = "under-explored in TNBC"
    elif count <= 50:
        category = "moderately studied in TNBC"
    else:
        category = "well-studied in TNBC"

    time.sleep(0.4)  # be polite to NCBI's unauthenticated rate limit (3 req/s)
    return {"tnbc_specific_paper_count": count, "category": category, "query_used": query}


def score_candidate(gene_symbol):
    print(f"\n=== {gene_symbol} ===")
    stats = get_stats_signal(gene_symbol)
    if stats is None:
        return None

    print("  Fetching druggability from ChEMBL...")
    drug = get_druggability(gene_symbol)
    print("  Fetching literature-gap signal from PubMed...")
    lit = get_literature_gap(gene_symbol)

    scorecard = {
        "gene": gene_symbol,
        "statistical_signal": stats,
        "druggability": drug,
        "literature_gap": lit,
    }

    print(f"  Subtype-significant: {stats['subtype_significant']} (p={stats['subtype_kruskal_p']:.4f})")
    print(f"  Survival-significant: {stats['survival_significant']} (p={stats['survival_logrank_p']:.4f})")
    if drug.get("chembl_target_found"):
        print(f"  ChEMBL target: {drug['chembl_target_name']} ({drug['chembl_target_id']}), "
              f"{drug['n_potent_ligands']} potent ligands, mechanisms: {drug['known_action_types']}")
    else:
        print("  No ChEMBL target entry found.")
    if "tnbc_specific_paper_count" in lit:
        print(f"  TNBC-specific literature: {lit['tnbc_specific_paper_count']} papers ({lit['category']})")

    return scorecard


if __name__ == "__main__":
    genes = sys.argv[1:]
    if not genes:
        print("Usage: python3 candidate_scoring.py GENE1 [GENE2 ...]")
        sys.exit(1)

    scorecards = []
    for gene in genes:
        sc = score_candidate(gene)
        if sc:
            scorecards.append(sc)

    out_path = f"{RESULTS_DIR}/candidate_scorecard.json"
    # Merge with any existing scorecard entries instead of clobbering them
    existing = []
    try:
        with open(out_path) as f:
            existing = json.load(f)
    except FileNotFoundError:
        pass
    existing_genes = {e["gene"] for e in existing}
    merged = [e for e in existing if e["gene"] not in {s["gene"] for s in scorecards}] + scorecards

    with open(out_path, "w") as f:
        json.dump(merged, f, indent=2)
    print(f"\nSaved {out_path} ({len(merged)} genes total)")
