"""
Cleavr batch screen -- run locally, after source_candidates.py.

Runs every gene in ../results/candidate_gene_list.json through the same
pipeline GNRHR/AR/CXCR4 went through (subtype comparison, survival split,
ChEMBL druggability, PubMed literature-gap), and -- critically -- applies
Benjamini-Hochberg FDR correction across ALL the Kruskal-Wallis p-values
together.

WHY THE FDR STEP MATTERS: testing one gene, p<0.05 is a reasonable bar.
Testing 15+ genes at once, on average one of them will hit p<0.05 by pure
chance even if nothing real is going on -- that's the multiple-testing
problem flagged in docs/ROADMAP.md item 4. BH-FDR correction adjusts the
p-values so "significant after correction" actually means something at
this scale. A gene that's significant before correction but not after is
reported honestly, not dropped or hidden.

HOW TO RUN:
    cd ~/Desktop/tnbc-project/code
    python3 screen_candidates.py

OUTPUT:
    ../results/screen_results.json   -- full detail per gene
    ../results/screen_results.csv    -- flat summary table
"""
import csv
import json

import biomarker_pipeline as bp
import candidate_scoring as cs

RESULTS_DIR = "../results"


def benjamini_hochberg(pvals):
    """Standard BH step-up procedure. Returns adjusted p-values in the same
    order as the input. Implemented directly (not via a library) so the
    correction is easy to audit line by line."""
    m = len(pvals)
    indexed = sorted(range(m), key=lambda i: pvals[i])
    adjusted = [0.0] * m
    prev_min = 1.0
    for rank, idx in zip(range(m, 0, -1), reversed(indexed)):
        # rank counts down from m to 1 as we walk the sorted list in reverse
        val = pvals[idx] * m / rank
        prev_min = min(prev_min, val)
        adjusted[idx] = min(prev_min, 1.0)
    return adjusted


if __name__ == "__main__":
    with open(f"{RESULTS_DIR}/candidate_gene_list.json") as f:
        candidate_list = json.load(f)
    genes = [c["gene"] for c in candidate_list["candidates"]]
    print(f"Screening {len(genes)} candidates: {genes}\n")

    print("Building TNBC cohort (once, shared across all genes)...")
    cohort = bp.build_tnbc_cohort(extra_genes=genes)

    screened = []
    for gene in genes:
        print(f"\n--- {gene} ---")
        stats_result = bp.analyze_gene(gene, cohort)
        if stats_result is None:
            print(f"  SKIPPED: {gene} not in expression matrix")
            continue

        print("  Fetching druggability (ChEMBL) + literature gap (PubMed)...")
        drug = cs.get_druggability(gene)
        lit = cs.get_literature_gap(gene)

        screened.append({
            "gene": gene,
            "kruskal_p": stats_result["kruskal_wallis"]["p"],
            "survival_logrank_p": stats_result["survival_logrank_p"],
            "n_samples": stats_result["n_tnbc_with_expression"],
            "druggability": drug,
            "literature_gap": lit,
        })

    # FDR correction across every Kruskal-Wallis p-value from this screen
    kruskal_ps = [s["kruskal_p"] for s in screened]
    adjusted = benjamini_hochberg(kruskal_ps)
    for s, adj_p in zip(screened, adjusted):
        s["kruskal_p_fdr_adjusted"] = adj_p
        s["subtype_significant_raw"] = s["kruskal_p"] < 0.05
        s["subtype_significant_after_fdr"] = adj_p < 0.05

    # Sort: FDR-significant first, then by literature gap (under-explored first)
    gap_order = {"under-explored in TNBC": 0, "moderately studied in TNBC": 1, "well-studied in TNBC": 2}
    screened.sort(key=lambda s: (
        not s["subtype_significant_after_fdr"],
        gap_order.get(s["literature_gap"].get("category"), 3),
        s["kruskal_p_fdr_adjusted"],
    ))

    with open(f"{RESULTS_DIR}/screen_results.json", "w") as f:
        json.dump(screened, f, indent=2)

    with open(f"{RESULTS_DIR}/screen_results.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["gene", "kruskal_p", "kruskal_p_fdr_adjusted", "subtype_sig_raw",
                          "subtype_sig_after_fdr", "survival_logrank_p", "n_potent_ligands",
                          "chembl_target_found", "tnbc_paper_count", "literature_category"])
        for s in screened:
            writer.writerow([
                s["gene"], f"{s['kruskal_p']:.4g}", f"{s['kruskal_p_fdr_adjusted']:.4g}",
                s["subtype_significant_raw"], s["subtype_significant_after_fdr"],
                f"{s['survival_logrank_p']:.4g}",
                s["druggability"].get("n_potent_ligands", ""),
                s["druggability"].get("chembl_target_found", ""),
                s["literature_gap"].get("tnbc_specific_paper_count", ""),
                s["literature_gap"].get("category", ""),
            ])

    print(f"\n=== SCREEN COMPLETE: {len(screened)} genes ===")
    print(f"Significant before FDR correction: {sum(s['subtype_significant_raw'] for s in screened)}")
    print(f"Significant AFTER FDR correction:  {sum(s['subtype_significant_after_fdr'] for s in screened)}")
    print(f"\nSaved {RESULTS_DIR}/screen_results.json and {RESULTS_DIR}/screen_results.csv")
