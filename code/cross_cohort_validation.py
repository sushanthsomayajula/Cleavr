"""
Cleavr cross-cohort validation -- run locally, after fetch_metabric.py.

The TCGA screen (screen_candidates.py) found 8 genes significant after FDR
correction, in a single cohort of 115 TNBC samples. That's honest but
limited: with one cohort you can't tell a real biological signal apart from
a TCGA-specific batch effect or normalization artifact. This script re-runs
the same subtype-comparison logic (Kruskal-Wallis across Lehmann subtypes,
log-rank survival split) on METABRIC -- a different cohort (2509 breast
cancer patients, 320 TNBC by ER-/PR-/HER2-), profiled on a different
platform (Illumina microarray, not RNA-seq V2 RSEM). A gene that replicates
across both a different cohort AND a different expression platform is a much
stronger candidate than one that only shows up in TCGA.

METHODOLOGY NOTE: this reuses the same marker-gene proxy for Lehmann
subtype (BL1/BL2/M/LAR) as biomarker_pipeline.py, applied to METABRIC's own
z-scored marker expression -- it does NOT reuse TCGA's subtype assignments,
since subtype has to be derived independently within each cohort's own
expression distribution.

HOW TO RUN:
    cd ~/Desktop/tnbc-project/code
    python3 fetch_metabric.py          # once, to pull ../brca_metabric/*.csv
    python3 cross_cohort_validation.py

OUTPUT:
    ../results/cross_cohort_results.json
    ../results/cross_cohort_results.csv
"""
import csv
import json
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import StandardScaler
from lifelines.statistics import logrank_test

METABRIC_DIR = "../brca_metabric"
RESULTS_DIR = "../results"

MARKER_GENES = ["CCNE1", "CDC6", "EGFR", "MET", "VIM", "ZEB1", "AR", "FOXA1"]
SUBTYPE_ORDER = ["BL1", "BL2", "M", "LAR"]
GENES_TO_VALIDATE = ["PRIM2", "POLD2", "POLD1", "RRM1", "TYMS", "RRM2", "POLE", "FTO"]


def build_metabric_tnbc_cohort():
    clin = pd.read_csv(f"{METABRIC_DIR}/clinical_patient.csv")
    tnbc_patients = clin[
        (clin["ER_STATUS"] == "Negative") &
        (clin["PR_STATUS"] == "Negative") &
        (clin["HER2_STATUS"] == "Negative")
    ][["patientId", "OS_STATUS", "OS_MONTHS"]].copy()
    print(f"METABRIC TNBC patients (ER-/PR-/HER2-): {len(tnbc_patients)}")

    expr = pd.read_csv(f"{METABRIC_DIR}/expression.csv")
    tnbc_expr = expr[expr["patientId"].isin(tnbc_patients["patientId"])].copy()
    tnbc_expr = tnbc_expr.drop_duplicates(subset="patientId")
    print(f"METABRIC TNBC samples with expression data: {len(tnbc_expr)}")

    needed_cols = MARKER_GENES + GENES_TO_VALIDATE
    missing = [g for g in needed_cols if g not in tnbc_expr.columns]
    if missing:
        print(f"WARNING - genes not found in METABRIC expression matrix: {missing}")
    tnbc_expr = tnbc_expr.dropna(subset=[g for g in MARKER_GENES if g in tnbc_expr.columns])

    # Same marker-gene proxy for Lehmann subtype as biomarker_pipeline.py,
    # computed independently within METABRIC's own expression distribution.
    scaler = StandardScaler()
    scaled = scaler.fit_transform(tnbc_expr[MARKER_GENES])
    scaled_df = pd.DataFrame(scaled, columns=MARKER_GENES, index=tnbc_expr.index)
    scaled_df["BL1_score"] = scaled_df[["CCNE1", "CDC6"]].mean(axis=1)
    scaled_df["BL2_score"] = scaled_df[["EGFR", "MET"]].mean(axis=1)
    scaled_df["M_score"] = scaled_df[["VIM", "ZEB1"]].mean(axis=1)
    scaled_df["LAR_score"] = scaled_df[["AR", "FOXA1"]].mean(axis=1)
    score_cols = ["BL1_score", "BL2_score", "M_score", "LAR_score"]
    tnbc_expr["subtype"] = scaled_df[score_cols].idxmax(axis=1).str.replace("_score", "")
    print(tnbc_expr["subtype"].value_counts())

    final_df = tnbc_expr.merge(tnbc_patients, on="patientId", how="left")
    final_df["OS_MONTHS"] = pd.to_numeric(final_df["OS_MONTHS"], errors="coerce")
    return final_df


def validate_gene(gene_symbol, cohort_df):
    if gene_symbol not in cohort_df.columns:
        return None
    df = cohort_df.copy()
    df[gene_symbol] = pd.to_numeric(df[gene_symbol], errors="coerce")
    df = df.dropna(subset=[gene_symbol])

    groups = [df.loc[df["subtype"] == s, gene_symbol].values for s in SUBTYPE_ORDER]
    groups = [g for g in groups if len(g) > 0]
    if len(groups) < 2:
        return None
    kw_stat, kw_p = stats.kruskal(*groups)

    median_val = df[gene_symbol].median()
    df["expr_group"] = np.where(df[gene_symbol] >= median_val, "high", "low")
    surv_df = df.dropna(subset=["OS_MONTHS", "OS_STATUS"]).copy()
    surv_df["event"] = surv_df["OS_STATUS"].astype(str).str.startswith("1").astype(int)

    lr_p = None
    if surv_df["expr_group"].nunique() == 2:
        lr = logrank_test(
            surv_df.loc[surv_df["expr_group"] == "high", "OS_MONTHS"],
            surv_df.loc[surv_df["expr_group"] == "low", "OS_MONTHS"],
            event_observed_A=surv_df.loc[surv_df["expr_group"] == "high", "event"],
            event_observed_B=surv_df.loc[surv_df["expr_group"] == "low", "event"],
        )
        lr_p = float(lr.p_value)

    return {
        "gene": gene_symbol,
        "metabric_n": int(len(df)),
        "metabric_kruskal_p": float(kw_p),
        "metabric_subtype_significant": bool(kw_p < 0.05),
        "metabric_survival_logrank_p": lr_p,
        "metabric_subtype_medians": {s: float(np.median(g)) for s, g in zip(SUBTYPE_ORDER, groups)} if len(groups) == 4 else None,
    }


if __name__ == "__main__":
    with open(f"{RESULTS_DIR}/screen_results.json") as f:
        tcga_results = {r["gene"]: r for r in json.load(f)}

    cohort = build_metabric_tnbc_cohort()

    combined = []
    for gene in GENES_TO_VALIDATE:
        print(f"\n--- {gene} ---")
        m = validate_gene(gene, cohort)
        if m is None:
            print(f"  SKIPPED: {gene} not usable in METABRIC cohort")
            continue
        t = tcga_results.get(gene, {})
        replicated = bool(m["metabric_subtype_significant"] and t.get("subtype_significant_after_fdr"))
        print(f"  TCGA:     p_fdr={t.get('kruskal_p_fdr_adjusted', 'NA')}, sig={t.get('subtype_significant_after_fdr')}")
        print(f"  METABRIC: p={m['metabric_kruskal_p']:.4g}, sig={m['metabric_subtype_significant']}")
        print(f"  Replicated in both cohorts: {replicated}")
        combined.append({
            "gene": gene,
            "tcga_kruskal_p_fdr_adjusted": t.get("kruskal_p_fdr_adjusted"),
            "tcga_subtype_significant_after_fdr": t.get("subtype_significant_after_fdr"),
            "tcga_n": t.get("n_samples"),
            "metabric_kruskal_p": m["metabric_kruskal_p"],
            "metabric_subtype_significant": m["metabric_subtype_significant"],
            "metabric_n": m["metabric_n"],
            "metabric_survival_logrank_p": m["metabric_survival_logrank_p"],
            "replicated_in_both_cohorts": replicated,
        })

    combined.sort(key=lambda r: (not r["replicated_in_both_cohorts"], r["metabric_kruskal_p"]))

    with open(f"{RESULTS_DIR}/cross_cohort_results.json", "w") as f:
        json.dump(combined, f, indent=2)

    with open(f"{RESULTS_DIR}/cross_cohort_results.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["gene", "tcga_p_fdr", "tcga_sig", "tcga_n",
                          "metabric_p", "metabric_sig", "metabric_n",
                          "metabric_survival_p", "replicated_in_both"])
        for r in combined:
            writer.writerow([
                r["gene"], r["tcga_kruskal_p_fdr_adjusted"], r["tcga_subtype_significant_after_fdr"],
                r["tcga_n"], f"{r['metabric_kruskal_p']:.4g}", r["metabric_subtype_significant"],
                r["metabric_n"], r["metabric_survival_logrank_p"], r["replicated_in_both_cohorts"],
            ])

    n_replicated = sum(r["replicated_in_both_cohorts"] for r in combined)
    print(f"\n=== CROSS-COHORT VALIDATION COMPLETE: {len(combined)} genes checked ===")
    print(f"Replicated (significant in BOTH TCGA-after-FDR and METABRIC): {n_replicated}")
    print(f"\nSaved {RESULTS_DIR}/cross_cohort_results.json and {RESULTS_DIR}/cross_cohort_results.csv")
