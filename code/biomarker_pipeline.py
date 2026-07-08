"""
Cleavr TNBC biomarker pipeline -- run locally.

Formerly run_gnrhr_analysis.py, hardcoded to one gene (GNRHR). This is the
refactored version: cohort construction and subtype assignment run once,
and any gene symbol can be tested against that cohort via analyze_gene().

HOW TO RUN:
1. Open Terminal
2. cd into the folder with brca_tcga/ inside, e.g.:
     cd ~/Desktop/tnbc-project/code
3. Install anything missing (safe to run even if some are already installed):
     pip3 install pandas numpy scipy matplotlib seaborn scikit-learn lifelines --break-system-packages
4. Run it, e.g.:
     python3 biomarker_pipeline.py GNRHR AR
   (defaults to GNRHR alone if you don't pass any gene symbols)

OUTPUT (written into ../results/<gene>/, one folder per gene):
  - <gene>_by_subtype.png       <- expression by TNBC subtype
  - <gene>_survival_km.png      <- survival curve, high vs low expression
  - <gene>_analysis_results.json <- all the stats numbers
  - <gene>_tnbc_final.csv       <- full per-patient table

Note on filenames: the original GNRHR-only script wrote analysis_results.json
and gnrhr_tnbc_final.csv (no "gnrhr_" prefix on the JSON). Now that this
runs multiple genes, every output is prefixed with the gene name so runs
don't overwrite each other.
"""
import json
import os
import subprocess
import sys
from itertools import combinations


def ensure(pkg, import_name=None):
    try:
        __import__(import_name or pkg)
    except ImportError:
        print(f"Installing {pkg} ...")
        subprocess.run([sys.executable, "-m", "pip", "install", pkg, "--break-system-packages"], check=True)


for pkg, imp in [("pandas", None), ("numpy", None), ("scipy", None),
                  ("matplotlib", None), ("seaborn", None),
                  ("scikit-learn", "sklearn"), ("lifelines", None)]:
    ensure(pkg, imp)

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import StandardScaler
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

STUDY_DIR = "../brca_tcga"
RESULTS_DIR = "../results"

# Fixed set used only to assign Lehmann subtype (BL1/BL2/M/LAR). This does
# NOT change based on which gene you're testing -- it's the coordinate
# system the target gene gets measured against.
MARKER_GENES = ["CCNE1", "CDC6", "EGFR", "MET", "VIM", "ZEB1", "AR", "FOXA1"]
SUBTYPE_ORDER = ["BL1", "BL2", "M", "LAR"]


def build_tnbc_cohort(extra_genes=None):
    """Build the TNBC cohort (ER-/PR-/HER2-) with Lehmann subtype assigned.
    Gene-independent -- runs once no matter how many genes you go on to test.
    extra_genes pulls additional gene columns into the same expression-matrix
    read so build_tnbc_cohort + analyze_gene only touch the big file once."""
    genes_needed = set(MARKER_GENES) | set(extra_genes or [])

    # 1. Clinical data -> TNBC cohort
    clin_raw = pd.read_csv(f"{STUDY_DIR}/data_clinical_patient.txt", sep="\t", header=None, skiprows=4)
    header = clin_raw.iloc[0].tolist()
    clin = clin_raw.iloc[1:].copy()
    clin.columns = header
    clin = clin.reset_index(drop=True)

    tnbc_patients = clin[
        (clin["ER_STATUS_BY_IHC"] == "Negative") &
        (clin["PR_STATUS_BY_IHC"] == "Negative") &
        (clin["IHC_HER2"] == "Negative")
    ][["PATIENT_ID", "OS_STATUS", "OS_MONTHS"]].copy()
    tnbc_patients.columns = ["patientId", "OS_STATUS", "OS_MONTHS"]
    print(f"TNBC patients (ER-/PR-/HER2-): {len(tnbc_patients)}")

    # 2. Expression matrix -> marker genes + whatever's being tested
    expr = pd.read_csv(f"{STUDY_DIR}/data_mrna_seq_v2_rsem.txt", sep="\t")
    expr_sub = expr[expr["Hugo_Symbol"].isin(genes_needed)].copy()
    expr_sub = expr_sub.drop(columns=["Entrez_Gene_Id"]).set_index("Hugo_Symbol")

    expr_t = expr_sub.T
    expr_t["patientId"] = expr_t.index.str[:12]
    expr_t["sampleId"] = expr_t.index
    expr_t = expr_t[expr_t["sampleId"].str.endswith("-01")]
    expr_t = expr_t.drop_duplicates(subset="patientId")

    missing = genes_needed - set(expr_sub.index)
    if missing:
        print(f"WARNING - genes not found in matrix: {missing}")

    tnbc_expr = expr_t[expr_t["patientId"].isin(tnbc_patients["patientId"])].copy()
    for g in genes_needed:
        if g in tnbc_expr.columns:
            tnbc_expr[g] = pd.to_numeric(tnbc_expr[g], errors="coerce")
    print(f"TNBC samples with expression data: {len(tnbc_expr)}")

    # 3. Marker-gene proxy for Lehmann subtype (BL1 / BL2 / M / LAR)
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


def analyze_gene(gene_symbol, cohort_df):
    """Run the per-gene analysis (subtype comparison, figure, survival split)
    against an already-built cohort. Writes <gene>_*.png/json/csv into
    ../results/ and returns the results dict."""
    if gene_symbol not in cohort_df.columns:
        print(f"SKIPPING {gene_symbol}: not found in expression matrix")
        return None

    gene_dir = f"{RESULTS_DIR}/{gene_symbol.lower()}"
    os.makedirs(gene_dir, exist_ok=True)

    df = cohort_df.copy()
    df[gene_symbol] = pd.to_numeric(df[gene_symbol], errors="coerce")

    groups = [df.loc[df["subtype"] == s, gene_symbol].dropna().values for s in SUBTYPE_ORDER]
    kw_stat, kw_p = stats.kruskal(*groups)
    print(f"\n[{gene_symbol}] Kruskal-Wallis across subtypes: H={kw_stat:.3f}, p={kw_p:.4f}")

    group_medians, group_means, group_n = {}, {}, {}
    for s, g in zip(SUBTYPE_ORDER, groups):
        group_medians[s] = float(np.median(g))
        group_means[s] = float(np.mean(g))
        group_n[s] = int(len(g))
        print(f"  {s}: n={len(g)}, median={np.median(g):.2f}, mean={np.mean(g):.2f}")

    pairwise = {}
    for a, b in combinations(SUBTYPE_ORDER, 2):
        ga = df.loc[df["subtype"] == a, gene_symbol].dropna().values
        gb = df.loc[df["subtype"] == b, gene_symbol].dropna().values
        u, p = stats.mannwhitneyu(ga, gb, alternative="two-sided")
        pairwise[f"{a}_vs_{b}"] = {"U": float(u), "p": float(p)}
        print(f"  {a} vs {b}: U={u:.1f}, p={p:.4f}")

    # Figure: expression by subtype
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.violinplot(data=df, x="subtype", y=gene_symbol, order=SUBTYPE_ORDER, palette="muted", ax=ax, cut=0)
    sns.stripplot(data=df, x="subtype", y=gene_symbol, order=SUBTYPE_ORDER, color="black", alpha=0.35, size=3, ax=ax)
    ax.set_xlabel("TNBC Subtype (Lehmann, marker-gene proxy)", fontsize=12)
    ax.set_ylabel(f"{gene_symbol} Expression (RNA-seq V2 RSEM)", fontsize=12)
    ax.set_title(f"{gene_symbol} Expression Across TNBC Subtypes", fontsize=13)
    ax.text(0.02, 0.98, f"Kruskal-Wallis p = {kw_p:.3f}", transform=ax.transAxes,
            va="top", fontsize=10, style="italic")
    plt.tight_layout()
    plt.savefig(f"{gene_dir}/{gene_symbol.lower()}_by_subtype.png", dpi=150)
    plt.close()
    print(f"Saved {gene_symbol.lower()}/{gene_symbol.lower()}_by_subtype.png")

    # Survival: high vs low (median split)
    median_val = df[gene_symbol].median()
    df["expr_group"] = np.where(df[gene_symbol] >= median_val, "high", "low")

    surv_df = df.dropna(subset=["OS_MONTHS", "OS_STATUS"]).copy()
    surv_df["event"] = surv_df["OS_STATUS"].astype(str).str.startswith("1").astype(int)

    os_high = surv_df.loc[surv_df["expr_group"] == "high", "OS_MONTHS"]
    os_low = surv_df.loc[surv_df["expr_group"] == "low", "OS_MONTHS"]
    os_u, os_p = stats.mannwhitneyu(os_high, os_low, alternative="two-sided")

    lr = logrank_test(
        surv_df.loc[surv_df["expr_group"] == "high", "OS_MONTHS"],
        surv_df.loc[surv_df["expr_group"] == "low", "OS_MONTHS"],
        event_observed_A=surv_df.loc[surv_df["expr_group"] == "high", "event"],
        event_observed_B=surv_df.loc[surv_df["expr_group"] == "low", "event"],
    )
    print(f"n evaluable for survival: {len(surv_df)} (high={len(os_high)}, low={len(os_low)})")
    print(f"Mann-Whitney on OS months, high vs low {gene_symbol}: p={os_p:.4f}")
    print(f"Log-rank test p={lr.p_value:.4f}")

    fig, ax = plt.subplots(figsize=(7, 5))
    kmf = KaplanMeierFitter()
    for label, mask in [(f"{gene_symbol} high", surv_df["expr_group"] == "high"),
                         (f"{gene_symbol} low", surv_df["expr_group"] == "low")]:
        kmf.fit(surv_df.loc[mask, "OS_MONTHS"], surv_df.loc[mask, "event"], label=label)
        kmf.plot_survival_function(ax=ax)
    ax.set_xlabel("Months")
    ax.set_ylabel("Overall survival probability")
    ax.set_title(f"TNBC overall survival by {gene_symbol} expression (log-rank p={lr.p_value:.3f})")
    plt.tight_layout()
    plt.savefig(f"{gene_dir}/{gene_symbol.lower()}_survival_km.png", dpi=150)
    plt.close()
    print(f"Saved {gene_symbol.lower()}/{gene_symbol.lower()}_survival_km.png")

    results = {
        "gene": gene_symbol,
        "n_tnbc_with_expression": int(len(df)),
        "subtype_counts": df["subtype"].value_counts().to_dict(),
        "kruskal_wallis": {"H": float(kw_stat), "p": float(kw_p)},
        "group_medians": group_medians,
        "group_means": group_means,
        "group_n": group_n,
        "pairwise_mannwhitney": pairwise,
        "median_split": float(median_val),
        "survival_mannwhitney_p": float(os_p),
        "survival_logrank_p": float(lr.p_value),
        "n_survival_evaluable": int(len(surv_df)),
    }
    with open(f"{gene_dir}/{gene_symbol.lower()}_analysis_results.json", "w") as f:
        json.dump(results, f, indent=2)

    out_df = df.drop(columns=["expr_group"], errors="ignore")
    out_df.to_csv(f"{gene_dir}/{gene_symbol.lower()}_tnbc_final.csv", index=False)

    print(f"\n=== DONE: {gene_symbol} ===")
    return results


if __name__ == "__main__":
    genes_to_run = sys.argv[1:] or ["GNRHR"]
    cohort = build_tnbc_cohort(extra_genes=genes_to_run)

    all_results = {}
    for gene in genes_to_run:
        result = analyze_gene(gene, cohort)
        if result is not None:
            all_results[gene] = result

    print("\n=== ALL GENES DONE ===")
    print(json.dumps(all_results, indent=2))
