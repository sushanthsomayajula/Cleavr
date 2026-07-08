"""
GNRHR (LHRH receptor) expression across TNBC subtypes -- run locally.

My Cowork sandbox is stuck (infra issue on Anthropic's side), so this is written
to run directly on your Mac with the Python setup you already have (same one that
ran the original notebook -- scikit-learn already installed there).

WHAT CHANGED FROM THE ORIGINAL NOTEBOOK:
The original notebook pulled GNRH1 (the ligand gene). Hu's nanoparticles target the
RECEPTOR on the tumor cell surface, encoded by GNRHR. This script analyzes GNRHR
instead, using your already-downloaded TCGA study bundle (brca_tcga/) directly --
no network calls, so it can't fail on API issues.

HOW TO RUN:
1. Open Terminal
2. cd into the folder you selected in Cowork (the one with brca_tcga/ inside), e.g.:
     cd /path/to/tnbc-project
   (Finder > right-click the folder > "New Terminal at Folder" is the easiest way)
3. Install anything missing (safe to run even if some are already installed):
     pip3 install pandas numpy scipy matplotlib seaborn scikit-learn lifelines --break-system-packages
4. Run it:
     python3 run_gnrhr_analysis.py

OUTPUT (written into this same folder):
  - gnrhr_by_subtype.png     <- main figure for the interview
  - gnrhr_survival_km.png    <- survival curve, GNRHR-high vs GNRHR-low
  - analysis_results.json    <- all the stats numbers
  - gnrhr_tnbc_final.csv     <- full per-patient table

Once this finishes, tell Claude in Cowork it's done -- it can read the JSON/CSV
directly from the shared folder (that part doesn't need the sandbox) and will
write up the one-page interview summary from the real numbers.
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


for pkg, imp in [("pandas", None), ("numpy", None), ("scipy", None),
                  ("matplotlib", None), ("seaborn", None),
                  ("scikit-learn", "sklearn"), ("lifelines", None)]:
    ensure(pkg, imp)

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import StandardScaler

STUDY_DIR = "../brca_tcga"

# ---------------------------------------------------------------
# 1. Clinical data -> TNBC cohort (ER-/PR-/HER2-)
# ---------------------------------------------------------------
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

# ---------------------------------------------------------------
# 2. Expression matrix -> GNRHR + 8 Lehmann marker genes
# ---------------------------------------------------------------
genes_needed = ["GNRHR", "CCNE1", "CDC6", "EGFR", "MET", "VIM", "ZEB1", "AR", "FOXA1"]

expr = pd.read_csv(f"{STUDY_DIR}/data_mrna_seq_v2_rsem.txt", sep="\t")
expr_sub = expr[expr["Hugo_Symbol"].isin(genes_needed)].copy()
expr_sub = expr_sub.drop(columns=["Entrez_Gene_Id"]).set_index("Hugo_Symbol")

expr_t = expr_sub.T
expr_t["patientId"] = expr_t.index.str[:12]
expr_t["sampleId"] = expr_t.index
expr_t = expr_t[expr_t["sampleId"].str.endswith("-01")]
expr_t = expr_t.drop_duplicates(subset="patientId")

missing = set(genes_needed) - set(expr_sub.index)
if missing:
    print(f"WARNING - genes not found in matrix: {missing}")

tnbc_expr = expr_t[expr_t["patientId"].isin(tnbc_patients["patientId"])].copy()
for g in genes_needed:
    tnbc_expr[g] = pd.to_numeric(tnbc_expr[g], errors="coerce")
print(f"TNBC samples with expression data: {len(tnbc_expr)}")

# ---------------------------------------------------------------
# 3. Marker-gene proxy for Lehmann subtype (BL1 / BL2 / M / LAR)
# ---------------------------------------------------------------
marker_genes = ["CCNE1", "CDC6", "EGFR", "MET", "VIM", "ZEB1", "AR", "FOXA1"]
scaler = StandardScaler()
scaled = scaler.fit_transform(tnbc_expr[marker_genes])
scaled_df = pd.DataFrame(scaled, columns=marker_genes, index=tnbc_expr.index)

scaled_df["BL1_score"] = scaled_df[["CCNE1", "CDC6"]].mean(axis=1)
scaled_df["BL2_score"] = scaled_df[["EGFR", "MET"]].mean(axis=1)
scaled_df["M_score"] = scaled_df[["VIM", "ZEB1"]].mean(axis=1)
scaled_df["LAR_score"] = scaled_df[["AR", "FOXA1"]].mean(axis=1)

score_cols = ["BL1_score", "BL2_score", "M_score", "LAR_score"]
tnbc_expr["subtype"] = scaled_df[score_cols].idxmax(axis=1).str.replace("_score", "")
print(tnbc_expr["subtype"].value_counts())

# ---------------------------------------------------------------
# 4. Merge survival, run stats
# ---------------------------------------------------------------
final_df = tnbc_expr.merge(tnbc_patients, on="patientId", how="left")
final_df["GNRHR"] = pd.to_numeric(final_df["GNRHR"], errors="coerce")
final_df["OS_MONTHS"] = pd.to_numeric(final_df["OS_MONTHS"], errors="coerce")

order = ["BL1", "BL2", "M", "LAR"]
groups = [final_df.loc[final_df["subtype"] == s, "GNRHR"].dropna().values for s in order]
kw_stat, kw_p = stats.kruskal(*groups)
print(f"\nKruskal-Wallis across subtypes: H={kw_stat:.3f}, p={kw_p:.4f}")

group_medians, group_means, group_n = {}, {}, {}
for s, g in zip(order, groups):
    group_medians[s] = float(np.median(g))
    group_means[s] = float(np.mean(g))
    group_n[s] = int(len(g))
    print(f"  {s}: n={len(g)}, median={np.median(g):.2f}, mean={np.mean(g):.2f}")

from itertools import combinations
pairwise = {}
for a, b in combinations(order, 2):
    ga = final_df.loc[final_df["subtype"] == a, "GNRHR"].dropna().values
    gb = final_df.loc[final_df["subtype"] == b, "GNRHR"].dropna().values
    u, p = stats.mannwhitneyu(ga, gb, alternative="two-sided")
    pairwise[f"{a}_vs_{b}"] = {"U": float(u), "p": float(p)}
    print(f"  {a} vs {b}: U={u:.1f}, p={p:.4f}")

# ---------------------------------------------------------------
# 5. Figure
# ---------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

fig, ax = plt.subplots(figsize=(8, 5))
sns.violinplot(data=final_df, x="subtype", y="GNRHR", order=order, palette="muted", ax=ax, cut=0)
sns.stripplot(data=final_df, x="subtype", y="GNRHR", order=order, color="black", alpha=0.35, size=3, ax=ax)
ax.set_xlabel("TNBC Subtype (Lehmann, marker-gene proxy)", fontsize=12)
ax.set_ylabel("GNRHR Expression (RNA-seq V2 RSEM)", fontsize=12)
ax.set_title("GNRHR (LHRH Receptor) Expression Across TNBC Subtypes", fontsize=13)
ax.text(0.02, 0.98, f"Kruskal-Wallis p = {kw_p:.3f}", transform=ax.transAxes,
        va="top", fontsize=10, style="italic")
plt.tight_layout()
plt.savefig("../results/gnrhr_by_subtype.png", dpi=150)
plt.close()
print("\nSaved gnrhr_by_subtype.png")

# ---------------------------------------------------------------
# 6. Survival: GNRHR-high vs GNRHR-low (median split)
# ---------------------------------------------------------------
median_gnrhr = final_df["GNRHR"].median()
final_df["gnrhr_group"] = np.where(final_df["GNRHR"] >= median_gnrhr, "high", "low")

surv_df = final_df.dropna(subset=["OS_MONTHS", "OS_STATUS"]).copy()
surv_df["event"] = surv_df["OS_STATUS"].astype(str).str.startswith("1").astype(int)

os_high = surv_df.loc[surv_df["gnrhr_group"] == "high", "OS_MONTHS"]
os_low = surv_df.loc[surv_df["gnrhr_group"] == "low", "OS_MONTHS"]
os_u, os_p = stats.mannwhitneyu(os_high, os_low, alternative="two-sided")

from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

lr = logrank_test(
    surv_df.loc[surv_df["gnrhr_group"] == "high", "OS_MONTHS"],
    surv_df.loc[surv_df["gnrhr_group"] == "low", "OS_MONTHS"],
    event_observed_A=surv_df.loc[surv_df["gnrhr_group"] == "high", "event"],
    event_observed_B=surv_df.loc[surv_df["gnrhr_group"] == "low", "event"],
)
print(f"\nn evaluable for survival: {len(surv_df)} (high={len(os_high)}, low={len(os_low)})")
print(f"Mann-Whitney on OS months, high vs low GNRHR: p={os_p:.4f}")
print(f"Log-rank test p={lr.p_value:.4f}")

fig, ax = plt.subplots(figsize=(7, 5))
kmf = KaplanMeierFitter()
for label, mask in [("GNRHR high", surv_df["gnrhr_group"] == "high"),
                     ("GNRHR low", surv_df["gnrhr_group"] == "low")]:
    kmf.fit(surv_df.loc[mask, "OS_MONTHS"], surv_df.loc[mask, "event"], label=label)
    kmf.plot_survival_function(ax=ax)
ax.set_xlabel("Months")
ax.set_ylabel("Overall survival probability")
ax.set_title(f"TNBC overall survival by GNRHR expression (log-rank p={lr.p_value:.3f})")
plt.tight_layout()
plt.savefig("../results/gnrhr_survival_km.png", dpi=150)
plt.close()
print("Saved gnrhr_survival_km.png")

# ---------------------------------------------------------------
# 7. Save results
# ---------------------------------------------------------------
results = {
    "n_tnbc_patients": int(len(tnbc_patients)),
    "n_tnbc_with_expression": int(len(final_df)),
    "subtype_counts": final_df["subtype"].value_counts().to_dict(),
    "kruskal_wallis": {"H": float(kw_stat), "p": float(kw_p)},
    "group_medians": group_medians,
    "group_means": group_means,
    "group_n": group_n,
    "pairwise_mannwhitney": pairwise,
    "median_gnrhr_split": float(median_gnrhr),
    "survival_mannwhitney_p": float(os_p),
    "survival_logrank_p": float(lr.p_value),
    "n_survival_evaluable": int(len(surv_df)),
}
with open("../results/analysis_results.json", "w") as f:
    json.dump(results, f, indent=2)

final_df.to_csv("../results/gnrhr_tnbc_final.csv", index=False)

print("\n=== DONE ===")
print(json.dumps(results, indent=2))
