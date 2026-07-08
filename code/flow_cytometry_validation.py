"""
Cleavr flow cytometry validation -- run locally, whenever new flow
cytometry data exists.

Cleavr's screen and cross-cohort validation are computational: RNA-level
expression differences across TNBC subtypes, tested in silico against
public sequencing data. Flow cytometry checks whether that pattern holds
at the protein level, in real TNBC cell lines -- which is what actually
determines whether a receptor is usable for nanoparticle targeting. RNA
and protein levels do not always agree, so this comparison is a real check,
not a formality.

WHY CELL LINES NEED THEIR OWN SUBTYPE TABLE: flow cytometry measures cell
lines, not patient samples, so there is no Lehmann subtype label to read
off a clinical record the way biomarker_pipeline.py does for TCGA/METABRIC.
The subtype assignments below are pulled directly from the literature, not
guessed:
  - Lehmann BD et al. "Identification of human triple-negative breast
    cancer subtypes and preclinical models for selection of targeted
    therapies." J Clin Invest. 2011;121(7):2750-67. PMID 21633166.
    (original TNBCtype algorithm, applied to a panel of TNBC cell lines)
  - Espinosa Fernandez JR et al. "Identification of triple-negative breast
    cancer cell lines classified under the same molecular subtype using
    different molecular characterization techniques: Implications for
    translational research." PLoS One. 2020;15(4):e0231953. PMID 32353087.
    (cross-checked the same cell lines against a second, revised algorithm
    -- TNBCtype-IM -- in vitro AND in vivo as xenografts, across 5
    independent datasets. The six cell lines flagged "stable" below were
    concordant on every single check; the rest were only classified once.)

HOW TO RUN:
    cd ~/Desktop/tnbc-project/code
    python3 flow_cytometry_validation.py
(First run with no data yet creates an empty template and stops -- fill it
in after each experiment, then re-run.)

INPUT (hand-maintained, append a row after every experiment):
    ../flow_cytometry/measurements.csv
    columns: date, cell_line, gene, mfi, isotype_mfi, pct_positive, replicate, notes

OUTPUT:
    ../results/flow_cytometry_validation.json
    ../results/flow_cytometry_validation.csv
"""
import csv
import json
import os
import re
from collections import defaultdict

FLOW_DIR = "../flow_cytometry"
RESULTS_DIR = "../results"
INPUT_CSV = f"{FLOW_DIR}/measurements.csv"

# Cell line -> Lehmann subtype (BL1/BL2/M/LAR only; MSL and IM cell lines
# are excluded since Cleavr's pipeline only models the 4-subtype scheme
# used for TCGA/METABRIC). See module docstring for the two source papers.
CELL_LINE_SUBTYPE = {
    "HCC70":    {"subtype": "BL2", "stable": True},
    "SUM149PT": {"subtype": "BL2", "stable": True},
    "HCC1806":  {"subtype": "BL2", "stable": True},
    "BT549":    {"subtype": "M",   "stable": True},
    "MDAMB453": {"subtype": "LAR", "stable": True},
    "HCC2157":  {"subtype": "BL1", "stable": True},
    "MDAMB468": {"subtype": "BL1", "stable": False},
    "HCC1937":  {"subtype": "BL1", "stable": False},
    "HCC3153":  {"subtype": "BL1", "stable": False},
    "SUM185PE": {"subtype": "LAR", "stable": False},
}


def normalize_cell_line(name):
    """Strip hyphens/spaces/case so 'MDA-MB-453', 'mda mb 453', and
    'MDAMB453' all match the same lookup key."""
    return re.sub(r"[^A-Z0-9]", "", name.upper())


NORMALIZED_LOOKUP = {normalize_cell_line(k): v for k, v in CELL_LINE_SUBTYPE.items()}


def load_measurements():
    if not os.path.exists(INPUT_CSV):
        return None
    rows = []
    with open(INPUT_CSV) as f:
        for row in csv.DictReader(f):
            if not row.get("cell_line") or not row.get("gene"):
                continue
            rows.append(row)
    return rows


def load_rna_pattern(gene):
    """Pull the TCGA subtype medians already computed by biomarker_pipeline.py
    for this gene, if it's been run for that gene yet."""
    path = f"{RESULTS_DIR}/{gene.lower()}/{gene.lower()}_analysis_results.json"
    if not os.path.exists(path):
        return None
    with open(path) as f:
        r = json.load(f)
    return r.get("group_medians")


if __name__ == "__main__":
    rows = load_measurements()

    if rows is None:
        os.makedirs(FLOW_DIR, exist_ok=True)
        with open(INPUT_CSV, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["date", "cell_line", "gene", "mfi", "isotype_mfi",
                              "pct_positive", "replicate", "notes"])
        print(f"No measurements yet -- created an empty template at {INPUT_CSV}")
        print("Add one row per experiment (cell line x gene x replicate), then re-run this script.")
        print(f"Recognized cell lines: {sorted(CELL_LINE_SUBTYPE)}")
        raise SystemExit(0)

    if len(rows) == 0:
        print(f"{INPUT_CSV} exists but has no data rows yet. Add measurements and re-run.")
        raise SystemExit(0)

    print(f"Loaded {len(rows)} measurement rows from {INPUT_CSV}")

    by_gene_cellline = defaultdict(list)
    unmapped_cell_lines = set()
    for row in rows:
        cl_key = normalize_cell_line(row["cell_line"])
        gene = row["gene"].strip().upper()
        try:
            mfi = float(row["mfi"])
            isotype = float(row["isotype_mfi"])
        except (ValueError, KeyError):
            print(f"  SKIPPING row (bad mfi/isotype_mfi): {row}")
            continue
        if isotype <= 0:
            print(f"  SKIPPING row (isotype_mfi must be > 0): {row}")
            continue
        signal = mfi / isotype  # fold-over-isotype, the actual protein-level readout
        info = NORMALIZED_LOOKUP.get(cl_key)
        if info is None:
            unmapped_cell_lines.add(row["cell_line"])
            continue
        by_gene_cellline[(gene, row["cell_line"], info["subtype"], info["stable"])].append(signal)

    if unmapped_cell_lines:
        print(f"\nWARNING - cell line(s) not in the subtype lookup table, excluded: {sorted(unmapped_cell_lines)}")
        print("Add them to CELL_LINE_SUBTYPE in this script (with a literature citation) to include them.")

    per_gene = defaultdict(lambda: defaultdict(list))
    for (gene, cell_line, subtype, stable), signals in by_gene_cellline.items():
        mean_signal = sum(signals) / len(signals)
        per_gene[gene][subtype].append({
            "cell_line": cell_line, "mean_signal": mean_signal,
            "n_replicates": len(signals), "stable_classification": stable,
        })

    report = []
    for gene, subtype_data in per_gene.items():
        print(f"\n=== {gene} ===")
        protein_by_subtype = {}
        for subtype, lines in subtype_data.items():
            avg = sum(l["mean_signal"] for l in lines) / len(lines)
            protein_by_subtype[subtype] = avg
            for l in lines:
                note = "" if l["stable_classification"] else " (single-algorithm classification, less certain)"
                print(f"  {subtype} / {l['cell_line']}: signal={l['mean_signal']:.2f} (n={l['n_replicates']} replicates){note}")

        rna_medians = load_rna_pattern(gene)
        entry = {
            "gene": gene,
            "protein_signal_by_subtype": protein_by_subtype,
            "n_subtypes_with_data": len(protein_by_subtype),
            "rna_medians_tcga": rna_medians,
        }

        if rna_medians is None:
            print(f"  No RNA-level result found for {gene} (run biomarker_pipeline.py {gene} first). Protein data recorded, no comparison yet.")
            entry["comparison"] = None
        elif len(protein_by_subtype) < 2:
            print(f"  Only {len(protein_by_subtype)} subtype(s) have flow data so far; need at least 2 to compare a pattern.")
            entry["comparison"] = "insufficient_subtypes"
        else:
            common_subtypes = sorted(set(protein_by_subtype) & set(rna_medians))
            protein_rank = sorted(common_subtypes, key=lambda s: protein_by_subtype[s])
            rna_rank = sorted(common_subtypes, key=lambda s: rna_medians[s])
            concordant = protein_rank == rna_rank
            print(f"  Protein-level ranking (low to high): {protein_rank}")
            print(f"  RNA-level ranking (TCGA, low to high): {rna_rank}")
            print(f"  Rankings match exactly: {concordant}")
            entry["comparison"] = {
                "subtypes_compared": common_subtypes,
                "protein_rank_low_to_high": protein_rank,
                "rna_rank_low_to_high": rna_rank,
                "exact_rank_match": concordant,
            }
        report.append(entry)

    with open(f"{RESULTS_DIR}/flow_cytometry_validation.json", "w") as f:
        json.dump(report, f, indent=2)

    with open(f"{RESULTS_DIR}/flow_cytometry_validation.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["gene", "n_subtypes_with_data", "exact_rank_match", "protein_rank_low_to_high", "rna_rank_low_to_high"])
        for e in report:
            comp = e["comparison"]
            if isinstance(comp, dict):
                writer.writerow([e["gene"], e["n_subtypes_with_data"], comp["exact_rank_match"],
                                  ">".join(comp["protein_rank_low_to_high"]), ">".join(comp["rna_rank_low_to_high"])])
            else:
                writer.writerow([e["gene"], e["n_subtypes_with_data"], comp, "", ""])

    print(f"\nSaved {RESULTS_DIR}/flow_cytometry_validation.json and .csv")
