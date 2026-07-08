"""
Cleavr absolute receptor quantification -- run locally, before or alongside
flow_cytometry_validation.py.

flow_cytometry_validation.py compares a relative signal (MFI / isotype MFI)
across subtypes. That tells you which cell line has more receptor than
another, but not an actual number. This script does the real thing:
converts raw fluorescence into an estimated antibody-binding capacity (ABC),
essentially receptors per cell, using calibration beads with a certified,
known number of binding sites -- the standard method in quantitative flow
cytometry (e.g. Bangs Labs Quantum Simply Cellular, BD Quantibrite).

METHOD: different bead kits recommend different curve shapes (some linear,
some log-log). Rather than assume one, this fits both a linear regression
(ABC = a*MFI + b) and a log-log regression (log10(ABC) = a*log10(MFI) + b)
per calibration run and keeps whichever fits better (higher R^2), reporting
that choice and the R^2 so it's auditable, not just trusted. Both fits are
implemented directly (ordinary least squares) rather than via a stats
library, for the same auditability reason the BH-FDR correction in
screen_candidates.py is implemented by hand.

WHY PER-DATE CALIBRATION: PMT voltage and instrument settings can vary
between runs, so a calibration curve from one day doesn't necessarily apply
to another. Calibration beads should be run on the same day, same settings,
as the samples they calibrate -- this script fits one curve per date in
calibration_beads.csv and only applies it to that date's measurements.

HOW TO RUN:
    cd ~/Desktop/tnbc-project/code
    python3 receptor_quantification.py
(First run with no calibration data yet creates an empty template and
stops -- fill it in, then re-run. Needs at least 3 bead populations per
date to fit a curve; kits typically ship 4-6.)

INPUT:
    ../flow_cytometry/calibration_beads.csv
    columns: date, bead_lot, population, molecules_per_bead, mfi
    ../flow_cytometry/measurements.csv  (same file flow_cytometry_validation.py reads)

OUTPUT:
    ../results/receptor_quantification.json
    ../results/receptor_quantification.csv
    (flow_cytometry_validation.py automatically picks this up on its next
    run and uses net_receptors_per_cell instead of the raw MFI ratio,
    wherever a matching date/cell_line/gene/replicate exists.)
"""
import csv
import json
import math
import os
from collections import defaultdict

FLOW_DIR = "../flow_cytometry"
RESULTS_DIR = "../results"
BEADS_CSV = f"{FLOW_DIR}/calibration_beads.csv"
MEASUREMENTS_CSV = f"{FLOW_DIR}/measurements.csv"

MIN_BEAD_POPULATIONS = 3


def ols_fit(x, y):
    """Ordinary least squares, y = slope*x + intercept. Returns (slope,
    intercept, r_squared). Implemented directly, not via numpy/scipy, so
    the math is easy to audit line by line."""
    n = len(x)
    mean_x = sum(x) / n
    mean_y = sum(y) / n
    num = sum((xi - mean_x) * (yi - mean_y) for xi, yi in zip(x, y))
    den = sum((xi - mean_x) ** 2 for xi in x)
    if den == 0:
        return None
    slope = num / den
    intercept = mean_y - slope * mean_x
    y_pred = [slope * xi + intercept for xi in x]
    ss_res = sum((yi - ypi) ** 2 for yi, ypi in zip(y, y_pred))
    ss_tot = sum((yi - mean_y) ** 2 for yi in y)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return slope, intercept, r2


def fit_calibration_curve(points):
    """points: list of (mfi, molecules_per_bead). Tries both a linear fit
    and a log-log fit, keeps whichever has the higher R^2. Returns a dict
    describing the chosen curve, or None if there aren't enough points."""
    if len(points) < MIN_BEAD_POPULATIONS:
        return None

    mfis = [p[0] for p in points]
    abcs = [p[1] for p in points]

    linear = ols_fit(mfis, abcs)
    candidates = []
    if linear:
        slope, intercept, r2 = linear
        candidates.append({"method": "linear", "slope": slope, "intercept": intercept, "r_squared": r2})

    if all(v > 0 for v in mfis) and all(v > 0 for v in abcs):
        log_mfis = [math.log10(v) for v in mfis]
        log_abcs = [math.log10(v) for v in abcs]
        loglog = ols_fit(log_mfis, log_abcs)
        if loglog:
            slope, intercept, r2 = loglog
            candidates.append({"method": "log-log", "slope": slope, "intercept": intercept, "r_squared": r2})

    if not candidates:
        return None
    return max(candidates, key=lambda c: c["r_squared"])


def apply_curve(curve, mfi):
    if mfi is None or mfi <= 0:
        return None
    if curve["method"] == "linear":
        return curve["slope"] * mfi + curve["intercept"]
    log_abc = curve["slope"] * math.log10(mfi) + curve["intercept"]
    return 10 ** log_abc


def load_bead_data():
    if not os.path.exists(BEADS_CSV):
        return None
    by_date = defaultdict(list)
    with open(BEADS_CSV) as f:
        for row in csv.DictReader(f):
            if not row.get("date") or not row.get("mfi") or not row.get("molecules_per_bead"):
                continue
            try:
                mfi = float(row["mfi"])
                abc = float(row["molecules_per_bead"])
            except ValueError:
                continue
            by_date[row["date"]].append((mfi, abc))
    return by_date


def load_measurement_rows():
    if not os.path.exists(MEASUREMENTS_CSV):
        return []
    with open(MEASUREMENTS_CSV) as f:
        return [row for row in csv.DictReader(f) if row.get("cell_line") and row.get("gene")]


if __name__ == "__main__":
    by_date = load_bead_data()

    if by_date is None:
        os.makedirs(FLOW_DIR, exist_ok=True)
        with open(BEADS_CSV, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["date", "bead_lot", "population", "molecules_per_bead", "mfi"])
        print(f"No calibration bead data yet -- created an empty template at {BEADS_CSV}")
        print("Add one row per bead population per calibration run (need >= 3 populations per date), then re-run.")
        print("molecules_per_bead comes from the bead kit's lot-specific certificate of analysis, not measured.")
        raise SystemExit(0)

    if not by_date:
        print(f"{BEADS_CSV} exists but has no usable rows yet. Add bead measurements and re-run.")
        raise SystemExit(0)

    curves = {}
    for date, points in by_date.items():
        curve = fit_calibration_curve(points)
        if curve is None:
            print(f"WARNING - {date}: only {len(points)} bead population(s), need >= {MIN_BEAD_POPULATIONS}. Skipping this date.")
            continue
        curves[date] = curve
        print(f"{date}: {curve['method']} calibration curve, R^2 = {curve['r_squared']:.4f}, "
              f"from {len(points)} bead populations")
        if curve["r_squared"] < 0.95:
            print(f"  WARNING - R^2 below 0.95, this curve is noisier than kit manufacturers typically expect. Re-run the calibration if possible.")

    measurement_rows = load_measurement_rows()
    if not measurement_rows:
        print(f"\nNo rows in {MEASUREMENTS_CSV} yet. Calibration curves are fit and saved; "
              f"quantification will run once measurement data exists.")

    results = []
    for row in measurement_rows:
        date = row.get("date", "")
        curve = curves.get(date)
        if curve is None:
            continue
        try:
            mfi = float(row["mfi"])
            isotype_mfi = float(row["isotype_mfi"])
        except (ValueError, KeyError):
            continue

        abc_target = apply_curve(curve, mfi)
        abc_isotype = apply_curve(curve, isotype_mfi)
        if abc_target is None or abc_isotype is None:
            continue
        net = abc_target - abc_isotype

        results.append({
            "date": date,
            "cell_line": row["cell_line"],
            "gene": row["gene"].strip().upper(),
            "replicate": row.get("replicate", ""),
            "mfi": mfi,
            "isotype_mfi": isotype_mfi,
            "abc_target": abc_target,
            "abc_isotype": abc_isotype,
            "net_receptors_per_cell": net,
            "calibration_method": curve["method"],
            "calibration_r_squared": curve["r_squared"],
        })

    with open(f"{RESULTS_DIR}/receptor_quantification.json", "w") as f:
        json.dump(results, f, indent=2)

    with open(f"{RESULTS_DIR}/receptor_quantification.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["date", "cell_line", "gene", "replicate", "mfi", "isotype_mfi",
                          "abc_target", "abc_isotype", "net_receptors_per_cell",
                          "calibration_method", "calibration_r_squared"])
        for r in results:
            writer.writerow([r["date"], r["cell_line"], r["gene"], r["replicate"],
                              f"{r['mfi']:.1f}", f"{r['isotype_mfi']:.1f}",
                              f"{r['abc_target']:.0f}", f"{r['abc_isotype']:.0f}",
                              f"{r['net_receptors_per_cell']:.0f}",
                              r["calibration_method"], f"{r['calibration_r_squared']:.4f}"])

    print(f"\n=== QUANTIFICATION COMPLETE: {len(results)} measurement rows converted ===")
    for r in results:
        print(f"  {r['date']} {r['cell_line']} {r['gene']}: ~{r['net_receptors_per_cell']:.0f} receptors/cell "
              f"({r['calibration_method']} fit, R^2={r['calibration_r_squared']:.3f})")
    print(f"\nSaved {RESULTS_DIR}/receptor_quantification.json and .csv")
    print("Run flow_cytometry_validation.py next -- it will use these numbers automatically.")
