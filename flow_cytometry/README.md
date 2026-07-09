# flow_cytometry/

Raw wet-lab flow cytometry measurements, hand-entered after each experiment. This is the protein-level validation step for Cleavr's computational predictions (RNA-level expression differences across TNBC subtypes).

Two tools work off this data: `flow_cytometry_validation.py` compares a relative signal (MFI/isotype fold-change) across subtypes, `receptor_quantification.py` converts that into an absolute receptors-per-cell estimate using a calibration bead standard curve. The relative tool needs only `measurements.csv`; the quantification tool also needs `calibration_beads.csv` from the same experiment day.

For entering and reviewing data without editing CSVs directly, `code/qflow_tool.html` is a self-contained browser tool (no server, no dependencies, open the file directly) that does both steps in one place: editable data-entry tables, calibration curve fitting, ABC conversion, and the subtype comparison, plus CSV import/export using the exact same column schema as the two scripts above, so files move freely between the browser tool and the command line. It also fixes a gap in `flow_cytometry_validation.py`: that script only ever uses the raw MFI/isotype ratio and never actually picks up the calibrated receptor count, even though this file used to say it would. The browser tool uses the calibrated number whenever a calibration curve exists for that date, falling back to the raw ratio otherwise.

## How to use

1. Run calibration beads on the flow cytometer the same day as your samples, same instrument settings. Add one row per bead population (need at least 3, kits typically ship 4-6) to `calibration_beads.csv`, using the molecules-per-bead value from the kit's lot-specific certificate of analysis, not a measured value.
2. After each flow cytometry run, add one row per cell line x gene x replicate to `measurements.csv`.
3. Run `python3 ../code/receptor_quantification.py` from `code/`. It fits a calibration curve (linear or log-log, whichever fits better) per date and converts each sample's MFI into a background-subtracted receptors-per-cell estimate. Output lands in `../results/receptor_quantification.json` and `.csv`.
4. Run `python3 ../code/flow_cytometry_validation.py`. It normalizes MFI against the isotype control, groups by subtype, and compares the protein-level pattern against Cleavr's existing RNA-level result for that gene (if `biomarker_pipeline.py` has already been run for it). Note: this script always uses the raw MFI/isotype ratio, even on dates with a fitted calibration curve; `code/qflow_tool.html` uses the calibrated receptor count instead whenever one is available.

## Columns

**`measurements.csv`**

| Column | Meaning |
|---|---|
| `date` | Experiment date |
| `cell_line` | Must be one of the cell lines in the lookup table below (case/hyphen-insensitive) |
| `gene` | Gene symbol, matching what was run through `biomarker_pipeline.py` |
| `mfi` | Median fluorescence intensity, target antibody |
| `isotype_mfi` | Median fluorescence intensity, isotype control (must be > 0) |
| `pct_positive` | Percent positive cells (recorded, not currently used in the comparison) |
| `replicate` | Replicate number, for averaging |
| `notes` | Free text |

**`calibration_beads.csv`**

| Column | Meaning |
|---|---|
| `date` | Calibration run date, must match the `measurements.csv` date it calibrates |
| `bead_lot` | Bead kit lot number |
| `population` | Bead population label (e.g. peak1, peak2) |
| `molecules_per_bead` | Certified antibody-binding capacity for this population, from the kit's certificate of analysis |
| `mfi` | Measured MFI for this bead population |

## Why cell lines need their own subtype table

Flow cytometry uses cell lines, not patient samples, so there's no Lehmann subtype label to read off a clinical record the way `biomarker_pipeline.py` does for TCGA/METABRIC. The mapping in `code/flow_cytometry_validation.py` is pulled from two papers, not guessed:

- Lehmann BD et al. J Clin Invest. 2011;121(7):2750-67. PMID 21633166 (original TNBCtype algorithm applied to a cell line panel)
- Espinosa Fernandez JR et al. PLoS One. 2020;15(4):e0231953. PMID 32353087 (cross-checked against a revised algorithm, in vitro and in vivo, across 5 datasets)

| Cell line | Subtype | Confidence |
|---|---|---|
| HCC70 | BL2 | Stable (concordant on every check) |
| SUM149PT | BL2 | Stable |
| HCC1806 | BL2 | Stable |
| BT549 | M | Stable |
| MDA-MB-453 | LAR | Stable |
| HCC2157 | BL1 | Stable |
| MDA-MB-468 | BL1 | Single-algorithm call, less certain |
| HCC1937 | BL1 | Single-algorithm call, less certain |
| HCC3153 | BL1 | Single-algorithm call, less certain |
| SUM185PE | LAR | Single-algorithm call, less certain |

Cell lines classified as MSL or IM subtype in the literature are excluded, since Cleavr's pipeline only models the 4-subtype BL1/BL2/M/LAR scheme. A cell line not in this table gets excluded from the comparison with a warning printed at run time, rather than silently dropped or guessed at.
