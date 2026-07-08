# flow_cytometry/

Raw wet-lab flow cytometry measurements, hand-entered after each experiment. This is the protein-level validation step for Cleavr's computational predictions (RNA-level expression differences across TNBC subtypes).

## How to use

1. After each flow cytometry run, add one row per cell line x gene x replicate to `measurements.csv`.
2. Run `python3 ../code/flow_cytometry_validation.py` from `code/`. It normalizes MFI against the isotype control, groups by subtype, and compares the protein-level pattern against Cleavr's existing RNA-level result for that gene (if `biomarker_pipeline.py` has already been run for it).
3. Output lands in `../results/flow_cytometry_validation.json` and `.csv`.

## Columns

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
