# Cleavr

Systematic biomarker discovery for triple-negative breast cancer (TNBC), subtype by subtype.

**Live site:** https://cleavr.bio

Cleavr tests candidate TNBC biomarkers against public, patient-level sequencing data, broken out by molecular subtype, and cross-checks every claim against the primary literature before publishing it. Every number is independently re-derived from raw source files in a second pass. Null results are reported as directly as positive ones.

This repo currently holds two studies. **Case Study 01** asks whether GNRHR (the LHRH/GnRH receptor, the target of an active TNBC nanoparticle-targeting strategy) expresses differently across TNBC's molecular subtypes: no statistically significant difference (Kruskal-Wallis p = 0.71), and no relationship to survival (log-rank p = 0.98) in this cohort (n = 116). **Case Study 02** is a systematic 15-gene screen with FDR correction, cross-validated against a second cohort (METABRIC): 7 of 8 FDR-significant genes replicated. See `website/studies.html` for both.

## Structure

```
Cleavr/
├── README.md                  <- this file
├── index.html                  <- redirects to website/index.html (GitHub Pages entry point)
├── website/
│   ├── index.html               <- home / project pitch
│   ├── studies.html              <- index of every case study and screen
│   ├── gnrhr.html                <- Case Study 01: full write-up + figures
│   ├── explorer.html             <- interactive tool: pick a subtype, see expression data + a live 3D receptor structure
│   ├── tools.html                <- Tools section: index of the open lab tools (linked from the main nav)
│   ├── qflow.html                <- qFlow: quantitative flow cytometry workbench (6 modules, browser-local)
│   └── README.md                 <- short note about this subfolder
├── code/
│   ├── biomarker_pipeline.py   <- run this to reproduce/extend the analysis for any gene
│   ├── source_candidates.py    <- systematic candidate-gene sourcing (Open Targets)
│   ├── candidate_scoring.py    <- druggability (ChEMBL) + literature-gap (PubMed) scoring
│   ├── screen_candidates.py    <- batch screen with Benjamini-Hochberg FDR correction
│   ├── fetch_metabric.py       <- pulls the METABRIC cohort (cBioPortal API) for cross-cohort validation
│   ├── cross_cohort_validation.py  <- re-tests TCGA-significant genes against METABRIC
│   ├── flow_cytometry_validation.py  <- compares protein-level (flow) signal against the RNA-level result
│   ├── receptor_quantification.py    <- converts flow MFI into absolute receptors/cell via bead calibration
│   ├── qflow_tool.html        <- dev copy of the qFlow workbench, kept byte-identical to website/qflow.html
│   └── tnbc-analysis.ipynb     <- original Case Study 01 analysis as a notebook, with narration
├── flow_cytometry/
│   ├── measurements.csv        <- hand-entered flow cytometry data, one row per cell line x gene x replicate
│   ├── calibration_beads.csv   <- bead standard data for absolute quantification
│   └── README.md               <- column docs + cell-line-to-subtype table, sourced from literature
├── results/
│   ├── <gene>/                 <- one folder per gene tested (e.g. results/gnrhr/, results/fto/), each with
│   │                              <gene>_by_subtype.png, <gene>_survival_km.png, <gene>_tnbc_final.csv,
│   │                              <gene>_analysis_results.json
│   ├── candidate_gene_list.json    <- systematically sourced candidate list
│   ├── candidate_scorecard.json    <- druggability + literature-gap scores
│   ├── screen_results.json/.csv    <- full 15-gene screen, FDR-corrected
│   └── cross_cohort_results.json/.csv  <- TCGA vs. METABRIC replication results
├── test/
│   ├── qflow_core.test.js     <- extracts qFlow's compute core and checks it (96 assertions)
│   ├── py_parity.py           <- runs the real code/ Python functions for cross-language parity
│   ├── test_data.json         <- shared calibration/measurement dataset for both sides
│   └── README.md              <- how to run the suite
├── docs/
│   ├── GNRHR_TNBC_Summary.txt  <- plain text write-up
│   ├── GNRHR_TNBC_Summary.md   <- same write-up, Markdown formatted
│   ├── ROADMAP.md              <- technical roadmap, what's done and what's next
│   └── UNDERSTANDING_THE_FINDINGS.txt  <- plain-language walkthrough of every concept and number
├── brca_metabric/                <- gitignored; METABRIC cohort data, re-fetchable via code/fetch_metabric.py
└── brca_tcga/                    <- gitignored; full raw TCGA download (~1.3 GB), not pushed
```

The large raw cohort data (`brca_tcga/`, ~1.3 GB, and `brca_metabric/`) is intentionally excluded from this repo (see `.gitignore`): it's re-downloadable from cBioPortal and doesn't belong in git. `code/biomarker_pipeline.py` and `code/tnbc-analysis.ipynb` expect `../brca_tcga` to exist locally if you want to reproduce the analysis from scratch:

```
cd code
python3 biomarker_pipeline.py GNRHR AR
```

Pass any gene symbol(s) present in the expression matrix; the cohort and subtype assignment are built once and reused, only the target gene changes.

## Methodology notes (see `docs/` for the full write-up)

- TNBC cohort (n=116) = ER-negative, PR-negative, HER2-negative by IHC. 115 had usable RNA-seq data.
- Subtype (BL1/BL2/M/LAR) is a **marker-gene proxy**, not the validated Lehmann classifier: this TCGA study doesn't ship real subtype calls.
- Kruskal-Wallis across subtypes: p = 0.71 (not significant). Pairwise comparisons remain non-significant even after Bonferroni correction (all p = 1.0).
- Survival analysis (log-rank p = 0.98) is based on only 18 deaths out of 115 patients (15.7% event rate), underpowered, treat as inconclusive rather than a confident null.
- All numbers above were independently re-derived from the raw TCGA files in a second pass and cross-checked against the saved results; they match exactly.

## Roadmap

**Done:** parameterized pipeline (`code/biomarker_pipeline.py`, any gene, not just GNRHR); systematic candidate sourcing from Open Targets (`code/source_candidates.py`); batch screening with Benjamini-Hochberg FDR correction across all genes tested together (`code/screen_candidates.py`); druggability (ChEMBL) and TNBC-specific literature-gap (PubMed) scoring per candidate (`code/candidate_scoring.py`); cross-cohort validation against METABRIC (`code/fetch_metabric.py`, `code/cross_cohort_validation.py`), 7 of 8 FDR-significant genes replicated; flow cytometry wet-lab validation path, both relative signal (`code/flow_cytometry_validation.py`) and absolute receptors-per-cell via bead calibration (`code/receptor_quantification.py`). A public **Tools** section (`website/tools.html`, linked from the main nav) now hosts **qFlow** (`website/qflow.html`), a browser-local quantitative-flow workbench covering the receptor-quantification family: calibration/ABC, MESF/ERF/PE unit conversion, %positive gating, stain/separation index, spillover compensation, and FCS 3.0/3.1 parsing. Its compute core is checked by an automated suite (`test/`, 96 assertions) that verifies numeric parity with the Python scripts. See `docs/ROADMAP.md` for the full writeup, including honest caveats on the screen and validation.

**Still open:**
- Actually running the first real flow cytometry experiment on the top candidate (FTO): the tooling is built and tested against synthetic data, but no real measurements exist yet.
- Literature cross-referencing, deeper: surface the actual relevant papers per candidate, not just a count.
- A queryable biomarker index across genes and subtypes.

## Related public tools

Cleavr doesn't try to replace general-purpose portals like [cBioPortal](https://www.cbioportal.org/), [GEPIA2](http://gepia2.cancer-pku.cn/), [UALCAN](https://ualcan.path.uab.edu/), the [GTEx Portal](https://gtexportal.org/), or [UCSC Xena](https://xenabrowser.net/). Case Study 01 is a narrow, single-question companion analysis with every number independently verified; the roadmap above is what turns this into something broader.

## License

No license yet: code and data here are visible for reference, but not licensed for reuse. Open to revisiting this as the project matures.

## Contact

Sushanth Somayajula: [sushanthsomayajula@gmail.com](mailto:sushanthsomayajula@gmail.com)
