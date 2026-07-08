# Cleavr

Systematic biomarker discovery for triple-negative breast cancer (TNBC), subtype by subtype.

**Live site:** https://sushanthsomayajula.github.io/Cleavr/

Cleavr tests candidate TNBC biomarkers against public, patient-level sequencing data, broken out by molecular subtype, and cross-checks every claim against the primary literature before publishing it. Every number is independently re-derived from raw source files in a second pass. Null results are reported as directly as positive ones.

This repo currently holds **Case Study 01**: does GNRHR (the LHRH/GnRH receptor — the target of an active TNBC nanoparticle-targeting strategy) express differently across TNBC's molecular subtypes? Result: no statistically significant difference (Kruskal-Wallis p = 0.71), and no relationship to survival (log-rank p = 0.98) in this cohort (n = 116).

## Structure

```
Cleavr/
├── README.md                  <- this file
├── index.html                  <- redirects to website/index.html (GitHub Pages entry point)
├── website/
│   ├── index.html               <- home / project pitch
│   ├── gnrhr.html                <- Case Study 01: full write-up + figures
│   ├── explorer.html             <- interactive tool: pick a subtype, see expression data + a live 3D receptor structure
│   └── README.md                 <- short note about this subfolder
├── code/
│   ├── biomarker_pipeline.py   <- run this to reproduce/extend the analysis for any gene
│   └── tnbc-analysis.ipynb     <- same analysis as a notebook (Jupyter), with narration
├── results/
│   ├── gnrhr_by_subtype.png    <- main figure (GNRHR)
│   ├── gnrhr_survival_km.png   <- survival figure (GNRHR)
│   ├── gnrhr_tnbc_final.csv    <- one row per patient, all variables (GNRHR)
│   ├── gnrhr_analysis_results.json   <- every stat, machine-readable (GNRHR)
│   ├── ghrhr_tnbc_final.csv    <- side investigation into a different, unrelated gene (GHRHR) — see docs
│   └── ar_*                    <- same set of outputs for AR (2nd gene, validates the pipeline works generically)
├── docs/
│   ├── GNRHR_TNBC_Summary.txt  <- plain text write-up
│   ├── GNRHR_TNBC_Summary.md   <- same write-up, Markdown formatted
│   └── UNDERSTANDING_THE_FINDINGS.txt  <- plain-language walkthrough of every concept and number
├── data/                        <- gitignored; small reference files only, not pushed
└── brca_tcga/                   <- gitignored; full raw TCGA download (~1.3 GB), not pushed
```

The large raw TCGA data (`brca_tcga/`, ~1.3 GB) and `data/` are intentionally excluded from this repo (see `.gitignore`) — they're re-downloadable from cBioPortal and don't belong in git. `code/biomarker_pipeline.py` and `code/tnbc-analysis.ipynb` expect `../brca_tcga` to exist locally if you want to reproduce the analysis from scratch:

```
cd code
python3 biomarker_pipeline.py GNRHR AR
```

Pass any gene symbol(s) present in the expression matrix — the cohort and subtype assignment are built once and reused, only the target gene changes.

## Methodology notes (see `docs/` for the full write-up)

- TNBC cohort (n=116) = ER-negative, PR-negative, HER2-negative by IHC. 115 had usable RNA-seq data.
- Subtype (BL1/BL2/M/LAR) is a **marker-gene proxy**, not the validated Lehmann classifier — this TCGA study doesn't ship real subtype calls.
- Kruskal-Wallis across subtypes: p = 0.71 (not significant). Pairwise comparisons remain non-significant even after Bonferroni correction (all p = 1.0).
- Survival analysis (log-rank p = 0.98) is based on only 18 deaths out of 115 patients (15.7% event rate) — underpowered, treat as inconclusive rather than a confident null.
- All numbers above were independently re-derived from the raw TCGA files in a second pass and cross-checked against the saved results — they match exactly.

## Roadmap

**Done:** parameterized pipeline (`code/biomarker_pipeline.py`, any gene, not just GNRHR); systematic candidate sourcing from Open Targets (`code/source_candidates.py`); batch screening with Benjamini–Hochberg FDR correction across all genes tested together (`code/screen_candidates.py`); druggability (ChEMBL) and TNBC-specific literature-gap (PubMed) scoring per candidate (`code/candidate_scoring.py`). See `docs/ROADMAP.md` for the full writeup, including honest caveats on the first 15-gene screen.

**Still open:**
- Cross-cohort validation against a second dataset (METABRIC) so results don't rest on one cohort.
- Literature cross-referencing, deeper — surface the actual relevant papers per candidate, not just a count.
- A queryable biomarker index across genes and subtypes.

## Related public tools

Cleavr doesn't try to replace general-purpose portals like [cBioPortal](https://www.cbioportal.org/), [GEPIA2](http://gepia2.cancer-pku.cn/), [UALCAN](https://ualcan.path.uab.edu/), the [GTEx Portal](https://gtexportal.org/), or [UCSC Xena](https://xenabrowser.net/). Case Study 01 is a narrow, single-question companion analysis with every number independently verified; the roadmap above is what turns this into something broader.

## License

No license yet — code and data here are visible for reference, but not licensed for reuse. Open to revisiting this as the project matures.

## Contact

Sushanth Somayajula — [sushanthsomayajula@gmail.com](mailto:sushanthsomayajula@gmail.com)
