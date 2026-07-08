# Cleavr

Analysis of GNRHR (LHRH receptor) expression across TNBC molecular subtypes, using TCGA data. Public-facing write-up and interactive explorer live in `website/`.

## Folder structure

```
tnbc-project/
├── README.md                  <- this file
├── website/
│   ├── index.html              <- open this first. Full write-up + both figures.
│   ├── explorer.html           <- interactive tool: pick a subtype, see expression data + a live 3D receptor structure
│   └── README.md                <- short note about this subfolder
├── code/
│   ├── run_gnrhr_analysis.py  <- run this to reproduce everything from scratch
│   └── tnbc-analysis.ipynb    <- same analysis as a notebook (Jupyter), with narration
├── results/
│   ├── gnrhr_by_subtype.png   <- main figure
│   ├── gnrhr_survival_km.png  <- survival figure
│   ├── gnrhr_tnbc_final.csv   <- one row per patient, all variables
│   ├── ghrhr_tnbc_final.csv   <- side investigation into a different, unrelated gene (GHRHR) — see docs
│   └── analysis_results.json  <- every stat, machine-readable
├── docs/
│   ├── GNRHR_TNBC_Summary.txt  <- plain text write-up
│   ├── GNRHR_TNBC_Summary.md   <- same write-up, Markdown formatted
│   └── UNDERSTANDING_THE_FINDINGS.txt  <- plain-language walkthrough of every concept and number
├── data/
│   └── raw/                    <- see note below, mostly not duplicated here on purpose
└── brca_tcga/                  <- full raw TCGA download (~1.3 GB)
```

## Where to start

Open `website/index.html` in a browser for the full write-up, or `website/explorer.html` for the interactive tool. Both are static files — double-click to open, no server needed. They link to each other and to the files in `results/`, `code/`, and `docs/` using paths one level up (`../`), since they live inside `website/`.

## Where the raw data actually lives

The full raw TCGA download (`brca_tcga/`, ~1.3 GB) stays at the project root — it was **not** duplicated into `data/raw/` because copying 1.3 GB of raw data around a project is bad practice, not good organization. `data/raw/brca_tcga/` only holds the small reference files (clinical metadata, case lists, license) for quick lookup.

`code/run_gnrhr_analysis.py` and `code/tnbc-analysis.ipynb` both point at `../brca_tcga` and save their outputs to `../results/`, so run them from inside `code/`:

```
cd code
python3 run_gnrhr_analysis.py
```

## Methodology notes (see docs/ for full write-up)

- TNBC cohort (n=116) = ER-negative, PR-negative, HER2-negative by IHC. 115 had usable RNA-seq data (one patient has no primary-tumor RNA-seq sample at all in TCGA — not a bug).
- Subtype (BL1/BL2/M/LAR) is a **marker-gene proxy**, not the validated Lehmann classifier — TCGA's Firehose Legacy study doesn't ship real subtype calls.
- Kruskal-Wallis across subtypes: p = 0.71 (not significant). Pairwise comparisons remain non-significant even after Bonferroni correction (all p = 1.0).
- Survival analysis (log-rank p = 0.98) is based on only 18 deaths out of 115 patients (15.7% event rate) — underpowered, treat as inconclusive rather than a confident null.
- All numbers above were independently re-derived from the raw TCGA files and cross-checked against the saved results — they match exactly.

## Related public tools

This project doesn't try to replace general-purpose portals like [cBioPortal](https://www.cbioportal.org/), [GEPIA2](http://gepia2.cancer-pku.cn/), [UALCAN](https://ualcan.path.uab.edu/), the [GTEx Portal](https://gtexportal.org/), or [UCSC Xena](https://xenabrowser.net/) — it's a narrower, purpose-built companion analysis answering one specific question (GNRHR expression across TNBC subtypes, tied to a specific nanoparticle-targeting mechanism) with every number independently verified.

## Cleanup note

Files that were superseded loose copies at the project root (from earlier passes of this analysis) may still exist alongside these organized folders — anything starting with `._` (macOS metadata) or `.ipynb_checkpoints/` is safe to delete by hand. Claude cannot delete files itself.
