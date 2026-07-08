# Cleavr: technical roadmap

Working notes on what turns Case Study 01 (one gene, one dataset) into an actual biomarker-discovery pipeline. Ordered roughly by what unlocks the next step.

## 1. Parameterize the pipeline

**Done.** `code/biomarker_pipeline.py` (formerly `run_gnrhr_analysis.py`) now splits into `build_tnbc_cohort()` (cohort construction + marker-gene subtype assignment, gene-independent, runs once) and `analyze_gene(gene_symbol, cohort_df)` (subtype comparison, figure, survival split, JSON/CSV output, runs per gene). Validated by running AR (LAR-subtype marker) alongside GNRHR through the same pipeline; GNRHR numbers matched the original script exactly.

## 2. Cross-cohort validation

**Done.** `code/fetch_metabric.py` pulls METABRIC clinical + expression data directly from the public cBioPortal REST API (320 TNBC patients by ER-/PR-/HER2-, all with expression data, nearly 3x the TCGA TNBC cohort of 115). `code/cross_cohort_validation.py` independently re-derives Lehmann subtype within METABRIC's own expression distribution (same marker-gene proxy, not reused from TCGA) and re-tests the 8 genes that were significant after FDR correction in the TCGA screen. Result: 7 of 8 replicated (significant in both cohorts, on two different platforms: RNA-seq V2 RSEM vs. Illumina microarray). RRM2 did not replicate (p=0.23 in METABRIC) despite being significant in TCGA, and is reported as such rather than dropped. FTO, the standout candidate from the original screen, replicated (METABRIC p=0.0007), and the direction of its subtype pattern (BL1 lowest) held in both cohorts, though the exact subtype ranking wasn't identical between cohorts. Output: `results/cross_cohort_results.json`. Now published on the site as Case Study 02, alongside GNRHR, via `website/studies.html`.

## 3. A real candidate gene list

GNRHR was chosen because of a specific external motivation (Dr. Hu's nanoparticle work), not systematically. A real biomarker-discovery tool needs a principled way to generate candidates instead of testing genes one at a time as they come up in conversation. Options worth evaluating: the druggable genome (ChEMBL target list), Open Targets' TNBC-associated target scores, DepMap dependency data for TNBC cell lines, or a literature-mining pass over recent TNBC targeted-therapy papers to extract candidate receptor/target genes automatically.

**Partial progress:** `code/candidate_scoring.py` now takes any gene already run through `biomarker_pipeline.py` and pulls two more axes live from public APIs: druggability (ChEMBL: does a known ligand/drug exist, how potent, what mechanism) and literature gap (PubMed: how many TNBC-specific papers already exist on this gene, independent of how well-studied the gene is generally). It deliberately reports these three axes (statistical signal, druggability, literature gap) separately rather than collapsing them into one composite score, since a single number would overstate the precision of the underlying evidence. Ran on GNRHR/AR/CXCR4 as a proof of concept; output in `results/candidate_scorecard.json`. Still ad hoc: genes are still hand-picked one at a time, not generated from a systematic candidate list. That's the next piece.

## 4. Multiple-testing discipline at scale

Testing one gene, a Kruskal-Wallis + six pairwise comparisons with Bonferroni correction is enough. Testing dozens or hundreds of candidate genes changes the statistics problem entirely: screening many genes for "significant" subtype differences will produce false positives by chance alone. This needs a genome-wide-style multiple-testing correction (Benjamini-Hochberg FDR, not just Bonferroni within one gene's tests) built into the pipeline from the start, not bolted on after the fact. This is the single easiest way to accidentally p-hack a "discovery" pipeline, worth being paranoid about before showing results to anyone.

## 5. Literature cross-referencing, automated

Every claim in Case Study 01 was manually checked against PubMed/primary sources; that doesn't scale past one gene. Worth building a lightweight step that, given a gene symbol, pulls recent relevant literature (PubMed/Europe PMC API) and surfaces it alongside the statistical result, so a claim like "expression differs by subtype" is never presented without what's already published about that gene in TNBC.

## 6. From static site to real backend

The current site is precomputed static HTML, fine for one gene, not for "look up any gene on demand." Once there's a working parameterized pipeline (step 1) and it's fast enough, the natural next step is a small backend (FastAPI or similar) that runs the pipeline on request and a frontend that queries it, instead of hand-building a page per gene. Doesn't need to happen before there's a second or third case study worth serving dynamically.

## Immediate next steps (rough order)

1. ~~Refactor `run_gnrhr_analysis.py` into a `gene -> results` function.~~ Done: see `code/biomarker_pipeline.py`.
2. ~~Pick one more gene and run it through the refactored pipeline.~~ Done: AR, see `results/ar/ar_analysis_results.json`.
3. ~~Pull METABRIC and get the two cohorts' clinical fields harmonized enough to run the same pipeline on both.~~ Done: see `code/fetch_metabric.py` and `code/cross_cohort_validation.py`.
4. ~~Once a few more genes exist: build the FDR-corrected multi-gene comparison view.~~ Done: see `code/screen_candidates.py`.
5. ~~Decide what from the cross-cohort result is ready to publish on the site.~~ Done: Case Study 02 published, see `website/studies.html`. FTO replicating across cohort and platform is the strongest single finding so far, but "strongest so far" isn't the same as "validated"; still no wet-lab confirmation.
6. Next: wet-lab flow cytometry validation of FTO at the protein level.
