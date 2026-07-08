# Cleavr — technical roadmap

Working notes on what turns Case Study 01 (one gene, one dataset) into an actual biomarker-discovery pipeline. Ordered roughly by what unlocks the next step.

## 1. Parameterize the pipeline

`code/run_gnrhr_analysis.py` currently hardcodes GNRHR. The first real engineering step is refactoring it into something that takes a gene symbol (or list of symbols) as input and produces the same outputs (subtype comparison, survival split, figures, JSON stats) for any gene, not just this one. Concretely: pull the TNBC cohort construction, marker-gene subtype assignment, and survival-split logic into a reusable module; keep only "which gene(s) to test" as the variable input.

## 2. Cross-cohort validation

Right now every finding rests on one TCGA study (Firehose Legacy BRCA). The next credibility step is running the same pipeline against a second, independent cohort — METABRIC is the natural choice: also public via cBioPortal, large TNBC-relevant sample size, different sequencing platform. A finding that replicates across both TCGA and METABRIC is much stronger evidence than either alone. Needs: harmonizing clinical field names (survival status/time, receptor status fields differ between studies), and deciding how to handle the fact that METABRIC uses microarray, not RNA-seq (different expression units — z-score or rank-based comparison, not raw value comparison, across cohorts).

## 3. A real candidate gene list

GNRHR was chosen because of a specific external motivation (Dr. Hu's nanoparticle work), not systematically. A real biomarker-discovery tool needs a principled way to generate candidates instead of testing genes one at a time as they come up in conversation. Options worth evaluating: the druggable genome (ChEMBL target list), Open Targets' TNBC-associated target scores, DepMap dependency data for TNBC cell lines, or a literature-mining pass over recent TNBC targeted-therapy papers to extract candidate receptor/target genes automatically.

## 4. Multiple-testing discipline at scale

Testing one gene, a Kruskal-Wallis + six pairwise comparisons with Bonferroni correction is enough. Testing dozens or hundreds of candidate genes changes the statistics problem entirely — screening many genes for "significant" subtype differences will produce false positives by chance alone. This needs a genome-wide-style multiple-testing correction (Benjamini-Hochberg FDR, not just Bonferroni within one gene's tests) built into the pipeline from the start, not bolted on after the fact. This is the single easiest way to accidentally p-hack a "discovery" pipeline — worth being paranoid about before showing results to anyone.

## 5. Literature cross-referencing, automated

Every claim in Case Study 01 was manually checked against PubMed/primary sources — that doesn't scale past one gene. Worth building a lightweight step that, given a gene symbol, pulls recent relevant literature (PubMed/Europe PMC API) and surfaces it alongside the statistical result, so a claim like "expression differs by subtype" is never presented without what's already published about that gene in TNBC.

## 6. From static site to real backend

The current site is precomputed static HTML — fine for one gene, not for "look up any gene on demand." Once there's a working parameterized pipeline (step 1) and it's fast enough, the natural next step is a small backend (FastAPI or similar) that runs the pipeline on request and a frontend that queries it, instead of hand-building a page per gene. Doesn't need to happen before there's a second or third case study worth serving dynamically.

## Immediate next steps (rough order)

1. Refactor `run_gnrhr_analysis.py` into a `gene -> results` function.
2. Pick one more gene (ideally one with a clear TNBC subtype-specific hypothesis in the literature) and run it through the refactored pipeline as a real test of the parameterization, not just GNRHR again.
3. Pull METABRIC and get the two cohorts' clinical fields harmonized enough to run the same pipeline on both.
4. Only after 2–3 genes exist: build the FDR-corrected multi-gene comparison view.
