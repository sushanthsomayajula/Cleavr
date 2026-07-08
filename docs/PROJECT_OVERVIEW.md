# Cleavr — Project Overview

A plain-language guide to everything that's been built, why each piece exists, and how it all fits together. Written so you can explain any part of this project to someone else — or to yourself, six months from now — without having to reconstruct the reasoning.

---

## 1. What Cleavr actually is, right now

Cleavr is, honestly, two things layered on top of each other: one real, verified scientific case study, and a professional public container built around it. Both are genuine. Neither should be oversold as more than it is.

The case study: does GNRHR (the LHRH/GnRH receptor) express differently across triple-negative breast cancer's molecular subtypes? We tested this against real patient data and found no statistically significant difference — a legitimate, well-caveated null result.

The container: a live site at `cleavr.bio`, a public GitHub repository, a written roadmap, and a name and brand identity. This is the professional presentation layer — real and useful, but it is not itself evidence that a biomarker-discovery *company* works. That evidence would come from doing this same rigorous process across many genes and multiple datasets, which hasn't happened yet.

---

## 2. The science — Case Study 01 (GNRHR)

### What we asked
TNBC isn't one disease biologically — the Lehmann classification (2011) splits it into molecular subtypes (BL1, BL2, M, MSL, IM, LAR) with different underlying biology. LHRH-conjugated nanoparticles are a real, published targeting strategy for TNBC that works by binding the GNRHR receptor on the tumor cell surface. Nobody had checked whether that receptor is expressed consistently across TNBC's subtypes — if it isn't, a receptor-targeted therapy could work well in one subtype and poorly in another.

### How we answered it
We used TCGA (The Cancer Genome Atlas) breast cancer data, pulled via cBioPortal — a public repository of tumor genomic data used throughout cancer research. We defined a TNBC cohort (116 patients, ER-negative/PR-negative/HER2-negative), approximated molecular subtype using marker genes (since this particular TCGA study doesn't ship official Lehmann subtype calls), and compared GNRHR expression across those subtype groups using standard non-parametric statistics (Kruskal-Wallis, Mann-Whitney), then checked whether expression level predicted survival.

### What we found
No statistically significant difference in GNRHR expression across subtypes (p = 0.71), and no relationship between GNRHR expression and survival (p = 0.98). LAR trended numerically highest, but not significantly so.

### Why you can trust this specific result
This is the part worth understanding, not just believing: every number was independently recomputed in a second pass and matched exactly. We caught and fixed a real error early on — an earlier pass had accidentally measured GNRH1 (the hormone/ligand gene) instead of GNRHR (the receptor gene) — two similarly-named genes with completely different biological roles. We also verified the 3D structure model against the original PDB file's authoritative residue-numbering records (called DBREF records) rather than trusting a secondary description, and caught a numbering error before it shipped that would have mislabeled real receptor as bacterial scaffold protein in the 3D viewer. Full plain-language walkthrough of every statistical concept: `docs/UNDERSTANDING_THE_FINDINGS.txt`.

---

## 3. The engineering — how the site actually works

**It's a static site. No backend, no database, no server-side code.** Every page is a plain HTML file with embedded CSS and JavaScript. This was a deliberate choice, not a limitation: the site doesn't need to compute anything live — all the statistics were computed once, saved as JSON/CSV, and the pages just read and display that saved data. Static sites are free to host, essentially unbreakable, and load instantly. The tradeoff, and it's an important one for the roadmap: a static site can't answer a *new* question on demand (like "show me gene X instead"). That requires a real backend, which is explicitly item 6 on the roadmap, not something built yet.

**Three pages:**
- `index.html` (Home) — the company pitch: problem, approach, current status, roadmap, contact.
- `gnrhr.html` (GNRHR Study) — the full scientific write-up, the one described in section 2 above.
- `explorer.html` (Explorer) — an interactive tool. Pick a subtype, and JavaScript re-renders a chart and stats panel using data that's embedded directly in the page (no server call). The one exception: the 3D receptor structure is fetched live, in your browser, from RCSB (the official protein structure database) — that's the only part of the site that needs an active internet connection beyond just loading the page.

**How `cleavr.bio` actually routes to the site:** Domain names don't point anywhere by themselves — DNS (Domain Name System) is the lookup table that maps a name like `cleavr.bio` to actual server addresses. We added four **A records** (each an IP address) pointing the bare domain at GitHub's servers, and one **CNAME record** (a "this name is really just an alias for that name" record) pointing `www.cleavr.bio` at your GitHub Pages URL. A file named `CNAME` in the repo itself tells GitHub which custom domain to expect and serve. Once all of that matched up, GitHub issued an HTTPS certificate automatically, which is why the site loads securely.

---

## 4. The business and branding side

**Why "Cleavr":** wordplay on "cleavage" (breast) and "cleaver" (cutting/destroying cancer) — a real, defensible piece of brand reasoning, not just a random name.

**What we found before committing:** a real, currently-operating company already uses the identical name (`cleavr.io`, a server-hosting platform) — different industry, but it means Cleavr won't own its own name in search results, and there's enough thematic overlap with the roadmap's eventual "software platform" direction to be worth real trademark diligence later, not just a Google search. You weighed that honestly and chose to keep the name anyway. That's a legitimate call — plenty of companies share a name with an unrelated small company — but it's a known, accepted tradeoff now, not an unknown risk.

**The domain:** `cleavr.bio` was purchased on Namecheap for the (promotional) price of $6.18/year — worth checking the renewal price before next year, since intro pricing on new TLDs is common and often much higher on renewal.

**On calling yourself "Founder":** legitimate to use — it describes what you're doing, not a legal status, and doesn't require incorporation. The thing that actually protects your credibility is keeping whatever's *next to* the title honest about stage ("early-stage," pointing to what's actually built) rather than the title itself.

---

## 5. The habit underneath all of this: verify before you publish

This is worth naming explicitly, because it's the actual transferable skill here, more than any single finding. Multiple times in this project, a claim got checked against a primary source instead of taken on faith, and the checking changed the outcome:

- **GNRH1 vs. GNRHR** — caught and corrected before any numbers were finalized.
- **PDB residue numbering** — the 3D structure's fusion-domain coloring was verified against the file's own authoritative records rather than a secondary description, catching an error that would have mislabeled real receptor as bacterial scaffold.
- **The GHRH-R claim** — when a "quick Google search" suggested GHRH-R (a completely different gene from GNRHR) was highest in the LAR subtype, we tested it directly against our own data rather than accepting the search result, and found the opposite of what was claimed.
- **The AI-search-summary citation mismatch** — a pasted AI summary attached a real citation to a claim the actual paper didn't support. Caught by reading what the paper actually studied, not just trusting the citation was doing what it claimed to do.
- **Dr. Hu's citation, today** — pulled her actual papers from PubMed rather than reconstructing citations from memory, and corrected a misheard "State Department" into her real, verified NC State University affiliation before it went anywhere public.

If Cleavr becomes a real biomarker-discovery tool, this — the discipline of not publishing a claim until it's been checked against a primary source — is the actual product, more than any individual result.

---

## 6. What's real vs. what's still a gap

**Real, done, and solid:**
- One fully verified scientific case study
- A professional, publicly accessible site on a real domain
- Clean GitHub repo with real documentation
- A named, specific technical roadmap

**Not yet real — the honest gap:**
- The analysis pipeline still only runs on one gene (GNRHR); it isn't generalized yet
- Everything rests on one dataset (TCGA); nothing has been cross-validated against a second cohort
- There's no automated literature-checking — every verification so far has been done by hand, which doesn't scale past one gene
- No multiple-testing correction framework exists yet for screening many genes at once (this becomes essential, not optional, the moment you test gene #2 and beyond — see `docs/ROADMAP.md` for why)

None of that is a criticism — it's just where the project actually is, stated plainly so the next session can pick up exactly where it left off instead of guessing.

---

## 7. Where everything lives

```
Cleavr/
├── README.md                   ← technical/repo overview
├── docs/
│   ├── PROJECT_OVERVIEW.md      ← this file
│   ├── ROADMAP.md                ← the actual technical next-steps plan
│   ├── UNDERSTANDING_THE_FINDINGS.txt  ← every statistical concept, explained plainly
│   └── GNRHR_TNBC_Summary.md/.txt      ← condensed one-page write-up
├── website/                     ← the live site's source (index/gnrhr/explorer .html)
├── code/                        ← the actual analysis (pipeline, scoring, screening, cross-cohort validation)
├── results/                     ← every number and figure the site displays, one folder per gene tested
└── brca_tcga/, brca_metabric/   ← raw cohort data (not in the public repo — too large, re-fetchable)
```

## 8. A short glossary, for anything above that wasn't obvious

- **TCGA** — The Cancer Genome Atlas, a large public database of tumor genomic data.
- **cBioPortal** — the web tool most researchers use to actually pull TCGA data.
- **RSEM** — the algorithm/unit used to measure how much a gene is being expressed.
- **PDB** — Protein Data Bank, the public archive of solved 3D protein structures.
- **DNS, A record, CNAME** — the internet's address book; see section 3 above.
- **GitHub Pages** — free static-site hosting built into GitHub, used here to serve the actual site.
- **p-value / Kruskal-Wallis / Mann-Whitney / log-rank** — statistical tests; fully explained in `UNDERSTANDING_THE_FINDINGS.txt`.
- **FDR (false discovery rate) correction** — the statistical safeguard needed once you're testing many genes at once instead of one, to avoid finding "significant" results that are really just chance. Not yet built — see the roadmap.

---

*Next session: pick up at `docs/ROADMAP.md`, item 1 — refactoring the pipeline to accept any gene, not just GNRHR.*
