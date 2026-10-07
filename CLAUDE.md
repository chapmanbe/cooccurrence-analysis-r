# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

`CooccurrenceAnalysis` (R) is a domain-neutral package for co-occurrence patterns in binary (presence/absence) event data. It is an R port of the Julia package of the same name (`../cooccurrence-analysis`), which remains the **reference implementation**: every ported function is checked against it, and where the two disagree the Julia behavior is the specification unless a divergence is documented. The data model is neutral: **records** (`id`) accumulate **events**, each carrying one **item**, and records may belong to **groups**. Nothing in this package names an application domain, and nothing should.

## What is ported

Pairwise association statistics (Fisher or chi-square, with BH, Bonferroni or Holm), the significance-filtered co-occurrence network, Leiden (default) or Louvain communities via igraph, per-node metrics, group-stratified networks, the HDP Bernoulli mixture, frequent itemsets and association rules (enumerated by `arules`), and their plots. **Not yet ported:** the flat-K Bernoulli mixtures (including the Turing MCMC/ADVI fits), timing filters, `label_propagation`, and `compare_networks`. Those raise an error naming the Julia function rather than silently doing something else; keep it that way. A new port adds the function, pairs every Julia testset for it with a testthat test, and only then removes the raise. The refactor plan that delegates to validated R packages is `docs/specs/2026-10-06-r-ecosystem-refactor.md`.

## Commands

```bash
# Dependencies (no renv here; a consuming project pins versions in its own lockfile)
Rscript -e 'install.packages(c("igraph", "ggplot2", "testthat", "pkgload", "lintr", "arules"))'
# Optional cross-check dependency (GitHub-only for current R); its test skips without it
Rscript -e 'remotes::install_github("griffithdan/cooccur")'

# Tests; the suite builds its own synthetic fixtures, so no data file is needed
Rscript -e 'testthat::test_local(".")'

# Lint, with the rules in .lintr
Rscript -e 'lintr::lint_package()'
```

Develop with `pkgload::load_all(".")`; there is no build step. `NAMESPACE` is **maintained by hand**, not by roxygen2: add every new export there yourself, or callers fail with "could not find function". There is no `man/`; the roxygen comments are documentation in the source.

## Data contract

Functions take an **event table**: a data.frame with one row per (record, item) event and columns `id` and `item` (plus a group column for stratified analysis). A record holding an item twice counts once; every count is a count of records. Records with fewer than `min_items` distinct items (default 2) cannot express co-occurrence and are excluded.

## Architecture

One file per method family in `R/`: `network_construction.R` (2x2 statistics, pairwise tests, network building), `network_metrics.R` (modularity, communities via igraph, centrality), `network_group_stratification.R` (one network per group), `network_visualization.R` (ggplot2), `mining.R` (frequent itemsets and rules via `arules`), `pipeline.R` (`run_network_pipeline`), `utils.R`. Statistics return plain data.frames and lists; plots return ggplot objects and never save.

## Parity rules (each was a real bug)

- **Rounding** goes through `round_digits()`, never `round(x, digits)`. R's decimal rounding disagrees with Julia's on ties (`round(0.15, 1)`), which changes published tables.
- **Fisher tests are one-sided** (upper tail), and the **BH family is exactly the pairs that co-occur** in at least one record, not every possible pair. Widening the family changes every q-value.
- **`odds_ratio` applies the Haldane correction by default** (+0.5 to all four cells when any cell is zero), as Julia does; the raw ratio is `Inf`.
- **Sort in C-locale byte order** (`.sort_c`); tests pin `LC_COLLATE = "C"`. Group results are named and ordered by group value, never by position.

## Divergences from the Julia reference

Each is a Julia bug the port fixes deliberately; the paired test is marked `R-only`.

- **Chi-square on a zero margin raises.** Julia returns `NaN` for a pair whose item is in every record, and the `NaN` then turns every q-value in the family into `NaN`.
- **Communities use full Leiden (default) or Louvain from igraph.** Julia runs only Louvain's first phase (local moving, no aggregation), which stalls below the optimum: on a ring of 20 triangles it stops at Q = 0.70 where aggregation reaches 0.77. Partitions therefore differ from Julia's, and `community_method`'s default is now `"leiden"`. Both algorithms are randomized, so `detect_communities()`, `stratified_network_analysis()` and `run_network_pipeline()` require `seed` (a whole number, applied locally, or `NULL` for the caller's stream). Expected partitions in tests come from fixtures unambiguous by construction, not from Julia.
- **`min_count` in mining is used as a count.** Julia converts it to a fraction and RuleMiner back with `ceil(min_support * n)`, which can land one above it (7 of 100 gives 8).

## Conventions

- **Anything the caller knows and the package cannot is an argument without a guessable default.** `exclusive_items` (items possible in only one group) has no default; `NULL` is the explicit "none", its keys must be group values, and an item may belong to one group only. A list that matches none of the data's items raises: a mismatch there silently disables the exclusion. `seed` likewise has no default: a randomized result needs the caller to choose reproducibility.
- **Labels are arguments, not vocabulary.** Plot titles and axis labels are parameters with neutral defaults (`"Item"`, `"items"`); callers pass their own domain nouns. Never hardcode a domain term.
- **Plots match the Julia figures' encoding**: gray70 edges at alpha 0.6, nodes colored by community in Makie's seven `wong_colors` (no black), recycled by community label. What a figure *records* (panel order, colors) is read back from `ggplot2::ggplot_build()`, not rebuilt from the arguments, so a test of the record tests the drawing.
- **No silent failures.** Tests run with partial-match warnings on (`tests/testthat/setup.R`); use `vapply`, not `sapply`; no `suppressWarnings` or `try`. Raise on NA ids, NA groups, unknown keys, and duplicate items rather than dropping them.
- Every change to numerical behavior needs a test whose expected values come from **running the Julia reference**, not from reasoning about what it should return.
