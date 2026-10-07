# Spec: delegate to validated R packages

**Status:** phases 0-2 and R1 implemented · **Date:** 2026-10-06 · **Scope:** `CooccurrenceAnalysis` (R)

## 1. Why

This package is a port of the Julia `CooccurrenceAnalysis`, and the port copied
Julia's habit of hand-writing algorithms that Julia's ecosystem lacks. R has
validated, published implementations of several of them. Hand-written code we
don't need is code we have to test and keep correct ourselves. A review on
2026-10-06 (see `tests/testthat/test-r-ecosystem-crosscheck.R`) found:

- The 2×2 statistics already delegate to base R (`stats::phyper`, `stats::p.adjust`)
  and agree with `cooccur` (Veech 2013) and `stats::fisher.test` to rounding.
- `modularity_q()` agrees with `igraph::modularity()` to 1e-10 but is hand-written.
- `.louvain()` is **phase 1 only** (local moving, no aggregation), mirroring
  Julia's `network_metrics.jl:180`. It is not full Louvain.
- Record/item incidence is built twice by hand (`.record_counts()` in
  `network_construction.R` and `record_item_matrix()` in `hdp_wrappers.R`),
  and the HDP carries its own sparse-matrix helpers (`.sparse_binary`,
  `.sx_times`, `.sx_crossprod`).
- Not-yet-ported features (association rules, chi-square, Bonferroni/Holm,
  label propagation, flat-K mixtures) all have mature R implementations.
- No R package implements the HDP Bernoulli mixture or the group-stratified
  network with exclusive items. Those stay as they are.

## 2. Governing rule

The Julia package remains the reference (CLAUDE.md). This refactor adds one rule
for when to delegate:

> **Delegate when the R function computes the same quantity under this package's
> conventions.** Keep the exported signature and return shape, and keep this
> package's input validation in front of the call (NA ids, unknown keys, the BH
> family, Haldane correction, C-locale ordering, `round_digits()`). The existing
> Julia-derived tests must pass **unchanged**.
>
> **When the R function computes a different quantity, even a better one,
> delegating is a divergence, not a refactor.** It needs a decision (§6), a new
> option or a documented change of default, and an entry in a "Divergences from
> the Julia reference" section of CLAUDE.md and README.

Every delegation also gets a cross-check test against the package it delegates
to, alongside the Julia-derived test.

## 3. Inventory and actions

| ID | Component | Now | R package | Action | Julia parity |
|----|-----------|-----|-----------|--------|--------------|
| R1 | Community detection | Hand-written phase-1 local moving | `igraph::cluster_louvain` / `cluster_leiden` | Replace (**D1 decided**: treat the phase-1 limit as a bug) | **Diverges, deliberately** |
| R2 | `modularity_q()` | Hand-written | `igraph::modularity` | Delegate | Same (verified) |
| R3 | Incidence and pair counts | Two hand-built copies with tabulate loops | `Matrix::sparseMatrix` + `Matrix::crossprod` | One shared internal helper | Same (integers) |
| R4 | HDP sparse helpers | `.sparse_binary`, `.sx_times`, `.sx_crossprod` | `Matrix` | Replace | Same up to floating-point summation order |
| R5 | Fisher test, BH | `stats::phyper`, `stats::p.adjust` | — | **No change** | — |
| R6 | Chi-square; Bonferroni/Holm | Raise "not ported" | `stats::chisq.test`, `stats::p.adjust` | Port as thin wrappers | Same (verify) |
| R7 | Label propagation | Raises | `igraph::cluster_label_prop` | Port only if **D3** says yes | Both are nondeterministic |
| R8 | Association rules | Not ported | `arules::apriori` / `eclat` | Port as a wrapper | Same counts (verify) |
| R9 | Flat-K Bernoulli mixture | Not ported | `flexmix` (`FLXMCmvbinary`), `poLCA` | Per **D2** | Diverges unless hand-ported |
| R10 | Centrality, layout | Already `igraph` | — | **No change** | — |
| R11 | Network plots | Hand-built ggplot | `ggraph` | **No change** (see §7) | — |
| R12 | HDP | Hand-written CAVI | none found | **No change**; add flat-K sanity test | — |

### R1 Community detection

`.louvain()` runs local moving to convergence and stops. It never aggregates
communities into super-nodes, so it can stall at a partition that full Louvain
or Leiden would improve. On the planted-partition test graph it matches igraph
to within 0.02 Q; on larger or less clean graphs the gap can grow.

Verified facts that shape the design:

- `igraph::cluster_louvain` is **stochastic**: on the Zachary graph, two
  different seeds gave different memberships (igraph 2.3.4). A wrapper needs a
  `seed` argument with no default, or the caller sets the seed explicitly. The
  wrapper should use `withr`-style local seeding, or save and restore
  `.Random.seed`, so it doesn't change the caller's RNG stream.
- igraph exposes no gain-tie count. Return `n_gain_ties = NA_integer_` for the
  igraph methods, and say so in the documentation.

**Decision (D1):** the phase-1-only implementation is a bug in the reference,
not a specification. Julia's own comment ("sufficient for small graphs") makes
it a shortcut, not a methodological choice, and nothing in either package
depends on stopping before aggregation. The port replaces it with the full
algorithm and records the change as a divergence.

Changes:

- `detect_communities(net, method = "leiden", seed, resolution = 1)` calls
  `igraph::cluster_leiden(g, objective_function = "modularity",
  weights = E(g)$weight, resolution = resolution, n_iterations = -1)`.
  igraph's Leiden defaults to the CPM objective and two iterations, so both
  arguments must be set. Leiden is the default because it guarantees
  connected communities, which Louvain doesn't. `method = "louvain"` calls
  `igraph::cluster_louvain(g, weights = E(g)$weight, resolution = resolution)`.
  Any other method still raises. Because the default changes the meaning of
  the old `"louvain"` default, any call that relied on it implicitly now gets
  Leiden; the divergence entry has to say so.
- **`seed` has no default.** Results depend on the seed (verified above), and
  under the package's convention a value that changes results and that the
  package can't choose is an argument without a guessable default. Seed the
  call locally and restore the caller's `.Random.seed` afterward.
- Delete `.louvain()`, `.modularity_gain()`, `.neighbors()` (if nothing else
  uses it), the `max_iter` argument, and the `n_gain_ties` return field. That
  field existed only to trace mismatches against Julia's tie-breaking, and
  that comparison no longer applies.
- The rest of the return shape is unchanged: assignments renumbered `1..K` by
  first appearance in `net$items` order (igraph's raw labels are arbitrary), `communities` sorted in C-locale order, and `modularity`
  computed by `modularity_q()` (not taken from igraph) so it is reported the
  same way for every method.
- `run_network_pipeline()` and `stratified_network_analysis()` take `seed` and
  pass it, along with `method`, through. Stratified analysis uses the same seed
  for every group, so a group's partition doesn't depend on how many groups
  come before it.

Tests:

- Expected partitions can no longer come from running Julia. Rewrite the
  partition assertions in `test-network.R`, `test-network-stratification.R`,
  and `test-validation.R` (52 references to communities or pipelines in
  total; each needs review, though most won't change) against fixtures whose
  community structure is unambiguous by construction: disconnected cliques,
  or planted blocks with no between-block edges. Then the expected partition
  can be worked out on paper and doesn't depend on the seed.
- Add a regression test for the bug: a graph where local moving alone stalls
  below the optimum and aggregation improves it (for example, a ring of
  cliques), asserting a Q above the phase-1 result.
- Keep Julia-derived expectations for everything upstream of partitioning
  (edges, weights, q-values), since those don't change.
- Replace the CLAUDE.md parity rule on Louvain tie-breaking with the
  divergence entry, and open an issue (or fix) against the Julia package so
  the reference converges on full Louvain too.

### R2 `modularity_q()`

Reimplement as `igraph::modularity(g, membership, weights = E(g)$weight)` on the
graph built from `adjacency`. **Keep the special case**: igraph returns `NaN` for
an edgeless graph (verified), while this package returns `0`. Handle the
empty and edgeless cases before calling igraph. `.modularity_gain()` stays only
if `"local_moving"` stays.

### R3 One incidence helper

Add an internal `.incidence(event_df, min_items, group_col = NULL)` that returns:

- `x`: a sparse logical records × items `Matrix` (`lgCMatrix`), deduplicated,
  with columns in C-locale item order and rows in first-appearance id order
  (as now);
- `ids`, `group` (if requested), and `n_item`.

`.record_counts()` becomes `pair_counts <- as.matrix(Matrix::crossprod(x))` with
the diagonal used for `n_item`. `record_item_matrix()` calls the same helper and
returns a dense logical matrix as it does now, so its exported contract does
not change. Validation (missing columns, NA ids and groups) moves into the
helper so it is written once.

### R4 HDP sparse helpers

Replace `.sparse_binary`, `.sx_times`, and `.sx_crossprod` with `Matrix`
products. Summation order may change, so check `test-hdp-known-answers.R` at
its current tolerances. If a known answer moves outside tolerance, **stop and
report**; do not loosen the tolerance to make it pass.

### R6 Chi-square, Bonferroni, Holm

- `test = "chisq"` → `stats::chisq.test(ct, correct = FALSE)$p.value`. Julia's
  `HypothesisTests.ChisqTest` on a 2×2 table is Pearson without continuity
  correction; **verify by running the Julia reference** before writing expected
  values. Decide what happens with expected cell counts below 5: R warns, and
  the no-silent-failures rule bans `suppressWarnings`. Either raise, or
  record the condition in a column, matching whatever Julia does.
- `correction = "bonferroni" | "holm"` → `stats::p.adjust`. The family is still
  exactly the pairs that co-occur.
- Each option pairs with its Julia testset (CLAUDE.md) before its raise is removed.

### R8 Association rules

`mine_frequent_itemsets()` and `mine_association_rules()` wrap `arules`
(`apriori`, or `eclat` for itemsets). Map the Julia arguments: `min_support` →
`support`; `min_count` → `support = min_count / N`; `max_length` → `maxlen`.
Return the Julia columns (`Itemset`/`LHS`, `RHS`, `Support`, `Confidence`,
`Coverage`, `Lift`, `N`, `Length`) as a plain data.frame, sorted in C-locale
order. To verify against the Julia reference (which uses RuleMiner):
whether rules with an empty LHS are emitted (`arules` emits `{} => X` unless
`minlen = 2`), the counting at the support boundary (≥ vs >), and how ties are
ordered. `arules` goes in `Suggests`, guarded by `requireNamespace()` with an
error that names the package.

### R9 Flat-K Bernoulli mixture

Julia's `fit_bernoulli_mixture` is EM with multiple restarts, empirical-Bayes
Beta priors (MAP rather than ML), and BIC-based `select_K`. `flexmix` and
`poLCA` fit maximum-likelihood EM, which is a different estimator. See **D2**.
The Turing MCMC and ADVI variants (`bayesian_turing.jl`) are out of scope for
this refactor.

### R12 HDP sanity test

Add one cross-check: on data simulated from K = 3 Bernoulli components with one
group, an HDP fit with small `gamma` should recover the components as well as
`flexmix` with K = 3 does. Use the adjusted Rand index against the truth (above
a threshold set in review) and matched item probabilities within a tolerance
set in review. `flexmix` goes in `Suggests`, and the test uses
`skip_if_not_installed`.

## 4. Dependencies

- **Imports:** add `Matrix` (a recommended package that ships with R).
- **Suggests:** `arules`, `flexmix`, `cooccur` (already added; GitHub-only for
  current R, so either add `Remotes: griffithdan/cooccur` or rely on the skip).
- NAMESPACE is maintained by hand: add `importFrom` or use `pkg::` calls, and
  export every new public function.

## 5. Phases and acceptance

Each phase is its own branch and merge. "All tests pass" means
`testthat::test_local(".")` with zero failures and zero warnings, and
`lintr::lint_package()` clean.

1. **Phase 0, cross-check baseline.** Commit `test-r-ecosystem-crosscheck.R`
   and the `Suggests` change. *Done when:* committed, and the suite passes.
2. **Phase 1, internal delegation, no behavior change (R2, R3, R4).** *Done
   when:* every existing test passes **unmodified**; the cross-checks still
   pass; `git diff tests/` shows only new tests.
3. **Phase 2, new options via wrappers (R6, R8, and R7 if D3 says yes).** *Done
   when:* each option has a testthat test whose expected values come from
   running the Julia reference, plus a cross-check against the wrapped R
   function; the "not ported" raises are removed only for options that pass
   both.
4. **Phase 3, divergences (R1, R9, R12).** R1 is decided and can start once
   Phase 1 merges. *Done when:* the local-moving code is gone; partition tests
   use structure that is unambiguous by construction, plus the
   aggregation-regression test; CLAUDE.md and README carry the Divergences
   section; and D2 is decided.
   *R1 done 2026-10-06. R12 done 2026-10-07 (`test-hdp-flexmix.R`). R9 closed by D2 = skip.*

## 6. Decisions needed

- **D1: Community method. Decided 2026-10-06:** Julia's phase-1-only Louvain
  is a bug; replace it with igraph. The default is `"leiden"` (modularity
  objective), and `"louvain"` is available as an option (§R1). Before merging, check whether
  any consuming analysis has reported partitions from the old method, since
  those numbers will change.
- **D2: Flat-K mixture.** (a) Hand-port Julia's MAP-EM for parity; (b) wrap
  `flexmix` and document that it is the ML estimator; (c) don't port it, since
  the HDP subsumes it for the analyses that use this package.
  **Recommendation: (c)** unless a consumer needs flat-K output, then (a).
  **Decided 2026-10-07: (c).** R9 stays a raise. R12's flexmix cross-check found the
  HDP and a flat K = 3 ML fit agree on well-separated data (ARI 0.975-0.995, identical).
  Not tested: overlapping components or sparse items, where MAP and ML could differ.
- **D3: Label propagation.** Port it as an `igraph::cluster_label_prop`
  wrapper with a required seed, or keep the raise. CLAUDE.md records that the
  analysis never uses it. **Recommendation: keep the raise.**

## 7. Non-goals

- **No runtime dependency on `cooccur`.** It is GitHub-only for current R and
  reports two-sided results; it stays a test-time cross-check.
- **No `ggraph`.** The plot-encoding tests read layers back through
  `ggplot2::ggplot_build()`. `ggraph` would change the layer structure and add
  a dependency without changing what is drawn.
- **No change** to the HDP algorithm, group stratification, the exclusive-items
  contract, rounding, or sort order.
