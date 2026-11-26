# Code review of f849c40..HEAD (2026-11-25)

## Why
The last reviews were `6a09b00` and `f849c40`, both on 2026-09-25. Since then, 18 commits
(about 4,800 lines across 35 files) have had no review: group stratification and
`exclusive_items`, the HDP mixture and multi-seed fitting, chi-square and Bonferroni/Holm,
`arules` mining, igraph Leiden/Louvain, and the seed and random-number handling.

## Plan: two independent reviews, then merge
1. **Claude:** `/code-review high f849c40..HEAD`. **Done**; findings below.
2. **Gemini, through `agy`:** **Done** (2026-11-25); merged below under "Merged findings".
   Raw output is in `.mmreview/20261125-agy-raw.txt` (gitignored). The notes that follow
   record why it was delayed. Originally **not run yet.** It must use Brian's Google account sign-in
   (Gemini Pro subscription), **never** the `gemini` CLI, which is set to `gemini-api-key`
   and bills per token. The Claude Code sandbox blocked `agy`'s Google hosts, and the
   auto-mode check denied allowing them. The settings change is being made in the Tools
   project; the brief is `2026-11-25-agy-sandbox-brief.md`.
3. **Merge:** combine both into one list, marking which findings both models raised and
   which only one did. Then fix, with tests first.

### Running the Gemini review
The prompt is `docs/reviews/2026-11-25-agy-review-prompt.md`. It is read-only, reads
`CLAUDE.md` and `git diff f849c40..HEAD` itself, and returns a JSON array. From the repo root:

    agy -p "$(cat docs/reviews/2026-11-25-agy-review-prompt.md)" --mode plan --print-timeout 0 < /dev/null > docs/reviews/2026-11-25-agy-review.json

Claude runs it once the sandbox allows it. Until then, Brian runs it with a `!` prefix or in
an ordinary terminal. If `agy` reports it is signed out, run `agy` once interactively.

## Claude findings (none fixed yet)

### Bugs
1. **`R/hdp_wrappers.R:160`:** `.init_for_records` matches records by
   `sprintf("%.17g", as.numeric(id))`, but ids can be any type.
   - Character ids (`"P001"`) all become `NA`, and the call fails with the misleading
     "init has a record more than once".
   - `"01"` and `"1"` collide.
   - Factor ids compare level codes.
   - No test calls `hdp_clustering()` with `init`.
2. **`R/hdp_multiseed.R:149`:** `hdp_fit_seeds` restores `.Random.seed` only if the caller
   had one. The `RNGkind()` call in its `on.exit` then reseeds from the last fit, leaving a
   fixed `.Random.seed` behind.
   - In a fresh session, `runif(1)` after `hdp_fit_seeds(c(1,2,3), f)` is the same every time.
   - Callers using `sample.kind = "Rounding"` get a warning.
   - `.with_seed()` in `utils.R` handles this correctly by removing `.Random.seed`.
3. **`R/mining.R:73`:** `arules::apriori` is called without `maxtime = 0`. Its default
   5-second limit stops enumeration on large data. It only warns, and our functions return
   the partial set as complete, which breaks the no-silent-failures rule.
4. **`R/mining.R:72`:** with `min_count` above the number of records, the support is
   `(min_n - 0.5)/nrow(x) > 1` and arules fails with an obscure validation error. It should
   return the empty table.
5. **`R/hdp_multiseed.R:138`:** seed validation and random-number save/restore are a separate
   copy of `.check_seed` / `.with_seed`. It accepts `3e9` and `Inf`, which `as.integer()`
   turns into `NA` (every fit then fails in `set.seed`). The pinned RNG kinds are also defined
   twice (`.HDP_RNG_KIND`), so the copies can drift apart.
6. **`R/network_visualization.R:244`:** uses `grDevices::col2rgb` / `rgb`, but grDevices is
   not listed under Imports in `DESCRIPTION`, so `R CMD check` will flag it.

### Cleanups
7. **`R/network_construction.R:303`:** the `cooccurrence_network` object does not record
   `test` or `correction`, although its documentation says it holds the parameters that
   produced it.
8. **`R/mining.R:64`:** `.itemset_counts` recounts every itemset with a dense `rowSums`.
   `arules::quality(found)$count` and `support(..., type = "absolute")` already give exact
   counts.
9. **British spellings:**
   - in `hdp_visualization.R:87,129,171,203` and `network_visualization.R:134-141,216-221`
   - the exported column `edge_colour` (from `plot_network_comparison`) is part of the
     public API, so renaming it **needs Brian's decision**
   - ggplot2's own `colour =` aesthetic is excluded.
10. **Three copies** of the event-table-to-incidence code (`.record_counts`, `.transactions`,
    `record_item_matrix`), each with its own checks and messages. `hdp_clustering` also checks
    `group_col` twice.

### Assessment
No parity break in the core pairwise statistics. The real problems are in the HDP layer and
mining. Proposed order once merged with Gemini: bugs 1–6 (each a failing test first), then
7, 8, 10, and 9 after the `edge_colour` decision.

## Merged findings (Claude + Gemini via `agy`, 2026-11-25)

Gemini (`agy --mode plan`) received `CLAUDE.md` and the full diff
`f849c40..HEAD -- R/ tests/ NAMESPACE DESCRIPTION` inside its prompt
(`.mmreview/20261125-agy-prompt.md`) and returned 20 findings. Codex did not run.
**Verified** means Claude reproduced the finding in R or confirmed it by reading the code.

### Both models raised (highest signal)
All six of Claude's bugs. Gemini found each one independently. **All six fixed 2026-11-25, each with a test that failed first** (`test-hdp.R`, `test-hdp-multiseed.R`, `test-mining.R`, new `test-package-metadata.R`):
- **B1** `R/hdp_wrappers.R:160`: `.init_for_records` coerces ids with `as.numeric`.
- **B2** `R/hdp_multiseed.R:149-153`: `hdp_fit_seeds` leaves `.Random.seed` behind when the
  caller had none.
- **B3** `R/mining.R:73`: no `maxtime = 0`, so the result is silently truncated at 5 s.
- **B4** `R/mining.R:72`: `min_count` above the number of records gives support > 1, and
  arules fails.
- **B5** `R/hdp_multiseed.R:144`: seed validation accepts `Inf` and `2^31`, which become `NA`.
- **B6** `R/network_visualization.R:242`: grDevices is used but not listed in `DESCRIPTION`
  Imports.

### Gemini only: bugs
G1–G10 **fixed 2026-11-25**, each with a test that failed first.

- **G1** `R/hdp_multiseed.R:25`: `adjusted_rand_index` computes `na * max(ib)` in integer
  arithmetic and allocates a dense `na × nb` table. **Verified:**
  `adjusted_rand_index(1:50000, 1:50000)` gives "NAs produced by integer overflow". It needs
  more than ~46k distinct labels, so it is unlikely with cluster labels, but the fix is cheap:
  count pairs with `table()` or a sparse key.
- **G2** `R/network_metrics.R:29`: `modularity_q` maps `NA` assignments to a real community
  through `match()`. **Verified:** two triangles joined by an edge, with
  `c(1,1,1,NA,NA,NA)`, return 0.357 instead of raising.
- **G3** `R/hdp_multiseed.R:236`: `hdp_pairwise_ari` validates fits by name, so a duplicate
  name checks element 1 twice and never checks element 2. A raw `hdp_bernoulli_result` has
  no `record_assignments`, so `identical(NULL, NULL)` passes. **Verified by reading.**
- **G4** `R/network_group_stratification.R:118-128`: rows of excluded items are dropped
  before any `NA`-id check, so an `NA` id on such a row disappears silently. **Verified by
  reading.**
- **G5** `R/hdp_wrappers.R:214`: `hdp_active_clusters` on a raw `hdp_bernoulli_result`
  returns `integer(0)`; it should raise a class error.
- **G6** `R/network_visualization.R:145`: `plot_network_comparison` misses duplicate names
  and `NA` names in `strata`. With duplicates, the second stratum is silently ignored.
- **G7** `R/hdp_multiseed.R:320`: `hdp_cluster_recurrence(reference = 1.5 or TRUE)` is
  coerced to 1 instead of raising.
- **G8** `R/hdp_clustering.R:316`: `init` containing `Inf` passes validation, and the
  responsibilities become NaN.
- **G9** `R/hdp_clustering.R:350`: a `+Inf` ELBO would win, contrary to the comment.
  Theoretical, but `is.finite()` costs nothing.
- **G10** `R/network_metrics.R:20` (`Inf` weights pass validation) and
  `R/hdp_multiseed.R:214` (`rel_tol = Inf` passes). Low severity.

### Gemini only: tests that cannot fail
T1–T3 **fixed 2026-11-25**. Each strengthened test now fails under the break the finding named: a bad aesthetic, a hardcoded method and seed, or an ignored resolution.

- **T1** `tests/testthat/test-hdp-visualization.R:10`: the three HDP plots are checked only
  with `expect_s3_class`, never with `ggplot_build()`, so a bad aesthetic still passes.
- **T2** `tests/testthat/test-communities.R:140`: the seed and method passthrough test uses
  `make_test_event_df()`, whose partition is the same under any method or seed. Hardcoding
  either one inside the pipeline would still pass.
- **T3** `tests/testthat/test-communities.R:116`: no test shows that a non-default
  `resolution` changes the partition.

### Claude only
Cleanups 7–10 above: test and correction missing from the network object, the duplicate
itemset recount, British spellings and `edge_colour`, and the three copies of the incidence
code.

Cleanups 7, 8 and 10 **done 2026-11-25**:
- **7:** the network stores `test` and `correction`. `cooccurrence_network_from_edge_data()`
  now requires both, since `edge_data` does not record them and the rebuild must stay exact.
  This breaks existing calls to it.
- **8:** itemset and rule counts come from arules (`quality()$count`, and
  `support(type = "absolute")` for rule left-hand sides). Adding 1 to any of the three
  fails between 3 and 11 mining tests.
- **10:** `.event_incidence()` and `.multi_item_matrix()` in `utils.R` replace the three
  copies, and `.check_group_col()` replaces the group checks in `record_item_matrix` and
  `stratified_network_analysis` and the duplicate in `hdp_clustering`. Every entry point now
  raises the same messages; a new test in `test-validation.R` pins them.

Cleanup 9 waits on the `edge_colour` decision.

### Assessment after the merge
The two models agree on every bug Claude found, and Gemini adds two more that were
reproduced (G1, G2) and several silent-failure gaps of the same kind (G3–G8). No finding
touches the parity rules for the core pairwise statistics. Proposed fix order, each with a
failing test first: B1–B6, G2, G3, G4, G1, then G5–G10, T1–T3, and the cleanups. Rename
`edge_colour` only after Brian decides.
