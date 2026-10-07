# Phase 3: Leiden via igraph Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the hand-written phase-1-only Louvain with igraph's full Leiden (default) and Louvain, with a required, locally applied seed, and bring the downstream CancerClustering project along.

**Architecture:** `detect_communities()` builds an igraph graph from `net$adjacency` and calls `igraph::cluster_leiden()` (modularity objective) or `igraph::cluster_louvain()` inside a helper that seeds R's RNG and restores the caller's stream afterward. Labels are renumbered by first appearance in `net$items` order, so they depend only on the partition. `stratified_network_analysis()` and `run_network_pipeline()` take `seed` and pass it through. The old `.louvain()`, `.modularity_gain()` and `.neighbors()` are deleted, along with `max_iter` and the `n_gain_ties` field.

**Tech Stack:** R, igraph (already in Imports), testthat 3e, lintr.

**Spec:** `docs/specs/2026-10-06-r-ecosystem-refactor.md` (§R1, D1, Phase 3). Read §2 (governing rule) and §R1 before starting.

**Scope:** this plan covers R1 only. The spec's other Phase 3 items wait: R9 (flat-K mixture) needs decision D2, and R12 (the HDP sanity test against `flexmix`) needs tolerances agreed with Brian. Task 3 records both as open.

**Downstream:** Task 4 makes the matching changes in CancerClustering on its own branch. The two repos merge together (Task 5), so CancerClustering is never broken against this package's `main`.

**One deliberate change to §R1:** the spec keeps "renumbered `1..K` by ascending raw label". igraph's raw labels are arbitrary, so this plan numbers labels by **first appearance in `net$items` order** instead. That makes the labels, and the community colors in the plots, a function of the partition alone. Task 3 updates the spec to match.

## Global Constraints

- Default method is `"leiden"`; `"louvain"` is the only other method; `"label_propagation"` keeps its "exists in the Julia reference but is not ported" raise; anything else raises via `.check_choice()`.
- Leiden is always called with `objective_function = "modularity"` and `n_iterations = -1`. igraph's defaults (CPM, 2 iterations) are wrong for this package.
- `seed` has **no default**. A missing `seed` raises. `seed = NULL` means "draw from the caller's RNG stream" (explicit, not a default). A whole number seeds locally and the caller's `.Random.seed` is restored, or removed if it didn't exist.
- `modularity` in the result is always `modularity_q()` of the returned partition, never igraph's own value.
- Return shape: `assignments` (integer, named by item, `net$items` order), `communities` (named `"1".."K"`, members C-locale sorted), `modularity`, `n_communities`. No `n_gain_ties`.
- No `suppressWarnings`, no `try`, `vapply` not `sapply` (CLAUDE.md "No silent failures").
- `NAMESPACE` is maintained by hand. No exports are added or removed in this plan.
- Expected partitions **cannot** come from running Julia (this is a deliberate divergence). They come from fixtures whose partition is unambiguous by construction, each verified under seeds 1..100 while this plan was written.
- Every task ends with `Rscript -e 'testthat::test_local(".", reporter = "check")'` showing `FAIL 0 | WARN 0` and `Rscript -e 'cat(length(lintr::lint_package()), "lints\n")'` showing `0 lints`.

## Facts verified while writing this plan (igraph 2.3.4)

- Same seed gives the same partition for both methods; different seeds can differ.
- `cluster_leiden(objective_function = "modularity", weights = ...)` returns a partition whose `igraph::modularity()` equals `modularity_q()`, so the weights are used.
- **`cluster_leiden(..., n_iterations = -1)` errors on a graph with no edges** ("missing value where TRUE/FALSE needed"). `detect_communities()` must handle the edgeless case before calling igraph.
- Isolated vertices alongside edges are fine (each becomes its own community).
- `two_triangles()` (bridge weight 0.1) and two disjoint triangles give partition `111222` under all 100 seeds, for both methods. A triangle gives `111`.
- Ring of 20 triangles: the old phase-1 code reaches Q = 0.70 (20 communities). Leiden reaches 0.77–0.775 and Louvain 0.77–0.775 under seeds 1..20. The partition itself varies by seed (K = 8–10), so tests assert Q, not membership.
- `set.seed(2^31)` errors; negative whole numbers are accepted.

## Review Focus

1. **Edgeless graph with vertices:** igraph Leiden crashes. Expected: one community per vertex, modularity 0, for both methods. Test in Task 2.
2. **Caller has no `.Random.seed`** (a fresh session): calling with a whole-number seed must not leave one behind. Test in Task 1.
3. **Out-of-range or malformed seed** (`2^31`, `1.5`, `NA`, `c(1, 2)`, `"1"`): it must raise this package's message, not `set.seed()`'s coercion warning. Test in Task 1.
4. **Disconnected components:** two disjoint triangles give 2 communities under every seed. Test in Task 2.
5. **Empty network (0 items) inside stratified analysis:** returns the empty result without touching igraph, whatever the seed. Test in Task 2.

---

### Task 0: Pre-flight (human checks, no code)

- [ ] **Step 1: Branch from main**

```bash
cd /Users/brainchapman/Code/Julia/cooccurrence-analysis-r
git checkout main && git checkout -b refactor/3-leiden
```

- [ ] **Step 2: Know the downstream consumer.** `~/Code/Julia/CancerClustering` uses this package through `rpkg/SEERClustering` and its pixi tasks. Those tasks run `pkgload::load_all("../cooccurrence-analysis-r")`, so CancerClustering picks up **whatever branch is checked out here**. From Task 2 on, CancerClustering's R tests break until Task 4 lands. Don't run CancerClustering's pipeline between Task 2 and Task 4 except as Task 4 directs.

- [ ] **Step 3: Pin the igraph floor.** The code passes `resolution =` to `cluster_leiden()`. Confirm the argument exists in the installed version:

Run: `Rscript -e 'cat(packageVersion("igraph"), "resolution" %in% names(formals(igraph::cluster_leiden)), "\n")'`
Expected: `2.3.4 TRUE` (or later). In `DESCRIPTION`, change `    igraph,` to `    igraph (>= 2.3.4),`, using the installed version unless you have verified that an older release has the `resolution` argument.

---

### Task 1: Seed helpers

**Files:**
- Modify: `R/utils.R` (append)
- Create: `tests/testthat/test-communities.R`

**Interfaces:**
- Produces: `.check_seed(seed, caller)` raises unless `seed` is `NULL` or one whole number in `[-.Machine$integer.max, .Machine$integer.max]`, and returns `invisible(seed)`. `.with_seed(seed, expr)` evaluates `expr` with the RNG seeded, restores the caller's stream, and returns `expr`'s value; with `seed = NULL` it evaluates `expr` on the caller's stream.

- [ ] **Step 1: Write the failing tests** in a new `tests/testthat/test-communities.R`:

```r
# Community detection via igraph (Leiden default, Louvain option). These
# replace the Julia-parity tests of the old phase-1-only Louvain, a
# documented divergence (CLAUDE.md, "Divergences from the Julia reference").
# Expected partitions come from fixtures that are unambiguous by
# construction, not from the Julia reference.

test_that("R-only: .with_seed seeds locally and restores the caller's stream", {
  set.seed(42)
  expected_next <- stats::runif(1)
  set.seed(42)
  a <- .with_seed(1L, stats::runif(3))
  expect_identical(stats::runif(1), expected_next)
  expect_identical(.with_seed(1L, stats::runif(3)), a)
})

test_that("R-only: .with_seed(NULL) draws from the caller's stream", {
  set.seed(3)
  a <- .with_seed(NULL, stats::runif(2))
  set.seed(3)
  expect_identical(stats::runif(2), a)
})

test_that("R-only: .with_seed leaves no .Random.seed behind in a fresh session", {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) saved <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(if (had) assign(".Random.seed", saved, envir = globalenv()))
  if (had) rm(".Random.seed", envir = globalenv())
  .with_seed(1L, stats::runif(1))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("R-only: .check_seed accepts NULL and whole numbers, raises otherwise", {
  expect_silent(.check_seed(NULL, "f"))
  expect_silent(.check_seed(7, "f"))
  expect_silent(.check_seed(-5L, "f"))
  for (bad in list(1.5, NA, NA_integer_, c(1, 2), "1", 2^31, Inf)) {
    expect_error(.check_seed(bad, "f"), "f: seed must be NULL or one whole number")
  }
})
```

- [ ] **Step 2: Run them to confirm they fail**

Run: `Rscript -e 'testthat::test_local(".", filter = "communities", reporter = "check")'`
Expected: FAIL, `could not find function ".with_seed"`.

- [ ] **Step 3: Implement.** Append to `R/utils.R`:

```r
#' Raise unless `seed` is NULL or one whole number `set.seed()` accepts
#' @keywords internal
.check_seed <- function(seed, caller) {
  ok <- is.null(seed) ||
    (is.numeric(seed) && length(seed) == 1L && is.finite(seed) && seed == round(seed) &&
       abs(seed) <= .Machine$integer.max)
  if (!ok) {
    stop(sprintf(paste0("%s: seed must be NULL or one whole number within the integer range ",
                        "(NULL draws from the current RNG stream)"), caller), call. = FALSE)
  }
  invisible(seed)
}

#' Evaluate `expr` with the RNG seeded to `seed`, then restore the caller's
#' stream (or remove `.Random.seed` if the caller had none), so a seeded call
#' never shifts the caller's later random draws. `seed = NULL` evaluates
#' `expr` on the caller's stream, which it advances.
#' @keywords internal
.with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  env <- globalenv()
  had <- exists(".Random.seed", envir = env, inherits = FALSE)
  if (had) saved <- get(".Random.seed", envir = env, inherits = FALSE)
  on.exit(if (had) assign(".Random.seed", saved, envir = env) else
    if (exists(".Random.seed", envir = env, inherits = FALSE)) rm(".Random.seed", envir = env))
  set.seed(seed)
  expr
}
```

`expr` is lazy: it is evaluated at the final `expr`, after `set.seed()`. Don't touch `expr` earlier.

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `Rscript -e 'testthat::test_local(".", filter = "communities", reporter = "check")'`
Expected: `FAIL 0 | WARN 0`.

- [ ] **Step 5: Full suite and lint, then commit**

```bash
Rscript -e 'testthat::test_local(".", reporter = "check")'
Rscript -e 'cat(length(lintr::lint_package()), "lints\n")'
git add R/utils.R tests/testthat/test-communities.R DESCRIPTION
git commit -m "Seed helpers for reproducible igraph community detection"
```

---

### Task 2: `detect_communities()` on igraph, seed threaded through callers

The suite must stay green at the commit, and a seedless `detect_communities()` breaks every caller. So the function, its two callers, the deleted internals and all the affected tests change together in this task.

**Files:**
- Modify: `R/network_metrics.R:33-166` (delete `.neighbors`, `.modularity_gain`, `.louvain`; rewrite `detect_communities`)
- Modify: `R/network_group_stratification.R:89-127` (`community_method`, `seed`)
- Modify: `R/pipeline.R` (`community_method`, `seed`)
- Modify: `tests/testthat/helper-fixtures.R` (add `ring_of_cliques`)
- Modify: `tests/testthat/test-communities.R` (append)
- Modify: `tests/testthat/test-network.R`, `tests/testthat/test-validation.R`, `tests/testthat/test-network-stratification.R`, `tests/testthat/test-r-ecosystem-crosscheck.R`

**Interfaces:**
- Consumes: `.check_seed(seed, caller)` and `.with_seed(seed, expr)` from Task 1; `communities_from_assignments(net, assignments)` and `modularity_q()` (unchanged).
- Produces:
  - `detect_communities(net, method = "leiden", seed, resolution = 1)`;
  - `stratified_network_analysis(event_df, exclusive_items, group_col = "Group", weight_metric = "lift", min_count = 30L, alpha = 0.05, test = "fisher", correction = "bh", community_method = "leiden", timing_filter = "all", seed)`;
  - `run_network_pipeline(event_df, weight_metric = "lift", min_count = 30L, alpha = 0.05, community_method = "leiden", stratify_by_group = TRUE, test = "fisher", correction = "bh", seed)`.

  `seed` goes **last** in each signature, so existing positional calls keep working.

- [ ] **Step 1: Add the ring fixture** to the end of `tests/testthat/helper-fixtures.R`:

```r
#' Ring of `n_cliques` complete graphs of `k` nodes, each joined to the next
#' by one unit edge. Phase-1-only local moving stops at one community per
#' clique; full Louvain/Leiden aggregates neighbors and scores higher.
ring_of_cliques <- function(n_cliques, k) {
  n <- n_cliques * k
  adj <- matrix(0, n, n)
  for (c in seq_len(n_cliques) - 1L) {
    idx <- c * k + seq_len(k)
    adj[idx, idx] <- 1
    nxt <- (c + 1L) %% n_cliques * k + 1L
    adj[c * k + k, nxt] <- adj[nxt, c * k + k] <- 1
  }
  diag(adj) <- 0
  adj
}
```

`net_from_adjacency()` and `two_triangles()` live at the top of `test-network.R`, and testthat files don't share scope. **Move both** into `helper-fixtures.R` unchanged (cut from `test-network.R` lines 9–25, paste at the end of the helper), so `test-communities.R` can use them.

- [ ] **Step 2: Write the failing tests.** Append to `tests/testthat/test-communities.R`:

```r
canon <- function(m) match(m, unique(m))

test_that("R-only: detect_communities finds two bridged triangles under any seed", {
  net <- net_from_adjacency(two_triangles())
  for (method in c("leiden", "louvain")) {
    for (s in 1:5) {
      comm <- detect_communities(net, method = method, seed = s)
      expect_identical(unname(comm[["assignments"]]), c(1L, 1L, 1L, 2L, 2L, 2L))
      expect_identical(comm[["communities"]],
                       list(`1` = c("n1", "n2", "n3"), `2` = c("n4", "n5", "n6")))
      expect_identical(names(comm), c("assignments", "communities", "modularity",
                                      "n_communities"))
      expect_equal(comm[["modularity"]], modularity_q(two_triangles(), c(1, 1, 1, 2, 2, 2)))
    }
  }
})

test_that("R-only: disconnected components are separate communities", {
  adj <- two_triangles()
  adj[3, 4] <- adj[4, 3] <- 0
  for (s in 1:5) {
    expect_identical(unname(detect_communities(net_from_adjacency(adj), seed = s)[["assignments"]]),
                     c(1L, 1L, 1L, 2L, 2L, 2L))
  }
})

test_that("R-only: full Leiden/Louvain beat the reference's phase-1 local moving", {
  # 20 triangles in a ring. Phase-1-only local moving (the Julia reference,
  # and this package before the change) stops at Q = 0.70 with one community
  # per triangle; aggregation reached Q >= 0.77 under every seed in 1..20.
  net <- net_from_adjacency(ring_of_cliques(20L, 3L))
  for (method in c("leiden", "louvain")) {
    for (s in 1:5) {
      comm <- detect_communities(net, method = method, seed = s)
      expect_gt(comm[["modularity"]], 0.75)
      expect_lt(comm[["n_communities"]], 20L)
    }
  }
})

test_that("R-only: labels are numbered by first appearance in item order", {
  net <- net_from_adjacency(ring_of_cliques(20L, 3L))
  for (s in 1:10) {
    a <- unname(detect_communities(net, seed = s)[["assignments"]])
    expect_identical(unique(a), seq_len(max(a)))
  }
})

test_that("R-only: same seed, same partition; seed = NULL follows the caller's stream", {
  net <- net_from_adjacency(ring_of_cliques(20L, 3L))
  expect_identical(detect_communities(net, seed = 11L), detect_communities(net, seed = 11L))
  set.seed(5)
  a <- detect_communities(net, seed = NULL)
  set.seed(5)
  expect_identical(detect_communities(net, seed = NULL), a)
})

test_that("R-only: seed is required, and bad methods and resolutions raise", {
  net <- net_from_adjacency(two_triangles())
  expect_error(detect_communities(net), "seed is required")
  expect_error(detect_communities(net, seed = 1.5), "seed must be NULL or one whole number")
  expect_error(detect_communities(net, method = "walktrap", seed = 1L), "unsupported method")
  expect_error(detect_communities(net, method = "label_propagation", seed = 1L),
               "label_propagation' exists in the Julia reference but is not ported")
  for (bad in list(0, -1, NA_real_, c(1, 2), "1")) {
    expect_error(detect_communities(net, seed = 1L, resolution = bad), "resolution")
  }
})

test_that("R-only: edgeless and empty networks need no igraph call", {
  # igraph's cluster_leiden(n_iterations = -1) errors on a graph with no edges.
  for (method in c("leiden", "louvain")) {
    comm <- detect_communities(net_from_adjacency(matrix(0, 3, 3)), method = method, seed = 1L)
    expect_identical(unname(comm[["assignments"]]), 1:3)
    expect_identical(comm[["modularity"]], 0)
  }
  empty <- net_from_adjacency(matrix(0, 0, 0), character(0))
  expect_identical(detect_communities(empty, seed = 1L)[["n_communities"]], 0L)
  expect_identical(detect_communities(empty, seed = NULL)[["communities"]], list())
})

test_that("R-only: an isolated vertex beside edges gets its own community", {
  adj <- matrix(0, 3, 3)
  adj[1, 2] <- adj[2, 1] <- 1
  comm <- detect_communities(net_from_adjacency(adj), seed = 1L)
  expect_identical(unname(comm[["assignments"]]), c(1L, 1L, 2L))
})

test_that("R-only: stratified analysis and the pipeline pass seed and method through", {
  ev <- make_test_event_df()
  res <- stratified_network_analysis(ev, exclusive_items = NULL, min_count = 1L, alpha = 1,
                                     community_method = "louvain", seed = 3L)
  for (g in names(res)) {
    expect_identical(res[[g]][["communities"]],
                     detect_communities(res[[g]][["network"]], method = "louvain", seed = 3L))
  }
  expect_error(stratified_network_analysis(ev, exclusive_items = NULL), "seed is required")
  piped <- run_network_pipeline(ev, min_count = 1L, alpha = 1, stratify_by_group = FALSE,
                                seed = 3L)
  expect_identical(piped[["communities"]], detect_communities(piped[["network"]], seed = 3L))
  expect_error(run_network_pipeline(ev, stratify_by_group = FALSE), "seed is required")
})
```

`net_from_adjacency(matrix(0, 0, 0), character(0))` must produce a valid empty network. If `dimnames<-` fails on a 0×0 matrix with `list(character(0), character(0))`, build the empty network inline in that test as `structure(list(items = character(0), adjacency = matrix(0, 0, 0)), class = "cooccurrence_network")`.

- [ ] **Step 3: Run them to confirm they fail**

Run: `Rscript -e 'testthat::test_local(".", filter = "communities", reporter = "check")'`
Expected: FAIL. Most fail on `unused argument (seed = ...)`; the stratified test fails on `seed`.

- [ ] **Step 4: Rewrite `detect_communities`.** In `R/network_metrics.R`, delete everything from the `#' Neighbors of vertex` roxygen block (line 33) through the end of `detect_communities` (line 166). That removes `.neighbors`, `.modularity_gain`, `.louvain` and the old `detect_communities`. Insert:

```r
#' Detect communities in a co-occurrence network
#'
#' `method = "leiden"` (default) runs [igraph::cluster_leiden()] with the
#' modularity objective, iterated to convergence; `"louvain"` runs
#' [igraph::cluster_louvain()]. Both use the edge weights and `resolution`.
#' Leiden is the default because it guarantees connected communities, which
#' Louvain does not. The Julia reference's `"label_propagation"` is not
#' ported; asking for it, or anything else, raises.
#'
#' **Divergence from the Julia reference:** Julia runs only Louvain's first
#' phase (local moving, no aggregation), which can stall below the optimum,
#' so partitions differ from Julia's by design (CLAUDE.md, "Divergences").
#'
#' Both algorithms are randomized. `seed` has no default: pass a whole number
#' for a reproducible partition (applied locally; the caller's RNG stream is
#' restored), or `NULL` to draw from the current stream. A graph with no
#' edges gets one community per vertex without calling igraph.
#'
#' Returns a list: `assignments` (integer, named by item, in `net$items`
#' order, labels `1..K` numbered by first appearance), `communities` (a list
#' named `"1"`, `"2"`, ..., each the sorted member items), `modularity`
#' ([modularity_q()] of the partition), and `n_communities`.
#' @export
detect_communities <- function(net, method = "leiden", seed, resolution = 1) {
  stopifnot(inherits(net, "cooccurrence_network"))
  if (identical(method, "label_propagation")) {
    stop(paste0("method = 'label_propagation' exists in the Julia reference but is not ",
                "ported to R: it is nondeterministic, and the analysis uses 'leiden'"),
         call. = FALSE)
  }
  .check_choice(method, c("leiden", "louvain"), "method")
  if (missing(seed)) {
    stop(paste0("detect_communities: seed is required: pass a whole number for a ",
                "reproducible partition, or NULL to draw from the current RNG stream"),
         call. = FALSE)
  }
  .check_seed(seed, "detect_communities")
  if (!(is.numeric(resolution) && length(resolution) == 1L && is.finite(resolution) &&
          resolution > 0)) {
    stop("detect_communities: resolution must be one positive number", call. = FALSE)
  }
  items <- net[["items"]]
  if (length(items) == 0L) {
    return(list(assignments = stats::setNames(integer(0), character(0)),
                communities = list(), modularity = 0, n_communities = 0L))
  }
  adj <- net[["adjacency"]]
  raw <- if (all(adj == 0)) {
    seq_along(items)
  } else {
    g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE,
                                             diag = FALSE)
    w <- igraph::E(g)$weight
    fit <- .with_seed(seed, if (identical(method, "leiden")) {
      igraph::cluster_leiden(g, objective_function = "modularity", weights = w,
                             resolution = resolution, n_iterations = -1)
    } else {
      igraph::cluster_louvain(g, weights = w, resolution = resolution)
    })
    as.integer(igraph::membership(fit))
  }
  communities_from_assignments(net, stats::setNames(match(raw, unique(raw)), items))
}
```

- [ ] **Step 5: Thread `seed` and the new default through the callers.**

In `R/network_group_stratification.R`:
- In the signature, change `community_method = "louvain", timing_filter = "all") {` to `community_method = "leiden", timing_filter = "all", seed) {`.
- Directly after the `if (missing(exclusive_items)) { ... }` block, add:

```r
  if (missing(seed)) {
    stop(paste0("stratified_network_analysis: seed is required: pass a whole number for ",
                "reproducible partitions (every group uses it), or NULL to draw from the ",
                "current RNG stream"), call. = FALSE)
  }
  .check_seed(seed, "stratified_network_analysis")
```

- Change `communities = detect_communities(net, method = community_method),` to `communities = detect_communities(net, method = community_method, seed = seed),`.
- In its roxygen block, add the line: `#' \`seed\` is passed to [detect_communities()] unchanged for every group, so a group's partition does not depend on which groups precede it.`

In `R/pipeline.R`:
- Change `alpha = 0.05, community_method = "louvain",` to `alpha = 0.05, community_method = "leiden",`, and `stratify_by_group = TRUE, test = "fisher", correction = "bh") {` to `stratify_by_group = TRUE, test = "fisher", correction = "bh", seed) {`.
- After the existing `if (!identical(stratify_by_group, FALSE)) { ... }` block (so the `compare_networks` raise still comes first), add:

```r
  if (missing(seed)) {
    stop(paste0("run_network_pipeline: seed is required: pass a whole number for a ",
                "reproducible partition, or NULL to draw from the current RNG stream"),
         call. = FALSE)
  }
```

- Change `communities = detect_communities(net, method = community_method),` to `communities = detect_communities(net, method = community_method, seed = seed),`.
- In the roxygen block, extend the last sentence: `` `test` and `correction` go to [build_cooccurrence_network()]; `community_method` and `seed` to [detect_communities()]. ``

- [ ] **Step 6: Update the existing tests.** Change exactly these; no other existing assertion changes.

`tests/testthat/test-network.R`:
- **Line 96, "detect_communities: louvain":** rename it to `"detect_communities: leiden and louvain"`. Wrap the body in `for (method in c("leiden", "louvain")) { ... }`, call `detect_communities(net, method = method, seed = 1L)`, and keep the four assertions.
- **Line 105, "louvain modularity gain matches delta Q (C1, C2, T6)":** delete the whole test. It tests the deleted `.modularity_gain` and `.louvain`; the ring-of-cliques test replaces its coverage.
- **Line 130, "renumbering, names, and unsupported methods (R)":** change `detect_communities(net_from_adjacency(two_triangles()))` to `detect_communities(net_from_adjacency(two_triangles()), seed = 1L)`. In the `label_propagation` `expect_error`, add `seed = 1L`.
- **Line 142, communities_from_assignments:** `comm <- detect_communities(net, seed = 1L)`, and replace `comm[setdiff(names(comm), "n_gain_ties")]` with `comm`.
- **Lines 191 and 197 (plots):** `detect_communities(net, seed = 1L)`.

`tests/testthat/test-validation.R`:
- **"Louvain counts exact-gain ties (R)" (line 179):** delete the whole test.
- **"network_integration" (line 165):** in the first `run_network_pipeline(...)` call (the one with `stratify_by_group = FALSE`), add `seed = 1L`. Leave the two `expect_error(... "compare_networks")` calls unchanged; they raise before the seed check.
- **"network with Bonferroni and Holm corrections" (around line 105):** add `seed = 1L` to the `run_network_pipeline(...)` call.
- **"detect_communities: label_propagation":** unchanged. The method check raises before the seed check.

`tests/testthat/test-network-stratification.R`:
- Add `seed = 1L` to every `stratified_network_analysis(` call (20 calls). Apply it mechanically and check the count:

```bash
python3 - <<'EOF'
p = "tests/testthat/test-network-stratification.R"
s = open(p).read()
n = s.count("stratified_network_analysis(")
s = s.replace("stratified_network_analysis(", "stratified_network_analysis(seed = 1L, ")
open(p, "w").write(s)
print(n, "calls seeded")
EOF
```

  Expected output: `20 calls seeded`. A named argument first is valid R. The `exclusive_items` check runs before the seed check, so `"exclusive_items has no default: forgetting it errors (R)"` still gets its expected error. Re-read each `expect_error(stratified_network_analysis(seed = 1L, ...` and confirm its regex still targets the intended error.
- **`two_triangle_stratum()` (line ~244):** `detect_communities(net, seed = 1L)`.
- **Line 213** `expect_error(run_network_pipeline(make_test_event_df()), "compare_networks")`: unchanged.

`tests/testthat/test-r-ecosystem-crosscheck.R`: replace the test `"detect_communities recovers planted blocks and is competitive with igraph Louvain"` with:

```r
test_that("detect_communities delegates to igraph: same seed, same partition", {
  skip_if_not_installed("igraph")
  pn <- planted_net()
  g <- igraph::graph_from_adjacency_matrix(pn$net$adjacency, mode = "undirected",
                                           weighted = TRUE, diag = FALSE)
  w <- igraph::E(g)$weight
  canon <- function(m) match(m, unique(m))
  set.seed(7L)
  ref_louvain <- igraph::membership(igraph::cluster_louvain(g, weights = w))
  set.seed(7L)
  ref_leiden <- igraph::membership(igraph::cluster_leiden(g, objective_function = "modularity",
                                                          weights = w, n_iterations = -1))
  louvain <- detect_communities(pn$net, method = "louvain", seed = 7L)
  leiden <- detect_communities(pn$net, seed = 7L)
  expect_identical(unname(louvain$assignments), canon(as.integer(ref_louvain)))
  expect_identical(unname(leiden$assignments), canon(as.integer(ref_leiden)))
  for (fit in list(louvain, leiden)) {
    expect_gt(adjusted_rand(unname(fit$assignments), pn$truth), 0.9)
    expect_equal(fit$modularity, igraph::modularity(g, membership = unname(fit$assignments),
                                                    weights = w), tolerance = 1e-10)
  }
})
```

- [ ] **Step 7: Confirm nothing references the deleted code**

Run: `grep -rnE "\.louvain|modularity_gain|\.neighbors|n_gain_ties|max_iter" R tests`
Expected: no matches, apart from `max_iter` inside `R/hdp_*.R` and the HDP tests, which belong to the HDP and stay.

- [ ] **Step 8: Run the full suite and lint**

```bash
Rscript -e 'testthat::test_local(".", reporter = "check")'
Rscript -e 'cat(length(lintr::lint_package()), "lints\n")'
```

Expected: `FAIL 0 | WARN 0` and `0 lints`. If a stratification test fails on an `expect_error` regex, the seed insertion changed which check fires first. Fix the call, not the regex.

- [ ] **Step 9: Commit**

```bash
git add R tests DESCRIPTION
git commit -m "Leiden (default) and Louvain via igraph replace phase-1-only Louvain"
```

---

### Task 3: Documentation of the divergence

**Files:**
- Modify: `CLAUDE.md`
- Modify: `docs/specs/2026-10-06-r-ecosystem-refactor.md`
- Create: `docs/julia-issue-full-louvain.md` (draft only)

**Interfaces:**
- Consumes: the behavior shipped in Task 2.

- [ ] **Step 1: CLAUDE.md.** Make four edits:
  - In "What is ported", replace `Louvain communities` with `Leiden (default) or Louvain communities via igraph`.
  - In "Architecture", replace `` `network_metrics.R` (modularity, Louvain, centrality) `` with `` `network_metrics.R` (modularity, communities via igraph, centrality) ``.
  - In "Parity rules", delete the whole `- **Louvain tie-breaking differs from Julia**: ...` bullet.
  - In "Divergences from the Julia reference", append:

```markdown
- **Communities use full Leiden (default) or Louvain from igraph.** Julia runs only Louvain's first phase (local moving, no aggregation), which stalls below the optimum: on a ring of 20 triangles it stops at Q = 0.70 where aggregation reaches 0.77. Partitions therefore differ from Julia's, and `community_method`'s default is now `"leiden"`. Both algorithms are randomized, so `detect_communities()`, `stratified_network_analysis()` and `run_network_pipeline()` require `seed` (a whole number, applied locally, or `NULL` for the caller's stream). Expected partitions in tests come from fixtures unambiguous by construction, not from Julia.
```

  Also add `seed` to the "Anything the caller knows..." convention bullet: after the `exclusive_items` sentence, add `` `seed` likewise has no default: a randomized result needs the caller to choose reproducibility. ``

- [ ] **Step 2: Spec status.** In `docs/specs/2026-10-06-r-ecosystem-refactor.md`:
  - change the status line to `**Status:** phases 0-2 and R1 implemented · ...`;
  - under Phase 3, add `*R1 done 2026-10-06. Open: R9 (needs D2) and R12 (needs tolerances).*`;
  - in §R1, replace `assignments renumbered `1..K` by ascending raw label` with `assignments renumbered `1..K` by first appearance in `net$items` order (igraph's raw labels are arbitrary)`.

- [ ] **Step 3: Draft (do not file) the Julia issue.** Create `docs/julia-issue-full-louvain.md`:

```markdown
# Louvain stops after phase 1 (no aggregation)

`network_metrics.jl:180` documents "Simple Louvain community detection (phase 1 only — sufficient for small graphs)". Without aggregation, local moving can stall below the modularity optimum. On a ring of 20 triangles (each joined to the next by one unit edge), phase 1 stops at one community per triangle, Q = 0.70; full Louvain or Leiden reaches Q ≈ 0.77. The R port (`cooccurrence-analysis-r`) now uses igraph's Leiden (default) and Louvain and records this as a divergence. Proposal: add the aggregation phase, or call a full implementation, so the two packages agree again on method if not on exact partitions (both algorithms are randomized).
```

Filing it on GitHub (`chapmanbe/cooccurrence-analysis`) is outward-facing. **Ask Brian before filing**, and don't run `gh issue create` without his yes.

- [ ] **Step 4: Lint, run the suite, commit**

```bash
Rscript -e 'testthat::test_local(".", reporter = "check")'
git add CLAUDE.md docs
git commit -m "Document the Leiden divergence from the Julia reference"
```

---

### Task 4: CancerClustering follows the change

Work in `~/Code/Julia/CancerClustering`. **Read its `AGENTS.md` "Hard rules" first.** These apply throughout:
- rule 1: no SEER data leaves the machine, so print only aggregates of the kind the manuscript already publishes (modularity, community counts, ARI), never counts per site or per patient;
- rule 3: only Snakemake rules write to `results/`;
- rule 4: never type a number into the manuscript;
- rule 6: every artifact's rule lists its outputs;
- rule 9: set and record seeds explicitly.

**Files:**
- Modify: `rpkg/SEERClustering/R/sex_anatomy.R:97-115` (`seer_stratified_network_analysis`)
- Modify: `scripts/fit_network.R:48,60-64`, `scripts/fit_stratified_network.R:53-56,73-74`
- Modify: `port.smk` and `Snakefile` (outputs named `-louvain.json`)
- Modify: `contract/port_tolerances.yaml` (network_v2 and network_stratified_v2 blocks); in outcome B also `contract/compare_port.py` and `contract/test_compare_port.py`
- Modify: `rpkg/SEERClustering/tests/testthat/*` (calls that now need `seed`)
- Modify: `PORTING.md` ("Notes")

**Interfaces:**
- Consumes: `detect_communities(net, method = "leiden", seed, resolution = 1)` and `stratified_network_analysis(..., community_method = "leiden", ..., seed)` from Task 2.
- Produces: `seer_stratified_network_analysis(event_df, ontology, exclusive_items = sex_exclusive_sites(ontology), weight_metric = "phi", min_count = 30L, alpha = 0.05, test = "fisher", correction = "bh", community_method = "leiden", seed)`, with `seed` last and no default.

- [ ] **Step 1: Branch**

```bash
cd ~/Code/Julia/CancerClustering
git checkout main && git checkout -b port/leiden-communities
```

- [ ] **Step 2: Measure against the frozen Julia partitions (decision gate input).** This reads only the frozen Julia bundles under `data/derived/port_baselines/` and never loads raw SEER data. Save the script to the session scratchpad, **not** `results/` (rule 3):

```r
# Run from ~/Code/Julia/CancerClustering.
pkgload::load_all("../cooccurrence-analysis-r", quiet = TRUE)
pkgload::load_all("rpkg/SEERClustering", export_all = FALSE, quiet = TRUE)
ari <- function(a, b) {
  tab <- table(a, b)
  c2 <- function(x) x * (x - 1) / 2
  sn <- sum(c2(tab)); sa <- sum(c2(rowSums(tab))); sb <- sum(c2(colSums(tab)))
  expected <- sa * sb / c2(length(a))
  denom <- 0.5 * (sa + sb) - expected
  if (denom == 0) return(1)
  (sn - expected) / denom
}
measure <- function(dir, label) {
  b <- read_network_bundle(dir)
  p <- b$params
  prev <- stats::setNames(b$nodes$prevalence, b$nodes$item)
  net <- cooccurrence_network_from_edge_data(as.data.frame(b$edges), prev, b$summary$n_records,
                                             p$weight_metric, p$min_count, p$alpha)
  julia <- stats::setNames(b$nodes$community, b$nodes$item)[net$items]
  for (m in c("leiden", "louvain")) {
    fit <- detect_communities(net, method = m, seed = 2026L)
    cat(sprintf("%-34s %-7s ARI=%.3f  Q_julia=%.4f  Q_r=%.4f  K_julia=%d  K_r=%d\n", label, m,
                ari(julia, fit$assignments), b$summary$modularity, fit$modularity,
                b$summary$n_communities, fit$n_communities))
  }
}
for (onto in c("summary", "intermediate")) {
  measure(sprintf("data/derived/port_baselines/network_v2_%s/bundle", onto),
          sprintf("network_v2_%s", onto))
  for (g in c("Female", "Male")) {
    measure(sprintf("data/derived/port_baselines/network_stratified_v2_%s/bundle/%s", onto, g),
            sprintf("stratified_v2_%s/%s", onto, g))
  }
}
```

Expected: 12 lines (6 networks × 2 methods). If a baseline directory is missing, stop and tell Brian which `contract/freeze_*_baseline.py` needs running. Don't run it yourself.

- [ ] **Step 3: Brian decides (stop and ask).** Show him the 12 lines and ask for two decisions:
  1. **Method for CancerClustering.** Recommend `"leiden"`, the package default, unless the table shows Louvain matching Julia where Leiden doesn't and Brian prefers continuity.
  2. **Outcome.** **A:** the chosen method has ARI = 1.000 for all six networks, so the partition checks against Julia stay as they are. **B:** some network has ARI < 1 and Q_r ≥ Q_julia everywhere, so those partition checks become a modularity floor (Brian signs off, as the `port_tolerances.yaml` header allows). **C:** some Q_r < Q_julia. Stop: full Leiden/Louvain scoring below phase 1 needs investigating before anything changes.

  Record his answers in the commit message of Step 9.

- [ ] **Step 4: The wrapper takes and passes `seed`.** In `rpkg/SEERClustering/R/sex_anatomy.R`:
  - change `community_method = "louvain") {` to `community_method = "leiden", seed) {`, using Brian's method if he chose differently;
  - after the `if (missing(ontology)) { ... }` block, add:

```r
  if (missing(seed)) {
    stop(paste0("seer_stratified_network_analysis() needs seed =: community detection is ",
                "randomized; pass a whole number (the analysis uses 2026)"), call. = FALSE)
  }
```

  - change `community_method = community_method)` to `community_method = community_method, seed = seed)`;
  - in the roxygen block, add `` #' `seed` goes to `CooccurrenceAnalysis::stratified_network_analysis()` unchanged (every sex uses it); the analysis scripts pass 2026 (AGENTS.md rule 9). ``

  Then find the tests that call it and add `seed = 2026L` to each:

Run: `grep -rn "seer_stratified_network_analysis(\|detect_communities(\|stratified_network_analysis(" rpkg/SEERClustering/tests rpkg/SEERClustering/R scripts`
Every call in `R/` and `scripts/` that reaches `detect_communities` must pass a seed. Add one test to the wrapper's test file:

```r
test_that("seer_stratified_network_analysis requires seed", {
  ev <- make_seer_fixture()  # use the fixture this test file already builds
  expect_error(seer_stratified_network_analysis(ev, ontology = "summary"), "needs seed")
})
```

Use the existing fixture constructor's actual name from that file, not `make_seer_fixture`. The comment marks it.

- [ ] **Step 5: The scripts use the chosen method and record the seed.**
  - `scripts/fit_network.R:48`: `comm <- detect_communities(net, method = "leiden", seed = 2026L)`.
  - Replace lines 60–64 (the `-louvain.json` writer and its comment) with:

```r
  # Community detection is randomized (igraph Leiden); the method and seed
  # are recorded beside the bundle (AGENTS.md rule 9).
  jsonlite::write_json(list(method = "leiden", seed = 2026L),
                       paste0(output, "-communities.json"), auto_unbox = TRUE, pretty = TRUE)
```

  - `scripts/fit_stratified_network.R:53-56`: pass `community_method = "leiden", seed = 2026L` in the `seer_stratified_network_analysis(...)` call. Replace lines 73–74 (the ties JSON) with the same `-communities.json` writer as above, and update the header comment on line 11 from `-louvain.json` to `-communities.json`.
  - In `port.smk` and `Snakefile`, rename every output and path ending `-louvain.json` to `-communities.json` (rule 6). Check: `grep -rn "louvain.json" port.smk Snakefile scripts contract` returns nothing.

- [ ] **Step 6: Contract: drop the tie diagnostic, record the method.** In `contract/port_tolerances.yaml`, replace every `louvain_gain_ties` entry (kind `r_diagnostic`, `keys: [n_gain_ties]`) and its comment block with:

```yaml
  # Not a comparison: the community method and seed R used (igraph Leiden
  # since CooccurrenceAnalysis phase 3; Julia runs phase-1-only Louvain).
  - name: community_method
    kind: r_diagnostic
    r: <same directory as before>/<step>-communities.json
    keys: [method, seed]
```

  For each entry, `r:` is the old path with `-louvain.json` changed to `-communities.json`. In the file header, change `Louvain identical partition (min_ari: 1.0) or Brian's sign-off` to `community partition identical to Julia's (min_ari: 1.0), or under Brian's sign-off a modularity floor (R's full Leiden must score at least Julia's phase-1 Louvain)`.

- [ ] **Step 7 (outcome B only): Modularity floor for the networks Brian signed off.**
  1. Add a failing test to `contract/test_compare_port.py`, importing the new function alongside the file's existing imports from `compare_port`:

```python
def test_modularity_floor():
    j = pd.DataFrame({"n_records": [10], "modularity": [0.70]})
    assert compare_modularity_floor(j, j.assign(modularity=[0.77]))["passed"]
    assert compare_modularity_floor(j, j.assign(modularity=[0.70]))["passed"]
    assert not compare_modularity_floor(j, j.assign(modularity=[0.69]))["passed"]
```

  2. Run `pixi run test-contract` and confirm it fails on the import.
  3. Add to `contract/compare_port.py`, after `compare_partitions`:

```python
def compare_modularity_floor(j: pd.DataFrame, r: pd.DataFrame, column="modularity",
                             abs_tol=1e-9) -> dict:
    """R's partition must score at least Julia's modularity. CooccurrenceAnalysis
    (R) runs full Leiden/Louvain from igraph, which aggregates where the Julia
    reference's phase-1-only Louvain stops, so partitions may differ but Q may
    not fall."""
    qj, qr = float(j[column].iloc[0]), float(r[column].iloc[0])
    return {"passed": qr >= qj - abs_tol, "detail": {"julia": qj, "r": qr}}
```

     In `_check`, before the final `raise`, add:

```python
    if spec["kind"] == "modularity_floor":
        return compare_modularity_floor(j, r, spec.get("column", "modularity"),
                                        spec.get("abs_tol", 1e-9))
```

  4. Run `pixi run test-contract` and confirm it passes.
  5. In `port_tolerances.yaml`, for **each network whose ARI < 1** (the Female/Male entries are separate):
     - replace its `kind: partition` entry with `- name: modularity_floor`, `kind: modularity_floor`, and `julia:`/`r:` pointing at that network's `summary.arrow` (the paths of its `summary` table entry);
     - delete its `summary` table entry, since `modularity` and `n_communities` there are partition-dependent;
     - add `columns:` to its `nodes` and `nodes_csv` table entries listing every column except `community` (read the names with `head -1` on the CSV, or `pyarrow.feather.read_table(path).column_names`; headers carry no data);
     - restrict its `scalars` entry from the glob to an explicit `keys:` list that leaves out `*_modularity` and `*_n_communities`;
     - for its figure-data entries (`fig_*`), add `columns:` leaving out any community or node-color column.

     Each edit gets a one-line YAML comment: `# partition-dependent; floor signed off by Brian <date>`.

- [ ] **Step 8: Notes for the manuscript (no prose edits).** The manuscript currently reports Julia's numbers, so nothing published changes until the R pipeline replaces the Julia rules. Add a dated entry to `PORTING.md` under "Notes":

```markdown
- 2026-10-06: CooccurrenceAnalysis (R) now detects communities with igraph Leiden (seed 2026), not the Julia reference's phase-1-only Louvain. Outcome <A|B> of the phase-3 measurement. When the R network rules replace the Julia ones: revise the methods text that names Louvain (manuscript/_manuscript/index.qmd:102; supplement.qmd:65 and :86, "implemented directly") -- Brian's prose, his edit -- and, under outcome B, expect network_v2_*_modularity / *_n_communities in manuscript/_variables.yml and the sex-comparison ARI to change (they flow through merge_results; never hand-edit them, AGENTS.md rule 4).
```

- [ ] **Step 9: Verify and commit**

```bash
pixi run test-r-cooccurrence
pixi run test-r-seer
pixi run test-contract
```

All pass. Then run the port checks for both network steps at both vocabularies, as AGENTS.md "Commands" describes (rules `check_network_v2` and `check_network_stratified_v2`). All PASS. A failing check outside the entries edited in Steps 6–7 is a real regression: stop and report it.

```bash
git add rpkg scripts port.smk Snakefile contract PORTING.md
git commit -m "Follow CooccurrenceAnalysis to igraph Leiden communities (seed 2026)"
```

The message body records Brian's Step 3 decisions and the 12-line measurement (aggregates only).

---

### Task 5: Whole-branch verification (both repos)

- [ ] **Step 1: Fresh full run.** `Rscript -e 'testthat::test_local(".", reporter = "check")'` gives `FAIL 0 | WARN 0 | SKIP 0`. A skip means `cooccur` or `arules` is missing; report it rather than ignoring it.
- [ ] **Step 2: Diff audit.** Run `git diff main --stat` and `git diff main -- tests/`. Confirm that the only deleted tests are the two named in Task 2 Step 6, and that no existing `expect_*` assertion was loosened. Adding `seed =` and the `n_gain_ties` setdiff removal are the only allowed edits to surviving tests.
- [ ] **Step 3: CancerClustering against this branch.** With `refactor/3-leiden` checked out here, run `pixi run test` in `~/Code/Julia/CancerClustering` on `port/leiden-communities`. All pass.
- [ ] **Step 4: Report to Brian:**
  - the test counts in both repos;
  - the two deleted tests;
  - his Task 4 decisions and the measurement;
  - the Julia issue draft awaiting his decision.

  Merge both branches to their `main`s **together**, and only with his go-ahead: CancerClustering first, then this repo, back to back, so neither `main` ever runs against an incompatible sibling.
