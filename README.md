# CooccurrenceAnalysis (R)

Domain-neutral analysis of co-occurrence in binary (presence/absence) event data:

- pairwise association statistics: one-sided Fisher exact or Pearson chi-square, with
  Benjamini-Hochberg, Bonferroni or Holm correction;
- significance-filtered co-occurrence networks, with Leiden (default) or Louvain communities
  from igraph and per-node metrics;
- group-stratified networks, with items that can occur in only one group excluded from the
  others;
- a hierarchical Dirichlet process (HDP) Bernoulli mixture fit jointly across groups by
  truncated stick-breaking coordinate-ascent variational inference, with multi-seed fitting
  and stability summaries;
- frequent itemsets and association rules, enumerated by `arules`;
- ggplot2 plots of all of the above.

Every method here was published by others (see [References](#references)); the package makes
no statistical innovation of its own. One part is not a published algorithm as written: the
update scheme of the HDP mixture's variational inference, which adapts Wang, Paisley & Blei
(2011) with two approximations (see [Known limitations](#known-limitations)).

The data model has no domain vocabulary. **Records** (`id`) accumulate **events**, each
carrying one **item**, and records may belong to **groups**. Plot titles and axis labels are
arguments, so callers supply their own nouns.

## Status

Version 0.1.0, not on CRAN, and no rendered documentation site exists. The documentation is
the roxygen comments in `R/`; there is no `man/` directory, so `?function` help pages are not
installed. Expect breaking changes before 1.0: two were made on 2026-11-25
(`cooccurrence_network_from_edge_data()` now requires `test` and `correction`, and the
`edge_colour` column returned by `plot_network_comparison()` is now `edge_color`).

## Installation

```r
# install.packages("remotes")
remotes::install_github("chapmanbe/cooccurrence-analysis-r")
install.packages("arules")  # needed for itemsets and rules
```

## Example

Functions take an event table with one row per (record, item) event. A record holding an
item twice counts once, and records with fewer than `min_items` distinct items (default 2)
are excluded.

```r
library(CooccurrenceAnalysis)

events <- data.frame(
  id   = rep(1:40, each = 2),
  item = c(rep(c("a", "b"), 15), rep(c("c", "d"), 15), rep(c("a", "c"), 10))
)
pairs <- compute_pairwise_associations(events)            # every co-occurring pair
net <- build_cooccurrence_network(events, min_count = 5L)
comm <- detect_communities(net, seed = 1L)               # seed is required
plot_cooccurrence_network(net, communities = comm)
```

`seed` has no default anywhere a result is randomized: pass a whole number for a
reproducible result, or `NULL` to use the caller's random-number stream.

## How this package was developed

### Origin: a port from Julia

The package began on 2026-09-25 as an R port of the Julia package `CooccurrenceAnalysis`.
That package was itself split, on 2026-07-09, from the analysis project where its code had
been written since March 2026, so that the algorithms could live in a reusable,
domain-neutral package. The port went in steps (pairwise statistics and
networks, group stratification, the HDP mixture), each written test first. Early tests took
their expected values from running the Julia code.

A refactor (2026-10-06/07, `docs/specs/2026-10-06-r-ecosystem-refactor.md`) then replaced
hand-written algorithms with established R implementations where one exists: igraph for
community detection and modularity, `stats` for the test statistics and p-value correction,
and `arules` for itemset enumeration. The HDP mixture and the group-stratified network have
no R equivalent and remain this package's own code.

### The package is now its own specification

The port is finished, and Julia is no longer treated as the reference. The Julia code was
never reviewed the way this package has been, and porting it turned up four of its bugs,
each fixed here and each tested by a test marked `R-only`:

- Chi-square on a zero margin returns `NaN` in Julia, which turns every q-value in the family
  into `NaN`. Here it raises.
- Julia's Louvain runs only the first phase (local moving, no aggregation) and stalls below
  the optimum: Q = 0.70 against 0.77 on a ring of 20 triangles. Here communities come from
  igraph's full Leiden or Louvain, so partitions differ from Julia's.
- Julia's mining converts `min_count` to a fraction and back, which can land one record
  above it (7 of 100 gives 8). Here it is used as a count.
- The adjusted Rand index of fewer than two elements is 0 in Julia and 1 here, as the
  definition gives.

Some of Julia's choices are kept deliberately, because results already computed with the
Julia package, including a manuscript's, used them: `round_digits()` rounds as Julia's
`round` does, Fisher tests are one-sided, the BH family is exactly the pairs that co-occur,
and odds ratios get the Haldane correction. The Julia history records no reasoning for the
one-sided test or the family beyond the code itself.

Not ported: the flat-K Bernoulli mixtures (deliberately, since the HDP covers them), timing
filters, `label_propagation` and `compare_networks`. Calling any of them raises an error
naming the Julia function.

### Use of AI coding assistants

Brian Chapman directed the work and made the design and modeling decisions. Most of the code,
tests and documentation were written with AI coding assistants in Claude Code: 36 of the 39
commits through 2026-11-25 carry a `Co-Authored-By` line for Claude (Opus 5.5 or Sonnet 5.5).
The Julia package was written the same way: 26 of its 28 commits carry a Claude co-author
line, as does the first version of its network code in the parent project.
The working rules the assistants followed are in `CLAUDE.md`. They include writing tests
first, raising on missing or invalid input rather than dropping it, and requiring an
independent source for every expected value.

## How it is tested

```r
testthat::test_local(".")   # the suite builds its own synthetic fixtures; no data needed
lintr::lint_package()
```

As of 2026-11-25 the suite has 169 tests (858 expectations), and all pass with no warnings.
Partial-match warnings are turned on during tests, and no test suppresses a warning.

**Where expected values come from.** A test that checks a number must take it from a source
independent of this package, in this order of preference:

1. a fixture with one right answer by construction, or a value computed by hand from the
   definition;
2. an established implementation: `stats::fisher.test`, `stats::chisq.test`,
   `stats::p.adjust`, `cooccur` (Veech 2013), `flexmix`, igraph or `arules`;
3. the Julia package, only alongside one of the above.

An audit on 2026-11-25 (`docs/reviews/2026-11-25-oracle-audit.md`) found 10 of the then 150
tests taking their values from Julia alone. All but the HDP known-answer traces now also check
against independent code in `tests/testthat/helper-fixtures.R`. For the HDP, an evaluator
written separately from the model's definition (`tests/testthat/helper-hdp-elbo.R`)
recomputes the evidence lower bound (ELBO) at every iteration of the four known-answer fits
and agrees with both this package and Julia. The HDP's recovery of planted clusters is also
compared against `flexmix`.

**Checking that tests can fail.** Each test added or strengthened on 2026-11-25 was confirmed
by deliberately breaking the code and watching it fail. Examples: chi-square computed with N - 1, lift with
n + 1, Holm swapped for Hochberg, a stick-breaking prior changed by 1%, a hardcoded seed, an
ignored `resolution`, and record counts off by one.

The `cooccur` cross-check needs the GitHub version (`remotes::install_github("griffithdan/cooccur")`)
and is skipped without it; `flexmix` and `arules` tests skip if those packages are missing.

## Code review

Two review passes on 2026-09-25 led to fixes in commits `6a09b00` and `f849c40`. The code
written after them (about 4,800 lines across 35 files) had a two-model review on 2026-11-25,
recorded in `docs/reviews/2026-11-25-code-review.md`:

- **Claude** reviewed the diff with Claude Code's `/code-review` and found six bugs and four
  cleanups.
- **Gemini**, run separately through Google's Antigravity CLI, received the project rules and
  the full diff and returned 20 findings. It found all six of Claude's bugs independently and
  added ten more problems and three tests that could not fail.
- A third reviewer (OpenAI Codex) was planned but did not run.

Each finding was verified before being fixed, by reproducing it in R or by reading the code.
Every bug fix began with a test that failed, and the merged list records which model raised
each finding. All 23 findings (6 shared bugs, 10 Gemini-only bugs, 3 tests that could not
fail, and 4 cleanups) were fixed by 2026-11-25.

## Known limitations

- **The HDP updates are an approximation.** The model is the standard truncated
  stick-breaking HDP (Teh et al. 2006; Sethuraman 1994), with empirical-Bayes Beta priors on
  the items following Ye, Zhang & Nie (2018). The Julia code describes its algorithm as "in
  the spirit of" Wang, Paisley & Blei (2011, §3.1) batch CAVI, but departs from it on
  purpose: it uses a single-level shared-atom representation with no table indicators, so
  there is no q(c_jt). Two update steps are not exact coordinate ascent:
  - the global-stick update pools raw soft counts across groups rather than table counts;
  - the group-stick prior, Beta(α·b_k, α·Σ_{l>k} b_l), uses the global weights b at their
    variational mean in place of an expectation with no closed form.

  The ELBO formula itself is checked independently, but because of these two steps the ELBO
  can fall between iterations, and it does in the known-answer traces. Neither step was
  meant as new method. The design plan (2026-04-29, in the parent project) set out to port
  Wang & Blei's published algorithm ("not doing research on the algorithm, just porting a
  known-good method"). It wrote the plug-in into the update as the standard form, and it
  called for auxiliary counts in the global-stick update; the implementation pooled raw soft
  counts instead. The plan required the ELBO to rise monotonically, but the first tests checked
  only that it was finite. A review in July 2026 flagged the pooling and the falling ELBO and
  chose to document it rather than change the updates. The plug-in was never flagged, and the
  Julia comment calls that math exact. The scheme has not been compared with another
  published truncated HDP scheme, such as Liang et al. (2007), nor with another HDP
  implementation.
  Treat HDP results as approximate posterior summaries, and check their stability across
  seeds (`hdp_fit_seeds()`).
- **The HDP known-answer traces still rest on Julia** for the full sequence of updates; only
  their ELBO values are checked independently.
- `round_digits()` deliberately reproduces Julia's built-in `round`, which the Julia package
  used for every table and label (`round_digits(0.15, 1)` is 0.2, although 0.15 is stored as
  0.1499…). Whether any published table depends on this is not recorded.
- No `R CMD check` result is recorded here.

## References

Checked on 2026-11-25 against Crossref by DOI, or against the publisher's page where there is
no DOI. Holm (1979) is confirmed only from secondary sources.

**Association statistics and multiple testing**
- Pearson, K. (1900). On the criterion that a given system of deviations from the probable in
  the case of a correlated system of variables is such that it can be reasonably supposed to
  have arisen from random sampling. *Philosophical Magazine* 50(302), 157–175.
  https://doi.org/10.1080/14786440009463897
- Fisher, R. A. (1935). The logic of inductive inference. *Journal of the Royal Statistical
  Society* 98(1), 39–54. https://doi.org/10.2307/2342435
- Haldane, J. B. S. (1956). The estimation and significance of the logarithm of a ratio of
  frequencies. *Annals of Human Genetics* 20(4), 309–311.
  https://doi.org/10.1111/j.1469-1809.1955.tb01285.x
- Anscombe, F. J. (1956). On estimating binomial response relations. *Biometrika* 43(3–4),
  461–464. https://doi.org/10.1093/biomet/43.3-4.461
- Holm, S. (1979). A simple sequentially rejective multiple test procedure. *Scandinavian
  Journal of Statistics* 6(2), 65–70.
- Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and
  powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B*
  57(1), 289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

**Networks and communities**
- Hidalgo, C. A., Blumm, N., Barabási, A.-L., & Christakis, N. A. (2009). A dynamic network
  approach for the study of human phenotypes. *PLoS Computational Biology* 5(4), e1000353.
  https://doi.org/10.1371/journal.pcbi.1000353
- Newman, M. E. J., & Girvan, M. (2004). Finding and evaluating community structure in
  networks. *Physical Review E* 69(2), 026113. https://doi.org/10.1103/PhysRevE.69.026113
- Blondel, V. D., Guillaume, J.-L., Lambiotte, R., & Lefebvre, E. (2008). Fast unfolding of
  communities in large networks. *Journal of Statistical Mechanics* 2008(10), P10008.
  https://doi.org/10.1088/1742-5468/2008/10/P10008
- Traag, V. A., Waltman, L., & van Eck, N. J. (2019). From Louvain to Leiden: guaranteeing
  well-connected communities. *Scientific Reports* 9, 5233.
  https://doi.org/10.1038/s41598-019-41695-z

**Mixture model and clustering comparison**
- Sethuraman, J. (1994). A constructive definition of Dirichlet priors. *Statistica Sinica*
  4(2), 639–650.
- Teh, Y. W., Jordan, M. I., Beal, M. J., & Blei, D. M. (2006). Hierarchical Dirichlet
  processes. *Journal of the American Statistical Association* 101(476), 1566–1581.
  https://doi.org/10.1198/016214506000000302
- Blei, D. M., & Jordan, M. I. (2006). Variational inference for Dirichlet process mixtures.
  *Bayesian Analysis* 1(1), 121–143. https://doi.org/10.1214/06-BA104
- Liang, P., Petrov, S., Jordan, M., & Klein, D. (2007). The infinite PCFG using hierarchical
  Dirichlet processes. *Proceedings of EMNLP-CoNLL 2007*, 688–697.
  https://aclanthology.org/D07-1072/
- Wang, C., Paisley, J., & Blei, D. M. (2011). Online variational inference for the
  hierarchical Dirichlet process. *Proceedings of AISTATS 2011*, PMLR 15, 752–760.
  https://proceedings.mlr.press/v15/wang11a.html
- Ye, M., Zhang, P., & Nie, L. (2018). Clustering sparse binary data with hierarchical
  Bayesian Bernoulli mixture model. *Computational Statistics & Data Analysis* 123, 32–49.
  https://doi.org/10.1016/j.csda.2018.01.020
- Hubert, L., & Arabie, P. (1985). Comparing partitions. *Journal of Classification* 2(1),
  193–218. https://doi.org/10.1007/BF01908075

**Itemsets and association rules**
- Agrawal, R., & Srikant, R. (1994). Fast algorithms for mining association rules.
  *Proceedings of VLDB 1994*, 487–499.
- Brin, S., Motwani, R., Ullman, J. D., & Tsur, S. (1997). Dynamic itemset counting and
  implication rules for market basket data. *Proceedings of SIGMOD '97*, 255–264.
  https://doi.org/10.1145/253260.253325

**R packages used or checked against**
- Csárdi, G., & Nepusz, T. igraph (R package). https://igraph.org
- Hahsler, M., Grün, B., & Hornik, K. (2005). arules: a computational environment for mining
  association rules and frequent item sets. *Journal of Statistical Software* 14(15).
  https://doi.org/10.18637/jss.v014.i15
- Leisch, F. (2004). FlexMix: a general framework for finite mixture models and latent class
  regression in R. *Journal of Statistical Software* 11(8).
  https://doi.org/10.18637/jss.v011.i08
- Griffith, D. M., Veech, J. A., & Marsh, C. J. (2016). cooccur: probabilistic species
  co-occurrence analysis in R. *Journal of Statistical Software* 69(Code Snippet 2).
  https://doi.org/10.18637/jss.v069.c02
- Veech, J. A. (2013). A probabilistic model for analysing species co-occurrence. *Global
  Ecology and Biogeography* 22(2), 252–260. https://doi.org/10.1111/j.1466-8238.2012.00789.x

## License

MIT; see `LICENSE.md`.
