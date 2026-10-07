# Knowledge Digest: R-ecosystem refactor, what is done and what is left
*Written 2026-10-07 at the close of the session that ran Phases 0-3. Read this once; work from it.*

> **Kickoff:** Read this file, then (1) decide D2 (flat-K Bernoulli mixture: port, wrap, or skip) and, if D2 says
> anything but "skip," do R9; (2) write the R12 sanity test (HDP vs a flat-K fit); (3) before filing the Julia
> issue, run the *Julia* reference on the ring-of-20-triangles fixture; (4) do not keep the sentence "Leiden
> guarantees connected communities" in the docs unless you have checked it against the Leiden paper; it is
> [ASSERTED] below. Do not merge, push, or file anything without the user's go-ahead.

## Objective
Finish the refactor in `docs/specs/2026-10-06-r-ecosystem-refactor.md` so this package delegates to validated R
packages wherever they compute the same quantity. Phases 0-3 (items R1-R4, R6, R8) are done and merged to local
`main`. Remaining: R9, R12, the Julia issue, the push, and the manuscript prose. "Done" = the spec's Phase 3 box
ticked for R9/R12, both repos pushed (CancerClustering first), and the Julia issue filed or deliberately dropped.

## Central argument
The Julia package is the reference only where it is right. Where an R package gives the same quantity under this
package's conventions, delegate; where Julia is wrong (phase-1-only Louvain, NaN chi-square, `min_count` off by
one), fix it in R, record the divergence in CLAUDE.md, and mark the paired test `R-only`. Never loosen a test to
make a delegation pass.

## Decisions & framing (do not relitigate)
- **Julia's phase-1-only Louvain is a bug**, not a spec (user, 2026-10-06). Replaced by igraph.
- **Default community method is Leiden** (`method = "leiden"`, modularity objective, `n_iterations = -1`);
  `"louvain"` is the option. igraph's Leiden defaults (CPM objective, 2 iterations) are wrong for this package.
- **`seed` has no default** in `detect_communities`, `stratified_network_analysis`, `run_network_pipeline`
  (`NULL` = draw from the caller's stream). `.with_seed()` pins Mersenne-Twister/Inversion/Rejection so a seed
  reproduces whatever `RNGkind` the caller set, and restores the caller's stream.
- **Labels are numbered by first appearance in `net$items` order**, because igraph's raw labels are arbitrary.
- **Raw community label ids are not compared exactly against Julia anywhere** in CancerClustering's contract
  (Julia's follow its internal labels and cannot be reproduced); partitions are checked by ARI, or by a modularity
  floor at the `intermediate` vocabulary (user signed off 2026-10-06).
- **CancerClustering runs community detection on the `--seed` it records** (20260101 from its Snakefile).
- **`cooccur` is a test-time cross-check only**, never a runtime dependency (GitHub-only for current R).
- **No `ggraph`**; plot-encoding tests read `ggplot_build()` layers.
- **Merge order across repos: CancerClustering first, then this repo**, because CancerClustering loads this
  repo's working copy (`pkgload::load_all("../cooccurrence-analysis-r")`).

## Structure (what exists now)
`R/`: `network_construction.R` (2x2 stats, chi-square/Fisher, BH/Bonferroni/Holm), `network_metrics.R` (igraph
communities, modularity), `mining.R` (arules-enumerated itemsets/rules, counts computed here), `utils.R` (`.incidence`,
`.with_seed`, `.check_seed`), HDP files untouched. Plans/specs: `docs/specs/`, `docs/superpowers/plans/`.

## Sources
Literature and package pages. Existence was checked by web search this session; nothing below was read in full.
| Source | Locator | Status | Supports |
|---|---|---|---|
| cooccur: Probabilistic Species Co-Occurrence Analysis in R (J Stat Softw) | https://www.jstatsoft.org/v069/c02 | [VERIFIED] page located | model "essentially equivalent to Fisher's exact test" (search snippet) |
| arules (J Stat Softw) | https://jstatsoft.org/v14/i15 | [VERIFIED] page located | arules wraps Borgelt's Apriori and Eclat |
| poLCA (J Stat Softw) | https://jstatsoft.org/v42/i10 | [VERIFIED] page located | EM latent-class models for polytomous data |
| BayesBinMix (arXiv / R Journal) | https://arxiv.org/pdf/1609.06960 | [VERIFIED] page located | Bayesian MCMC mixtures of multivariate Bernoulli, unknown K |
| flexmix `FLXMCmvbinary` | named in a search result, no page opened | [ASSERTED] | the Bernoulli mixture driver in flexmix |
| Traag, Waltman & van Eck, the Leiden algorithm paper | no locator found this session | [ASSERTED] | "Leiden guarantees connected communities; Louvain can return disconnected ones" |

## Key claims
Verified by running something this session (not citations):
- `cooccur` `p_gt` equals this package's one-sided Fisher `p_value` to the 5 places `cooccur` rounds to — high, verified.
- `modularity_q` equals `igraph::modularity` to 1e-10 on weighted graphs — high, verified.
- Julia chi-square on a zero margin returns NaN and turns every BH q-value NaN — high, verified by running Julia.
- Julia `min_count = 7` of 100 becomes threshold 8 via `ceil(7/100*100)` in RuleMiner — high, verified by running `ceil` and reading `core_utils.jl:9`.
- igraph `cluster_leiden(n_iterations = -1)` errors on an edgeless graph; same seed gives same partition — high, verified (igraph 2.3.3 and 2.3.4).
- At the `summary` vocabulary all three CancerClustering networks reproduce Julia's partition (ARI 1.000); at `intermediate` ARI is 0.80-0.93 and R's Q is higher in all three (0.7069/0.6773/0.7008 vs Julia 0.6860/0.6631/0.6878) — high, verified against the frozen baselines.
Not verified:
- "Leiden guarantees connected communities" (in `R/network_metrics.R:39` and spec D1) — [ASSERTED], from memory.
- "Phase 1 stops at Q = 0.70 on a ring of 20 triangles" as a statement about **Julia** (in `docs/julia-issue-full-louvain.md`) — [ASSERTED]: it was measured on this port's old `.louvain()`, which mirrors Julia but is not Julia. Run Julia before filing.
- `flexmix`/`poLCA` fit maximum-likelihood EM while Julia's `fit_bernoulli_mixture` uses empirical-Bayes priors (MAP) — [ASSERTED]; Julia's source was read for function names, flexmix's docs were not.

## Constraints
Per CLAUDE.md: American English, no domain nouns in the package, no `suppressWarnings`/`try`, `vapply` not `sapply`,
`NAMESPACE` hand-maintained, every numerical change needs a test whose expected values come from **running Julia**
(except documented divergences). CancerClustering AGENTS.md: no SEER data leaves the machine; only Snakemake rules
write `results/`; never type a number into the manuscript. Tool note: the shell hook rewrites `head`/`grep` into
summarizers; use `grep -F`, Python, or R to read exact lines.

## Non-goals / rejected framings
- Keeping Julia's partitions for parity — rejected: they come from a bug.
- Porting association-rule enumeration by hand — rejected: `arules` exists (done as a wrapper).
- A runtime dependency on `cooccur` or `ggraph` — rejected (see Decisions).
- Hand-editing `manuscript/_variables.yml` — forbidden (AGENTS.md rule 4); values flow through `merge_results`.

## Open questions for the drafting session
- [ ] **D2:** port Julia's flat-K MAP-EM, wrap `flexmix`, or skip (the HDP subsumes it)? Recommendation was skip unless a consumer needs flat-K output.
- [ ] **R12:** agree tolerances (ARI floor, matched item-probability tolerance) for an HDP-vs-`flexmix` sanity test on K=3 simulated data.
- [ ] Verify the Leiden connectedness claim against the paper, or soften the roxygen and spec wording.
- [ ] Run the ring-of-20-triangles fixture through **Julia** before filing `docs/julia-issue-full-louvain.md`; then ask the user whether to file it.
- [ ] Push `main` in both repos (CancerClustering first); both are ahead of `origin/main` by earlier local commits too. The user decides.
- [ ] When the R network rules replace the Julia ones in CancerClustering: revise the methods text naming Louvain (`manuscript/_manuscript/index.qmd:102`, `supplement.qmd:65` and `:86`), which is the user's prose, and expect `network_v2_intermediate_modularity`, `_n_communities` and the sex-comparison ARI to change.
- [ ] Deferred minors: `run_network_pipeline` should call `.check_seed()` up front; `modularity_floor` is one-sided (a collapse to fewer communities at equal Q would pass).
- [ ] igraph floor `>= 2.3.3` was tested only under CancerClustering's renv; `arules` and `cooccur` are absent from that renv, so mining tests and the cooccur cross-check skip there.
- [ ] Any [ASSERTED] item above still to verify before it is cited or filed.

## Suggested drafting model
Strong model for D2 and the Leiden-claim check (judgment). Sonnet is fine for R12 once tolerances are fixed in
writing, and for filing mechanics.
