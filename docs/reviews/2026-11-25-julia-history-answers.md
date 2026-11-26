# Answers from the Julia `CooccurrenceAnalysis` history

Written 2026-11-25 from `/Users/brainchapman/Code/Julia/cooccurrence-analysis` at `8865bc9`,
and from its parent analysis project (abbreviated CC below; its name and the application
domain are redacted here). No code was changed in either repository.

**Where the history lives.** This package was split off from CC's old `clustering/`
directory. Its first commit, `3b918d4` (2026-07-09), is "a byte-faithful rename of the old
`clustering/` tree" (`CC/notes/post_migration_plan.md`). The reason recorded for the split,
paraphrased to keep the domain out: the parent project's domain vocabulary interfered with
AI-assisted refactoring, so the algorithms needed a domain-neutral repository, and the split
also keeps the parent project's source data out of a reusable package. CC's git history (323 commits, from 2026-03-17) covers the period before
the split; its `clustering/` code was removed in `7400ed8` (2026-07-29). The commit numbers in
section 3 of the brief come from this repository.

## 1. Where the HDP CAVI update scheme came from

**The design was intended as standard, not new.** `CC/notes/hdp_mixture_plan.md`, from
`26abf33` (2026-04-29), says: "Wang & Blei 2011 is the reference … We're not doing research on
the algorithm, just porting a known-good method." It also cites Beall et al. 2024 (*Stat Med*)
and Teh et al. 2006.

**Both steps trace back to that plan and its implementation:**
- **β plug-in:** the plan writes the plug-in directly into the update (`:295-297`,
  γ' = α·E_q[β_k] + Σr). It is presented as the standard form, not as an approximation.
- **Pooled soft counts:** the plan said to "Pass auxiliary counts from per-group sticks.
  Standard derivation, see Wang & Blei 2011 §3" (`:304`). The code in `26abf33` pools raw
  soft counts instead, with only the comment "treats all-group assignments as exchangeable".

**Falling ELBO:** the plan says the ELBO "must monotonically increase … if it ever decreases,
there's a bug" (`:455`). `26abf33`'s commit message claims "ELBO monotone verified", but its
tests check only `isfinite(result.elbo)`. An adversarial review flagged the soft-count
pooling and the non-monotone ELBO (`CC/notes/review_revisions_plan.md:26-30`, 2026-07-07). It
decided to "*document* it, not re-engineer the sampler". **The plug-in was never flagged.**
Liang et al. 2007 is cited nowhere.

## 2. Louvain phase 1: deliberate or an oversight?

**It was deliberate from the start.** The first version, CC `300ca01` (2026-03-18, Claude Opus
4.6 co-author), already says "Simple Louvain community detection (phase 1 only — sufficient for
small graphs)." No reasoning for "sufficient" is recorded.

It survived unchanged into this package (`src/network_metrics.jl:180`). The 2026-07-09 review
fixed only the gain bug (C1/C2, `314f1bb`).

**Comparison:** yes, after the R port. CC `fc37f99` (2026-10-07) measured R's igraph Leiden
against the frozen Julia results:
- **Summary vocabulary:** ARI 1.000, so the partitions match.
- **Intermediate vocabulary:** ARI 0.80-0.93, and Leiden has higher modularity in every case
  (0.7069 vs 0.6860 for the cohort, 0.6773 vs 0.6631 Female, 0.7008 vs 0.6878 Male).

Brian chose Leiden on 2026-10-06. CC `f589c10` changed the manuscript methods to Leiden, cited
Traag et al. 2019, and disclosed the earlier phase-1 Louvain (`CC/manuscript/supplement.qmd:65`).

## 3. Were the four R-found bugs known?

1. **Chi-square on a zero margin returns NaN.** Not known. `test_association`
   (`src/statistical_validation.jl:108-110`) passes the table to `ChisqTest` with no margin
   guard. The only `:chisq` test (`test/runtests.jl:299`) uses non-degenerate items. `:fisher`
   is the default, so the bug only reaches callers who choose `:chisq`.
2. **Louvain phase 1 only.** Known as a limitation since 2026-03-18 (see question 2).
3. **`min_count` round-trip.** Not known. `src/mining.jl:31-32,72-73` converts with
   `min_count / nrow`. The test (`runtests.jl:232-240`) uses counts of 10 and 2 against a
   threshold of 5, so it cannot catch an off-by-one.
4. **ARI of fewer than 2 elements.** Partly known. C11 (`1c32ac9`) made `compare_networks`
   return `missing` when fewer than 2 items are shared (`network_group_stratification.jl:188-191`),
   but `_adjusted_rand_index` itself still has `n < 2 && return 0.0` (`:115`).

No TODO, FIXME, issue or review in either repository mentions bugs 1, 3 or 4.

## 4. Why one-sided Fisher, and why the BH family of co-occurring pairs?

- **One-sided Fisher:** `src/statistical_validation.jl:107`,
  `tail=:right  # One-sided: testing for positive association`. That comment is the only
  reasoning recorded.
- **The BH family:** `network_construction.jl:109-110` ("Skip pairs with zero co-occurrence")
  removes pairs with n11 = 0 before BH. The rule was already in the pristine tree as
  `observed == 0 && continue`, and `8795423` kept it. No reason is recorded.
- **Related record:** a docstring on the ARM path (`statistical_validation.jl:183-192`) calls
  the BH family there selection-conditioned and "optimistic". It says the manuscript's
  confirmatory test rests on `known_pair_recovery.jl`, which is not in this repository.

I did not search CC's git history for when these two choices were introduced.

## 5. The tie rule in `round_digits`

The Julia package has no `round_digits`: `git log --all -S round_digits` finds nothing. Every
rounding call is Base `round(x, digits=d)`, used for display. Base rounding scales, rounds ties
to even, then unscales. On Julia 1.12.4, `round(0.15, digits=1)` is 0.2, `round(0.25, digits=1)`
is 0.2 and `round(0.35, digits=1)` is 0.4.

The R rule therefore reproduces Julia's default `round`. No record shows it was chosen to match
an earlier table. I did not trace which of CC's tables depend on it.

## 6. Which published results came from the Julia package, and from which commit?

**Nothing has been published in a journal.** CC's manuscript is "Draft for circulation — not yet
a journal manuscript" (`CC/manuscript/index.qmd:21`). It began as `collaborator_draft.md`,
shared with the domain collaborator. "Published" in `d7e92ff` and in the CC labbook ("backing a
published retraction", `labbook/.capture/2026-07-30.facts.md:229`) means that circulated draft.

**Which commits:**
- **v1 results** (deleted in CC `4fa4ca1`, 2026-10-06): produced before 2026-07-29 by CC's
  in-tree `clustering/`, which used the buggy Louvain gain.
- **The frozen Julia v2 baselines the R port checks against:** their `fit.meta.json` records only
  CC's commit (for example `b26998a`), not the Julia package's commit. The network fit finished
  at 16:26 CDT on 2026-09-24 (`created_utc` 21:26Z), before `3f0b229`. I infer from that timing
  that the Julia code was at `d7e92ff`; this is not recorded. `src/` changed only in comments
  after that commit.

**Effect of the fixes:** Louvain was already handled in the manuscript methods (question 2).
One §4.2 figure "needs re-running" (CC labbook, same line).

## 7. References in the Julia code and docs, compared with brief §5

**Cited in the Julia package:** Teh et al. 2006, Wang, Paisley & Blei 2011, Sethuraman 1994,
"Sudderth 2006", and Ye, Zhang & Nie 2018 (*CSDA* 123:32-49, `bayesian_clustering.jl:83,128`).
CC's HDP plan also cites Beall et al. 2024 (*Stat Med*).

**Problems:**

- **Wrong:** `hdp_clustering.jl:490` cites "Wang & Blei 2011". This is the original April
  citation (CC `26abf33:3`, and the plan). Paisley was added to the header at `:5`, but the
  docstring at `:490` was not updated. CC `notes/refactoring_plan.md:190` records the Wang,
  Paisley & Blei paper being added to `references.bib`.
- **Cited for a different method:** Wang, Paisley & Blei is credited for an algorithm the code
  says it departs from. "Sudderth 2006" has no full reference.
- **Missing from brief §5:**
  - Ye, Zhang & Nie 2018 (I did not check this citation)
  - Beall et al. 2024
  - FP-Growth (Han, Pei & Yin 2000; `mining.jl:36`)
  - label propagation (Raghavan et al. 2007), which the CC supplement (`:67`) still names
  - BIC (Schwarz 1978)
  - EM (Dempster et al. 1977)
  - NUTS and ADVI
- **Uncited:** Newman, Haldane-Anscombe, BH and Holm are named in the Julia code but given no
  citation.
