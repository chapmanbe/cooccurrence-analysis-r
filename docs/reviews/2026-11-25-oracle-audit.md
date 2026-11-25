# Oracle audit: which tests rest on Julia alone (2026-11-25)

## Why
CLAUDE.md now treats Julia as one cross-check, not the specification. A test whose expected
values come only from Julia passes if Julia and R share a bug. This audit lists those tests.

## Result
Of 150 `test_that` blocks: 10 take their values from Julia alone, 60 from an independent
source, 80 check no values (errors, classes, shapes). No test pairs a Julia value with an
independent one.

An agent recomputed each Julia-only value that has a practical oracle, and every one matched
to 1e-12. So far, no shared bug has turned up in these values.

| Test | Pinned | Independent oracle to add |
|---|---|---|
| `test-validation.R:10` | odds ratios, raw and Haldane | closed form `ad/bc`, `(a+.5)(d+.5)/((b+.5)(c+.5))` |
| `test-validation.R:55` | chi-square p, BH q, edges | `chisq.test(correct = FALSE)` per pair, then `p.adjust("BH")` |
| `test-validation.R:91` | Bonferroni, Holm q | `fisher_greater()` per pair, then `p.adjust` |
| `test-mining.R:11`, `:21`, `:41`, `:51` | itemsets, rules, metrics | hand counts from the 20 records, or `arules::apriori` called directly |
| `test-hdp-multiseed.R:49` | ARI values | the existing `adjusted_rand()` helper |
| `test-network.R:140` | `round_digits` ties | `round(x * 10^d) / 10^d`. The tie rule is a convention, chosen deliberately |
| `test-hdp-known-answers.R:39` | full CAVI traces (4 cases) | none for the whole trace. Partial: an ELBO evaluator written from the definition, and a fixed-point check |

## Open questions the audit raised
1. **The HDP global-stick update is an approximation inherited from Julia.** It pools raw
   soft counts rather than table counts (`R/hdp_clustering.R:119`), so the ELBO is still a
   lower bound but can fall, and the pinned traces do fall. This is documented, not a bug,
   but no one has checked the approximation outside Julia. An independent ELBO evaluator would
   test the ELBO formula. Whether the approximation is acceptable is a modeling decision.
2. **ARI for one record:** `adjusted_rand_index(1, 1)` is pinned to 0 in
   `test-hdp-multiseed.R:66`, while the `adjusted_rand()` helper returns 1. The suite holds
   both conventions, so one should be chosen.
3. **`round_digits` ties:** 0.15 is stored as 0.1499…, so correct rounding gives 0.1. Julia
   and R both give 0.2, by design, to keep published tables unchanged. The only source for this
   choice is Julia.
4. **Chi-square fixture:** expected cells are below 5, so `chisq.test` warns there. The
   values are right, but the case is weak.

## Agent's full table

Counts (150 `test_that` blocks; the known-answer loop is one source block, 4 runtime cases):
JULIA-ONLY 10 · JULIA+INDEPENDENT 0 · INDEPENDENT 60 · NO-VALUES 80.

No test mixes Julia values with an in-test independent check. Method-level independent
coverage exists on *other* fixtures (test-r-ecosystem-crosscheck.R: cooccur, fisher.test,
p.adjust, igraph; test-validation.R:71 chisq.test), but not on the pinned fixture values.

"Verified" = I recomputed the pinned value with the named oracle during this audit (matches to <1e-12).

## JULIA-ONLY

| file | line | test | values pinned | suggested independent oracle |
|---|---|---|---|---|
| test-validation.R | 10 | odds_ratio zero-cell handling (C5) | 33.333..., 121, 1.2857..., 101.667, 67.857, 0.06294 | Closed form in-test: `ad/bc`, Haldane `(a+.5)(d+.5)/((b+.5)(c+.5))` (verified) |
| test-validation.R | 55 | test_association with chisq | 4 p-values, 4 BH q-values, edge set on make_test_event_df | `stats::chisq.test(tab, correct=FALSE)` per pair + `p.adjust("BH")` (verified; note expected cells < 5, chisq.test warns) |
| test-validation.R | 91 | network with Bonferroni and Holm | 4 Bonferroni, 4 Holm q-values | `fisher_greater()` per pair + `p.adjust("bonferroni"/"holm")` (verified) |
| test-mining.R | 11 | mine_frequent_itemsets | 4 itemsets, N 10/5/3/2, support .5/.25/.15/.1 | Hand count from fixture doc (20 multi-item records), or `arules::apriori(target="frequent itemsets")` directly |
| test-mining.R | 21 | thresholds are inclusive | itemset sets at min_count 5/3/2, min_support .25 | Hand count; inclusiveness is a spec choice (>=) |
| test-mining.R | 41 | mine_association_rules | 7 rules; Item01=>GB1 conf 3/13 excluded | Hand count, or `arules::apriori(target="rules")` + `quality()` |
| test-mining.R | 51 | rule metrics match the reference exactly | 8 rules: N, confidence, coverage, lift, support | Closed form from counts: supp=n_ab/20, conf=n_ab/n_a, lift=conf/(n_b/20) (verified) |
| test-hdp-multiseed.R | 49 | adjusted_rand_index matches the Julia reference | 0.4444, 0.2381, 0, 1, 0.7232 | `adjusted_rand()` in helper-fixtures.R (verified) or `mclust::adjustedRandIndex` |
| test-hdp-known-answers.R | 39 | known answer from the reference: 4 cases | full ELBO traces, n_iter, theta, pi_mean, beta_mean, assignments, effective_k | Partial only: (a) independent first-principles ELBO evaluator at returned params; (b) fixed-point check: one more coordinate update leaves params unchanged; (c) ELBO monotonicity under exact CAVI. Full trace: none practical (deterministic but only Julia implements it) |
| test-network.R | 140 | round_digits matches the Julia reference on ties | 0.15->0.2, 1.015->1.01, 2.5->2 | Closed form: `round(x * 10^d) / 10^d` (IEC half-even at 0 digits); the rule itself is a convention, not a truth |

## JULIA+INDEPENDENT

None.

## INDEPENDENT

| file | line | test | source |
|---|---|---|---|
| test-bh-family.R | 24 | only co-occurring pairs are tested | fixture by construction |
| test-bh-family.R | 32 | p_value is one-sided Fisher | fisher_greater() (stats::fisher.test) |
| test-bh-family.R | 45 | p_adjusted is BH over tested pairs | stats::p.adjust |
| test-bh-family.R | 61 | widening the family loses an edge | fixture by construction |
| test-bh-family.R | 71 | edge filter reads q | fixture by construction |
| test-bh-family.R | 95 | depleted pair never an edge | fixture by construction |
| test-communities.R | 53 | two bridged triangles under any seed | unambiguous fixture; modularity_q |
| test-communities.R | 68 | disconnected components | unambiguous fixture |
| test-communities.R | 77 | full Leiden/Louvain beat phase-1 | ring of cliques, thresholds (Julia Q=0.70 cited, not asserted) |
| test-communities.R | 121 | edgeless and empty networks | trivial closed form |
| test-communities.R | 133 | isolated vertex own community | unambiguous fixture |
| test-hdp-flexmix.R | 26 | HDP vs flexmix | flexmix + planted truth ARI |
| test-hdp-multiseed.R | 62 | ARI R-only edge cases | closed form / convention |
| test-hdp-multiseed.R | 71 | solve_lsap brute force | brute-force enumeration |
| test-hdp-multiseed.R | 85 | solve_lsap greedy trap | hand-computed |
| test-hdp-multiseed.R | 103 | solve_lsap rectangular | brute force |
| test-hdp-multiseed.R | 121 | hdp_best_elbo ties | hand |
| test-hdp-multiseed.R | 131 | hdp_best_elbo tolerance | hand |
| test-hdp-multiseed.R | 195 | ARI without overflow | closed form (identity) |
| test-hdp-multiseed.R | 274 | hdp_match_clusters permutation | by construction |
| test-hdp-multiseed.R | 290 | hdp_match_clusters optimal | hand cosines |
| test-hdp-multiseed.R | 301 | hdp_cluster_recurrence | hand |
| test-hdp-visualization.R | 63 | profile item order | hand |
| test-hdp.R | 49 | empirical_bayes_priors | closed form (mean*conc + floor): 6/6, 2/10, 1/11 |
| test-hdp.R | 65 | 3-cluster planted recovery | ARI vs planted truth (threshold) |
| test-hdp.R | 149 | hdp_clustering on event_df | n_records_per_group from fixture arithmetic |
| test-hdp.R | 187 | categorization planted labels | planted |
| test-hdp.R | 212 | C9 regression | planted blocks, thresholds |
| test-hdp.R | 227 | native 3-group categorization | planted |
| test-hdp.R | 264 | stopping rule per-iteration | synthetic trace, hand |
| test-hdp.R | 296 | joint fit shares atoms | planted |
| test-hdp.R | 307 | record_item_matrix | hand |
| test-hdp.R | 326 | Beta KL and stick means | numerical integration / closed form |
| test-hdp.R | 343 | init matched by id | hand |
| test-mining.R | 31 | min_count as a count | by construction |
| test-mining.R | 73 | three-item rules | closed form from counts |
| test-network-stratification.R | 45 | exclusion records | fixture arithmetic |
| test-network-stratification.R | 57 | C-locale order | definition |
| test-network-stratification.R | 70 | NULL = no exclusion | fixture arithmetic (36) |
| test-network-stratification.R | 89 | excluded item dropped | fixture arithmetic (35, 28) |
| test-network-stratification.R | 104 | X-Y known answer | fisher_greater() |
| test-network-stratification.R | 161 | shared fixture | by construction |
| test-network-stratification.R | 192 | 3 groups | by construction |
| test-network-stratification.R | 221 | plot_network_comparison panels | gray70 = #B3B3B3 |
| test-network-stratification.R | 253 | draws as Julia: alpha, colour | Wong (2011) palette, spec values |
| test-network-stratification.R | 268 | wong colours | Wong (2011) published palette |
| test-network.R | 7 | phi_coefficient | closed form |
| test-network.R | 24 | columns, order, rounding | closed form from counts |
| test-network.R | 51 | vertices, weights, prevalence | fixture arithmetic |
| test-network.R | 86 | renumbering, names | unambiguous fixture |
| test-network.R | 116 | metrics on two triangles | hand-computed |
| test-network.R | 133 | top_by_centrality | hand |
| test-network.R | 164 | modularity_q edgeless | convention (vs igraph NaN) |
| test-r-ecosystem-crosscheck.R | 24 | agree with cooccur | cooccur |
| test-r-ecosystem-crosscheck.R | 51 | agree with fisher.test | stats::fisher.test, p.adjust |
| test-r-ecosystem-crosscheck.R | 82 | modularity_q vs igraph | igraph::modularity |
| test-r-ecosystem-crosscheck.R | 96 | detect_communities vs igraph | igraph cluster_louvain/leiden |
| test-validation.R | 32 | Haldane-corrected pairwise OR | closed form |
| test-validation.R | 71 | chisq vs chisq.test | stats::chisq.test |
| test-validation.R | 110 | adjust_pvalues | hand-computable Bonferroni/Holm |

## NO-VALUES (80)

test-bh-family.R: 53, 81, 104, 111, 118
test-communities.R: 7, 16, 23, 33, 42, 92, 100, 109, 140, 155
test-hdp-known-answers.R: 54 (self-check of the pinned Julia trace against two stopping rules), 74, 86
test-hdp-multiseed.R: 144, 165, 183, 203, 217, 232, 238, 252
test-hdp-visualization.R: 9, 13, 17, 21, 26, 33, 40, 56
test-hdp.R: 90, 109, 123, 129, 140, 175, 198, 207, 245, 287, 320, 362
test-mining.R: 86, 97, 106, 115, 125
test-network-stratification.R: 34, 65, 77, 120, 127, 135, 145, 152, 217, 274
test-network.R: 14, 40, 66, 75, 98, 107, 146, 152, 158
test-package-metadata.R: 3
test-validation.R: 42, 81, 121, 131, 140, 147, 154, 159, 165

## Flags

1. **HDP known answers: the ELBO is not monotone.** Exact CAVI never lowers the ELBO. The pinned
   traces fall: flat_a1.3_g0.7 drops from -284.79097 (iter 12) to -284.79137 (iter 16), and
   eb_q3_p29 drops from -242.4831 (iter 9) to -242.5299 (iter 13) before rising again. This is
   expected only if some update is not an exact coordinate maximizer (for example a point-estimate
   or truncated update of the top-level stick beta). Otherwise the ELBO formula or an update has a
   bug that Julia and R share. The test file's own comment treats the dip as a feature ("a
   best-ELBO stopping rule sees no improvement"). This most needs an independent ELBO evaluator.
2. **The ARI convention for n = 1 differs.** `adjusted_rand_index(1, 1)` is pinned to 0
   (test-hdp-multiseed.R:66), but the helper `adjusted_rand()` returns 1 when nc2 = 0. Both are
   conventions, and the suite holds both.
3. **The round_digits tie rule is Julia's convention, not correct rounding.** 0.15 is stored as
   0.1499..., so the correctly rounded answer is 0.1, while Julia and R both give 0.2. This is
   deliberate, but its only source is Julia.
4. **chisq values on make_test_event_df rest on tables with expected cells below 5.** The numbers
   are right (verified), but the test as a statistic is questionable there.
