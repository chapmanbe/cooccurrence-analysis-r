# Brief for the Julia `CooccurrenceAnalysis` repository: what the R port found, and what its history should answer

Written 2026-11-25 from the R package's repository
(`/Users/brainchapman/Code/Julia/cooccurrence-analysis-r`, pushed to
`github.com/chapmanbe/cooccurrence-analysis-r` at `e083f5c`). Read this once, then use the
Julia repository's own code, docs and git history (`git log`, `git log -S`, `git blame`) to
answer the questions in section 4.

**Ground rules for the session that reads this:**
- Answer from evidence: a file and line, a commit hash, a doc, or a cited paper. Where the
  history doesn't say, write "not recorded". Don't reconstruct intent.
- Don't change the Julia code. This is research; fixes are a separate decision for Brian.
- Write the answers to `julia-history-answers.md`, next to this file, one section per
  question, each answer under 150 words with its evidence.

## 1. Background

The R package began on 2026-09-25 as a port of this Julia package, with tests whose expected
values came from running Julia. The port is now finished, and the R package no longer treats
Julia as its reference: Julia had never been reviewed the way the R code has, and porting
turned up four Julia bugs (section 3). R's working rules are in its `CLAUDE.md`.

36 of the R repository's 39 commits were written with Claude (Claude Code); Brian Chapman
directed the work and made the design and modeling decisions.

## 2. How the R package was checked

- **Independent oracles.** Every expected value in a test must come from a fixture with one
  right answer, a hand calculation from the definition, or an established implementation
  (`stats::fisher.test`, `chisq.test`, `p.adjust`, `cooccur`, `flexmix`, igraph, `arules`).
  Julia can serve only as a cross-check alongside one of those.
- **Oracle audit, 2026-11-25** (`docs/reviews/2026-11-25-oracle-audit.md` in the R repo):
  10 of 150 R tests took their values from Julia alone. Every Julia value that had a practical
  oracle was recomputed independently and matched to 1e-12, so no shared bug turned up in
  them. The remaining Julia-only values are the four HDP known-answer CAVI traces.
- **Independent ELBO evaluator** (`tests/testthat/helper-hdp-elbo.R`): written from the
  model's definition, it reproduces R's and Julia's ELBO at every iteration of those four
  traces (1e-10 against R, 1e-9 against Julia). So **the ELBO formula is right**. The ELBO
  still falls between some iterations, which means the **updates** are not exact coordinate
  ascent on that bound (section 4, question 1).
- **Two-model code review, 2026-11-25** (`docs/reviews/2026-11-25-code-review.md`): Claude
  (`/code-review`) and Gemini (through Google's Antigravity CLI) reviewed the R code
  independently. Gemini found all six of Claude's bugs and added 13 more findings; Codex was
  planned but did not run. All 23 findings were fixed, each bug fix starting with a failing
  test. Several of the R bugs (seed handling, `arules` limits, id coercion) are R-specific and
  say nothing about Julia.

## 3. Julia bugs the port found

Each is fixed in R and covered by a test marked `R-only`.

1. **Chi-square on a zero margin returns `NaN`**, when an item is in every record. The `NaN`
   then makes every BH q-value in the family `NaN`. R raises instead.
2. **Louvain is phase 1 only** (local moving, no aggregation; the R notes cite
   `network_metrics.jl:180`). It stalls below the optimum: Q = 0.70 on a ring of 20 triangles,
   where full Louvain reaches 0.77. R uses igraph's Leiden (default) or Louvain, so partitions
   differ from Julia's.
3. **`min_count` in mining is round-tripped**: converted to a fraction, then back by RuleMiner
   with `ceil(min_support * n)`, which can land one above it (7 of 100 gives 8).
4. **The adjusted Rand index of fewer than two elements is 0.** The definition gives 1, since
   no pair can disagree.

R keeps several Julia choices deliberately, because published results depend on them. Fisher
tests are one-sided; the BH family is exactly the pairs that co-occur in at least one record;
odds ratios get the Haldane correction (+0.5 to all four cells when any cell is zero); and
`round_digits()` rounds ties the Julia way (`round_digits(0.15, 1)` is 0.2, although 0.15 is
stored as 0.1499…).

## 4. Questions for the Julia history

1. **Where did the HDP CAVI update scheme come from?** Two steps are not exact coordinate
   ascent:
   - the global-stick update pools raw soft counts across groups rather than table counts
     (R: `R/hdp_clustering.R`, near line 119);
   - the group-stick prior Beta(α·b_k, α·Σ_{l>k} b_l) uses the global weights b at their
     variational mean, a plug-in for an expectation with no closed form.

   Find the Julia code for both and any comment, doc, commit message or reference explaining
   them. Was the scheme taken from a paper (for example Liang et al. 2007 or Wang, Paisley &
   Blei 2011, below), derived here, or written ad hoc? Was the falling ELBO ever noticed?
2. **Was Louvain phase 1 a deliberate simplification or an oversight?** Was it ever compared
   with a full Louvain or Leiden implementation?
3. **Were any of the four bugs in section 3 known?** Look for issues, TODOs, comments or
   tests that mention them.
4. **Why one-sided Fisher, and why the BH family of co-occurring pairs?** Where were those
   decisions made, and with what reasoning?
5. **The tie rule in `round_digits`**: was matching some earlier table the reason, or was it
   Julia's default `round` behavior that later became fixed because results were published?
6. **Which published results came from the Julia package**, and from which commit? The fixes
   in section 3, especially the Louvain one, could change them.
7. **What references does the Julia code or its docs cite?** Compare them against section 5
   and report any that are missing, wrong, or cited for a method the code doesn't implement.

## 5. Published sources for the methods (checked 2026-11-25)

Neither package makes a statistical innovation: both implement methods published by others.
Each citation below was checked against Crossref by DOI, or against the publisher's page where
there is no DOI.

| Method | Citation | Checked against |
|---|---|---|
| Pearson chi-square | Pearson, K. (1900). On the criterion that a given system of deviations…. *Philosophical Magazine* 50(302):157–175. doi:10.1080/14786440009463897 | Crossref |
| Fisher exact test | Fisher, R. A. (1935). The logic of inductive inference. *JRSS* 98(1):39–54 (the discussion follows at 55–82). doi:10.2307/2342435 | Crossref; pages from the University of Adelaide repository |
| Benjamini-Hochberg | Benjamini, Y. & Hochberg, Y. (1995). Controlling the false discovery rate…. *JRSS B* 57(1):289–300. doi:10.1111/j.2517-6161.1995.tb02031.x | Crossref |
| Holm | Holm, S. (1979). A simple sequentially rejective multiple test procedure. *Scand. J. Statist.* 6(2):65–70. JSTOR 4615733 | secondary sources only; the DOI is not in Crossref |
| Haldane-Anscombe correction | Haldane, J. B. S. (1956). The estimation and significance of the logarithm of a ratio of frequencies. *Ann. Hum. Genet.* 20(4):309–311. doi:10.1111/j.1469-1809.1955.tb01285.x · Anscombe, F. J. (1956). On estimating binomial response relations. *Biometrika* 43(3–4):461–464. doi:10.1093/biomet/43.3-4.461 | Crossref |
| Lift (as "interest") | Brin, S., Motwani, R., Ullman, J. D. & Tsur, S. (1997). Dynamic itemset counting and implication rules for market basket data. *SIGMOD '97*, 255–264. doi:10.1145/253260.253325 | Crossref |
| Disease co-occurrence networks | Hidalgo, C. A., Blumm, N., Barabási, A.-L. & Christakis, N. A. (2009). A dynamic network approach for the study of human phenotypes. *PLoS Comput. Biol.* 5(4):e1000353. doi:10.1371/journal.pcbi.1000353 | Crossref |
| Modularity | Newman, M. E. J. & Girvan, M. (2004). Finding and evaluating community structure in networks. *Phys. Rev. E* 69(2):026113. doi:10.1103/PhysRevE.69.026113 | Crossref |
| Louvain | Blondel, V. D., Guillaume, J.-L., Lambiotte, R. & Lefebvre, E. (2008). Fast unfolding of communities in large networks. *J. Stat. Mech.* P10008. doi:10.1088/1742-5468/2008/10/P10008 | Crossref |
| Leiden | Traag, V. A., Waltman, L. & van Eck, N. J. (2019). From Louvain to Leiden: guaranteeing well-connected communities. *Sci. Rep.* 9:5233. doi:10.1038/s41598-019-41695-z | Crossref (article number not returned) |
| Stick-breaking | Sethuraman, J. (1994). A constructive definition of Dirichlet priors. *Statistica Sinica* 4(2):639–650 | journal's own page |
| HDP | Teh, Y. W., Jordan, M. I., Beal, M. J. & Blei, D. M. (2006). Hierarchical Dirichlet processes. *JASA* 101(476):1566–1581. doi:10.1198/016214506000000302 | Crossref |
| Variational DP mixtures | Blei, D. M. & Jordan, M. I. (2006). Variational inference for Dirichlet process mixtures. *Bayesian Analysis* 1(1):121–143. doi:10.1214/06-BA104 | Crossref |
| Truncated HDP variational inference (candidates for question 1) | Liang, P., Petrov, S., Jordan, M. & Klein, D. (2007). The infinite PCFG using hierarchical Dirichlet processes. *EMNLP-CoNLL 2007*, 688–697 · Wang, C., Paisley, J. & Blei, D. M. (2011). Online variational inference for the hierarchical Dirichlet process. *AISTATS 2011*, PMLR 15:752–760 | ACL Anthology; PMLR |
| Apriori | Agrawal, R. & Srikant, R. (1994). Fast algorithms for mining association rules. *VLDB 1994*, 487–499 | DBLP via search; the VLDB PDF refused direct download |
| Adjusted Rand index | Hubert, L. & Arabie, P. (1985). Comparing partitions. *J. Classification* 2(1):193–218. doi:10.1007/BF01908075 | Crossref |

These are the R packages R delegates to or checks against: Griffith, Veech & Marsh (2016),
`cooccur`, *J. Stat. Softw.* 69(c02), doi:10.18637/jss.v069.c02, built on Veech (2013),
*Global Ecol. Biogeogr.* 22(2):252–260, doi:10.1111/j.1466-8238.2012.00789.x; Hahsler, Grün &
Hornik (2005), `arules`, *J. Stat. Softw.* 14(15), doi:10.18637/jss.v014.i15; Leisch (2004),
`flexmix`, *J. Stat. Softw.* 11(8), doi:10.18637/jss.v011.i08.

Whether Liang et al. or Wang et al. match the Julia HDP updates is **unverified**. They are
listed as the closest published forms to check against, not as the scheme's source.
