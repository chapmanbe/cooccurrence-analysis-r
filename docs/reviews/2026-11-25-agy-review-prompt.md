You are an adversarial code reviewer. Do NOT edit any file. Read only.

Repository: the current directory, an R package `CooccurrenceAnalysis` that ports a Julia
package of the same name. Read `CLAUDE.md` first: its "Parity rules" and "Divergences from
the Julia reference" sections are the specification. The Julia reference, if you need it,
is at `../cooccurrence-analysis`.

Scope: the changes in `git diff f849c40..HEAD -- R/ tests/ NAMESPACE DESCRIPTION`
(run that command yourself). This code has had no review. It covers group-stratified
networks and `exclusive_items`, the HDP Bernoulli mixture and multi-seed fitting
(best-ELBO selection, ARI, LSAP cluster matching), chi-square and Bonferroni/Holm,
rule mining via `arules`, Leiden/Louvain via igraph, and seed/RNG handling.

Hunt for, in priority order:
1. Numerical or statistical errors: wrong formulas, off-by-one, wrong tail, wrong
   multiple-testing family, rounding not via `round_digits()`, ELBO/variational updates
   that are mathematically wrong, label-switching mistakes in cluster matching.
2. Silent failures: NA dropped instead of raised, `sapply`, `suppressWarnings`, `try`,
   partial matching, inputs that should raise but return something plausible.
3. Reproducibility: seeds that leak into or clobber the caller's RNG stream, results that
   depend on locale or on igraph version.
4. Tests that cannot fail: expected values derived from the R code itself rather than from
   the Julia reference or a fixture unambiguous by construction; assertions too loose to
   catch the bug they claim to guard.
5. Violations of the CLAUDE.md conventions (hardcoded domain terms, missing NAMESPACE
   exports, plots that save).

Do not report style nits. For every finding, give a concrete failing input or state.
If you are unsure, say so rather than overstate.

Output: a JSON array only, no prose before or after. Each element:
{"file": "R/....R", "line": 123, "severity": "high|medium|low",
 "category": "numerical|silent-failure|reproducibility|test-weakness|convention",
 "summary": "one sentence", "failure_scenario": "concrete input -> wrong result",
 "confidence": "high|medium|low"}
Return [] if you find nothing.
