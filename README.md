# When Does the Normal Approximation Break Down?

An interactive Shiny app for BIOL 3P96 (Biostatistics) at Brock University.

## What this app does

Many familiar statistical methods for counts quietly rely on a smooth,
bell-shaped normal curve, even though counts are whole numbers that can never
fall below zero. Students are usually given rules of thumb for when this
shortcut is safe:

- **Binomial counts:** use the normal approximation only if np ≥ 5 and
  n(1 − p) ≥ 5 (at least 5 expected successes and 5 expected failures).
- **Chi-square tests:** no expected count below 1, and no more than 20% of
  expected counts below 5.

This app shows where those rules come from, what goes wrong when they are
broken, and why the two rules are really the same idea. Each category in a
chi-square test is itself a binomial count, so the chi-square rules are the
binomial problem happening inside every cell of the table.

Two tabs build the idea step by step:

1. **Binomial counts** — Bars show the true binomial distribution with the
   normal curve laid on top. When the expected count is squeezed against zero
   (or against n), part of the normal curve spills onto impossible counts,
   shaded in orange, and the real distribution becomes lopsided. A second plot
   shows 95% confidence intervals from 50 simulated studies, using either the
   Wald interval (estimate ± 1.96 SE) or the Agresti–Coull interval (add about
   2 successes and 2 failures first). The exact long-run coverage of each
   method is reported, so students can see the Wald interval fall well short
   of 95% when np is small, and see a simple fix rescue it.
2. **Chi-square counts** — You set the total sample size, the number of
   categories, and how many of them are rare. The expected counts are shown
   colour-coded against the 5 and 1 thresholds. The test is then run on
   thousands of simulated studies in which the null hypothesis is true, and the
   histogram of chi-square statistics is compared with the theoretical
   chi-square curve. The false positive rate is reported against the nominal
   5%. A second view shows the rarest category's count on its own, using the
   same plot as Tab 1.

## Key ideas illustrated

- A normal curve keeps about 95% of its area within 2 SDs of its centre. When
  np = 5, the centre sits at least √5 ≈ 2.2 SDs above zero, which is where the
  "5" comes from
- These rules are rules of thumb, not laws: the approximation improves
  gradually, and some textbooks use 10 instead of 5
- The Wald confidence interval for a proportion can badly undercover when np
  is small, and collapses to zero width when a study sees no successes
- The chi-square statistic only follows the chi-square curve when each
  category's count is roughly bell-shaped, which needs large enough expected
  counts
- A category with an expected count below 1 acts like a switch: one or two
  chance individuals can decide the outcome of the whole test
- The rules apply to *expected* counts, not observed counts

## How to use

1. On Tab 1, move the n and p sliders. Watch the bars pile up against an edge
   and the orange area grow as np or n(1 − p) drops below 5. Press
   **Resample** to draw 50 new studies, and switch between **Wald** and
   **Agresti–Coull** to compare the same studies under each method.
2. On Tab 2, change the sample size, number of categories, and number of rare
   categories. Press **Resample** to rerun the simulated studies. Switch the
   right-hand plot to **The rarest category on its own** to see the link back
   to Tab 1.

## Learning goals

- Explain why a count that is close to zero cannot be well described by a
  normal curve
- Connect the np ≥ 5 rule to the "within 2 SDs" property of the normal curve
- Describe what goes wrong with a confidence interval or test when the
  approximation fails
- Explain why the chi-square expected-count rules follow from the binomial
  rule
- Identify alternatives when the rules are not met: exact methods such as
  `binom.test()`, adjusted intervals such as Agresti–Coull, combining rare
  categories, or simulated p-values with
  `chisq.test(..., simulate.p.value = TRUE)`

## Notes

- Chi-square statistics are calculated without a continuity correction.
- Agresti–Coull intervals are shown untruncated, so lower limits can dip
  slightly below 0 at small np. Most software truncates these at 0.
- With 5000 simulated studies in Tab 2, differences of up to about half a
  percent from 5% in the false positive rate are just simulation noise.

## References

- Agresti, A. & Coull, B.A. (1998). Approximate is better than "exact" for
  interval estimation of binomial proportions. *The American Statistician*
  52:119–126.
- Brown, L.D., Cai, T.T. & DasGupta, A. (2001). Interval estimation for a
  binomial proportion. *Statistical Science* 16:101–133.
- Cochran, W.G. (1954). Some methods for strengthening the common χ² tests.
  *Biometrics* 10:417–451.

## Course context

Developed for BIOL 3P96 — Biostatistics, Brock University.
Built with R and Shiny (base R graphics only).
