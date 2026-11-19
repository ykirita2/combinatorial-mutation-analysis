# Combinatorial epistasis analysis for evoCas12j

Code and input data for the combinatorial-library analyses reported in
*Complete combinatorial mutagenesis resolves antagonistic epistasis and yields
a hypercompact Cas12j genome editor* (Toyoshima et al.).

Two complete combinatorial libraries were measured: every combination of the
11 substitutions selected in the N-terminal-segment scan (2,048 genotypes, of
which 1,823 passed the read-count threshold) and of the 10 substitutions
selected in the C-terminal-segment scan (1,024 genotypes, 1,021 retained).
The scripts here take those measurements and produce every statistic and
figure panel derived from them.

## Contents

```
data/
  combinatorial_library_N.csv    N-segment library: genotype, replicate scores, read count
  combinatorial_library_C.csv    C-segment library: the same
  substitution_assignment.csv    library column -> substitution, scanning score
src/
  epistasis_analysis.py          all analyses; writes results/
  build_fig3.py                  main-text Figure 3
  build_figS3.py                 Supplementary Figure 2
results/                         output of epistasis_analysis.py
```

## Running

```
python src/epistasis_analysis.py --data data --out results
```

Python 3.13 with numpy, pandas and scipy. Runtime is a few minutes, dominated
by the strategy simulation. Every random draw is seeded (`SEED_STRATEGY`,
`SEED_CV` in the script), so repeated runs give identical output.

## What each output contains

| file | content |
| --- | --- |
| `library_summary.csv` | read-count filtering, replicate structure, the per-score standard error used for error propagation, and two independent estimators of it |
| `pairwise_interaction_detectability.csv` | for each of the 55 and 45 pairs: the four-genotype contrast on the parental background, the same contrast averaged over all measurable backgrounds, the number of backgrounds, whether either is resolved against measurement noise, and a chi-square test of whether the spread across backgrounds exceeds noise (Benjamini-Hochberg corrected) |
| `interaction_vs_background_load.csv` | mean contrast grouped by the number of further substitutions carried by the background |
| `local_optima.csv` | local optima under single-substitution moves and the basin of attraction of each |
| `substitution_background_effects.csv` | effect of each substitution on the parental genotype and across all backgrounds in which it could be measured |
| `model_order_cv.csv` | 10-fold cross-validated coefficient of determination for additive, pairwise and third-order models |
| `coefficient_transfer.csv` | correlation between pairwise coefficients fitted to genotypes carrying at most two substitutions and to the complete library |
| `strategy_benchmark.csv` | in silico comparison of seven engineering strategies, each with the same final validation step |

## Definitions

A pairwise interaction is the four-genotype contrast (a double-mutant cycle)

```
epsilon = y_ij - y_i - y_j + y_0
```

evaluated on every background in which all four genotypes were measured.
Positive values indicate synergy. Because the contrast is a sum of four
measurements, its standard error is twice that of a single score; the standard
error of the background-averaged value is that quantity divided by the square
root of the number of backgrounds. A pair is called resolved when the absolute
value of the estimate exceeds twice its standard error.

The strategy comparison treats each measured landscape as ground truth, adds
independent noise of the estimated measurement error once, runs each strategy,
and then applies one common selection step: the ten best candidates measured
by that strategy are re-measured independently and the genotype with the
highest average of the two measurements is selected. This step is required for
comparability, because taking the maximum observed score over a large measured
set favours genotypes whose noise happened to be positive.

## Figures

`build_fig3.py` and `build_figS3.py` reproduce the published panels. They read
the tables written by `epistasis_analysis.py`; see the header of each script
for the variables it expects.

## Data

The underlying sequencing data are deposited in the Sequence Read Archive under
BioProject PRJNA1514873. The cryo-EM map and model are deposited as EMD-83117
and PDB 44ZO.
