# Robust likelihood ratio test under measurement errors, version 1.0

The implementation in R software performs a robust hypothesis test
comparing two samples, which takes into account measurement errors.
This is a very unique approach to testing under the presence
of measurement errors. The presented code performs a complete analysis,
including the computation of the test statistic for a given dataset, the 
critical value, and evaluates the power of the test. Computing
the critical value as well as the power is very tedious but feasible
with this software for moderate sample sizes.

Feel free to use or modify the code.

## Requirements

No additional packages for R software are needed.

## Usage

* The analysis must be performed in this order, starting with computing the test statistic, then the critical value, and finally power of the test:
TestStatistic.R, CriticalValue.R, Power.R.

## Authors
  * Jan Kalina, The Czech Academy of Sciences, Institute of Computer Science
  * Michel Broniatowski,  Université Pierre et Marie Curie (Sorbonne Université)
  * Jana Jurečková, Faculty of Mathematics and Physics, Charles University & The Czech Academy of Sciences, Institute of Information Theory and Automation

## Contact

Do not hesitate to contact us (kalina@cs.cas.cz) or write an Issue.

## How to cite

Please consider citing the following:

Broniatowski M, Jurečková J, Kalina J (2018): Likelihood ratio testing under measurement errors. Entropy 2018, 20, Article 966.

## Acknowledgement

This work was supported by the Czech Science Foundation grants GA19-05704S and GA18-01137S.