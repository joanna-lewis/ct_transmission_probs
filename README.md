# Per-partnership transmission probabilities for _Chlamydia trachomatis_ infection: Evidence synthesis of population-based survey data

_Joanna Lewis<sup>(1,2)</sup>, Peter J. White<sup>(1,3)</sup> and Malcolm J. Price<sup>(4,5,6)</sup>_

_<sup>(1)</sup>National Institute for Health Research Health Protection Research Unit in Modelling Methodology and Medical Research Council Centre for Global Infectious Disease Analysis, Imperial College London School of Public Health, London, UK_

_<sup>(2)</sup>Population, Policy and Practice, UCL Great Ormond Street Institute of Child Health, 30 Guilford Street, London, UK_

_<sup>(3)</sup>Modelling and Economics Unit, National Infection Service, Public Health England, London, UK_

_<sup>(4)</sup>NIHR Birmingham Biomedical Research Centre, University Hospitals Birmingham NHS Foundation Trust and University of Birmingham, UK_

_<sup>(5)</sup>University Hospitals Birmingham NHS Foundation Trust, UK_

_<sup>(6)</sup>Institute of Applied Health Research, University of Birmingham, UK_

This repository contains R and Stan code to estimate the per-partnership transmission probabilities of _Chlamydia trachomatis_.

There are three Stan models:
* `transmission_probability_reparam.stan` is the model used for the main analysis.
* `transmission_probability_reparam_different_means.stan` is a model allowing the mean number of partners per year to differ between men and women (sensitivity analysis).
* `transmission_probability_reparam_sens_analysis.stan` is a sensitivity analysis to investigate the assumption of random mixing between men and women in the population.

The R scripts apply each of the models to two datasets: the Natsal-2 study from England and the three NHANES surveys conducted in the US between 2009 and 2013. There are six scripts (three models times two datasets), identified by their filenames which contain both the name of the stan model and the name of the dataset.

Finally, the script labelled `run_transmission_probability_reparam_natsal2_nonocon.R` contains a sensitivity analysis in which the model was applied to data on the estimated number of new partners without a condom, rather than the reported total number of new partners (regardless of condom use).