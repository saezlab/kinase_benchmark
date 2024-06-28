# Systematic comparison of methods for kinase activity estimation

<!-- badges: start -->
<!-- badges: end -->

## Overview
In this study, we present a flexible framework to assess different combinations of computational algorithms and kinase-substrate libraries 
for the inference of kinase activities. For the benchmark, we use a set of kinase perturbation experiments to evaluate which combination
is able to recapitulate the perturbed kinases from the phosphoproteomics data.

If you want to test your own method try out our package [benchmarKIN](https://github.com/saezlab/benchmarKIN) and check out the 
[documentation](https://benchmarkin.readthedocs.io/). 

## Kinase substrate libraries
We have included the following kinase-substrate libraries:
- PhosphoSitePlus
- PTMsigDB
- Omnipath
- Gold Standard set of GPS 6.0
- iKiP-db
- NetworKIN

Additionally have tested the combination with predicted targets including the Kinase Library and Phosformer.

## Methods
We have included the following methods for the comparison:
- fgsea
- KARP
- Kologomorov-Smirnov
- linear model (RoKAI, decoupler)
- mean
- normalised mean
- median
- multivarite linear model
- PC1
- PTM-SEA
- sum
- upper quantile
- viper
- Wilcox
- z-score (KSEA, RoKAI)
