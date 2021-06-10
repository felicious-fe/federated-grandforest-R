Federated Grand Forest
============

A federated graph-guided Random Forest algorithm

## Introduction

Grand Forest is a graph-guided Random Forest algorithm, integrating secondary graph-structured data in order guide the feature selection towards interacting features. While it can be used for prediction, the main purpose of Grand Forest is descriptive, as it provides an efficient way of discovering highly informative subnetworks. Grand Forest is based on [ranger](https://github.com/imbs-hl/ranger). This version of the federated GrandForest implements the ability to merge individually trained forests with the same model parameters on different data sets.

## Installation

To install the latest development version of the R package run:

```R
devtools::install_github("felicious-fe/federated-grandforest-R")
```

## Reporting bugs

If you find any bugs, or if you experience any crashes, please report them at <https://github.com/felicious-fe/federated-grandforest-R>.
