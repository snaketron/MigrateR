# cellvel
Quantification of cell migration speed with Bayesian probability models

## Overview 
Cell imaging allows us to profile the migration speeds of cells under various
treatment groups (e.g. chemical compounds at different doses) and biological 
conditions (e.g. cancer vs. healthy cells). These experiments are very 
expensive, and this implies that few biological replicates can be screened 
over the course of an experiment. Meanwhile, the cell migration readouts 
are often noisy and prone to batch effects. As a result of these factors, we 
get highly uncertain estimates of the migration speeds in our cell samples,
and these challenges are compounded by use of inadequate statistical methods. 

cellvel implements hierarchical Bayesian models to describe noisy and sparse 
cell migration data. It accounts for potential batch effects, and provides 
robust quantitative estimates of the cell speed in each sample, while also
performing differential cell migration analysis between different treatment.
cellel implements functions that facilitate planning of future experiments, 
providing answer to questions such as: "how many replicates and cells should
we screen in order to capture a given effect ?".

## How to use cellvel

scBubbletree is an R-package

https://github.com/snaketron/cellvel

To install this package, start R and enter:

```r
library("devtools")
devtools::install_github("snaketron/cellvel")
```

Case studies are provided in the directory /vignettes

## Workflow & output 

![alt text](inst/extdata/logo.png)
