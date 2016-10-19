# mixChecker: An R package for how (not) to over-interpret STRUCTURE/ADMIXTURE bar plots

## Summary

In [A tutorial on how (not) to over-interpret STRUCTURE/ADMIXTURE bar plots](http://biorxiv.org/content/early/2016/07/28/066431), Daniel Falush, Lucy Lan Dorp and Daniel Lawson discuss how mixture solutions to understand population genetics data may be misinterpreted. This is the R package that performs the analysis in that paper.

## Abstract

Genetic clustering algorithms, implemented in popular programs such as STRUCTURE and ADMIXTURE, have been used extensively in the characterisation of individuals and populations based on genetic data. A successful example is reconstruction of the genetic history of African Americans who are a product of recent admixture between highly differentiated populations. Histories can also be reconstructed using the same procedure for groups which do not have admixture in their recent history, where recent genetic drift is strong or that deviate in other ways from the underlying inference model. Unfortunately, such histories can be misleading. We have implemented an approach to assessing the goodness of fit of the model using the ancestry “palettes” estimated by CHROMOPAINTER and apply it to both simulated and real examples. Combining these complementary analyses with additional methods that are designed to test specific hypothesis allows a richer and more robust analysis of recent demographic history based on genetic data.

## What you need

To replicate these analyses you need to:

* Run [ADMIXTURE](https://www.genetics.ucla.edu/software/admixture/)/[STRUCTURE](http://pritchardlab.stanford.edu/structure.html) or similar to obtain a mixture solution
* Run [FineSTRUCTURE](www.paintmychromosomes.com) or similar to obtain a Similarity matrix
* Be familiar enough withR to load these data matrices
* Follow the examples in this package

## How to install this package

```R
install.packages("devtools")
library(devtools)
install_github("danjlawson/mixChecker")
```

## Important places to look for inline help

```R
?mixCecker # For general instructions
?compareMixtureToData # For the main function that processes the data
?compareMixtureToDataDirect # For more manual reordering etc
?mixturePlot # For plotting
```R

