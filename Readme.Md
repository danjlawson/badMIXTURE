# badMIXTURE: An R package for how (not) to over-interpret STRUCTURE/ADMIXTURE bar plots

## Summary

In [A tutorial on how (not) to over-interpret STRUCTURE/ADMIXTURE bar plots](http://biorxiv.org/content/early/2016/07/28/066431), Daniel Falush, Lucy van Dorp and Daniel Lawson discuss how mixture solutions to understand population genetics data may be misinterpreted. This is the R package that performs the analysis in that paper.

## Abstract

Genetic clustering algorithms, implemented in popular programs such as STRUCTURE and ADMIXTURE, have been used extensively in the characterisation of individuals and populations based on genetic data. A successful example is reconstruction of the genetic history of African Americans who are a product of recent admixture between highly differentiated populations. Histories can also be reconstructed using the same procedure for groups which do not have admixture in their recent history, where recent genetic drift is strong or that deviate in other ways from the underlying inference model. Unfortunately, such histories can be misleading. We have implemented an approach to assessing the goodness of fit of the model using the ancestry “palettes” estimated by CHROMOPAINTER and apply it to both simulated and real examples. Combining these complementary analyses with additional methods that are designed to test specific hypothesis allows a richer and more robust analysis of recent demographic history based on genetic data.

## What you need

To replicate these analyses you need to:

* Run [ADMIXTURE](https://www.genetics.ucla.edu/software/admixture/)/[STRUCTURE](http://pritchardlab.stanford.edu/structure.html) or similar to obtain a mixture solution
* Run [mixPainter](www.paintmychromosomes.com) (a fork of finestructure version 2) to use chromosome painting to compare all samples to all populations correctly. Note that you can do this yourself with chromopainter, if you know what you are doing, though there are many pitfalls.
* Be familiar enough with R to load these data matrices
* Follow the [example from plink data](https://github.com/danjlawson/badMIXTUREexample)
* Follow the examples in this package from chromopainter data (which are based on the above)

## How to install this package

```R
install.packages("devtools")
library(devtools)
install_github("danjlawson/badMIXTURE")
```

## Important places to look for inline help

```R
?badMIXTURE # For general instructions
?compareMixtureToData # For the main function that processes the data
?compareMixtureToDataDirect # For more manual reordering etc
?mixturePlot # For plotting
```


# EXAMPLE:

This example is not complete...

```
library("badMIXTURE")
## READ THE "Recent" Scenario 
recent_ariids=read.table("data/Recent_admix.ids",as.is=T)

recent_ariQ=readQ("data/Recent_admix.pruned.11.Q",recent_ariids,poporder)
recent_ariU=readChunkcounts("data/Recent_admix.pruned_ul_cp_pop_unlinked.chunkcounts.out")
recent_ariL=readChunkcounts("data/Recent_admix.pruned_cp_pop_linked.chunkcounts.out")

## READ THE "Remnants" Scenario 
remnants_ariids=read.table("data/Remnants_admix.ids",as.is=T)

remnants_ariQ=readQ("data/Remnants_admix.pruned.11.Q",remnants_ariids,poporder)
remnants_ariU=readChunkcounts("data/Remnants_admix.pruned_ul_cp_pop_unlinked.chunkcounts.out")
remnants_ariL=readChunkcounts("data/Remnants_admix.pruned_cp_pop_linked.chunkcounts.out")

## READ THE "Marginalisation" Scenario 
marginalisation_ariids=read.table("data/Marginalisation_admix.ids",as.is=T)

marginalisation_ariQ=readQ("data/Marginalisation_admix.pruned.11.Q",remnants_ariids,poporder)
marginalisation_ariU=readChunkcounts("data/Marginalisation_admix.pruned_ul_cp_pop_unlinked.chunkcounts.out")
marginalisation_ariL=readChunkcounts("data/Marginalisation_admix.pruned_cp_pop_linked.chunkcounts.out")


###########################
## Colours 
## For the ancestral populations
mycols=rep("",11)
mycols[1:3]=c("#00008899","#00880099","#BB00BB99")
mycols[4:11]="lightgrey"

## For the observable populations
mycols2=rep("",12)
mycols2[1:4]=c("#000088FF","#333333FF","#008800FF","#BB00BBFF")
mycols2[5:12]="lightgrey"

###########################
## Apply badMIXTURE

recent_admU<-compareMixtureToData(recent_ariQ,recent_ariU,ids=recent_ariids,
                           mycols=mycols,mycols2=mycols2,
                           poporder=poporder,remself=remself)
recent_admL<-compareMixtureToData(recent_ariQ,recent_ariL,ids=recent_ariids,
                           mycols=mycols,mycols2=mycols2,
                           poporder=poporder,remself=remself)

remnants_admU<-compareMixtureToData(remnants_ariQ,remnants_ariU,ids=remnants_ariids,
                           mycols=mycols,mycols2=mycols2,
                           poporder=poporder,remself=remself)
remnants_admL<-compareMixtureToData(remnants_ariQ,remnants_ariL,ids=remnants_ariids,
                           mycols=mycols,mycols2=mycols2,
                           poporder=poporder,remself=remself)

marginalisation_admU<-compareMixtureToData(marginalisation_ariQ,marginalisation_ariU,ids=marginalisation_ariids,
                           mycols=mycols,mycols2=mycols2,
                           poporder=poporder,remself=remself)
marginalisation_admL<-compareMixtureToData(marginalisation_ariQ,marginalisation_ariL,ids=marginalisation_ariids,
                           mycols=mycols,mycols2=mycols2,
                           poporder=poporder,remself=remself)

## Get the position of the rightmost individual that we are interested in
txmaxes=c(cumsum(sapply(recent_admU$poplist,length))[4],
          cumsum(sapply(remnants_admU$poplist,length))[4],
          cumsum(sapply(marginalisation_admU$poplist,length))[4])+9 # 9 = 3 *3 as there are 3 gaps, each of 3 individuals
          
rmunlinked=0.0005
pdf("AriSims_Test_Unlinked.pdf",height=10,width=6)
mixturePlot(recent_admU,main="Recent Unlinked",residual.max=rmunlinked,xlim=c(0,txmaxes[1]))
mixturePlot(remnants_admU,main="Remnants Unlinked",residual.max=rmunlinked,xlim=c(0,txmaxes[2]))
mixturePlot(marginalisation_admU,main="Marginalisation Unlinked",residual.max=rmunlinked,xlim=c(0,txmaxes[3]))
dev.off()

rmlinked=0.002
pdf("AriSims_Test.pdf",height=10,width=6)
mixturePlot(recent_admL,main="Recent Linked",residual.max=rmlinked,xlim=c(0,txmaxes[1]))
mixturePlot(remnants_admL,main="Remnants Linked",residual.max=rmlinked,xlim=c(0,txmaxes[2]))
mixturePlot(marginalisation_admL,main="Marginalisation Linked",residual.max=rmlinked,xlim=c(0,txmaxes[3]))
dev.off()
```
