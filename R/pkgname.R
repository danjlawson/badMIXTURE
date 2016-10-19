#' @title Compare admixture estimates from STRUCTURE-like programs to palettes
#'
#' @description
#' mixChecker is a packge to check whether a mixture is a good description of data. This is thought of and described in terms of genetic data, but it is appropriate for other types of mixture.
#'
#' The key concept is that we can compare the similarity of a number of objects (individuals) to a number of reference points, here thought of as the mean of clusters of the individuals. It is not important whether these represent anything - they can even be random clusterings - although the better chosen they are, the more powerful the test is.
#'
#' Given a mixture, we expect the similarities to fit the mixture solution. This is explicitly proven for some special cases in genetics, but is likely to hold for many similarities, for example, covariance between locations in PCA (Prinicple Components Analysis).
#' 
#' We therefore need "data" (the set of similarities to each of the clusters, which can be computed from similarities to all individuals), a "mixture solution" to test, and a record of which individuals make up which clusters.
#'
#' The functions you're likely to need from \pkg{mixChecker} are first \code{\link{compareMixtureToData}}, which performs the calculations, and
#' \code{\link{mixturePlot}} which displays the results.  The example section for these functions contain several cases to get you started. See also \code{?arisim} for information about the included data, and \code{arisimsmall} for a subset of this dataset that is easy to experiment on.
#'
#' If you want to reorder, look at \code{\link{compareMixtureToDataDirect}} and the example in it.
#' 
#' @references Falush, Van-Dorp & Lawson \url{http://biorxiv.org/content/early/2016/07/28/066431}
"_PACKAGE"
#> [1] "_PACKAGE"
