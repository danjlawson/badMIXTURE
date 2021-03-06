% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compareMixtureToData.R
\name{compareMixtureToDataDirect}
\alias{compareMixtureToDataDirect}
\title{Compare a mixture solution to some data, without reordering or changing any of the data.}
\usage{
compareMixtureToDataDirect(mix, data.NbyP, poplist, remself = 0,
  ancestral = "solve", adm = NULL, data.PbyP = NULL, tdend = NULL,
  mycols = NULL, mycols2 = NULL, gap = 3)
}
\arguments{
\item{mix}{A proposed mixture solution of dimension N by K; for example, as generated by STRUCTURE or ADMIXTURE.}

\item{data.NbyP}{A matrix containing the data of dimension N by P, for example, as generated by ChromoPainter. the data must reflect similarity to P clusterings of the data.}

\item{poplist}{A list of length P, with elements containing the rownames of the above two matrices which were used in constructing the data column p in data.NbyP (a transformation of the ids give to \code{\link{compareMixtureToData}}).}

\item{remself}{Number of iterations that "self-copying" (cluster specific sharing of drift) is removed. Set to 10 to essentially remove all self-copying.}

\item{ancestral}{Ancestral model, i.e. how the latent clusters are defined. There are two options. "mixture": meaning that the populations are defined as a mixture of the individuals that comprise them according to the admixture model. This is not advisable because it allows ancestry from different true latent popluations to affect inference of others. "solve": find the definition of the clusters that best explains the data by root-mean-sqaure-distance. Default: "solve"}

\item{adm}{The result of a call to compareMixtureToData. This is used to collect the plotting properties below (poplist, tdend, mycols, mycols2, gap)}

\item{data.PbyP}{=NULL, # Needed for colours only}

\item{tdend}{Dendrogram of the populations to determine plot order. Default: NULL, meaning create this from the data}

\item{mycols}{The colour for each of the K ancestral populations; Default: NULL, meaning use rainbow(K). Can be modified in the returned object instead.}

\item{mycols2}{The colour for each of the P clusters; Default: NULL, meaning determined using \code{\link{rgbDistCols}} so that similar clusters have similar colours. Can be modified in the returned object instead.}

\item{gap}{The spacing between populations, relative to the spacing between individuals which is 1. Default: 3}
}
\value{
An object of class admixfs ; see \code{\link{compareMixtureToData}}.
}
\description{
This function takes a mixture solution, a data matrix, and a mapping of data observations into clusters. It then predicts what it expects the mixture solutions to look like.

This function assumes that the data are ordered according to the clusters, and the clusters are ordered according to the palette. Arranging this is non-trivial so you may wish to use the version of this code that does this for you, \code{\link{compareMixtureToData}}). You can then easily call this function by reordering rows/columns of the data and the mixture, as shown in the example.
}
\examples{
\dontrun{
data(arisimsmall)
## Example where we reorder the populations manually
adm0<-compareMixtureToData(arisimsmall$mixture,arisimsmall$data,arisimsmall$fam)
mypoplist=adm0$poplist[paste0("Pop",c(13,5:7,1:4,9,11,12))]
mymix=adm0$mix[unlist(mypoplist),]
mydata=adm0$data.NbyP[unlist(mypoplist),names(mypoplist)]
adm<-compareMixtureToDataDirect(mymix,mydata,mypoplist)
gpref=mixturePlot(adm,cex.names=0.5,residual.scale=0.2,height.cluster=.4)
}
}
\keyword{admixture}
\keyword{chromopainter}
\keyword{finestructure}
\keyword{mixture}
\keyword{structure}
