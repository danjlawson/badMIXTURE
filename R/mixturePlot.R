
###############################
#' @title Plot a mixture fit vs the data.
#' @description
#' This function takes a mixture solution as returned by \code{\link{compareMixtureToData}}) or \code{\link{compareMixtureToDataDirect}}) and plots it in a variety of ways.
#'
#' Focus your attention on only the first two arguments: \code{adc} which is the processed mixture solution) and \code{show} which is how you choose the plots to display. The remaining arguments affect only the presentation of the plot.
#'
#' It is essential to understand that there are two key structures being explored simultaneously in these plots. The first is the P-dimensional clustering of the data, which defines the similarities; each of the N data points has a similarity to these P clusters. A representation in this space is called a palette. The second is the K-dimensional set of admixture weights. Each of the K "ancestral" or latent variables also has a P dimensional palette. Further, each of the N data points has a K dimensional "admixture" breakdown, seen as a mixture of the K ancestral palettes.
#'
#' The variety of plots available each explore different aspects of the fit of the mixture solution to the data.  We assume that the individuals are ordered by cluster and plot clusters as separated. These must match the palette in order.
#' 
#' @section Plot types
#' @description 
#' A variety of plot types for "show" are available:
#' \itemize{
#'   \item Admixture : the Admixture profile, shown via the classic barplot. 
#'   \item Cluster : The assignment of individuals into the P clusters.
#'   \item Raw : The palettes, as observed in the data.
#'   \item Pred : The predicted palettes.
#'   \item Ancestral : The palettes of the ancestral populations.... See ancestralshow.
#'   \item Residual : The matrix of residuals. See resid.cols.
#'   \item Same : The part of the palette that was predicted and is present.
#'   \item Omni : A visualisation of the raw palette and its fit to the data. The palette is split into three parts: a) the predicted palette that is present in the data, shown immediately above the x-axis up to the black line. b) the additional palette that is present in the data but missing in the prediction, shown above the black line up to 1. c) the predicted palette that is not present in the data; this is shown beneath the x axis. Therefore everything shown above the x-axis represents the data, and everything below the black line represents the prediction.
#' }
#' 
#' There are also some other plots which are more experimental:
#' \itemize{
#'   \item Over : A comparison of the predicted palettes to the observed palette; the predicted palette is shown on top, with "over-prediction" shown in black. The "missing" part of the prediction is show at the bottom.
#'   \item Over2 : As "Over", but without "blacking out" the section of the prediction that was overpredicted.
##'   \item Under : 
##'   \item Under2 :
#' }
#'
#' @section Advanced usage
#' @description 
#' 
#' If you want to produce high quality plots, you may wish to control the plotting very carefully. It is still helpful to use these functions as they ensure that the x-axis aligns; however you might want to display your own annotation.
#'
#' This is possible by calling with a single "show" argument, and setting  layout=FALSE. Then the plot will be drawn in the next available plotting area. You can disable all text overlay and add this yourself. See examples.
#' 
#' @param adc A mixture solution with data, as returned by \code{\link{compareMixtureToData}}) or \code{\link{compareMixtureToDataDirect}}).
#' @param show A character vector of plot types to create. See plots below.
#' @param mainname Title for the main plot
#' @param uselayout Whether we call layout, or allow the user to set up the plot environment. Default: TRUE
#' @param popy height for the cluster names; repeated (default: 0.5)
#' @param height.cluster Height of the "cluster" plot relative to the others
#' @param cex Overall scale; see par. Default: 1.0
#' @param cex.main As described in par. Default: 1.0
#' @param cex.axis As described in par. Default: 1.0
#' @param cex.poplab cex for the population labels. Default: 0.6
#' @param cex.residlab cex for the residuals. Default: 0.6
#' @param cex.residleg cex for the legend in residual plots. Default: 0.6
#' @param cex.names cex for the names ...? Default: 0.2
#' @param cex.ancestral cex for the ancestral ...? Default: 1.0
#' @param marx Margin in the x direction. Default: c(3,1)
#' @param mary.admixture Margin in the y direction of the admixture plot. Default: c(2,3)
#' @param mary.cluster Margin in the y direction of the cluster plot. Default: c(1,1)
#' @param mary.diff Margin in the y direction of the difference, diff & omni plots. Default: c(2,2)
#' @param mary.ancestral Margin in the y direction of the ancestral plot. Default: c(1,2)
#' @param mary.residual Margin in the y direction of the residual plot. Default: c(1,1)
#' @param mary.raw Margin in the y direction of the raw plot. Default: c(1,1)
#' @param omni.lwd Line width for the omni plot. Default: 2.0
#' @param omni.ylim Limits of the omni plot. default: NULL, meaning the range of the data
#' @param line.above Line where text is written above the plots with mtext; Default: 0
#' @param line.below Line where text is written below the plots with mtext; Default: 0
#' @param xlim xlim for the plots; all plots are restricted to this range. Note that you must calculate the spaces between populations; see below. (Default: NULL, meaning the range of the data)
#' @param showtitles whether to plot any titles. Default: TRUE
#' @param ancestralspace Spacing between ancestral vectors (they are each of width 1). Default: 0.5
#' @param ancestralshow A vector of which of the ancestral populations to show. Default: NULL, meaning all of them.
#' @param residual.max Maximum value for the resiidual plots. Default:NULL, meaning the range of the data.
#' @param residual.scale Scale of the residuals in the AdmixtureResidualMatrix, as compared to the ..? Defualt: 0.1.
#' @param residual.showscale Whether to show the scale in the residual matrix plot. Default: TRUE
#' @param ancestral.mean Whether to show the global mean in the Ancestral plot. Default: TRUE
#' @param ancestral.names A vector of names for the ancestral populations in the Ancestral plot. Default: NULL, meaning label them A1..K.
#' @param ancestral.shownames Whether to show the names in the ancestral plot. Default: TRUE.
#' @param cluster.srt String rotation in the Cluster plot, as passed to par(srt). Default: 0.
#' @param cluster.names Whether to show the names in the Cluster plot. Default: TRUE.
#' @param resid.cols Colours for the Residual Matrix. Default: \code{\link{MakeColorRYWGB}})(), which makes one direction Yellow->Red, and the other Green->Blue, with white around zero.
#' @param gpref A list of the above preferences, passed to the helper functions that do the plotting
#' @keywords mixture chromopainter admixture structure finestructure

#' @return gpref A list of the graphical preferences that have been specified, so you can call this function again and get consistent output.
#' @export
#' @examples
#' \dontrun{
#' ## Straightforward example:
#' data(arisim_remnants)
#' adm<-compareMixtureToData(arisim_remnants$mixture,
#'                           arisim_remnants$data,arisim_remnants$ids)
#' gpref=mixturePlot(adm)
#'
#' ## Example showing individual names clearly
#' data(arisimsmall)
#' adm<-compareMixtureToData(arisimsmall$mixture,arisimsmall$data,arisimsmall$ids)
#' 
#' gpref=mixturePlot(adm,
#'                   c("Admixture","Cluster","Raw","Residual"),
#'                   cex.names=0.5,residual.scale=0.2,height.cluster=.4)
#'
#' ## Complex example using layout and data subsetting
#' data(arisimsmall)
#' adm<-compareMixtureToData(arisimsmall$mixture,arisimsmall$data,arisimsmall$ids)
#' adm2=adm
#' adm2$mycols2[]="lightgrey"
#' adm2$mycols2[4:7]=c("darkred","orange","green","blue")
#' 
#' layout(matrix(c(1,2,3,4,5,5,5,5),nrow=4,ncol=2),heights=c(1,0.4,1,1))
#' gpref=mixturePlot(adm2,
#'                   c("Admixture","Cluster","Raw","Residual"),
#'                   mainname=paste0("Admixture, Marginalization, Direct"),
#'                   mary.admixture=c(0,1),
#'                   cex.names=0,
#'             ancestralshow=c(1,2,4), 
#'             residual.max=NULL,
#'             residual.scale=0.2,
#'             residual.showscale=TRUE, 
#'             ancestral.mean=FALSE,
#'             xlim=c(24,56),
#'             uselayout=FALSE)
#' gpref=mixturePlot(adm2,"Ancestral",gpref=gpref,uselayout=FALSE)
#' }

mixturePlot<-function(adc,
                      show=c("Admixture","Cluster","Residual"),
                      mainname="Admixture",
                      uselayout=TRUE, 
                      popy=0.5,
                      height.cluster=0.2,
                      cex=1,
                      cex.main=1,
                      cex.axis=1,
                      cex.poplab=0.6,
                      cex.residlab=0.6,
                      cex.residleg=0.6,
                      cex.names=0.2,
                      cex.ancestral=1,
                      marx=c(3,1),
                      mary.admixture=c(2,3),
                      mary.ancestral=c(1,2),
                      mary.cluster=c(1,1),
                      mary.residual=c(1,1),
                      mary.diff=c(2,2),
                      mary.raw=c(1,1),
                      omni.lwd=2,
                      omni.ylim=NULL,
                      line.above=0,
                      line.below=0,
                      xlim=NULL,
                      showtitles=TRUE,
                      ancestralshow=NULL,
                      ancestralspace=0.5,
                      residual.max=NULL,
                      residual.scale=0.1,
                      residual.showscale=TRUE,
                      ancestral.mean=TRUE,
                      ancestral.names=NULL,
                      ancestral.shownames=TRUE,
                      cluster.srt=0,
                      cluster.names=TRUE,
                      resid.cols=MakeColorRYWGB(),
                      gpref=NULL)
{
    
    theights<-rep(1,length(show))
    theights[show=="Cluster"]<-height.cluster
    if(uselayout) layout(matrix(1:length(show),nrow=length(show)),heights=theights)
    if(is.null(xlim)) xlim=adc$xrange
    if(is.null(gpref)){
        gpref<-list(cex=cex,
                cex.main=cex.main,
                cex.axis=cex.axis,
                cex.poplab=cex.poplab,
                cex.names=cex.names,
                cex.residlab=cex.residlab,
                cex.residleg=cex.residleg,
                cex.ancestral=cex.ancestral,
                omni.lwd=omni.lwd,
                omni.ylim=omni.ylim,
                cluster.names=cluster.names,
                cluster.srt=cluster.srt,
                ancestral.mean=ancestral.mean,
                ancestral.names=ancestral.names,
                ancestral.shownames=ancestral.shownames,
                residual.max=residual.max,
                residual.scale=residual.scale,
                residual.showscale=residual.showscale,
                showtitles=showtitles,
                mary.raw=mary.raw,
                mary.admixture=mary.admixture,
                mary.cluster=mary.cluster,
                mary.ancestral=mary.ancestral,
                mary.diff=mary.diff,
                mary.residual=mary.residual,
                line.above=line.above,
                line.below=line.below,
                mainname=mainname,
                marx=marx,
                popy=popy,
                ancestralshow=ancestralshow,
                ancestralspace=ancestralspace,
                resid.cols=resid.cols,
                xlim=xlim)
    }
    
    for(s in show){
        if(s=="Admixture"){
            plotAdmixtureCoancestryAd(adc,gpref)
        }else if(s=="Cluster"){
            plotAdmixtureCoancestryCluster(adc,gpref)
        }else if(s=="Over"){
            plotAdmixtureCoancestryOver(adc,"over",TRUE,
                                        gpref)
        }else if(s=="Over2"){
            plotAdmixtureCoancestryOver(adc,"over",FALSE,
                                        gpref)
        }else if(s=="Under"){
            plotAdmixtureCoancestryOver(adc,"under",TRUE,
                                        gpref)
        }else if(s=="Under2"){
            plotAdmixtureCoancestryOver(adc,"under",FALSE,
                                        gpref)
        }else if(s=="Raw"){
            plotAdmixtureCoancestryRaw(adc,"data.NbyP",gpref,name="Coancestry")
        }else if(s=="Pred"){
            plotAdmixtureCoancestryRaw(adc,"pred.NbyPatK",gpref,name="Predicted")
        }else if(s=="Same"){
            plotAdmixtureCoancestryRaw(adc,"same.KbyP",gpref,name="Predicted and Present")
        }else if(s=="Omni"){
            plotAdmixtureCoancestryOmni(adc,gpref)
        }else if(s=="Mismatch"){
            plotAdmixtureCoancestryMismatch(adc,gpref)
        }else if(s=="Ancestral"){
            plotAdmixtureCoancestryAncestral(adc,gpref)
        }else if(s=="Residual"){
            plotAdmixtureResidualMatrix(adc,gpref)
        }else{
            stop(paste("Unreckognised plot type :",s))
        }
    }
    return(gpref)
}


#############################################
## SPECIFIC PLOT FUNCTIONS

plotAdmixtureCoancestryAd<-function(adc,
                                    gpref) {

    if(gpref$cex.names>0) {axisnames=TRUE
    }else {axisnames=FALSE ; gpref$cex.names=1}
    par(mar=c(gpref$mary.admixture[1],gpref$marx[1],
              gpref$mary.admixture[2],gpref$marx[2]),
        cex=gpref$cex)
    barplot(t(adc$mix),col=adc$mycols,border=NA,
            space=adc$tspace,las=2,main="",
            cex.main=gpref$cex.main,axes=FALSE,axisnames=axisnames,
            cex.names=gpref$cex.names,xlim=gpref$xlim,xpd=FALSE)
    if(!is.null(gpref$xlim)) {
        par(xpd=FALSE)
        rect(gpref$xlim[2],0,gpref$xlim[2]+100,1,col="white",border=NA)
        rect(gpref$xlim[1]-100,0,gpref$xlim[1],1,col="white",border=NA)
        par(xpd=TRUE)
    }
    
    axis(2,las=1,cex.axis=gpref$cex.axis)
    if(gpref$showtitle)  mtext(paste0(gpref$mainname," K=",adc$K), 3,
                               line=gpref$line.above, cex=gpref$cex.main)
}
                                    
plotAdmixtureCoancestryCluster<-function(adc,
                                         gpref) {
    par(mar=c(gpref$mary.cluster[1],gpref$marx[1],
              gpref$mary.cluster[2],gpref$marx[2]),cex=gpref$cex)
    barplot(t(adc$selfmatrix),col=adc$mycols2,border=NA,
            space=adc$tspace,las=2,
            main="",axes=FALSE,axisnames=FALSE,xlim=gpref$xlim,xpd=FALSE)
    if(gpref$cluster.names) text(adc$popxcentres,gpref$popy,names(adc$poplist),
         adj=0.5,cex=gpref$cex.poplab,srt=gpref$cluster.srt)
    if(!is.null(gpref$xlim)) {
        rect(gpref$xlim[2],0,gpref$xlim[2]+100,1,col="white",border=NA)
        rect(gpref$xlim[1]-100,0,gpref$xlim[1],1,col="white",border=NA)
    }
}

plotAdmixtureCoancestryRaw<-function(adc,val="data.NbyP",
                                         gpref,name="Coancestry") {
    par(mar=c(gpref$mary.raw[1],gpref$marx[1],
              gpref$mary.raw[2],gpref$marx[2]),cex=gpref$cex)
    plotmat=adc[[val]]
    barplot(t(plotmat/rowSums(plotmat)),col=adc$mycols2,border=NA,
            space=adc$tspace,las=2,
            main="",axes=FALSE,axisnames=FALSE,xlim=gpref$xlim,xpd=FALSE)
    if(gpref$showtitles) mtext(paste0(name,", P=",adc$P),3,line=gpref$line.above,cex=gpref$cex.main)
    if(!is.null(gpref$xlim)) {
        rect(gpref$xlim[2],0,gpref$xlim[2]+100,1,col="white",border=NA)
        rect(gpref$xlim[1]-100,0,gpref$xlim[1],1,col="white",border=NA)
    }
    axis(2,las=1,cex.axis=gpref$cex.axis)
}

plotAdmixtureCoancestryAncestral<-function(adc,
                                           gpref) {
    par(mar=c(gpref$mary.ancestral[1],gpref$marx[1],
              gpref$mary.ancestral[2],gpref$marx[2]),cex=gpref$cex)
    if(gpref$ancestral.mean){
        tmp<-rbind(adc$meanpainting.KbyP,adc$coancestry.KbyP)
    }else{
        tmp<-adc$coancestry.KbyP
    }
    tmp[tmp<0]<-0
    tmp<-tmp/rowSums(tmp)
    if(is.null(gpref$ancestralshow)) {
        show=1:dim(tmp)[1]
    }else {
        show=gpref$ancestralshow
        if(gpref$ancestral.mean)  show=c(1,show+1)
    }    
    barplot(t(tmp[show,]),col=adc$mycols2,border=NA,
            space=gpref$ancestralspace,las=2,ylim=c(-0.2,1),
            main="",axes=FALSE,axisnames=FALSE)
    if(gpref$showtitles) mtext(paste0("Ancestral, P=",adc$P," K=",adc$K),3,line=gpref$line.above,cex=gpref$cex.main)

    xleft<-gpref$ancestralspace + (0:length(show))*(1+gpref$ancestralspace)
    tcols=adc$mycols
    tborder=rep(NA,length=adc$K)
    if(gpref$ancestral.mean) {
        tcols<-c("white",tcols)
        tborder=c(1,tborder)
    }
    tcols=tcols[show]
    tborder<-tborder[show]
    if(is.null(gpref$ancestral.names)){
        tnames=paste0("A",1:adc$K)
    }else tnames=gpref$ancestral.names
    if(gpref$ancestral.mean) tnames<-c("Mean",tnames)[show]
    for(i in 1:length(show)) {
        rect(xleft[i],-0.2,xleft[i]+1,-0.05,col=tcols[i],border=tborder[i])
        if(gpref$ancestral.shownames) text(xleft[i]+0.5,-0.125,tnames[i],adj=0.5,cex=gpref$cex.ancestral)
    }
    axis(2,las=1,cex.axis=gpref$cex.axis,at=c(0,0.2,0.4,0.6,0.8,1))

}


plotAdmixtureCoancestryOver<-function(adc,
                                      dir="over",
                                      blackout=TRUE,
                                      gpref) {
    if(dir=="over"){
        if(blackout) {
            top=adc$meandiff.KbyP.over2A
        }else{
            top=adc$meandiff.KbyP.over
        }
        toptext="excess over"
        bottom=adc$diff.KbyP.under
    }else if(dir=="under"){
        if(blackout) {
            top=adc$meandiff.KbyP.under2A
        }else{
            top=adc$meandiff.KbyP.under
        }
        toptext="missing under"
        bottom=adc$diff.KbyP.over
    }else stop(paste("Unreckognised direction:",dir))
    if(blackout) {
        topcolor=adc$mycols2A
    }else{
        topcolor=adc$mycols2
    }
    par(mar=c(gpref$mary.diff[1],gpref$marx[1],
              gpref$mary.diff[2],gpref$marx[2]),cex=gpref$cex)
    barplot(t(top),col=topcolor,border=NA,
        space=adc$tspace,las=2,
        main="",
        cex.main=gpref$cex.main,
        ylim=c(-max(rowSums(bottom)),max(rowSums(top))),
        axes=FALSE,
        axisnames=FALSE,xlim=gpref$xlim,xpd=FALSE)

    axis(2,las=1,cex.axis=gpref$cex.axis)
    barplot(t(-bottom),col=adc$mycols2,border=NA,
            space=adc$tspace,
            las=2,main="",cex.main=gpref$cex.main,add=TRUE,
            axes=FALSE,
            axisnames=FALSE)
    lines(adc$xrange,c(0,0))
    if(gpref$showtitles)     mtext(paste0("P=",adc$P," Predicted Coancestry matrix (",toptext," mean, ",dir,"prediction=",format(sum(adc$dist.KbyP.predfail),digits=4),")"),
          3,line=gpref$line.above,cex=gpref$cex.main)
    if(gpref$showtitles)     mtext(paste("Additional difference to data, total distance =",format(sum(adc$dist.KbyP),digits=4)),
          1,line=gpref$line.below,cex=gpref$cex.main)
}


plotAdmixtureCoancestryMismatch<-function(adc,
                                      gpref) { 
    topcolor=adc$mycols2
    top=adc$diff.KbyP.over
    bottom= adc$diff.KbyP.under
    par(mar=c(gpref$mary.diff[1],gpref$marx[1],
              gpref$mary.diff[2],gpref$marx[2]),cex=gpref$cex)
    barplot(t(top),col=topcolor,border=NA,
        space=adc$tspace,las=2,
        main="",
        cex.main=gpref$cex.main,
        ylim=c(-max(rowSums(bottom)),max(rowSums(top))),
        axes=FALSE,
        axisnames=FALSE,xlim=gpref$xlim,xpd=FALSE)
    axis(2,las=1,cex.axis=gpref$cex.axis)
    barplot(t(-bottom),col=adc$mycols2,border=NA,
            space=adc$tspace,
            las=2,main="",cex.main=gpref$cex.main,add=TRUE,
            axes=FALSE,
            axisnames=FALSE)
    lines(adc$xrange,c(0,0))
    if(gpref$showtitles)     mtext(paste0("P = ",adc$P," overprediction=",format(sum(top),digits=4),")"),
          3,line=gpref$line.above,cex=gpref$cex.main)
    if(gpref$showtitles)     mtext(paste0("P = ",adc$P," underprediction=",format(sum(bottom),digits=4),")"),
          1,line=gpref$line.below,cex=gpref$cex.main)
}

plotAdmixtureCoancestryOmni<-function(adc,
                                      gpref) { 

    topcolor=c(adc$mycols2,adc$mycols2)
    top=cbind(adc$same.KbyP,adc$diff.KbyP.under)
    bottom= adc$diff.KbyP.over
    if(!is.null(gpref$omni.ylim)) {
        ylim=gpref$omni.ylim
    }else{
        ylim=c(-max(rowSums(bottom)),max(rowSums(top)))
    }
    
    par(mar=c(gpref$mary.diff[1],gpref$marx[1],
              gpref$mary.diff[2],gpref$marx[2]),cex=gpref$cex)
    barplot(t(top),col=topcolor,border=NA,
            space=adc$tspace,las=2,
            main="",
            cex.main=gpref$cex.main,
            ylim=ylim,
            axes=FALSE,
            axisnames=FALSE,xlim=gpref$xlim,xpd=FALSE)

    tx0<-cumsum(1+adc$tspace)
    tx<-ty<-rep(NA,max(tx0))
    tx[tx0]<-tx0
    ty[tx0]<-rowSums(adc$same.KbyP)
    tx<-c(NA,tx)
    ty<-c(NA,ty)
    tmp<-c(which(is.na(ty)),Inf)
    tw<-tmp[which(diff(tmp)>1)]
    tx[tw]<-tx[tw+1]-1
    ty[tw]<-ty[tw+1]
    lines(tx,ty,type="S",lwd=gpref$omni.lwd)

    axis(2,las=1,cex.axis=gpref$cex.axis)
    barplot(t(-bottom),col=adc$mycols2,border=NA,
            space=adc$tspace,
            las=2,main="",cex.main=gpref$cex.main,add=TRUE,
            axes=FALSE,
            axisnames=FALSE,xpd=FALSE)
    if(min(ylim)<0) lines(adc$xrange,c(0,0))
    if(gpref$showtitles) {
        mtext(paste0("P = ",adc$P," matching | extra to data"),
              3,line=gpref$line.above,cex=gpref$cex.main)
        mtext(paste0("P = ",adc$P," overprediction = ",
                     format(sum(bottom),digits=4),")"),
              1,line=gpref$line.below,cex=gpref$cex.main)
    }
    if(!is.null(gpref$xlim)) {
        rect(gpref$xlim[2],0,gpref$xlim[2]+100,1,col="white",border=NA)
        rect(gpref$xlim[1]-100,0,gpref$xlim[1],1,col="white",border=NA)
    }

}


plotAdmixtureResidualMatrix<-function(adc,
                                      gpref) {
    # default resid.cols are produced by brewer.pal(11,"Spectral")
    xhack<-adc$tspace
#    xhack[xhack>0]<-1
    xhack<-cumsum(xhack + 1) # Indices of the data ; left out are the indices of the spaces
    
    tdata<-adc$data.NbyP - adc$pred.NbyPatK

    plotdata<-matrix(NA,nrow=max(xhack),ncol=dim(tdata)[2])

    if(is.null(gpref$residual.max)) {
        zlim=c(-1,1)*max(abs(as.numeric(tdata)))
    }else{
        zlim=c(-1,1)*gpref$residual.max
    }
    tdata[tdata>zlim[2]]<-zlim[2]
    tdata[tdata<zlim[1]]<-zlim[1]
    plotdata[xhack,]<-tdata

    par(mar=c(gpref$mary.residual[1]+2,gpref$marx[1],
              gpref$mary.residual[2],gpref$marx[2]),cex=gpref$cex)
    barplot(t(adc$selfmatrix),col="white",border=NA,ylim=c(-gpref$residual.scale,1),
            space=adc$tspace,las=2,
            main="",axes=FALSE,axisnames=FALSE,xlim=gpref$xlim,xpd=FALSE)
    if(gpref$showtitle) mtext("Residuals (data - prediction)",
          3,line=gpref$line.above,cex=gpref$cex.main)

    ty<-dim(plotdata)[2]
    ybounds<-seq(1/ty,1-1/ty,length.out=ty+1)
    image((1:max(xhack))-0.5,
          ybounds,
          plotdata,col=gpref$resid.cols,axes=FALSE,xlab="",
          ylab="",zlim=zlim,add=TRUE)

    if(!is.null(gpref$xlim)) {
        par(xpd=FALSE)
        rect(gpref$xlim[2],0,gpref$xlim[2]+100,1,col="white",border=NA)
        rect(gpref$xlim[1]-100,0,gpref$xlim[1],1,col="white",border=NA)
        par(xpd=TRUE)
    }

    ## colour scale
    plotmin<-min(gpref$xlim)
    plotmax<-max(gpref$xlim)
    plotdist<-plotmax-plotmin
    colindex<-(matrix(seq(zlim[1],zlim[2],
                          length.out=plotdist),
                      ncol=1,nrow=plotdist))
    ## Plot the scale
    if(gpref$residual.showscale){
        image((plotmin:plotmax)-0.5,
              c(-0.5,-1/ty),
              colindex,col=gpref$resid.cols,axes=FALSE,xlab="",
              ylab="",zlim=zlim,add=TRUE)
        axis(1,labels=format(seq(zlim[1],zlim[2],length.out=5),digits=2),
             at=seq(plotmin+0.5,plotmax-0.5,length.out=5),cex.axis=gpref$cex.residleg)
    }
    par(xpd=TRUE)

    ## Draw the boxes
    for(i in 1:adc$P){
        rect(plotmin-plotdist/10,ybounds[i],
             plotmin-plotdist/100,ybounds[i+1],
             col=adc$mycols2[i],border="white")
    }
    if(gpref$cluster.names) text(plotmin-plotdist/50,
         (ybounds[-1] +ybounds[-length(ybounds)])/2,
         colnames(adc$data.NbyP),
         adj=1,cex=gpref$cex.residlab)
    par(xpd=FALSE)

}

