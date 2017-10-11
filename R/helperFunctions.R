#' @title Find and test mixPainter binaries
#'
#' @description
#' badMIXTURE uses chromosome painting to test whether a given mixture is
#' a good explanation of the data. mixPainter is a version of
#' finestructure/ChromoPainter that makes this simple.
#'
#' This function gives you the command-line location of a version of mixPainter
#' on your system. It looks for "mixPainter" in your path, or asks you to download
#' the mixPainter zip package and put the results somewhere it can find.
#' Note that it
#' is possible that no versions work on your system; currently we do support:
#' * linux (gcc 5.3.0, compiled on Scientific Linux release 6.4)
#' * Mac (g++-5 5.2.0, OSX Sierra 10.12)
#' We will make every attempt to provide a working version for other mac or
#' linux versions but can make no guarantee that we can get it working for you.
#' In that case you would be able to script your own interface with chromopainter.
#'
#' It will also print the version information of mixPainter.
#' 
#' @return a character string of the location of mixPainter
#' @export
mixPainter<-function(){

    tmp=suppressWarnings(system("which mixPainter",intern=T))
    if(length(tmp)>0) {
        version=system("mixPainter -V")
        print("Found mixPainter in your path")
        return(tmp)
    }
    
    endings=c(".mac",".linux","")
    attempt0=sapply(endings,function(x){
        tmp=system.file("bin",paste0("mixPainter",x), package = "badMIXTURE")
        if(nchar(tmp)>=1) return(TRUE)
        FALSE
    })
    if(all(!attempt0)){
        stop(paste("You don't have mixPainter. Try downloading it from https://people.maths.bris.ac.uk/~madjl/finestructure/mixPainter.zip and putting the contents into",
                   system.file("bin", package = "badMIXTURE")))
    }
    endings=endings[attempt0]
    attempt=sapply(endings,function(x){
        cmd=paste(system.file("bin",paste0("mixPainter",x), package = "badMIXTURE"),"-h")
        tmp=suppressWarnings(system(cmd,intern=T,ignore.stderr=TRUE))
        if(length(tmp)>1) return(TRUE)
        FALSE
    })
    if(all(attempt==FALSE)) stop("ERROR: Cannot find a mixPainter executable that works on your system!")
    ret=    system.file("bin",paste0("mixPainter",
                             names(which(attempt)[1])),
                package="badMIXTURE")

    version=system(paste(ret,"-V"))
    ret
}


#' @title Aggregate a matrix by taking column means
#'
#' @description
#' Take means over columns in a matrix, by grouping all columns listed in poplist
#' e.g. if mat is M*N matrix and poplist is length K, returns a M*K matrix
#' the names of poplist are used to assign names to the returned matrix
#' 
#' @param mat A matrix
#' @param poplist list of which members of the original mat are in each of the new columns
#' @export
matColMeans<-function(mat,poplist)
{
	res<-matrix(0,nrow=dim(mat)[1],ncol=length(poplist))
	colindex<-lapply(poplist,function(x){which(colnames(mat)%in%x)})
	res<-t(apply(mat,1,function(x){
		sapply(colindex,function(y){mean(x[y])})
	}))
	colnames(res)<-names(poplist)
	rownames(res)<-rownames(mat)
	res
}


#' @title Aggregate a matrix by taking column sums
#'
#' @description
#' Sums over columns in a matrix, by grouping all columns listed in poplist
#' e.g. if mat is M*N matrix and poplist is length K, returns a M*K matrix
#' the names of poplist are used to assign names to the returned matrix
#' 
#' @param mat A matrix of dimension M by N
#' @param poplist list of which members (column identifiers) of the original mat are in each of the new columns, as returned by \code{\link{idsToList}}.
#' @return A matrix on dimension M by K.
#' @export
matColSums<-function(mat,poplist)
{
	res<-matrix(0,nrow=dim(mat)[1],ncol=length(poplist))
	colindex<-lapply(poplist,function(x){which(colnames(mat)%in%x)})
	res<-t(apply(mat,1,function(x){
		sapply(colindex,function(y){sum(x[y])})
	}))
#        if(dim(res)[1]==1) res<-t(res)
	colnames(res)<-names(poplist)
	rownames(res)<-rownames(mat)
	res
}



#' @title Convert an ID table into a population list
#'
#' @description
#'
#' Takes a data frame converting individual identifiers into cluster ids, and returns this same structure. The cluster ids (and individual ids within each cluster id) appear in the same order they appeared in the data.
#' 
#' @param ids A data frame containing 2 or more columns. The first are individual identifiers, the second are population identifiers. Any further columns are ignored. NB: Chromopainter excludes individuals in the third column that have a 0; you are advised to do this before passing to this function.
#' @return A list of populations, named  by the population identifiers. Each is a chartacter vector of the individual identifiers in each population.
#' @export
idsToList<-function(ids){
    tpops<-unique(ids[,2])
    poplistunsrt<-lapply(tpops,function(x){
        ids[ids[,2]==x,1]
    })
    names(poplistunsrt)<-tpops
    poplistunsrt
}

#' @title Reorder a Q matrix into the order that we 
#'
#' @description
#'
#' A function to reorder the data to make the populations appear in the same order as the latent clusters. Takes an admixture matrix and reorders the rows according to a desired vector of population names. It then reorders the columns so that the latent clusters are ordered similarly to the labels.
#' 
#' @param Q A matrix of admixture estimates for each individual
#' @param ids A data frame containing 2 or more columns. The first are individual identifiers, the second are population identifiers. Any further columns are ignored. NB: Chromopainter excludes individuals in the third column that have a 0; you are advised to do this before passing to this function.
#' @param poporder A list of the order that we want populations to appear in. If NULL, populations will be ordered as they appear in the data
#' @param thresh A threshold for the admixture memberships when counting how to order columns of Q
#' @return A reordered matrix of the Q values
#' @export
reorderQ=function(Q,ids,poporder=NULL, thresh=0.95){
    ##
    if(is.null(poporder)) poporder=unique(ids[,2])
#    if(!(all(poporder%in%unique(ids[,2])))) stop("Not all populations in poporder are represented in ids!")
    if(!(all(unique(ids[,2])%in%poporder))) {
        warning("Not all ids are represented in the populations provided in poporder!")
        poporder
    }
    if(!(all(rownames(Q)%in%ids[,1]))) stop("Not all individuals in Q are found in ids!")
    if(!(all(ids[,1]%in%rownames(Q)))) stop("Not all individuals in ids are found in Q!")

    mapstatelist<-lapply(poporder,function(x){
        ids[(ids[,2]==x),1]
    })
    names(mapstatelist)<-poporder
    Q<-Q[unlist(mapstatelist),]
    tq=apply(Q,2,function(x){
        tw=which(x > thresh)
        if(length(tw)==0) tw=which.max(x)
        x=mean(tw)
    })
    Q<-Q[,order(tq)]
    Q
}


#' @title Load a chunkcounts matrix from chromopainter
#'
#' @description
#'
#' Read chromopainter output and present it in the form that badMIXTURE needs; i.e. a matrix normalised to have rowsums = 1, which is the default behaviour
#' 
#' @param filename Name of the chunkcounts file to be loaded
#' @param skip number of lines to skip at the top of the file. Set to 1 if there is a CFACTOR header.
#' @return A matrix of normalised chunkcounts
#' @seealso readQ
#' @examples
#' \dontrun{
#' recent_ariL=readChunkcounts((system.file("extdata",
#' "Recent_admix.pruned_cp_pop_linked.chunkcounts.out", package = "badMIXTURE")))
#' }
#' @export
readChunkcounts=function(filename,skip=1){
    x=as.matrix(utils::read.table(filename,header=T,row.names=1,skip=skip))
    x=x/rowSums(x)
    x
}

#' @title Load an ADMIXTURE Q file
#'
#' @description
#'
#' Read ADMIXTURE output, and assign row names to it. Optionally, reorder it (if poporder is specified)
#' 
#' @param filename Name of the ADMIXTURE Q file to read
#' @param ids An ids dataframe; individual ids in column 1, population ids in column 2
#' @param poporder The order that we we want populations to be displayed in. Should correspond to the unique labels in the ids file. No reordering done if not provided
#' @param thresh A threshold for counting an individual as a member of a population, for use in reordering columns of Q to match the order in poporder
#' @examples
#' \dontrun{
#' ## Admixture Q matrices don't have row names. Read an ID file for this:
#' recent_ariids=read.table(system.file("extdata",
#'   "Recent_admix.ids", package = "badMIXTURE"),as.is=T)
#' ## Now read the Q matrix
#' recent_ariQ=readQ(system.file("extdata",
#'   "Recent_admix.pruned.11.Q", package = "badMIXTURE"),recent_ariids)
#'
#' ## Now we want to reorder the individuals, specifying the order that
#' we would like our populations to appear
#' poporder=paste0("Pop",c(7,5,6,13,4,8,9,11:12,1:3))
#' recent_ariQ=readQ(system.file("extdata",
#'   "Recent_admix.pruned.11.Q", package = "badMIXTURE"),recent_ariids,poporder)
#' }
#' @return A matrix of normalised chunkcounts
#' @seealso reorderQ
#' @export
readQ=function(filename,ids,poporder=NULL,thresh=0.95){
    x=as.matrix(utils::read.table(filename))
    rownames(x)=ids[,1]
    if(!is.null(poporder)){
        x=reorderQ(x,ids,poporder,thresh=thresh)
    }
    x
}


###########################################
## COLORS


#' @title Make a colour palette moving from Yellow to Red to Purple
#'
#' @description
#' Make a colour palette moving from Yellow to Red to Purple.
#' 
#' @seealso MakeColorRYWGB, MakeColorWYRPB
#' @param colby The step size of the difference between successive colours. Each colour transition is characterised by moving one or move RGB colours from 0 to 1 (or vice versa) in steps of colby. Default: 0.05
#' @param final A colour to place at the end, expressed as an RGB numeric vector. Default: NULL, meaning place no additional colour at the end.
#' @export
MakeColorYRP<-function(colby=0.05,final=NULL){
tmp<-c(grDevices::rgb(1,seq(1,0,-colby),0),grDevices::rgb(1,0,seq(colby,1,colby)),grDevices::rgb(seq(1-colby,0,-colby),0,1.0))
	if(is.null(final)) return(tmp)
	c(tmp,grDevices::rgb(final[1],final[2],final[3]))
}

#' @title Make a colour palette moving from White to Yellow to Red to Purple to Blue
#'
#' @description
#' Make a colour palette moving from White to Yellow to Red to Purple to Blue.
#' 
#' @seealso MakeColorYRP, MakeColorRYWGB
#' @param colby The step size of the difference between successive colours. Each colour transition is characterised by moving one or move RGB colours from 0 to 1 (or vice versa) in steps of colby. Default: 0.05
#' @export
MakeColorWYRPB<-function(colby=0.05){
    if(length(colby)<5)colby<-rep(colby,each=5)
    c(grDevices::rgb(1,1,seq(1,0,-colby[1])),
      grDevices::rgb(1,seq(1,0,-colby[2]),0),
      grDevices::rgb(1,0,seq(colby[3],1,colby[3])),
      grDevices::rgb(seq(1-colby[4],0,-colby[4]),0,1.0),
      grDevices::rgb(0,0,seq(1-colby[5],0,-colby[5])))
}

#' @title Make a colour palette moving from Red to Yellow to White to Green to Blue
#'
#' @description
#' Make a colour palette moving from Red to Yellow to White to Green to Blue.
#' 
#' @seealso MakeColorYRP, MakeColorWYRPB
#' @param colby The step size of the difference between successive colours. Each colour transition is characterised by moving one or move RGB colours from 0 to 1 (or vice versa) in steps of colby. Default: 0.05
#' @export
MakeColorRYWGB<-function(colby=0.05){
    if(length(colby)<6)colby<-rep(colby,each=6)
    c(grDevices::rgb(seq(0.5,1,colby[1]/2),seq(0,0,colby[1]),0),
      grDevices::rgb(1,seq(0,1,colby[2]),0),
      grDevices::rgb(1,1,seq(colby[3],1,colby[3])),
      grDevices::rgb(seq(1,colby[4],-colby[4]),1,1),
      grDevices::rgb(0,seq(1,colby[5],-colby[5]),1),
      grDevices::rgb(0,0,seq(1,0.5+colby[6],-colby[6]/2)))
}


