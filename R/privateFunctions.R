#' @import grDevices
#' @import graphics
#' @import stats

#' @title  Relabel a dendrogram based on IDs
#'
#' @description
#' Look at each tip of popdend, and find all individuals in the ids structure that have this population. Replace labels with a ; separated character string of the individuals in that population
#' 
#' @param popdend A dendrogram
#' @param ids An N by 3 data frame consisting of: column 1: row names (for both the data and the mix). column 2: the cluster menbership that created the groups in dataraw (with the column names in dataraw as the values). column 3: inclusion (0 for absent, 1 for present).
#' @param relabel A function that takes the names in ids[,1] and returns the names that will be found in popdend.

popDendRelabelMembers<-function(popdend,ids,relabel=function(x)x){
    tdend<-dendrapply(popdend,function(x){
        if(is.leaf(x)){
            attr(x,"label")<-paste(
                relabel(ids[ids[,2]==attr(x,"label"),1]),
                collapse=";")
        }
        x
    })
    tdend
}


#' @title Aggregate a matrix of dimension M to dimension K by taking the best cut of a dendrogram
#'
#' @description
#' Cut a tree for the rows of a matrix to get K tips
#' Make the popdata that has this new set of rows by taking row means (or whatever combine is set to do)
#' 
#' @param popdataraw A matrix
#' @param popdend A dendrogram relating the rows of popdataraw
#' @param K A height at which to cut the dendrogram
#' @param combine A function to combine columns of popdataraw
#' @param simplify Whether to simplify the labels of the dendrogram after it has been cut
aggregrateDataForK<-function(popdataraw,popdend,K, combine=matColMeans,simplify=TRUE){
        uch<-uniqueCutHeights(popdend)
        tcutdend<-cut(popdend,uch[as.character(K)])
        tcutlabels<-lapply(tcutdend$lower,labels)
        if(simplify){
            tcutlabels<-strsplit(as.character(tcutlabels),";")
            tcutlabels<-vapply(tcutlabels,
                               function(x)strsplit(x,";")[[1]],FUN.VALUE="character")
        }else{
            names(tcutlabels)=sapply(tcutlabels,
                               function(x)paste(x,collapse=";"))
        }

        popdataraw.cut<-t(combine(t(popdataraw),tcutlabels))
        popdataraw.cut
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

#' @title Find all of the unique heights of a dendrogram at which it can be cut
#'
#' @description
#' Obtain all cuts of a dendrogram, return this in a useful format
#' 
#' @param tdend A dendrogram object
uniqueCutHeights<-function(tdend)
{
    theights<-dendrogramHeights(tdend)
    theights<-theights-c(0,diff(theights)/2)
    names(theights)<-rev(1+1:length(theights))
    theights
}

#' @title Find all of the unique heights of a dendrogram at which it can be cut
#'
#' @description
#' Obtain all cuts of a dendrogram
#' 
#' @param tdend A dendrogram object
dendrogramHeights<-function(tdend){
    ### NOTE: better implementation in dendextend: heights_per_k.dendrogram(dend)
    nodeHeight<-function(x){
        if(is.leaf(x)) {
            return(attr(x,"height"))#attr(x,"height"))
        }else{
            return(c(attr(x,"height"),nodeHeight(x[[1]]),nodeHeight(x[[2]])))
        }
    }
    ret=sort(unique(unlist(sapply(tdend,nodeHeight))))
    c(ret,(attr(tdend,"height")+ret[length(ret)])/2)
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

###########################################
## COLORS


#' #' @title Make a colour palette moving from Yellow to Red to Purple
#'
#' @description
#' Make a colour palette moving from Yellow to Red to Purple.
#' 
#' 
#' @param colby The step size of the difference between successive colours. Each colour transition is characterised by moving one or move RGB colours from 0 to 1 (or vice versa) in steps of colby. Default: 0.05
#' @param final A colour to place at the end, expressed as an RGB numeric vector. Default: NULL, meaning place no additional colour at the end.
MakeColorYRP<-function(colby=0.05,final=NULL){
    require(grDevices)
tmp<-c(rgb(1,seq(1,0,-colby),0),rgb(1,0,seq(colby,1,colby)),rgb(seq(1-colby,0,-colby),0,1.0))
	if(is.null(final)) return(tmp)
	c(tmp,rgb(final[1],final[2],final[3]))
}

#' @title Make a colour palette moving from White to Yellow to Red to Purple to Blue
#'
#' @description
#' Make a colour palette moving from White to Yellow to Red to Purple to Blue.
#' 
#' 
#' @param colby The step size of the difference between successive colours. Each colour transition is characterised by moving one or move RGB colours from 0 to 1 (or vice versa) in steps of colby. Default: 0.05
MakeColorWYRPB<-function(colby=0.05){
    ## makes white/yellow/red/purple/black colour scheme, adding rgb(final) if not null
    require(grDevices)
    if(length(colby)<5)colby<-rep(colby,each=5)
    c(rgb(1,1,seq(1,0,-colby[1])),
      rgb(1,seq(1,0,-colby[2]),0),
      rgb(1,0,seq(colby[3],1,colby[3])),
      rgb(seq(1-colby[4],0,-colby[4]),0,1.0),
      rgb(0,0,seq(1-colby[5],0,-colby[5])))
}

#' @title Make a colour palette moving from Red to Yellow to White to Green to Blue
#'
#' @description
#' Make a colour palette moving from Red to Yellow to White to Green to Blue.
#' 
#' 
#' @param colby The step size of the difference between successive colours. Each colour transition is characterised by moving one or move RGB colours from 0 to 1 (or vice versa) in steps of colby. Default: 0.05
MakeColorRYWGB<-function(colby=0.05){
    ## makes white/yellow/red/purple/black colour scheme, adding rgb(final) if not null
    if(length(colby)<6)colby<-rep(colby,each=6)
    c(rgb(seq(0.5,1,colby[1]/2),seq(0,0,colby[1]),0),
      rgb(1,seq(0,1,colby[2]),0),
      rgb(1,1,seq(colby[3],1,colby[3])),
      rgb(seq(1,colby[4],-colby[4]),1,1),
      rgb(0,seq(1,colby[5],-colby[5]),1),
      rgb(0,0,seq(1,0.5+colby[6],-colby[6]/2)))
}



#' @title Obtain a colour palette for a matrix using multi-dimensional scaling
#'
#' @description
#' Makes a colour for each row of the data by embedding it in 4 dimensions, "r,g,b,alpha"
#' 
#' @param mydata data for which colours are required; each row is an observation
#' @param colorder Mapping of MDS directions to RGBalpha channels;allows relabelling of colour directions and changing polarity (if -ve)
#' @param colmax the maximum a colour can be (for preventing white)
#' @param alphamin the minimum the alpha can be (for preventing colours from being too faded)
rgbDistCols<-function(mydata, colorder=c(1,2,3,4),colmax=0.8, alphamin=0.5){
    ##     require("MASS") # Only requred for the isoMDS implementation, disabled
    d <- dist(mydata) # euclidean distances between the rows
    fit <- cmdscale(d,eig=TRUE, k=4) # k is the number of dim
#    fit <- isoMDS(d, k=4) # k is the number of dim
    emb<-as.matrix(fit$points)
    row.names(emb)<-NULL
    emb<-apply(emb,2,function(x){
        (order(x)-1)/(length(x)-1)
    })
    mymap<-function(x,sign=1){
        if(sign>0) return(x)
        return(1-x)
    }
    return(rgb(mymap(emb[,abs(colorder[1])],colorder[1])*colmax,
               mymap(emb[,abs(colorder[2])],colorder[2])*colmax,
               mymap(emb[,abs(colorder[3])],colorder[3])*colmax,
               (mymap(emb[,abs(colorder[4])],colorder[4]))*(1-alphamin)+alphamin))
}
