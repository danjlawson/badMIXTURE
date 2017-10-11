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
