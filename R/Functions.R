#' Calculate FPKM/RPKM
#'
#' Function to calculate fragments/reads per kilobase per million using a gene count matrix
#'
#' @param mat data frame count matrix with rows as genes and columns as samples
#' @param geneLength data frame with two columns: GeneId corresponding to rownames(mat) and geneLength with the exon lengths of each gene
#' @param libSize numeric vector of length ncol(mat) containing (in order) the library sizes for each sample in mat. Leave NULL to use colSums(mat).
#' @return data frame with FPKM/RPKM normalized gene counts
#' @export

fpkm<-function(mat, geneLength, libSize=NULL){
  message("Calculating per million scale factors...")
  scales<-c()
  if(is.null(libSize)){
    scales<-colSums(mat)/1000000
  } else {
    scales<-libSize/1000000
  }
  message("Calculating Reads per million (RPM)...")
  #pb<-txtProgressBar(min=0, max=length(scales), style=3)
  #for(i in 1:length(scales)){
  #  mat[,i]<-mat[,i]/scales[i]
    #setTxtProgressBar(pb, i)
  #}
  mat<-sweep(mat, 2, scales, "/")
  mat<-mat[order(rownames(mat)),]
  geneLength<-geneLength[order(geneLength$GeneId),]
  if(!all.equal(rownames(mat), geneLength$GeneId)){
    stop("Check genes in count matrix and gene length match")
  }
  message("Normalizing to gene length to get RPKM/FPKM...")
  #pb<-txtProgressBar(min=0, max=nrow(geneLength), style=3)
  #for(i in 1:nrow(geneLength)){
  #  mat[i,]<-mat[i,]/(geneLength$Length[i]/1000)
  #  setTxtProgressBar(pb, i)
  #}

  mat<- 1000 * sweep(mat, 1, geneLength$Length, "/")

  return(mat)
}

#' Calculate TPM
#'
#' Function to calculate Transcripts per million using a gene count matrix
#'
#' @param mat data frame count matrix with rows as genes and columns as samples
#' @param geneLength data frame with two columns: GeneId corresponding to rownames(mat) and geneLength with the exon lengths of each gene
#' @param libSize numeric vector of length ncol(mat) containing (in order) the library sizes for each sample in mat. Leave NULL to use colSums(mat).
#' @return data frame with TPM normalized gene counts
#' @export

tpm<-function(mat, geneLength, libSize=NULL){
  message("Calculting per million scale factors...")
  scales<-c()
  if(is.null(libSize)){
    scales<-colSums(mat)/1000000
  } else {
    scales<-libSize/1000000
  }

  mat<-mat[order(rownames(mat)),]
  geneLength<-geneLength[order(geneLength$GeneId),]
  if(nrow(mat)>nrow(geneLength)){
    mat<-mat[which(rownames(mat) %in% geneLength$GeneId),]
  }
  if(nrow(mat)<nrow(geneLength)){
    geneLength<-geneLength[which(geneLength$GeneId %in% rownames(mat)),]
  }
  if(!all.equal(rownames(mat), geneLength$GeneId)){
    stop("Check genes in count matrix and gene length match")
  }
  message("Calculating reads per kilobase (RPK)...")
  #pb<-txtProgressBar(min=0, max=nrow(geneLength), style=3)
  #for(i in 1:nrow(geneLength)){
  #  mat[i,]<-mat[i,]/(geneLength$Length[i]/1000)
  #  setTxtProgressBar(pb, i)
  #}
  mat <- 1000 * sweep(mat, 1, geneLength$Length, "/")
  message("Calculating TPM...")
  #pb<-txtProgressBar(min=0, max=length(scales),style = 3)
  #for(i in 1:length(scales)){
  #  mat[,i]<-mat[,i]/scales[i]
  #  setTxtProgressBar(pb, i)
  #}
  mat <- sweep(mat, 2, scales, "/")
  return(mat)
}

fpkm2tpm<-function(mat, libSize=NULL){
  if(is.null(libSize)){
    return(sweep(mat, 2, colSums(mat)/1000000, "/"))
  } else {
    return(sweep(mat, 2, libSize/1000000, "/"))
  }
}

fpkm2raw<-function(mat, geneLength, libSize){
  mat<-mat[order(rownames(mat)),]
  geneLength<-geneLength[order(geneLength$GeneId),]
  if(nrow(mat)>nrow(geneLength)){
    mat<-mat[which(rownames(mat) %in% geneLength$GeneId),]
  }
  if(nrow(mat)<nrow(geneLength)){
    geneLength<-geneLength[which(geneLength$GeneId %in% rownames(mat)),]
  }
  if(!all.equal(rownames(mat), geneLength$GeneId)){
    stop("Check genes in count matrix and gene length match")
  }
  mat<-sweep((mat/1000), 1, geneLength$Length, "*")
  scales<-libSize/1000000
  mat<-sweep(mat, 2, scales, "*")
  return(mat)
}

tpm2raw<-function(mat, geneLength, libSize){
  mat<-mat[order(rownames(mat)),]
  geneLength<-geneLength[order(geneLength$GeneId),]
  if(nrow(mat)>nrow(geneLength)){
    mat<-mat[which(rownames(mat) %in% geneLength$GeneId),]
  }
  if(nrow(mat)<nrow(geneLength)){
    geneLength<-geneLength[which(geneLength$GeneId %in% rownames(mat)),]
  }
  if(!all.equal(rownames(mat), geneLength$GeneId)){
    stop("Check genes in count matrix and gene length match")
  }
  scales<-libSize/1000000
  mat<-sweep(mat,2, scales, "*")
  mat<-sweep((mat/1000), 1, geneLength$Length, "*")
  return(mat)
}
