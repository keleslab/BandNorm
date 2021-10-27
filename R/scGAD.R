#' Performing scGAD for each single-cell data
#'
#' This function allows you to calculate GAD value for each gene in each cell.
#' @param path A path to the single-cell Hi-C data. The data format should be the same as in the bandnorm hic_df input.
#' @param hic_df You can also load dataframe containing all Hi-C data to here, but it is not recommended, since in scGAD, we are dealing with high resolution matrices, and it will consume a lot of memory if we load it directly.
#' @param genes A data frame containing 5 columns: chrmomsome, start, end, strand, gene name.
#' @param cores Number of cores used for parallel running. Default is 4.
#' @param threads Number of threads for fread function, default is 8.
#' @param depthNorm Whether to normalize the sequencing depth effect. Default is TRUE.
#' @export
#' @examples
#' #

scGAD = function(path = NULL, hic_df = NULL, genes, depthNorm = TRUE, cores = 4, threads = 8){
  setDTthreads(threads)
  discardCounts = max(genes$s2 - genes$s1)
  chr = genes$chr
  s1_low <- ifelse(genes$strand == "+", genes$s1 - 1000, genes$s1)
  s2 <- ifelse(genes$strand == "-", genes$s2, genes$s2 + 1000)
  if (is.null(hic_df)){
    cl <- makeCluster(cores[1])
    registerDoParallel(cl)
    names = list.files(path)
    paths = list.files(path, full.names = TRUE)
    output = foreach(k=1:length(names), .packages=c("dplyr", "data.table", "matrixStats"), .combine = 'cbind') %dopar% {
      cell = fread(paths[k], select = c(1, 2, 4, 5))
      colnames(cell) = c("V1", "V2", "V4", "v5")
      cell = cell[abs(V4 - V2) <= discardCounts]
      setkey(cell, V1)
      
      gad_score = rep(NA, nrow(genes))
      idx = list()
      pchr = chr[1]
      temp = cell[J(pchr)]
      for (i in 1:nrow(genes)){
        cchr = chr[i]
        if (cchr != pchr) {
          temp = cell[J(cchr)]
        }
        
        gad_score[i] = sum2(.subset2(temp, 4)[temp$V2 %between% c(s1_low[i], s2[i]) &
                                                temp$V4 %between% c(s1_low[i], s2[i])])
        pchr = chr[i]
      }
      gad_score
    }
  } else{
    cl <- makeCluster(cores[1])
    registerDoParallel(cl)
    hic_df = setDT(hic_df)
    names = unique(hic_df$cell)
    output = foreach(k=1:length(names), .packages=c("dplyr", "data.table", "matrixStats"), .combine = 'cbind') %dopar% {
      setkey(hic_df, cell)
      tempcell = hic_df[J(names[k])]
      tempcell = tempcell[abs(binA - binB) <= discardCounts]
      setkey(tempcell, chrom)
      gad_score = rep(NA, nrow(genes))
      idx = list()
      pchr = chr[1]
      temp = tempcell[J(pchr)]
      for (i in 1:nrow(genes)){
        cchr = chr[i]
        if (cchr != pchr) {
          temp = tempcell[J(cchr)]
        }
        
        gad_score[i] = sum2(.subset2(temp, 5)[temp$binA %between% c(s1_low[i], s2[i]) &
                                                temp$binB %between% c(s1_low[i], s2[i])])
        pchr = chr[i]
      }
      gad_score
    }
  }
  
  colnames(output) = names
  rownames(output) = genes$gene_name
  output = output[rowSums(output) > 0, ]
  output = output[!is.na(rowSums(output)), ]
  
  if (depthNorm) {
    output = t(t(output) / colSums(output)) * 1e04
  }
  GAD = (output - rowMeans(output))/sqrt(rowVars(output))
  GAD
}


#' Performing Projection for scGAD on other single-cell data
#'
#' This function allows you to project scGAD value on various assays.
#' @param DataList A list of matrices containing all assays. Each matrix should have a name that corresponds to the assay. 
#' @param doNorm A vector of boolean. Each entries determines whether each matrix should be normalized.
#' @export
#' @examples
#' #

runProjection = function(DataList, doNorm, cellTypeList){
  geneList = lapply(DataList, rownames)
  common_gene = Reduce(intersect, geneList)
  selectGene = function(assay, genes = common_gene){
    assay[genes, ]
  }
  DataList = lapply(DataList, selectGene)
  nameAssays = names(DataList)
  GADList = list()
  for (i in 1:length(nameAssays)){
    s = CreateSeuratObject(DataList[[i]])
    if (doNorm[i]){
      s = NormalizeData(object = s)
      all.genes <- rownames(s)
      s <- ScaleData(s, features = all.genes)
    }else {
      s@assays$RNA@scale.data = DataList[[i]]
    }
    Idents(s) = cellTypeList[[i]]
    s$method = nameAssays[i]
    DefaultAssay(s) = "RNA"
    GADList[[i]] = s
  }
  
  GADList <- lapply(X = GADList, FUN = function(x) {
    x <- FindVariableFeatures(x, selection.method = "mean.var.plot", nfeatures = 2500)
  })
  features <- SelectIntegrationFeatures(object.list = GADList)
  suppressWarnings(getAnchors <- FindIntegrationAnchors(object.list = GADList, anchor.features = features,
                                                        k.filter = 50))
  suppressWarnings(getCombined <- IntegrateData(anchorset = getAnchors))

  DefaultAssay(getCombined) <- "integrated"

  suppressWarnings(getCombined <- ScaleData(getCombined, verbose = FALSE))
  getCombined <- RunPCA(getCombined, npcs = 5, verbose = FALSE)
  getCombined <- RunUMAP(getCombined, reduction = "pca", dims = 1:5)
  getCombined
}

