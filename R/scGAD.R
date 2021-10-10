#' Performing scGAD for each single-cell data
#'
#' This function allows you to calculate GAD value for each gene in each cell.
#' @param path A path to the single-cell Hi-C data. The data format should be the same as bandnorm input.
#' @param hic_df You can also load dataframe containing all Hi-C data to here, but it is not recommended, since in scGAD, we are dealing with high resolution matrices, and it will consume a lot of memory if we load it directly.
#' @param genes A data frame containing 4 columns: chrmomsome, start, end, gene name.
#' @param cores Number of cores used for parallel running. Default is 10.
#' @param depthNorm Whether to normalize the sequencing depth effect. Default is TRUE.
#' @export
#' @examples
#' #

scGAD = function(path = NULL, hic_df = NULL, genes, depthNorm = TRUE, cores = 25){
  discardCounts = max(genes$s2 - genes$s1)

  if (is.null(hic_df)){
    cl <- makeCluster(cores[1])
    registerDoParallel(cl)
    names = list.files(path)
    paths = list.files(path, full.names = TRUE)
    output = foreach(k=1:length(names), .packages=c("dplyr", "data.table"), .combine = 'cbind') %dopar% {
      cell = fread(paths[k])
      cell = cell[abs(cell$V4 - cell$V2) <= discardCounts, ]
      gad_score = rep(NA, nrow(genes))
      for (i in 1:nrow(genes)){
        temp = cell[V1 == genes[i, ]$chr, ]
        gad_score[i] = sum(temp[V2 >= genes[i, ]$s1 - 10000 &
                                  V2 <= genes[i, ]$s2 &
                                  V4 >= genes[i, ]$s1 - 10000 &
                                  V4 <= genes[i, ]$s2, ]$V5)
      }
      gad_score
    }
  } else{
    names = unique(hic_df$cell)
    output = c()
    for (i in 1:length(names)){
      tempcell = hic_df[hic_df$cell == names[i], ]
      tempcell = tempcell[abs(tempcell$binA - tempcell$binB) <= discardCounts, ]
      gad_score = rep(NA, nrow(genes))
      for (i in 1:nrow(genes)){
        temp = tempcell[tempcell$chrom == genes[i, ]$chr, ]
        gad_score[i] = sum(temp[temp$binA >= genes[i, ]$s1 - 10000 &
                                  temp$binA <= genes[i, ]$s2 &
                                  temp$binB >= genes[i, ]$s1 - 10000 &
                                  temp$binB <= genes[i, ]$s2, ]$count)
      }
      output = cbind(output, gad_score)
    }
  }

  colnames(output) = names
  rownames(output) = genes$gene_name
  
  if (depthNorm) {
    output = t(t(output) / colSums(output)) * 1e04
  }
  GAD = (output - rowMeans(output))/sqrt(rowVars(output))
  GAD
}




