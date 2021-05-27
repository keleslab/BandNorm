#' Performing GAD for each single-cell data
#'
#' This function allows you to calculate GAD value for each gene in each cell.
#' @param path A path to the single-cell Hi-C data. The data format should be the same as bandnorm input.
#' @param genes A data frame containing 4 columns: chrmomsome, start, end, gene name.
#' @param cores Number of cores used for parallel running. Default is 10.
#' @param output The path for outputting scGAD. Each file contains two columns: gene name and gad value.
#' @export
#' @examples
#' #

scGAD = function(path, genes, cores = 10, output){
  colnames(genes) = c("chr", "s1", "s2", "gene_name")
  cl <- makeCluster(cores[1])
  registerDoParallel(cl)
  names = list.files(path)
  output = foreach(k=1:length(names), .packages=c("dplyr", "data.table")) %dopar% {
    paths = list.files(path, full.names = TRUE)
    cell = fread(paths[k])
    gad_score = matrix(0, nrow = nrow(genes), ncol = 2)
    for (i in 1:nrow(genes)){
      temp = cell[V1 == genes[i, ]$chr, ]
      gene_pos = temp[V2 >= genes[i, ]$s1 &
                        V2 <= genes[i, ]$s2 &
                        V4 >= genes[i, ]$s1 &
                        V4 <= genes[i, ]$s2, ]

      up = temp[V2 >= genes[i, ]$s2 &
                  V2 <= (genes[i, ]$s2 + genes[i, ]$s2 - genes[i, ]$s1 + 1) &
                  V4 >= genes[i, ]$s2 &
                  V4 <= (genes[i, ]$s2 + genes[i, ]$s2 - genes[i, ]$s1 + 1), ]

      down = temp[V2 >= (genes[i, ]$s1 - genes[i, ]$s2 + genes[i, ]$s1) &
                    V2 <= genes[i, ]$s1 &
                    V4 >= (genes[i, ]$s1 - genes[i, ]$s2 + genes[i, ]$s1) &
                    V4 <= genes[i, ]$s1, ]

      gad_score[i, ] = c(genes[i, ]$gene_name, ((2 * sum(gene_pos$V5)) / (sum(up$V5) + sum(down$V5))))
    }
    write.table(gad_score, file = paste(output, "/", names[k], sep = ""),
                row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  }
}

