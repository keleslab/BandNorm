#' Performing scGAD for each single-cell data
#'
#' This function allows you to calculate GAD value for each gene in each cell.
#' @param path A path to the single-cell Hi-C data. The data format should be the same as in the bandnorm hic_df input.
#' @param hic_df You can also load dataframe containing all Hi-C data to here, but it is not recommended, since in scGAD, we are dealing with high resolution matrices, and it will consume a lot of memory if we load it directly.
#' @param genes A data frame containing 5 columns: chrmomsome, start, end, strand, gene name.
#' @param cores Number of cores used for parallel running. Default is 4.
#' @param threads Number of threads for fread function, default is 8.
#' @param depthNorm Whether to normalize the sequencing depth effect. Default is TRUE.
#' @param binPair Use bin pair or valid reads. The former is faster since it is already binned, but the latter will be more accurate. Default is TRUE.
#' @param res The resolution of the data. Used only when binPair = TRUE. Default is 10000.
#' @param format The format of the valid pairs. "short", "medium", "long", "4DN" are supported. See https://github.com/aidenlab/juicer/wiki/Pre.
#' @export
#' @examples
#' #

scGAD = function(path = NULL, hic_df = NULL, genes, depthNorm = TRUE, cores = 4, threads = 8, binPair = TRUE, format = "short", res = 10000){
  setDTthreads(threads)
  discardCounts = max(genes$s2 - genes$s1)
  genes$s1 <- ifelse(genes$strand == "+", genes$s1 - 1000, genes$s1)
  genes$s2 <- ifelse(genes$strand == "-", genes$s2, genes$s2 + 1000)
  colnames(genes) = c("chr", "start", "end", "strand", "names")
  genes = makeGRangesFromDataFrame(genes, keep.extra.columns=TRUE)

  if (is.null(hic_df)){
    if (binPair){
      names = list.files(path)
      paths = list.files(path, full.names = TRUE)
      getCount = function(k){
        cell = fread(paths[k], select = c(1, 2, 4, 5))
        colnames(cell) = c("V1", "V2", "V4", "V5")
        cell = cell[abs(V4 - V2) <= discardCounts]
        GInt = GenomicInteractions(GRanges(cell$V1,
                                           IRanges(cell$V2, width = res)),
                                   GRanges(cell$V1,
                                           IRanges(cell$V4, width = res)),
                                   counts = cell$V5)
        one <- overlapsAny(anchorOne(GInt), genes)
        two <- overlapsAny(anchorTwo(GInt), genes)
        x.valid <- GInt[one & two]
        hits <- list()
        hits$one <- findOverlaps(anchorOne(x.valid), genes, select = "first")
        hits$two <- findOverlaps(anchorTwo(x.valid), genes, select = "first")
        counts = data.table(reads = x.valid$counts[hits[[1]] == hits[[2]]], pos = hits$one[hits$one == hits$two])
        tabulated <- unique(counts$pos)
        counts <- setDT(counts)[,.(reads = sum(reads)), by = 'pos']$reads
        dat = data.table(names = genes[unique(tabulated)]$names, counts = counts)
        colnames(dat) = c("names", names[k])
        dat
      }
      cl <- makeCluster(cores)
      clusterEvalQ(cl, {library(data.table)
        library(GenomicInteractions)
      })
      output <- parLapply(cl, 1:length(names), getCount)
      output = suppressMessages(Reduce(full_join, output))
      output[is.na(output)] = 0
    }
    else {
      names = list.files(path)
      paths = list.files(path, full.names = TRUE)
      getCount = function(k){
        cell = fread(paths[k])
        if (format == "short"){
          cell = cell[cell$V2 == cell$V6, ]
          cell = cell[, c(2, 3, 7)]
        }else if (format == "medium"){
          cell = cell[cell$V3 == cell$V7, ]
          cell = cell[, c(3, 4, 8)]
        }else if (format == "long"){
          cell = cell[cell$V2 == cell$V6, ]
          cell = cell[, c(2, 3, 7)]
        }else if (format == "4DN"){
          cell = cell[cell$V2 == cell$V4, ]
          cell = cell[, c(2, 3, 5)]
        }
        colnames(cell) = c("V1", "V2", "V4")
        cell = cell[abs(V4 - V2) <= discardCounts]
        GInt = GenomicInteractions(GRanges(cell$V1,
                                           IRanges(cell$V2, width = res)),
                                   GRanges(cell$V1,
                                           IRanges(cell$V4, width = res)),
                                   counts = 1)
        one <- overlapsAny(anchorOne(GInt), genes)
        two <- overlapsAny(anchorTwo(GInt), genes)
        x.valid <- GInt[one & two]
        hits <- list()
        hits$one <- findOverlaps(anchorOne(x.valid), genes, select = "first")
        hits$two <- findOverlaps(anchorTwo(x.valid), genes, select = "first")
        counts = data.table(reads = x.valid$counts[hits[[1]] == hits[[2]]], pos = hits$one[hits$one == hits$two])
        tabulated <- unique(counts$pos)
        counts <- setDT(counts)[,.(reads = sum(reads)), by = 'pos']$reads
        dat = data.table(names = genes[unique(tabulated)]$names, counts = counts)
        colnames(dat) = c("names", names[k])
        dat
      }
      cl <- makeCluster(cores)
      clusterEvalQ(cl, {library(data.table)
        library(GenomicInteractions)
      })
      output <- parLapply(cl, 1:length(names), getCount)
      output = suppressMessages(Reduce(full_join, output))
      output[is.na(output)] = 0
    }
  } else{
    hic_df = setDT(hic_df)
    names = unique(hic_df$cell)
    getCount = function(k){
      library(data.table)
      cell = hic_df[hic_df$cell == names[k], ]
      GInt = GenomicInteractions(GRanges(cell$chrom,
                                         IRanges(cell$binA, width = res)),
                                 GRanges(cell$chrom,
                                         IRanges(cell$binB, width = res)),
                                 counts = cell$count)
      one <- overlapsAny(anchorOne(GInt), genes)
      two <- overlapsAny(anchorTwo(GInt), genes)
      x.valid <- GInt[one & two]
      hits <- list()
      hits$one <- findOverlaps(anchorOne(x.valid), genes, select = "first")
      hits$two <- findOverlaps(anchorTwo(x.valid), genes, select = "first")
      counts = data.table(reads = x.valid$counts[hits[[1]] == hits[[2]]], pos = hits$one[hits$one == hits$two])
      tabulated <- unique(counts$pos)
      counts <- aggregate(reads ~ pos, data = counts, FUN = sum)$reads
      dat = data.table(names = genes[unique(tabulated)]$names, counts = counts)
      colnames(dat) = c("names", names[k])
      dat
    }
    cl <- makeCluster(cores)
    clusterEvalQ(cl, {library(data.table)
                      library(GenomicInteractions)
                     })
    output <- parLapply(cl, 1:length(names), getCount)
    output = suppressMessages(Reduce(full_join, output))
    output[is.na(output)] = 0
  }
  finNames = output$names
  output = as.matrix(output[, -1])
  rownames(output) = finNames
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
  print(1)
  features <- SelectIntegrationFeatures(object.list = GADList)
  suppressWarnings(getAnchors <- FindIntegrationAnchors(object.list = GADList, anchor.features = features,
                                                        k.filter = 50))
  suppressWarnings(getCombined <- IntegrateData(anchorset = getAnchors, k.weight = 80))

  DefaultAssay(getCombined) <- "integrated"

  suppressWarnings(getCombined <- ScaleData(getCombined, verbose = FALSE))
  getCombined <- RunPCA(getCombined, npcs = 15, verbose = FALSE)
  getCombined <- RunUMAP(getCombined, reduction = "pca", dims = 1:5)
  getCombined
}

