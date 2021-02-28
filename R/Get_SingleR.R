#' Title
#'
#' @param seurat.obj seurat.obj
#' @param species the species of the seurat.obj
#' @param types the types of annotation, 'main' or 'fine'
#'
#' @importFrom dplyr tbl_df select
#' @importFrom SingleR SingleR
#' @importFrom Seurat as.SingleCellExperiment
#' @return the seurat,obj including the celltype annotation
#' @export
#'
#' @examples
#' \dontrun{
#' input.file <- system.file('extdata','pbmc_1k.RDS',package = 'LIANLABDATA')
#' pbmc_1k <- readRDS(input.file)
#' seurat.obj <- Get_singleR( pbmc_1K, species = 'h',types = 'main')
#' }
Get_singleR = function(seurat.obj,species=c('m','h'),types = c('main','fine')) {

  HPCA.se <- BED.se <- IMG.se <- MR.se <- cluster <- NULL
  path = getwd()

  seurat.obj@meta.data$cell.type <- Idents(seurat.obj)
  test <- as.SingleCellExperiment(seurat.obj)

    if (species=='h') {
    input.file <- system.file('extdata','HumanPrimaryCellAtlasData.RData',package = 'LIANLAB')
    load(input.file)
    hp <- HPCA.se
    input.file <- system.file('extdata','BlueprintEncodeData.RData',package = 'LIANLAB')
    load(input.file)
    bp <- BED.se

  }else{
    input.file <- system.file('extdata','ImmGenData.RData',package = 'LIANLAB')
    load(input.file)
    hp <- IMG.se
    input.file <- system.file('extdata','MouseRNAseqData.RData',package = 'LIANLAB')
    load(input.file)
    bp <- MR.se
  }

  if (types == 'main') {
    Anno <- SingleR(test = test,
                    ref = list(HP = hp, BP = bp),
                    labels = list(hp$label.main, bp$label.main),
                    method = "cluster",
                    # de.method="wilcox",
                    # genes="de",
                    de.n = 100,
                    clusters = test$cell.type)
  }else{
    Anno <- SingleR(test = test,
                    ref = list(HP = hp, BP = bp),
                    labels = list(hp$label.fine, bp$label.fine),
                    method = "cluster",
                    # de.method="wilcox",
                    # genes="de",
                    de.n = 100,
                    clusters = test$cell.type)
  }


  #提取需要的细胞分类信息
  Anno$cluster <- rownames(Anno)
  fin <- Anno %>% tbl_df() %>% select(cluster,labels)
  print('Annotation information')
  print(table(Anno$labels))

  #将细胞注释信息重新添加到Seurat对象中去
  new.cluster.ids <- fin$labels
  names(new.cluster.ids) <- levels(seurat.obj)
  seurat.obj <- RenameIdents(seurat.obj, new.cluster.ids)
  setwd(path)
  return(seurat.obj)

}
