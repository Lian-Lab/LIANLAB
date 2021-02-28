#' To get seurat.obj easily
#'
#' @param object gene expression matrix
#' @param minfeature Include genes where at least this many features are detected.
#' @param mincell Include cells where at least this many features are detected.
#' @param species the species of the matrix
#' @param minFeatureRNA the min cut off the genes count
#' @param maxFeatureRNA the max cut off the genes count
#' @param Percent.mt the max cut off of the percetnage of the mitochondrial gene
#' @param is.tumor whether the sample is tumor
#' @param projectname the name of the project
#' @param dim.FindNeighbors number of dims to find neighbors
#' @param dims.TSNE number of dims to t-SNE
#'
#' @return the seurat.obj including PCA, t-SNE
#' @export
#'
#' @examples
#' \dontrun{
#' input.file <- system.file('extdata','pbmc_small.csv',package = 'LIANLAB')
#' pbmc_small <- read.csv(input.file, header = T, row.names = 1)
#'
#' seurat.obj = CreateSeuratObject(counts = pbmc_small, min.features = 200,min.cells = 10)
#' dim(Abseq.seurat)
#' See_mRNA(Abseq.seurat)
#' data[["percent.mt"]] <- PercentageFeatureSet(Abseq.seurat, pattern = "^MT.")  # MT(human)/mt(mouse)
#' VlnPlot(Abseq.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#' table(Abseq.seurat[["percent.mt"]]<10)
#'
#' seurat.obj = Get_Seurat(pbmc_small, minfeature = 200, mincell=10,species='h')
#' }
#'
Get_Seurat = function(object,minfeature=200,mincell=10,species=c('h','m'),minFeatureRNA=NULL,maxFeatureRNA=NULL,Percent.mt=NULL,is.tumor=F,projectname = NULL,dim.FindNeighbors=1:30,dims.TSNE=1:30){
  nFeature_RNA <- percent.mt <- NULL
  if (is.null(minFeatureRNA)) {
    minFeatureRNA <- 200
  }
  if (is.null(maxFeatureRNA)) {
    maxFeatureRNA <- 4000
  }
  if (is.tumor==F) {
    if (is.null(Percent.mt)) {
      Percent.mt <- 25
    }
  }else{
    if (is.null(Percent.mt)) {
      Percent.mt <- 100
    }
  }

  if (is.null(projectname)) {
    projectname <- 'SCseq'
  }
  minFeatureRNAs <- minFeatureRNA
  maxFeatureRNAs <- maxFeatureRNA
  Percent.mts <- Percent.mt


  object <-  CreateSeuratObject(counts = object,project = projectname, min.features = minfeature, min.cells = mincell)
  print(dim(object))

  if (species=='m') {
    object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^mt")
  }else{
    object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT")
  }


  VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

  object <- subset(object, subset = ((nFeature_RNA > minFeatureRNAs) & (nFeature_RNA < maxFeatureRNAs) & (percent.mt < Percent.mts)) )
  object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)

  all.genes <- rownames(object)
  object <- ScaleData(object, features = all.genes)
  object <- RunPCA(object, features = VariableFeatures(object = object))
  DimPlot(object = object, reduction = "pca")

  ElbowPlot(object = object,ndims = 50)
  object <- FindNeighbors(object = object, dims = dim.FindNeighbors)
  object <- FindClusters(object = object, resolution = 0.40)
  object <- RunTSNE(object, dims = dims.TSNE)
  object <- RunUMAP(object, reduction = "pca", dims = dims.TSNE)

  return(object)

}

#' To get seurat.obj for integrated
#'
#' @param object gene expression matrix
#' @param minfeature Include genes where at least this many features are detected.
#' @param mincell Include cells where at least this many features are detected.
#' @param species the species of the matrix
#' @param minFeatureRNA the min cut off the genes count
#' @param maxFeatureRNA the max cut off the genes count
#' @param Percent.mt the max cut off of the percetnage of the mitochondrial gene
#' @param is.tumor whether the sample is tumor
#' @param projectname the name of the project
#'
#' @importFrom Seurat DimPlot ElbowPlot FindNeighbors FindClusters RunTSNE RunUMAP VariableFeatures NormalizeData ScaleData RunPCA CreateSeuratObject PercentageFeatureSet
#' @return the seurat.obj
#' @export
#'
#' @examples
#' \dontrun{
#' input.file <- system.file('extdata','pbmc_small.csv',package = 'LIANLABDATA')
#' pbmc_small <- read.csv(input.file, header = T, row.names = 1)
#'
#' seurat.obj = CreateSeuratObject(counts = pbmc_small, min.features = 200,min.cells = 10)
#' dim(Abseq.seurat)
#' See_mRNA(Abseq.seurat)
#' data[["percent.mt"]] <- PercentageFeatureSet(Abseq.seurat, pattern = "^MT.")  # MT(human)/mt(mouse)
#' VlnPlot(Abseq.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#' table(Abseq.seurat[["percent.mt"]]<10)
#'
#' seurat.obj = Get_for_integrated(pbmc_small, minfeature = 200, mincell = 10, species = 'h')
#' }
#'
Get_for_integrated = function(object, minfeature = 200,mincell = 10,species = c('h','m'),minFeatureRNA=NULL,maxFeatureRNA=NULL,is.tumor=F,Percent.mt=NULL,projectname = NULL){

  nFeature_RNA <- percent.mt <- NULL
  if (length(minFeatureRNA)==0) {
    minFeatureRNA <- 200
  }
  if (length(maxFeatureRNA)==0) {
    maxFeatureRNA <- 4000
  }
  if (is.tumor==F) {
    if (is.null(Percent.mt)) {
      Percent.mt <- 25
    }
  }else{
    if (is.null(Percent.mt)) {
      Percent.mt <- 100
    }
  }
  if (length(projectname)==0) {
    projectname <- 'SCseq'
  }
  minFeatureRNAs <- minFeatureRNA
  maxFeatureRNAs <- maxFeatureRNA
  Percent.mts <- Percent.mt


  object <- CreateSeuratObject(counts = object,project = projectname, min.features = minfeature,min.cells = mincell)
  print(dim(object))
  if (species=='m') {
    object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^mt.")
  }else{
    object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT.")
  }


  VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

  object <- subset(object, subset = ((nFeature_RNA > minFeatureRNAs) & (nFeature_RNA < maxFeatureRNAs) & (percent.mt < Percent.mts)) )
  ##标准化
  object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)

  return(object)

}


#' To integrate different batch seurat.obj together
#'
#' @param umi.list the seurat list
#' @importFrom Seurat FindIntegrationAnchors IntegrateData
#' @return the integrated seurat.obj
#' @export
#'
#' @examples
#' \dontrun{
#' seurat.int = Get_integrated(list(seurat.obj1,seurat.obj2,seurat.obj3))
#' }
#'
Get_integrated = function(umi.list) {

  if (requireNamespace("future", quietly = TRUE)) {

    future::plan("multiprocess", workers = 1)
    umi <- FindIntegrationAnchors(object.list = umi.list, dims = 1:40)

    future::plan("multiprocess", workers = 15)
    umi.integrated <- IntegrateData(anchorset = umi, dims = 1:30)
    SeuratObject::DefaultAssay(umi.integrated) <- "integrated"
    future::plan("multiprocess", workers = 1)
    umi.integrated <- ScaleData(umi.integrated, verbose = T)
    future::plan("multiprocess", workers = 10)
    umi.integrated <- RunPCA(umi.integrated, verbose = T)

    DimPlot(object = umi.integrated, reduction = "pca")

    ElbowPlot(object = umi.integrated, ndims = 50)
    umi.integrated <- FindNeighbors(object = umi.integrated, dims = 1:30)
    umi.integrated <- FindClusters(object = umi.integrated, resolution = 0.6)
    colnames(umi.integrated@meta.data)

    ###tsne
    umi.integrated <- RunTSNE(object = umi.integrated, dims = 1:30,check_duplicates = FALSE)
    umi.integrated <- RunUMAP(umi.integrated, reduction = "pca", dims = 1:30)
    future::plan("multiprocess", workers = 1)

  }else{
    umi <- FindIntegrationAnchors(object.list = umi.list, dims = 1:40)

    umi.integrated <- IntegrateData(anchorset = umi, dims = 1:30)
    SeuratObject::DefaultAssay(umi.integrated) <- "integrated"
    umi.integrated <- ScaleData(umi.integrated, verbose = T)
    umi.integrated <- RunPCA(umi.integrated, verbose = T)

    DimPlot(object = umi.integrated, reduction = "pca")

    ElbowPlot(object = umi.integrated,ndims = 50)
    umi.integrated <- FindNeighbors(object = umi.integrated, dims = 1:30)
    umi.integrated <- FindClusters(object = umi.integrated, resolution = 0.6)

    ###tsne
    umi.integrated <- RunTSNE(object = umi.integrated, dims = 1:30,check_duplicates = FALSE)
    umi.integrated <- RunUMAP(umi.integrated, reduction = "pca", dims = 1:30)
  }

  return(umi.integrated)

}

