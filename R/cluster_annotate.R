#' annotate the celltype for cluster in the seurat.obj
#'
#' @param seurat_object annotation of the seurat.obj
#' @param cluster_markers marker genes of celltypes
#' @param filename the name of the generated file
#' @param width width of the figure
#' @param height height of the figure
#' @param n_col ncol for the vlnplot
#'
#' @importFrom Seurat RenameIdents VlnPlot
#' @importFrom SeuratObject Idents
#' @return a seurat.obj with annotation
#' @export
#'
#' @examples
#' \dontrun{
#' data(common_cluster,package="LIANLAB")
#' input.file <- system.file('extdata','pbmc_1k.RDS',package = 'LIANLABDATA')
#' pbmc_1k <- readRDS(input.file)
#' pbmc_1k <- cluster_annotate( seurat_object = pbmc_1k, cluster_markers = common_cluster)
#' }
cluster_annotate = function(seurat_object,cluster_markers = NULL,filename = "",width = 4,height = 2.5,n_col = NULL){
  if(is.null(cluster_markers)){
    stop("There is no cluster markers supply !!")
  }
  annotation_target <- names(cluster_markers)
  clustermean <- NULL
  seurat_violin <- seurat_object
  for (i in 1:length(cluster_markers)) {
    choose <- subset(seurat_object,features=cluster_markers[[i]])
    df <- as.data.frame(choose@assays$RNA@data)
    score <- apply(df, 2, mean)
    clustermean <- rbind(clustermean,tapply(score,Idents(choose),mean))
    row.names(clustermean)[nrow(clustermean)] <- annotation_target[i]
    seurat_violin$score <- as.numeric(score)
    colnames(seurat_violin@meta.data)[ncol(seurat_violin@meta.data)] <- annotation_target[i]
  }

  annotation <- row.names(clustermean)[apply(clustermean, 2, which.max)]

  for (i in 1:length(annotation)) {
    j <- as.numeric(apply(clustermean, 2, which.max)[i])
    if(clustermean[j,i] < 0.5){
      annotation[i] <- "Unknown"
    }
  }
  if(is.null(n_col)){
    n_col <- ceiling(30/length(levels(seurat_object)))
  }
  n_row <- ceiling(length(annotation_target)/n_col)
  new.cluster.ids  <-  annotation
  names(new.cluster.ids) <- levels(seurat_object)
  seurat_object <- RenameIdents(seurat_object,new.cluster.ids)

  pdf(paste0("annotation_violin_",filename,".pdf"),width = width*n_col,height = height*n_row)
  p <- VlnPlot(seurat_violin,features = annotation_target,pt.size = 0,ncol = n_col)
  print(p)
  dev.off()

  return(seurat_object)
}
