#' To draw the average expression of the cluster
#'
#' @param seurat.obj the seurat object
#' @param genes the genes you interested
#' @param split whether to split the heatmap
#' @param colorm the data.frame including clusters(first col) and colors(second col)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' input.file <- system.file('extdata','pbmc_1k.RDS',package = 'LIANLABDATA')
#' pbmc_1k <- readRDS(input.file)
#' genes = c('CD8A','CD3','CD4')
#' Get_ave_heatmap(pbmc_1k,genes,split=T)
#' }
Get_ave_heatmap = function(seurat.obj,genes,split=c(T,F),colorm=NULL,order_col = NULL,cluster_row = F,cluster_col = F){

  type <- NULL
  umi2 = AverageExpression(seurat.obj,assays = 'RNA')

  head(umi2)
  count_data = umi2$RNA
  count_data = count_data[rowSums(count_data)>0,]

  count_data = as.data.frame(t(scale(t(count_data))))
  count_data = count_data[which(rownames(count_data)%in%unique(genes)),]

  annotation_col = data.frame(row.names = levels(seurat.obj),type = levels(seurat.obj))
  annotation_col$type = as.factor(annotation_col$type)
  if (is.null(colorm)) {
    colorful <- LIANLAB::colorful
    colors = colorful[["colors"]]
    colsss = as.character(colors[1:length(levels(seurat.obj))])
    names(colsss) = levels(seurat.obj)

  }else{
    rownames(colorm) = colorm[,1]
    colorm = colorm[rownames(colorm)%in%levels(seurat.obj),]
    colsss = as.character(colorm[,2])
    names(colsss) = colorm[,1]
  }
  if (!is.null(order_col)) {
  count_data <- count_data[,order_col]
  }
  if (split==T) {
    p1 <-  pheatmap(count_data,
                    cluster_rows = cluster_row,
                    color = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")),
                    annotation_colors = list("type"=colsss),
                    annotation_col = annotation_col,
                    cluster_cols = cluster_col,
                    treeheight_row = 0,
                    border_color = 'white',
                    # cellwidth = 30,cellheight = 10,
                    show_rownames=T,
                    cutree_cols = ncol(count_data),
                    # main = i,
                    angle_col = 45,
                    scale = 'none',
                    gaps_col = rep(1:ncol(count_data)),
                    fontsize_row = 8,fontsize_col = 8
                    # fontface="bold",
                    # fontfamily= "Arval"
    )
  }else{
    p1 <- pheatmap(count_data,
                   cluster_rows = cluster_row,
                   color = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")),
                   annotation_colors = list("type"=colsss),
                   annotation_col = annotation_col,
                   cluster_cols = cluster_col,
                   treeheight_row = 0,
                   border_color = 'white',
                   # cellwidth = 30,cellheight = 10,
                   show_rownames=T,
                   # main = i,
                   angle_col = 45,
                   scale = 'none',
                   fontsize_row = 8,fontsize_col = 8
                   # fontface="bold",
                   # fontfamily= "Arval"
    )
  }
  print(p1)


}
