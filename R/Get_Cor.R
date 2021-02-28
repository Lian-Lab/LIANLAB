#' Compare the correlations between clusters in one/two seurat.obj, only in mouse or human
#'
#' @param seurat.1 seurat.obj1
#' @param seurat.2 seurat.obj2, could be NULL
#' @param species.1 the species of the seuart.obj1, 'mouse' (default)
#' @param gene To Compare the correlations in the genes you interested
#' @param species.2 the species of the seuart.obj2 (default)
#' @param method One of "pearson" (default), "kendall", or "spearman": can be abbreviated.
#'
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @return the data.frame of the correlations between cluster in two seurat.obj
#' @export
#'
#' @examples
#' \dontrun{
#' input.file <- system.file('extdata','pbmc_1k.RDS',package = 'LIANLABDATA')
#' pbmc_1k <- readRDS(input.file)
#' Get_Cor(pbmc_1k,species.1 = 'mouse')
#' }
Get_Cor = function(seurat.1,seurat.2=NULL,species.1='mouse',gene=NULL,species.2='mouse',method='pearson' ){

  if (!requireNamespace(c('corrplot','stats'), quietly = TRUE)) {
    stop("Please install the package 'corrplot' and 'stats' ")
  }

  celltpye <- NULL

  if (is.null(seurat.2)) {
    #平均矩阵
    seurat.1_data = as.data.frame(seurat.1@assays$RNA@data)
    if (!is.null(gene)) {
      seurat.1_data = seurat.1_data[gene,]
    }
    seurat.1_ave = data.frame(row.names = rownames(seurat.1_data))
    seurat.1$celltpye = seurat.1@active.ident
    metadata = seurat.1@meta.data

    for (i in levels(seurat.1@active.ident)) {

      cell_list = rownames(subset(metadata,subset=celltpye==i))
      seurat.1_ave[,i] = apply(seurat.1_data[,which(colnames(seurat.1_data)%in%cell_list)], 1 , 'mean')

    }

    seurat.1_ave = t(stats::na.omit(t(seurat.1_ave)))
    #相关性分析
    cor1=cor(as.matrix(seurat.1_ave),method = method)
    write.csv(cor1,'correlation.csv')

  }else{
    #对象1
    seurat.1_data = as.data.frame(seurat.1@assays$RNA@data)
    if (!is.null(gene)) {
      seurat.1_data = seurat.1_data[gene,]
    }

    seurat.1_ave = data.frame(row.names = rownames(seurat.1_data))
    seurat.1$celltpye = seurat.1@active.ident
    metadata = seurat.1@meta.data

    for (i in levels(seurat.1@active.ident)) {

      cell_list = rownames(subset(metadata,subset=celltpye==i))
      length(cell_list)
      seurat.1_ave[,i] = apply(seurat.1_data[,which(colnames(seurat.1_data)%in%cell_list)], 1 , 'mean')

    }
    colnames(seurat.1_ave) = paste('seurat.1',colnames(seurat.1_ave))
    colname.1 = colnames(seurat.1_ave)
    seurat.1_ave = t(stats::na.omit(t(seurat.1_ave)))

    #对象2
    seurat.2_data = as.data.frame(seurat.2@assays$RNA@data)
    if (!is.null(gene)) {
      seurat.2_data = seurat.2_data[gene,]
    }
    seurat.2_ave = data.frame(row.names = rownames(seurat.2_data))
    seurat.2$celltpye = seurat.2@active.ident
    metadata = seurat.2@meta.data

    for (i in levels(seurat.2@active.ident)) {

      cell_list = rownames(subset(metadata,subset=celltpye==i))
      length(cell_list)
      seurat.2_ave[,i] = apply(seurat.2_data[,which(colnames(seurat.2_data)%in%cell_list)], 1 , 'mean')

    }
    colnames(seurat.2_ave) = paste('seurat.2',colnames(seurat.2_ave))
    colname.2 = colnames(seurat.2_ave)
    seurat.1_ave = t(stats::na.omit(t(seurat.1_ave)))

    #基因转换
    if(species.1 != species.2){

      if (species.1!="mouse") {

        seurat.2_ave <-  Get_HM_gene_change(seurat.2_ave,species = 'm')
      }else{
        seurat.1_ave <-  Get_HM_gene_change(seurat.1_ave,species = 'm')
      }
    }

    seurat.1_ave = as.data.frame(seurat.1_ave)
    seurat.2_ave = as.data.frame(seurat.2_ave)
    #整合矩阵
    seurat.1_ave$gene = rownames(seurat.1_ave)
    seurat.2_ave$gene = rownames(seurat.2_ave)


    two_ave = merge(seurat.1_ave,seurat.2_ave,by='gene',all=F)

    rownames(two_ave) = two_ave$gene

    two_ave = two_ave[,-which(colnames(two_ave)%in%c('gene'))]
    two_ave[1,1]


    cor=cor(as.matrix(two_ave),method = method)
    cor1 = cor[colname.1,colname.2]
    write.csv(cor1,'correlation.csv')

  }

  #作图
  png(paste0("correlation.png"),width = 40*ncol(cor1)+200, height = 40*nrow(cor1)+100)
  print(corrplot::corrplot(corr = cor1,cl.lim=c(-1,1),tl.cex=1,tl.col='black',
                 col=rev(brewer.pal(n = 10, name = "RdBu"))))
  print(corrplot::corrplot(corr = cor1,add=TRUE, type="lower", method="number",
                 order="original", diag=FALSE,tl.pos="n", cl.pos="n"))
  dev.off()


  png("correlation_pheatmap.png", width = 40*ncol(cor1)+200, height = 40*nrow(cor1)+100)
  p <- pheatmap(cor1,
                treeheight_col = 0,
                treeheight_row = 0,
                cluster_cols = TRUE,cluster_rows = TRUE,
                color =  rev(brewer.pal(n = 11, name = "RdBu")),
                clustering_distance_rows = "correlation",
                clustering_distance_cols = "correlation",
                border_color = 'gray80',
                main= "correlation",fontsize_col = 13,fontsize_row = 13)
  print(p)
  dev.off()

  print(p)
  return(cor1)

}
