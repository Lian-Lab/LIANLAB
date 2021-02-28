#' To get basic drawing, including t-SNE, split figure and featureplot
#'
#' @param object seurat.obj
#' @param filename the name of the generated file
#' @param genenames the genes you interested
#' @param pt.size size of the point
#' @param split_by split the figure by
#' @param width the width of the figure
#' @param height the height of the figure
#' @param ncol the col of the figure
#'
#' @export
#'
#' @examples
#' \dontrun{
#' input.file <- system.file('extdata','pbmc_small.RDS',package = 'LIANLAB')
#' load(input.file)
#' Get_BASIC_drawing(pbmc_small,genenames=c('CD8A','CD4'))
#' }
Get_BASIC_drawing = function(object,filename=NULL,genenames,pt.size = 0.6,split_by ='samples',width = 7, height =6,ncol = 2){

  colorful <- LIANLAB::colorful
  colors = colorful[["colors"]]

  pdf(paste0(filename,"_tSNE",'_lable',".pdf"),width = width, height = height)
  g <- DimPlot(object = object, reduction = 'tsne', label = T,cols=colors,pt.size = 0.6)
  print(g)
  dev.off()

  pdf(paste0(filename,"_tSNE",'_nolable',".pdf"),width = width, height = height)
  g <- DimPlot(object = object, reduction = 'tsne', label = F,cols=colors,pt.size = 0.6)
  print(g)
  dev.off()

  pdf(paste0(filename,"_tSNE",'sample_lable',".pdf"),width = width*ncol, height = height * ceiling(length(unique(object@meta.data[,split_by]))/ncol))
  g <- DimPlot(object = object, reduction = 'tsne',split.by = split_by, label = T,ncol = ncol,cols=colors,pt.size = 0.6)
  print(g)
  dev.off()

  pdf(paste0(filename,"_tSNE",'sample_nolable',".pdf"),width = width*ncol, height = height * ceiling(length(unique(object@meta.data[,split_by]))/ncol))
  g= DimPlot(object = object, reduction = 'tsne',split.by = split_by, label = F,ncol = ncol,cols=colors,pt.size = 0.6)
  print(g)
  dev.off()

  myfeatureplot(object,genenames,filename = filename)

}


#' To get the top DE gene of cluster
#'
#' @param seurat_object seurat.obj
#' @param filename the name of the generated file
#' @param thresh.use Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells
#' @param species the species of the seurat.obj
#' @param all.markers the result of the FindAllMarkers()
#' @param colors the colors of clusters
#'
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom Seurat FindAllMarkers DoHeatmap DotPlot RotatedAxis NoLegend
#'
#' @return a list including two characters,top5 and top20
#' @export
#'
#' @examples
#' \dontrun{
#' input.file <- system.file('extdata','pbmc_1k.RDS',package = 'LIANLABDATA')
#' pbmc_1k <- readRDS(input.file)
#' top_genes <- myfindmarkers(pbmc_1k)
#' }
myfindmarkers = function(seurat_object,filename=NULL,thresh.use=0.25,species="human",all.markers=NULL,colors=NULL){

  cluster <- avg_log2FC <- gene <- NULL
  if(is.null(all.markers)){
    all.markers <- FindAllMarkers(object = seurat_object, only.pos = TRUE, min.pct = 0.25,
                                  logfc.threshold = thresh.use)
  }
  top20 <- all.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)
  top50 <- all.markers %>% group_by(cluster) %>% top_n(50, avg_log2FC)

  write.csv(top20, paste0(filename,"_clustertop20_gene.csv"), row.names = F)
  write.csv(top50, paste0(filename,"_clustertop50_gene.csv"), row.names = F)
  if(species == "huamn"){
    common_markers=c("CD3D","CD3E","CD8A","CD8B","CD4","CD14")
  }else{
    common_markers=c("Cd3d","Cd3e","Cd8a","Cd8b1","Cd4","Cd14")
  }

  if(length(grep(pattern = "TotalSeqC",rownames(seurat_object)))>0){
    antibody_name <- rownames(seurat_object)[grep(pattern = "TotalSeqC",rownames(seurat_object))]
    common_markers <- c(common_markers,antibody_name)
  }
  non_common <- setdiff(unique(as.character(all.markers$gene)),common_markers)
  all.markers <- subset(all.markers,gene%in%non_common)
  top5 <- all.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
  top20 <- all.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)
  heatgene <- unique(as.character(top20$gene))

  if(nrow(seurat_object@assays$RNA@scale.data)==0){
    SeuratObject::DefaultAssay(seurat_object)="RNA"
    seurat_object <- ScaleData(seurat_object)
  }
  ave_cluster <- AverageExpression(seurat_object,assays = "RNA",slot = "data",return.seurat = T)
  if (is.null(colors)) {
    colorful <- LIANLAB::colorful
    colors = colorful[["colors"]]
  }
  pdf(paste0("heatmap_",filename,"_top20.pdf"),width =15, height = 40)
  f1 <- DoHeatmap(ave_cluster, features = heatgene,
               label = TRUE,size=3,assay = "RNA",lines.width = 1,
               group.colors = colors,draw.lines = F)+
    scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))
  print(f1)
  dev.off()
  pdf(paste0("heatmap_",filename,"_top5.pdf"),width =15, height = 12)
  f2 <- DoHeatmap(object = ave_cluster, features = unique(as.character(top5$gene)),
               label = TRUE,size=3,group.colors = colors,draw.lines = F)+
    scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))
  print(f2)
  dev.off()

  png(paste0("heatmap_",filename,"_top20.png"),width = 1500, height = 4000)
  print(f1)
  dev.off()

  png(paste0("heatmap_",filename,"_top5.png"),width = 750/1.5, height = 600/1.5)
  print(f2)
  dev.off()

  pdf(paste0("Dotplot_",filename,"_nolegend.pdf"),width =20, height = 5.18)
  f3 <- DotPlot(object = seurat_object, features = rev(x = unique(top5$gene)), cols = c("blue","firebrick2"),
             dot.scale = 5) + RotatedAxis()+NoLegend()
  print(f3)
  dev.off()

  pdf(paste0("total_heatmap_",filename,"_top5.pdf"),width =15, height = 12)
  f4 <- DoHeatmap(object = seurat_object, features = unique(top5$gene),
               label = TRUE,size=3,group.colors = colors,draw.lines = T)+
    scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))
  print(f4)
  dev.off()

  pdf(paste0("total_heatmap_",filename,"_top20.pdf"),width =15, height = 40)
  f5 <- DoHeatmap(object = seurat_object, features = unique(top20$gene),
               label = TRUE,size=3,group.colors = colors,draw.lines = T)+
    scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))
  print(f5)
  dev.off()

  export <- list(top5=unique(as.character(top5$gene)),top20=unique(as.character(top20$gene)))

  return(export)
}
