#' Featureplot for genes and split by samples, group and so on.
#'
#' @param seurat.obj seurat.obj to draw
#' @param genenames genenames to draw a featureplot
#' @param filename the name of the generated file
#' @param spilt_by whether choose to spilt.by in meta.data
#'
#' @importFrom grDevices pdf
#' @importFrom Seurat FeaturePlot
#' @export
#'
#' @examples
#' \dontrun{
#' input.file <- system.file('extdata','pbmc_1k.RDS',package = 'LIANLABDATA')
#' pbmc_1k <- readRDS(input.file)
#' genenames <- c("MS4A1", "CD79B", "CD79A", "HLA-DRA", "TCL1A", "HLA-DQB1", "HVCN1")
#' myfeatureplot(pbmc_1k, genenames)
#' myfeatureplot(pbmc_1k, genenames, spilt_by = 'group')
#' }
myfeatureplot = function(seurat.obj, genenames, filename = NULL, spilt_by = NULL ){
  if (is.null(spilt_by)) {
    marker_no = length(genenames)
    featureplot_no = marker_no/3
    featureplot_count = 1

    while(featureplot_no > 1) {
      g=FeaturePlot(object = seurat.obj, features = genenames[(3*(featureplot_count-1)+1):(3*featureplot_count)],cols = c("lightblue1","lightgreen","darkorange","firebrick1"), pt.size = 0.5,ncol = 3)

      pdf(paste0("Featureplot_",filename,"_",as.character(featureplot_count),".pdf"),width = 21, height = 14)
      print(g)
      dev.off()
      png(paste0("Featureplot_",filename,"_",as.character(featureplot_count),".png"),width = 2100, height = 1400)
      print(g)
      dev.off()
      featureplot_no=featureplot_no-1
      featureplot_count=featureplot_count+1
    }
    if(featureplot_no > 0){

      g=FeaturePlot(object = seurat.obj, features = genenames[(3*(featureplot_count-1)+1):marker_no],cols = c("lightblue1","lightgreen","darkorange","firebrick1"), pt.size = 0.5,ncol=3)

      pdf(paste0("Featureplot_",filename,"_",as.character(featureplot_count),".pdf"),width = 7*length(featureplot_no), height = 6)
      print(g)
      dev.off()
      png(paste0("Featureplot_",filename,"_",as.character(featureplot_count),".png"),width = 700*length(featureplot_no), height = 600)
      print(g)
      dev.off()
    }
  }else {
    for (i in genenames) {
      num_spilt = unique(seurat.obj@meta.data[,spilt_by])
      n_col = ceiling(sqrt(num_spilt))
      n_row = round(num_spilt/n_col)

      g=FeaturePlot(object = seurat.obj, features = i,cols = c("lightblue1","lightgreen","darkorange","firebrick1"), pt.size = 0.5,ncol=n_col)

      pdf(paste0("Featureplot_",filename,"_",as.character(featureplot_count),".pdf"),width = 7*n_col, height = 6*n_row)
      print(g)
      dev.off()

      png(paste0("Featureplot_",filename,"_",as.character(featureplot_count),".png"),width = 700*n_col, height = 600*n_row)
      print(g)
      dev.off()
    }


  }
}

