#' to show genes like FACS with colors
#'
#' @param seurat_object seurat.obj
#' @param genes the genes to show
#' @param filename the name of the generated file
#' @param color_by to color the cluster, default 'seurat cluster'
#'
#' @export
#'
#' @examples
#' \dontrun{
#' input.file <- system.file('extdata','pbmc_1k.RDS',package = 'LIANLABDATA')
#' pbmc_1k <- readRDS(input.file)
#' Get_Scatterplot1(pbmc_1k,genes=c('CD8A','CD4'))
#' }
Get_Scatterplot1 = function(seurat_object, genes=c('x','y'),filename=NULL,color_by=NULL){
  cluster <- NULL
  x= genes[1]
  y = genes[2]
  list1 = seurat_object@assays$RNA@data[genes,]
  list1 = t(as.matrix(list1))

  if (is.null(color_by)) {
    Ident = as.data.frame(seurat_object@active.ident)
  }else{
    Ident = as.data.frame(seurat_object@meta.data[,color_by])
  }


  Ident = cbind(Ident,list1[,1:2])
  head(Ident)
  colnames(Ident)[1]= 'cluster'
  xlab = colnames(Ident)[2]
  ylab = colnames(Ident)[3]
  colnames(Ident)[2]= 'x'
  colnames(Ident)[3]= 'y'

  p <- ggplot(data = Ident, mapping = aes(x = x , y = y,colour =cluster)) +
    geom_point() + theme_bw() +
    theme(axis.title.x = element_text(face = 'bold',size = 12),
          axis.title.y = element_text(face = 'bold',size = 12),
          axis.text.x = element_text(face = 'bold',size = 12),
          axis.text.y = element_text(face = 'bold',size = 12)) +
    xlab(paste(xlab)) + ylab(paste(ylab))

  pdf(paste0("scatterplot_",filename,".pdf"),width =9, height = 8)
  print(p)
  dev.off()

}


#' To show the corralations between two genes
#'
#' @param mat.object expression matrix as genes in row, samples or cells in column
#' @param seurat_object it also support a seurat_object instead of expression matrix
#' @param log the mat.object whether had log
#' @param genes the genes to show
#' @param filename the name of the generated file
#' @param type Smoothing method (function) to use, accepts either NULL(default) or a character vector, e.g. "lm", "glm", "gam", "loess" or a function
#'
#' @importFrom ggplot2 geom_smooth geom_rug
#' @importFrom ggpubr stat_cor
#'
#' @export
#'
#' @examples
#' \dontrun{
#' input.file <- system.file('extdata','pbmc_1k.RDS',package = 'LIANLABDATA')
#' pbmc_1k <- readRDS(input.file)
#' Get_Scatterplot1(pbmc_1k,genes=c('CD8A','CD4'))
#' }
Get_Scatterplot2 = function(mat.object=NULL,seurat_object=NULL,log=F, genes=c('x','y'),filename=NULL,type=NULL){
  if (is.null(mat.object)) {

    if (is.null(seurat_object)) {

      stop("No correct expression matrix supplied!")

    }else{

    x= genes[1]
    y = genes[2]
    list1 = seurat_object@assays$RNA@data[genes,]

    list1 = t(as.matrix(list1))
    list1 = as.data.frame(list1)

    xlab = colnames(list1)[1]
    ylab = colnames(list1)[2]
    colnames(list1)[1]= 'x'
    colnames(list1)[2]= 'y'

    if (!is.null(type)) {
      p <- ggplot(data = list1,aes(x = x , y = y)) +
        geom_point(color='#9e4aa8') + geom_smooth(color='#fbc08b',size=2.0,method = type) + geom_rug(color='#7fc97e') +
        theme_bw() +theme(axis.text.x = element_text(face = 'bold',size = 12), axis.text.y = element_text(face = 'bold',size = 12)) +
        theme(axis.title.x = element_text(face = 'bold',size = 12),axis.title.y = element_text(face = 'bold',size = 12)) +
        theme(panel.border=element_blank()) +theme(axis.ticks=element_blank())+
        xlab(paste(xlab)) + ylab(paste(ylab))
    }else{
      p <- ggplot(data = list1,aes(x = x , y = y)) +
        geom_point(color='#9e4aa8') + geom_smooth(color='#fbc08b',size=2.0) + geom_rug(color='#7fc97e') +
        theme_bw() +theme(axis.text.x = element_text(face = 'bold',size = 12), axis.text.y = element_text(face = 'bold',size = 12)) +
        theme(axis.title.x = element_text(face = 'bold',size = 12),axis.title.y = element_text(face = 'bold',size = 12)) +
        theme(panel.border=element_blank()) +theme(axis.ticks=element_blank())+
        xlab(paste(xlab)) + ylab(paste(ylab))
    }

    pdf(paste0("scatterplot_",filename,".pdf"),width =8, height = 8)
    print(p)
    dev.off()

    p <- p + stat_cor(data=list1, method = "spearman")
    pdf(paste0("scatterplot_spearman",filename,".pdf"),width =8, height = 8)
    print(p)
    dev.off()

    }

  }else{

    if (log==F) {mat.object = log2(mat.object+1)}

    x= genes[1]
    y = genes[2]

    list1 = mat.object[genes,]

    list1 = t(as.matrix(list1))
    list1 = as.data.frame(list1)

    xlab = colnames(list1)[1]
    ylab = colnames(list1)[2]
    colnames(list1)[1]= 'x'
    colnames(list1)[2]= 'y'


    p <- ggplot(data = list1,aes(x = x , y = y)) +
      geom_point(color='#9e4aa8') + geom_smooth(color='#fbc08b',size=2.0) + geom_rug(color='#7fc97e') +
      theme_bw() +theme(axis.text.x = element_text(face = 'bold',size = 12), axis.text.y = element_text(face = 'bold',size = 12)) +
      theme(axis.title.x = element_text(face = 'bold',size = 12),axis.title.y = element_text(face = 'bold',size = 12)) +
      theme(panel.border=element_blank()) +theme(axis.ticks=element_blank())+
      xlab(paste(xlab)) + ylab(paste(ylab))

    pdf(paste0("scatterplot_",filename,".pdf"),width =8, height = 8)
    print(p)
    dev.off()

    p <- p + stat_cor(data=list1, method = "spearman")
    pdf(paste0("scatterplot_spearman",filename,".pdf"),width =8, height = 8)
    print(p)
    dev.off()
  }

}
