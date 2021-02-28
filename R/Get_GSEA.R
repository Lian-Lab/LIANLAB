#' Title
#'
#' @param object DEG data.frame or matrix, including col named 'avg_logFC'
#' @param species species of DEG
#' @param category geneset in GSEA, including 'c2'(default),'c3','c4','c5','c6','c7','hallmark'
#' @param filename the name of the generated file
#' @param p.valueCutoff cut off the p value
#'
#' @importFrom clusterProfiler bitr enricher GSEA
#' @importFrom dplyr select
#' @importFrom enrichplot gseaplot2 gseaplot
#' @importFrom utils data
#' @return result of the gsea
#' @export
#'
#' @examples
#' \dontrun{
#' input.file <- system.file('extdata','DEG.csv',package = 'LIANLABDATA')
#' markers <- read.csv(input.file,header = T, row.names = 1)
#' result <- Get_GSEA( markers,species = 'm',category='c2')
#' }
Get_GSEA = function(object,species=c('m','h'),category=NULL,filename=NULL,p.valueCutoff=0.05){

  geneset <- symbol <- mouse.symbol <- msigdf.human <- msigdf.mouse <- category_code <- colorful <- NULL
  if (is.null(category)) {
    category = 'c2'
  }
  genelist=data.frame(rownames(object),object$avg_logFC)
  colnames(genelist)[1]='SYMBOL'
  colnames(genelist)[2]='avg_logFC'


  if (species=='h') {
    msigdf.human <- msigdf::msigdf.human
    c2 <- msigdf.human %>% filter(category_code == category) %>% select(geneset, symbol ) %>% as.data.frame
  }else{
    msigdf.mouse <- msigdf::msigdf.mouse
    c2 <- msigdf.mouse %>% filter(category_code == category) %>% select(geneset, mouse.symbol) %>% as.data.frame
  }
  colorful <- LIANLAB::colorful
  colors = colorful[["colors"]]

  de <- genelist[,1]
  de = as.character(de)

  if (species=='h') {
    de = bitr(de,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = 'org.Hs.eg.db')
  }else{
    de = bitr(de,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = 'org.Mm.eg.db')
  }

  x <- enricher(de, TERM2GENE = c2, pvalueCutoff = p.valueCutoff)
  write.csv(x@result,paste0(filename,'_','enricher.csv'))

  de1 <- as.numeric(genelist$avg_logFC)
  names(de1) = genelist$SYMBOL
  de1 = sort(de1,decreasing = T)

  y <- GSEA(de1, TERM2GENE = c2, verbose=FALSE, pvalueCutoff = p.valueCutoff)
  write.csv(y@result,paste0(filename,'_','GSEA.csv'))

  png(paste0("GSEAplot2_",y@result$ID[1],".png"),width =600, height = 500)
  p <- gseaplot2(y, 1,title = y@result$ID[1],base_size=10,subplots = 1:4)
  print(p)
  dev.off()

  png(paste0("GSEAplot2_pvalue",y@result$ID[1],".png"),width =600, height = 500)
  p <- gseaplot2(y, 1,title = y@result$ID[1],base_size=10,subplots = 1:4,pvalue_table =TRUE)
  print(p)
  dev.off()

  png(paste0("GSEAplot2_top5.png"),width =600, height = 500)
  p <- gseaplot2(y, 1:5,base_size=10,subplots = 1:4,color =colors)
  print(p)
  dev.off()

  png(paste0("gseaplot_",y@result$ID[1],".png"),width =600, height = 500)
  p <- gseaplot(y, 1,base_size=5,title =y@result$ID[1], subplots = 1:5,color.line ='firebrick',color='#DAB546')
  print(p)
  dev.off()

  return(y)
}
