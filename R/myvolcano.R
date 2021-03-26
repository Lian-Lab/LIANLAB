#' Volcano for DE genes.
#'
#' @param df DE gene matrix, including cols named avg_log2FC and p_val_adj
#' @param gene.plot number of tags for genes
#' @param logfc.cutoff cutoff of logFC
#' @param p.cutoff cutoff of p value
#'
#' @importFrom dplyr %>% filter top_n group_by
#' @importFrom grDevices dev.off png
#' @importFrom stats setNames
#' @importFrom ggplot2 ggplot aes geom_point labs theme element_text unit element_line element_blank scale_color_manual element_rect
#' @importFrom ggrepel geom_text_repel
#'
#' @return DE genes
#' @export
#'
#' @examples
#' \dontrun{
#' input.file <- system.file('extdata','DEG.csv',package = 'LIANLABDATA')
#' DEG_matrix <- read.csv(input.file, header = TRUE, row.names = 1)
#' DEG <- myvolcano(df = DEG_matrix, gene.plot = 10, logfc.cutoff = 0.5, p.cutoff = 0.01)
#' }

myvolcano=function(df, gene.plot = 10, logfc.cutoff = 0.5, p.cutoff = 0.01){

  avg_log2FC <- p_val_adj <- group <- NULL

  df[,'color'] = rep('#BEBEBE',times=nrow(df))
  df[,'group'] = rep('#Other',times=nrow(df))
  df$label = rownames(df)
  df$gene = rownames(df)

  up = df %>% filter(avg_log2FC > logfc.cutoff & (p_val_adj < p.cutoff))
  down = df %>% filter(avg_log2FC < -logfc.cutoff & (p_val_adj < p.cutoff))

  df[up$gene,'color'] = '#FF0000'
  df[up$gene,'group'] = 'up'
  df[down$gene,'color'] = '#0000FF'
  df[down$gene,'group'] = 'down'

  df$p_val_adj = -log10(df$p_val_adj)
  df = df[order(df$group,-df$p_val_adj),]
  top = df %>% group_by(group) %>% top_n(n = gene.plot, wt = p_val_adj)

  top = top %>% filter(group !='Other')
  df[!df$gene %in% (top$gene),'label'] = ''

  p <- ggplot(df, aes(avg_log2FC, p_val_adj,color = group)) + geom_point(size=4)+
    labs(x='logFC',y='-Log10 P')+
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_text(margin = unit(c(4, 0, 0, 0), "mm"),size = 18,face = 'bold'),
      axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm"),size = 18,face = 'bold'),
      axis.line = element_line(colour = 'black',size = 1),
      axis.text = element_text(size = 16,face = 'bold',colour = 'black'),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_line(colour = 'grey',size=0.6,linetype = 2),
      axis.ticks = element_line(colour = 'black',size = 1),
      legend.title = element_blank(),
      legend.text = element_text(size = 16,face = 'bold',colour = 'black'),
      legend.position = 'top',legend.key.size = unit(1.5, 'lines')) +
    scale_color_manual(values = setNames(df$color, df$group)) +
    geom_text_repel(aes(label = df$label),color= 'black',fontface = 'bold',size = 6,box.padding = 0.5,point.padding = 0.5,segment.color = 'grey50')


  png(filename = paste0('volcano_',gene.plot,'.png'),width = 700,height = 600)
  print(p)
  dev.off()

  print(p)

  return(top$gene)
}
