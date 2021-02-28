#' Batch/Specify cell population for GO and KEGG analyse
#'
#' @param marker result of FindAllMarkers(), including cols named avg_logFC, p_val_adj, gene
#' @param spe_cluster the clusters you interested
#' @param species m/h
#' @param logfc.cutoff cut off the logFC
#' @param p.cutoff cut off the p value
#'
#' @export
#'
#' @examples
#' \dontrun{
#' input.file <- system.file('extdata','DEG.csv',package = 'LIANLABDATA')
#' markers <- read.csv(input.file,header = T, row.names = 1)
#' Get_GK1( markers_top,species = 'm')
#' Get_GK1( markers_top,spe_cluster = 'Naive CD4 T',species = 'm')
#' }
Get_GK1 = function(marker,spe_cluster = NULL,species = c('m','h'),logfc.cutoff=0.1,p.cutoff=0.01){
  cluster <- avg_logFC <- p_val_adj <- gene <- NULL
  if ( !is.null(spe_cluster)) {
    if(length(spe_cluster)==1){
      clustermarker = marker[which(marker$cluster == spe_cluster),]
      up = clustermarker %>% filter(avg_logFC > logfc.cutoff & p_val_adj < p.cutoff)
      down = clustermarker %>% filter(avg_logFC < -logfc.cutoff & p_val_adj < p.cutoff)
      up_gene=up$gene
      down_gene=down$gene
      if (length(up_gene)==0) {print('There is no up gene')}
      if (length(down_gene)==0) {print('There is no down gene')}
      if (length(up_gene)==0 & length(down_gene)==0) {
        print(paste0('cluster_',spe_cluster,' is wrong.'))
      }else{
        myGO_KEGG(species = species,filename = paste0('cluster_',spe_cluster), up_gene = up_gene,down_gene = down_gene)
      }

    }else{
      for (j in spe_cluster ) {
        clustermarker = marker[which(marker$cluster == j),]
        up = clustermarker %>% filter(avg_logFC > logfc.cutoff & p_val_adj < p.cutoff)
        down = clustermarker %>% filter(avg_logFC < -logfc.cutoff & p_val_adj < p.cutoff)
        up_gene=up$gene
        down_gene=down$gene
        if (length(up_gene)==0) {print('There is no up gene')}
        if (length(down_gene)==0) {print('There is no down gene')}
        if (length(up_gene)==0 & length(down_gene)==0) {
          print(paste0('cluster_',j,' is wrong.'))

        }else{
          myGO_KEGG(species = species,filename = paste0('cluster_',j), up_gene = up_gene,down_gene = down_gene)

        }

      }
    }
  }else{
    for (i in  unique(marker$cluster)) {
      clustermarker = marker[which(marker$cluster == i),]
      up = clustermarker %>% filter(avg_logFC > logfc.cutoff & p_val_adj < p.cutoff)
      down = clustermarker %>% filter(avg_logFC < -logfc.cutoff & p_val_adj < p.cutoff)

      up_gene=up$gene
      down_gene=down$gene

      if (length(up_gene)==0) {print('There is no up gene')}
      if (length(down_gene)==0) {print('There is no down gene')}
      if (length(up_gene)==0 & length(down_gene)==0) {
        print(paste0('cluster_',i,' is wrong.'))
      }else{
        myGO_KEGG(species = species,filename = paste0('cluster_',i), up_gene = up_gene,down_gene = down_gene)
      }
    }
  }
}



#' GO and KEGG analyse of the result between two clusters
#'
#' @param marker result of FindMarkers() or other,including cols named avg_logFC, p_val_adj, gene
#' @param species m/h
#' @param logfc.cutoff cut off the logFC
#' @param p.cutoff cut off the p value
#'
#' @export
#'
#' @examples
#' \dontrun{
#' input.file <- system.file('extdata','DEG.csv',package = 'LIANLABDATA')
#' markers <- read.csv(input.file,header = T, row.names = 1)
#' Get_GK2( markers_top,species = 'm')
#' }
Get_GK2 = function(marker,species = c('m','h'),logfc.cutoff=0.1,p.cutoff=0.01){
  avg_logFC <- p_val_adj <- gene <- NULL
  marker$gene = rownames(marker)
  up = marker %>% filter(avg_logFC > logfc.cutoff & p_val_adj < p.cutoff)
  down = marker %>% filter(avg_logFC < -logfc.cutoff & p_val_adj < p.cutoff)
  up_gene=up$gene
  down_gene=down$gene
  if (length(up_gene)==0) {print('There is no up gene')}
  if (length(down_gene)==0) {print('There is no down gene')}
  myGO_KEGG(species = species,filename = paste0('GK_'), up_gene = up_gene,down_gene = down_gene)

}


#' Integrated GO/KEGG pathways gene
#'
#' @param DATA the GO/KEGG result, including cols named geneID and Description
#' @param keyword the keywood of the pathways you interested
#'
#' @importFrom stringr str_split
#' @return the genes in the pathways you choose
#' @export
#'
#' @examples
#' \dontrun{
#' input.file <- system.file('extdata','GO_enrich.csv',package = 'LIANLAB')
#' markers <- read.csv(input.file,header = T, row.names = 1)
#' Get_GK_genes(markers)
#' Get_GK_genes(markers, keyword = c('DNA'))
#' }
Get_GK_genes = function(DATA,keyword=NULL){
  geneID <- Description <- NULL

  row_list = list()
  if (!is.null(keyword)) {
    j =1
    for (i in keyword) {
      ROW = grepl(pattern = paste0('*',i),DATA$Description)
      row_list[[j]] = rownames(DATA[ROW,])
      j = j+1
    }
    nROWs = row_list[[1]]
    for (i in row_list) {
      nROWs = union(nROWs,i)
    }

    DATA = DATA[nROWs,]
    write.csv(DATA,'select.csv')
  }

  if (nrow(DATA)!=0) {
    gene = data.frame(DATA$geneID)
    gene_lists = list()
    i = 1
    while(i<nrow(gene)+1){
      gene_lists[[i]] = as.character(unlist(str_split(gene[i,1], '/', n = Inf, simplify = FALSE)))
      i = i + 1
    }
    GENES = gene_lists[[1]]
    for (i in gene_lists) {
      GENES = union(GENES,i)
    }
    return(GENES)
  }else{
    print('There is no set of genes that satisfy the keywords')
  }

}
