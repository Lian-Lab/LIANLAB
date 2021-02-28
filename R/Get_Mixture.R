#' To change the genes from human/mouse to mouse/human
#'
#' @param Objects the matrix to change gene
#' @param species species of the matrix genes
#' @importFrom biomaRt getLDS
#' @return The genes after changed
#' @export
#'
#' @examples
#' \dontrun{
#' input.file <- system.file('extdata','pbmc_small.csv',package = 'LIANLABDATA')
#' matrix.obj <- read.csv(input.file,header = T, row.names = 1)
#' matrix.obj <- Get_HM_gene_change( matrix.obj, species = 'h')
#' character.obj <- c('Cd4','Cd8a','Cd8b1','Cd3e')
#' character.obj <- Get_HM_gene_change( character.obj, species = 'm')
#' }
Get_HM_gene_change = function(Objects,species=c('h','m')){
  input.file <- system.file('extdata','human_gene_ensembl.rds',package = 'LIANLAB')
  human <- readRDS(input.file)

  input.file <- system.file('extdata','mouse_gene_ensembl.rds',package = 'LIANLAB')
  mouse <- readRDS(input.file)

if (class(Objects)=="data.frame"|class(Objects)=="matrix") {
  if(species=='m'){
    gene_mouse=row.names(Objects)
    print("Converting mouse gene symbols to human...")
    gene_human=getLDS(attributes = c("mgi_symbol"),
                      filters = "mgi_symbol",
                      values=gene_mouse,mart = mouse,
                      attributesL = c("hgnc_symbol"),
                      martL = human,uniqueRows = T)
    gene_human=gene_human[!duplicated(gene_human$MGI.symbol),]
    gene_human=gene_human[!duplicated(gene_human$HGNC.symbol),]
    Objects=Objects[gene_human$MGI.symbol,]
    row.names(Objects)=gene_human$HGNC.symbol

  }else{
    gene_human=row.names(Objects)
    print("Converting human gene symbols to mouse...")
    gene_mouse=getLDS(attributes = c("hgnc_symbol"),
                      filters = "hgnc_symbol",
                      values=gene_human,mart = human,
                      attributesL = c("mgi_symbol"),
                      martL = mouse,uniqueRows = T)
    gene_mouse=gene_mouse[!duplicated(gene_mouse$HGNC.symbol),]
    gene_mouse=gene_mouse[!duplicated(gene_mouse$MGI.symbol),]
    Objects=Objects[gene_mouse$HGNC.symbol,]
    row.names(Objects)=gene_mouse$MGI.symbol

  }

}else {
  if(species=='m'){
    gene_mouse= Objects
    print("Converting mouse gene symbols to human...")
    gene_human=getLDS(attributes = c("mgi_symbol"),
                      filters = "mgi_symbol",
                      values=gene_mouse,mart = mouse,
                      attributesL = c("hgnc_symbol"),
                      martL = human,uniqueRows = T)
    gene_human=gene_human[!duplicated(gene_human$MGI.symbol),]
    gene_human=gene_human[!duplicated(gene_human$HGNC.symbol),]
    Objects = gene_human$HGNC.symbol

  }else{
    gene_human=Objects
    print("Converting human gene symbols to mouse...")
    gene_mouse=getLDS(attributes = c("hgnc_symbol"),
                      filters = "hgnc_symbol",
                      values=gene_human,mart = human,
                      attributesL = c("mgi_symbol"),
                      martL = mouse,uniqueRows = T)
    gene_mouse=gene_mouse[!duplicated(gene_mouse$HGNC.symbol),]
    gene_mouse=gene_mouse[!duplicated(gene_mouse$MGI.symbol),]
    Objects = gene_mouse$MGI.symbol
  }
}

  return(Objects)
}




#' Get DE genes in ave_logFC
#'
#' @param markers the result of FindAllMarkers(), including the col named "cluster"
#' @param m the number of the genes in top
#'
#' @return a data.frame has the genes in top you choose
#' @export
#'
#' @examples
#' \dontrun{
#' input.file <- system.file('extdata','DEG.csv',package = 'LIANLABDATA')
#' markers <- read.csv(input.file,header = T, row.names = 1)
#' markers_top <- top_m( markers, m = 20)
#' }
top_m = function(markers,m){
  cluster <- avg_logFC <- p_val_adj <- avg_logFC <- NULL
  markers = subset(markers,subset = p_val_adj<0.05&avg_logFC>0.5)

  dim(markers)
  markers = markers[order(markers$cluster,markers$avg_logFC,decreasing = T),]

  mar2 = data.frame()
  for (i in levels(markers$cluster)) {
    mar = subset(markers,subset=cluster==i)
    if (nrow(mar) < m |nrow(mar) == m) {
      mar1 = as.data.frame(mar)
    }else{mar1 = as.data.frame(mar[1:m, ])}
    mar2 = rbind(mar2, mar1)
  }
  return(mar2)
}


#' Get DE genes in ave_logFC
#'
#' @param marker the result of DE genes of a cluster or the result between the two clusters
#' @param m the number of the genes in top and down
#'
#' @return a data.frame has the genes in top/down you choose
#' @export
#'
#' @examples
#' \dontrun{
#' input.file <- system.file('extdata','DEG.csv',package = 'LIANLABDATA')
#' markers <- read.csv(input.file,header = T, row.names = 1)
#' markers <- subset(markers,subset = cluster == 'Naive CD4 T')
#' markers_top <- top_m( markers, m = 20)
#' }
top_up_down_m = function(marker,m){
  p_val_adj <- avg_logFC <- NULL
  marker = subset(marker,subset=p_val_adj<0.05&abs(avg_logFC)>0.5)
  marker = marker[order(marker$avg_logFC,decreasing = T),]
  marker_up = marker[1:m,]
  marker_down = marker[(nrow(marker)-m):nrow(marker),]
  marker = rbind(marker_up,marker_down)
  marker$gene = rownames(marker)
  return(marker)
}


#' To see the distribution of the number of genes
#'
#' @param seurat.obj the seurat.obj  to cut
#' @importFrom ggplot2 geom_histogram
#' @export
#'
#' @examples
#' \dontrun{
#' input.file <- system.file('extdata','pbmc_small.RDS',package = 'LIANLABDATA')
#' load(input.file)
#' See_mRNA(pbmc_small)
#' }
See_mRNA = function(seurat.obj){
  nFeature_RNA <- NULL
  ggplot(seurat.obj@meta.data,aes(nFeature_RNA))+geom_histogram()
}


#' Find out if the gene is in the object
#'
#' @param Object seurat.obj/matrix (row)/data.frame (row)/character
#' @param genes the gene you interested
#'
#' @return the genes in the object
#' @export
#'
#' @examples
#' \dontrun{
#' input.file <- system.file('extdata','pbmc_1k.RDS',package = 'LIANLABDATA')
#' pbmc_1k <- readRDS(input.file)
#' genes = c('CD8A','CD3','CD4')
#' Get_Genes(pbmc_small,genes)
#'
#' input.file <- system.file('extdata','DEG.csv',package = 'LIANLABDATA')
#' markers <- read.csv(input.file,header = T, row.names = 1)
#' Get_Genes(markers,genes)
#'
#' all_genes = c('CD8A','CD3','CD8A1')
#' Get_Genes(all_genes,genes)
#' }
Get_Genes <- function(Object,genes){
  if (class(Object)=="Seurat") {
    inside = intersect(genes,rownames(Object@assays$RNA@counts))

    #求向量x与向量y中不同的元素(只取x中不同的元素)
    outside = setdiff(genes, inside)
    print("=======================================")
    if (length(inside)!=0) {
      for (i in inside) {
        print(paste(i,'is in the Object'))
      }
      print(paste('sum is',length(inside)))

      print("=======================================")
      print("=======================================")

      for (i in outside) {
        print(paste(i,'is not in the Object'))
      }
      print(paste('sum is',length(outside)))

    }else{print('There is no gene in the Object')}

    print("=======================================")

  }else{
    if (class(Object)=='data.frame'|class(Object)=='matrix') {
      inside = rownames(Object[which(rownames(Object)%in%genes),])

      #求向量x与向量y中不同的元素(只取x中不同的元素)
      outside = setdiff(genes, inside)
      print("=======================================")
      if (length(inside)!=0) {
        for (i in inside) {
          print(paste(i,'is in the object'))
        }
        print(paste('sum is',length(inside)))

        print("=======================================")
        print("=======================================")

        for (i in outside) {
          print(paste(i,'is not in the object'))
        }
        print(paste('sum is',length(outside)))

      }else{print('There is no gene in the object')}
      print("=======================================")
    }else{
      if (class(Object)=='character') {
        inside = Object[which(Object%in%genes)]

        #求向量x与向量y中不同的元素(只取x中不同的元素)
        outside = setdiff(genes, inside)
        print("=======================================")
        if (length(inside)!=0) {
          for (i in inside) {
            print(paste(i,'is in the object'))
          }
          print(paste('sum is',length(inside)))

          print("=======================================")
          print("=======================================")

          for (i in outside) {
            print(paste(i,'is not in the object'))
          }
          print(paste('sum is',length(outside)))

        }else{print('There is no gene in the object')}
        print("=======================================")
      }
    }
  }
  return(inside)
}


#' This function is used to replace the cluster idents in target
#'
#' @param target_object the target seurat.obj
#' @param replace_object the replace seurat.obj, including the new annotation
#'
#' @return the seurat.obj had new annotation
#' @export
#'
#' @examples
#' \dontrun{
#' target_object = myreplace_ident(target_object,replace_object)
#' }
#'
myreplace_ident=function(target_object,replace_object){
  celltype <- NULL
  target_object$celltype = as.character(target_object@active.ident)
  replace_object$celltype = as.character(replace_object@active.ident)

  for (i in levels(replace_object@active.ident)) {
    target_object$celltype[rownames(target_object@meta.data)%in%rownames(subset(replace_object@meta.data, celltype==i))] = i
  }

  target_object@active.ident = as.factor(target_object$celltype)

  return(target_object)
}


#' Citeseq divide HTOtag
#'
#' @param path pathways cellranger exported data
#' @param filename the name of the generated file
#' @param positive.quantile The quantile of inferred 'negative' distribution for each hashtag - over which the cell is considered 'positive'. Default is 0.99
#' @param width the width of the figure
#' @param height the height of the figure
#'
#' @importFrom Seurat HTODemux Read10X CreateAssayObject RidgePlot FeatureScatter
#' @return seurat.obj including
#' @export
#'
HTO_get_singlet=function(path,filename=NULL,positive.quantile=0.99,width=7,height=6){

  hash.ID <- NULL
  #import cellranger exported data and create seurat object
  data <- Read10X(data.dir = path)
  object <- CreateSeuratObject(counts = data$`Gene Expression`)

  #cell-hash classification
  object[["Antibody"]] <- CreateAssayObject(counts = data$`Antibody Capture`)
  object <- NormalizeData(object,assay="Antibody",normalization.method = "CLR")
  object <- HTODemux(object,assay = "Antibody",positive.quantile = positive.quantile)
  fre <- as.data.frame(table(object$hash.ID)/ncol(object))
  count <- as.data.frame(table(object$hash.ID))
  merge <- merge(fre,count,all=T,by="Var1")
  colnames(merge) <- c("Class","Frequency","Cell_count")
  write.csv(merge,paste0(filename,"_Antibody_classification.csv"),row.names = F)

  SeuratObject::Idents(object) <- "hash.ID"
  pdf(paste0(filename,"_HTO_Ridgeplot.pdf"), width = width, height = height)
  p <- RidgePlot(object, assay = "Antibody", features = rownames(object[["Antibody"]]), ncol = 1)
  print(p)
  dev.off()
  pdf(paste0(filename,"_HTO_scatter.pdf"), width = width, height = height)
  p <- FeatureScatter(object, feature1 = "Sampletag1", feature2 = "Sampletag2")
  print(p)
  dev.off()
  pdf(paste0(filename,"_HTO_violin.pdf"), width = width, height = height)
  p <- VlnPlot(object, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
  print(p)
  dev.off()

  #remove doublet and negative cells
  SeuratObject::Idents(object) <- "Antibody_classification.global"
  object <- subset(object, idents = "Singlet")

  return(object)
}



#' To test the package need
#'
#' @param package_list the character of the test package
#'
#' @export
#'
#' @examples
#' \dontrun{
#' test_package(c('Seurat','DOSE'))
#' }
test_package = function(package_list){
  i <- 1
  while(i < length(package_list) + 1){
    if (!requireNamespace(package_list[i], quietly = TRUE)) {
      print(paste("Please install the package",package_list[i]))
    }
    i = i + 1
  }
    if (i == length(package_list)) {
      stop('Please install the package.')
    }
}
