#' This function is used to perform GSVA analysis
#'
#' @details This function is used to perform GSVA analysis, including pheatmap part which exhibit
#'   the GSVA score of selected genesets in every samples or cells barplot part which exhibit the
#'   top20 different expressed genesets in two group
#' @param em expression matrix as genes in row, samples or cells in column
#' @param seurat_object it also support a seurat_object instead of expression matrix
#' @param gspath the correct path of your geneset file(.gmt)
#' @param species species of the seurat.obj or matrix
#' @param text.size the size of the text
#' @param topn number of the top pathways
#' @param min.sz minimum size of the resulting gene sets
#' @param filename the name of the generated file
#' @param clusters to choose the clusters you interested
#' @param pathname to choose the pathways you interested
#' @param width the width of the figure
#' @param height the heigh of the figure
#'
#' @importFrom GSEABase getGmt
#' @importFrom Seurat AverageExpression
#' @importFrom GSVA gsva
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @export
#'
#' @examples
#' \dontrun{
#' input.file <- system.file('extdata','pbmc_1k.RDS',package = 'LIANLABDATA')
#' pbmc_1k <- readRDS(input.file)
#' GSVA_pheatmap(pbmc_1k,species="mouse")
#' }
GSVA_pheatmap=function(em=NULL,seurat_object=NULL,
                       gspath=NULL,
                       species="human",text.size=2,topn=20,min.sz=2,filename="geneset",
                       clusters=NULL,pathname=NULL,
                       width = 6,height = 5){

  test_package(c("GSEABase",'GSVA'))

  RNA <- clusters <- NULL
  if(is.null(em)){
    if(is.null(seurat_object)){
      stop("No correct expression matrix supplied!")
    }else{
      if(!is.null(clusters)){
        seurat_object <- subset(seurat_object,idents=clusters)
      }
      em <- AverageExpression(seurat_object)$RNA
    }
  }else{
    if(!is.null(seurat_object)){
      warning("There are two data supplied,only expression matrix used!")
      em = em
    }else{
      em = em
    }
  }

  if(species=="mouse"){
    em <- Get_HM_gene_change(em,species = 'm')
  }
  if (is.null(gspath)) {
    gspath <- system.file('extdata','h.all.v7.1.symbols.gmt',package = 'LIANLAB')
  }
  geneset <- getGmt(gspath)
  escore <- gsva(as.matrix(em), min.sz=min.sz, max.sz=2000,geneset,verbose=FALSE, parallel.sz=8)
  write.csv(escore,paste0(filename,'_gsva_escore.csv'))

  if(!is.null(pathname)){
    escore = escore[pathname,]
  }
  escore_scale <- as.data.frame(t(scale(as.data.frame(t(escore)))))
  escore_scale = escore_scale[1:topn,]

  pdf(paste0("GSVA_heatmap_",filename,".pdf"),width = 1.15*width,height = height)
  p <- pheatmap(escore,fontsize_row = 4,cluster_rows = T,cluster_cols = T,
             treeheight_row = 0,treeheight_col = 0,
             color = rev(brewer.pal(n=11,name = "RdBu")))
  print(p)
  dev.off()
}


#' This function is used to perform GSVA analysis between two clusters
#'
#' @param em expression matrix as genes in row, samples or cells in column
#' @param seurat_object it also support a seurat_object instead of expression matrix
#' @param gspath the correct path of your geneset file(.gmt)
#' @param species "mouse"(default),'human'
#' @param group a vector of two group,if em is a expression data,it need supply
#' @param adjPvalueCutoff cutof the P.adjust
#' @param logFCcutoff cutoff the logFC
#' @param text.size the size of the text
#' @param topn number of the top pathways
#' @param min.sz minimum size of the resulting gene sets
#' @param width the width of the figure
#' @param height the height of the figure
#' @param filename the name of the generated file
#' @param ylimit limit of the y aes
#'
#' @importFrom utils head tail
#' @importFrom stats model.matrix
#' @importFrom ggplot2 geom_bar scale_fill_manual ylab coord_flip geom_text scale_y_continuous
#' @importFrom limma lmFit eBayes topTable decideTests
#' @return seurat.obj with the score of the pathways
#' @export
#'
#' @examples
#' \dontrun{
#' input.file <- system.file('extdata','pbmc_1k.RDS',package = 'LIANLABDATA')
#' pbmc_1k <- readRDS(input.file)
#' GSVA_pheatmap( pbmc_1k,species="mouse",group=c('1','2'))
#' }
GSVA_difference=function(em=NULL,seurat_object=NULL,
                         gspath=NULL,
                         species="mouse",group=NULL,adjPvalueCutoff=0.001,
                         logFCcutoff=0.1,text.size=2,topn=20,min.sz=2,
                         width = 6,
                         height = 5,filename=NULL,ylimit=10){
  test_package(c("GSEABase",'GSVA','limma'))
  pathname <- NULL
  if(is.null(em)){
    if(is.null(seurat_object)){
      stop("No correct expression matrix supplied!")
    }else{

      g1=subset(seurat_object,idents=group[1])
      g2=subset(seurat_object,idents=group[2])
      merge_em=cbind(as.data.frame(g1@assays$RNA@data),as.data.frame(g2@assays$RNA@data))
      group_merge=factor(c(rep(group[1],ncol(g1)),rep(group[2],ncol(g2))),
                         levels = group)
    }
  }else{
    if(!is.null(seurat_object)){
      warning("There are two data supplied,only expression matrix used!")
      merge_em=em
    }else{
      merge_em=em
    }
    group_merge=factor(c(rep("group1",group[1]),rep("group2",group[2])),levels = c("group1","group2"))
  }

  if(species=="mouse"){
    merge_em <- Get_HM_gene_change(merge_em,species = 'm')

  }
  if (is.null(gspath)) {
    gspath <- system.file('extdata','h.all.v7.1.symbols.gmt',package = 'LIANLAB')
  }
  geneset=getGmt(gspath)
  gs_es <- gsva(as.matrix(merge_em), geneset, min.sz=min.sz, max.sz=2000, verbose=FALSE, parallel.sz=12)
  adjPvalueCutoff <- adjPvalueCutoff
  logFCcutoff <- logFCcutoff
  design <- model.matrix(~ factor(group_merge))
  vsname=paste0(levels(group_merge)[2],"vs",levels(group_merge)[1])
  colnames(design) <- c(levels(group_merge)[1],vsname)
  fit <- lmFit(gs_es, design)
  fit <- eBayes(fit)
  allGeneSets <- topTable(fit, coef=vsname,number=Inf)
  DEgeneSets <- topTable(fit, coef=vsname,number=Inf,p.value = adjPvalueCutoff, lfc = logFCcutoff,adjust.method = "BH")
  res <- decideTests(fit, p.value=adjPvalueCutoff)

  if(nrow(DEgeneSets)==0){
    stop("No geneset through the filter!")
  }
  DEgeneSets=DEgeneSets[order(DEgeneSets$t,decreasing = F),]
  write.csv(DEgeneSets,paste0(vsname,"_DEgenesets.csv"))

  count_up_gs=length(which(DEgeneSets$t<0))
  count_down_gs=length(which(DEgeneSets$t>0))
  if(nrow(DEgeneSets)>=topn){
    if(count_up_gs>=round(topn/2)&count_down_gs>=round(topn/2)){
      DEgeneSets_top=rbind(head(DEgeneSets,round(topn/2)),tail(DEgeneSets,round(topn/2)))
    }else if(count_up_gs<round(topn/2)){
      DEgeneSets_top=rbind(head(DEgeneSets,topn-count_up_gs),tail(DEgeneSets,count_up_gs))
    }else{
      DEgeneSets_top=rbind(head(DEgeneSets,count_down_gs),tail(DEgeneSets,topn-count_down_gs))
    }
    DEgeneSets_top$pathname=factor(row.names(DEgeneSets_top),levels = row.names(DEgeneSets_top))
    group_top=factor(c(rep("down",length(which(DEgeneSets_top$t<0))),
                       rep("up",length(which(DEgeneSets_top$t>0)))),levels = c("down","up"))
    DEgeneSets_top$group=group_top
  }else{
    DEgeneSets_top=DEgeneSets
    DEgeneSets_top$pathname=factor(row.names(DEgeneSets_top),levels = row.names(DEgeneSets_top))
    group_top=factor(c(rep("down",length(which(DEgeneSets_top$t<0))),
                       rep("up",length(which(DEgeneSets_top$t>0)))),levels = c("down","up"))
    DEgeneSets_top$group=group_top
  }
  count_up=length(which(DEgeneSets_top$group=="up"))
  count_down=length(which(DEgeneSets_top$group=="down"))
  ###different geneset plot-----
  if(count_up>0&count_down>0){
    pdf(paste0(filename,"GSVA_",vsname,"_DEgeneset.pdf"),width = width,height = height)
    g1 <- ggplot(DEgeneSets_top,aes(pathname,t,fill=group))+geom_bar(stat = "identity")+
      scale_fill_manual(values = c("lightskyblue","hotpink2"))+
      theme(panel.border = element_blank(),panel.grid.major = element_blank(),
            panel.background = element_blank(),
            panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
            axis.text.x = element_text(vjust = 0.5, hjust = 0.5))+
      theme(axis.title.y =  element_blank())+
      theme(axis.ticks.y =  element_blank())+
      theme(axis.line.y =   element_blank())+
      theme(axis.text.y =   element_blank())+
      ylab("t value of GSVA score")+
      coord_flip()+
      geom_text(aes(y=c(rep(0.2,nrow(DEgeneSets_top))),label=c(as.character(pathname[1:count_down]),
                     rep("",count_up)),
                    vjust = 0.5,hjust=0),size=text.size) +
      geom_text(aes(y=c(rep(-0.2,nrow(DEgeneSets_top))),label=c(rep("",count_down),
                    as.character(pathname[(count_down+1):nrow(DEgeneSets_top)])),
                    vjust = 0.5,hjust=1),size=text.size)
    print(g1)
    dev.off()
  }else if(count_up>0&count_down==0){
    warning("There are no down regulated genesets found under these filter conditions,
          you could try smaller logFCcutoff value...")
    pdf(paste0(filename,"GSVA_",vsname,"_DEgeneset.pdf"),width = width,height = height)
    g1 <- ggplot(DEgeneSets_top,aes(pathname,t,fill=group))+geom_bar(stat = "identity")+
      scale_fill_manual(values = c("hotpink2"))+
      theme(panel.border = element_blank(),panel.grid.major = element_blank(),
            panel.background = element_blank(),
            panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
            axis.text.x = element_text(vjust = 0.5, hjust = 0.5))+
      theme(axis.title.y = element_blank())+
      theme(axis.ticks.y = element_blank())+
      theme(axis.line.y = element_blank())+
      theme(axis.text.y = element_blank())+
      ylab("t value of GSVA score")+
      coord_flip()+
      geom_text(aes(y=c(rep(-0.2,nrow(DEgeneSets_top))),label=pathname,
                    vjust = 0.5,hjust=1),size=2)+
      scale_y_continuous(limits = c(-ylimit,10*ceiling(max(DEgeneSets_top$t)/10)))
    print(g1)
    dev.off()
  }else if(count_up==0&count_down>0){
    warning("There are no up regulated genesets found under these filter conditions,
          you could try smaller logFCcutoff value...")
    pdf(paste0(filename,"GSVA_",vsname,"_DEgeneset.pdf"),width = width,height = height)
    g1 <- ggplot(DEgeneSets_top,aes(pathname,t,fill=group))+geom_bar(stat = "identity")+
      scale_fill_manual(values = c("lightskyblue"))+
      theme(panel.border = element_blank(),panel.grid.major = element_blank(),
            panel.background = element_blank(),
            panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
            axis.text.x = element_text(vjust = 0.5, hjust = 0.5))+
      theme(axis.title.y =  element_blank())+
      theme(axis.ticks.y =  element_blank())+
      theme(axis.line.y =   element_blank())+
      theme(axis.text.y =   element_blank())+
      ylab("t value of GSVA score")+
      coord_flip()+
      geom_text(aes(y=c(rep(0.2,nrow(DEgeneSets_top))),label=pathname,
                    vjust = 0.5,hjust=0),size=2)+
      scale_y_continuous(limits = c(10*floor(min(DEgeneSets_top$t)/10),ylimit))
    print(g1)
    dev.off()
  }else{
    stop("There are no diffrent genesets found under these filter conditions,
          please try smaller logFCcutoff value...")
  }
  if(!is.null(seurat_object)){
    seurat_object@neighbors=c(list(diff_gsname=as.character(DEgeneSets_top$pathname)),
                              seurat_object@neighbors)
    return(seurat_object)
  }
}




#' Insert pathway score into the seurat.obj
#'
#' @param seurat_object it also support a seurat_object instead of expression matrix
#' @param gspath the correct path of your geneset file(.gmt)
#' @param species "mouse"(default),'human'
#' @param min.sz minimum size of the resulting gene sets
#' @param parallel.sz number of threads of execution to use when doing the calculations in parallel
#'
#' @importFrom Seurat FindVariableFeatures as.sparse
#' @return seurat.obj with the score of the pathways
#' @export
#'
#' @examples
#' \dontrun{
#' input.file <- system.file('extdata','pbmc_1k.RDS',package = 'LIANLABDATA')
#' pbmc_1k <- readRDS(input.file)
#' GSVA_pheatmap(pbmc_1k,species="mouse")
#' }
GSVA_score_insert=function(seurat_object=NULL,
                           gspath=NULL,
                           species="mouse",min.sz = 2,parallel.sz = 12){
  test_package(c("GSEABase",'GSVA'))
  if(is.null(seurat_object)){
    stop("No seurat object supplied!")
  }else{
    seurat_object2=FindVariableFeatures(seurat_object)
    seurat_object2=subset(seurat_object2,
                          features=seurat_object2@assays$RNA@var.features)
    total_matrix=as.data.frame(seurat_object2@assays$RNA@data)
  }

  if(species=="mouse"){
    total_matrix = Get_HM_gene_change(total_matrix, species = 'm')

  }
  if (is.null(gspath)) {
    gspath <- system.file('extdata','h.all.v7.1.symbols.gmt',package = 'LIANLAB')
  }
  geneset=getGmt(gspath)
  total_es <- gsva(as.matrix(total_matrix), geneset,
                   min.sz = min.sz, max.sz = 2000, verbose = FALSE, parallel.sz = parallel.sz)
  total_es=as.data.frame(t(scale(as.data.frame(t(total_es)))))

  new_data_slot=rbind(as.data.frame(seurat_object@assays$RNA@data),total_es)
  seurat_object@assays$RNA@data=as.sparse(new_data_slot)

  new_counts_slot=rbind(as.data.frame(seurat_object@assays$RNA@counts),total_es)
  seurat_object@assays$RNA@counts=as.sparse(new_counts_slot)

  new_scale_slot=rbind(as.data.frame(seurat_object@assays$RNA@scale.data),total_es)
  seurat_object@assays$RNA@scale.data=as.matrix(new_scale_slot)

  seurat_object@neighbors=c(list(pathname=row.names(total_es),GSVA_matrix=total_es),
                            seurat_object@neighbors)
  write.csv(total_es,"Inserted_gsva_data.csv")
  write.csv(data.frame(pathname=row.names(total_es)),"Inserted_geneset_name.csv")
  return(seurat_object)
}



#' To find geneset markers of each cluster
#'
#' @param em expression matrix as genes in row, samples or cells in column
#' @param seurat_object it also support a seurat_object instead of expression matrix
#' @param gspath the correct path of your geneset file(.gmt)
#' @param species species of the seurat.obj or matrix
#' @param text.size the size of the text
#' @param topn number of the top pathways
#' @param min.sz minimum size of the resulting gene sets
#' @param filename the name of the generated file
#' @param adjPvalueCutoff cutof the P.adjust
#' @param logFCcutoff cutof the logFC
#' @param parallel.sz number of threads of execution to use when doing the calculations in parallel.
#'
#' @return seurat.obj with the marker_gs
#' @export
#'
#' @examples
#' \dontrun{
#' input.file <- system.file('extdata','pbmc_1k.RDS',package = 'LIANLABDATA')
#' pbmc_1k <- readRDS(input.file)
#' GSVA_pheatmap( seurat_object = pbmc_1k,species="mouse")
#' }
GSVA_find_gs_markers=function(em=NULL,seurat_object=NULL,
                              gspath=NULL,
                              species="mouse",adjPvalueCutoff=0.001,filename=NULL,
                              logFCcutoff=0.1,text.size=2,topn=20,min.sz=2,
                              parallel.sz=12){
  test_package(c("GSEABase",'GSVA'))
  group1 <- group2 <- NULL
  marker_gs <- data.frame()
  for (i in levels(seurat_object)) {
    if(is.null(em)){
      if(is.null(seurat_object)){
        stop("No correct expression matrix supplied!")
      }else{
        print(paste0("####-------",i," cluster------#####"))
        group1=i
        g1=subset(seurat_object,idents=group1)
        g2=subset(seurat_object,idents=setdiff(levels(seurat_object),group1))
        merge_em=cbind(as.data.frame(g1@assays$RNA@data),as.data.frame(g2@assays$RNA@data))
        group_merge=factor(c(rep("group1",ncol(g1)),rep("group2",ncol(g2))),
                           levels = c("group1","group2"))
      }
    }else{
      if(!is.null(seurat_object)){
        warning("There are two data supplied,only expression matrix used!")
        merge_em <- em
      }else{
        merge_em <- em
      }
      group_merge=factor(c(rep("group1",group1),rep("group2",group2)),levels = c("group1","group2"))
    }

    if(species=="mouse"){
      merge_em <- Get_HM_gene_change(merge_em,species = 'm')
    }
    if (is.null(gspath)) {
      gspath <- system.file('extdata','h.all.v7.1.symbols.gmt',package = 'LIANLAB')
    }
    geneset <- getGmt(gspath)
    gs_es <- gsva(as.matrix(merge_em), geneset, min.sz=min.sz,
                  max.sz=2000, verbose=FALSE, parallel.sz=parallel.sz)
    adjPvalueCutoff <- adjPvalueCutoff
    logFCcutoff <- logFCcutoff
    design <- model.matrix(~ factor(group_merge))
    vsname=paste0("group2","vs","group1")
    colnames(design) <- c("group1",vsname)
    fit <- lmFit(gs_es, design)
    fit <- eBayes(fit)
    allGeneSets <- topTable(fit, coef=vsname,number=Inf)
    DEgeneSets <- topTable(fit, coef=vsname,number=Inf,p.value=adjPvalueCutoff, lfc = logFCcutoff,adjust.method = "BH")
    res <- decideTests(fit, p.value=adjPvalueCutoff)
    #head(DEgeneSets)
    if(nrow(DEgeneSets)==0){
      stop("No geneset through the filter!")
    }
    DEgeneSets=DEgeneSets[order(DEgeneSets$t,decreasing = F),]
    count_up_gs=length(which(DEgeneSets$t>0))

    if(count_up_gs>=topn){
      DEgeneSets_top_name=head(row.names(DEgeneSets),topn)
    }else{
      DEgeneSets_top_name=head(row.names(DEgeneSets),count_up_gs)
    }
    marker_gs_i=data.frame(cluster=c(rep(i,length(DEgeneSets_top_name))),geneset=DEgeneSets_top_name)
    marker_gs=rbind(marker_gs,marker_gs_i)
  }
  write.csv(marker_gs,paste0("marker_geneset_",filename,".csv"))
  if(!is.null(seurat_object)){
    seurat_object@neighbors=c(list(marker_gs=marker_gs),seurat_object@neighbors)
    return(seurat_object)
  }
}

