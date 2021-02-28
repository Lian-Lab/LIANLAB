#' GO and KEGG analyse for the DEG
#'
#' @param species species of genes
#' @param up_gene character of up genes
#' @param down_gene character of down genes
#' @param png.width width of the figure in png
#' @param png.height height of the figure in png
#' @param pdf.width width of the figure in pdf
#' @param pdf.height height of the figure in pdf
#' @param show number of the pathway to show
#' @param filename the name of the generated file
#'
#' @importFrom utils write.csv
#' @importFrom clusterProfiler bitr enrichGO simplify enrichKEGG
#' @importFrom enrichplot dotplot cnetplot
#' @export
#'
#' @examples
#' \dontrun{
#' input.file <- system.file('extdata','DEG.csv',package = 'LIANLAB')
#' DEG_matrix <- read.csv(input.file, header = TRUE, row.names = 1)
#' up_gene <- rownames(subset(DEG_matrix,avg_logFC > 0))
#' down_gene <- rownames(subset(DEG_matrix,avg_logFC < 0))
#' DEG <- myGO_KEGG(species="h",up_gene,down_gene)
#' }
myGO_KEGG=function(species="m",up_gene,down_gene,png.width=1000,png.height=500,pdf.width=700,pdf.height=5,show=20,filename=NULL){

  test_package(c("org.Hs.eg.db",'org.Mm.eg.db','graphics','utils','DOSE'))

  if(species=="h"|species=="human"){
    organism_db="org.Hs.eg.db"
    organism="hsa"
  }else{
    organism_db="org.Mm.eg.db"
    organism="mmu"
  }

  #upgene
  up_geneid <- bitr(up_gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=organism_db)
  up_geneid <- up_geneid$ENTREZID
  go_pvalue <- 0.05
  up_ego <- enrichGO(gene=up_geneid,OrgDb = organism_db,
                     keyType = 'ENTREZID', ont="BP",pvalueCutoff=go_pvalue,readable=TRUE)
  n <- 0
  while(length(up_ego)==0 & n < 10){
    go_pvalue <- go_pvalue+0.05
    up_ego <- enrichGO(gene=up_geneid,OrgDb = organism_db,
                       keyType = 'ENTREZID', ont="BP",pvalueCutoff=go_pvalue,readable=TRUE)
    n=n+1
  }

  n <- 0
  if(length(up_ego)==0){
    print("There were no UP GO terms enriched,please check your data!!")
  }else{
    while (nrow(up_ego)<5&n<10){
      go_pvalue <- go_pvalue+0.05
      up_ego <- enrichGO(gene=up_geneid,OrgDb = organism_db,
                         keyType = 'ENTREZID', ont="BP",pvalueCutoff=go_pvalue,readable=TRUE)
      n=n+1
    }
  }

  #downgene
  down_geneid <- bitr(down_gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=organism_db)
  down_geneid <- down_geneid$ENTREZID
  go_pvalue <- 0.05
  down_ego <- enrichGO(gene = down_geneid,OrgDb = organism_db,
                       keyType = 'ENTREZID', ont="BP",pvalueCutoff=go_pvalue,readable=TRUE)

  n <- 0
  while(length(down_ego)==0&n<10){                                                                                                                                                                                                                 go_pvalue=go_pvalue+0.05
  down_ego <- enrichGO(gene=down_geneid,OrgDb = organism_db,
                       keyType = 'ENTREZID', ont="BP",pvalueCutoff=go_pvalue,readable=TRUE)
  n=n+1
  }
  n <- 0
  if(length(down_ego)==0){
    print("There were no DOWN GO terms enriched,please check your data!!")
  }else{
    while (nrow(down_ego)<5&n<10){
      go_pvalue=go_pvalue+0.05
      down_ego <- enrichGO(gene=down_geneid,OrgDb = organism_db,
                           keyType = 'ENTREZID', ont="BP",pvalueCutoff=go_pvalue,readable=TRUE)
      n=n+1
    }
  }

  up_ego <- simplify(up_ego)
  down_ego <- simplify(down_ego)
  utils::write.csv(up_ego,paste0(filename,"_GO_UP_result.csv"))
  utils::write.csv(down_ego,paste0(filename,"GO_DOWN_result.csv"))

  if (nrow(up_ego)>1) {
    png(paste0(filename,"_GO_UP_top20_barplot.png"), width = png.width, height = png.height)
    b1=graphics::barplot(up_ego, showCategory=show,title=paste0("GO_","_UP")) #条形图
    print(b1)
    dev.off()
    png(paste0(filename,"_GO_UP_top20_dotplot.png"), width = png.width, height = png.height)
    d1=dotplot(up_ego,showCategory=show,title=paste0("GO_","_UP"))
    print(d1)
    dev.off()
    pdf(paste0(filename,"_GO_UP_top20_barplot.pdf"), width = pdf.width, height = pdf.height)
    b2=graphics::barplot(up_ego, showCategory=show,title=paste0("GO_","_UP")) #条形图
    print(b2)
    dev.off()
    pdf(paste0(filename,"_GO_UP_top20_dotplot.pdf"), width = pdf.width, height = pdf.height)
    d2=dotplot(up_ego,showCategory=show,title=paste0("GO_","_UP"))
    print(d2)
    dev.off()
  }

  if (nrow(down_ego)>1) {
    pdf(paste0(filename,"_GO_down_top20_barplot.pdf"), width = pdf.width, height = pdf.height)
    b3=graphics::barplot(down_ego, showCategory=show,title=paste0("GO_","_DOWN")) #条形图
    print(b3)
    dev.off()
    pdf(paste0(filename,"_GO_down_top20_dotplot.pdf"), width = pdf.width, height = pdf.height)
    d3=dotplot(down_ego,showCategory=show,title=paste0("GO_","_DOWN"))
    print(d3)
    dev.off()
    png(paste0(filename,"_GO_down_top20_barplot.png"), width = png.width, height = png.height)
    b4=graphics::barplot(down_ego, showCategory=show,title=paste0("GO_","_DOWN")) #条形图
    print(b4)
    dev.off()
    png(paste0(filename,"_GO_down_top20_dotplot.png"), width = png.width, height = png.height)
    d4=dotplot(down_ego,showCategory=show,title=paste0("GO_","_DOWN"))
    print(d4)
    dev.off()
  }



  kegg_pvalue <- 0.05
  up_ekk <- enrichKEGG(gene=up_geneid,keyType = 'ncbi-geneid',organism = organism, pAdjustMethod = "BH",
                       pvalueCutoff=kegg_pvalue)
  n <- 0
  while (nrow(up_ekk) < 5 & n < 10){
    kegg_pvalue <- kegg_pvalue+0.05
    up_ekk <- enrichKEGG(gene=up_geneid,keyType = 'ncbi-geneid',organism = organism, pAdjustMethod = "BH",
                         pvalueCutoff=kegg_pvalue)
    n = n+1
  }
  up_symbol_ekk <- DOSE::setReadable(up_ekk, OrgDb = organism_db,keyType = "ENTREZID")

  kegg_pvalue <- 0.05
  down_ekk <- enrichKEGG(gene=down_geneid,keyType = 'ncbi-geneid',organism = organism, pAdjustMethod = "BH",
                         pvalueCutoff = kegg_pvalue)
  n <- 0
  while (nrow(down_ekk)<5 & n<10){
    kegg_pvalue <- kegg_pvalue+0.05
    down_ekk <- enrichKEGG(gene=down_geneid,keyType = 'ncbi-geneid',organism = organism, pAdjustMethod = "BH",
                           pvalueCutoff=kegg_pvalue)
    n=n+1
  }
  down_symbol_ekk <- DOSE::setReadable(down_ekk, OrgDb = organism_db,keyType = "ENTREZID")
  utils::write.csv(up_symbol_ekk,paste0(filename,"_KEGG_UP_result.csv"))
  utils::write.csv(down_symbol_ekk,paste0(filename,"_KEGG_DOWN_result.csv"))

  if (nrow(up_symbol_ekk)>1) {
    png(paste0(filename,"_KEGG_up_cnetplot.png"), width = png.width, height = png.height)
    k1=cnetplot(up_symbol_ekk, categorySize="pvalue")
    print(k1)
    dev.off()
    pdf(paste0(filename,"_KEGG_up_cnetplot.pdf"), width = pdf.width, height = pdf.height)
    k2=cnetplot(up_symbol_ekk, categorySize="pvalue",showCategory=10)
    print(k2)
    dev.off()
  }

  if (nrow(down_symbol_ekk)>1) {
    png(paste0(filename,"_KEGG_down_cnetplot.png"), width = png.width, height = png.height)
    k3=cnetplot(down_symbol_ekk, categorySize="pvalue",showCategory=10)
    print(k3)
    dev.off()
    pdf(paste0(filename,"_KEGG_down_cnetplot.pdf"), width = pdf.width, height = pdf.height)
    k4=cnetplot(down_symbol_ekk, categorySize="pvalue",showCategory=10)
    print(k4)
    dev.off()
  }

  if (nrow(up_ekk)>1) {
    png(paste0(filename,"_KEGG_UP_dotplot.png"), width = png.width, height = png.height)
    k5=dotplot(up_ekk,showCategory = show,title = paste0("KEGG_","_UP" ))
    print(k5)
    dev.off()
    pdf(paste0(filename,"_KEGG_UP_dotplot.pdf"), width = pdf.width, height = pdf.height)
    k6=dotplot(up_ekk,showCategory = show,title = paste0("KEGG_","_UP" ))
    print(k6)
    dev.off()
    png(paste0(filename,"_KEGG_UP_barplot.png"), width = png.width, height = png.height)
    k7=graphics::barplot(up_ekk,showCategory = show,title = paste0("KEGG_","_UP" ))
    print(k7)
    dev.off()
    pdf(paste0(filename,"_KEGG_UP_barplot.pdf"), width = pdf.width, height = pdf.height)
    k8=graphics::barplot(up_ekk,showCategory =show,title = paste0("KEGG_","_UP" ))
    print(k8)
    dev.off()
  }

  if (nrow(down_ekk)>1) {
    png(paste0(filename,"_KEGG_DOWN_dotplot.png"), width = png.width, height = png.height)
    k9=dotplot(down_ekk,showCategory = show,title = paste0("KEGG_","_DOWN" ))
    print(k9)
    dev.off()
    pdf(paste0(filename,"_KEGG_DOWN_dotplot.pdf"), width = pdf.width, height = pdf.height)
    k10=dotplot(down_ekk,showCategory = show,title = paste0("KEGG_","_DOWN" ))
    print(k10)
    dev.off()
    png(paste0(filename,"_KEGG_DOWN_barplot.png"), width = png.width, height = png.height)
    k11=graphics::barplot(down_ekk,showCategory = show,title = paste0("KEGG_","_DOWN" ))
    print(k11)
    dev.off()
    pdf(paste0(filename,"_KEGG_DOWN_barplot.pdf"), width = pdf.width, height = pdf.height)
    k12=graphics::barplot(down_ekk,showCategory = show,title = paste0("KEGG_","_DOWN" ))
    print(k12)
    dev.off()
  }




}
