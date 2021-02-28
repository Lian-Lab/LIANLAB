#' To analyse the seurat.obj easily
#'
#' @param seurat_object seurat.obj
#' @param gene_ann gene information rownames are genenames, and a col named 'gene_short_name'
#' @param n_cell barcode of the mixed cells, default NULL
#' @param num_expressed the lowest gene expression, default
#' @param batch_if batch factor, default NULL
#' @param remove_batch whether to remove the batch effect caused by the batch factor
#' @param dispersion_empiricals dispersion empiricals
#' @param mean_expressions mean expressions
#' @param ordering_genes marks genes that will be used for clustering in subsequent calls to clusterCells
#'
#' @importFrom monocle newCellDataSet detectGenes differentialGeneTest setOrderingFilter plot_ordering_genes reduceDimension plot_cell_trajectory dispersionTable orderCells
#' @importFrom VGAM negbinomial.size
#' @importFrom BiocGenerics estimateSizeFactors estimateDispersions
#' @return a list, including cds and cluster DE genes
#' @export
#'
#' @examples
#' \dontrun{
#' input.file <- system.file('extdata','pbmc_1k.RDS',package = 'LIANLABDATA')
#' pbmc_1k <- readRDS(input.file)
#' cds.list <- Get_monocle(pbmc_1k)
#' }
Get_monocle = function(seurat_object,gene_ann=NULL,n_cell=NULL,
                       num_expressed=1,batch_if=NULL,remove_batch = F,
                       dispersion_empiricals=0.1,mean_expressions=0.01,ordering_genes=NULL){
  test_package(c("monocle",'VGAM','BiocGenerics'))
  num_cells_expressed <- mean_expression <- dispersion_empirical <- NULL
  ###Monocle2 三件套
  #表达矩阵
  ct <- seurat_object@assays$RNA@data
  sample_ann <- seurat_object@meta.data
  if (is.null(gene_ann)) {
    ###基因信息
    gene_ann <- data.frame(
      "gene_short_name" = row.names(ct),
      row.names = row.names(ct)
    )
  }

  ###转换成AnnotationDataframe对象
  pd <- new("AnnotatedDataFrame",
            data=sample_ann)
  fd <- new("AnnotatedDataFrame",
            data=gene_ann)

  # 构建CDS对象
  #expressionFamily ，选择数据分布:
  #FPKM/TPM 值是log-正态分布的；UMIs和原始count值用负二项分布模拟的效果更好
  #负二项分布有两种方法，这里选用了negbinomial.size，另外一种negbinomial稍微更准确一点，但速度大打折扣，它主要针对非常小的数据集
  sc_cds <- newCellDataSet(
    as.matrix(ct),
    phenoData = pd,
    featureData =fd,
    expressionFamily = negbinomial.size(),
    lowerDetectionLimit=1)

  #####质控过滤
  cds <- sc_cds

  ##剔除杂细胞
  if (!is.null(n_cell)) {
    cds <- cds[,-n_cell]
  }


  # 设置一个基因表达量的过滤阈值，结果会在cds@featureData@data中新增一列num_cells_expressed，记录这个基因在多少细胞中有表达
  cds <- detectGenes(cds, min_expr = 0.1)
  # 结果保存在cds@featureData@data
  print(head(cds@featureData@data))
  ##基因过滤
  #选出表达大于5的
  expressed_genes <- row.names(subset(cds@featureData@data,
                                      num_cells_expressed >= num_expressed))
  cds <- cds[expressed_genes,]
  #####计算Size factors和"dispersion" values，
  # Size factors可以帮助我们对不同细胞中mRNA表达的差异进行归一化处理
  # "dispersion" values将有助于后续的基因差异表达分析。
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)

  disp_table <- dispersionTable(cds) # 挑有差异的
  unsup_clustering_genes <- subset(disp_table, mean_expression >= mean_expressions) # 挑表达量不太低的
  unsup_clustering_genes <- subset(unsup_clustering_genes, dispersion_empirical >= dispersion_empiricals) # 挑表达量不太低的

  cluster_DEG_genes <- differentialGeneTest(cds[as.character(unsup_clustering_genes$gene_id),],fullModelFormulaStr = '~seurat_clusters',cores = 6)
  dim(cluster_DEG_genes)
  cluster_DEG_gene <- rownames(cluster_DEG_genes)[order(cluster_DEG_genes$qval)][1:2000]

  cds <- setOrderingFilter(cds, as.character(cluster_DEG_gene))  # 准备聚类基因名单
  print(length(cluster_DEG_genes))
  print(plot_ordering_genes(cds))

  if (!is.null(ordering_genes)) {
    cds <- setOrderingFilter(cds, ordering_genes)
  }
  ###降维
  # 默认使用DDRTree的方法
  if (remove_batch==F) {
    # 进行降维
    cds <- reduceDimension(cds, max_components = 2,
                           reduction_method = 'DDRTree')
  }else{
    ##去除批次的影响  修改residualModelFormulaStr
    i = 2
    batch_ifs = batch_ifs = paste0("~",batch_if[1])
    while(i < length(batch_if)+1){
      batch_ifs = paste(batch_ifs,"+",batch_if[i])
      i = i + 1
    }
    cds <- reduceDimension(cds, max_components = 2,
                           reduction_method = 'DDRTree',
                           residualModelFormulaStr = as.character(batch_ifs) )
  }
  ###细胞排序
  cds <- orderCells(cds)
  ##可视化
  p <- plot_cell_trajectory(cds, color_by = "Pseudotime")
  print(p)

  mono.cds = list('monocle_object'=cds,'cluster_DEG_genes'=cluster_DEG_genes)
  return(mono.cds)

}
