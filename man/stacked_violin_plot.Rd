% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/myviolin_plot.R
\name{stacked_violin_plot}
\alias{stacked_violin_plot}
\title{This function is used to plot the stacked violinplot of single cell RNAseq data}
\usage{
stacked_violin_plot(
  gene,
  seurat_object,
  cluster = NULL,
  limits.max = 7,
  width = 13,
  height = 10.3,
  flip = T,
  filename = "",
  text.size = 10,
  Mean = T,
  col = NULL
)
}
\arguments{
\item{gene}{a vector of the genes name to plot}

\item{seurat_object}{the seurat object of your data}

\item{cluster}{choose to cluster you interested}

\item{limits.max}{numeric,the max legend value to regulate the color levels}

\item{width}{the width of the figure}

\item{height}{the height of the figure}

\item{flip}{logical,set the direction of coordinate axis}

\item{filename}{the name of the generated file}

\item{text.size}{numeric,the text size of gene names}

\item{Mean}{whether choose to mean}

\item{col}{col of the figure}
}
\description{
This function is used to plot the stacked violinplot of single cell RNAseq data
}
\details{
This function is used to plot the stacked violinplot of single cell RNAseq data or any
seurat object for specified genes expressed on celltype clusters.
}
\examples{
\dontrun{
input.file <- system.file('extdata','pbmc_1k.RDS',package = 'LIANLABDATA')
pbmc_1k <- readRDS(input.file)
gene=c('CD8A','CD4')
}
}
