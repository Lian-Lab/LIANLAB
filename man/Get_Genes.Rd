% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Get_Mixture.R
\name{Get_Genes}
\alias{Get_Genes}
\title{Find out if the gene is in the object}
\usage{
Get_Genes(Object, genes)
}
\arguments{
\item{Object}{seurat.obj/matrix (row)/data.frame (row)/character}

\item{genes}{the gene you interested}
}
\value{
the genes in the object
}
\description{
Find out if the gene is in the object
}
\examples{
\dontrun{
input.file <- system.file('extdata','pbmc_1k.RDS',package = 'LIANLABDATA')
pbmc_1k <- readRDS(input.file)
genes = c('CD8A','CD3','CD4')
Get_Genes(pbmc_small,genes)

input.file <- system.file('extdata','DEG.csv',package = 'LIANLABDATA')
markers <- read.csv(input.file,header = T, row.names = 1)
Get_Genes(markers,genes)

all_genes = c('CD8A','CD3','CD8A1')
Get_Genes(all_genes,genes)
}
}
