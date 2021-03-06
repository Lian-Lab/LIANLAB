% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Get_Percentplot.R
\name{Get_percentplot_split}
\alias{Get_percentplot_split}
\title{Get percentage of seurat.obj, which has more than one samples}
\usage{
Get_percentplot_split(
  seurat_object,
  filename = NULL,
  order = "increasing",
  split_by = NULL,
  colorm = NULL
)
}
\arguments{
\item{seurat_object}{seurat.obj, including cols could split}

\item{filename}{the name of the generated file}

\item{order}{the order of the samples in figure}

\item{split_by}{split by the seurat.obj by cols, default 'samples'}

\item{colorm}{the colors of the clusters}
}
\description{
Get percentage of seurat.obj, which has more than one samples
}
\examples{
\dontrun{
input.file <- system.file('extdata','pbmc_1k.RDS',package = 'LIANLABDATA')
pbmc_1k <- readRDS(input.file)
Get_percentplot_split(pbmc_1k)
}
}
