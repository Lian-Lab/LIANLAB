% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Get_GSEA.R
\name{Get_GSEA}
\alias{Get_GSEA}
\title{Title}
\usage{
Get_GSEA(
  object,
  species = c("m", "h"),
  category = NULL,
  filename = NULL,
  p.valueCutoff = 0.05
)
}
\arguments{
\item{object}{DEG data.frame or matrix, including col named 'avg_logFC'}

\item{species}{species of DEG}

\item{category}{geneset in GSEA, including 'c2'(default),'c3','c4','c5','c6','c7','hallmark'}

\item{filename}{the name of the generated file}

\item{p.valueCutoff}{cut off the p value}
}
\value{
result of the gsea
}
\description{
Title
}
\examples{
\dontrun{
input.file <- system.file('extdata','DEG.csv',package = 'LIANLABDATA')
markers <- read.csv(input.file,header = T, row.names = 1)
result <- Get_GSEA( markers,species = 'm',category='c2')
}
}
