% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Get_Seurat.R
\name{Get_Seurat}
\alias{Get_Seurat}
\title{To get seurat.obj easily}
\usage{
Get_Seurat(
  object,
  minfeature = 200,
  mincell = 10,
  species = c("h", "m"),
  minFeatureRNA = NULL,
  maxFeatureRNA = NULL,
  Percent.mt = NULL,
  is.tumor = F,
  projectname = NULL,
  dim.FindNeighbors = 1:30,
  dims.TSNE = 1:30
)
}
\arguments{
\item{object}{gene expression matrix}

\item{minfeature}{Include genes where at least this many features are detected.}

\item{mincell}{Include cells where at least this many features are detected.}

\item{species}{the species of the matrix}

\item{minFeatureRNA}{the min cut off the genes count}

\item{maxFeatureRNA}{the max cut off the genes count}

\item{Percent.mt}{the max cut off of the percetnage of the mitochondrial gene}

\item{is.tumor}{whether the sample is tumor}

\item{projectname}{the name of the project}

\item{dim.FindNeighbors}{number of dims to find neighbors}

\item{dims.TSNE}{number of dims to t-SNE}
}
\value{
the seurat.obj including PCA, t-SNE
}
\description{
To get seurat.obj easily
}
\examples{
\dontrun{
input.file <- system.file('extdata','pbmc_small.csv',package = 'LIANLAB')
pbmc_small <- read.csv(input.file, header = T, row.names = 1)

seurat.obj = CreateSeuratObject(counts = pbmc_small, min.features = 200,min.cells = 10)
dim(Abseq.seurat)
See_mRNA(Abseq.seurat)
data[["percent.mt"]] <- PercentageFeatureSet(Abseq.seurat, pattern = "^MT.")  # MT(human)/mt(mouse)
VlnPlot(Abseq.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
table(Abseq.seurat[["percent.mt"]]<10)

seurat.obj = Get_Seurat(pbmc_small, minfeature = 200, mincell=10,species='h')
}

}
