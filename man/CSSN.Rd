\name{CSSN}
\alias{CSSNEst}
\alias{CSSNPredict}
\title{Recover spatially-varying cell-specific gene co-expression networks for single-cell spatial expression data and predict cell-specific gene co-expression networks on tissue positions where cells are not captured.}
\usage{
CSSNEst(X, cell.info, nu, d = 0.1, m.info = 70, is.scale = TRUE)
CSSNPredict(GN, cell.info, miss.indx, m.info = 70)
}
\arguments{
\item{X}{spatial single-cell RNA-sequencing data, where rows represent genes and columns represent cells. Matrix values need to be non-negative and continuous.}

\item{cell.info}{cell information matrix (n,3), where rows stands for cells, the first column represents cell type label of cells, the second and third columns are Centroid_X Centroid_Y coordinates respectively.}

\item{nu}{n dimension vector, the degree of freedom of the prior Inverse-Wishart distribution of \eqn{\Sigma_i}.}

\item{d}{numeric threshold parameter, the range is 0~1. The default value is 0.1.}

\item{m.info}{numeric, cell density parameter. The default is 70.}

\item{is.scale}{bool, if TRUE, the variance of gene expression data in each cell type will be scaled to 1. The default is TRUE.}

\item{GN}{estimated gene co-expression networks, an array of (G,G,n).}

\item{miss.indx}{coordinates of missing cells, a (miss.num, 2) matrix, where miss.num is the number of missing cells.}

}
\value{
CSSNEst returns an array with (G,G,miss.num) of gene co-expression networks for each missing cell and CSSNEst returns an R list including the following information.
\item{Correlation Matrix}{Gene correlation matrix of n cells, an array with dimension (G,G,n), where G is the number of genes and n is the number of cells.}
\item{Gene Networks}{Gene co-expression networks of n cells, an array with dimension (G,G,n), where G is the number of genes and n is the number of cells.}
}
\description{
The function CSSNEst implements a two-step algorithm to estimate cell-specific gene co-expression networks, the function CSSNPredict implements a two-step algorithm to estimate cell-specific gene co-expression networks.
}
\examples{
library(CSSN)
library(ggplot2)
library(pheatmap)
data(example_data)

# gene number
G <- nrow(X)
# cell number
n <- ncol(X)

#---- set spatial pattern manually----
pal <- c(rgb(221, 160, 221, maxColorValue = 255), 
        rgb(0, 206, 209, maxColorValue = 255))
pal <- setNames(pal, c("1", "2"))

#-----Cell's Spatial Pattern------
cell.type <- as.vector(cell.info[,1])
gg <- ggplot(cell.info, aes(x = X, y = Y, col = as.factor(cell.type), shape = as.factor(cell.type)))
pl <- gg + geom_point(size = 2.5) +
 scale_color_manual(values = c(pal[1], pal[2])) +
 theme_bw()+
 theme(legend.text=element_text(size=20),
       axis.title.x=element_text(size=16),
       axis.title.y=element_text(size=16),
       axis.text.x = element_text(size = 12,face = "bold"),
       axis.text.y = element_text(size = 12,face = "bold")
 ) + labs(x = "H", y = "L") +
 guides(color = guide_legend(title = "Cell Type",
                             title.theme = element_text(size = 25),
                             override.aes = list(size = 5)
 ),
 shape = guide_legend(title = "Cell Type",
                      title.theme = element_text(size = 25),
                      override.aes = list(size = 5)))
ggsave("cell spatial.png", pl, width = 9, height = 12)


#----run CSSNEst--------
nu <- rep(2*G, n)
Result <- CSSNEst(X, cell.info, nu = nu, d = 0.1, m.info = 70, is.scale = TRUE)
Sparse.Corr <- Result$`Gene Networks`

#-----The first five cell's estimated gene co-expression networks-----
colors_func <- colorRampPalette(c('white', "black"))
colors <- colors_func(2)
filename <- paste0("Est_", 1:5, ".png")
for(i in 1:10){
 p2 <- pheatmap(Sparse.Corr[,,i],
                color = colors,
                legend_breaks = c(0,1),
                cluster_cols = F, cluster_rows = F,
                show_rownames = F, show_colnames = F,
                width = 3.3, height = 2.8,
                filename = filename[i]
                
 )
}
# Prediction
set.seed(1)
miss.num <- 5
miss.x <- runif(miss.num, min(cell.info[,2]), max(cell.info[,2]))
miss.y <- runif(miss.num, min(cell.info[,3]), max(cell.info[,3]))
miss.indx <- cbind(miss.x, miss.y)
pre <- CSSNPredict(Sparse.Corr, cell.info, miss.indx)

}
\references{

}
