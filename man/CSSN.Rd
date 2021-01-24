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
library(MASS)
library(Matrix)
library(CholWishart)
library(pheatmap)
library(ggplot2)

# set seed for reproducible
set.seed(20210120)
# gene number
G <- 10
# cell-type number
K <- 2
# radius of square neighborhood
r <- 50
# the area of whole tissue cryosection
L <- 100 # length, micrometre
H <- 50 # width, micrometre
area  <- L * H

#---- Step1 generate spatial coordinates of cells ----
# Poisson process
lambda <- c(2,4) / 10
# 1. divide the area R into n.region(100) sub-regions
sub.num <- 10
n.region <- sub.num * sub.num
N <- matrix(NA, K, n.region)
for(i in 1:n.region){
 N[, i] <- rpois(K, lambda)
}
# total cell numbers in one tissue cryosection 
n <- sum(N)
# Random and independently place N[,i] points in n.region sub-regions.
N.i <- colSums(N)
xlab <- seq(0, H, length.out = sub.num + 1)
ylab <- seq(0, L, length.out =sub.num + 1)
centroid.x <- NULL
centroid.y <- NULL
for(i in 1:sub.num){
 for(j in 1:sub.num){
   centroid.x <- c(centroid.x, runif(N.i[(i-1) * sub.num + j], xlab[i], xlab[i+1]))
   centroid.y <- c(centroid.y, runif(N.i[(i-1) * sub.num + j], ylab[j], ylab[j+1]))
 }
}

#---- set spatial pattern manually----
pal <- c(rgb(221, 160, 221, maxColorValue = 255), 
        rgb(0, 206, 209, maxColorValue = 255))
pal <- setNames(pal, c("1", "2"))
cell.type <- NULL
ct <- 1:K
for(i in 1:n.region){
 cell.type <- c(cell.type, rep(ct, N[,i]))
}

cell.type[1:10] <- sample(1:K, 10, prob = c(0.2, 0.8), replace = TRUE)

t1 <- which(centroid.x < 20 & centroid.x > 0 & centroid.y < 100 & centroid.y > 50)
t2 <- which(centroid.x < 50 & centroid.x > 20 & centroid.y < 50 & centroid.y > 0)

cell.type[t1] <- sample(c(1,2), length(t1), prob = c(0.1, 0.9), replace = TRUE)
cell.type[t2] <- sample(c(1,2), length(t2), prob = c(0.8, 0.2), replace = TRUE)

cell.info <- data.frame(cell.type, centroid.x, centroid.y)
colnames(cell.info) <- c("CT", "X", "Y")

#-----Cell's Spatial Pattern------
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

#---- Step 2 Generate true Sigma.k (Sparse + Positive define) ----
# K blocks
Indicator <- function(i, j, number){
 if(abs(i - j) == number){
   return(1)
 }else{
   return(0)
 }
}
Sigma.k <- array(0, dim = c(G, G, K))
sub.piece <- G / K
subSgm.1 <- matrix(0, sub.piece, sub.piece)
subSgm.2 <- matrix(0, sub.piece, sub.piece)

rho <- 0.7
for(i in 1:sub.piece){
 for(j in 1:sub.piece){
   subSgm.1[i, j] <- rho^(abs(i - j))
   subSgm.2[i, j] <- -0.3 * Indicator(i, j, 1) + 1.3 * Indicator(i, j, 0)
 }
}

diag(subSgm.1) <- diag(subSgm.1) + 0.5
Sigma.k[,,1] <- as.matrix(bdiag(list(subSgm.1, subSgm.2)))
Sigma.k[,,2] <- as.matrix(bdiag(list(subSgm.2, subSgm.1)))

nu <- rep(G + 5, n)
NeiFind <- function(i, r){
 ind.i <- as.numeric(cell.info[i,2:3])
 nei.inx <- which(abs(cell.info[,2] - ind.i[1]) < r & abs(cell.info[,3] - ind.i[2]) < r)
 #Remove the i-th cell
 cell.index <- nei.inx[nei.inx != i]
 cell.type <- cell.info[cell.index, 1]
 # Return a matrix with index of neighborhood cells and corresponding cell type.
 return(cbind(cell.index, cell.type))
}


ExpSigma <- function(i, r, nu.i){
 nei.mat <- data.frame(NeiFind(i, r))
 colnames(nei.mat) <- c("cell.index", "cell.type")
 ni <- nrow(nei.mat)
 if(ni == 0){ 
   Lambda.i <- (nu.i - G - 1) * Sigma.k[,,cell.type[[i]]]
 }else{
   cell.label <- as.integer(names(table(nei.mat$cell.type)))
   nei.nk <- as.numeric(table(nei.mat$cell.type))
   weight <- nei.nk / ni
   tmp <- 0
   for(j in 1:length(cell.label)){
     tmp <- tmp + Sigma.k[,,cell.label[j]] * weight[j]
   }
   Lambda.i <- (nu.i - G - 1) * tmp 
 }
}

Sigma.i <- array(NA, dim = c(G, G, n))
X <- matrix(NA, G, n)
Lambda <- array(NA, dim = c(G, G, n))
c.thre <- 0.5
Corr.true <- array(NA, dim = c(G, G, n))
for(i in 1:n){
 Lambda[,,i] <- ExpSigma(i, r, nu[i])
 Sigma.i[,,i] <- rInvWishart(1, nu[i], Lambda[,,i])[,,1]
 Sigma.i[,,i][abs(Sigma.i[,,i]) < c.thre] <- 0
 # ensure Sigma_i are positive definite
 diag(Sigma.i[,,i]) <- diag(Sigma.i[,,i]) + 5
 Corr.true[,,i] <- diag(diag(Sigma.i[,,i])^(-0.5)) \%*\% Sigma.i[,,i] \%*\% diag(diag(Sigma.i[,,i])^(-0.5))
 X[,i] <- mvrnorm(1, mu = rep(0, G), Sigma = Sigma.i[,,i])
}


#----run CSSNEst--------
library(CSSN)
nu <- rep(2*G, n)
Result <- CSSNEst(X, cell.info, nu = nu, d = 0.1, m.info = 70, is.scale = TRUE)
Sparse.Corr <- Result$`Gene Networks`
Corr.true[Corr.true != 0 ] <- 1

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
