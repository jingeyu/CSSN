## CSTS: bayesian joint modeling of single-cell expression data and bulk spatial transcriptomic data

The R package CSTS implements the two-step method proposed by Jinge Yu and Xiangyu Luo (2021+) that recovers spatially-varying cell-specific gene co-expression networks for single-cell spatial expression data. In this package, we employed a two-step algorithm to perform efficient thresholded Bayesian posterior estimates for gene co-expression networks in each cell and predict gene co-expression networks where cells are missing. CSTS can be installed in commonly used operating systems including Windows, Linux and Mac OS. 


## Prerequisites and Installation

1. R version >= 3.6.
2. R packages: pheatmap (>= 1.0.12), igraph (>= 1.2.6), plotly (>= 4.9.2.1)
3. Install the package CSTS.

```
devtools::install_github("jingeyu/CSTS")
```


## Example Code
Following shows an example that generates data and runs the main function "CSTSEst" and "CSTSPredict" in our package. 

``` {r, eval=FALSE}
library(CSTS)
#############################################
#Data generation
#############################################
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
  # Remove the i-th cell
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
  Corr.true[,,i] <- diag(diag(Sigma.i[,,i])^(-0.5)) %*% Sigma.i[,,i] %*% diag(diag(Sigma.i[,,i])^(-0.5))
  X[,i] <- mvrnorm(1, mu = rep(0, G), Sigma = Sigma.i[,,i])
}

#----run CSTS--------
library("CSTS")
nu <- rep(2*G, n)
Result <- CSTSEst(X, cell.info, nu = nu, d = 0.1, m.info = 70, is.scale = TRUE)
Sparse.Corr <- Result$`Gene Networks`
Corr.true[Corr.true != 0 ] <- 1

error <- 0
for(i in 1:n){
  error <- c(error, sum(abs(Corr.true[,,i][upper.tri(Corr.true[,,i])] - Sparse.Corr[,,i][upper.tri(Sparse.Corr[,,i])])))
}

# Error:
print(sum(error)/n)

```

``` {r, eval=FALSE}
example(CSTSEst)
set.seed(1)
miss.num <- 5
miss.x <- runif(miss.num, min(cell.info[,2]), max(cell.info[,2]))
miss.y <- runif(miss.num, min(cell.info[,3]), max(cell.info[,3]))
miss.indx <- cbind(miss.x, miss.y)
library(CSTS)
pre <- CSTSPredict(Sparse.Corr, cell.info, miss.indx)
```
or you can simply run
``` {r, eval=FALSE}
library(CSTS)
example(CSTSEst)
example(CSTSPredict)
```

## Remarks
* If you have any questions regarding this package, please contact Jinge Yu at yjgruc@ruc.edu.cn.
