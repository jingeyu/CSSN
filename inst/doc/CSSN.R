## ----setup, echo=FALSE, results="hide"----------------------------------------
knitr::opts_chunk$set(tidy = FALSE,
                      cache = FALSE,
                      dev = "png",
                      message = FALSE, error = FALSE, warning = TRUE)

## ---- loadLibraries-----------------------------------------------------------
library(CSSN)
library(stringr)
library(ggplot2)
library(pheatmap)

## ----echo = TRUE--------------------------------------------------------------
# example_data includes two matrices, X and cell.info.
dir <- system.file(package = 'CSSN') #directory for the example data
load(paste0(dir,"/example_data.RData"))
# X is the single-cell spatial expression data, continuous.
dim(X)
# cell.info is the matrix contains each cell's cell type information and their 2-D coordinates.
dim(cell.info)

# Look at the expression data.
head(X[,1:5])
# Look at the cell information matrix.
head(cell.info)

# Cell type number
table(cell.info[,1])
# number of cells
n <- ncol(X)
# number of genes
G <- nrow(X)
print(c(n, G))

## ----Visualiztion, warning=FALSE----------------------------------------------
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
 ) + labs(x = "X", y = "Y") +
 guides(color = guide_legend(title = "Cell Type",
                             title.theme = element_text(size = 25),
                             override.aes = list(size = 5)
 ),
 shape = guide_legend(title = "Cell Type",
                      title.theme = element_text(size = 25),
                      override.aes = list(size = 5)))
pl

## ----Two-step-----------------------------------------------------------------
# Two-step algorithm
nu <- rep(2*G, n)
# Output all cells' gene co-expression networks.
Result <- CSSNEst(X, cell.info, nu = nu, d = 0.1, m.info = 70, is.scale = TRUE, is.all = TRUE)
# Result is a list of length n, where n is the number of cell.

#-----The first cell's estimated gene co-expression networks-----
colors_func <- colorRampPalette(c('white', "black"))
colors <- colors_func(2)
par(mfrow=c(1,2))
pheatmap(Result[[1]],
                color = colors,
                legend_breaks = c(0,1),
                cluster_cols = F, cluster_rows = F,
                show_rownames = F, show_colnames = F,
                width = 2.8, height = 2.8
                # filename = filename[i]
)

## -----------------------------------------------------------------------------
indx.cell <- c(1,3,7,10)
result <- CSSNEst(X, cell.info, nu = nu, d = 0.1, m.info = 70, is.scale = TRUE, is.all = FALSE, indx.cell = indx.cell, output.corr = TRUE)
result$`Correlation Matrix`[[1]]
result$`Gene Networks`[[1]]

## ----Prediction---------------------------------------------------------------
# generate coordinates of missing cells
set.seed(1)
miss.num <- 5
miss.x <- runif(miss.num, min(cell.info[,2]), max(cell.info[,2]))
miss.y <- runif(miss.num, min(cell.info[,3]), max(cell.info[,3]))
miss.indx <- cbind(miss.x, miss.y)

# Apply prediction algorithm
pre <- CSSNPredict(Result, cell.info, miss.indx)
length(pre)
dim(pre[[1]])

