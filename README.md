## CSSN: Recovering spatially-varying cell-specific gene co-expression networks for single-cell spatial expression data

The R package CSSN implements the two-step method proposed by Jinge Yu and Xiangyu Luo (2021+) that recovers spatially-varying cell-specific gene co-expression networks for single-cell spatial expression data. In this package, we employed a two-step algorithm to perform efficient thresholded Bayesian posterior estimates for gene co-expression networks in each cell and predict gene co-expression networks where cells are missing. CSSN can be installed in commonly used operating systems including Windows, Linux and Mac OS. 


## Prerequisites and Installation

1. R version >= 3.6.
2. R packages: pheatmap (>= 1.0.12), stats (>= 4.0.3)
3. Install the package CSSN.

```
devtools::install_github("jingeyu/CSSN")
```


## Example Code
Following shows an example that generates data and runs the main function "CSSNEst" and "CSSNPredict" in our package. 

``` {r, eval=FALSE}
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
```
or you can simply run
``` {r, eval=FALSE}
library(CSSN)
example(CSSNEst)
```

## Remarks
* If you have any questions regarding this package, please contact Jinge Yu at yjgruc@ruc.edu.cn.

