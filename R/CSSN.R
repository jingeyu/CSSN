CSSNEst <- function(X, cell.info, nu, d = 0.1,
                    m.info = 70, is.scale = TRUE,
                    is.all = TRUE, indx.cell = NULL, 
                    output.corr = FALSE){
  library(Matrix)
  #---- Step0 Preprocess of Data ------
  # if there is NaNs in the dataset
  if(sum(is.na(X)) > 0){
    stop("Please remove NaNs in the dataset!!!")
  }

  if(d < 0 & d > 1){
    stop("Please set threshold parameter d between 0 and 1!")
  }

  if(d < 0 & d > 1){
    stop("Please set threshold parameter d between 0 and 1!")
  }

  # Cell number:
  n <- ncol(X)
  # Gene number
  G <- nrow(X)
  # Notice that cell types in cell.info are factor
  cell.type <- as.vector(cell.info[, 1])
  # Cell Type number
  K <- length(unique(cell.type))
  # Transform factor cell types to numeric ones.
  ct <- names(table(cell.type))
  for(i in 1:length(ct)){
    cell.type <- gsub(ct[i], i, cell.type)
  }
  cell.type <- as.numeric(cell.type)
  ind.cell.type <- list()

  #---- Step1 Centralization and Scaling ----
  for(k in 1:K){
    tmp <- which(cell.type == k)
    X[, tmp] <- X[, tmp] - rowMeans(X[,tmp])
    if(is.scale == TRUE){
      X[,tmp] <- X[, tmp] / apply(X[, tmp], 1, sd)
    }
    ind.cell.type[[k]] <- tmp
  }

  NA.num <- sum(is.na(X))
  if(NA.num > 0){
    stop("Please remove genes having zero values across one cell-type!")
  }

  cell.info <- as.data.frame(cell.info)
  colnames(cell.info) <- c("CT", "X", "Y")
  #---- Step2 Estimate Sigma_k ----
  Sgm.hat <- array(NA, dim = c(G,G,K))
  for(k in 1:K){
    nk <- length(ind.cell.type[[k]])
    Sgm.hat[,,k] <- X[, ind.cell.type[[k]]] %*% t(X[, ind.cell.type[[k]]]) / (nk - 1)
  }


  #---- Step3 Find the Neighbors of Cell i ----
  #--- Set neighborhood radius r ---
  # according to average density of cells in the whole tissue slice ---

  # Area of cell tissue slice
  W <- max(cell.info[,2]) - min(cell.info[,2])
  L <- max(cell.info[,3]) - min(cell.info[,3])

  S <- L * W
  # Square neighborhood radius
  r <- sqrt(S * m.info / n) / 2

  NeiFind <- function(i){
    ind.i <- as.numeric(cell.info[i,2:3])
    nei.inx <- which(abs(cell.info[,2] - ind.i[1]) < r & abs(cell.info[,3] - ind.i[2]) < r)
    # Remove the i-th cell
    cell.index <- nei.inx[nei.inx != i]
    cell.type.nei <- cell.type[cell.index]
    # Return a matrix with index of neighborhood cells and corresponding cell type.
    return(cbind(cell.index, cell.type.nei))
  }


  #---- Step4 Estimate Posterior Expectation of Sigma_i ----
  SgmEst <- function(nu.i, i, K, Sgm.hat){
    nei.mat <- data.frame(NeiFind(i))
    colnames(nei.mat) <- c("cell.index", "cell.type")

    ni <- nrow(nei.mat)
    if(ni == 0){
      Lambda.i <- (nu.i - G - 1) * Sgm.hat[,,cell.type[[i]]]
    }else{
      cell.label <- as.integer(names(table(nei.mat$cell.type)))
      nei.nk <- as.numeric(table(nei.mat$cell.type))
      weight <- nei.nk / ni
      tmp <- 0
      for(j in 1:length(cell.label)){
        tmp <- tmp + Sgm.hat[,,cell.label[j]] * weight[j]
      }
      Lambda.i <- (nu.i - G - 1) * tmp
    }
    # mu is known to be 0.
    # posterior mean
    S.i <- X[,i] %*% t(X[,i])
    Sgm.i.hat <- (S.i + Lambda.i) / (nu.i - G)

    return(Sgm.i.hat)
  }

  if(is.all == TRUE){
    if(output.corr == TRUE){
      Corr.ori <- list()
      Corr.sparse  <- list()
      for(i in 1:n){
        tmp <- SgmEst(nu[i], i, K, Sgm.hat)
        tmp <- cov2cor(tmp)
        Corr.ori[[i]] <- tmp
        tmp[abs(tmp) < 0.1] <- 0
        tmp[tmp != 0] <- 1
        Corr.sparse[[i]] <- Matrix(data = tmp, sparse = TRUE)
      }
      return(list("Correlation Matrix" = Corr.ori,
                  "Gene Networks" = Corr.sparse))
    } else{
      Corr.sparse <- list()
      for(i in 1:n){
        tmp <- SgmEst(nu[i], i, K, Sgm.hat)
        tmp <- cov2cor(tmp)
        tmp[abs(tmp) < 0.1] <- 0
        tmp[tmp != 0] <- 1
        Corr.sparse[[i]] <- Matrix(data = tmp, sparse = TRUE)
      }
      return(Corr.sparse)
    }
   
  } else{
    n2 <- length(indx.cell)
    if(output.corr == TRUE){
      Corr.ori <- list()
      Corr.sparse  <- list()
      for(i in 1:n2){
        tmp <- SgmEst(nu[indx.cell[i]], indx.cell[i], K, Sgm.hat)
        tmp <- cov2cor(tmp)
        Corr.ori[[i]] <- tmp
        tmp[abs(tmp) < 0.1] <- 0
        tmp[tmp != 0] <- 1
        Corr.sparse[[i]] <- Matrix(data = tmp, sparse = TRUE)
      }
      return(list("Correlation Matrix" = Corr.ori,
                  "Gene Networks" = Corr.sparse))
    } else{
      Corr.sparse  <- list()
      for(i in 1:n2){
        tmp <- SgmEst(nu[indx.cell[i]], indx.cell[i], K, Sgm.hat)
        tmp <- cov2cor(tmp)
        tmp[abs(tmp) < 0.1] <- 0
        tmp[tmp != 0] <- 1
        Corr.sparse[[i]] <- Matrix(data = tmp, sparse = TRUE)
      }
      return(Corr.sparse)
  }

  }
}


CSSNPredict <- function(GN, cell.info, miss.indx, m.info = 70){
  G <- nrow(GN[[1]])
  miss.num <- nrow(miss.indx)
  n <- nrow(cell.info)
  n1 <- length(GN)
  # if the input gene co-expression matrix is of all cells.
  if(n != n1){
    stop("Please input the gene co-expression matrix of all cells!")
  }
  
  if(min(miss.indx[,1]) < min(cell.info[,2]) | min(miss.indx[,2]) < min(cell.info[,3]) | max(miss.indx[,1]) > max(cell.info[,2]) | max(miss.indx[,1]) > max(cell.info[,2])){
    stop("Please remove cells out of the range of tissue slice!")
  }
  if(sum(miss.x %in% cell.info[,2]) >0 & sum(miss.y %in% cell.info[,3]) > 0){
    stop("Please remove cells already in the tissue slice!")
  }

  W <- max(cell.info[,2]) - min(cell.info[,2])
  L <- max(cell.info[,3]) - min(cell.info[,3])
  S <- L * W
  r <- sqrt(S * m.info / n) / 2

  NeiFind <- function(miss.indx.cell){
    nei.indx <- which(abs(cell.info[,2] - miss.indx.cell[1]) < r & abs(cell.info[,3] - miss.indx.cell[2]) < r)
    cell.type.nei <- cell.type[nei.indx]
    return(cbind(nei.indx, cell.type.nei))
  }

  est.miss <- list()
  for(i in 1:miss.num){
    miss.nei <- NeiFind(miss.indx[i,])
    tmp <- matrix(0, G, G)
    for(j in 1:nrow(miss.nei)){
      tmp <- tmp + GN[[miss.nei[j,1]]]
    }
    tmp <- tmp / miss.num
    tmp[tmp < 0.5] <- 0
    tmp[tmp >= 0.5] <- 1
    est.miss[[i]] <- tmp
  }
  return(est.miss)
}
