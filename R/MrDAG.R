# 135 characters #####################################################################################################################
#' @title MrDAG: Mendelian randomization for multiple outcomes and exposures with Bayesian Directed Acyclic Graphs (DAGs) 
#' exploration and causal effects estimation
#' @description Markov chain Monte Carlo (MCMC) implementation of Bayesian multivariable, multi-response Mendelian randomization (MR) 
#' model for summary-level data with Directed Acyclic Graphs (DAGs) exploration and causal effects estimation
#'
#' @param data Number of observations (IVs in a summary-level MR design) times the number of traits (both outcomes and exposures). 
#' There is a restriction in the order of appearance of the traits in the \code{data} matrix: The group of outcomes ought to appear 
#' first, followed by the group of exposures. See also \code{MrDAGcheck} argument
#' @param niter Number of MCMC iterations (excluding burn-in)
#' @param burnin Number of MCMC iterations to be performed during burn-in
#' @param thin Parameter that defines how often the MCMC output should be stored, i.e., at every thin-th iteration
#' @param w Prior probability of edge inclusion (\code{0.05} default) between each pair of nodes (vertices) in the graph
#' @param a Degrees of freedom of the Wishart prior distribution on the precision matrix, i.e., the inverse of the covariance matrix 
#' between the traits (\code{m} default, where \code{m} is the total number of outcomes and exposures)
#' @param U (Proportional to the) expected value of the Wishart prior distribution on the precision matrix, i.e., the inverse of the
#' covariance matrix between the traits. Specifically, proportional to an m-dimensional diagonal matrix as default, where \code{m} is 
#' the total number of outcomes and exposures
#' @param tempmax Annealing parameter \code{T} used to facilitate the convergence of the MCMC algorithm to the target distribution 
#' (\code{10} default). Temperature 1/T increases linearly during the burn-in until T=1 at the end of the burn-in
#' @param MrDAGcheck List object that contains the indices of the traits in \code{data} that are defined as outcomes and exposures.
#' If \code{NULL} (default), no partition of the traits between outcomes and exposures is specified and \code{MrDAG} algorithm 
#' performs structure learning between all traits without constraints. Note that if so, this procedure does not correspond to an MR 
#' analysis
#' @param filename Name of the file for \code{MrDAG} output object (\code{"MrDAG_object"} default). If \code{NULL}, the ouput is not
#' saved
#' @param filepath Directory where \code{MrDAG} object is saved. If the directory does not exist, \code{MrDAG} algorithm will create 
#' it by using the current working directory as the root directory. If the directory is not specified (\code{NULL} default), 
#' \code{MrDAG} object is saved in the current working directory
#' @param fastMCMC If logical \code{TRUE} (default), the first valid DAG proposed by randomly selecting either to add, delete or 
#' swap a directed edge from the current DAG is used in the proposal distribution. Otherwise, all possible DAGs which differ from the 
#' current one by adding, deleting or swapping a directed edge are generated and, then, a DAG is sampled at random for the proposal 
#' distribution
#' @param savememory If logical \code{TRUE} (default), the visited graph and the posterior draws from the modified Cholesky 
#' decomposition (L,D) (\insertCite{Zuber2024;textual}{MrDAG} and \insertCite{Castelletti2021;textual}{MrDAG}) are stored as a list 
#' object, otherwise they are stored as an array
#' @param seed Seed to be used in the initialisation of the MCMC algorithm (\code{3112021} default)
#'
#' @details
#' For details regarding the model and the algorithm, see \insertCite{Zuber2024;textual}{MrDAG}
#'
#' @export
#'
#' @return The value returned is a list object \code{list(graphs, L,D, logmarglik, validpropMrDAG, acceptpropDAG, timeMrDAG, 
#' hyperpar, samplerpar, opt)}
#' \itemize{
#'   \item{\code{graphs}}{ Explored DAGs that belong to learned Completed Partially DAGs (CPDAGs) 
#'         (\insertCite{Chickering2002;textual}{MrDAG}) or Essential Graphs (EGs) (\insertCite{Andersson1997;textual}{MrDAG}). The 
#'         class of \code{graphs} depends on the \code{savememory} option. The number of explored DAGs corresponds to the number of 
#'         thinned (\code{thin}) MCMC iterations (excluding burn-in) }
#'   \item{\code{L}}{ Posterior samples of the lower triangular matrix of the modified Cholesky decomposition (L,D) 
#'         (\insertCite{Zuber2024;textual}{MrDAG} and \insertCite{Castelletti2021;textual}{MrDAG}) }
#'   \item{\code{D}}{ Posterior samples of the diagonal matrix of the modified Cholesky decomposition (L,D) 
#'         (\insertCite{Zuber2024;textual}{MrDAG} and \insertCite{Castelletti2021;textual}{MrDAG}) }
#'   \item{\code{logmarglik}}{ Log-marginal likelihood of explored DAGs which belong to the Markov Equivalent Classes whose unique 
#'         representative chain graphs are the EGs learned during MCMC iterations (including burn-in) without thinning }
#'   \item{\code{validpropMrDAG}}{ If \code{MrDAGcheck} is different from \code{NULL}, the proportion of proposed DAGs that comply 
#'         with the partial ordering (\insertCite{Perkovic2017;textual}{MrDAG} implied by the partition of the traits between 
#'         exposures and outcomes with directed causal effects only from the former to the latter }
#'   \item{\code{acceptpropDAG}}{ Proportion of proposed DAGs that are accepted by the Metropolis-Hastings ratio after 
#'         burn-in }
#'   \item{\code{timeMrDAG}}{ Time in minutes employed by \code{MrDAG} algorithm to analyse the data. Considerable gains are obtained 
#'         by specifying \code{fastMCMC = TRUE} option }
#'   \item{\code{hyperpar}}{ List of the hyper-parameters \code{list(w, a, U, tempmax)} and, if specified, \code{MrDAGcheck} list }
#'   \item{\code{samplerpar}}{ List of parameters used in the MCMC algorithm \code{list(niter, burnin, thin)} }
#'   \item{\code{opt}}{ List of options used in \code{MrDAG} algorithm \code{list(fastMCMC, savememory, seed)} } }
#'
#' @references
#' \insertAllCited{}
#'
#'
#' @examples
# 100 characters ##################################################################################
#' # Example: Analysis of lifestyle and behavioural exposures that might impact mental health 
#' # phenotypes. 708 independent Instrumental Variables (IVs) were selected to be associated at 
#' # genome-wide significance with six lifestyle and behavioural traits after pruning or clumping,
#' # which are considered exposures of the risk of seven mental health outcome phenotypes
#'
#' # After loading the data set
#'
#' data(LBT2MD_data)
#'
#' # the indices of the traits that define the outcomes and exposures are provided in a list object
#'
#' MrDAGcheck <- NULL
#' MrDAGcheck$Y_idx <- 1 : 7    # Mental health phenotypes
#' MrDAGcheck$X_idx <- 8 : 13   # Lifestyle and behavioural traits
#'
#' # MrDAG algorithm is run to generate 1,000 posterior samples of all unknowns
#'
#' output <- MrDAG(data = LBT2MD_data, 
#'                 niter = 5000, burnin = 2500, thin = 5, tempmax = 20, w = 0.01, 
#'                 MrDAGcheck = MrDAGcheck, filename = NULL)
#'
#'
# @importFrom ggpattern geom_tile_pattern scale_pattern_manual
# @importFrom ggplot2 aes element_blank element_text geom_hline geom_text geom_vline ggplot guides rel scale_fill_gradient2 theme 
#             unit xlab ylab
# @importFrom gRbase is.DAG
# @importFrom grDevices dev.off jpeg pdf png tiff
# @importFrom mvtnorm rmvnorm
# @importFrom pcalg dag2cpdag
# @importFrom qgraph qgraph
# @importFrom reshape2 melt
#' @importFrom Rdpack reprompt
# @importFrom stats quantile rbinom rgamma runif
# @importFrom utils capture.output globalVariables setTxtProgressBar txtProgressBar


########## MrDAG ##########

MrDAG <- function(data, 
                  niter, burnin, thin, 
                  w = 0.05, a = NULL, U = NULL, tempmax = 10, 
                  MrDAGcheck = NULL, 
                  filename = "MrDAG_object", filepath = NULL, 
                  fastMCMC = TRUE, savememory = FALSE, 
                  seed = 31122021)
{
  cat("##\n")
  cat("## Bayesian Mendelian randomization with Directed Acyclic Graphs exploration\n")
  cat("## Version: 0.1.0\n")
  cat("##\n")
  cat("## Copyright (C) 2024 L. Bottolo and V. Zuber\n")
  cat("##\n")
  cat("## The model selected is:\n")
  cat("##\n")
  
  if (is.null(MrDAGcheck))
  {
    
    if (fastMCMC == FALSE)
    {
      
      if (savememory == FALSE)
      {
        cat("## Directed Acyclic Graphs exploration\n")
      } else {
        cat("## Directed Acyclic Graphs exploration with memory saving\n")
      }
    
    } else {
      
      if (savememory == FALSE)
      {
        cat("## Directed Acyclic Graphs exploration with fast MCMC\n")
      } else {
        cat("## Directed Acyclic Graphs exploration with fast MCMC and memory saving\n")
      }
    }
  
  } else {
    
    if (fastMCMC == FALSE)
    {
      
      if (savememory == FALSE)
      {
        cat("## Mendelian randomization with Directed Acyclic Graphs exploration\n")
      } else {
        cat("## Mendelian randomization with Directed Acyclic Graphs exploration with memory saving\n")
      }
    
    } else {
      
      if (savememory == FALSE)
      {
        cat("## Mendelian randomization with Directed Acyclic Graphs exploration with fast MCMC\n")
      } else {
        cat("## Mendelian randomization with Directed Acyclic Graphs exploration with fast MCMC and memory saving\n")
      }
    
    }
  
  }
  
  cat("##\n")
  
  start_time <- Sys.time()
  set.seed(seed)
  
  X <- data
  n <- dim(X)[1]
  m <- dim(X)[2]
  
  X <- scale(X, center = TRUE, scale = FALSE)
  
  if (is.null(colnames(X)))
  {
    
    if (is.null(MrDAGcheck))
    {
      colnames(X) <- paste0("EO", seq(1 : m))
    } else {
      colnames(X) <- c(paste0("O", seq(1 : length(MrDAGcheck$Y_idx))), 
                       paste0("E", seq(1 : length(MrDAGcheck$X_idx))))
    }
  
  } else {
    colnames_YX <- toupper(substr(colnames(X), 1, 3))
    
    if (length(colnames_YX) != length(unique(colnames_YX)))
    {
      idx <- which(duplicated(colnames_YX))
      idx <- unique(c(idx, match(colnames_YX[idx], colnames_YX)))
      colnames_YX[idx] <- toupper(colnames(X)[idx])
    }
    
    colnames(X) <- colnames_YX
  }
  
  if (is.null(rownames(X)))
  {
    rownames(X) <- paste0("IV", seq(1 : n))
    rownames_X <- rownames(X)
  }
  
  if (is.null(a))
  {
    a <- m
  }
  
  if (is.null(U))
  {
    U <- diag(1, m) / n
  }
  
  hyperpar <- NULL
  hyperpar$w <- w
  hyperpar$a <- a
  hyperpar$U <- U
  rownames(hyperpar$U) <- colnames_YX
  colnames(hyperpar$U) <- colnames_YX
  hyperpar$tempmax <- tempmax
  hyperpar$MrDAGcheck <- MrDAGcheck
  
  samplerpar <- NULL
  samplerpar$niter <- niter
  samplerpar$burnin <- burnin
  samplerpar$thin <- thin
  
  tXX <- t(X)%*%X
  
  niter_TMP <- niter
  niter <- burnin + niter
  iterCounter <- 1
  
  validpropMrDAG <- 0
  acceptpropDAG <- 0
  
  ########## Storing ##########
  
  if (savememory == TRUE)
  {
    graphs <- vector("double", niter_TMP / thin)
    L <- vector("double", niter_TMP / thin)
    D <- vector("double", niter_TMP / thin)
    names(graphs) <- paste0("iter", seq(1, niter_TMP, thin))
    names(L) <- paste0("iter", seq(1, niter_TMP, thin))
    names(D) <- paste0("iter", seq(1, niter_TMP, thin))
  } else {
    graphs <- array(0, dim = c(m, m, niter_TMP / thin), dimnames = list(colnames_YX, colnames_YX, 
                                                                        paste0("iter", seq(1, niter_TMP, thin))))
    L <- array(0, dim = c(m, m, niter_TMP / thin), dimnames = list(colnames_YX, colnames_YX, 
                                                                   paste0("iter", seq(1, niter_TMP, thin))))
    D <- array(0, dim = c(m, m, niter_TMP / thin), dimnames = list(colnames_YX, colnames_YX, 
                                                                   paste0("iter", seq(1, niter_TMP, thin))))
  }
  
  Temp <- c(seq(tempmax, 1, length.out = burnin), rep(1, niter_TMP))
  logmarglik <- vector("double", niter)
  names(logmarglik) <- paste0("iter", seq(1, niter, 1))
  
  currentDAG <- matrix(0, ncol = m, nrow = m, dimnames = list(colnames_YX, colnames_YX))
  PBText <- utils::txtProgressBar(min = 2, max = niter, style = 3)
  
  ########## MCMC iterations ##########
  
  for (iter in 1 : niter)
  {
    proposed <- proposeDAG(currentDAG, MrDAGcheck, fastMCMC)
    
    if (iter > burnin & !is.null(MrDAGcheck))
    {
      validpropMrDAG <- validpropMrDAG + proposed$validpropMrDAG
    }
    
    acceptReject <- acceptRejectDAG(tXX, currentDAG, 
                                    proposed$proposedDAG, proposed$opNode, 
                                    proposed$opType, proposed$currentOpCard, proposed$proposedOpCard, 
                                    n, hyperpar$a, hyperpar$U, hyperpar$w, Temp[iter])
    
    if (acceptReject$isAccepted == TRUE)
    {
      currentDAG <- proposed$proposedDAG
      logmarglik[iter] <- acceptReject$logmarglik
      
      if (iter > burnin)
      {
        acceptpropDAG <- acceptpropDAG + 1
      }
    
    } else {
      
      if (iter == 1)
      {
        
        for (node in 1 : nrow(U))
        {
          logmarglik <- logmarglik + mlnodeDAGWishart(node, currentDAG, tXX, n, a, U)
        }
      
      } else {
        logmarglik[iter] <- logmarglik[iter - 1]
      }
    }
    
    postPar <- rDAGWishart(1, currentDAG, a + n, U + tXX)
    
    if (iter > burnin & iter %% thin == 0)
    {
      
      if (savememory == TRUE) {
        graphs[iterCounter] <- bd_encode(currentDAG)
        L[iterCounter] <- bd_encode(postPar$L)
        D[iterCounter] <- bd_encode(postPar$D)
      } else {
        graphs[, , iterCounter] <- currentDAG
        L[, , iterCounter] <- postPar$L
        D[, , iterCounter] <- postPar$D
      }
      
      iterCounter <- iterCounter + 1
    }
    
    utils::setTxtProgressBar(PBText, iter)
    close(PBText)
  }
  
  cat("\n")
  
  if (!is.null(MrDAGcheck))
  {
    validpropMrDAG <- validpropMrDAG / niter_TMP
  } else {
    validpropMrDAG <- NULL
  }
  
  acceptpropDAG <- acceptpropDAG / niter_TMP
  
  opt <- list(fastMCMC = fastMCMC, savememory = savememory, seed = seed)
  end_time <- Sys.time()
  timeMrDAG <- as.numeric(difftime(time1 = end_time, time2 = start_time, units = "min"))
  timeMrDAG <- paste0(round(timeMrDAG, digits = 3), "m")
  
  output <- list(graphs = graphs, L = L,D = D, 
                 logmarglik = logmarglik, validpropMrDAG = validpropMrDAG, acceptpropDAG = acceptpropDAG, 
                 timeMrDAG = timeMrDAG, 
                 hyperpar = hyperpar, 
                 samplerpar = samplerpar, 
                 opt = opt)
  
  if (!is.null(filename))
  {
    saveas <- paste(filename, "RData", sep = ".")
    
    if (is.null(filepath))
    {
      filepath <- getwd()
    } else {
      
      if (!dir.exists(filepath))
      {
        dir.create(filepath, recursive = TRUE)
        cat("\nThe directory", paste0("\"", filepath, "\""), "does not exist. It is now created.\n")
      }
    
    }
    
    if (.Platform$OS.type == "windows")
    {
      filepath <- paste0(normalizePath(filepath), "\\")
    } else {
      filepath <- paste0(normalizePath(filepath), "/")
    }
    
    save(output, file = paste0(filepath, saveas))
  }
  
  return(output)
}

########## acceptRejectDAG ##########

acceptRejectDAG <- function(tXX, currentDAG, proposedDAG, 
                            node, opType, currentOpCard, proposedOpCard, 
                            n, a, U, w, Temp)
{
  
  if (sum(abs((currentDAG - proposedDAG))) > 0)
  {
    logpriorRatios <- c(log(w / (1 - w)), log((1 - w) / w), log(1))
    logpriorRatio <- logpriorRatios[opType]
    logproposalRatio <- log(currentOpCard) - log(proposedOpCard)
    
    if (opType != 3) {
      currentlogmarglik <- mlnodeDAGWishart(node, currentDAG, tXX, n, a, U)
      proposedlogmarglik <- mlnodeDAGWishart(node, proposedDAG, tXX, n, a, U)
    } else {
      currentlogmarglik <- mlnodeDAGWishart(node[1], currentDAG, tXX, n, a, U) + 
                           mlnodeDAGWishart(node[2], currentDAG, tXX, n, a, U)
      proposedlogmarglik <- mlnodeDAGWishart(node[1], proposedDAG, tXX, n, a, U) + 
                            mlnodeDAGWishart(node[2], proposedDAG, tXX, n, a, U)
    }
    
    acceptRatio <- min(0, 1 / Temp * (proposedlogmarglik - currentlogmarglik + logpriorRatio + logproposalRatio))
    isAccepted <- log(stats::runif(1)) < acceptRatio 
    
    if (isAccepted == TRUE)
    {
       logmarglik <- 0
       
       for (node in 1 : nrow(U))
       {
          logmarglik <- logmarglik + mlnodeDAGWishart(node, proposedDAG, tXX, n, a, U)
       }
    
    } else {
       logmarglik <- NULL
    }
  
  } else {
    isAccepted <- FALSE
    logmarglik <- NULL
  }
  
  return(list(isAccepted = isAccepted, logmarglik = logmarglik))
}

########## bd_decode ##########

bd_decode <- function(string, separator = ";") 
{
  vec4mat <- as.numeric(strsplit(string, separator)[[1]])
  q <- length(vec4mat)
  matrix(vec4mat, ncol = q)
  
  return(vec4mat)
}

########## bd_encode ##########

bd_encode <- function(matrix, separator = ";")
{
  matrix <- paste(matrix, collapse = separator)
  
  return(matrix)
}

########## causal_effect ##########

causal_effect <- function(response, target, L,D)
{
  y <- response
  L_I <- L
  L_I[, target] <- 0
  diag(L_I) <- 1
  Sigma_I <- solve(t(L_I)) %*% D %*% solve(L_I)
  effect <- sapply(target, function(x) Sigma_I[x, y] / Sigma_I[x, x])
  
  return(effect)
}

########## fa ##########

fa <- function(node, DAG)
{
  # DAG <- as(DAG, "matrix")
  pa <- which(DAG[, node] != 0)
  fa <- c(node, pa)
  
  return(fa)
}

########## mlnodeDAGWishart ##########

mlnodeDAGWishart <- function(node, DAG, tXX, n, a, U)
{
  j <- node
  pa <- pa(j, DAG)
  m <- ncol(tXX)
  a.star <- (a + length(pa) - m + 1)
  Upost <- U + tXX
  
  if (length(pa) == 0)
  {
    U_jj <- U[j, j]
    Upost_jj <- Upost[j, j]
    
    # prior.normcost <- 1 / gamma(a.star) * (1 / 2) * U_jj ^a.star
    # post.normcost <- 1 / gamma(a.star + n / 2) * (1 / 2) * Upost_jj ^(a.star + n / 2)
    prior.normcost <- -lgamma(a.star / 2) + a.star / 2 * log(U_jj / 2)
    post.normcost <- -lgamma(a.star / 2 + n / 2) + (a.star / 2 + n / 2) * log(Upost_jj / 2)
  } else {
    U_paj.j <- U[pa, j]
    U_jj <- U[j, j] - t(U_paj.j) %*% solve(U[pa,pa]) %*% U_paj.j
    Upost_paj.j <- Upost[pa, j]
    Upost_jj <- Upost[j, j] - t(Upost_paj.j) %*% solve(Upost[pa, pa]) %*% Upost_paj.j
    
    # prior.normcost <- 1 / gamma(a.star) * (1 / 2) * U_jj ^a.star * det(as.matrix(U[pa, pa])) ^(1 / 2)
    # post.normcost <- 1 / gamma(a.star + n / 2) * (1 / 2) * (Upost_jj) ^(a.star + n / 2) * det((as.matrix(Upost[pa, pa]))) ^(1 / 2)
    prior.normcost <- -lgamma(a.star / 2) + a.star / 2 * log(U_jj / 2) + 0.5 * log(det(as.matrix(U[pa, pa])))
    post.normcost <- -lgamma(a.star / 2 + n / 2) + (a.star / 2 + n / 2) * log(Upost_jj / 2) + 0.5 * log(det(as.matrix(Upost[pa, pa])))
  }
  
  # nodemll <- -n / 2 * log(2 * pi) + log(prior.normcost) - log(post.normcost)
  nodelml <- -n / 2 * log(2 * pi) + prior.normcost - post.normcost
  
  return(nodelml)
}

########## opcard ##########

opcard <- function(DAG)
{
  A <- DAG
  q <- ncol(A)
  A_na <- A
  diag(A_na) <- NA
  
  id_set <- c()
  dd_set <- c()
  rd_set <- c()
  
  set_id <- which(A_na == 0, TRUE)
  
  if(length(set_id) != 0)
  {
    id_set <- cbind(1, set_id)
  }
  
  set_dd <- which(A_na == 1, TRUE)
  
  if(length(set_dd != 0)){
    dd_set = cbind(2, set_dd)
  }
  
  set_rd <- which(A_na == 1, TRUE)
  
  if(length(set_rd != 0))
  {
    rd_set <- cbind(3, set_rd)
  }
  
  O <- rbind(id_set, dd_set, rd_set)
  op.cardvec <- vector(length = nrow(O))
  
  for (i in 1 : nrow(O))
  {
    op.cardvec[i] <- gRbase::is.DAG(operation(O[i,1], DAG, O[i, 2 : 3]))
  }
  
  op.card <- sum(op.cardvec)
  
  return(op.card)
}

########## operation ##########

operation <- function(op, A, nodes)
{
  x <- nodes[1]
  y <- nodes[2]
  
  if (op == 1)
  {
    A[x, y] <- 1
  }
  
  if (op == 2)
  {
    A[x, y] <- 0
  }
  
  if (op == 3)
  {
    A[x, y] <- 0
    A[y, x] <- 1
  }
  
  return(A)
}

########## pa ##########

pa <- function(node, DAG)
{
  # DAG <- as(DAG, "matrix")
  pa <- which(DAG[, node] != 0)
  
  return(pa)
}

########## proposeDAG ##########

proposeDAG <- function(DAG, MRCheck, fastMCMC)
{
  A <- DAG
  A_na <- A
  diag(A_na) <- NA
  validpropMrDAG <- NA
  
  ## Define the set of possible operations
  ## The cardinality of O will change depending on how many edges are present in the DAG
  id_set <- c()
  dd_set <- c()
  rd_set <- c()
  
  ## Set of nodes for id
  set_id <- which(A_na == 0, TRUE)
  if (length(set_id) != 0)
  {
    id_set <- cbind(1, set_id)
  }
  
  ## Set of nodes for dd
  set_dd <- which(A_na == 1, TRUE)
  if (length(set_dd != 0))
  {
    dd_set <- cbind(2, set_dd)
  }
  
  ## Set of nodes for rd
  set_rd <- which(A_na == 1, TRUE)
  if (length(set_rd != 0))
  {
    rd_set <- cbind(3, set_rd)
  }
  
  O <- rbind(id_set, dd_set, rd_set)
  
  ## Sample one random operation and verify it produces a DAG
  if (fastMCMC == FALSE)
  {
    proposedOpCardVec <- vector(length = nrow(O))
    
    for (i in 1 : nrow(O))
    {
      proposedOpCardVec[i] <- gRbase::is.DAG(operation(O[i, 1], DAG, O[i, 2 : 3]))
    }
    
    proposedOpCard <- sum(proposedOpCardVec)
    i <- sample(which(proposedOpCardVec), 1)
    ANext <- operation(O[i, 1], A, O[i, 2 : 3])
    
    if (!is.null(MRCheck))
    {
      
      if (sum(ANext[MRCheck$Y_idx, MRCheck$X_idx]) != 0)
      {
        ANext <- A
        validpropMrDAG <- 0
      } else {
        validpropMrDAG <- 1
      }
    
    }
    
    currentOpCard <- opcard(ANext)
  } else {
    repeat
    {
      i <- sample(nrow(O), 1)
      ANext <- operation(O[i, 1], A, O[i, 2 : 3])
      verify <- gRbase::is.DAG(ANext)
      
      if (verify == TRUE)
      {
        
        if (!is.null(MRCheck))
        {
          
          if (sum(ANext[MRCheck$Y_idx, MRCheck$X_idx]) == 0)
          {
            validpropMrDAG <- 1
            break
          } else {
            ANext <- A
            validpropMrDAG <- 0
            break
          }
        
        } else {
          break
        }
      
      }
    
    }
    
    proposedOpCard <- nrow(O)
    currentOpCard <- nrow(O)
  }
  
  opType <- O[i, 1]
  if (opType == 3) 
  {
    opNode <- O[i, -1]
  } else {
    opNode <- O[i, 3]
  }
  
  return(list(proposedDAG = ANext, opType = opType, opNode = opNode, 
              currentOpCard = currentOpCard, proposedOpCard = proposedOpCard, 
              validpropMrDAG = validpropMrDAG))
}

########## rDAG ##########

rDAG <- function(m, w)
{
  DAG <- matrix(0, nrow = m, ncol = m)
  DAG[lower.tri(DAG)] <- stats::rbinom(n = m * (m - 1) / 2, size = 1, prob = w)
  
  return(DAG)
}

########## rDAGWishart ##########

rDAGWishart <- function(n, DAG, a, U)
{
  m <- ncol(DAG)
  ajs <- sapply(1 : m, function(j) { a + sum(DAG[, j] == 1) - m + 1 })
  L.array <- array(0, dim = c(m, m, n))
  D.array <- array(0, dim = c(m, m, n))
  
  for (i in 1 : n)
  {
    params <- lapply(1 : m, function(j) { rnodeDAGWishart(j, DAG, ajs[j], U) })
    sigmas <- sapply(1 : m, function(x) { params[[x]]$sigmaj })
    L <- lapply(1 : m, function(x) { params[[x]]$Lj })
    D.array[, , i] <- diag(sigmas)
    
    for (j in 1 : m)
    {
      whc <- which(DAG[, j] == 1)
      L.array[whc, j, i] <- as.numeric(L[[j]])
    }
    
    diag(L.array[, , i]) <- 1
  }
  
  if (n == 1)
  {
    D.array <- D.array[, , 1]
    L.array <- L.array[, , 1]
  }
  
  return(list(D = D.array, L = L.array))
}

########## rnodeDAGWishart ##########

rnodeDAGWishart <- function(node, DAG, aj, U)
{
  j <- node
  pa <- pa(j, DAG)
  out <- list(sigmaj = 0, Lj = 0)

  if (length(pa) == 0)
  {
    U_jj <- U[j, j]
    out$sigmaj <- stats::rgamma(1, shape = aj / 2, rate = U_jj / 2) ^(-1)
  } else {
    U_paj.j <- U[pa, j]
    invU_papa <- solve(U[pa, pa])
    U_jj <- U[j, j] - t(U_paj.j) %*% invU_papa %*% U_paj.j
    out$sigmaj <- stats::rgamma(1, shape = aj / 2, rate = U_jj / 2) ^(-1)
    out$Lj <- mvtnorm::rmvnorm(1, -invU_papa %*% U_paj.j, out$sigmaj * invU_papa)
  }
  
  return(out)
}

