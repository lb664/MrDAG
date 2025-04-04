# 135 characters #####################################################################################################################
#' @title Estimation of the Posterior Probability of Edge Inclusion
#' @description Estimation of the Posterior Probability of Edge Inclusion (PPEI) based on the Directed Acyclic Graphs (DAGs) explored
#' by \code{MrDAG} algorithm
#'
#' @param output Output produced by \code{MrDAG} algorithm
#' @param ord Indices of the traits in \code{data} that define the outcomes and the exposures. They specify the order of appearance 
#' of the traits when printing the PPEI. If \code{ord = NULL} (default), the order is the same as \code{MrDAGcheck} in \code{MrDAG} 
#' output. Note that despite \code{ord} ordering, the outcomes always precede the exposures
#'
#' @export
#'
#' @return The value returned is a list object \code{list(edgeProb, group, ord)}
#' \itemize{
#'   \item{\code{edgeProb}}{ Estimate of the PPEI between each trait. Feedback loop of the same trait is not allowed (NA). If a
#'          partition of the traits between outcomes and exposures is specified, the PPEIs between outcomes and exposures are not
#'          considered (NA) }
#'   \item{\code{group}}{ If \code{MrDAGcheck} in \code{MrDAG} algorithm contains the indices of the traits in \code{data} that 
#'         define the outcomes and the exposures, \code{group} coincides with this list }
#'   \item{\code{ord}}{ Indices of the traits in \code{data} that specify the outcomes and the exposures. It might differ from 
#'         \code{group} if \code{ord} has been specified and it is different from \code{MrDAGcheck}} }
#'
# @references
# \insertAllCited{}
#'
#'
#' @examples
# 100 characters ##################################################################################
#' # Example: Estimation of the Posterior Probabilities of Edge Inclusion (PPEIs) of lifestyle and 
#' # behavioural traits, and mental health phenotypes. 708 independent Instrumental Variables 
#' # (IVs) were selected to be associated at genome-wide significance with six lifestyle and 
#' # behavioural traits after pruning or clumping which are considered exposures of the risk of 
#' # seven mental health outcome phenotypes
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
#'                 niter = 7500, burnin = 2500, thin = 5, tempMax = 20, pp = 0.01, 
#'                 MrDAGcheck = MrDAGcheck, fileName = NULL)
#'
#' # Finally, PPEIs are calculated and presented with lifestyle and behavioural traits first, 
#' # followed by mental health phenotypes
#'
#' ord <- c(8 : 13, 1 : 7)
#' PPEI <- get_edgeprob(output, ord = ord)


########## get_edgeprob ##########

get_edgeprob <- function(output, ord = NULL)
{
  niter <- output$samplerPar$niter
  burnin <- output$samplerPar$burnin
  thin <- output$samplerPar$thin
  niter_TMP <- (niter - burnin) / thin
  q <- length(output$hyperPar$MrDAGcheck$Y_idx)
  p <- length(output$hyperPar$MrDAGcheck$X_idx)  
  colnames_YX <- colnames(output$hyperPar$U)
  
  if (is.null(dim(output$graphs)))
  {
    m <- sqrt(length(bd_decode(output$graphs[1])))
  } else {
    m <- dim(output$graphs)[1]
  }
  
  if (!is.null(output$hyperPar$MrDAGcheck))
  {
    group <- NULL
    group$Outcomes <- output$hyperPar$MrDAGcheck$Y_idx
    group$Exposures <- output$hyperPar$MrDAGcheck$X_idx
  } else {
    group <- NULL
  }
  
  if (!is.null(ord))
  {
    
    if (!is.null(output$hyperPar$MrDAGcheck))
    {
      Y_ord <- ord[ord <= q]
      X_ord <- ord[ord > q]
    } else {
      Y_ord <- ord
      X_ord <- NULL
    }
  
  } else {
    
    if (!is.null(output$hyperPar$MrDAGcheck))
    {
      Y_ord <- output$hyperPar$MrDAGcheck$Y_idx
      X_ord <- output$hyperPar$MrDAGcheck$X_idx
    } else {
      Y_ord <- seq(1, m, 1)
      X_ord <- NULL
    }
  
  }
  
  if (is.null(dim(output$graphs)))
  {
    graphs <- array(dim = c(m, m, niter_TMP), dimnames = list(colnames_YX, colnames_YX, 
                                                              paste0("iter", seq(burnin + thin, niter, thin))))
    
    for (iter in 1 : niter_TMP)
    {
      graphs[, , iter] <- bd_decode(output$graphs[iter])
    }
  
  } else {
    graphs <- output$graphs
  }
  
  edgeProb <- apply(graphs, c(1, 2), mean)
  diag(edgeProb) <- NA
  edgeProb[1 : q, q + (1 : p)] <- NA
  colnames(edgeProb) <- colnames_YX
  rownames(edgeProb) <- colnames_YX
  
  if (!is.null(output$hyperPar$MrDAGcheck))
  {
    edgeProb <- edgeProb[c(Y_ord, X_ord), c(Y_ord, X_ord)]
  } else {
    edgeProb <- edgeProb[Y_ord, Y_ord]
  }
  
  return(list(edgeProb = edgeProb, group = group, ord = ord))
}
