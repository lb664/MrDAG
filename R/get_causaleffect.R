# 135 characters #####################################################################################################################
#' @title Estimation of the (average) causal effects under intervention on an exposure
#' @description Estimation of the (average) causal effect under intervention on a trait (target) and measured on another one 
#' (response) based on the Directed Acyclic Graphs (DAGs) explored by \code{MrDAG} algorithm
#'
#' @param output Output produced by \code{MrDAG} algorithm
#' @param response Trait (response) where the effect of the intervention on another trait (target) is measured
#' @param target Trait (target) under intervention
#' @param BMA If \code{TRUE}, Bayesian Model Averaging of the estimated (average) causal effect across all explored DAGs in the 
#' visited Completed Partially DAG (CPDAG) (\insertCite{Chickering2002;textual}{MrDAG}) or Essential Graph (EG) 
#' (\insertCite{Andersson1997;textual}{MrDAG}) is performed
#' @param CI Level (\code{0.95} default) of the credible interval (CI) of the (average) causal effect. It is calculated based on 
#' suitable quantiles of the estimated (average) causal effect across all explored DAGs
#'
#' @export
#'
#' @return The value returned is a list object \code{list(causaleffect, causaleffect_LL, causaleffect_UL, group, BMA, CI)}
#' \itemize{
#'   \item{\code{causaleffect}}{ Estimate of the (average) causal effect under intervention on an trait (target) on another 
#'         one (response) based on the DAGs explored by \code{MrDAG} algorithm. Feedback loop causal effect of the same trait is not 
#'         allowed (NA) }
#'   \item{\code{causaleffect_LL}}{ Lower limit of the \code{CI} of the (average) causal effect. It is calculated as the 
#'         (1-\code{CIboot})/2 % quantile of the (average) causal effect across all explored DAGs }
#'   \item{\code{causaleffect_UL}}{ Upper limit of the \code{CI} of the (average) causal effect. It is calculated as the 
#'         1-\[(1-\code{CIboot})/2\] % quantile of the (average) causal effect across all explored DAGs }
#'   \item{\code{group}}{ If \code{MrDAGcheck} in \code{MrDAG} algorithm contains the indices of the traits in \code{data} that 
#'         define the outcomes and the exposures, \code{group} coincides with this list }
#'   \item{\code{BMA}}{ Logical option }
#'   \item{\code{CI}}{ Level of the credible interval option } }
#'
#' @references
#' \insertAllCited{}
#'
#'
#' @examples
# 100 characters ##################################################################################
#' # Example: Estimation of the (average) causal effects under intervention on lifestyle and 
#' # behavioural exposures, and measured on mental health phenotypes. 708 independent Instrumental 
#' # Variables (IVs) were selected to be associated at genome-wide significance with six lifestyle 
#' # and behavioural traits after pruning or clumping which are considered exposures of the risk 
#' # of seven mental health outcome phenotypes
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
#' # Finally, the (average) causal effects and credible intervals (CI) of the intervention on ALC 
#' # and measured on SCZ are estimated
#'
#' causaleffects <- get_causaleffect(output, 7, 8)   # Intervention on ALC and measured on SCZ


########## get_causaleffect ##########

get_causaleffect <- function(output, response, target, BMA = TRUE, CI = 0.95)
{
  niter <- output$samplerpar$niter
  burnin <- output$samplerpar$burnin
  thin <- output$samplerpar$thin
  niter_TMP <- (niter - burnin) / thin
  q <- length(output$hyperpar$MrDAGcheck$Y_idx)
  p <- length(output$hyperpar$MrDAGcheck$X_idx)
  colnames_YX <- colnames(output$hyperpar$U) 
  
  if (is.null(dim(output$graphs)))
  {
    m <- sqrt(length(bd_decode(output$graphs[1])))
  } else {
    m <- dim(output$graphs)[1]
  }
  
  if (!is.null(output$hyperpar$MrDAGcheck))
  {
    group <- NULL
    group$Outcomes <- output$hyperpar$MrDAGcheck$Y_idx
    group$Exposures <- output$hyperpar$MrDAGcheck$X_idx
  } else {
    group <- NULL
  }
  
  if (is.null(dim(output$graphs)))
  {
    graphs <- array(dim = c(m, m, niter_TMP), dimnames = list(colnames_YX, colnames_YX, 
                                                              paste0("iter", seq(burnin + thin, niter, thin))))
    L <- array(dim = c(m, m, niter_TMP), dimnames = list(colnames_YX, colnames_YX, 
                                                         paste0("iter", seq(burnin + thin, niter, thin))))
    D <- array(dim = c(m, m, niter_TMP), dimnames = list(colnames_YX, colnames_YX, 
                                                         paste0("iter", seq(burnin + thin, niter, thin))))
    
    for (iter in 1 : niter_TMP)
    {
      graphs[, , iter] <- bd_decode(output$graphs[iter])
      L[, , iter] <- bd_decode(output$L[iter])
      D[, , iter] <- bd_decode(output$D[iter])
    }
  
  } else {
    graphs <- output$graphs
    L <- output$L
    D <- output$D
  }
  
  causaleffect <- matrix(0, ncol = length(target), nrow = niter_TMP)
  colnames(causaleffect) <- paste("Causal effect under intervention on", colnames_YX[target], "and measured on", colnames_YX[response], sep = " ")
  rownames(causaleffect) <- paste0("iter", seq(burnin + thin, niter, thin))
  
  for (iter in 1 : niter_TMP)
  {
    causaleffect[iter, ] <- causal_effect(response, target, L[, , iter], D[, , iter])
  }
  
  if (BMA == TRUE)
  {
    causaleffect_LL <- apply(causaleffect, 2, stats::quantile, (1 - CI) / 2)
    causaleffect_UL <- apply(causaleffect, 2, stats::quantile, 1 - ((1 - CI) / 2))
    causaleffect <- colMeans(causaleffect)
    
    return(list(causaleffect = causaleffect, causaleffect_LL = causaleffect_LL, causaleffect_UL = causaleffect_UL, group = group, BMA = BMA, CI = CI))
  } else {
    return(list(causaleffect = causaleffect, group = group, BMA = BMA, CI = CI))
  }
  
}
