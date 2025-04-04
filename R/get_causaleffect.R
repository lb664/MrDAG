# 135 characters #####################################################################################################################
#' @title Estimation of the (average) marginal causal effects under intervention on an exposure
#' @description Estimation of the (average) marginal causal effect under intervention on a trait (target) and measured on another one 
#' (response) based on the Directed Acyclic Graphs (DAGs) explored by \code{MrDAG} algorithm
#'
#' @param output Output produced by \code{MrDAG} algorithm
#' @param response Trait (response) where the effect of the intervention on another trait (target) is measured. The index refers to
#' position of the trait (response) in the output produced by \code{MrDAG} algorithm
#' @param target Trait (target) under intervention. The index refers to position of the trait (response) in the output produced by 
#' \code{MrDAG} algorithm
#' @param BMA If \code{TRUE}, Bayesian Model Averaging of the estimated (average) causal effect across all explored DAGs in the 
#' visited Completed Partially DAG (CPDAG) (\insertCite{Chickering2002;textual}{MrDAG}) or Essential Graph (EG) 
#' (\insertCite{Andersson1997;textual}{MrDAG}) is performed
#' @param CI Level (\code{0.95} default) of the credible interval (CI) of the (average) causal effect. It is calculated based on 
#' suitable quantiles of the estimated (average) causal effect across all explored DAGs
#'
#' @export
#'
#' @return The value returned is a list object \code{list(causalEffect, causalEffect_LL, causalEffect_UL, group, BMA, CI)}
#' \itemize{
#'   \item{\code{causalEffect}}{ Estimate of the (average) causal effect under intervention on an trait (target) on another 
#'         one (response) based on the DAGs explored by \code{MrDAG} algorithm. Feedback loop causal effect of the same trait is not 
#'         allowed (NA). Note that the estimated causal effect corresponds to the marginal causal effect of a target on a response 
#'         *without* considering the effects of the others traits }
#'   \item{\code{causalEffect_LL}}{ Lower limit of the \code{CI} of the (average) causal effect. It is calculated as the 
#'         (1-\code{CIboot})/2 % quantile of the (average) causal effect across all explored DAGs }
#'   \item{\code{causalEffect_UL}}{ Upper limit of the \code{CI} of the (average) causal effect. It is calculated as the 
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
#'                 niter = 7500, burnin = 2500, thin = 5, tempMax = 20, pp = 0.01, 
#'                 MrDAGcheck = MrDAGcheck, fileName = NULL)
#'
#' # Finally, the posterior (average) causal effects and credible intervals (CI) of the intervention
#' # on ALC and measured on SCZ are estimated
#'
#' causalEffect <- get_causaleffect(output, 7, 8)   # Intervention on ALC(8) and measured on SCZ(7)
#'
#' causalEffect <- get_causaleffect(output, 7, c(4, 8))   # Intervention on ALC(8) and measured on 
#'                                                        # SCZ(7) after conditioning on BD(4)


########## get_causaleffect ##########

get_causaleffect <- function(output, response, target, BMA = TRUE, CI = 0.95)
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
  
  causalEffect <- matrix(0, ncol = length(target), nrow = niter_TMP)
  colnames(causalEffect) <- paste("Causal effect under intervention on", colnames_YX[target], "and measured on", colnames_YX[response], sep = " ")
  rownames(causalEffect) <- paste0("iter", seq(burnin + thin, niter, thin))
  
  for (iter in 1 : niter_TMP)
  {
    causalEffect[iter, ] <- causal_effect(response, target, L[, , iter], D[, , iter])
  }
  
  if (BMA == TRUE)
  {
    causalEffect_LL <- apply(causalEffect, 2, stats::quantile, (1 - CI) / 2)
    causalEffect_UL <- apply(causalEffect, 2, stats::quantile, 1 - ((1 - CI) / 2))
    causalEffect <- colMeans(causalEffect)
    
    return(list(causalEffect = causalEffect, causalEffect_LL = causalEffect_LL, causalEffect_UL = causalEffect_UL, group = group, BMA = BMA, CI = CI))
  } else {
    return(list(causalEffect = causalEffect, group = group, BMA = BMA, CI = CI))
  }
  
}
