# 135 characters #####################################################################################################################
#' @title Estimation of the (average) causal effects under intervention on the exposures
#' @description Estimation of the (average) causal effects under intervention on the exposures (targets) and measured on the outcomes 
#' (responses) based on the Directed Acyclic Graphs (DAGs) explored by \code{MrDAG} algorithm
#'
#' @param output Output produced by \code{MrDAG} algorithm
#' @param ord Indices of the traits in \code{data} that define the outcomes and the exposures. They specify the order of appearance 
#' of the traits when printing the (average) causal effect. If \code{ord = NULL}, the order is the same as \code{MrDAGcheck} in 
#' \code{MrDAG} ouput
#' @param BMA If \code{TRUE}, Bayesian Model Averaging of the estimated (average) causal effect across all explored DAGs in the 
#' visited Completed Partially DAG (CPDAG) (\insertCite{Chickering2002;textual}{MrDAG}) or Essential Graph (EG) 
#' (\insertCite{Andersson1997;textual}{MrDAG}) is performed
#' @param CI Level (\code{0.95} default) of the credible interval (CI) of the (average) causal effect. It is calculated based on 
#' suitable quantiles of the estimated (average) causal effect across all explored DAGs
#' @param progress Logical value set as \code{FALSE} by default to print on the screen the progress of the causal effects estimation
#' @param plot Logical value (default \code{FALSE}). If \code{TRUE}, \code{get_causaleffects} generates the output to be used in a 
#' plot
#'
#' @export
#'
#' @return The value returned is a list object \code{list(causalEffects, causalEffects_LL, causalEffects_UL, group, ord, BMA, CI)}
#' \itemize{
#'   \item{\code{causalEffects}}{ Estimate of the (average) causal effects under intervention on the exposures based on the DAGs 
#'         explored by \code{MrDAG} algorithm. Feedback loop causal effect of the same outcome is not allowed (NA) }
#'   \item{\code{causalEffects_LL}}{ Lower limit of the \code{CI} of the (average) causal effects. It is calculated as the 
#'         (1-\code{CIboot})/2 % quantile of the (average) causal effect across all explored DAGs }
#'   \item{\code{causalEffects_UL}}{ Upper limit of the \code{CI} of the (average) causal effects. It is calculated as the 
#'         1-\[(1-\code{CIboot})/2\] % quantile of the (average) causal effect across all explored DAGs }
#'   \item{\code{group}}{ If \code{MrDAGcheck} in \code{MrDAG} algorithm contains the indices of the traits in \code{data} that 
#'         define the outcomes and the exposures, \code{group} coincides with this list }
#'   \item{\code{ord}}{ Indices of the traits in \code{data} that specify the outcomes and the exposures. It might differ from 
#'         \code{group} if \code{ord} has been specified and it is different from \code{MrDAGcheck} }
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
#' # Finally, the posterior (average) causal effects and 90% credible intervals are estimated
#'
#' ord <- c(8 : 13, 1 : 7)
#' causalEffects <- get_causaleffects(output, ord = ord, CI = 0.90, progress = TRUE)


########## get_causaleffects ##########

get_causaleffects <- function(output, ord = NULL, BMA = TRUE, CI = 0.95, progress = FALSE, plot = FALSE)
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
  
  if (!is.null(output$hyperPar$MrDAGcheck))
  {
    dimnames <- list(colnames_YX, colnames_YX[output$hyperPar$MrDAGcheck$Y_idx])
    
    if (BMA == TRUE)
    {
      causalEffects <- matrix(NA, nrow = m, ncol = q, dimnames = dimnames)
      causalEffects_LL <- matrix(NA, nrow = m, ncol = q, dimnames = dimnames)
      causalEffects_UL <- matrix(NA, nrow = m, ncol = q, dimnames = dimnames)
    } else {
      causalEffects <- array(NA, dim = c((niter - burnin) / thin, m, q), dimnames = list(paste0("iter", seq(burnin + thin, niter, thin)), 
                                                                                         colnames_YX, colnames_YX[output$hyperPar$MrDAGcheck$Y_idx])) 
    }
    
    for (k in 1 : q)
    {
      response_idx <- k
      target_idx <- (1 : m)[-k]
      causalEffect_TMP <- get_causaleffect(output, response = response_idx, target = target_idx, BMA = BMA, CI = CI)
      
      if (BMA == TRUE)
      {
        causalEffects[, k] <- append(causalEffect_TMP$causalEffect, NA, after = (k - 1))
        causalEffects_LL[, k] <- append(causalEffect_TMP$causalEffect_LL, NA, after = (k - 1))
        causalEffects_UL[, k] <- append(causalEffect_TMP$causalEffect_UL, NA, after = (k - 1))
      } else {
        causalEffects[, target_idx, k] <- causalEffect_TMP$causalEffect
      }
      
      if (progress == TRUE)
      {
        cat("Estimation (average) causal effect under intervention on exposure X and measured on outcome", colnames_YX[k], "\n")
      } else {
        
        if (k == 1)
        {
          cat("Estimation (average) causal effects under intervention on the exposures\n")
        }
      
      }
    
    }
    
    if (plot == TRUE)
    {
      
      if (BMA == TRUE)
      {
        causalEffects <- causalEffects[X_ord, Y_ord]
        causalEffects_LL <- causalEffects_LL[X_ord, Y_ord]
        causalEffects_UL <- causalEffects_UL[X_ord, Y_ord]
      } else {
        causalEffects <- causalEffects[, X_ord, Y_ord]
      }
    
    } else {
      
      if (BMA == TRUE)
      {
        causalEffects <- causalEffects[X_ord, Y_ord]
        causalEffects_LL <- causalEffects_LL[X_ord, Y_ord]
        causalEffects_UL <- causalEffects_UL[X_ord, Y_ord]
      } else {
        causalEffects <- causalEffects[, X_ord, Y_ord]
      }
    
    }
  
  } else {
    dimnames <- list(colnames_YX, colnames_YX)
    
    if (BMA == TRUE)
    {
      causalEffects <- matrix(NA, nrow = m, ncol = m, dimnames = dimnames)
      causalEffects_LL <- matrix(NA, nrow = m, ncol = m, dimnames = dimnames)
      causalEffects_UL <- matrix(NA, nrow = m, ncol = m, dimnames = dimnames)
    } else {
      causalEffects <- array(NA, dim = c((niter - burnin) / thin, m, m), dimnames = list(paste0("iter", seq(burnin + thin, niter, thin)), 
                                                                                        colnames_YX, colnames_YX)) 
    }
    
    for (k in 1 : m)
    {
      response_idx <- k
      target_idx <- (1 : m)[-k]
      causalEffect_TMP <- get_causaleffect(output, response = response_idx, target = target_idx, BMA = BMA, CI = CI)
      
      if (BMA == TRUE)
      {
        causalEffects[, k] <- append(causalEffect_TMP$causalEffect, NA, after = (k - 1))
        causalEffects_LL[, k] <- append(causalEffect_TMP$causalEffect_LL, NA, after = (k - 1))
        causalEffects_UL[, k] <- append(causalEffect_TMP$causalEffect_UL, NA, after = (k - 1))
      } else {
        causalEffects[, target_idx, k] <- causalEffect_TMP$causalEffect
      }
      
      if (progress == TRUE)
      {
        cat("Estimation (average) causal effects under intervention on a trait and measured on trait", colnames_YX[k], "\n")
      } else {
        
        if (k == 1)
        {
          cat("Estimation (average) causal effects under intervention on the traits\n")
        }
      
      }
    
    }
    
    if (BMA == TRUE)
    {
      causalEffects <- causalEffects[Y_ord, Y_ord]
      causalEffects_LL <- causalEffects_LL[Y_ord, Y_ord]
      causalEffects_UL <- causalEffects_UL[Y_ord, Y_ord] 
    } else {
      causalEffects <- causalEffects[, Y_ord, Y_ord]
    }
  
  }
  
  if (BMA == TRUE)
  {
    return(list(causalEffects = causalEffects, causalEffects_LL = causalEffects_LL, causalEffects_UL = causalEffects_UL, group = group, ord = ord, BMA = BMA, CI = CI))
  } else {
    return(list(causalEffects = causalEffects, group = group, ord = ord, BMA = BMA, CI = CI))
  }

}
