# 135 characters #####################################################################################################################
#' @title MrDAG data set: Lifestyle and behavioural exposures that might impact mental health phenotypes
#' @description The data set contains lifestyle and behavioural traits that are considered exposures of the risk of mental health 
#' phenotypes. As outcomes, seven mental health phenotypes are considered, including (in alphabetic order) attention deficit 
#' hyperactivity disorder (ADHD), anorexia nervosa (AN), autism spectrum disorder (ASD), bipolar disorder (BD), cognition (COG),  
#' major depressive disorder (MDD) and schizophrenia (SCZ). As exposures, six lifestyle and behavioural traits that have 
#' previously been investigated for their protective/risk effects on mental health are considered, including (in alphabetic order) 
#' alcohol consumption (ALC), education (in years) (EDU), leisure screen time (LST), physical activity (PA), lifetime smoking index 
#' (SM) and sleep duration (SP)
#'
#' @docType data
#'
#' @format A data frame consisting of 708 independent Instrumental Variables (IVs) selected to be associated at genome-wide 
#' significance with the exposures after pruning or clumping. For details, see \insertCite{Zuber2025;textual}{MrDAG}
#'
#' @references
#' \insertAllCited{}
#'
#'
#' @examples
# 100 characters ##################################################################################
#' # Example:
#'
#' data(LBT2MD_data)
#' head(LBT2MD_data)
"LBT2MD_data"
