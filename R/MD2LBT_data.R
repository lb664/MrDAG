# 135 characters #####################################################################################################################
#' @title MrDAG data set: Mental health phenotypes that might impact lifestyle and behavioural traits
#' @description The data set contains mental health phenotypes that are considered exposures of the risk of lifestyle and behavioural 
#' traits. As outcomes, six lifestyle and behavioural traits are considered, including (in alphabetic order) alcohol consumption 
#' (ALC), education (in years) (EDU), leisure screen time (LST), physical activity (PA), lifetime smoking index (SM) and sleep 
#' duration (SP). As exposures, seven mental health phenotypes are considered, including (in alphabetic order) attention deficit 
#' hyperactivity disorder (ADHD), anorexia nervosa (AN), autism spectrum disorder (ASD), bipolar disorder (BD), cognition (COG), 
#' major depressive disorder (MDD) and schizophrenia (SCZ)
#'
#' @docType data
#'
#' @format A data frame consisting of 470 independent Instrumental Variables (IVs) selected to be associated at genome-wide 
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
#' data(MD2LBT_data)
#' head(MD2LBT_data)
"MD2LBT_data"
