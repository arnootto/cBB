#' Mule Deer Survival Data
#'
#' A dataset containing survival and mortality data for mule deer fawns in Colorado, Idaho, and Montana from 1981 to 1996.
#'
#' @format A data frame with 26 rows and 10 columns:
#' \describe{
#'   \item{Year}{The year of observation, ranging from 1981 to 1996.}
#'   \item{State}{The state where the data was collected: Colorado, Idaho, or Montana.}
#'   \item{Radiocollared_fawns}{Number of radio-collared fawns monitored.}
#'    \item{Mortalities}{Number of fawns that died during the observation period.}
#'   \item{Predation_n}{Number of mortalities attributed to predation.}
#'   \item{Predation_percent}{Percentage of mortalities due to predation.}
#'   \item{Winter_malnutrition_n}{Number of mortalities attributed to winter malnutrition.}
#'   \item{Winter_malnutrition_percent}{Percentage of mortalities due to winter malnutrition.}
#'   \item{Other_mortality_n}{Number of mortalities due to other causes.}
#'   \item{Other_mortality_percent}{Percentage of mortalities due to other causes.}
#' }
#' @source Unsworth, J. W., Pac, D. F., White, G. C., & Bartmann, R. M. (1999). Mule deer survival in Colorado, Idaho, and Montana. The Journal of Wildlife Management, 63(1), 315-326.
#' @examples
#' data(MuleDeer)
#' summary(MuleDeer)
#' @docType data
#' @name MuleDeer
#' @usage data(MuleDeer)
#' @keywords datasets
#' @format A data frame with 26 rows and 10 variables
NULL
