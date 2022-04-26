#' @title EcM Fungal Data with Environmental Variables
#'
#' @description Ectomycorrhizal (EcM) fungal community and environmental data were excerpted from Gong et al. (2022). Sampling was conducted along the elevation gradient (2,900 m to 4,500 m) in the eastern slope of Baima Snow Mountain during both dry (November 2017) and wet (August 2018) seasons after studying community compositions of the host plants of EcM fungi.
#'
#' @docType data
#' @keywords datasets
#' @name baima
#' @usage data(baima.fun)
#' @usage data(baima.env)
#'
#' @format 263 samples from the root tips were excerpted. There are two linked data sets: \code{baima.fun}, a data frame containing 3,099 amplicon sequence variants (ASVs) of root associated EcM fungi; \code{baima.env}, a data frame containing 14 environmental variables.
#' @format The fields in the environmental data are:
#' \describe{
#' \item{environmental_medium}{Sample type}
#' \item{latitude}{Latitude}
#' \item{lontitude}{Lontitude}
#' \item{altitude}{Elevation (m)}
#' \item{season}{Season, a factor with levels \code{dry} and \code{wet}}
#' \item{em.GR}{Richness of EcM plant at genus level}
#' \item{em.abun}{The number of individuals of each EcM genus}
#' \item{sea.MT}{Dry-season and wet-season mean temperature}
#' \item{pH}{Soil pH}
#' \item{TP}{Total phosphorus (g/kg)}
#' \item{TK}{Total potassium (g/kg)}
#' \item{AN}{Alkaline-hydrolysable nitrogen (mg/kg)}
#' \item{AP}{Available phosphorus (mg/kg)}
#' \item{AK}{Available potassium (mg/kg)}
#' }
#'
#' @references Gong S, Feng B, Jian S P, et al. Elevation Matters More than Season in Shaping the Heterogeneity of Soil and Root Associated Ectomycorrhizal Fungal Community. Microbiology spectrum, 2022, 10(1): e01950-21.
#'
#' @examples help(baima)
NULL
