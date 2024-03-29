% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/barplot_hp.R
\name{barplot_hp}
\alias{barplot_hp}
\title{Visualization of HP Using Column Diagram}
\usage{
barplot_hp(
  x,
  order.var = TRUE,
  decreasing.var = TRUE,
  cutoff = -1,
  col.fill = "valid",
  col.color = NULL,
  col.width = 0.6,
  show.effect = TRUE,
  effect.cex = 2.7,
  title.cex = 10,
  axis.cex = 8
)
}
\arguments{
\item{x}{A \code{\link{rdacca.hp}} object, which contains the output of HP from canonical analysis.}

\item{order.var}{The predictors in the matrix layout should be ordered by. Default is \code{TRUE}, which orders the predictors by their effect values. IF \code{FALSE}, sort by the order of predictors in input data.}

\item{decreasing.var}{If \code{order.var=TRUE}, how the predictors should be ordered. Default is \code{TRUE}, from greatest to least.}

\item{cutoff}{Effect values below \code{cutoff} will not be displayed, default is \code{-1}. Note: Negative values due to adjustment of R-squared mean negligible contributions, but they are included in the computation of the total contribution of each predictor category.}

\item{col.fill}{How the bars should be colored. Options include \code{"valid"} (according to the validity of effects) or \code{"vars"} (color by predictors), default is \code{"valid"}.}

\item{col.color}{Color of bars.}

\item{col.width}{Width of bars, default is \code{0.6}.}

\item{show.effect}{Show the effect values above bars, default is \code{TRUE}.}

\item{effect.cex}{Font size of the effect values, default is \code{2.7}.}

\item{title.cex}{Font size of axis titles, default is \code{10}.}

\item{axis.cex}{Font size of axis labels, default is \code{8}.}
}
\value{
\itemize{Returns a ggplot2.}
}
\description{
Visualization of individual effects in hierarchical partitioning (HP) using column diagram.
}
\details{
This function is used to visualize the object of \code{\link{rdacca.hp}} (Lai et al. 2022), which calculates the individual effects of predictor variables or groups of predictor variables in canonical analysis based on HP.
}
\examples{
## A simple example of partial dbRDA
data(baima.fun)
data(baima.env)

# Bray-Curtis index was used to calculate community composition dissimilarity
baima.fun.bray <- vegan::vegdist(baima.fun, method = "bray")

# Quantify the individual effects of soil properties on EcM fungal community composition
soil <- baima.env[c("pH", "TP", "TK", "AN", "AP", "AK")]
baima.soil.vp <- rdacca.hp::rdacca.hp(baima.fun.bray, soil, method = "dbRDA", type = "adjR2")

# Plot individual effects
barplot_hp(baima.soil.vp, col.fill = "var", 
 col.color = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69"))
}
\references{
Lai J., Zou Y., Zhang J., Peres-Neto P. (2022) Generalizing hierarchical and variation partitioning in multiple regression and canonical analyses using the rdacca.hp R package. Methods in Ecology and Evolution, 13:782-788.
}
