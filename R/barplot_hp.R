#' @title Visualization of HP Using Column Diagram
#'
#' @description Visualization of individual effects in hierarchical partitioning (HP) using column diagram.
#'
#' @param rdacca A \code{\link{rdacca.hp}} object, which contains the output of HP from canonical analysis.
#' @param order.var The predictors in the matrix layout should be ordered by. Default is \code{TRUE}, which orders the predictors by their effect values. IF \code{FALSE}, sort by the order of predictors in input data.
#' @param decreasing.var If \code{order.var=TRUE}, how the predictors should be ordered. Default is \code{TRUE}, from greatest to least.
#' @param cutoff Effect values below \code{cutoff} will not be displayed, default is \code{-1}. Note: Negative values due to adjustment of R-squared mean negligible contributions, but they are included in the computation of the total contribution of each predictor category.
#' @param col.fill How the bars should be colored. Options include \code{"valid"} (according to the validity of effects) or \code{"vars"} (color by predictors), default is \code{"valid"}.
#' @param col.color Color of bars.
#' @param col.width Width of bars, default is \code{0.6}.
#' @param show.effect Show the effect values above bars, default is \code{TRUE}.
#' @param effect.cex Font size of the effect values, default is \code{2.7}.
#' @param title.cex Font size of axis titles, default is \code{10}.
#' @param axis.cex Font size of axis labels, default is \code{8}.
#'
#' @details This function is used to visualize the object of \code{\link{rdacca.hp}} (Lai et al. 2022), which calculates the individual effects of predictor variables or groups of predictor variables in canonical analysis based on HP.
#'
#' @return \itemize{Returns a ggplot2.}
#'
#' @references Lai J, Zou Y, Zhang J, et al. Generalizing hierarchical and variation partitioning in multiple regression and canonical analyses using the rdacca.hp R package. Methods in Ecology and Evolution, 2022.
#'
#' @export
#' @examples
#' \dontrun{
#' require(rdacca.hp)
#' require(ggplot2)
#'
#' ## A simple example of partial dbRDA
#' data(baima.fun)
#' data(baima.env)
#'
#' # Bray-Curtis index was used to calculate community composition dissimilarity
#' baima.fun.bray <- vegdist(baima.fun, method = "bray")
#' 
#' # Quantify the individual effects of soil properties on EcM fungal community composition
#' soil <- baima.env[c("pH", "TP", "TK", "AN", "AP", "AK")]
#' baima.soil.vp <- rdacca.hp(baima.fun.bray, soil, method = "dbRDA", type = "adjR2")
#' 
#' # Plot individual effects
#' barplot_hp(baima.soil.vp, col.fill = "var", 
#'    col.color = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69"))
#' }

barplot_hp <- function(rdacca, order.var = TRUE, decreasing.var = TRUE, cutoff = -1, col.fill = "valid", col.color = NULL, col.width = 0.6, show.effect = TRUE, effect.cex = 2.7, title.cex = 10, axis.cex = 8) {
    Constrained <- 100*rdacca$Total_explained_variation
    Hier.part <- as.data.frame(rdacca$Hier.part)
    Hier.part$Var <- rownames(Hier.part)
    Hier.part <- Hier.part[which(Hier.part$Individual >= cutoff), ]
    Hier.part$Individual <- 100*Hier.part$Individual
    if (order.var)
        Hier.part <- Hier.part[order(Hier.part$Individual, decreasing = !decreasing.var), ]
    Hier.part$Var <- factor(Hier.part$Var, levels = Hier.part$Var)
    if (col.fill == "valid")
        Hier.part$valid <- apply(Hier.part[3], 1, function(x) ifelse(x <= 0, "0", "1"))
    
    p.hp <- ggplot(Hier.part, aes_string(x = "Var", y = "Individual")) +
        theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.line.x = element_line(color = "black"),
            axis.text = element_text(color = "black", size = axis.cex),
            axis.ticks = element_line(color = "black"),
            axis.ticks.y = element_blank(),
            axis.title = element_text(color = "black", size = title.cex),
            plot.title = element_text(hjust = 0.5, size = title.cex),
            legend.position = "none") +
            coord_flip() +
            scale_y_continuous(expand = expansion(mult = c(ifelse(min(Hier.part$Individual) < 0, 0.3, 0), 0.3))) +
            labs(y = "Individual (%)", x = "", title = paste("Constrained:", Constrained, "%", "  ", "Residual", 100-Constrained, "%"))
    
    if (col.fill == "valid") {
        if (is.null(col.color))
            col.color <- colorRampPalette(c("#EEEEEE", "#000000"))(2)
        p.hp <- p.hp + 
            geom_col(aes_string(fill = "valid"), width = col.width) +
            scale_fill_manual(values = col.color, limits = c("0", "1"))
    }
    else if (col.fill == "var") {
        if (is.null(col.color))
            col.color <- colorRampPalette(c("#EEEEEE", "#000000"))(length(levels(Hier.part$Var)))
        p.hp <- p.hp +
            geom_col(aes_string(fill = "Var"), width = col.width) +
            scale_fill_manual(values = col.color)
    }
    
    if (show.effect)
        p.hp$layers[[2]] <- geom_text(aes_string(label = "Individual", hjust = ifelse(Hier.part$Individual < 0, 1.2, -0.2)), color = "black", size = effect.cex)
    p.hp <- p.hp + geom_hline(yintercept = 0)
    return(p.hp)
}
