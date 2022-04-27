# UpSetVP

<!-- badges: start -->
<!-- badges: end -->

Visualization of variance decomposition analysis (VPA) and hierarchical partitioning (HP) with unlimited number of predictor variables (or matrices of predictors) using UpSet matrix layout (Conway et al. 2017; Lex et al. 2014).

## Installation

Install the released version of `UpSetVP` from
[GitHub](https://github.com/LiuXYh/UpSetVP) with:

``` r
# install.packages('devtools')
devtools::install_github('LiuXYh/UpSetVP', force = TRUE)
```

## A Simple Example

Load packages.

``` r
library(rdacca.hp)
library(ggplot2)
library(patchwork)
library(UpSetVP)
```

Ectomycorrhizal (EcM) fungal community and environmental data were excerpted from Gong et al. (2022).

``` r
data(baima.fun)
data(baima.env)
```

Quantify the relative importance of individual soil properties (pH, TP, TK, AN, AP, AK) on the composition of EcM fungal community by using partial dbRDA.

``` r
# Bray-Curtis index was used to calculate community composition dissimilarity
baima.fun.bray <- vegdist(baima.fun, method = 'bray')

# VPA and HP by using rdacca.hp package (Lai et al. 2022)
soil <- baima.env[c('pH', 'TP', 'TK', 'AN', 'AP', 'AK')]
baima.soil.vp <- rdacca.hp(baima.fun.bray, soil, method = 'dbRDA', var.part = TRUE, type = 'adjR2')

# Plot unique, common, as well as individual effects
upset_vp(baima.soil.vp, plot.hp = TRUE)
```

<img src="man/figures/1.png" height="80%" width="80%" />

``` r
# Only plot individual effects
barplot_hp(baima.soil.vp, col.fill = 'var')
```

<img src="man/figures/2.png" height="40%" width="40%" />

##
The relative importance of groups of environmental factors on EcM fungal community composition.<br>
Environmental factors including elevation, season, space (dbMEM1 and dbMEM2), host (em.GR and em.abun), climate (sea.MT), and soil (pH, TP, TK, AN, AP, and AK).

``` r
# Distance-based Moran's eigenvector maps (dbMEM) was used to extract spatial relationships
space.dbmem <- adespatial::dbmem(baima.env[c('latitude', 'lontitude')])

# VPA and HP by using rdacca.hp package
env.list <- list(
    elevation = baima.env['altitude'],
    season = baima.env['season'],
    space = data.frame(space.dbmem)[1:2],
    host = baima.env[c('em.GR', 'em.abun')],
    climate = baima.env['sea.MT'],
    soil = baima.env[c('pH', 'TP', 'TK', 'AN', 'AP', 'AK')]
)
baima.env.vp <- rdacca.hp(baima.fun.bray, env.list, method = 'dbRDA', var.part = TRUE, type = 'adjR2')

# Plot unique, common, as well as individual effects
upset_vp(baima.env.vp, plot.hp = TRUE, order.part = 'degree')
```

<img src="man/figures/3.png" height="80%" width="80%" />

``` r
# Only plot individual effects
barplot_hp(baima.env.vp, col.fill = 'var', col.color = c('#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69'))
```

<img src="man/figures/4.png" height="40%" width="40%" />

## References

Conway J R, Lex A, Gehlenborg N. UpSetR: an R package for the visualization of intersecting sets and their properties. Bioinformatics, 2017, 33(18): 2938-2940.<br>
Gong S, Feng B, Jian S P, et al. Elevation Matters More than Season in Shaping the Heterogeneity of Soil and Root Associated Ectomycorrhizal Fungal Community. Microbiology spectrum, 2022, 10(1): e01950-21.<br>
Lai J, Zou Y, Zhang J, et al. Generalizing hierarchical and variation partitioning in multiple regression and canonical analyses using the rdacca. hp R package. Methods in Ecology and Evolution, 2022.<br>
Lex A, Gehlenborg N, Strobelt H, et al. UpSet: visualization of intersecting sets. IEEE transactions on visualization and computer graphics, 2014, 20(12): 1983-1992.
