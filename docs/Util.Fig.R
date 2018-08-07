require(grid)
require(gridExtra)
require(gtable)
require(ggplot2)

match.widths.grob <- function(g.list) {

    max.width <- g.list[[1]]$widths[2:7]

    for(j in 2:length(g.list)) {
        max.width <- grid::unit.pmax(max.width, g.list[[j]]$widths[2:7])
    }

    for(j in 1:length(g.list)) {
        g.list[[j]]$widths[2:7] <- as.list(max.width)
    }
    return(g.list)
}

match.widths <- function(p.list) {
    g.list <- lapply(p.list, ggplotGrob)
    return(match.widths.grob(g.list))
}

grid.vcat <- function(p.list, ...) {
    g.list <- match.widths(p.list)
    ret <- grid.arrange(grobs = g.list, ncol = 1,  ...)
    return(ret)
}

match.heights.grob <- function(g.list, stretch = TRUE)  {
    max.height <- g.list[[1]]$heights[2:7]

    if(stretch) {
        for(j in 2:length(g.list)) {
            max.height <- grid::unit.pmax(max.height, g.list[[j]]$heights[2:7])
        }
    }

    for(j in 1:length(g.list)) {
        g.list[[j]]$heights[2:7] <- as.list(max.height)
    }

    return(g.list)
}

match.heights <- function(p.list, stretch = FALSE) {
    g.list <- lapply(p.list, ggplotGrob)
    return(match.heights.grob(g.list, stretch))
}

grid.hcat <- function(p.list, ...) {
    g.list <- match.heights(p.list, stretch = TRUE)
    ret <- grid.arrange(grobs = g.list, nrow = 1,  ...)
    return(ret)
}

################################################################
match.2by2 <- function(p1, p2, p3, p4) {

    g1 <- ggplotGrob(p1)
    g2 <- ggplotGrob(p2)
    g3 <- ggplotGrob(p3)
    g4 <- ggplotGrob(p4)

    gg <- match.widths.grob(list(g1, g3))
    g1 <- gg[[1]]
    g3 <- gg[[2]]

    gg <- match.widths.grob(list(g2, g4))
    g2 <- gg[[1]]
    g4 <- gg[[2]]

    gg <- match.heights.grob(list(g1, g2))
    g1 <- gg[[1]]
    g2 <- gg[[2]]

    gg <- match.heights.grob(list(g3, g4))
    g3 <- gg[[1]]
    g4 <- gg[[2]]
    
    list(g1, g2, g3, g4)
}

################################################################

row.order <- function(mat) {
    require(cba)
    require(proxy)

    if(nrow(mat) < 3) {
        return(1:nrow(mat))
    }

    D <- proxy::dist(mat, method = function(a,b) 1 - cor(a,b, method = 'spearman'))
    D[!is.finite(D)] <- 0
    h.out <- hclust(D)
    o.out <- cba::order.optimal(D, h.out$merge)
    return(o.out$order)
}

gg.plot <- function(...) {
    require(ggplot2)
    ggplot(...) + theme_bw() +
        theme(plot.background = element_blank(),
              panel.background = element_blank(),
              strip.background = element_blank(),
              legend.background = element_blank(),
              axis.text = element_text(size = 8),
              axis.title = element_text(size = 8),
              legend.key.width = unit(.5, 'lines'),
              legend.title = element_text(size = 6),
              legend.text = element_text(size = 6))    
}

order.pair <- function(pair.tab) {

    require(tidyr)
    require(dplyr)
    M <- pair.tab %>% tidyr::spread(key = col, value = weight, fill = 0)
    rr <- M[, 1] %>% unlist()
    cc <- colnames(M)[-1] %>% unlist()

    log.msg('Built the Mat: %d x %d', nrow(M), ncol(M))
    ro <- row.order(M %>% dplyr::select(-row) %>% as.matrix())
    log.msg('Sort the rows: %d', length(ro))

    co <- order(apply(M[ro, -1], 2, which.max), decreasing = TRUE)

    ## co <- row.order(t(M %>% dplyr::select(-row) %>% as.matrix()))
    ## log.msg('Sort the columns: %d', length(co))

    list(rows = rr[ro], cols = cc[co], M = M)
}

