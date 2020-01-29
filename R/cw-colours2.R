#' set colours by category, define tol*qualitative vectors

library(RColorBrewer)
source("R/cw-palette.R")
## tol5qualitative= c("#332288","#88CCEE","#117733","#DDCC77","#CC6677")
## tol6qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677","#AA4499")
## tol7qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677","#AA4499")

## redblue<-colorRampPalette(c("red","white","blue"))
##     magblue<-colorRampPalette(c("magenta","white","blue"))
##     grnvi<-colorRampPalette(c("#00441b","gray90","#40004b"))
## grnmag<-colorRampPalette(c("#00441b","white","magenta"))

## brn <- "#543005"; grn <- "#003c30"
## brngrn <- colorRampPalette(c("#543005","white","#003c30"))

## # name=trait, entry=category
## ## https://www.r-bloggers.com/the-paul-tol-21-color-salute/
## tol5qualitative= c("#332288","#88CCEE","#117733","#DDCC77","#CC6677")
## tol6qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677","#AA4499")
## tol7qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677","#AA4499")
## study2col <- structure(tol7qualitative, names=unique(trait2study))

colored_labelled_boxes <- function (colors, labels, dend, rowLabels = NULL, cex.rowLabels = 0.9, 
    add = TRUE, y_scale, y_shift, text_shift = 1, sort_by_labels_order = TRUE, 
    horiz = FALSE, ...) {
    n_colors <- if (is.null(dim(colors))) 
        length(colors)
    else nrow(colors)
    n_groups <- if (is.null(dim(colors))) 
        1
    else ncol(colors)
    if (!missing(dend)) {
        if (is.hclust(dend)) 
            dend <- as.dendrogram(dend)
        if (!is.dendrogram(dend)) 
            stop("'dend' should be a dendrogram.")
        dend_labels <- labels(dend)
        dend_order <- order.dendrogram(dend)
    }
    else {
        dend_labels <- rep("W", n_colors)
        dend_order <- seq_len(n_colors)
    }
    if (!sort_by_labels_order) 
        dend_order <- seq_len(n_colors)
    if (!horiz) {
        if (missing(y_shift)) 
            y_shift <- -dendextend:::max_labels_height(dend_labels) + par("usr")[3L] - 
                strheight("X")
        if (missing(y_scale)) 
            y_scale <- strheight("X") * n_groups
    }
    else {
        if (missing(y_shift)) 
            y_shift <- -(min(strwidth(dend_labels)) + par("usr")[2L] + 
                strwidth("X")) * 1.2
        if (missing(y_scale)) 
            y_scale <- strwidth("X") * n_groups * 1.2
    }
    y_shift <- y_shift - y_scale
    colors <- as.matrix(colors)
    dimC <- dim(colors)
    if (is.null(rowLabels) & (length(dimnames(colors)[[2]]) == 
        dimC[2])) 
        rowLabels = names(as.data.frame(colors))
    op <- options()
    pr <- par(no.readonly = TRUE)
    options(stringsAsFactors = FALSE)
    par(xpd = TRUE)
    if (length(dend_order) != dimC[1]) 
        stop("ERROR: length of colors vector not compatible with number of objects in the hierarchical tree.")
    C <- colors[dend_order, ]
    C <- as.matrix(C)
    L <- as.matrix(labels[dend_order, ])
    step <- 1/(n_colors - 1)
    ystep <- 1/n_groups * 1.2
    if (!add) {
        barplot(height = 1, col = "white", border = FALSE, space = 0, 
            axes = FALSE, ...)
    }
    charWidth <- strwidth("W")/2
    charHeight <- strheight("W")/2
    for (j in 1:n_groups) {
        ind <- (1:n_colors)
        xl <- (ind - 1.5) * step
        xr <- (ind - 0.5) * step 
        yb <- rep(ystep * (j - 1), n_colors) 
        yt <- rep(ystep * j, n_colors) 
        w <- step/10
        if (add) {
            xl <- dendextend:::rescale(xl, to = c(1 - 0.5, n_colors - 0.5)) + 0.1
            xr <- dendextend:::rescale(xl, to = c(1 + 0.5, n_colors + 0.5)) - 0.1
            yb <- (yb + ystep/10) * y_scale + y_shift 
            yt <- (yt - ystep/10) * y_scale + y_shift 
        }
        if (horiz) {
            rect(-yb, xl, -yt, xr, col = as.character(C[, j]), 
                border = as.character(C[, j]))
            text(-(yb+yt)/2, (xl+xr)/2, labels=L[,j])
            par(srt = 90)
            if (is.null(rowLabels)) {
                s <- as.character(j)
                text(s, pos = 1, offset = 0.5, y = charHeight * 
                  text_shift - dendextend:::rotated_str_dim(s)[2]/2, x = -(ystep * 
                  (j) * y_scale + y_shift), cex = cex.rowLabels)
            }
            else {
                s <- rowLabels[j]
                text(s, pos = 1, offset = 0.5, y = charHeight * 
                  text_shift - dendextend:::rotated_str_dim(s)[2]/2, x = -(ystep * 
                  (j) * y_scale + y_shift), cex = cex.rowLabels)
            }
        }
        else {
            rect(xl, yb, xr, yt, col = as.character(C[, j]), 
                border = as.character(C[, j]))
            text((xl+xr)/2, (yb+yt)/2, labels=L[,j])
            if (is.null(rowLabels)) {
                text(as.character(j), pos = 2, x = charWidth * 
                  text_shift, y = ystep * (j - 0.5) * y_scale + 
                  y_shift, cex = cex.rowLabels)
            }
            else {
                text(rowLabels[j], pos = 2, x = charWidth * text_shift, 
                  y = ystep * (j - 0.5) * y_scale + y_shift, 
                  cex = cex.rowLabels)
            }
        }
    }
    options(op)
    par(pr)
    return(invisible(C))
}

colored_boxes <- function (colors, dend, rowLabels = NULL, cex.rowLabels = 0.9, 
    add = TRUE, y_scale, y_shift, text_shift = 1, sort_by_labels_order = TRUE, 
    horiz = FALSE, ...) {
    n_colors <- if (is.null(dim(colors))) 
        length(colors)
    else nrow(colors)
    n_groups <- if (is.null(dim(colors))) 
        1
    else ncol(colors)
    if (!missing(dend)) {
        if (is.hclust(dend)) 
            dend <- as.dendrogram(dend)
        if (!is.dendrogram(dend)) 
            stop("'dend' should be a dendrogram.")
        dend_labels <- labels(dend)
        dend_order <- order.dendrogram(dend)
    }
    else {
        dend_labels <- rep("W", n_colors)
        dend_order <- seq_len(n_colors)
    }
    if (!sort_by_labels_order) 
        dend_order <- seq_len(n_colors)
    if (!horiz) {
        if (missing(y_shift)) 
            y_shift <- -dendextend:::max_labels_height(dend_labels) + par("usr")[3L] - 
                strheight("X")
        if (missing(y_scale)) 
            y_scale <- strheight("X") * n_groups
    }
    else {
        if (missing(y_shift)) 
            y_shift <- -(min(strwidth(dend_labels)) + par("usr")[2L] + 
                strwidth("X"))
        if (missing(y_scale)) 
            y_scale <- strwidth("X") * n_groups
    }
    y_shift <- y_shift - y_scale
    colors <- as.matrix(colors)
    dimC <- dim(colors)
    if (is.null(rowLabels) & (length(dimnames(colors)[[2]]) == 
        dimC[2])) 
        rowLabels = names(as.data.frame(colors))
    op <- options()
    pr <- par(no.readonly = TRUE)
    options(stringsAsFactors = FALSE)
    par(xpd = TRUE)
    if (length(dend_order) != dimC[1]) 
        stop("ERROR: length of colors vector not compatible with number of objects in the hierarchical tree.")
    C <- colors[dend_order, ]
    C <- as.matrix(C)
    step <- 1/(n_colors - 1)
    ystep <- 1/n_groups
    if (!add) {
        barplot(height = 1, col = "white", border = FALSE, space = 0, 
            axes = FALSE, ...)
    }
    charWidth <- strwidth("W")/2
    charHeight <- strheight("W")/2
    for (j in 1:n_groups) {
        ind <- (1:n_colors)
        xl <- (ind - 1.5) * step
        xr <- (ind - 0.5) * step 
        yb <- rep(ystep * (j - 1), n_colors) 
        yt <- rep(ystep * j, n_colors) 
        w <- step/10
        if (add) {
            xl <- dendextend:::rescale(xl, to = c(1 - 0.5, n_colors - 0.5)) + 0.1
            xr <- dendextend:::rescale(xl, to = c(1 + 0.5, n_colors + 0.5)) - 0.1
            yb <- (yb + ystep/10) * y_scale + y_shift 
            yt <- (yt - ystep/10) * y_scale + y_shift 
        }
        if (horiz) {
            rect(-yb, xl, -yt, xr, col = as.character(C[, j]), 
                border = as.character(C[, j]))
            par(srt = 90)
            if (is.null(rowLabels)) {
                s <- as.character(j)
                text(s, pos = 1, offset = 0.5, y = charHeight * 
                  text_shift - dendextend:::rotated_str_dim(s)[2]/2, x = -(ystep * 
                  (j) * y_scale + y_shift), cex = cex.rowLabels)
            }
            else {
                s <- rowLabels[j]
                text(s, pos = 1, offset = 0.5, y = charHeight * 
                  text_shift - dendextend:::rotated_str_dim(s)[2]/2, x = -(ystep * 
                  (j) * y_scale + y_shift), cex = cex.rowLabels)
            }
        }
        else {
            rect(xl, yb, xr, yt, col = as.character(C[, j]), 
                border = as.character(C[, j]))
            if (is.null(rowLabels)) {
                text(as.character(j), pos = 2, x = charWidth * 
                  text_shift, y = ystep * (j - 0.5) * y_scale + 
                  y_shift, cex = cex.rowLabels)
            }
            else {
                text(rowLabels[j], pos = 2, x = charWidth * text_shift, 
                  y = ystep * (j - 0.5) * y_scale + y_shift, 
                  cex = cex.rowLabels)
            }
        }
    }
    options(op)
    par(pr)
    return(invisible(C))
}

