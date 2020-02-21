#!/bin/R
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

###################################################################################################
## Load required packages
###################################################################################################
libs <- paste0(Sys.getenv("CONDA_PREFIX"),"/lib/R/library")
.libPaths(libs)
.libPaths()
library(RColorBrewer)
library(igraph)
library(ape)
library(phytools)


min0 <- function(x){
  min(x[x>0])
}

expand_hist_mids <- function(hist){
  binsize <- abs(hist$mids[2]-hist$mids[1])
  mids <- c(hist$mids-binsize,hist$mids[length(hist$mids)],hist$mids[length(hist$mids)]+binsize)
  return(mids)
}

expand_hist_counts <- function(hist){
  binsize <- abs(hist$mids[2]-hist$mids[1])
  counts <- c(0,hist$counts,0)
  return(counts)
}

plotdim <- function(x){
  if(x>1){
    xl <- round(sqrt(x)+0.5)
    yl <- ifelse(sqrt(x) %% 1<0.5,xl-1,xl)
    return(c(xl,yl))
  }else{
    return(c(1,1))
  }
}

countSpaces <- function(string){
  nchar(string) - nchar( sub(" *", "", string))
}

#graph_to_newick function by Christian Koenig
#https://gist.github.com/ChrKoenig
graph_to_newick <- function(graph, root){
  ##### Define callback functions for DFS
  # 1. function to be called whenever a vertex is visited
  f.in <- function(graph, data, extra) {
    curr_node = extra$names[data['vid']+1] # Get vertex name (Add 1 to index because igraph uses 0-based indexing)
    prev_node = extra$order.in[which(extra$order.in == curr_node)-1]
    next_node = extra$order.out[which(extra$order.out == curr_node)+1]
    if(length(extra$distances[prev_node]) == 0){ # first node / root
      cat("")
    } else{
      if(extra$distances[prev_node] < extra$distances[curr_node]){cat("(")}
    }
    FALSE # Do not terminate
  }
  
  # 2. function to be called whenever the subtree of a vertex is completed
  f.out <- function(graph, data, extra) {
    curr_node = extra$names[data['vid']+1] # Get vertex name (Add 1 to index because igraph uses 0-based indexing)
    prev_node = extra$order.in[which(extra$order.in == curr_node)-1]
    next_node = extra$order.out[which(extra$order.out == curr_node)+1]
    if(length(extra$distances[prev_node]) == 0){ # first node / root
      cat("")
    } else if(is.na(next_node)){ # last / root
      cat(curr_node)
      cat(";")
    } else{
      cat(curr_node)
      if(extra$distances[next_node] < extra$distances[curr_node]){cat(")")}
      if(extra$distances[next_node] == extra$distances[curr_node]){cat(",")}
      if(extra$distances[next_node] > extra$distances[curr_node]){cat(",")}
    }
    FALSE # Do not terminate
  }
  
  ##### Instructions
  # Obtain extra arguments (distance from root, order of node visits etc.) by running a test DFS without callback functions
  tmp_dfs = graph.dfs(graph, root = which(names(V(graph)) == root), dist = T, order = T, order.out = T) 
  # Organize results as extra arguments for the actual DFS
  extra = list(names  = names(V(graph)), order.in = names(tmp_dfs$order), order.out = names(tmp_dfs$order.out), 
               distances = tmp_dfs$dist, maxdist = max(tmp_dfs$dist))
  # Run DFS with callback functions
  newick_string = capture.output(graph.dfs(graph, root = which(names(V(graph)) == root), 
                                           dist = T, father = T, in.callback = f.in, out.callback = f.out, extra = extra))[1]
  # Tidy up
  newick_string = gsub(pattern = "\\$root", paste(root, ";", sep  = ""), newick_string)
  return(newick_string)
}

plot.phylo2 <- function (x, type = "phylogram", use.edge.length = TRUE, node.pos = NULL, 
                         show.tip.label = TRUE, show.node.label = FALSE, edge.color = "black", 
                         edge.width = 1, edge.lty = 1, font = 3, cex = par("cex"), 
                         adj = NULL, srt = 0, no.margin = FALSE, root.edge = FALSE, 
                         label.offset = 0, underscore = FALSE, x.lim = NULL, y.lim = NULL, 
                         direction = "rightwards", lab4ut = NULL, tip.color = "black", 
                         plot = TRUE, rotate.tree = 0, open.angle = 0, node.depth = 1, 
                         align.tip.label = FALSE, ...) 
{
  Ntip <- length(x$tip.label)
  if (Ntip < 2) {
    warning("found less than 2 tips in the tree")
    return(NULL)
  }
  .nodeHeight <- function(edge, Nedge, yy) .C(node_height, 
                                              as.integer(edge[, 1]), as.integer(edge[, 2]), as.integer(Nedge), 
                                              as.double(yy))[[4]]
  .nodeDepth <- function(Ntip, Nnode, edge, Nedge, node.depth) .C(node_depth, 
                                                                  as.integer(Ntip), as.integer(edge[, 1]), as.integer(edge[, 
                                                                                                                           2]), as.integer(Nedge), double(Ntip + Nnode), as.integer(node.depth))[[5]]
  .nodeDepthEdgelength <- function(Ntip, Nnode, edge, Nedge, 
                                   edge.length) .C(node_depth_edgelength, as.integer(edge[, 
                                                                                          1]), as.integer(edge[, 2]), as.integer(Nedge), as.double(edge.length), 
                                                   double(Ntip + Nnode))[[5]]
  Nedge <- dim(x$edge)[1]
  Nnode <- x$Nnode
  if (any(x$edge < 1) || any(x$edge > Ntip + Nnode)) 
    stop("tree badly conformed; cannot plot. Check the edge matrix.")
  ROOT <- Ntip + 1
  type <- match.arg(type, c("phylogram", "cladogram", "fan", 
                            "unrooted", "radial"))
  direction <- match.arg(direction, c("rightwards", "leftwards", 
                                      "upwards", "downwards"))
  if (is.null(x$edge.length)) {
    use.edge.length <- FALSE
  }
  else {
    if (use.edge.length && type != "radial") {
      tmp <- sum(is.na(x$edge.length))
      if (tmp) {
        warning(paste(tmp, "branch length(s) NA(s): branch lengths ignored in the plot"))
        use.edge.length <- FALSE
      }
    }
  }
  if (is.numeric(align.tip.label)) {
    align.tip.label.lty <- align.tip.label
    align.tip.label <- TRUE
  }
  else {
    if (align.tip.label) 
      align.tip.label.lty <- 3
  }
  if (align.tip.label) {
    if (type %in% c("unrooted", "radial") || !use.edge.length || 
        is.ultrametric(x)) 
      align.tip.label <- FALSE
  }
  if (type %in% c("unrooted", "radial") || !use.edge.length || 
      is.null(x$root.edge) || !x$root.edge) 
    root.edge <- FALSE
  phyloORclado <- type %in% c("phylogram", "cladogram")
  horizontal <- direction %in% c("rightwards", "leftwards")
  xe <- x$edge
  if (phyloORclado) {
    phyOrder <- attr(x, "order")
    if (is.null(phyOrder) || phyOrder != "cladewise") {
      x <- reorder(x)
      if (!identical(x$edge, xe)) {
        ereorder <- match(x$edge[, 2], xe[, 2])
        if (length(edge.color) > 1) {
          edge.color <- rep(edge.color, length.out = Nedge)
          edge.color <- edge.color[ereorder]
        }
        if (length(edge.width) > 1) {
          edge.width <- rep(edge.width, length.out = Nedge)
          edge.width <- edge.width[ereorder]
        }
        if (length(edge.lty) > 1) {
          edge.lty <- rep(edge.lty, length.out = Nedge)
          edge.lty <- edge.lty[ereorder]
        }
      }
    }
    yy <- numeric(Ntip + Nnode)
    TIPS <- x$edge[x$edge[, 2] <= Ntip, 2]
    yy[TIPS] <- 1:Ntip
  }
  z <- reorder(x, order = "postorder")
  if (phyloORclado) {
    if (is.null(node.pos)) 
      node.pos <- if (type == "cladogram" && !use.edge.length) 
        2
    else 1
    if (node.pos == 1) 
      yy <- .nodeHeight(z$edge, Nedge, yy)
    else {
      ans <- .C(node_height_clado, as.integer(Ntip), as.integer(z$edge[, 
                                                                       1]), as.integer(z$edge[, 2]), as.integer(Nedge), 
                double(Ntip + Nnode), as.double(yy))
      xx <- ans[[5]] - 1
      yy <- ans[[6]]
    }
    if (!use.edge.length) {
      if (node.pos != 2) 
        xx <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, 
                         node.depth) - 1
      xx <- max(xx) - xx
    }
    else {
      xx <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, 
                                 Nedge, z$edge.length)
    }
  }
  else {
    twopi <- 2 * pi
    rotate.tree <- twopi * rotate.tree/360
    if (type != "unrooted") {
      TIPS <- x$edge[which(x$edge[, 2] <= Ntip), 2]
      xx <- seq(0, twopi * (1 - 1/Ntip) - twopi * open.angle/360, 
                length.out = Ntip)
      theta <- double(Ntip)
      theta[TIPS] <- xx
      theta <- c(theta, numeric(Nnode))
    }
    switch(type, fan = {
      theta <- .nodeHeight(z$edge, Nedge, theta)
      if (use.edge.length) {
        r <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, 
                                  Nedge, z$edge.length)
      } else {
        r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, 
                        node.depth)
        r <- 1/r
      }
      theta <- theta + rotate.tree
      if (root.edge) r <- r + x$root.edge
      xx <- r * cos(theta)
      yy <- r * sin(theta)
    }, unrooted = {
      nb.sp <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, 
                          node.depth)
      XY <- if (use.edge.length) unrooted.xy(Ntip, Nnode, 
                                             z$edge, z$edge.length, nb.sp, rotate.tree) else unrooted.xy(Ntip, 
                                                                                                         Nnode, z$edge, rep(1, Nedge), nb.sp, rotate.tree)
      xx <- XY$M[, 1] - min(XY$M[, 1])
      yy <- XY$M[, 2] - min(XY$M[, 2])
    }, radial = {
      r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, node.depth)
      r[r == 1] <- 0
      r <- 1 - r/Ntip
      theta <- .nodeHeight(z$edge, Nedge, theta) + rotate.tree
      xx <- r * cos(theta)
      yy <- r * sin(theta)
    })
  }
  if (phyloORclado) {
    if (!horizontal) {
      tmp <- yy
      yy <- xx
      xx <- tmp - min(tmp) + 1
    }
    if (root.edge) {
      if (direction == "rightwards") 
        xx <- xx + x$root.edge
      if (direction == "upwards") 
        yy <- yy + x$root.edge
    }
  }
  if (no.margin) 
    par(mai = rep(0, 4))
  if (show.tip.label) 
    nchar.tip.label <- nchar(x$tip.label)
  max.yy <- max(yy)
  getLimit <- function(x, lab, sin, cex) {
    s <- strwidth(lab, "inches", cex = cex)
    if (any(s > sin)) 
      return(1.5 * max(x))
    Limit <- 0
    while (any(x > Limit)) {
      i <- which.max(x)
      alp <- x[i]/(sin - s[i])
      Limit <- x[i] + alp * s[i]
      x <- x + alp * s
    }
    Limit
  }
  if (is.null(x.lim)) {
    if (phyloORclado) {
      if (horizontal) {
        xx.tips <- xx[1:Ntip]
        if (show.tip.label) {
          pin1 <- par("pin")[1]
          tmp <- getLimit(xx.tips, x$tip.label, pin1, 
                          cex)
          tmp <- tmp + label.offset
        }
        else tmp <- max(xx.tips)
        x.lim <- c(0, tmp)
      }
      else x.lim <- c(1, Ntip)
    }
    else switch(type, fan = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy * 
                        cex)
        x.lim <- range(xx) + c(-offset, offset)
      } else x.lim <- range(xx)
    }, unrooted = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy * 
                        cex)
        x.lim <- c(0 - offset, max(xx) + offset)
      } else x.lim <- c(0, max(xx))
    }, radial = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.03 * cex)
        x.lim <- c(-1 - offset, 1 + offset)
      } else x.lim <- c(-1, 1)
    })
  }
  else if (length(x.lim) == 1) {
    x.lim <- c(0, x.lim)
    if (phyloORclado && !horizontal) 
      x.lim[1] <- 1
    if (type %in% c("fan", "unrooted") && show.tip.label) 
      x.lim[1] <- -max(nchar.tip.label * 0.018 * max.yy * 
                         cex)
    if (type == "radial") 
      x.lim[1] <- if (show.tip.label) 
        -1 - max(nchar.tip.label * 0.03 * cex)
    else -1
  }
  if (phyloORclado && direction == "leftwards") 
    xx <- x.lim[2] - xx
  if (is.null(y.lim)) {
    if (phyloORclado) {
      if (horizontal) 
        y.lim <- c(1, Ntip)
      else {
        pin2 <- par("pin")[2]
        yy.tips <- yy[1:Ntip]
        if (show.tip.label) {
          tmp <- getLimit(yy.tips, x$tip.label, pin2, 
                          cex)
          tmp <- tmp + label.offset
        }
        else tmp <- max(yy.tips)
        y.lim <- c(0, tmp)
      }
    }
    else switch(type, fan = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy * 
                        cex)
        y.lim <- c(min(yy) - offset, max.yy + offset)
      } else y.lim <- c(min(yy), max.yy)
    }, unrooted = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy * 
                        cex)
        y.lim <- c(0 - offset, max.yy + offset)
      } else y.lim <- c(0, max.yy)
    }, radial = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.03 * cex)
        y.lim <- c(-1 - offset, 1 + offset)
      } else y.lim <- c(-1, 1)
    })
  }
  else if (length(y.lim) == 1) {
    y.lim <- c(0, y.lim)
    if (phyloORclado && horizontal) 
      y.lim[1] <- 1
    if (type %in% c("fan", "unrooted") && show.tip.label) 
      y.lim[1] <- -max(nchar.tip.label * 0.018 * max.yy * 
                         cex)
    if (type == "radial") 
      y.lim[1] <- if (show.tip.label) 
        -1 - max(nchar.tip.label * 0.018 * max.yy * 
                   cex)
    else -1
  }
  if (phyloORclado && direction == "downwards") 
    yy <- y.lim[2] - yy
  if (phyloORclado && root.edge) {
    if (direction == "leftwards") 
      x.lim[2] <- x.lim[2] + x$root.edge
    if (direction == "downwards") 
      y.lim[2] <- y.lim[2] + x$root.edge
  }
  asp <- if (type %in% c("fan", "radial", "unrooted")) 
    1
  else NA
  plot.default(0, type = "n", xlim = x.lim, ylim = y.lim, 
               xlab = "", ylab = "", axes = FALSE, asp = asp, ...)
  if (plot) {
    if (is.null(adj)) 
      adj <- if (phyloORclado && direction == "leftwards") 
        1
    else 0
    if (phyloORclado && show.tip.label) {
      MAXSTRING <- max(strwidth(x$tip.label, cex = cex))
      loy <- 0
      if (direction == "rightwards") {
        lox <- label.offset + MAXSTRING * 1.05 * adj
      }
      if (direction == "leftwards") {
        lox <- -label.offset - MAXSTRING * 1.05 * (1 - 
                                                     adj)
      }
      if (!horizontal) {
        psr <- par("usr")
        MAXSTRING <- MAXSTRING * 1.09 * (psr[4] - psr[3])/(psr[2] - 
                                                             psr[1])
        loy <- label.offset + MAXSTRING * 1.05 * adj
        lox <- 0
        srt <- 90 + srt
        if (direction == "downwards") {
          loy <- -loy
          srt <- 180 + srt
        }
      }
    }
    if (type == "phylogram") {
      phylogram.plot(x$edge, Ntip, Nnode, xx, yy, horizontal, 
                     edge.color, edge.width, edge.lty)
    }
    else {
      if (type == "fan") {
        ereorder <- match(z$edge[, 2], x$edge[, 2])
        if (length(edge.color) > 1) {
          edge.color <- rep(edge.color, length.out = Nedge)
          edge.color <- edge.color[ereorder]
        }
        if (length(edge.width) > 1) {
          edge.width <- rep(edge.width, length.out = Nedge)
          edge.width <- edge.width[ereorder]
        }
        if (length(edge.lty) > 1) {
          edge.lty <- rep(edge.lty, length.out = Nedge)
          edge.lty <- edge.lty[ereorder]
        }
        outpoints <- circular.plot.print(z$edge, Ntip, Nnode, xx, yy, theta, 
                                         r, edge.color, edge.width, edge.lty)
      }
      else cladogram.plot(x$edge, xx, yy, edge.color, 
                          edge.width, edge.lty)
    }
    if (root.edge) {
      rootcol <- if (length(edge.color) == 1) 
        edge.color
      else "black"
      rootw <- if (length(edge.width) == 1) 
        edge.width
      else 1
      rootlty <- if (length(edge.lty) == 1) 
        edge.lty
      else 1
      if (type == "fan") {
        tmp <- polar2rect(x$root.edge, theta[ROOT])
        segments(0, 0, tmp$x, tmp$y, col = rootcol, 
                 lwd = rootw, lty = rootlty)
      }
      else {
        switch(direction, rightwards = segments(0, yy[ROOT], 
                                                x$root.edge, yy[ROOT], col = rootcol, lwd = rootw, 
                                                lty = rootlty), leftwards = segments(xx[ROOT], 
                                                                                     yy[ROOT], xx[ROOT] + x$root.edge, yy[ROOT], 
                                                                                     col = rootcol, lwd = rootw, lty = rootlty), 
               upwards = segments(xx[ROOT], 0, xx[ROOT], 
                                  x$root.edge, col = rootcol, lwd = rootw, 
                                  lty = rootlty), downwards = segments(xx[ROOT], 
                                                                       yy[ROOT], xx[ROOT], yy[ROOT] + x$root.edge, 
                                                                       col = rootcol, lwd = rootw, lty = rootlty))
      }
    }
    if (show.tip.label) {
      if (is.expression(x$tip.label)) 
        underscore <- TRUE
      if (!underscore) 
        x$tip.label <- gsub("_", " ", x$tip.label)
      if (phyloORclado) {
        if (align.tip.label) {
          xx.tmp <- switch(direction, rightwards = max(xx[1:Ntip]), 
                           leftwards = min(xx[1:Ntip]), upwards = xx[1:Ntip], 
                           downwards = xx[1:Ntip])
          yy.tmp <- switch(direction, rightwards = yy[1:Ntip], 
                           leftwards = yy[1:Ntip], upwards = max(yy[1:Ntip]), 
                           downwards = min(yy[1:Ntip]))
          segments(xx[1:Ntip], yy[1:Ntip], xx.tmp, yy.tmp, 
                   lty = align.tip.label.lty)
        }
        else {
          xx.tmp <- xx[1:Ntip]
          yy.tmp <- yy[1:Ntip]
        }
        text(xx.tmp + lox, yy.tmp + loy, x$tip.label, 
             adj = adj, font = font, srt = srt, cex = cex, 
             col = tip.color)
      }
      else {
        angle <- if (type == "unrooted") 
          XY$axe
        else atan2(yy[1:Ntip], xx[1:Ntip])
        lab4ut <- if (is.null(lab4ut)) {
          if (type == "unrooted") 
            "horizontal"
          else "axial"
        }
        else match.arg(lab4ut, c("horizontal", "axial"))
        xx.tips <- xx[1:Ntip]
        yy.tips <- yy[1:Ntip]
        if (label.offset) {
          xx.tips <- xx.tips + label.offset * cos(angle)
          yy.tips <- yy.tips + label.offset * sin(angle)
        }
        if (lab4ut == "horizontal") {
          y.adj <- x.adj <- numeric(Ntip)
          sel <- abs(angle) > 0.75 * pi
          x.adj[sel] <- -strwidth(x$tip.label)[sel] * 
            1.05
          sel <- abs(angle) > pi/4 & abs(angle) < 0.75 * 
            pi
          x.adj[sel] <- -strwidth(x$tip.label)[sel] * 
            (2 * abs(angle)[sel]/pi - 0.5)
          sel <- angle > pi/4 & angle < 0.75 * pi
          y.adj[sel] <- strheight(x$tip.label)[sel]/2
          sel <- angle < -pi/4 & angle > -0.75 * pi
          y.adj[sel] <- -strheight(x$tip.label)[sel] * 
            0.75
          text(xx.tips + x.adj * cex, yy.tips + y.adj * 
                 cex, x$tip.label, adj = c(adj, 0), font = font, 
               srt = srt, cex = cex, col = tip.color)
        }
        else {
          if (align.tip.label) {
            POL <- rect2polar(xx.tips, yy.tips)
            POL$r[] <- max(POL$r)
            REC <- polar2rect(POL$r, POL$angle)
            xx.tips <- REC$x
            yy.tips <- REC$y
            segments(xx[1:Ntip], yy[1:Ntip], xx.tips, 
                     yy.tips, lty = align.tip.label.lty)
          }
          if (type == "unrooted") {
            adj <- abs(angle) > pi/2
            angle <- angle * 180/pi
            angle[adj] <- angle[adj] - 180
            adj <- as.numeric(adj)
          }
          else {
            s <- xx.tips < 0
            angle <- angle * 180/pi
            angle[s] <- angle[s] + 180
            adj <- as.numeric(s)
          }
          font <- rep(font, length.out = Ntip)
          tip.color <- rep(tip.color, length.out = Ntip)
          cex <- rep(cex, length.out = Ntip)
          for (i in 1:Ntip) text(xx.tips[i], yy.tips[i], 
                                 x$tip.label[i], font = font[i], cex = cex[i], 
                                 srt = angle[i], adj = adj[i], col = tip.color[i])
        }
      }
    }
    if (show.node.label) 
      text(xx[ROOT:length(xx)] + label.offset, yy[ROOT:length(yy)], 
           x$node.label, adj = adj, font = font, srt = srt, 
           cex = cex)
  }
  L <- list(type = type, use.edge.length = use.edge.length, 
            node.pos = node.pos, node.depth = node.depth, show.tip.label = show.tip.label, 
            show.node.label = show.node.label, font = font, cex = cex, 
            adj = adj, srt = srt, no.margin = no.margin, label.offset = label.offset, 
            x.lim = x.lim, y.lim = y.lim, direction = direction, 
            tip.color = tip.color, Ntip = Ntip, Nnode = Nnode, root.time = x$root.time, 
            align.tip.label = align.tip.label)
  assign("last_plot.phylo", c(L, list(edge = xe, xx = xx, 
                                      yy = yy)), envir = .PlotPhyloEnv)
  invisible(L)
  outpoints
}

circular.plot.print <- function (edge, Ntip, Nnode, xx, yy, theta, r, edge.color, edge.width, 
                                 edge.lty) 
{
  r0 <- r[edge[, 1]]
  r1 <- r[edge[, 2]]
  theta0 <- theta[edge[, 2]]
  costheta0 <- cos(theta0)
  sintheta0 <- sin(theta0)
  x0 <- r0 * costheta0
  y0 <- r0 * sintheta0
  x1 <- r1 * costheta0
  y1 <- r1 * sintheta0
  segments(x0, y0, x1, y1, col = edge.color, lwd = edge.width, 
           lty = edge.lty)
  outm <- cbind(x0, y0, x1, y1)
  tmp <- which(diff(edge[, 1]) != 0)
  start <- c(1, tmp + 1)
  Nedge <- dim(edge)[1]
  end <- c(tmp, Nedge)
  foo <- function(edge.feat, default) {
    if (length(edge.feat) == 1) 
      return(as.list(rep(edge.feat, Nnode)))
    edge.feat <- rep(edge.feat, length.out = Nedge)
    feat.arc <- as.list(rep(default, Nnode))
    for (k in 1:Nnode) {
      tmp <- edge.feat[start[k]]
      if (tmp == edge.feat[end[k]]) {
        feat.arc[[k]] <- tmp
      }
      else {
        if (nodedegree[k] == 2) 
          feat.arc[[k]] <- rep(c(tmp, edge.feat[end[k]]), 
                               each = 50)
      }
    }
    feat.arc
  }
  nodedegree <- tabulate(edge[, 1L])[-seq_len(Ntip)]
  co <- foo(edge.color, "black")
  lw <- foo(edge.width, 1)
  ly <- foo(edge.lty, 1)
  for (k in 1:Nnode) {
    i <- start[k]
    j <- end[k]
    X <- rep(r[edge[i, 1]], 100)
    Y <- seq(theta[edge[i, 2]], theta[edge[j, 2]], length.out = 100)
    x <- X * cos(Y)
    y <- X * sin(Y)
    x0 <- x[-100]
    y0 <- y[-100]
    x1 <- x[-1]
    y1 <- y[-1]
    segments(x0, y0, x1, y1, col = co[[k]], lwd = lw[[k]], 
             lty = ly[[k]])
  }
  outm
}


###################################################################################################
## Files from snakemake
###################################################################################################

stats_image         <- snakemake@input #"Stats/all_stats.Rdata"
outplot             <- snakemake@output # 


###################################################################################################
## Load image 
###################################################################################################
print("Loading data.")
load(stats_image)

names(outplot) <- c(types,"sum")

krakRep <- list()
tax_tab <- list()
for(type in types){
  krakRep[[type]] <- read.delim(krakStats_files[[type]],stringsAsFactors=F,header=F,strip.white=F)
  
  tax_vec <- vector("character",100)
  tax_tab[[type]] <- list()
  for(i in 2:nrow(krakRep[[type]])){
    tax <- sub(" *","",krakRep[[type]]$V6[i])
    r <- krakRep[[type]]$V4[i]
    ftax <- paste0(r,"__",tax)
    l <- 1 + countSpaces(krakRep[[type]]$V6[i])/2
    tax_vec[l] <- ftax
    if(i==2){
      tot <- krakRep[[type]]$V2[i]
    }else{
    if(krakRep[[type]]$V2[i]>tot/2000 ){
      if(!"edgel" %in%ls()){
        edgel <- matrix(tax_vec[l-1:0],nrow=1,ncol=2)
      }else{
        edgel <- rbind(edgel,matrix(tax_vec[l-1:0],nrow=1,ncol=2))
      }
    }
    }
    stax <- paste(tax_vec[1:l],sep=";",collapse=";")
    if(! r %in% names(tax_tab[[type]])){
      tax_tab[[type]][[r]] <- data.frame("taxString"=stax,
                                 "reads"=krakRep[[type]]$V2[i],
                                 stringsAsFactors = F)
    }else{
      tax_tab[[type]][[r]] <- rbind(tax_tab[[type]][[r]],
                            data.frame("taxString"=stax,
                                       "reads"=krakRep[[type]]$V2[i],
                                       stringsAsFactors = F))
    }
  }
  krakRep[[type]]$newName <- paste0(krakRep[[type]]$V4,"__",
                                    gsub("[,/(/);:]","",sub("^ *","",krakRep[[type]]$V6)))
}
edgel <- unique(edgel)
gra <- graph_from_edgelist(as.matrix(edgel),directed = T)
V(gra)$name <- gsub("[,/(/);:]","",V(gra)$name)
df_tree = collapse.singles(read.newick(text = graph_to_newick(gra,"R__root")))
df_tree$edge.length <- rep(0.1,nrow(df_tree$edge))

###################################################################################################
## plot tree
###################################################################################################

relReads <- list()
for(type in types){
  plotfile <- outplot[type]
  relReads[[type]] <- log10(0.1 + sapply(c(df_tree$tip.label,df_tree$node.label)[reorder(df_tree, order = "postorder")$edge[,2]],
                                         function(x){
                                           if(x %in% krakRep[[type]]$newName){
                                             krakRep[[type]]$V2[krakRep[[type]]$newName==x]}else{
                                               0}})/max(krakRep[[type]]$V2[-1]))
  col <- colorRampPalette(brewer.pal(9,"YlGnBu")[-1])(256)[cut(relReads[[type]],256)]
  png(plotfile,width=17,height=17,units = "cm",res=300,pointsize=9)
  par("mar"=c(0.2,0.2,0.2,0.2))
  a <- plot.phylo2(df_tree,type="fan",show.tip.label = T,edge.color = "black",edge.width =0.3,
                   cex=0.3,label.offset = 0.02)
  points(a[,3:4],lwd=0.2,
         bg=col,pch=21,cex=0.75*(1+relReads[[type]])+ifelse(relReads[[type]]==-1,0,0.4))
  dev.off()
}

###################################################################################################
## plot mg vs mt for every level or histogram of mg or mt
###################################################################################################
plotfile <- outplot[3]
png(plotfile,width=17,height=17,units = "cm",res=300,pointsize=7)
par("mar"=c(2.9,4,0.5,0.5),tcl="-0.3",mgp=c(1.7,0.5,0),lwd=0.5,mfrow=plotdim(length(unique(unlist(c(sapply(tax_tab,names)))))),cex=1)
if(length(types)>1){
  for(rank in unique(unlist(c(sapply(tax_tab,names))))){
    if(rank %in% names(tax_tab[[1]]) & rank %in% names(tax_tab[[2]])){
      perTaxM <-merge(tax_tab[[1]][[rank]],tax_tab[[2]][[rank]],by=1,all=T,suffixes=paste0(".",types))
      perTaxM[is.na(perTaxM)] <- 0
      if(nrow(perTaxM)>0){
        plot(perTaxM$reads.mg,perTaxM$reads.mt,
             pch=16,cex=0.3,log="xy",las=1,
             xlim=c(0.09,max(sapply(tax_tab,function(x) max(x[[rank]][,2])))*1.2),xaxs="i",
             ylim=c(0.09,max(sapply(tax_tab,function(x) max(x[[rank]][,2])))*1.2),yaxs="i",
             ylab="",
             xlab="counts (metaG)",cex.lab=9/7,axes=F,
             col=brewer.pal(9,"Set1")[8])
        points(perTaxM$reads.mg[perTaxM$reads.mt==0],
               rep(0.1,length(which(perTaxM$reads.mt==0))),
               pch=16,cex=0.2,
               col=brewer.pal(9,"Set1")[8])
        points(rep(0.1,length(which(perTaxM$reads.mg==0))),
               perTaxM$reads.mt[perTaxM$reads.mg==0],
               pch=16,cex=0.2,
               col=brewer.pal(9,"Set1")[8])
        mtext("counts (metaT)",2,2.8,cex=9/7)
        mtext(paste0(rank,"__"),3,-1,cex=9/7)
        axis(1,
             at=c(0.09,10^c(0:10)),
             labels=c(0,format(10^c(0:10),scientific=T,digits = 2)),lwd=0.5,las=1)
        axis(2,at=c(0.09,10^c(0:10)),
             labels=c(0,format(10^c(0:10),scientific=T,digits = 2)),lwd=0.5,las=1)
        box(bty="o")
      }
    }else{
      ttyp <- ifelse(rank %in% names(tax_tab[[1]]),1,2)
      if(nrow(tax_tab[[ttyp]][[rank]])>0){
        depth_hist_log <- hist(log10(tax_tab[[ttyp]][[rank]][,2]),
                               breaks=seq(from=0,
                                          to=log10(max(tax_tab[[ttyp]][[rank]][,2])*1.2),
                                          length.out = 100),plot=F)
        plot(expand_hist_mids(depth_hist_log),
             expand_hist_counts(depth_hist_log),
             type="l",axes=F,ylab="",xlab=paste0("read counts ",ifelse(ttyp=="mg","(metaG)","(metaT)")),
             cex.lab=9/7,col=brewer.pal(9,"Set1")[8],
             xlim=c(expand_hist_mids(depth_hist_log)[1]-1.03,
                    max(expand_hist_mids(depth_hist_log))*1.03),xaxs="i",
             ylim=c(0,max(c(max(depth_hist_log$counts),
                            length(which(tax_tab[[ttyp]][[rank]][,2]==0))))*1.1),yaxs="i")
        mtext("counts",2,2.8,cex=9/7)
        mtext(paste0(rank,"__"),3,-1,cex=9/7)
        axis(1,at=c(0:10),
             labels = format(c(10^(0:10)),scientific=T),lwd=0.5,las=1)
        axis(2,lwd=0.5,las=1)
        box(bty="o")
      }
    }
  }
}else{
  for(rank in unique(unlist(c(sapply(tax_tab,names))))){
    if(nrow(tax_tab[[1]][[rank]])>0){
      depth_hist_log <- hist(log10(tax_tab[[1]][[rank]][,2]),
                             breaks=seq(from=0,
                                        to=log10(max(tax_tab[[1]][[rank]][,2])*1.2),
                                        length.out = 100),plot=F)
      plot(expand_hist_mids(depth_hist_log),
           expand_hist_counts(depth_hist_log),
           type="l",axes=F,ylab="",xlab=paste0("read counts ",ifelse(types[1]=="mg","(metaG)","(metaT)")),
           cex.lab=9/7,col=brewer.pal(9,"Set1")[8],
           xlim=c(expand_hist_mids(depth_hist_log)[1]-1.03,
                  max(expand_hist_mids(depth_hist_log))*1.03),xaxs="i",
           ylim=c(0,max(c(max(depth_hist_log$counts),
                          length(which(tax_tab[[1]][[rank]][,2]==0))))*1.1),yaxs="i")
      mtext("counts",2,2.8,cex=9/7)
      mtext(paste0(rank,"__"),3,-1,cex=9/7)
      axis(1,at=c(0:10),
           labels = format(c(10^(0:10)),scientific=T),lwd=0.5,las=1)
      axis(2,lwd=0.5,las=1)
      box(bty="o")
    }
  }
}
dev.off()

