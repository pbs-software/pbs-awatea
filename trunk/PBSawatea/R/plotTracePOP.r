# plotTracePOP.r - now adding running median, and taking off overall
#  median and lowess line. Using cquantile from cumuplot.
# trying to add in the MPD as a big circle for
#  trace plots. 20th Oct 2010 (20/10/2010!)

cquantile <- function(z, probs)  # cumulative quantile, from cumuplot
  {
  cquant <- matrix(0, nrow = length(z), length(probs))
  for (i in seq(along = z)) if (is.R())
    {
    cquant[i, ] <- quantile(z[1:i], probs = probs, names = FALSE)
    }
  else {
        cquant[i, ] <- quantile(z[1:i], probs = probs)
       }
  cquant <- as.data.frame(cquant)
  names(cquant) <- paste(formatC(100 * probs, format = "fg", 
      wid = 1, digits = 7), "%", sep = "")
  return(cquant)
}

# AME doing this, just do one prob at a time (so it returns a vector
#  not a matrix)
cquantile.vec <- function(z, prob)  # cumulative quantile of vector
  {                                 #  prob is a single number
  cquant <- rep(NA, length(z))
  if(length(prob) != 1) stop("length prob should be 1")
  for (i in 1:length(z))
    {
    cquant[i] <- quantile(z[1:i], probs = prob, names = FALSE)
    }
  return(cquant)
}


plotTracePOP = function (mcmc, axes = FALSE, same.limits = FALSE, between = list(x = axes, 
    y = axes), div = 1, span = 1/4, log = FALSE, base = 10, main = NULL, 
    xlab = NULL, ylab = NULL, cex.main = 1.2, cex.lab = 1, cex.strip = 0.8, 
    cex.axis = 0.8, las = 0, tck = 0.5, tick.number = 5, lty.trace = 1, 
    lwd.trace = 1, col.trace = "grey", lty.median = 1, lwd.median = 1, 
    col.median = "black", lty.quant = 2, lwd.quant = 1, col.quant = "black", 
    plot = TRUE, probs=c(0.025, 0.5, 0.975) ,...)  # AME probs
{
    panel.trace <- function(x, y, ...) {
        panel.xyplot(x, y, type = "l", lty = lty.trace, lwd = lwd.trace, 
            col = col.trace)
        if (any(is.finite(y)) && var(y) > 0) {
            # print(x)  # gives 1 2 3 ... 1000 for each parameter/yr
            # panel.xyplot(range(x), rep(median(y), 2), type = "l", 
            #  lty = lty.median, lwd = lwd.median, col = col.median)
            panel.xyplot(x, cquantile.vec(y, prob=0.025),
              type = "l", lty = lty.quant, lwd = lwd.quant,
              col = col.quant)
            panel.xyplot(x, cquantile.vec(y, prob=0.5),
              type = "l", lty = lty.median, lwd = lwd.median,
              col = col.median)
            panel.xyplot(x, cquantile.vec(y, prob=0.975),
              type = "l", lty = lty.quant, lwd = lwd.quant,
              col = col.quant)
            panel.xyplot(x[1], y[1], pch=19, col="red") # AME
            panel.xyplot(x[1], y[1], pch=1, col="black") 
                     # AME, based on plt.trace, assume x[1]=1
            # suppressWarnings(panel.loess(x, y, span = span,
            #  lty = lty.loess, lwd = lwd.loess, col =col.loess,...))
        }
    }
    relation <- if (same.limits) 
        "same"
    else "free"
    if (is.null(dim(mcmc))) {
        mcmc.name <- rev(as.character(substitute(mcmc)))[1]
        mcmc <- matrix(mcmc, dimnames = list(NULL, mcmc.name))
    }
    mcmc <- if (log) 
        log(mcmc/div, base = base)
    else mcmc/div
    mcmc <- as.data.frame(mcmc)
    n <- nrow(mcmc)
    p <- ncol(mcmc)
    x <- data.frame(Factor = ordered(rep(names(mcmc), each = n), 
        names(mcmc)), Draw = rep(1:n, p), Value = as.vector(as.matrix(mcmc)))
    require(grid, quiet = TRUE, warn = FALSE)
    require(lattice, quiet = TRUE, warn = FALSE)
    if (trellis.par.get()$background$col == "#909090") {
        for (d in dev.list()) dev.off()
        trellis.device(color = FALSE)
    }
    mymain <- list(label = main, cex = cex.main)
    myxlab <- list(label = xlab, cex = cex.lab)
    myylab <- list(label = ylab, cex = cex.lab)
    myrot <- switch(as.character(las), `0` = 90, `1` = 0, `2` = 0, 
        `3` = 90)
    myscales <- list(x = list(draw = FALSE), y = list(draw = axes, 
        relation = relation, cex = cex.axis, tck = tck, tick.number = tick.number, 
        rot = myrot))
    mystrip <- list(cex = cex.strip)
    graph <- xyplot(Value ~ Draw | Factor, panel = panel.trace, 
        data = x, as.table = TRUE, between = between, main = mymain, 
        xlab = myxlab, ylab = myylab, par.strip.text = mystrip, 
        scales = myscales, ...)
    if (plot) {
        print(graph)
        invisible(x)
    }
    else {
        invisible(graph)
    }
}


