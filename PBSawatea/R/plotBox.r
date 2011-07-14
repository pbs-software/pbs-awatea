#plotBox--------------------------------2011-07-13
# Modified boxplot with quantile whiskers.
#-----------------------------------------------RH
plotBox = function (x, ..., range=1.5, width=NULL, varwidth=FALSE, 
    notch=FALSE, outline=TRUE, names, plot=TRUE, border=par("fg"), 
    col=NULL, log="", pars=list(boxwex=0.8, staplewex=0.5, 
        outwex=0.5, whisklty=1 ), horizontal=FALSE, add=FALSE, at=NULL,
        quants=c(0.025,0.25,0.5,0.75,0.975), outliers=FALSE) 
{
    args <- list(x, ...)
    namedargs <- if (!is.null(attributes(args)$names)) 
        attributes(args)$names != ""
    else rep(FALSE, length.out = length(args))
    groups <- if (is.list(x)) 
        x
    else args[!namedargs]
    if (0L == (n <- length(groups))) 
        stop("invalid first argument")
    if (length(class(groups))) 
        groups <- unclass(groups)
    if (!missing(names)) 
        attr(groups, "names") <- names
    else {
        if (is.null(attr(groups, "names"))) 
            attr(groups, "names") <- 1L:n
        names <- attr(groups, "names")
    }
    cls <- sapply(groups, function(x) class(x)[1L])
    cl <- if (all(cls == cls[1L])) 
        cls[1L]
    else NULL
    for (i in 1L:n) groups[i] <- list(boxplot.stats(unclass(groups[[i]]), 
        range))
    stats <- matrix(0, nrow = 5L, ncol = n)
    conf <- matrix(0, nrow = 2L, ncol = n)
    ng <- out <- group <- numeric(0L)
    ct <- 1
    for (i in groups) {
        stats[, ct] <- i$stats
        conf[, ct] <- i$conf
        ng <- c(ng, i$n)
        if ((lo <- length(i$out))) {
            out <- c(out, i$out)
            group <- c(group, rep.int(ct, lo))
        }
        ct <- ct + 1
    }
    #----RH tweaks for qantile plots----
    stats = sapply(x,quantile,quants,na.rm=TRUE)
    if (!outliers) {
    	out = NULL; group = NULL }
    #-----------------------------------
    if (length(cl) && cl != "numeric") 
        oldClass(stats) <- cl
    z <- list(stats = stats, n = ng, conf = conf, out = out, 
        group = group, names = names)
#browser();return()
    if (plot) {
        if (is.null(pars$boxfill) && is.null(args$boxfill)) 
            pars$boxfill <- col
        do.call("bxp", c(list(z, notch = notch, width = width, 
            varwidth = varwidth, log = log, border = border, 
            pars = pars, outline = outline, horizontal = horizontal, 
            add = add, at = at), args[namedargs]))
        invisible(z)
    }
    else z
}
#------------------------------------------plotBox
#plotBox(xBox)