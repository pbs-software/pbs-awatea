# Taking cue from Roger Bivand's maptools:
.PBSawaEnv <- new.env(FALSE, parent=globalenv())  # be sure to exportPattern("^\\.PBS") in NAMESPACE

.onAttach <- function(lib,pkg)
{
	pkg_info = utils::sessionInfo( package="PBSawatea" )$otherPkgs$PBSawatea
	if( is.character( pkg_info$Packaged ) )
		pkg_date <- strsplit( pkg_info$Packaged, " " )[[1]][1]
	else
		pkg_date  <- date()
	
	userguide_path <- system.file( "doc/PBSawatea.pdf", package = "PBSawatea" )
	year <- substring(date(),nchar(date())-3,nchar(date()))

	packageStartupMessage("
-----------------------------------------------------------
PBS Awatea ", pkg_info$Version, " -- Copyright (C) 2011-",year," Fisheries and Oceans Canada

A rough guide 'PBSawatea.pdf' is located at 
", userguide_path, "

Packaged on ", pkg_date, "
Pacific Biological Station, Nanaimo

All available PBS packages can be found at
https://github.com/pbs-software

Aotearoa, six months in a leaky boat...
-----------------------------------------------------------

")
}
.onUnload <- function(libpath) {
	rm(.PBSawaEnv)
}

# No Visible Bindings
# ===================
if(getRversion() >= "2.15.1") utils::globalVariables(names=c(
	".findSquare",
	"addLabel","addLegend","assYrs",
	"B0.MCMC","B0.mpd","B0refLevels","B0refNames","blimYrs","boxpars","boxwidth","btarYrs",
	"Cnames","CAnames","clipVector","cordat","createFdir","currentMCMC","currentMSY","currentProj","currentRes","currYear",
	"delim",
	"evadat","expandGraph",
	"findPat","findPV","fval",
	"genMatrix","getYes","gfcode","global",
	"h.mpd",
	"istock",
	"J",
	"likdat","linguaFranca","lucent",
	"MAexp","mainTitle","MAobs","maxcol","maxgrad","mess","minCpueYr",
	"N","NCAset","Ncpue","npars","Nsurv",
	"obj",
	"packList","pad0","pardat","PBSawatea","plotBubbles","plt.ageResids","policy","priorBoundsList","priorDistList","priorInput",
	"qboxplot","quants3","quants5",
	"R0.mpd","refPointsList","refPointsHistList","refs","resetGraph","resFileList","rpType","rwtNo",
	"SAnames","sen.lab","Series","sexlab","show0","sigmaR","Snames","ssnames","startYear","stddat",
	"tcall","texThatVec","trevorMCMC",
	"ulimYrs","unpackList","use.Pnames","useCA","useSA","utarYrs",
	"variable","Vexp",
	"years",
	"z1z2"
	), package="PBSawatea")

## quantBox-----------------------------2018-04-03
##  Redefine boxplot to show quantiles (RH 150910)
##  http://r.789695.n4.nabble.com/Box-plot-with-5th-and-95th-percentiles-instead-of-1-5-IQR-problems-implementing-an-existing-solution-td3456123.html
##  Use PBStools solution without requiring PBStools.
##  This needs to be in `zzz.r' for package compilation but will repeat
##  in `plotFuns.r' for running code locally without loading package.
## ---------------------------------------------RH
local(envir=.PBSmodEnv,expr={
	myboxplot.stats <- function (x, coef=NULL, do.conf=TRUE, do.out=TRUE)
	{
		nna <- !is.na(x)
		n <- sum(nna)
		if (!exists("quants5"))
			quants5 = c(0.05,0.25,0.50,0.75,0.95)
		stats <- quantile(x, quants5, na.rm=TRUE) ## one day figure out how to make this dynamic
		iqr <- diff(stats[c(2, 4)])
		out <- x < stats[1] | x > stats[5]
		conf <- if (do.conf)
			stats[3] + c(-1.58, 1.58) * diff(stats[c(2, 4)])/sqrt(n)
		list(stats = stats, n = n, conf = conf, out = x[out & nna])
	}
	boxcode = deparse(boxplot.default)
	boxcode = gsub("boxplot\\.stats","tcall(myboxplot.stats)",boxcode)
	eval(parse(text=c("qboxplot=",boxcode)))
})

