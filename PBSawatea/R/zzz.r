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
	"addLabel","addLegend",
	"B0.MCMC","B0.mpd","B0refLevels","B0refNames","blimYrs","boxwidth","btarYrs",
	"Cnames","CAnames","clipVector","cordat","currentMCMC","currentMSY","currentProj","currentRes","currYear",
	"delim",
	"evadat","expandGraph",
	"findPat","fval",
	"genMatrix","getYes","gfcode","global",
	"h.mpd",
	"istock",
	"J",
	"likdat","lucent",
	"MAexp","mainTitle","MAobs","maxgrad","mess","minCpueYr",
	"N","Ncpue","npars","Nsurv",
	"obj",
	"packList","pad0","pardat","PBSawatea","plotBubbles","plt.ageResids","policy","priorBoundsList","priorDistList","priorInput",
	"quants3","quants5",
	"R0.mpd","refPointsList","refPointsHistList","refs","resetGraph","resFileList","rpType",
	"SAnames","sen.lab","Series","sexlab","show0","sigmaR","Snames","ssnames","startYear","stddat",
	"tcall","trevorMCMC",
	"ulimYrs","unpackList","use.Pnames","useCA","useSA","utarYrs",
	"variable","Vexp",
	"years",
	"z1z2"
	), package="PBSawatea")

