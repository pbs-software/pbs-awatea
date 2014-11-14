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
http://code.google.com/p/pbs-software/

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
	"B0.MCMC","B0.mpd","B0refLevels","B0refNames","blimYrs","boxwidth","btarYrs",
	"Cnames","CAnames","cordat","currentMCMC","currentMSY","currentProj","currentRes","currYear",
	"delim",
	"evadat",
	"fval",
	"gfcode","global",
	"h.mpd",
	"J",
	"likdat",
	"MAexp","mainTitle","MAobs","maxgrad","mess","minCpueYr",
	"N","npars","Nsurv",
	"obj",
	"pardat","PBSawatea","plt.ageResids","policy","priorBoundsList","priorDistList","priorInput",
	"R0.mpd","refPointsList","refPointsHistList","refs","resFileList","rpType",
	"Series","Snames","SAnames","sexlab","sigmaR","stddat",
	"trevorMCMC",
	"ulimYrs","use.Pnames","useCA","useSA","utarYrs",
	"variable","Vexp",
	"years",
	"z1z2"
	), package="PBSawatea")

