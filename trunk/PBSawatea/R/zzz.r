.onLoad <- function(lib,pkg)
{
	pkg_info = utils::sessionInfo( package="PBSawatea" )$otherPkgs$PBSawatea
	if( is.character( pkg_info$Packaged ) )
		pkg_date <- strsplit( pkg_info$Packaged, " " )[[1]][1]
	else
		pkg_date  <- date()
	
	userguide_path <- system.file( "doc/PBSawatea.pdf", package = "PBSawatea" )
	
	packageStartupMessage("
-----------------------------------------------------------
PBS Awatea ", pkg_info$Version, " -- Copyright (C) 2011-12 Fisheries and Oceans Canada

A rough guide 'PBSawatea.pdf' is located at 
", userguide_path, "

Packaged on ", pkg_date, "
Pacific Biological Station, Nanaimo

Aotearoa, six months in a leaky boat...
-----------------------------------------------------------

")
}
# No Visible Bindings
# ===================
if(getRversion() >= "2.15.1") utils::globalVariables(names=c(
	"B0.MCMC","B0.mpd","B0refLevels","B0refNames","boxwidth",
	"currentMCMC","currentMSY","currentProj","currentRes",
	"delim",
	"fval",
	"gfcode","global",
	"h.mpd",
	"mainTitle","maxgrad","mess","minCpueYr",
	"npars",
	"obj",
	"plt.ageResids","policy","priorBoundsList","priorDistList",
	"R0.mpd","refPointsList","refs","resFileList","rpType",
	"Series","sigmaR",
	"variable"),
	package="PBSawatea")
