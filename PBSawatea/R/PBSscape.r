##PBSscape-------------------------------2017-12-01
## Modified functions from Arni Magnussen's packages 'scape' and 'plotMCMC'
## ------------------------------------------------------------------------
##  plt.mpdGraphs   : wrapper function to call MPD graphing functions
##  plt.mcmcGraphs  : wrapper function to call MCMC graphing functions
#-------------------------------------------AME/RH
## Force this globally (without have to make a documentation file)
assign("quants3", c(0.05,0.50,0.95), envir=.GlobalEnv)
assign("quants5", c(0.05,0.25,0.50,0.75,0.95), envir=.GlobalEnv)

## Flush the cat down the console (change to '.flash.cat' to avoid conflict with function in PBStools)
## Note: `.flush.cat' already in PBStools but this package is not assumed to be loaded.
.flash.cat = function(...) { cat(...); flush.console() }

#---History---------------------------------------
# PBSscape.r - originally created from `ymrscape` for POP 2012.
#   Since then, it has been extensively revised by RH to handle
#   Rock Sole (single-sex model) and Silvergray Rockfish (2013).

# ymrScape.r - for ymr, just the function definitions here. Calls to
#   them will go into .Snw. Call this from .Snw. 23rd February 2011
# SEE popScapeRuns2.r for figures to do for both runs at once, as
#   went into the POP SAR. Wait until we have final results.

# popScape2.r - going to do all postscript files as better than .png.
#   Getting rid of menus and automatically doing all required figures
#   each time I run it, then automatically into latex to check, then
#   when all good put into the Word file with captions etc. Will be
#   way more efficient for me to work, rather than all the time it
#   takes opening and closing windows.
#   **Jump to end to load in MCMC and plot pairs plots and acf.

# popScape.r notes:
# Always choose one of the "Save all ** plots to PNG" options, 
#   as I'm not repeating the new figures. 
# Go through figures, see what else to add. Some will just need tweaking - see
#   notes on Appendix D printout.
# Checking that first value for MCMC chain is what's in currentRes as
#   the "MPD". Do get:
#   > currentRes$B[30, ]
#      Year     VB    SB       Y         U        R
#   30 1969 145424 88107 10382.4 0.0713941 17866.36
#   > currentMCMC$B[1, "1969"]
#   [1] 88107    # so this is correct (they call it SB then B)
#   But recruitment isn't - damn, maybe it's a year off again, 
#   like what I figured out for MCMC:
#   > currentMCMC$R[1, "1969"]
#   [1] 12773.1        # so disagrees
#   > currentMCMC$R[1, "1970"]
#   [1] 17866.4    # so it is a year off.
#   But recruits aren't actually plotted from MPD stuff, so think is okay.

# Plots like biomass.png use the MPD, but shouldn't they really use the median
#  of posterior (or mean), and maybe show credible intervals?? Check difference
#  between median of Bt and MPD.  Yes, they should, so doing that now


# popScape.r - developing further in assess7/, from version in
#  assess6/.  Had to use postscript for the recruitment marginal
#  posterior densities due to resizing. They will look fine when
#  the .doc is converted to .pdf.
#  Using the lattice approaches, not the three ....PJS.jpg sections
#   which are not lattice, and more fiddly to change. 18-Oct-10
# popScape.r - renaming carScapeAndy.r, and continuing to develop. Haven't
#  kept track of the many changes in the list below - could always do later
#  using CSDiff if really necessary. 1-Sep-10
# carScapeAndy.r - AME editing so that menu works for POP.
#  Search for AME for other changes. 12-Aug-10
# Added some plots: residuals by age w/o outliers, and residuals by age w/o outliers. May only work for save all MPD plots as doing quickly.   27th August 2010
# AME: all three types of plots of residuals of catch-at-age data,
#  with x-axis being year, age or year of birth, fill in NA's to give
#  white space for years/ages in which there are no data. 30th August
#  2010
# goto HERE - ask Michael for Trellis book.
# AME changing all wmf -> png
# For Paul: I fixed the problem with the loading in of .res files,
#  that menu now works. So first load in the .res file from the
#  main menu (I think this
#  could be simplified as it does automatically load in all .res
#  files in the directory, but for now we'll stick with this).
# For the MPD plots (option 2 on main menu): options 1-3 work.
#  For option 1 (Plot biomass, recruitment, catch),
#  we need to discuss whether to use the plotB or plotB2 function.
#  Currently plotB2 is used, for which below it says:
#  "function to accommodate PJS request not to show biomass
#  prior to fishery and survey indices period." Consequently the 
#  biomass estimates are only show for the later years.
#  option 2 just plots initial age structure and recruitments,
#  but by changing the options it can show more ages, so we can
#  maybe play with that.
#  option 4 (Plot commercial catch-at-age results) does not work, 
#  but the code mentions a bug that needs sorting.

#  option 5 now Plot survey catch-at-age results, following have
#   thus been incremented by 1.
#  option 6 (Plot abundance index) does not work, but the similar
#  command:
#    plotIndex(currentRes, what="s", xlim=c(1960, 2012))
#  does work for the survey indices, so I could just change the
#  the options to select something similar if you like. 
#  option 7 (All residual plots) doesn't give an error now.
#  options 8-10 plot results from each .res file in the directory
#  so that they can be easily compared - looks useful.
#
#
#--------------------------------------------------------------------#
# rsScape.r : Graphical analyses for rocksole Coleraine ouput.       #
# Developers: A.R. Kronlund, P.J. Starr                              #
# Modified by Allan C. Hicks to use output of Awatea forSPO7         #
# Required libraries: gdata, scape, plotMCMC, PBSmodelling          #
#                                                                    #
# Date Revised:                                                      #
# 08-Nov-05  Initial implementation.                                 #
# 10-Nov-05  Added generic residual plot plus stdRes calc functions. #
# 11-Nov-05  Added plot and save all functions.                      #
# 11-Nov-05  Added standardised age residuals a'la PJS doc.          #
# 12-Nov-05  Added res file selection function.                      #
# 12-Nov-05  Added MCMC load and selected plotting.                  #
# 13-Nov-05  Added biomass and projection plots. Probability tables. #
# 14-Nov-05  PJS adds a very small number of minor edits.            #
# 15-Nov-05  Addressed PJS improvements list.                        #
# 15-Nov-05  Added policy list at end of file so you don't have to   #
#            go hunting for it in the function code.                 #
# 15-Nov-05  Added console prompt to projection plots to step thru   #
#            multiple pages of plots depending on number of policies.#
# 17-Nov-05  Fixed age residuals, fixed multple fishery CA plots.    #
#            Fixed importCol to accommodate multiple gears.          #
#            Replace "importCol" with "importCol2" below.            #
# 21-Nov-05  Fixed bug in saving projection plots, Arni's scales to  #
#            trace plots.                                            #
# 21-Nov-05  Added quantile box plots for reconstruction-projections.#
# 22-Nov-05  Added observation-based reference points funciton.      #
# 23-Nov-05  Added import funs for PJS Delay Difference model output.#
# 24-Nov-05  Fixed bug in saved policyProjection plots.              #
# 25-Nov-05  Added new biomass trajectory plot to plot biomass only  #
#            over the period of tuning data.                         #
# 09-Dec-05  Added multi-panel plots for SSB/VB.                     #
# 09-Dec-05  Add PJS trace+mpd plots to rsScape.r.                   #
#            Correct bad recruits histogram and plot age 1's.        #
# 13-Dec-05  Finished PJS trace plots.                               #
# 16-Dec-05  Revised plots to not include the pre-1966 period.       #
#                                                                    #
# 16-Apr-06  Modified importCol2 to take Awatea output (ACH)         #
#            Added 0.5 to Year when plotting indices                 #
# 05-Sep-07  Minor modifications to plot age instead of length by    #                                    PJS                                       #
# 12-Aug-10  Commenting out mainMenu() to run from a script. AME.    #
# **-Aug-10  Many further changes specific to POP assessment. AME    #
# **-Sep-10   "                                                      #
# NOTE:                                                              #
# (1) For most plots there must be a "res" file loaded.              #
# (2) For MCMC plots you must have MCMC output loaded (*.pst files). #
# (3) For Projection plots and tables you must have MCMC and also    #
#     Projection output loaded.  See the MCMC menu.                  #
#                                                                    #
# To do:                                                             #
#                                                                    #
# (1) Check with Paul regarding calculations for age residuals.      #
#--------------------------------------------------------------------#
#                                                                    #
# Awatea res file has following structure (some elements may be      #
# missing dependent on model configuration and importCol details).   #
#                                                                    #
# N predicted numbers at age                                         #
# B predicted biomass, recruitment, and observed landings            #
# Sel predicted selectivity and observed maturity (age things)       #
# Dev predicted recruitment deviates from the stock-recruitment curve#
# CPUE, Survey commercial and survey abundance index and fit         #
# CAc, CAs commercial and survey C@A (catch at age) and fit          #
# CLc, CLs commercial and survey C@L (catch at length) and fit       #
# LA observed L@A and fit                                            #
#                                                                    #
# MCMC                                                               #
#                                                                    #
# The importProj function loads a list with elements "B" and "Y" for #
# biomass by catch policy and year, and catch by harvest policy and  #
# year.  The "B" element is itself a list of matrices with a matrix  #
# for each level of the catch.  This matrix has rows equal to the    #
# length of the chain and columns corresponding to projection years. #
# There are no plotting routines for these data.                     #
#--------------------------------------------------------------------#

#--------------------------------------------------------------------#
#                    Awatea Related Functions                     #
#--------------------------------------------------------------------#

# BUG FIX: Original importCol appears not to recognize multiple fishery series.

#plt.mpdGraphs--------------------------2018-04-17
# Plot the MPD graphs to encapsulated postscript files.
#-------------------------------------------------
# RH (2014-09-23)
#  Aded arguments `ptype',`pngres', and `ngear'
#  to render output figures in multiple formats
#  (only `eps' and `png' at the moment), and
#  to accommodate multiple gear types
#-------------------------------------------AME/RH
plt.mpdGraphs <- function(obj, save=FALSE, ssnames=paste("Ser",1:9,sep=""),
   ptypes=tcall(PBSawatea)$ptype, pngres=400, ngear=1,
   pchGear=seq(21,20+ngear,1), ltyGear=seq(1,ngear,1), 
   colGear=rep(c("black","blue"),ngear)[1:ngear])
{
	#AME some actually MCMC. # Doing as postscript now. # Taking some out for ymr.
	closeAllWin()

	# AME adding, plot exploitation rate, not writing new function:
	# RH modified to deal with multiple commercial gears

	B    = obj$B
	xlim = range(B$Year,na.rm=TRUE)
	U    = B[,grep("U",names(B)),drop=FALSE] # need to use `drop' argument for ngear=1
	ylim = range(U,na.rm=TRUE)
	for (p in ptypes) {
		if (p=="eps") postscript("exploit.eps", width=6.5, height=4.5, horizontal=FALSE,  paper="special")
		else if (p=="png") png("exploit.png", units="in", res=pngres, width=6.5, height=4.5)
		par(mfrow=c(1,1), mar=c(3,3.5,0.5,1), oma=c(0,0,0,0), mgp=c(2.4,0.5,0))
		plot(0,0,xlim=xlim,ylim=ylim,type="n",xlab="",ylab="Exploitation rate", las=1, cex.lab=1.2, tcl=-0.4)
		axis(1, at=intersect(seq(1900,3000,5),xlim[1]:xlim[2]), labels=FALSE, tcl=-0.2)
		mtext("Year", side=1, line=1.75, cex=1.2)
		abline(h=seq(0.05,0.95,0.05), lty=3, col="gainsboro")
		sapply(ngear:1,function(g,x,y){
			lines(x, y[,g], lty=ltyGear[g], col=colGear[g])
			exheat = list(lo = y[,g] <= 0.05, mid = y[,g] > 0.05 & y[,g] <= 0.10, hi = y[,g] > 0.10)
			for (j in 1:length(exheat)){
				zex = exheat[[j]]
				points(x[zex], y[,g][zex], cex=0.8, pch=pchGear[g], bg=switch(j,"gold","orange","red"), col=colGear[g])
			}
		}, x = B$Year, y = U)
		if (ngear>1) 
			legend("topleft",bty="n",lty=ltyGear,pch=pchGear,legend=Cnames,inset=0.05,seg.len=4,pt.bg="white",col=colGear)
			#legend("topleft",bty="n",lty=ltyGear,pch=pchGear,legend=paste0("gear ",1:ngear),inset=0.05,seg.len=4,pt.bg="white",col=colGear)
		if (p %in% c("eps","png")) dev.off()
	}

	## AME had added recruits.eps for POP, but that's MCMC, so moving
	##  to plt.mcmcGraphs, changing that filename to recruitsMCMC.eps
	##  and just adding here to do MPD for recruits.
	for (p in ptypes) {
		if (p=="eps") postscript("recruits.eps", width=6.5, height=4, horizontal=FALSE,  paper="special")
		else if (p=="png") png("recruits.png", units="in", res=pngres, width=6.5, height=4)
		par(mfrow=c(1,1), mar=c(3.25,3.5,1,1), oma=c(0,0,0,0), mgp=c(2,0.75,0))
		plot(obj$B$Year, obj$B$R, type="o", xlab="Year",
			ylab="Recruitment, Rt (1000s)", ylim=c(0, max(obj$B$R, na.rm=TRUE)))
		if (p %in% c("eps","png")) dev.off()
	}

	## Plot the selectivity.
	objRed = obj         ## Reduced object, just plotting Selectivity to age 20
	ageP = objRed$Sel$P; names(ageP)=objRed$Sel$Age
	selP = split(ageP,paste(objRed$Sel$Series,objRed$Sel$Sex,sep="."))
	xmax = max(as.numeric(sapply(selP,function(x){names(x[is.element(x,1)])[1]})),na.rm=TRUE) #maximum minimum age when P first hits 1
	if (is.na(xmax)) xmax = 20
	if (any(round(unlist(currentRes$extra$parameters[c("log_varRest","log_surveyvarR")]),5)!=100)) xmax=40 ## temporary fix
	objRed$Sel = objRed$Sel[objRed$Sel$Age <= xmax,]
	for (p in ptypes) {
		if (p=="eps")      postscript("selectivity.eps", width=6.5, height=4.5, horizontal=FALSE,  paper="special")
		else if (p=="png") png("selectivity.png", units="in", res=pngres, width=6.5, height=4.5)
		par(mfrow=c(1,1), mar=c(3.2,3.2,0.5,0.5), oma=c(0,0,0,0), mgp=c(2,0.75,0))
		plotSel( objRed, main=paste(mainTitle,"Selectivity"), xlim=c(0,xmax))
		if (p %in% c("eps","png")) dev.off()
	}

	## Plot the catch at age
	## NOTE: There is a bug in plotCA that prevents plotting multiple
	##       series given a list of character vectors in series.
	##         ACH: I'm not sure if this applies to CL

	# RH: Removed a bunch of commented text and code by AME
	#   Also transferred code manipulation for plotting MPD age fits
	#   into a new function called `plotAges' (located in `plotFuns.r')
	#   to handle both commercial and survey age fits.

	cyrs = sort(unique(obj[["CAc"]][["Year"]]))
	cgrp = split(cyrs,ceiling(seq_along(cyrs)/25))  ## RH (180417) changed maximum number of years on a page from 20 to 25
	if (length(cgrp)==1)
		plotAges(obj, what="c", maxcol=5, sexlab=sexlab, ptypes=ptypes, pngres=pngres, cex.lab=1.2, cex.axis=1, col.point=lucent("black",0.8), cex.points=0.4)
	else {
		for (i in 1:length(cgrp))
			plotAges(obj, what="c", maxcol=5, sexlab=sexlab, ptypes=ptypes, pngres=pngres, years=cgrp[[i]], set=LETTERS[i], cex.lab=1.2, cex.axis=1, col.point=lucent("black",0.8), cex.points=0.4)
	}
	plotAges(obj, what="s", maxcol=5, sexlab=sexlab, ptypes=ptypes, pngres=pngres, cex.lab=1.2, cex.axis=1, col.point=lucent("black",0.8), cex.points=0.4)

	## Plot the fishery index (CPUE data I think)
	## Plot the survey indices.

	## Now do on one plot as postscript:
	plotIndexNotLattice(obj, ssnames=ssnames, ptypes=ptypes, pngres=pngres ) # survey indices

	## Single plot of CPUE:
	plotCPUE(obj$CPUE, yLim=c(0, 200))

	# Plot standardised residuals. Now doing four plots on one page.
	# --------------------------------------------------------------
	#  Commercial.
	#  postscript("commAgeResids.eps", height=8.5, width=6.8, horizontal=FALSE,  paper="special")
	#  par(mai=c(0.45,0.55,0.1,0.1)) # JAE changed  for each figure
	#  par(omi=c(0,0,0.4,0))        # Outer margins of whole thing, inch
	#   Outliers don't get plotted, except for qq plot
	#  par(mfrow=c(4,1))
	# RH revamp because PJS wants this by sex (in addition to original irregardless of sex)
	objCAc = obj$CAc
	for (g in sort(unique(objCAc$Series))) { # treat each gear type separately
		objCAc.g = objCAc[is.element(objCAc$Series,g),]
		stdRes.CAc.g = stdRes.CA( objCAc.g )
		for (p in ptypes) {
			pname = paste0("commAgeResSer",g)
			if (p=="eps") postscript(paste0(pname,".eps"), width=6.5, height=8.5, horizontal=FALSE,  paper="special")
			else if (p=="png") png(paste0(pname,".png"), units="in", res=pngres, width=6.5, height=8.5)
			par(mfrow=c(4,1), mai=c(0.45,0.3,0.1,0.1), omi=c(0,0.25,0.4,0), mgp=c(2,0.75,0))
			plt.ageResidsPOP( stdRes.CAc.g, main="")
			mtext(CAnames[g],side=3,outer=TRUE,line=0.25,cex=1.5)
			mtext("Standardised Residuals",side=2,outer=TRUE,line=0,cex=1.5)
			plt.yearResidsPOP(stdRes.CAc.g)
			plt.cohortResids(stdRes.CAc.g)   # cohort resid, by year of birth
			plt.ageResidsqqPOP(stdRes.CAc.g)
			if (p %in% c("eps","png")) dev.off()
		}
		for (s in sort(unique(objCAc.g$Sex))) {
			objCAc.g.s = objCAc.g[is.element(objCAc.g$Sex,s),]
			stdRes.CAc.g.s = stdRes.CA( objCAc.g.s )
			for (p in ptypes) {
				pname = paste0("commAgeResSer",g,s)
				if (p=="eps") postscript(paste0(pname,".eps"), width=6.5, height=8.5, horizontal=FALSE,  paper="special")
				else if (p=="png") png(paste0(pname,".png"), units="in", res=pngres, width=6.5, height=8.5)
				par(mfrow=c(4,1), mai=c(0.45,0.3,0.1,0.1), omi=c(0,0.25,0.4,0), mgp=c(2,0.75,0))
				plt.ageResidsPOP( stdRes.CAc.g.s, main="" ) 
				mtext(paste0(CAnames[g]," - ",s),side=3,outer=TRUE,line=0.25,cex=1.5)
				mtext("Standardised Residuals",side=2,outer=TRUE,line=0,cex=1.5)
				plt.yearResidsPOP(stdRes.CAc.g.s)
				plt.cohortResids(stdRes.CAc.g.s)   # cohort resid, by year of birth
				plt.ageResidsqqPOP(stdRes.CAc.g.s)
				if (p %in% c("eps","png")) dev.off()
			}
		}
	}

	# And now for surveys.
	# AME adding - plotting the CA residuals for the two surveys:
	# RH revamp because PJS wants this by sex (in addition to original irregardless of sex)
	objCAs = obj$CAs
	for (g in sort(unique(objCAs$Series))) { # treat each survey separately (g = survey number)
		objCAs.g = objCAs[is.element(objCAs$Series,g),]
		stdRes.CAs.g = stdRes.CA( objCAs.g )
		for (p in ptypes) {
			pname = paste0("survAgeResSer",g)
			if (p=="eps") postscript(paste0(pname,".eps"), width=6.5, height=8.5, horizontal=FALSE,  paper="special")
			else if (p=="png") png(paste0(pname,".png"), units="in", res=pngres, width=6.5, height=8.5)
			par(mfrow=c(4,1), mai=c(0.45,0.3,0.1,0.1), omi=c(0,0.25,0.4,0), mgp=c(2,0.75,0))
			plt.ageResidsPOP( stdRes.CAs.g, main="")
			mtext(Snames[g],side=3,outer=TRUE,line=0.25,cex=1.5)  ## g indexes Snames not SAnames
			mtext("Standardised Residuals",side=2,outer=TRUE,line=0,cex=1.5)
			plt.yearResidsPOP(stdRes.CAs.g)
			plt.cohortResids(stdRes.CAs.g)   # cohort resid, by year of birth
			plt.ageResidsqqPOP(stdRes.CAs.g)
			if (p %in% c("eps","png")) dev.off()
		}
		for (s in sort(unique(objCAs.g$Sex))) {
			objCAs.g.s = objCAs.g[is.element(objCAs.g$Sex,s),]
			stdRes.CAs.g.s = stdRes.CA( objCAs.g.s )
			for (p in ptypes) {
				pname = paste0("survAgeResSer",g,s)
				if (p=="eps") postscript(paste0(pname,".eps"), width=6.5, height=8.5, horizontal=FALSE,  paper="special")
				else if (p=="png") png(paste0(pname,".png"), units="in", res=pngres, width=6.5, height=8.5)
				par(mfrow=c(4,1), mai=c(0.45,0.3,0.1,0.1), omi=c(0,0.25,0.4,0), mgp=c(2,0.75,0))
				plt.ageResidsPOP( stdRes.CAs.g.s, main="" ) 
				mtext(paste0(Snames[g]," - ",s),side=3,outer=TRUE,line=0.25,cex=1.5)  ## g indexes Snames not SAnames
				mtext("Standardised Residuals",side=2,outer=TRUE,line=0,cex=1.5)
				plt.yearResidsPOP(stdRes.CAs.g.s)
				plt.cohortResids(stdRes.CAs.g.s)   # cohort resid, by year of birth
				plt.ageResidsqqPOP(stdRes.CAs.g.s)
				if (p %in% c("eps","png")) dev.off()
			}
		}
	}
#	seriesList <- sort( unique( obj$CAs$Series) )  
#	nseries=length(seriesList)
#	surveyHeadName=if (!exists("tcall")) ssnames else tcall(PBSawatea)$Snames
#	for ( i in 1:nseries ) {   # POP no fits for survey 3
#		ii=seriesList[i]
#		stdRes.CA.CAs=stdRes.CA( obj$CAs[obj$CAs$Series == ii,] )
#		for (j in ptypes) {
#			pname = paste0("survAgeResSer", ii)
#			if (j=="eps") postscript(paste0(pname,".eps"), height=8.5, width=6.5, horizontal=FALSE,  paper="special")
#			else if (j=="png") png(paste0(pname,".png"), res=pngres, height=8.5*pngres, width=6.5*pngres)
#			par(mfrow=c(4,1), mai=c(0.45,0.55,0.1,0.1), omi=c(0,0,0.4,0), mgp=c(2,0.75,0))
#			plt.ageResidsPOP(stdRes.CA.CAs, main="" )
#			mtext(surveyHeadName[ii],side=3,outer=TRUE,line=0.25,cex=1.5)
#			plt.yearResidsPOP(stdRes.CA.CAs)
#			plt.cohortResids(stdRes.CA.CAs)   # cohort resid, by year of birth
#			plt.ageResidsqqPOP(stdRes.CA.CAs)
#			dev.off()
#		}
#	}

	## Plot observed and expected mean ages from commercial and survey C@A data.
	for (p in ptypes) {
		if (p=="eps") postscript("meanAge.eps", width=7, height=8.5, horizontal=FALSE,  paper="special")
		else if (p=="png") png("meanAge.png", units="in", res=pngres, width=7, height=8.5)
		plotMeanAge(obj=currentRes)
		if (is.element(p,c("eps","png"))) dev.off()
	}

	# Plot stock-recruitment function (based on MPD's)
	# xLimSR and yLimSR fixed here for YMR to have Run 26 and 27 figs
	#  on same scales. Use these first two values to scale to data:
	# xLimSR =c(0, max(obj$B$SB))
	# yLimSR=c(0, max(obj$B$R, na.rm=TRUE))
	#xLimSR=c(0, max(c(max(obj$B$SB),45000)))   # so it draw bigger if necessary
	#yLimSR=c(0, max(c(max(obj$B$R, na.rm=TRUE),55000)))
	xLimSR=c(0, 1.5*max(obj$B$SB,na.rm=TRUE))   # so it draw bigger if necessary
	xxx=(seq(0, xLimSR[2], length.out=100))
	yyy=srFun(xxx)
	yLimSR=c(0, 1.1*max(c(yyy,obj$B$R),na.rm=TRUE))

	for (p in ptypes) {
		if (p=="eps") postscript("stockRecruit.eps", width=6.5, height=4, horizontal=FALSE,  paper="special")
		else if (p=="png") png("stockRecruit.png", units="in", res=pngres, width=6.5, height=4)
		par(mfrow=c(1,1), mar=c(3.25,3.5,1,1), oma=c(0,0,0,0), mgp=c(2,0.75,0))
		plot(xxx, yyy, lwd=2, xlim=xLimSR, ylim=yLimSR, type="l",
			xlab=expression( paste("Spawning biomass ",  italic(B)[italic(t)-1], " (tonnes) in year ", italic(t), "-1", sep="") ),
			ylab=expression( paste("Recruitment ", italic(R)[italic(t)], " (1000s) in year ", italic(t), sep="") ) )
		text(obj$B[-length(years), "SB"], obj$B[-1, "R"], labels=substring(as.character(years), 3), cex=0.6, col="blue")
		if (p %in% c("eps","png")) dev.off()
	}

	##windows()
	##plt.lengthResids( stdRes.CL( obj$CLs ),
	##  main=paste("Survey",mainTitle,"Series",i) )
	##if ( save )
	##  savePlot( "surveyLengthResids", type="png" )
	## plt.idx( obj$CPUE,  main="Commercial Fishery",save="fishResids" )
	## plt.idx( obj$Survey,main="Survey",save="surveyResids" )
  closeAllWin()
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.mpdGraphs


#plt.mcmcGraphs-------------------------2014-09-23
# plt.mcmcGraphsAndySAR.r - Andy adding in wmf's for SAR. Based
#  on ymrScape.r. Just adding in the ones needed. 5th Feb 2012.
#-------------------------------------------------
# plt.mcmcGraphsFromPBSscape16oct12.r.   5th Feb 2012.
#-------------------------------------------------
#  AME editing (with *AME*) to give a policy option 
#  to be specified in run-master.Snw. 16th August 2012.
#-------------------------------------------------
# From PBSawatea from Rowan, dated 16th Oct 2012. Have checked
#  plt.mcmcGraphs() is identical to if I load PBSawatea and type
#  plt.mcmcGraphs (except formatting). Need to add in .wmf for
#  SAR - seemed to have that in ymrScape.r for YMR SAR, but the .wmf
#  commands are not in this one yet. So edit this, call it pltmcmcGraphsAndySAR.r
#  then send back to Rowan to put into package. May not have to edit any other files?
#-------------------------------------------------
# RH (2014-09-23)
#  Aded arguments `ptype',`pngres', and `ngear'
#  to render output figures in multiple formats
#  (only `eps' and `png' at the moment), and
#  to accommodate multiple gear types
#-------------------------------------------AME/AM
plt.mcmcGraphs <-
function (mcmcObj, projObj=NULL, mpdObj=NULL, save=FALSE, 
   ptypes=tcall(PBSawatea)$ptype, pngres=400, ngear=1,
   ylim.recruitsMCMC=NULL, ylim.exploitMCMC=NULL,
   ylim.VBcatch=NULL, ylim.BVBnorm=NULL,
   xlim.snail=NULL, ylim.snail=NULL,
   plotPolicies=names(projObj$Y[1:6]),
   onePolicy=names(projObj$Y[2]), mpd=list(),
   SAR.width=7.5, SAR.height=4, trevObj=NULL)
# plotPolicies is 6 policies projections to plot *AME*
# onePolicy is one to use for some figures *AME*
#*AME*xlim.pdfrec was =c(0, 200000). Put options for others
#  that will be useful if want to scale two model runs
#  to the same ylim. If NULL then fits ylim automatically.
{
	panel.cor <- function(x, y, digits=2, prefix="",...)
	{
		usr <- par("usr"); on.exit(par(usr))
		par(usr=c(0, 1, 0, 1))
		r <- abs(cor(x, y))
		txt <- format(c(r, 0.123456789), digits=digits)[1]
		txt <- paste(prefix, txt, sep="")
		text(0.5, 0.5, txt, cex=1.75)
	}
	panel.cor.small = eval(parse(text=sub("1\\.75", "1.25", deparse(panel.cor))))
#browser();return()

	for (p in ptypes) {
		if (p=="eps") postscript("recruitsMCMC.eps", width=6.25, height=4, horizontal=FALSE,  paper="special")
		else if (p=="png") png("recruitsMCMC.png", units="in", res=pngres, width=6.25, height=4)
		par(mfrow=c(1,1), mar=c(3,3,0.5,0.5), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
		plotRmcmcPOP(mcmcObj$R, yLim=ylim.recruitsMCMC) # *AME*
		if (p %in% c("eps","png")) dev.off()
	}

	for (p in ptypes) {
		if (p=="eps") postscript("exploitMCMC.eps", width=6.25, height=4*ngear, horizontal=FALSE,  paper="special")
		else if (p=="png") png("exploitMCMC.png", units="in", res=pngres, width=6.25, height=4*ngear)
		par(mfrow=c(ngear,1), mar=c(3,3,0.5,0.5), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
		for (g in 1:ngear) {
			gfile = mcmcObj$U[,grep(paste0("_",g),names(mcmcObj$U))]
			names(gfile) = substring(names(gfile),1,4)
			plotRmcmcPOP(gfile, yLab=paste0("Exploitation rate",ifelse(ngear>1,paste0(" - ",Cnames[g]),"")), yLim=ylim.exploitMCMC, yaxis.by=0.01)
		}
		if (p %in% c("eps","png")) dev.off()
	}

	for (p in ptypes) {
		if (p=="eps") postscript("pdfParameters.eps", width=6.25, height=7, horizontal=FALSE,  paper="special")
		else if (p=="png") png("pdfParameters.png", units="in", res=pngres, width=6.25, height=7)
		par(mfrow=c(1,1), mar=c(3,3,0.5,0.5), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
		plotDensPOPparsPrior(mcmcObj$P, lty.outer=2, between=list(x=0.3, y=0.2),mpd=mpd[["mpd.P"]])
		if (p %in% c("eps","png")) dev.off()
	}

	for (p in ptypes) {
		if (p=="eps") postscript("pdfBiomass%d.eps", width=6.5, height=8, horizontal=FALSE,  paper="special", onefile=FALSE)
		else if (p=="png") png("pdfBiomass%d.png", units="in", res=pngres, width=6.5, height=8)
		par(mfrow=c(1,1), mar=c(3,3,0.5,0.5), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
		plotDensPOP(mcmcObj$B/1000, xlab="Female spawning biomass, Bt (1000 t)", 
			between=list(x=0.2, y=0.2), ylab="Density", lwd.density=2, #panel.height=list(x=rep(1,5),unit="inches"), #*****Needs resolving
			same.limits=TRUE, layout=c(4,5), lty.outer=2, mpd=mpd[["mpd.B"]]/1000) 
		if (p %in% c("eps","png")) dev.off()
	}

	for (p in ptypes) {
		if (p=="eps") postscript("pdfRecruitment%d.eps", width=6.5, height=8, horizontal=FALSE,  paper="special", onefile=FALSE)
		else if (p=="png") png("pdfRecruitment%d.png", units="in", res=pngres, width=6.5, height=8)
		par(mfrow=c(1,1), mar=c(3,3,0.5,0.5), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
		plotDensPOP(mcmcObj$R/1000, xlab="Recruitment, Rt (1000s)", 
			between=list(x=0.2, y=0.2), ylab="Density", lwd.density=2,
			same.limits=TRUE, layout=c(4,5), lty.median=2, lty.outer=2, mpd=mpd[["mpd.R"]]/1000)
		if (p %in% c("eps","png")) dev.off()
	}

	for (p in ptypes) {
		if (p=="eps") postscript("pdfBiomass.eps", width=6.5, height=8, horizontal=FALSE,  paper="special")
		else if (p=="png") png("pdfBiomass.png", units="in", res=pngres, width=6.5, height=8)
		par(mfrow=c(1,1), mar=c(3,3,0.5,0.5), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
		plotDensPOP(mcmcObj$B[,getYrIdx(names(mcmcObj$B))]/1000, xlab="Female spawning biomass, Bt (1000 t)", 
			between=list(x=0.2, y=0.2), ylab="Density", lwd.density=2, #panel.height=list(x=rep(1,5),unit="inches"), #*****Needs resolving
			same.limits=TRUE, lty.outer=2, mpd=mpd[["mpd.B"]][getYrIdx(names(mcmcObj$B))]/1000) #, layout=c(0,length(getYrIdx(names(mcmcObj$B)))) )
		if (p %in% c("eps","png")) dev.off()
	}

	for (p in ptypes) {
		if (p=="eps") postscript("pdfRecruitment.eps", width=6.5, height=8, horizontal=FALSE,  paper="special")
		else if (p=="png") png("pdfRecruitment.png", units="in", res=pngres, width=6.5, height=8)
		par(mfrow=c(1,1), mar=c(3,3,0.5,0.5), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
		plotDensPOP(mcmcObj$R[,getYrIdx(names(mcmcObj$R))]/1000, xlab="Recruitment, Rt (1000s)", 
			between=list(x=0.2, y=0.2), ylab="Density", lwd.density=2,
			same.limits=TRUE, lty.median=2, lty.outer=2, mpd=mpd[["mpd.R"]][getYrIdx(names(mcmcObj$R))]/1000) #, layout=c(4,5))
		if (p %in% c("eps","png")) dev.off()
	}

	for (p in ptypes) {
		if (p=="eps") postscript("traceBiomass.eps", width=6.25, height=7, horizontal=FALSE,  paper="special")
		else if (p=="png") png("traceBiomass.png", units="in", res=pngres, width=6.25, height=7)
		par(mfrow=c(1,1), mar=c(3,3,0.5,0.5), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
		plotTracePOP(mcmcObj$B[,getYrIdx(names(mcmcObj$B))]/1000, axes=TRUE, between=list(x=0.2, y=0.2),
			xlab="Sample", ylab="Female spawning biomass, Bt (1000 t)", mpd=mpd[["mpd.B"]][getYrIdx(names(mcmcObj$B))]/1000)
		if (p %in% c("eps","png")) dev.off()
	}

	for (p in ptypes) {
		if (p=="eps") postscript("traceRecruits.eps", width=6.25, height=7, horizontal=FALSE,  paper="special")
		else if (p=="png") png("traceRecruits.png", units="in", res=pngres, width=6.25, height=7)
		par(mfrow=c(1,1), mar=c(3,3,0.5,0.5), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
		plotTracePOP(mcmcObj$R[,getYrIdx(names(mcmcObj$R))]/1000, axes=TRUE, between=list(x=0.2, y=0.2), 
			xlab="Sample", ylab="Recruitment, Rt (1000s)", mpd=mpd[["mpd.R"]][getYrIdx(names(mcmcObj$R))]/1000)
		if (p %in% c("eps","png")) dev.off()
	}

	for (p in ptypes) {
		if (p=="eps") postscript("traceParams.eps", width=6.25, height=7, horizontal=FALSE,  paper="special")
		else if (p=="png") png("traceParams.png", units="in", res=pngres, width=6.25, height=7)
		par(mfrow=c(1,1), mar=c(3,3,0.5,0.5), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
		idx <- apply(mcmcObj$P, 2, allEqual)
		plotTracePOP(mcmcObj$P[, !idx], axes=TRUE, between=list(x=0.2, y=0.2), 
			xlab="Sample", ylab="Parameter estimate", mpd=mpd[["mpd.P"]][!idx])
		if (p %in% c("eps","png")) dev.off()
	}

	for (p in ptypes) {
		if (p=="eps") postscript("splitChain.eps", width=6.25, height=7, horizontal=FALSE,  paper="special")
		else if (p=="png") png("splitChain.png", units="in", res=pngres, width=6.25, height=7)
		par(mfrow=c(1,1), mar=c(3,3,0.5,0.5), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
		#plotChains(mcmc=mcmcObj$P, axes=TRUE, pdisc=0, between=list(x=0, y=0), col.trace=c("green","red","blue"), xlab="Sample", ylab="Cumulative Frequency", cex.lab=1.5, yaxt="n")
		plotChains(mcmc=mcmcObj$P, axes=TRUE, pdisc=0, between=list(x=0, y=0), col.trace=c("red","green3","blue"), xlab="Parameter Value", ylab="Cumulative Frequency", cex.axis=1.3, cex.lab=1.4, yaxt="n")
		if (p %in% c("eps","png")) dev.off()
	}

	for (p in ptypes) {
		if (p=="eps") postscript("paramACFs.eps", width=8, height=8, horizontal=FALSE,  paper="special")
		else if (p=="png") png("paramACFs.png", width=8, height=8, units="in", res=pngres)
		plotACFs(currentMCMC, lag.max=60)
		if (p %in% c("eps","png")) dev.off()
	}

	for (p in ptypes) {
		if (p=="eps") postscript("VBcatch.eps", width=6.25, height=4.5*ngear, horizontal=FALSE,  paper="special")
		else if (p=="png") png("VBcatch.png", units="in", res=pngres, width=6.25, height=4.5*ngear)
		par(mfrow=c(ngear,1), mar=c(3,3,0.5,1), oma=c(0,ifelse(ngear>1,2,0),0,0), mgp=c(1.75,0.5,0))
		for (g in 1:ngear) {
			gfile = mcmcObj$VB[,grep(paste0("_",g),names(mcmcObj$VB))]
			names(gfile) = substring(names(gfile),1,4)
			plotVBcatch(gfile, mpdObj, gear=g, yLab=ifelse(ngear==1,"Catch and vulnerable biomass (t)",Cnames[g]), yLim=c(0,max(sapply(gfile,quantile,quants5[5]))),cex.lab=1.25)
			if (ngear>1) mtext("Catch and vulnerable biomass (t)",outer=TRUE,side=2,line=0.5,cex=1.5)
		}
		if (p %in% c("eps","png")) dev.off()
	}

#	win.metafile("SARVBcatch.wmf", width=0.8*SAR.width, height=0.8*SAR.height)
#	par(mai=c(0.72, 0.82, 0.2, 0.42), mgp=c(2,1,0))  
#	plotVBcatch( mcmcObj$VB, currentRes, yLim=ylim.VBcatch)
#	# mtext(SAR.main, side=3, font=1, cex=cex.main, line=line.main)
#	dev.off()
#
	for (p in ptypes) {
		if (p=="eps") postscript("BVBnorm.eps", width=6.25, height=5, horizontal=FALSE,  paper="special")
		else if (p=="png") png("BVBnorm.png", units="in", res=pngres, width=6.25, height=5)
		par(mfrow=c(1,1), mar=c(3,3,0.5,1), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
		plotBVBnorm(mcmcObj, xLeg=0.02, yLeg=0.2, yLim=ylim.BVBnorm, ngear=ngear, VB.col=c("blue","red"))
		if (p %in% c("eps","png")) dev.off()
	}

#    win.metafile("SARBVBnorm.wmf", height=SAR.height*0.7, width=SAR.width/2.3) # half width to do side-by-side
#    par(mai=c(0.62, 0.65, 0.05, 0.05), mgp=c(2,1,0))
#    plotBVBnorm(mcmcObj, xLeg=0.02, yLeg=0.25)
#    #, yLim=c(0, 1.1))  # yLim fixed for YMR11 submission, yLeg increased for SAR 
#    #  mtext(SAR.main, side=3, font=1, cex=cex.main, line=line.main)
#    dev.off()

#browser(); return()
	options(scipen=10)
	for (p in ptypes) {
		if (p=="eps") postscript("Bproj.eps", width=6.25, height=7, horizontal=FALSE, paper="special")
		else if (p=="png") png("Bproj.png", units="in", res=pngres, width=6.25, height=7)
		par(mfrow=c(1,1), mar=c(3,3,0.5,0.5), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
		plt.quantBio(mcmcObj$B, projObj$B, xyType="quantBox", policy=plotPolicies, save=FALSE)  # *AME* 
		if (p %in% c("eps","png")) dev.off()
	}

	for (p in ptypes) {
		if (p=="eps") postscript("Rproj.eps", width=6.25, height=7, horizontal=FALSE, paper="special")
		else if (p=="png") png("Rproj.png", units="in", res=pngres, width=6.25, height=7)
		par(mfrow=c(1,1), mar=c(3,3,0.5,0.5), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
		plt.quantBio(mcmcObj$R, projObj$R, xyType="quantBox", policy=plotPolicies, save=FALSE, yaxis.lab="Recruitment (1000s)")
		if (p %in% c("eps","png")) dev.off()
	}

	for (p in ptypes) {
		if (p=="eps") postscript("RprojOnePolicy.eps", width=6.25, height=5, horizontal=FALSE, paper="special")
		else if (p=="png") png("RprojOnePolicy.png", units="in", res=pngres, width=6.25, height=5)
		par(mfrow=c(1,1), mar=c(3,3,0.5,0.5), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
		plt.quantBioBB0(mcmcObj$R, projObj$R, xyType="quantBox", policy=onePolicy, 
			save=FALSE, xaxis.by=10, yaxis.lab="Recruitment (1000s)")    # *AME* (onePolicy)
		if (p %in% c("eps","png")) dev.off()
	}

	for (p in ptypes) {
		if (p=="eps") postscript("snail.eps", width=6.25, height=5, horizontal=FALSE, paper="special")
		else if (p=="png") png("snail.png", units="in", res=pngres, width=6.25, height=5)
		par(mfrow=c(1,1), mar=c(3,3.75,0.5,0.5), oma=c(0,0,0,0), mgp=c(2,0.5,0))
		plotSnail(mcmcObj$BoverBmsy, mcmcObj$UoverUmsy, p=quants3[c(1,3)], xLim=xlim.snail, yLim=ylim.snail, ngear=ngear, assYrs=2010) ## RSR in 5RF
		if (p %in% c("eps","png")) dev.off()
	}

	# Doing 6 pairs on a page:
	npr = 6
	nuP = length(use.Pnames)
	npp = ceiling(nuP/npr) # number of pairs plots
	for (i in 1:npp) {
		if (i<npp) ii = (1:npr)+(i-1)*npr
		else       ii = (nuP-npr+1):nuP
		for (p in ptypes) {
			pname = paste0("pairs",i)
			if (p=="eps") postscript(paste0(pname,".eps"), width=7, height=7, horizontal=FALSE,  paper="special")
			else if (p=="png") png(paste0(pname,".png"), units="in", res=pngres, width=7, height=7)
			par(mar=c(2,2,0.5,0.5), oma=c(0,0,0,0), mgp=c(2,0.75,0))
			pairs(mcmcObj$P[, ii], col="grey25", pch=20, cex=0.2, gap=0, lower.panel=panel.cor, cex.axis=1.5)
			if (p %in% c("eps","png")) dev.off()
		}
	}
	# Doing 1 pairs plot with all parameters
	npp = 1
	nuP = length(use.Pnames)
	npr = ceiling(nuP/npp)
	for (i in 1:npp) {
		if (i<npp) ii = (1:npr)+(i-1)*npr
		else       ii = (nuP-npr+1):nuP
		for (p in ptypes) {
			pname = "pairsPars"
			if (p=="eps") postscript(paste0(pname,".eps"), width=10, height=10, horizontal=FALSE,  paper="special")
			else if (p=="png") png(paste0(pname,".png"), units="in", res=pngres, width=10, height=10)
			par(mar=c(2,2,0.5,0.5), oma=c(0,0,0,0), mgp=c(2,0.75,0))
			pairs(mcmcObj$P[, ii], col="grey25", pch=20, cex=0.2, gap=0, lower.panel=panel.cor.small, cex.axis=1.25)
			if (p %in% c("eps","png")) dev.off()
		}
	}

	if (is.null(trevObj)) trevObj=trevorMCMC # created by Sweave code
	names(trevObj) = gsub("_","",names(trevObj))
	for (p in ptypes) {
		if (p=="eps") postscript("pairsMSY.eps", width=7, height=7, horizontal=FALSE,  paper="special")
		else if (p=="png") png("pairsMSY.png", units="in", res=pngres, width=7, height=7)
		par(mar=c(2,2,0.5,0.5), oma=c(0,0,0,0), mgp=c(2,0.75,0))
		pairs(trevObj, col="grey25", pch=20, cex=.2, gap=0, lower.panel=panel.cor, cex.axis=1.5)
		if (p %in% c("eps","png")) dev.off()
	}

	while(dev.cur() > 1)  dev.off()    # tidy up any remainingfrom the %d.eps
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.mcmcGraphs


load.allResFiles <- function( resList=NULL )  
{
	# Loads all Awatea "res" files in working directory into a list.
	if ( is.null(resList) )
		resList <- dir( path="../.", pattern="results.dat$" )
	#YMR AME path was "."
	# YMR pattern was ".res$", but need to go through Excel
	#  to get .res, and .dat has same numbers (see .Snw)
	# NOT called now by ymrrun1dos.Snw,just import results.dat

	result <- as.list( c(1:length(resList)) )
	names( result) <- resList

	for ( i in 1:length(result) )
		result[[i]] <- importRes( res.file=resList[i], Dev=TRUE, CPUE=TRUE,Survey=TRUE, CLc=TRUE, CLs=TRUE, CAs=TRUE, CAc=TRUE, extra=TRUE) # AME added CAc=TRUE, CAs=TRUE.
	result                             
}


#plotB2---------------------------------2011-08-31
# This is an alteration of Arni Magnussons "plotB" function to accommodate
# PJS request not to show biomass prior to fishery and survey indices period.
#----------------------------------------------AME
plotB2 <- function (model, what="d", series=NULL, years=NULL, axes=TRUE,
    div=1, legend="bottom", main="", xlab="", ylab="",
    cex.main=1.2, cex.legend=1, cex.lab=1, cex.axis=0.8,
    las=1, tck=c(1, what == "d")/2, tick.number=5, lty.grid=3,
    col.grid="white", pch=16, cex.points=0.8, col.points="black",
    lty.lines=1:3, lwd.lines=2, col.lines="black", ratio.bars=3,
    col.bars="grey", plot=TRUE, ...)
{
    panel.linebar <- function(x, y, bars, ...) {
        panel.abline(h=pretty(y, tick.number), lty=lty.grid,
            col=col.grid)
        panel.superpose(x, y, ...)
        panel.barchart(bars$Year, bars$Value, horizontal=FALSE,
            box.ratio=ratio.bars, col=col.bars)
    }
    panel.bar <- function(x, y, ...) {
        panel.abline(h=pretty(y, tick.number), lty=lty.grid,
            col=col.grid)
        panel.barchart(x, y, horizontal=FALSE, box.ratio=ratio.bars,
            col=col.bars)
    }
    if (class(model) != "scape")
        stop("The 'model' argument should be a scape object, not ",
            chartr(".", " ", class(model)), ".")
    what <- match.arg(what, c("d", "s", "l"))
    las <- as.numeric(las)
    x <- model$B
    x <- data.frame(Year=rep(x$Year, ncol(x) - 1), Series=rep(names(x)[-1],
        each=nrow(x)), Value=as.vector(as.matrix(x[, -1])))
    x$Value <- x$Value
    if (is.null(series))
        series <- unique(as.character(x$Series))
    if (is.null(years))
        years <- unique(x$Year)
    ok.series <- x$Series %in% series
    if (!any(ok.series))
        stop("Please check if the 'series' argument is correct.")
    ok.years <- x$Year %in% years
    if (!any(ok.years))
        stop("Please check if the 'years' argument is correct.")
    x <- x[ok.series & ok.years, ]


    Bframe <- x[x$Series %in% grep("B", series, value=TRUE),]
    Bframe$Series <- factor(Bframe$Series)

	## Find the first year where there are fishery CPUE data.
    cpueYear1 <- min( model$CPUE$Year[ !is.na( model$CPUE$Obs ) ] )

	## Set all SB and VB values to NA for years less than cpueYear1.
    Bframe$Value[ Bframe$Year < cpueYear1 ] <- NA

    Rframe <- x[x$Series == "R", ]
    Yframe <- x[x$Series == "Y", ]

    Bframe$Value <- Bframe$Value/div[1]
    Rframe$Value <- Rframe$Value/rep(div, length.out=2)[2]
    Yframe$Value <- Yframe$Value/div[1]

    #mess = c(
    #"require(grid, quietly=TRUE, warn.conflicts=FALSE)",
    #"require(lattice, quietly=TRUE, warn.conflicts=FALSE)"
    #)
    #eval(parse(text=mess))
    if (trellis.par.get()$background$col == "#909090") {
        for (d in dev.list()) dev.off()
        trellis.device(color=FALSE)
    }
    main <- rep(main, length.out=2)
    xlab <- rep(xlab, length.out=2)
    ylab <- rep(ylab, length.out=2)
    las <- rep(las, length.out=2)
    mymain <- list(label=main[1], cex=cex.main)
    myxlab <- list(label=xlab[1], cex=cex.lab)
    myylab <- list(label=ylab[1], cex=cex.lab)
    myrot <- switch(as.character(las[1]), "0"=list(x=list(rot=0),
        y=list(rot=90)), "1"=list(x=list(rot=0), y=list(rot=0)),
        "2"=list(x=list(rot=90), y=list(rot=0)), "3"=list(x=list(rot=90),
            y=list(rot=90)))
    myscales <- c(list(draw=axes, cex=cex.axis, tck=tck,
        tick.number=tick.number), myrot)
    lty.lines <- rep(lty.lines, length.out=nlevels(Bframe$Series))
    lwd.lines <- rep(lwd.lines, length.out=nlevels(Bframe$Series))
    col.lines <- rep(col.lines, length.out=nlevels(Bframe$Series))
    mykey <- list(space=legend, text=list(lab=levels(Bframe$Series),
        cex=cex.legend), lines=list(lty=lty.lines, lwd=lwd.lines,
        col=col.lines))
    if (what == "s") {
        graph <- xyplot(Rframe$Value ~ Bframe$Value[Bframe$Series ==
            "SB"], main=mymain, xlab=myxlab, ylab=myylab,
            scales=myscales, pch=pch, cex=cex.points, col=col.points,
            ...)
        graph$x.limits[1] <- 0
    }
    else if (what == "d" && nrow(Bframe) > 0) {
        graph <- xyplot(Value ~ Year, groups=Series, data=Bframe,
            panel=panel.linebar, type="l", bars=Yframe,
            main=mymain, xlab=myxlab, ylab=myylab, scales=myscales,
            key=mykey, lty=lty.lines, lwd=lwd.lines, col=col.lines,
            ...)
    }
    else {
        graph <- xyplot(Value ~ Year, data=Yframe, panel=panel.bar,
            main=mymain, xlab=myxlab, ylab=myylab, scales=myscales,
            ...)
    }
    graph$y.limits[1] <- 0
    if (plot) {
        print(graph)
        invisible(x)
    }
    else {
        invisible(graph)
    }
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotB2


#plotBmcmcPOP---------------------------2011-08-31
# AME writing plotBmcmcPOP(), to plot spawning biomass and vulnerable biomass
#  from posterior, so do as boxplot, and also the catch on same graph. So
#  combining some of plt.quantBio and plotB2. Don't need lattice, just one
#  figure, no panels..     Vulnerable biomass has no posterior saved, which
#  must be why it's not been done before. Hmmm.... still worth seeing spawning
#  though?
# Taking what's needed from plt.quantBio, this basically works:
#        plt.quantBio( currentMCMC$B, xyType=rpType ), though does 2x3 plots
# obj should be the specficic MCMC posterior by year (so just a data.frame),
#  e.g. currentMCMC$B.  currentRes1 is local currentRes.
#----------------------------------------------AME

#plotVBcatch----------------------------2014-09-22
# AME adding, based on plotBmcmcPOP (tweaking some)
#  currentMCMC$B.  currentRes1 is local currentRes.
#-------------------------------------------AME/RH
plotVBcatch=function(obj, currentRes1=currentRes,
   p = get("quants5"),
   xyType="quantBox",
   lineType=c(3,2,1,2,3),
   refLines=NULL, xLim=NULL, yLim=NULL,
   userPrompt=FALSE, save=TRUE, 
   xLab="Year",
   yLab="Catch and vulnerable biomass (t)",
   textLab=c("catch", "vulnerable"),
   yaxis.by=10000, tcl.val=-0.2,
   gear=1, ...)
   # xLab - x position for label, etc.
{
  # See plt.quantBio if want other xyTypes, as took out here:
  plt.qB <- function( obj, xyType="lines", new=TRUE, xLim, yLim, ... ) 
  {
    if ( new )
     plot( xLim,yLim, type="n", xlab=xLab, ylab=yLab, ... )
    
    yrs <- as.numeric(dimnames(obj)[[2]])

    # Quantile boxplots - assumes five quantiles.
    if ( xyType=="quantBox" )
    {
      delta <- 0.25
      # Draw the outer whiskers.
      segments( yrs,obj[1,], yrs,obj[5,], lty=1,col=1 )
      # Overlay the box.
      for ( i in 1:length(yrs) )
        rect( yrs[i]-delta,obj[2,i], yrs[i]+delta, obj[4,i], ... )
      # Add the median.
      segments( yrs-delta,obj[3,],yrs+delta,obj[3,],lty=1,col=1 )
    }
  }

  # Plot quantiles of biomass using the posterior densities.

  yrs1 <- NULL
  yrs2 <- NULL
  result1 <- NULL
  result2 <- NULL

  # Calculate the quantiles of the reconstructed biomass.
  result1 <- apply( obj,2,quantile,probs=p )
  yrs1 <- as.numeric(dimnames(result1)[[2]])

  if ( is.null(yLim) )
    {
      yLim <- c(0, max(c(max(result1), max(currentRes1$B$VB)))) #range(result1)
    }
  if ( is.null(xLim) )
    {
       xLim=range(yrs1)
    }

  # xLegPos=xLeg*diff(xLim)+xLim[1]    # position of xLeg
  # yLegPos=yLeg*diff(yLim)+yLim[1]

  plt.qB( result1,xLim=xLim,yLim=yLim, xyType=xyType, yaxt="n", ...)
  #points(currentRes1$B$Year, currentRes1$B$Y, type="h", lwd=3)   # catch -- RH: won't work when Ngear > 1
  Cgears  = currentRes1$B[,-1][,grep("Y",names(currentRes1$B[,-1])),drop=FALSE]
  Ctotal = apply(Cgears,1,sum)
  zpos = Cgears[,gear]>0
  points(currentRes1$B$Year[zpos], Cgears[,gear][zpos], type="h", lwd=3)   # gear-specific catch -- RH: use in case Ngear > 1
  # points(obj$B$Year, currentRes1$B$VB, type="p")
                          # was vuln biom MPD
  # text( xLab, yLab, textLab, pos=4, offset=0)   # Taking out
  #   as not really needed if give a decent caption
  axis(1, at=intersect(seq(1900,3000,5),xLim[1]:xLim[2]), tcl=tcl.val, labels=FALSE)
  axis(2, at=seq(0, yLim[2], by=yaxis.by), tcl=tcl.val, labels=FALSE)
  axis(2, at=pretty(yLim), labels=format(pretty(yLim),scientific=FALSE,big.mark=","))

  # legend(xLegPos, yLegPos, c("Vulnerable", "Spawning", "Catch"), bty="n")
  # points(xLegPos-2, yLegPos, type="p")
    # mtext( side=1, line=2, cex=1.0, "Year" )
    # mtext( side=2, line=2, cex=1.0, "Biomass" )
  # }
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotVBcatch
# RH -- following 3 lines for debugging only
#g=1; ngear=1
#test = currentMCMC$VB[,grep(paste0("_",g),names(currentMCMC$VB))]; names(test) = substring(names(test),1,4)
#plotVBcatch(test, currentRes, gear=g, yLab=ifelse(ngear==1,"Catch and vulnerable biomass (t)",Cnames[g]), yLim=c(0,max(sapply(test,quantile,quants5[5]))),cex.lab=1.25)

## plotBVBnorm--------------------------2018-05-16
## AME doing, tried in separate file, but then changed that to
##  lattice and wouldn't be good format for Arni's boxplots.
##  Based on plotVBcatch (tweaking some)
##  currentMCMC$B.  currentRes1 is local currentRes.
##  xLab - x position for label, etc.
## -----------------------------------------AME/RH
plotBVBnorm=function(mcmcObj,
   p = get("quants5"),
   xyType="quantBox",
   lineType=c(3,2,1,2,3),
   refLines=NULL, xLim=NULL, yLim=NULL,
   userPrompt=FALSE, save=TRUE, xLeg=0.05, yLeg=0.2,  # legend in relative (0,1) space
   yaxis.by=0.05, tcl.val=-0.2,
   B.col="black", VB.col="black", ngear=1, ...)
   # xLab - x position for label, etc.
{
	BVBlist = as.list(0:ngear); names(BVBlist)=c("Spawning Biomass",Cnames[1:ngear])
	BVBlist[[1]] = mcmcObj$B
	for (g in 1:ngear) {
		gfile = mcmcObj$VB[,grep(paste0("_",g),names(mcmcObj$VB))]
		names(gfile) = substring(names(gfile),1,4)
		BVBlist[[g+1]] = gfile
	}
   # Calculate medians to be plotted
	BVB0list = sapply(BVBlist,function(x){x/x[,1]},simplify=FALSE)  #B/B0 and VB/VB0 each chain
	BVB0med  = sapply(BVB0list,function(x){apply(x,2,median)},simplify=FALSE)  # median each year

	# Plot quantiles of biomass using the posterior densities.
	if ( is.null(yLim) )
		yLim = c(0,max(sapply(BVB0med,max)))
	yrs = sapply(BVB0med,function(x){as.numeric(names(x))},simplify=FALSE)
	if ( is.null(xLim) )
		xLim = range(sapply(yrs,range))
	VB.col = rep(VB.col,ngear)[1:ngear]

	plot(0, 0, xlim=xLim, ylim=yLim, type="n", xlab="Year",ylab="Biomass relative to unfished equilibrium",cex.lab=1.25)
	axis(1, at=intersect(seq(1900,3000,5),xLim[1]:xLim[2]), tcl=tcl.val, labels=FALSE)
	axis(2, at=seq(0, yLim[2], by=yaxis.by), tcl=tcl.val, labels=FALSE)
	points(yrs[[1]], BVB0med[[1]], type="p", col=B.col) 
	for (i in 1:ngear)
		points(yrs[[i+1]], BVB0med[[i+1]], type="l", col=VB.col[i])
	legtxt = bquote(paste(italic(B[t])/italic(B)[0], "    Spawning"))
	for (i in 1:ngear)
		legtxt = c(legtxt, bquote(paste(italic(V[t])/italic(V)[0], "    Vulnerable - ", .(tolower(Cnames[i])) ) ) )
	addLegend(xLeg,yLeg,legend=as.expression(legtxt),bty="n", lty=c(0,rep(1,ngear)), pch=c(1, rep(NA,ngear)), col=c(B.col, VB.col))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotBVBnorm


#plotRmcmcPOP---------------------------2014-09-19
# AME adding, plotting recruitment posteriors quantiles as one graph over time.
#  Already have the full posterior densities done.
#  Using plotBmcmcPOP as template, but will be simpler as no extra stuff. Prob
#   not simplifying down as much as could due to time constraints.
# Adding yLab and then using for exploitation plot also
#----------------------------------------------AME
plotRmcmcPOP=function(obj, 
   p = get("quants5"),
   xyType="quantBox",
   lineType=c(3,2,1,2,3),
   refLines=NULL, xLim=NULL, yLim=NULL,
   userPrompt=FALSE, save=TRUE, tcl.val=-0.2,
   yaxis.by=10000, yLab="Recruitment", ...)
{
  # See plt.quantBio if want other xyTypes, as took out here:
  plt.qB <- function( obj, xyType="lines", new=TRUE, xLim, yLim,... ) 
    {
    if ( new )
     plot( xLim,yLim, type="n", xlab="Year",ylab=yLab )

    yrs <- as.numeric(dimnames(obj)[[2]])

    # Quantile boxplots - assumes five quantiles.
    if ( xyType=="quantBox" )
    {
      delta <- 0.25
      # Draw the outer whiskers.
      segments( yrs,obj[1,], yrs,obj[5,], lty=1,col=1 )
      # Overlay the box.
      for ( i in 1:length(yrs) )
        rect( yrs[i]-delta,obj[2,i], yrs[i]+delta, obj[4,i],... )
      # Add the median.
      segments( yrs-delta,obj[3,],yrs+delta,obj[3,],lty=1,col=1 )
    }
  }
  # Plot quantiles of biomass using the posterior densities.
	yrs1 = yrs2 = result1 = result2 = NULL

  # Calculate the quantiles of the reconstructed biomass.
  result1 <- apply( obj,2,quantile,probs=p )
  yrs1 <- as.numeric(dimnames(result1)[[2]])

  if ( is.null(yLim) )
    {
      yLim <- range(result1)
    }
    if ( is.null(xLim) )
    {
      xLim=range(yrs1)
    }
  plt.qB( result1,xLim=xLim,yLim=yLim, xyType=xyType )
  axis(1, at=intersect(seq(1900,3000,5),xLim[1]:xLim[2]), tcl=tcl.val, labels=FALSE)
  axis(2, at=seq(0, yLim[2], by=yaxis.by), tcl=tcl.val, labels=FALSE)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotRmcmcPOP

#--------------------------------------------------------------------#
#                        Utility Functions                           #
#--------------------------------------------------------------------#
allEqual <- function(x)
{
  result <- all( x==x[1] )
  result
}

closeAllWin <- function()
{
  winList <- dev.list()
  if ( !is.null(winList) )
    for ( i in 1:length(winList) )
      dev.off()
}

graphics <- function( view="portrait" )
{
	if ( view=="portrait" )
		do.call("windows", list(width=8.5, height=11, record=TRUE))
	if ( view=="landscape" )
		do.call("windows", list(width=11, height=8.5, record=TRUE))
}

panLab <- function( x, y, txt, ... )
{
#  orig.par <- par()
  usr <- par( "usr" )
  par( usr=c(0,1,0,1) )
  text( x, y, txt, ... )
  par( usr=usr )
#  par( orig.par )
  return( NULL )
}

panLegend <- function( x, y, legTxt, ... )
{
#  orig.par <- par()
  usr <- par( "usr" )
  par( usr=c(0,1,0,1) )
  legend( x, y, legend=legTxt, ... )
  par( usr=usr )
#  par( orig.par )
  return( NULL )
}

#--------------------------------------------------------------------#
#                        Calculation Functions                       #
#--------------------------------------------------------------------#

#calc.projExpect------------------------2011-08-31
# Calculate the expectation of projection to reference.
# Compare refYears to projection years.
#----------------------------------------------AME
calc.projExpect <- function( obj, projObj, refYrs )
{
  policyList <- names(projObj)
  projYrs <- dimnames( projObj[[1]] )[[2]]
  refYears <- as.character(refYrs)

  nPolicies <- length(policyList)
  nProjYrs <- length(projYrs)
  nRefYrs <- length(refYears)

  # Final results.
  result <- as.list( c(1:nRefYrs) )
  names( result ) <- refYears

  # Loop over the reference years.
  for ( i in 1:nRefYrs )
  {
    # Create a results matrix for nPolicies (j), nProjYr (k).
    val <- matrix( NA, nrow=nPolicies,ncol=nProjYrs )
    dimnames( val ) <- list( policyList,projYrs )

    # Loop over the policies.
    for ( j in 1:nPolicies )
    {
      # This are the projection results for policy j.
      proj <- projObj[[j]]

      # Loop over the projection years.
      for ( k in 1:ncol(proj) )
      {
        tmp <- proj[,k] / obj[,refYears[i]]
        val[j,k] <- mean( tmp )
      }
    }
    result[[i]] <- val
  }
  cat( "\n\nExpectation of projection biomass > reference year:\n\n" )
  print( result )
  result
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~calc.projExpect


#calc.projExpect2-----------------------2011-08-31
# Calculate expectation (projection biomass / reference biomass).
#----------------------------------------------AME
calc.projExpect2 <- function( obj, projObj, refList )
{
  policyList <- names(projObj)
  nPolicies <- length(policyList)

  projYrs <- dimnames( projObj[[1]] )[[2]]
  nProjYrs <- length(projYrs)

  refYears <- refList$refYrs
  nRefYrs <- nrow( refYears )
  funs <- refList$funVec

  # Final results, each reference period stored as list element.
  result <- as.list( c(1:nRefYrs) )
  names( result ) <- dimnames(refYears)[[1]]

  # Loop over the reference years.
  for ( i in 1:nRefYrs )
  {
    # Create a results matrix for nPolicies (j), nProjYr (k).
    val <- matrix( NA, nrow=nPolicies,ncol=nProjYrs )
    dimnames( val ) <- list( policyList,projYrs )

    # Build reference years and coerce to character for indexing.
    period <- as.character( c( refYears[i,1]:refYears[i,2] ) )

    # Calculate the reference value for the performance measure.
    refVal <- calc.refVal( obj,period,funs[[i]] )

    # Loop over the catch level policies.
    for ( j in 1:nPolicies )
    {
      # These are the projection results for policy j.
      proj <- projObj[[j]]

      # Loop over the projection years.
      for ( k in 1:ncol(proj) )
      {
        tmp <- proj[,k] / refVal
        val[j,k] <- mean( tmp )
      }
    }
    result[[i]] <- val
  }
  cat( "\n\nExpectation of (projection biomass / reference period):\n\n" )
  print( data.frame( refYears,funs ) )
  cat( "\n" )
  print( result )
  result$refs <- refList
  result
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~calc.projExpect2


#calc.projProbs-------------------------2011-08-31
# Calculate the probability of being greater than refYears.
# Compare refYears to projection years.
#----------------------------------------------AME
calc.projProbs <- function( obj, projObj, refYrs )
{
  policyList <- names(projObj)
  projYrs <- dimnames( projObj[[1]] )[[2]]
  refYears <- as.character(refYrs)

  nPolicies <- length(policyList)
  nProjYrs <- length(projYrs)
  nRefYrs <- length(refYears)

  # Final results.
  result <- as.list( c(1:nRefYrs) )
  names( result ) <- refYears

  # Loop over the reference years.
  for ( i in 1:nRefYrs )
  {
    # Create a results matrix for nPolicies (j), nProjYr (k).
    val <- matrix( NA, nrow=nPolicies,ncol=nProjYrs )
    dimnames( val ) <- list( policyList,projYrs )

    # Loop over the policies.
    for ( j in 1:nPolicies )
    {
      # This are the projection results for policy j.
      proj <- projObj[[j]]

      # Loop over the projection years.
      for ( k in 1:ncol(proj) )
      {
        tmp <- proj[,k] > obj[,refYears[i]]
        val[j,k] <- sum( tmp ) / length(tmp)
      }
    }
    result[[i]] <- val
  }
  cat( "\n\nProbability of projection biomass > reference year:\n\n" )
  print( result )
  result
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~calc.projProbs


#calc.projProbs2------------------------2011-08-31
# Calculate the probability of being greater than refYears.
# Compare refYears to projection years.
#----------------------------------------------AME
calc.projProbs2 <- function( obj, projObj, refList )
{
  policyList <- names(projObj)
  nPolicies <- length(policyList)

  projYrs <- dimnames( projObj[[1]] )[[2]]
  nProjYrs <- length(projYrs)

  refYears <- refList$refYrs
  nRefYrs <- nrow( refYears )
  funs <- refList$funVec

  # Final results, each reference period stored as list element.
  result <- as.list( c(1:nRefYrs) )
  names( result ) <- dimnames(refYears)[[1]]

  # Loop over the reference years.
  for ( i in 1:nRefYrs )
  {
    # Create a results matrix for nPolicies (j), nProjYr (k).
    val <- matrix( NA, nrow=nPolicies,ncol=nProjYrs )
    dimnames( val ) <- list( policyList,projYrs )

    # Build reference years and coerce to character for indexing.
    period <- as.character( c( refYears[i,1]:refYears[i,2] ) )

    # Calculate the reference value for the performance measure.
    refVal <- calc.refVal( obj,period,funs[[i]] )

    # Loop over the catch level policies.
    for ( j in 1:nPolicies )
    {
      # These are the projection results for policy j.
      proj <- projObj[[j]]

      # Loop over the projection years.
      for ( k in 1:ncol(proj) )
      {
         tmp <- proj[,k] > refVal
         val[j,k] <- sum( tmp ) / length(tmp)
      }
    }
    result[[i]] <- val
  }
  cat( "\n\nProbability of projection biomass > reference period:\n\n" )
  print( data.frame( refYears,funs ) )
  cat( "\n" )
  print( result )
  result$refs <- refList
  result
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~calc.projProbs2


#calc.refProbs--------------------------2011-08-31
# Calculate the reference probabilities (basing on calc.projProbs2)
#----------------------------------------------AME
calc.refProbs <- function( projObj=currentProj$B, refPlist=refPointsList )  
{
  # refPlist will be a list, LRP, URP and Bmsy are defaults with values for each draw
  policyList <- names(projObj)
  nPolicies <- length(policyList)

  projYrs <- dimnames( projObj[[1]] )[[2]]
  nProjYrs <- length(projYrs)
  
  nRefPoints=length(names(refPlist))  # refYears <- refList$refYrs
    # nRefYrs <- nrow( refYears )
    # funs <- refList$funVec
  # browser(); return()
  # Final results, each reference point stored as list element.
  result <- as.list(1:nRefPoints)   # to get right number of list elements (1,2,3, but
                                    #  they'll get overwritten)
  names( result ) <- names(refPlist)

  # Loop over the reference points.
  for ( i in 1:nRefPoints )
  {
    # Create a results matrix for nPolicies (j), nProjYr (k).
    val <- matrix( NA, nrow=nPolicies,ncol=nProjYrs )
    dimnames( val ) <- list( policyList,projYrs )

    # Build reference years and coerce to character for indexing.
    # period <- as.character( c( refYears[i,1]:refYears[i,2] ) )

    # Calculate the reference value for the performance measure.
    # refVal <- calc.refVal( obj,period,funs[[i]] )
    refVal=refPlist[[i]]
    
    # Loop over the catch level policies.
    for ( j in 1:nPolicies )
    {
      # These are the projection results for policy j.
      proj <- projObj[[j]]

      # Loop over the projection years.
      for ( k in 1:ncol(proj) )
      {
         tmp <- proj[,k] > refVal
         val[j,k] <- sum( tmp ) / length(tmp)
      }
    }
    result[[i]] <- val
  }
  cat( "\n\nProbability of projection biomass > reference point:\n\n" )
  # print( data.frame( refPointsrefYears,funs ) )
  # cat( "\n" )
  print( result )
  # result$refs <- refPlist
  result
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~calc.refProbs


#calc.refProbsHist----------------------2013-09-23
# Calculate the reference probabilities (basing on calc.projProbs2)
#-----------------------------------------------RH
calc.refProbsHist <- function( projObj=currentProj$B, refPlist=refPointsHistList[c("blimHRP","btarHRP")] )  
{
	# refPlist will be a list, LRP, URP and Bmsy are defaults with values for each draw
	policyList=names(projObj)
	nPolicies =length(policyList)
	projYrs   =dimnames( projObj[[1]] )[[2]]
	nProjYrs  =length(projYrs)
	nRefPoints=length(names(refPlist))  # refYears <- refList$refYrs
	# Final results, each reference point stored as list element.
	result    =sapply(as.list(1:nRefPoints),function(x){NULL},simplify=FALSE) # to get right number of list elements, but they'll get overwritten
	names( result ) <- names(refPlist)
	# Loop over the reference points.
	for ( i in names(refPlist) )
	{
		# Retrieve the reference value for the performance measure.
		refVal=refPlist[[i]]
		if (is.null(refVal))  next
		# Create a results matrix for nPolicies (j), nProjYr (k).
		val <- matrix( NA, nrow=nPolicies,ncol=nProjYrs )
		dimnames( val ) <- list( policyList,projYrs )
		# Loop over the catch level policies.
		for ( j in 1:nPolicies )
		{
			# These are the projection results for policy j.
			proj <- projObj[[j]]
			# Loop over the projection years.
			for ( k in 1:ncol(proj) )
			{
				tmp <- proj[,k] > refVal
				val[j,k] <- sum( tmp ) / length(tmp)
			}
		}
		result[[i]] <- val
	}
	cat( "\n\nProbability of projection biomass > reference point:\n\n" )
	print( result )
	return(result)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~calc.refProbsHist


#calc.refVal----------------------------2011-08-31
# Calculates the reference value for performance measures.
# Returns a vector of length nrow(obj) reference values.
#----------------------------------------------AME
calc.refVal <- function( obj, refYrs, fun=mean )
{
  # Input:
  # obj: scape Biomass matrix with n rows and nYr columns.
  # refYrs: numeric years in reference period.
  # fun: The function to apply to reference period i.

  # Coerce obj to a matrix.
  obj <- as.matrix( obj )

  # Change reference years to character to allow named indexing.
  refYrs <- as.character( refYrs )

  # Extract relevant columns and apply function to rows.
  # Coerce to a matrix to accommodate single year reference.

  tmp <- matrix( obj[,refYrs],ncol=length(refYrs) )
  result <- apply( tmp,1,fun )
  result
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~calc.refVal


#getYrIdx-------------------------------2011-08-31
# Purpose is to return selected years for plotting.
# Default is to select 5 year increments.
#----------------------------------------------AME
getYrIdx <- function( yrNames,mod=5 )
{
  # Coerce to numerice and select the years modulo "mod".
  yrVals <- as.numeric( yrNames )
  idx <- yrVals %% mod==0

  # Select years from character vector yrNames.
  result <- yrNames[ idx ]
  result
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~getYrIdx


#out.pmTables---------------------------2011-08-31
# Write the decision tables to a comma delimited file.
#----------------------------------------------AME
out.pmTables <- function( obj, fileName="pm", dec=3 )
{
  nPM <- length( obj )
  for ( i in 1:(nPM-1) )
  {
    tmp <- obj[[i]]
    tmp <- format( round(tmp,digits=dec) )

    rowNames <- dimnames(tmp)[[1]]
    colNames <- dimnames(tmp)[[2]]
    catch <- as.numeric( rowNames )
    dimnames( tmp ) <- list( NULL,NULL )
    tmp <- cbind( catch, tmp )
    dimnames( tmp ) <- list( NULL, c("Catch",colNames) )
    write.table( tmp, quote=FALSE, row.names=FALSE,
      file=paste( fileName,i,".csv",sep=""),sep="," )
  }
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~out.pmTables


#stdRes.CA------------------------------2011-08-31
# This implements standardised residuals for the Awatea
# implementation of the Fournier robustified normal likelihood
# for proportions at length. Based on PJS summary of CASAL doc and ACH change to length.
#----------------------------------------------AME
stdRes.CA <- function( obj, trunc=3, myLab="Age Residuals", prt=TRUE )
{
  # Make a column for the standardised residuals.
  result <- cbind( obj,stdRes=rep(NA,nrow(obj)) )

  # Raw residuals.
  res <- obj$Obs - obj$Fit

  # Number of age bins.
  # QUESTION: Should this be from FIRST age to plus group?
  n <- length( unique(obj$Age) )

  # Kludgy, but loop through years.
  # Could reformat into matrix and vectorize.

  yrList <- sort( unique( obj$Year ) )
  for ( i in 1:length( yrList ) )
  {
    idx <- yrList[i]==obj$Year
    Z <- obj$Obs[idx]*(1.0-obj$Obs[idx]) + 0.1/n
    N <- min( obj$SS[idx],1000)
    SD <- sqrt( Z/N )
    result$stdRes[idx] <- res[idx]/SD
  }

  if ( prt )
  {
    sdRes <- sqrt( var( result$stdRes,na.rm=TRUE ) )
    sdTrunc <- ifelse( result$stdRes > trunc, trunc, result$stdRes )
    sdTrunc <- ifelse( result$stdRes < -trunc, -trunc, result$stdRes )
    sdResTrunc <- sqrt( var( sdTrunc,na.rm=TRUE ) )
    cat( "\n",myLab,"\n" )
    cat( "\n     Std. Dev. of standardised Residuals=",sdRes,"\n" )
    cat( "     Std. Dev. of Truncated standardised Residuals=",sdResTrunc,"\n" )
  }
  result
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~stdRes.CA


#stdRes.index---------------------------2011-08-31
# Compute the standardised residuals.
#----------------------------------------------AME
stdRes.index <- function( obj, label=NULL, prt=TRUE )
{
  res <- log(obj$Obs) - log(obj$Fit)
  stdRes <- res / obj$CV
  result <- cbind( obj,stdRes )

  if ( prt )
  {
    sdRes <- sqrt( var( stdRes,na.rm=TRUE ) )
    cat( "\n",label,"\n" )
    cat( "\n     Std. Dev. of standardised Residuals=",sdRes,"\n" )
  }
  result
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~stdRes.index


#--------------------------------------------------------------------#
#                         Plotting Functions                         #
#--------------------------------------------------------------------#

# AME introducing for YMR. Bubble plots of data:
# plt.bubbles=function( obj,
# For now just do in Sweave as quicker.

#plt.ageResidsPOP-----------------------2011-08-31
# AME changing for POP, just plotting age class resids
#  here, not qq-plot. Moving that to new function (so can do 4 on
#  page). See popScape.r for original.
#----------------------------------------------AME
plt.ageResidsPOP <- function( obj, ages=c(2,60), pct=c(5,25,50,75,95)
                             ,  main=NULL )
{
  # Input is the output from stdRes.CA
  # par( oma=c(2,1,1,1), mar=c(2,2,2,1), mfrow=c(2,1) )
  # Subset to required ages.
  obj <- obj[ (obj$Age >= ages[1]) & (obj$Age <=ages[2]), ]
  if( max(diff(sort(unique(obj$Age)))) > 1)
    {
      allAges=min(obj$Age):max(obj$Age)
      nodataAges=allAges[ !(allAges %in% obj$Age)]
      xx=split(c(obj$stdRes, rep(NA, length(nodataAges))), c(obj$Age, nodataAges))
      xpos <- boxplot( xx, whisklty=1, xlab="", ylab="", outline=FALSE )     #AME outline=FALSE removes outliers
    } else
    {            
    xpos <- boxplot( split( obj$stdRes, obj$Age ), whisklty=1, xlab="", ylab="", outline=FALSE )     #AME outline=FALSE removes outliers
    }
  abline( h=0, lty=2, col="red" )
  mtext( side=1, line=2, cex=0.8, "Age class" )
  # These wouldn't show up:
  # legend(x="topleft", "(a)", bty="n")
  # text( par("usr")[1] + 0.95*diff(par("usr")[1:2]),
  #     par("usr")[3] + 0.05*diff(par("usr")[3:4]), "(a)") 
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.ageResidsPOP


#plt.ageResidsqqPOP---------------------2011-08-31
# Plotting qq plot for age class resids.
#----------------------------------------------AME
plt.ageResidsqqPOP <- function( obj, ages=c(2,60),
                pct=c(5,25,50,75,95),  main=NULL )
{
  # Subset to required ages.
  obj <- obj[ (obj$Age >= ages[1]) & (obj$Age <=ages[2]), ]

  # Plot the q-normal plot of the standardised residuals 
  qqnorm( obj$stdRes,xlab="",ylab="",main="" )
  abline( a=0,b=1 )
  abline( h=quantile(obj$stdRes,p=pct/100,na.rm=TRUE),lty=2 )
  mtext( side=1, line=2, cex=0.8, "Theoretical quantiles" )

  #mtext( side=2, line=-1, cex=0.8, outer=TRUE, "Standardised Residuals" )
  if ( !is.null(main) )
    mtext( side=3, line=-0.5, cex=1.0, outer=TRUE, main )
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.ageResidsqqPOP


#plt.yearResidsPOP----------------------2011-08-31
# AME adding to plot age residuals by year. Is called for comm and survs.
#  fill.in=TRUE is to add the missing years for boxplot
#  ..POP does not do qq plot. See popScape.r for previous.
#----------------------------------------------AME
plt.yearResidsPOP <- function( obj, ages=c(2,60), pct=c(5,25,50,75,95),
                           main=NULL, fill.in=TRUE, ... )
{
  # Subset to required ages - still do as don't want age 1.
  obj <- obj[ (obj$Age >= ages[1]) & (obj$Age <=ages[2]), ]
  if(fill.in)
    {
      allYears=min(obj$Year):max(obj$Year)
      nodataYears=allYears[ !(allYears %in% obj$Year)]
      xx=split(c(obj$stdRes, rep(NA, length(nodataYears))), c(obj$Year, nodataYears))
      xpos <- boxplot( xx, whisklty=1, xlab="", ylab="",
          outline=FALSE, ... )     #AME outline=FALSE removes outliers
      # browser()
    } else
    {  
      xpos <- boxplot( split( obj$stdRes, obj$Year ), whisklty=1, xlab="", ylab="", outline=FALSE, ... )     #AME outline=FALSE removes outliers
    }
  abline( h=0, lty=2, col="red" )
  mtext( side=1, line=2, cex=0.8, "Year" )
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.yearResidsPOP


#plt.cohortResids-----------------------2011-08-31
# Plot age residuals by cohort.
#----------------------------------------------AME
plt.cohortResids <- function( obj, ages=c(2,59), pct=c(5,25,50,75,95),                     main=NULL )
{
  # Input is the CAc object from a Awatea res file. Ages to 59 as
  #  plus-age class will mess up year-of-birth calculation. Not automated.

  # par( oma=c(2,1,1,1), mar=c(2,2,2,1), mfrow=c(2,1) )

  # Subset to required ages - still do as don't want age 1 or 60
  #  for cohorts
  obj <- obj[ (obj$Age >= ages[1]) & (obj$Age <=ages[2]), ]
  # obj$stdRes has residuals for each age, year and both sexes. Need
  #  to assign a year of birth for each age as an extra column, then
  #  presumably just do the boxplot split using that.
  obj$birthyr=obj$Year - obj$Age
  if( max(diff(sort(unique(obj$birthyr)))) > 1)
    {
      allYears=min(obj$birthyr):max(obj$birthyr)
      nodataYears=allYears[ !(allYears %in% obj$birthyr)]
      xx=split(c(obj$stdRes, rep(NA, length(nodataYears))), c(obj$birthyr, nodataYears))
      xpos <- boxplot( xx, whisklty=1, xlab="", ylab="", outline=FALSE )     #AME outline=FALSE removes outliers
    } else
    {            
    xpos=boxplot( split( obj$stdRes, obj$birthyr ), whisklty=1, xlab="", ylab="", outline=FALSE )     #AME outline=FALSE removes outliers
    }
  abline( h=0, lty=2, col="red" )
  mtext( side=1, line=2, cex=0.8, "Year of birth" )

  # Plot the q-normal plot of the standardised residuals.
  # qqnorm( obj$stdRes,xlab="",ylab="",main="" )
  # abline( a=0,b=1 )
  # abline( h=quantile(obj$stdRes,p=pct/100,na.rm=TRUE),lty=2 )
  # mtext( side=1, line=2, cex=0.8, "Theoretical quantiles" )

  # mtext( side=2, line=0, cex=1.0, outer=TRUE, "Standardised Residuals" )
  # if ( !is.null(main) )
    mtext( side=3, line=-0.5, cex=1.0, outer=TRUE, main )

  # par( mfrow=c(1,1) )
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.cohortResids


#plt.allTraces--------------------------2011-08-31
# Plot all MCMC traces
#----------------------------------------------AME
plt.allTraces <- function( obj, bioYrList=NULL, recYrList=NULL, save=TRUE )
{
  # Input a Awatea MCMC object.

  plt.trace <- function( obj )
  {
    # Input "obj" is a vector of MCMC samples.
    # Produces one panel trace plot.

    nSample <- length( obj )
    plot( c(1:nSample), obj, type="n", axes=FALSE, xlab="", ylab="" )
    points( c(1:nSample),obj, cex=0.25, pch=16, col="darkgray" )

    # Plot MPD point (1st row).
    points( 1,obj[1], cex=2.0, pch=16, col="green" )
    points( 1,obj[1], cex=2.0, pch=1 )

    lines( lowess( c(1:nSample),obj,f=1/4), lty=1, lwd=1 )
    abline( h=mean(obj), lty=2 )
    axis( side=2 )
    box()
  }

  # Find the active parameters.
  #iPars <- apply( obj$P,2,allEqual )
  #ACH: changed allEqual because you want to find the params that are not all the same (the ones that are estimated)
  iPars <- !apply( obj$P,2,allEqual )
  nPars <- sum( iPars )

  # (1) Plot biomass traces, every 5th year by default (getYrIdx).
  if ( is.null(bioYrList) )
    bioYrList <- getYrIdx( names(obj$B) )
  nYrs <- length( bioYrList )

  graphics( view="landscape" )
  par( oma=c(2,2,1,1), mar=c(1,2,1,1), mfrow=c(3,max(round(nYrs/3),4)) )
  for ( i in 1:nYrs )
  {
    biomass <- obj$B[ ,bioYrList[i] ]
    plt.trace( biomass )
    panLab( 0.5,0.95, cex=1.0, bioYrList[i] )

    mtext( side=1, line=0.5, cex=1.0, outer=TRUE, "Sample" )
    mtext( side=2, line=0.5, cex=1.0, outer=TRUE, "Biomass" )
  }
  if ( save )
    savePlot( paste("biomassTracePJS"), type="jpg" )


  # (2) Plot recruitment traces, every 5th year by default (getYrIdx).

  if ( is.null(recYrList) )
    recYrList <- getYrIdx( names(obj$R) )
  nYrs <- length( recYrList )

  graphics( view="landscape" )
  par( oma=c(2,2,1,1), mar=c(1,2,1,1), mfrow=c(3,max(round(nYrs/3),4)) )
  for ( i in 1:nYrs )
  {
    recruits <- obj$R[ ,recYrList[i] ]
    plt.trace( recruits )
    panLab( 0.5,0.95, cex=1.0, recYrList[i] )

    mtext( side=1, line=0.5, cex=1.0, outer=TRUE, "Sample" )
    mtext( side=2, line=0.5, cex=1.0, outer=TRUE, "Recruitment" )
  }
  if ( save )
    savePlot( paste("recruitTracePJS"), type="jpg" )

  # Plot the active parameters.
  graphics( view="landscape" )
  par( oma=c(2,2,1,1), mar=c(1,2,1,1), mfrow=c(3,min(round(nPars/3),3)) )

  # Find the active parameters.
  iPars <- apply( obj$P,2,allEqual )
  activePars <- obj$P[,!iPars]
  parNames <- names( activePars )
  nPars <- ncol( activePars )

  for ( i in 1:nPars )
  {
    parVec <- activePars[ ,i]
    plt.trace( parVec )
    panLab( 0.5,0.95, cex=1.0, parNames[i] )

    mtext( side=1, line=0.5, cex=1.0, outer=TRUE, "Sample" )
  }
  if ( save )
    savePlot( paste("parTracePJS"), type="jpg" )

  par( mfrow=c(1,1) )
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.allTraces


#plt.expRate----------------------------2011-08-31
# Input an object from "load.allResFiles".
# Plot exploitation rate against year.
#----------------------------------------------AME
plt.expRate <- function( obj, yLim=c(0,0.5), xLim=c(1954,2005) )
{
  nPanels <- length(obj)

  par( oma=c(2,2,1,1), mar=c(2,2,1,1), mfrow=c(2,round(nPanels/2)) )

  for ( i in 1:nPanels )
  {
    res <- obj[[i]]
    year <- res$B$Year

    plot( year, res$B$U, type="n",
      xlab="", xlim=xLim, ylab="", ylim=yLim )
    lines( year, res$B$U, lwd=2, lty=1 )

    yrTicks <- as.numeric( year[ year %% 10==0 ] )

    axis( side=1, at=yrTicks )
    axis( side=2 )
    axis( side=4, labels=FALSE )
    box()

    panLab( 0.5,0.95, names(obj)[i] )

    mfg <- par( "mfg" )
    if ( mfg[1]==mfg[3] & mfg[2]==mfg[4] | i==nPanels )
    {
      mtext( side=1, line=0.5, cex=1.0, outer=TRUE, "Year" )
      mtext( side=2, line=0.5, cex=1.0, outer=TRUE, "Exploitation Rate" )
    }
  }

  par( mfrow=c(1,1) )
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.expRate


#plt.idx--------------------------------2017-12-01
# AME doing postscript for POP. Adapting to five surveys for YMR
#-------------------------------------------AME/RH
plt.idx <- function(obj, main="Residuals", save=NULL, ssnames=paste("Ser",1:9,sep=""),
   ptypes=tcall(PBSawatea)$ptype, pngres=400, ...)
{
	sType = substring(rev(as.character(substitute(obj)))[1],1,1)
	seriesList = sort( unique( obj$Series ) )
	nseries    = length(seriesList)
	surveyFigName = paste0(ifelse(sType=="S","surv",ifelse(sType=="C","cpue","unkn")),"Res",ssnames)
	surveyFigName = gsub(" ","",surveyFigName)
#browser();return()
	surveyHeadName=if (!exists("tcall")) ssnames else tcall(PBSawatea)[[paste0(sType,"names")]]

	for ( i in 1:nseries )
	{
		idx <- seriesList[i]==obj$Series
		result <- stdRes.index( obj[idx,], label=paste(main,"Series",i) )
		pname = surveyFigName[i]
		for (p in ptypes) {
			if (p=="eps") postscript(paste0(pname,".eps"), height=6, width=5, horizontal=FALSE,  paper="special")
			else if (p=="png") png(paste0(pname,".png"), units="in", res=pngres, height=6, width=5)
			par(mfrow=c(nseries,1), mar=c(1.75,2,2,0.5), oma=c(2,2,0,0), mgp=c(2,0.75,0))
			plt.stdResids( result, xLim=range(result[!is.na(result$Obs), ]$Year)) # restrict years for plot
			mtext( side=3, line=0, cex=1.0, outer=TRUE, surveyHeadName[i])
			if (p %in% c("eps","png")) dev.off()
		}
	}
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.idx


#plotIndexNotLattice--------------------2018-04-06
# Taking some of plt.idx, but doing plot.Index NOT as lattice
# obj = currentRes
# plotCI is now custom function in PBSawatea (not gplots)
#-------------------------------------------AME/RH
plotIndexNotLattice <- function(obj, main="", save=NULL,
   bar=1.96, ssnames=paste("Ser",1:9,sep=""),
   ptypes=tcall(PBSawatea)$ptype, pngres=400, ...)
{
	cvcol="slategrey"
	objSurv=obj$Survey; objCPUE=obj$CPUE
	seriesList <- sort( unique( objSurv$Series ) )   # sort is risky if not always in same order
	nseries=length(seriesList)
	surveyHeadName = if (!exists(".PBSmodEnv")) PBSawatea$Snames else tcall(PBSawatea)$Snames
	surveyHeadName = gsub("QC Sound","QCS",surveyHeadName)
	cvpro = tcall(PBSawatea)$cvpro
	if (is.null(cvpro) || all(cvpro==FALSE)) cvpro="unknown"

	# (1) Plot the survey indices 
	yrTicks=min(objSurv$Year):max(objSurv$Year)
	rc = .findSquare(nseries)
	for (p in ptypes) {
		if (p=="eps") postscript("survIndSer.eps", width=6.5, height=8.5, horizontal=FALSE,  paper="special")
		else if (p=="png") png("survIndSer.png", units="in", res=pngres, width=6.5, height=8.5)
		par(mfrow=c(rc[1],rc[2]), mar=c(2,2,1.25,0.5), oma=c(1.5,1.5,1,0.5), mgp=c(1.75,0.5,0))
		for ( i in 1:nseries ) {
			idx <- seriesList[i]==objSurv$Series
			seriesVals=objSurv[idx,]
			# seriesvals$Obs=seriesvals$Obs   # /q[i] - set to 1 anyway
			seriesVals$Hi <- seriesVals$Obs * exp(bar * seriesVals$CV)
			seriesVals$Lo <- seriesVals$Obs/exp(bar * seriesVals$CV)
			yearsnotNA=seriesVals[ !is.na(seriesVals$Obs), ]$Year
			# yearsPlot=min(yearsnotNA):max(yearsnotNA)
			xLim=range(yearsnotNA)
			if(i==1)
				xLimAll=xLim    # range to use in next plot
			xLimAll=range(xLim, xLimAll)    # range to use in next plot
			yLim=c(0, max(seriesVals$Hi, na.rm=TRUE))
			plotCI(seriesVals$Year, seriesVals$Obs, ui=seriesVals$Hi, li=seriesVals$Lo,
				xlim=xLim, ylim=yLim, xlab="", ylab="", gap=0, pch=19) # restrict years for plot, does error bars
			lines(seriesVals$Year, seriesVals$Fit, lwd=2)
			axis( side=1, at=yrTicks, tcl=-0.2, labels=FALSE )
			mtext( side=3, line=0.25, cex=ifelse(nseries>2,1.2,1.5), outer=FALSE, surveyHeadName[i]) #  outer=TRUE
			if (is.numeric(cvpro[i]) && round(cvpro[i],5)!=0)
				addLabel(0.95,0.95,paste("+ CV process error ",cvpro[i],sep=""),adj=c(1,1),cex=0.8,col=cvcol)
			if(i==nseries) {
				mtext(side=2, line=0, cex=1.5, outer=TRUE,"Relative biomass")
				mtext(side=1, line=0, cex=1.5, outer=TRUE, "Year")
			}
		}
		if (p %in% c("eps","png")) dev.off()
	}  # cex was 0.8 for POP
#browser();return()

	# (2) And again, but with the same year axis for each. Think will be instructive to see
	ymaxsurvIndSer3=0           # For ymaxsurvIndSer3.eps
	xLimmaxsurvIndSer3=NA
	XLIM=numeric()
	for (p in ptypes) {
		if (p=="eps") postscript("survIndSer2.eps", width=6.5, height=8.5, horizontal=FALSE,  paper="special")
		else if (p=="png") png("survIndSer2.png", units="in", res=pngres, width=6.5, height=8.5)
		par(mfrow=c(nseries,1), mar=c(0,2,0,0.5), oma=c(3.2,3,0.5,0), mgp=c(1.75,0.5,0))
		for ( i in 1:length(seriesList) ) {
			idx <- seriesList[i]==objSurv$Series
			yrTicks=as.numeric( objSurv$Year[idx])        ## all years
			#YrTicks=intersect(seq(1900,2100,10),yrTicks) ## decadal ticks (doesn't work for short survey in odd years)
			XLIM = range(c(XLIM,seriesVals$Year[!is.na(seriesVals$Obs)]))
			YrTicks=pretty(XLIM)
			seriesVals=objSurv[idx,]
			# seriesvals$Obs=seriesvals$Obs   # /q[i] - set to 1 anyway
			seriesVals$Hi <- seriesVals$Obs * exp(bar * seriesVals$CV)
			seriesVals$Lo <- seriesVals$Obs/exp(bar * seriesVals$CV)
			yearsnotNA=seriesVals[ !is.na(seriesVals$Obs), ]$Year
			# yearsPlot=min(yearsnotNA):max(yearsnotNA)
			# xLim=range(yearsnotNA)
			yLim=c(0, max(seriesVals$Hi, na.rm=TRUE))
			# For axis for survIndSer3.eps:
			ymaxsurvIndSer3=max(ymaxsurvIndSer3, max(seriesVals$Obs, na.rm=TRUE)/ mean(seriesVals$Obs, na.rm=TRUE) )
			xLimmaxsurvIndSer3=range(xLimmaxsurvIndSer3, yearsnotNA, na.rm=TRUE)     # setting xLimmaxsurvIndSer3=NA above
			#gplots::plotCI(seriesVals$Year, seriesVals$Obs, ui=seriesVals$Hi,
			plotCI(seriesVals$Year, seriesVals$Obs, ui=seriesVals$Hi, li=seriesVals$Lo, xaxt="n", xlim=xLimAll, ylim=yLim, xlab="", ylab="", gap=0, pch=19)
			# restrict years for plot, does error bars
			lines(seriesVals$Year, seriesVals$Fit, lwd=2, col="blue")
			axis( side=1, at=yrTicks, tcl=-0.2, labels=FALSE)
			axis( side=1, at=YrTicks, tcl=-0.4, labels=ifelse(i==nseries,TRUE,FALSE),cex.axis=1.2)
			#mtext(side=4, line=1.5, cex=0.8, outer=FALSE, paste0(strwrap(surveyHeadName[i],10),collapse="\n"))
			mtext(side=2, line=1.75, cex=ifelse(nseries>2,1,1.5), outer=FALSE, surveyHeadName[i], col="blue")
			if (is.numeric(cvpro[i]) && round(cvpro[i],5)!=0)
				addLabel(0.95,0.95,paste("+ CV process error ",cvpro[i],sep=""),adj=c(1,1),cex=.8+(.05*(nseries-1)),col=cvcol)
			if(i==nseries) {
				#axis( side=1, at=yrTicks, tcl=-0.2, labels=FALSE)
				#axis( side=1, at=YrTicks, tcl=-0.4, labels=TRUE)
				mtext(side=1, line=2, cex=1.5, outer=TRUE, "Year")
				mtext(side=2, line=1.2, cex=1.5, outer=TRUE, "Relative biomass")
			}
		}  # cex was 0.8 for POP
		if (p %in% c("eps","png")) dev.off()
	}
#browser();return()

	# (3) And again, but all series on same plot, normalised to their means to see trends. Maybe add CPUE also.
	# Calculate max of normalised surveys and CPUEs
	objsurv=objSurv[!is.na(objSurv$Obs),]
	norsurv=sapply(split(objsurv$Obs,objsurv$Series),function(x){xx=x[!is.na(x)]; if (length(xx)==0) 0 else xx/mean(xx)},simplify=FALSE)
	maxsurv=max(sapply(norsurv,max))
	yrssurv=range(objsurv$Year); yrsspan=yrssurv[1]:yrssurv[2]
	objcpue=objCPUE[is.element(objCPUE$Year,yrsspan),]
	norcpue=sapply(split(objcpue$Obs,objcpue$Series),function(x){xx=x[!is.na(x)]; if (length(xx)==0) 0 else xx/mean(xx)},simplify=FALSE)
	maxcpue=max(sapply(norcpue,max))

	for (p in ptypes) {
		if (p=="eps") postscript("survIndSer3.eps", width=6.5, height=6.5, horizontal=FALSE,  paper="special")
		else if (p=="png") png("survIndSer3.png", units="in", res=pngres, width=6.5, height=6.5)
		par(mfrow=c(1,1), mar=c(3,3,0.5,0.5), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
		#postscript("survIndSer3.eps", height=6.0, width=6.0, horizontal=FALSE,  paper="special")   # height was 6 for POP
		yrTicks=yrsspan
		#YrTicks=intersect(seq(1900,2100,10),yrTicks) ## decadal ticks (doesn't work for short survey in odd years)
		YrTicks=pretty(xLim)
		NSL=1:length(seriesList)
		CLRS=c("black","blue","red","green4","orange","purple","navy","salmon")
		clrs=rep(CLRS,max(NSL))[NSL]
		# Set up plot
		plot(NA, xlim=yrssurv, ylim=c(0, max(maxsurv,maxcpue)), xlab="Years", ylab="Survey indices normalised by means", xaxt="n")
		abline(h=1, col="grey")
		axis(1,at=yrTicks,tck=-0.01,labels=FALSE)
		axis(1,at=YrTicks,tck=-0.02,labels=TRUE)
		for ( i in NSL ) {
			idx <- is.element(objsurv$Series,seriesList[i])
			seriesVals=objsurv[idx,]
			yearsnotNA=seriesVals[ !is.na(seriesVals$Obs), ]$Year
			points(yearsnotNA, seriesVals$Obs[ !is.na(seriesVals$Obs)] / mean(seriesVals$Obs, na.rm=TRUE), pch=i, col=clrs[i], type="o")
		}
		legtxt = if (exists(".PBSmodEnv")) tcall(PBSawatea)$Snames else PBSawatea$Snames
		legtxt = gsub("QC Sound","QCS",legtxt)
		# Now draw on CPUE series also:
		nseries = length(seriesList)  #  Need nseries from surveys for col
		seriesListCPUE <- sort( unique( objcpue$Series ) )
		ncpue=length(seriesListCPUE)
		cpue = as.logical(obj$extra$likelihoods$CPUE)
		# sort risky if not always in same order.
		lty=rep(1,length(NSL))
		if(cpue) {
			lty=c(lty,rep(2,ncpue))
			NSL=c(NSL, (max(NSL)+1):(max(NSL)+ncpue))
			clrs=c(clrs,rep(CLRS,ncpue)[1:ncpue])
			legtxt=c(legtxt, gsub("Series","CPUE",seriesListCPUE))
			for ( i in 1:ncpue ) {
				idx <- is.element(objcpue$Series,seriesListCPUE[i])
				seriesVals=objcpue[idx,]
				yearsnotNA=seriesVals[ !is.na(seriesVals$Obs), ]$Year
				if (length(yearsnotNA)>0)
					points(yearsnotNA, seriesVals$Obs[ !is.na(seriesVals$Obs)] / mean(seriesVals$Obs, 
						na.rm=TRUE), pch=i+nseries, col=clrs[i+nseries], type="o", lty=2)
			}
		}
		addLegend(0.075, 0.975, bty="n", col=clrs, pch=NSL, legend=legtxt, cex=0.8)
		if (p %in% c("eps","png")) dev.off()
	}

	# (4) And again, but big figures for ease of viewing. 
	for ( i in 1:length(seriesList) ) {
		idx <- seriesList[i]==objSurv$Series
		seriesVals=objSurv[idx,]
		# seriesvals$Obs=seriesvals$Obs   # /q[i] - set to 1 anyway
		seriesVals$Hi <- seriesVals$Obs * exp(bar * seriesVals$CV)
		seriesVals$Lo <- seriesVals$Obs/exp(bar * seriesVals$CV)
		yearsnotNA=seriesVals[ !is.na(seriesVals$Obs), ]$Year
		# yearsPlot=min(yearsnotNA):max(yearsnotNA)
		xLim=range(yearsnotNA)
		yLim=c(0, max(seriesVals$Hi, na.rm=TRUE))
		for (p in ptypes) {
			if (p=="eps") postscript(paste0("survIndSer4-", i, ".eps"), width=6.5, height=6.5, horizontal=FALSE,  paper="special")
			else if (p=="png") png(paste0("survIndSer4-", i, ".png"), units="in", res=pngres, width=6.5, height=6.5)
			par(mfrow=c(1,1), mar=c(3,3.25,2,1), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
			plotCI(seriesVals$Year, seriesVals$Obs, ui=seriesVals$Hi, li=seriesVals$Lo, 
				xlim=xLim, ylim=yLim, xlab="Year", ylab="Relative biomass", gap=0, pch=19, cex.lab=1.5)
			lines(seriesVals$Year, seriesVals$Fit, lwd=2)
			axis( side=1, at=yrTicks, tcl=-0.2, labels=FALSE )
			mtext( side=3, line=0.25, cex=1.5, outer=FALSE, surveyHeadName[i]) #  outer=TRUE
			if (is.numeric(cvpro[i]) && round(cvpro[i],5)!=0)
				addLabel(0.95,0.95,paste("+ CV process error ",cvpro[i],sep=""),adj=c(1,0),cex=0.8,col=cvcol)
			if (p %in% c("eps","png")) dev.off()
		}
	}  # cex was 0.8 for POP
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotIndexNotLattice



#plt.numR-------------------------------2011-08-31
# Input an object from "load.allResFiles".
# Plot numbers at age at equilibrium.
# Plot recruitment age 1's.
#----------------------------------------------AME
plt.numR <- function( obj, minYr=NULL )
{
  nPanels <- length(obj)
  par( oma=c(2,2,1,1), mar=c(2,2,1,1), mfcol=c(2,2) )
  for ( i in 1:nPanels )
  {
    res <- obj[[i]]
    # Plot numbers at age for initial conditions.
    x <- res$N
    x <- x[ x$Year==min(x$Year), ]
    plot( x$Age, x$N, type="n", xlab="Age", ylab="", ylim=c(0,max(x$N)) )
    delta <- 0.4
    rect( x$Age-delta, 0, x$Age+delta, x$N, density=-1, col="gray" )
    axis( side=4, labels=FALSE )
    box()
    panLab( 0.5,0.95, names(obj)[i] )
    mfg <- par( "mfg" )
    if ( mfg[1]==1 & mfg[2]==1 )
      mtext( side=2, line=2.5, cex=1.0, "Initial numbers (000s)" )
    if ( is.null( minYr ) )
      minYr <- min(x$Year)
    # Age-1 recruits against year.
    x <- res$N
    x <- x[ x$Year<max(x$Year),]
    x <- x[ x$Age==min(x$Age), ]
    plot( x$Year, x$N, type="n", xlab="", ylab="", ylim=c(0,max(x$N)) )
    abline( v=seq(1955,2005,5 ), lty=3 )
    abline( v=seq(1950,2010,10), lty=2 )
    x <- x[ x$Year >= minYr, ]
    delta <- 0.35
    rect( x$Year-delta, 0, x$Year+delta, x$N, density=-1, col="gray" )
    mfg <- par( "mfg" )
    if ( mfg[1]==2 & mfg[2]==1 )
      mtext( side=2, line=2.5, cex=1.0, "Age 1's (000s)" )
    if ( mfg[1]==mfg[3] & mfg[2]==mfg[4] )
    {
      mtext( side=1, line=0, cex=1.0, outer=TRUE, "Year" )
    }
  }
  par( mfrow=c(1,1) )
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.numR


#plt.quantBioBB0------------------------2011-08-31
# This is for B/B0 for each run, with just one projection. Doing one
#  plot here instead of multiple, so taking some from plotBVBnorm
#  which did one for each run. Don't think used for biomasses.
#  Now using for single recruitment projection plot
#----------------------------------------------AME
plt.quantBioBB0 <- function( obj, projObj=NULL, policy=NULL,
                  p = get("quants5"),
                  xyType="lines",
                  lineType=c(3,2,1,2,3),
                  refLines=NULL, xLim=NULL, yLim=NULL,
                  userPrompt=FALSE, save=TRUE, main="", cex.main="",
                  tcl.val=-0.2, xaxis.by=1, yaxis.by=10000,
                  xaxis.lab="Year", yaxis.lab= "Spawning biomass" )
{
  plt.qB <- function( obj, xyType="lines", new=TRUE, xLim, yLim, line.col="black", med.col="black", ... )   #AME line.col="black", col is for filling rect med.col is for median
  {
    if ( new )
      plot( xLim,yLim, type="n", xlab=xaxis.lab,ylab=yaxis.lab )

    yrs <- as.numeric(dimnames(obj)[[2]])

    # Connect the quantiles with lines.
    if ( xyType=="lines" )
    {
      for ( i in 1:nrow(obj) )
      {
        # Plot reconstructed biomass.
        lines( yrs,obj[i,], lty=lineType[i],... )
      }
    }

    # ARK vertical line-dot plot.
    # Assumes that five quantiles requested, with median as one.
    if ( xyType=="lineDot" )
    {
      points( yrs,obj[2,], pch=3, cex=0.5,... )
      points( yrs,obj[3,], pch=1,... )
      points( yrs,obj[4,], pch=3, cex=0.5,... )
      segments( yrs,obj[1,], yrs,obj[5,], lty=1,... )
    }

    # Quantile boxplots - assumes five quantiles.
    if ( xyType=="quantBox" )
    {
      delta <- 0.25
      # Draw the outer whiskers.
      segments( yrs,obj[1,], yrs,obj[5,], lty=1, col=line.col) #, ... )
             #AME col=1 removed, col=line.col ... so can have red for projs
      # Overlay the box.
      for ( i in 1:length(yrs) )
        rect( yrs[i]-delta,obj[2,i], yrs[i]+delta, obj[4,i],
             border=line.col, col="white")#AME border,col=NA (empty)
      # Add the median.
      # segments( yrs-delta,obj[3,],yrs+delta,obj[3,],lty=1,col=line.col )   #AME black
      segments( yrs-delta,obj[3,],yrs+delta,obj[3,],lty=1,col=med.col )   #AME black
    }
  }

  # Plot quantiles of biomass using the posterior densities.
  # If proj!=NULL then add the projections for all policies.
  # If policy!=NULL than plot only the specified policy.

  # Just one plot now. Deleted some.

  yrs1 <- NULL
  yrs2 <- NULL
  result1 <- NULL
  result2 <- NULL

  # Calculate the quantiles of the reconstructed biomass.
  result1 <- apply( obj,2,quantile,probs=p )
  yrs1 <- as.numeric(dimnames(result1)[[2]])

  if ( is.null(yLim) )
    yLim <- range( result1 )

  # Reconstructed biomass.
  if ( is.null(projObj) )
  {
    plt.qB( result1,xLim=range(yrs1),yLim=yLim, xyType=xyType )
    # mtext( side=1, line=2, cex=1.0, "Year" )
    # mtext( side=2, line=2, cex=1.0, "Biomass" )
  }

  # Reconstructed biomass and projections.
  if ( !is.null(projObj) )
  {
    # Get the policies to be plotted.
    if ( is.null(policy) )
      policyList <- names( projObj )
    else
      policyList <- policy

    # Loop over the policies.
    result2 <- as.list(policyList)
    names( result2 ) <- policyList

    iPage <- 1
    nPolicies <- length(policyList)
    for ( j in 1:nPolicies )
    {
      # Calculate quantiles of the projected biomass for policy.
      pol <- policyList[j]
      result2[[j]] <- apply( projObj[[pol]],2,quantile,probs=p )

      # cat( "\n\nQuantiles of projection for policy=",
      #   policyList[j],"\n" )
      # print( result2[[j]] )

      yrs2 <- as.numeric( dimnames(result2[[j]])[[2]] )
      if ( is.null(xLim) )
        xLim <- range( yrs1,yrs2 )
      # yLim <- range( result1,result2[[j]] )  # AME to get same axes
      # yLim <- c(0,yLim[2])

      # Plot the quantiles of the biomass.
      if ( xyType=="quantBox" )
        plt.qB( result1, xyType=xyType, new=TRUE, xLim,yLim, col="black", line.col="black", med.col="red" )      
      else
        plt.qB( result1, xyType=xyType, new=TRUE, xLim,yLim, col="red") # line.col="red")

      if ( !is.null(refLines) )
        abline( v=refLines,lty=4 )

      # Plot the quantiles of the projected biomass.
      if ( xyType=="quantBox" )
        plt.qB( result2[[j]], xyType=xyType, new=FALSE, xLim,yLim, line.col="red", med.col="black")  # AME: col fills in box, I want to change line cols 
      else
        plt.qB( result2[[j]], xyType=xyType, new=FALSE, xLim,yLim, col="red" )

      #for ( i in 1:nrow(result2[[j]]) )
      #  lines( yrs2,result2[[j]][i,],lty=lineType[i],lwd=2,col=2 )

      abline( v=yrs2[1]-0.5, lty=2 )
      axis(1, at=seq(xLim[1], xLim[2], by=xaxis.by), tcl=tcl.val, labels=FALSE)
      axis(2, at=seq(0, yLim[2], by=yaxis.by), tcl=tcl.val, labels=FALSE)
      
      # mfg <- par( "mfg" )
      # if ( mfg[1]==1 & mfg[2]==1 )
      # {
      #  mtext( side=1, line=0.8, cex=1.0, outer=TRUE, "Year" )
      #                                 #AME line=0 changed
      #  mtext( side=2, line=0.8, cex=1.0, outer=TRUE,
      #        "Spawning biomass" ) # "   and Spawning
      #}
      #mtext( side=3, line=0.25, cex=0.8,
      #  paste( "Policy:",policyList[j]) )

      # Last panel on page or last policy.
      #if ( mfg[1]==mfg[3] & mfg[2]==mfg[4] | j==nPolicies )
      #{
      #  if ( save )
      #    savePlot( paste( "policyProj",iPage,sep=""),type="png" )
      #  iPage <- iPage + 1
      #
      #  if ( j < nPolicies )
      #  {
      #    windows()
      #    par( oma=c(2,2,1,1), mar=c(2,2,1.5,1),
      #         mfcol=c(nRow,nCol), userPrompt )
      #  }
      # }
    }
  }
  # par( mfrow=c(1,1) )

  val <- list( recon=result1, proj=result2 )
  val
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.quantBioBB0


#plt.ssbVbCatch-------------------------2011-08-31
# Plot the spawning stock and vulnerable biomass.
#----------------------------------------------AME
plt.ssbVbCatch <- function( obj, x1=1966, xLim=c(1954,2005), yLim=c(0,25000) )
{
  # Input an object from "load.allResFiles".
  # Plot spawning biomass, vulnerable biomass, catch against year.

  nPanels <- length(obj)

  par( oma=c(2,2,1,1), mar=c(2,2,1,1), mfrow=c(2,round(nPanels/2)) )

  for ( i in 1:nPanels )
  {
    res <- obj[[i]]
    year <- res$B$Year
    plot( year, res$B$SB, type="n", axes=FALSE,
      xlab="", xlim=xLim, ylab="", ylim=yLim )

    idx <- year >= x1
    lines( year[idx], res$B$SB[idx], lwd=2, lty=1 )
    lines( year[idx], res$B$VB[idx], lwd=2, lty=5 )

    yrTicks <- as.numeric( year[ year %% 10==0 ] )

    axis( side=1, at=yrTicks )
    axis( side=2 )
    axis( side=4, labels=FALSE )
    box()

    delta <- 0.35
    rect( year-delta, 0, year+delta, res$B$Y, density=-1, col="gray" )

    panLab( 0.5,0.95, names(obj)[i] )
    panLegend( 0.6, 0.95, legTxt=c("SB","VB"),
               ncol=1, lty=c(1,5), lwd=c(2,2), bty="n" )

    mfg <- par( "mfg" )
    if ( mfg[1]==mfg[3] & mfg[2]==mfg[4] | i==nPanels )
    {
      mtext( side=1, line=0.5, cex=1.0, outer=TRUE, "Year" )
      mtext( side=2, line=0.5, cex=1.0, outer=TRUE, "Metric tonnes" )
    }
  }
  par( mfrow=c(1,1) )
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.ssbVbCatch


#plt.stdResids--------------------------2011-12-01
# Plot standardised residuals (AME adding xlim's for POP).
#----------------------------------------------AME
plt.stdResids <- function( obj, pct=c(5,25,50,75,95),
                           main=NULL, yLim=NULL, xLim=xLim )
{
  # Input is a data.frame with columns "Year", "stdRes", "Fit".
  # cex.val=1.2       # size of text, if want to shrink down
  par( oma=c(1,1,2,1), mar=c(3,3,1,1), mfrow=c(3,1) )

  if ( is.null(yLim) )
    yLim <- range( obj$stdRes, na.rm=TRUE )
	ptcex=1; ptpch=19

  # Plot the standardised residuals against time.
  plot( obj$Year, obj$stdRes, type="n", xlab="", ylab="", ylim=yLim,
       xlim=xLim)
  abline( h=0, lty=2 )
  points( obj$Year, obj$stdRes, cex=ptcex, pch=ptpch) #, bg="orange" )
  mtext( side=1, line=2, cex=0.8, "Year" )

  # Plot the standardised residuals against predicted values.
  zin = !is.na(obj$Fit) & !is.na(obj$stdRes)
  plot( log(obj$Fit[zin]), obj$stdRes[zin], type="n", xlab="", ylab="", ylim=yLim )
  abline( h=0, lty=2 )
  points( log(obj$Fit[zin]), obj$stdRes[zin], cex=ptcex, pch=ptpch) #, bg="orange" )
  mtext( side=1, line=2, cex=0.8, "Predicted" )

  # Plot the q-normal plot of the standardised residuals.
  wiggle=qqnorm( obj$stdRes,xlab="",ylab="",main="" , pch=ptpch, cex=ptcex)
  abline( a=0,b=1 )
  abline( h=quantile(obj$stdRes,p=pct/100,na.rm=TRUE), lty=c(3,2,1,2,3), lwd=0.75 )
  points(wiggle, pch=ptpch, cex=ptcex) #, bg="orange")
  mtext( side=1, line=2, cex=0.8, "Theoretical quantiles" )

  mtext( side=2, line=-0.5, cex=0.8, outer=TRUE, "Standardised residuals" )
  if ( !is.null(main) )
    mtext( side=3, line=0, cex=1.0, outer=TRUE, main )
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.stdResids


#importMCMC.ddiff-----------------------2011-08-31
#  Import Functions for PJS Delay Difference Model
# Purpose is to make an plotMCMC object identical in format to
# the result of importMCMC from PJS delay difference model output.
# The difference is that B is biomass defined by ddiff model.
#----------------------------------------------AME
importMCMC.ddiff <- function()
{
  # Get the likelihood MCMCs. PJS includes other derived parameters here.
  L <- read.table( "mcmclike.csv", header=TRUE, sep="," )

  # Get the parameter MCMCs.
  params <-read.table( "mcmcparams.csv", header=TRUE, sep="," )
  rdevID <- grep( "rdev",names(params) )
  P <- params[ ,setdiff(c(1:ncol(params)),rdevID ) ]

  # Get the biomass MCMCs and strip "biom" off to leave only "yyyy".
  B <- read.table( "mcmcbiom.csv", header=TRUE, sep="," )
  names( B ) <- substring( names( B ),5 )

  # Get the recruitments.  They are already in the mcmcparams file
  # with headers "rdev_yyyy".
  R <- params[ ,rdevID ]

  # Now rename with year only as per plotMCMC projection object.
  names( R ) <- gsub( "rdev_","",names(R) )

  result <- list( L=L, P=P, B=B, R=R )
  result
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~importMCMC.ddiff


#importProj.ddiff-----------------------2011-08-31
# Purpose is to make an plotMCMC object identical in format to
# the result of importProj from PJS delay difference model output.
# The difference is that B is biomass defined by ddiff model.
#----------------------------------------------AME
importProj.ddiff <- function( yrVal="2006" )
{
  # Get the biomass projection - PJS does 1 year ahead projection.
  # The column "X" appears as the last column because trailing "," exist in the mcmcprojbiom.csv file.
  # Note also that "cat=" in csv file becomes "cat." in read.table.
  projbiom <- read.table( "mcmcprojbiom.csv",header=TRUE,sep="," )
  projNames <- gsub( "cat.","",names( projbiom ) )
  names( projbiom ) <- projNames

  # Remove the last column of NA values created by trailing ",".
  projbiom <- projbiom[ ,names(projbiom)!="X" ]
  projNames <- names( projbiom )

  # Convert the columns of projbiom to matrices as list elements.
  B <- as.list( 1:length(projNames) )
  names( B ) <- projNames
  for ( j in 1:length(projNames) )
  {
    B[[j]] <- matrix( projbiom[,j],ncol=1 )
    dimnames( B[[j]] ) <- list( NULL,yrVal )
  }
  # Get the harvest rate projection - PJS does 1 year ahead projection.
  # These appear to be F's rather than annual harvest rates?
  # The column "X" appears as the last column because trailing ","
  # exist in the mcmcprojbiom.csv file.
  # Note also that "cat=" in csv file becomes "cat." in read.table.
  projhr <- read.table( "mcmcprojhr.csv",header=TRUE,sep="," )
  projNames <- gsub( "cat.","",names( projhr ) )
  names( projhr ) <- projNames

  # Remove the last column of NA values created by trailing ",".
  projhr <- projhr[ ,names(projhr)!="X" ]
  projNames <- names( projhr )

  # Convert instantaneous F's to annual rates.
  projU <- (1.0-exp(-projhr))

  # Calculate yield (catch) at annual rate.
  catch <- projbiom * projU

  # Convert the columns of projbiom to matrices as list elements.
  Y <- as.list( 1:length(projNames) )
  names( Y ) <- projNames
  for ( j in 1:length(projNames) )
  {
    Y[[j]] <- matrix( catch[,j],ncol=1 )
    dimnames( Y[[j]] ) <- list( NULL,yrVal )
  }

  result <- list( B=B, Y=Y )
  result
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~importProj.ddiff


#msyCalc--------------------------------2012-10-16
# To load in MSY.out and calculated the MSY.
# Call this function with msy=msyCalc().
#----------------------------------------------AME
msyCalc=
function (dir=getwd(), error.rep=1) 
{
# rewriting msyCalc here to report the convergence numbers:
	control=readLines(paste(dir, "/Yields.ctl", sep=""))
	maxProj=as.numeric(control[3])
	tolerance=as.numeric(control[5])
	test=read.table(paste(dir, "/MSY.out", sep=""), header=TRUE, sep="\t")
	num.draws=dim(test)[1]
	nProjIndex=grep("nProj", names(test))
	nProjMat=test[, nProjIndex]   # matrix of number of projections
	nProjMatTF=rowSums(nProjMat > maxProj - 1)   # sum by row
	if (error.rep == 1 & (sum(nProjMatTF) > 0)) {
		stop(paste("Simulations reach maximum year for", sum(nProjMatTF), 
			"of the", num.draws * dim(nProjMat)[2], "simulations, so need to run for longer or reduce tolerance to reach equilibrium"))
	}
	yieldIndex=grep("Yield", names(test))
	yieldMat=test[, yieldIndex]
	uIndex=grep("U", names(test))
	uMat=test[, uIndex]
	VBIndex=grep("VB_", names(test))
	VBMat=test[, VBIndex]
	SBIndex=grep("SB_", names(test))
	SBMat=test[, SBIndex]
	imsy=apply(yieldMat, 1, which.max)
	if (error.rep == 1 & max(imsy) == dim(yieldMat)[2]) {
		stop("Need to use a higher max U to reach the MSY, for at least one draw")
	}
	msy=vector()
	umsy=vector()
	VBmsy=vector()
	Bmsy=vector()
	nProj=vector()
	for (i in 1:num.draws) {
		ind=imsy[i]
		msy[i]=yieldMat[i, ind]
		umsy[i]=uMat[i, ind]
		VBmsy[i]=VBMat[i, ind]
		Bmsy[i]=SBMat[i, ind]
		nProj[i]=nProjMat[i, ind]
	}
	return(list(yield=msy, u=umsy, VB=VBmsy, B=Bmsy, 
		nProj=nProj, uMin=uMat[,1], uMax=uMat[,dim(uMat)[2]], imsy=imsy, maxUind=rep(dim(yieldMat)[2], num.draws),  maxProj=rep(maxProj, num.draws), tolerance=rep(tolerance, num.draws), nProjMatTF=nProjMatTF))
	# imsy is index of tested u's for which you get umsy, so
	#  so if it's 1 or maxUind for an MCMC then that one 
	#  reached the bounds of u. Report that below in Sweave.
	# uMin, uMax are vectors of the min/max tested u's
	#  same for all MCMC samples, but need a vector to use
	#  sapply below (and same for the following variables).
	# maxProj is maximum number of projection years tried
	#  (from Yields.ctl file), so need to report if that's
	#  reached for any of the nProj. Again, need a vector.
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~msyCalc


#refPoints------------------------------2011-08-31
# Call from Sweave as  refPoints() or, in full:
# refPoints(currentMCMC, currentProj, currentMSY, refLevels=c(0.4, 0.8, 1))
#----------------------------------------------AME
refPoints <- function( mcmcObj=currentMCMC, projObj=currentProj,
                     msyObj=currentMSY, refLevels=c(0.4, 0.8, 1))
                     # refLevels are %age of msyObj
{
  refPlist=as.list(c("LRP", "URP", "Bmsy"))  # Can't have 0.4Bmsy
  names(refPlist)=c("LRP", "URP", "Bmsy")   # as numeric at start. '0.4Bmsy'
  for(i in 1:length(refLevels))
    {
    refPlist[[i]]=refLevels[i] * msyObj$B
    }
  return(refPlist)
}

refPointsB0 <- function( mcmcObj=currentMCMC, projObj=currentProj,
                     B0Obj=B0.MCMC, refLevels=B0refLevels,
                     refNames=B0refNames)
  {
  refPlist=as.list(refNames)
  names(refPlist)=c(refNames)
  for(i in 1:length(refLevels))
    {
    refPlist[[i]]=refLevels[i] * B0Obj
    }
  return(refPlist)
}

#refPointsHist--------------------------2013-09-27
# Call from Sweave as  refPointsHist(HRP.YRS=ROL.HRP.YRS)
# Originally implemented for Rock Sole 2013
#-----------------------------------------------RH
refPointsHist <- function( mcmcObj=currentMCMC, HRP.YRS)
   #blimYrs=1966:2005, btarYrs=1977:1985, ulimYrs=NULL, utarYrs=1966:2005
{
	unpackList(HRP.YRS)
	refPlist=list(blimHRP=NULL, btarHRP=NULL, ulimHRP=NULL, utarHRP=NULL)
	# Find the minimum B during the limit years
	if (!is.null(blimYrs))
		refPlist[["blimHRP"]]=apply(mcmcObj$B[,as.character(blimYrs),drop=FALSE],1,min)
	# Find the mean B during the target years
	if (!is.null(btarYrs))
		refPlist[["btarHRP"]]=apply(mcmcObj$B[,as.character(btarYrs),drop=FALSE],1,mean)
	# Find the minimum U during the limit years
	if (!is.null(ulimYrs))
		refPlist[["ulimHRP"]]=apply(mcmcObj$U[,findPat(ulimYrs,names(mcmcObj$U)),drop=FALSE],1,min) # RH: in case Ngear > 1
	# Find the mean U during the target years
	if (!is.null(utarYrs))
		refPlist[["utarHRP"]]=apply(mcmcObj$U[,findPat(utarYrs,names(mcmcObj$U)),drop=FALSE],1,mean) # RH: in case Ngear > 1
	return(refPlist)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~refPointsHist


#plotDensPOPparsPrior-------------------2011-08-31
# Adding the prior automatically.
#----------------------------------------------AME
plotDensPOPparsPrior <-
    function (mcmc, probs=get("quants3")[c(1,3)], points=FALSE, axes=TRUE, 
    same.limits=FALSE, between=list(x=axes, y=axes), 
    div=1, log=FALSE, base=10, main=NULL, xlab=NULL, 
    ylab=NULL, cex.main=1.2, cex.lab=1, cex.strip=0.8, 
    cex.axis=0.7, las=0, tck=0.5, tick.number=5, lty.density=1, 
    lwd.density=3, col.density="black", lty.median=2, lwd.median=1, 
    col.median="darkgrey", lty.outer=3, lwd.outer=1, col.outer="darkgrey", 
    pch="|", cex.points=1, col.points="black", plot=TRUE,
    MPD.height=0.04, mpd=mcmc[1,], ...)    # MPD.height, how far up to put MPD
{
    panel.dens <- function(x, ...) {
      if (any(is.finite(x)) && var(x) > 0) 
          panel.densityplot(x, lty=lty.density, lwd=lwd.density, 
                col.line=col.density, plot.points=points, 
                pch=pch, cex=cex.points, col=col.points,
                ...)
      else panel.densityplot(x, type="n", ...)

        panel.abline(v=quantile(x, probs=probs), lty=lty.outer, 
            lwd=lwd.outer, col=col.outer)
        panel.abline(v=median(x), lty=lty.median,
            lwd=lwd.median, 
            col=col.median)
        panel.xyplot(mpd[panel.number()], current.panel.limits()$ylim[2]*MPD.height,
            pch=19, col="red") # AME, MPD. 0.04 of way up
                                        #  assumes ..ylim[1]=0
        panel.xyplot(mpd[panel.number()], current.panel.limits()$ylim[2]*MPD.height,
            pch=1, col="black") #AME
        # panel.curve(priorDistList[[panel.number()]], min=-1, current.panel.limits()$xlim[1], current.panel.limits()$xlim[2], col="blue")
        # panel.curve(priorDistList[[1]], col="blue")
        #   panel.curve(priorDistList[[panel.number()]](x), from=max(priorBoundsList[[panel.number()]][1],
        #      current.panel.limits()$xlim[1]), to=min(priorBoundsList[[panel.number()]][2], 
        #      current.panel.limits()$xlim[2]), col="blue") # need the bounds, from max of lower bound and panel xlim[1], to min of upper bound and panel xlim[2]
        panel.curve(priorDistList[[panel.number()]]
            (x, priorInput[panel.number(), ] ),
            from=max(priorInput[panel.number(), 2] , 
              current.panel.limits()$xlim[1]),
            to=min(priorInput[panel.number(), 3],
              current.panel.limits()$xlim[2]),
            col="blue")
    }
    relation <- if (same.limits) 
        "same"
    else "free"
    if (is.null(dim(mcmc))) {
        mcmc.name <- rev(as.character(substitute(mcmc)))[1]
        mcmc <- matrix(mcmc, dimnames=list(NULL, mcmc.name))
    }
    mcmc <- if (log) 
        log(mcmc/div, base=base)
    else mcmc/div
    mcmc <- as.data.frame(mcmc)
    n <- nrow(mcmc)
    p <- ncol(mcmc)
    x <- data.frame(Factor=ordered(rep(names(mcmc), each=n), 
        names(mcmc)), Draw=rep(1:n, p), Value=as.vector(as.matrix(mcmc)))
    #mess = c(
    #"require(grid, quietly=TRUE, warn.conflicts=FALSE)",
    #"require(lattice, quietly=TRUE, warn.conflicts=FALSE)"
    #)
    #eval(parse(text=mess))
    if (trellis.par.get()$background$col == "#909090") {
        for (d in dev.list()) dev.off()
        trellis.device(color=FALSE)
    }
    mymain <- list(label=main, cex=cex.main)
    myxlab <- list(label=xlab, cex=cex.lab)
    myylab <- list(label=ylab, cex=cex.lab)
    myrot <- switch(as.character(las), `0`=0, `1`=0, `2`=90, 
        `3`=90)
    myscales <- list(y=list(draw=FALSE, relation="free"), 
        x=list(draw=axes, relation=relation, cex=cex.axis, 
            tck=tck, tick.number=tick.number, rot=myrot))
    mystrip <- list(cex=cex.strip)
    graph <- densityplot(~Value | Factor, panel=panel.dens, 
        data=x, as.table=TRUE, between=between, main=mymain, 
        xlab=myxlab, ylab=myylab, par.strip.text=mystrip, 
        scales=myscales, ...)
    if (!log) {
        if (is.list(graph$y.limits)) 
            graph$y.limits <- lapply(graph$y.limits, function(y) {
                y[1] <- 0
                return(y)
            })
        else graph$y.limits[1] <- 0
    }
    if (plot) {
        print(graph)
        invisible(x)
    }
    else {
        invisible(graph)
    }
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotDensPOPparsPrior


#findTarget-----------------------------2013-02-20
#  To derive decision tables for moving windows and find
#  the times to achieve recovery with given confidence.
#   Vmat  =matrix of projected B-values (MCMC projections x Year)
#   yrP   =user-specified projection years
#   yrG   =number of years for moving target window (e.g, 90y=3 YMR generations). Might not work for all possibilities.
#   ratio =recovery target ratio
#   target=recovery target values (e.g., B0, Bmsy).
#           =B0.MCMC for ratios of B0
#           =Bmsy.MCMC for ratios of Bmsy
#           =Bt.MCMC for moving window
#   conf  =confidence level required
#   plotit=logical to plot the probability of Bt/target
#   retVal=character name of object to return
#    retVal="N", look for the global object "Ttab" (number of years
#     to acheive target)
#    retVal="p.hi" gives global object "Ptab", a list of decision
#     tables where row is the catch option and column is the year
#     Values are probabilities of acheiving target.
#  (2012-02-20) 'xhi=x>=r' change to 'xhi=x>r'
#-----------------------------------------------RH
findTarget=function(Vmat, yrU=as.numeric(dimnames(Vmat)[[2]]), yrG=90, ratio=0.5, target=B0.MCMC,
    conf=0.95, plotit=FALSE, retVal="N") {
	
	# oldpar=par(no.readonly=TRUE);  on.exit(par(oldpar))
	yrA   =as.numeric(dimnames(Vmat)[[2]])   # years available
	yrP   =sort(intersect(yrA,yrU))          # years for proj
	yr0   =yrP[1]; yrN=rev(yrP)[1]         # 

	vmat=Vmat[,is.element(dimnames(Vmat)[[2]],as.character(yrP))]             # include only yrP years
	if (is.data.frame(target) || is.matrix(target)) {
		yrM  =yrP - yrG                                                        # moving target years
		yrM1 =intersect(as.numeric(dimnames(target)[[2]]),yrM)                 # available target years from MCMC
		if (length(yrM1)==0) {                                                   # projection not long enough for any overlap with 3 generations
			if (retVal=="N") return(NA)
			else {p.hi=rep(NA,length(yrP)); names(p.hi)=yrP }; return(p.hi) }
		yrMr =range(yrM1)                                                      # range of years to use from MCMC
		targM=target[,as.character(yrM1)]                                      # target data from MCMC
		yrM2 =setdiff(yrM,yrM1)                                                # missing target years (can occur before and after the MCMC years)
#browser(); return()
		if (length(yrM2)>0) {
			nrow=dim(target)[1]
			if (any(yrM2<yrMr[1])) {
				yrMo =yrM2[yrM2<yrMr[1]]                                         # years of data older than MCMCs
				ncol =length(yrMo)
				targ0=matrix(rep(target[,as.character(yrM1[1])],ncol),
					nrow=nrow, ncol=ncol, dimnames=list(1:nrow,yrMo))               # repeat B0 (first column)
				targM=cbind(as.data.frame(targ0),targM)                          # moving target
			}
			if (any(yrM2>yrMr[2])) {
				yrMn =yrM2[yrM2>yrMr[2]]                                         # years of data newer than MCMCs
				ncol =length(yrMn)
				targN=vmat[,as.character(yrMn)]                                  # start using projections
				targM=cbind(targM,targN)                                         # moving target
			}
		}
		rats=vmat/targM                                                        # matrix of ratios Bt/ moving target
	}
	else    # if it's a vector, so no moving window
		rats=apply(vmat,2,function(x,targ){x/targ},targ=target)                # matrix of ratios Bt/ target (B0 or Bmsy)
	p.hi=apply(rats,2,function(x,r){xhi=x>r; sum(xhi)/length(xhi)},r=ratio)  # vector of probabilities Bt/B0 > target ratio for each year.
	# p.hi can become each row of a decision table (AME checked
	#  the numbers for 0.4 Bmsy match my existing
	#  independent calculations). Need to save this for moving window.
#browser(); return()
	z.hi=p.hi >= conf                                                         # logical: is p.hi >= confidence limit specified

	if (all(z.hi))       yrT=yr0                      # all p.hi exceed the confidence level
	else if (!any(z.hi)) yrT=yrN                      # no  p.hi exceed the confidence level
	else {
		pdif=diff(p.hi)                                # one-year change in trend
		z1=diff(p.hi)>0                                # logical: trend increasing?
		z2=c(pdif[-1],FALSE)>0                         # logical: trend one period later increasing?
		z3=z.hi[-1]                                    # does the probability of equalling or exceeding the target ratio exceed the confidence level?
		z =z1 & z2 & z3                                # logical: potential years when target reached
		if (!any(z)) yrT=yrN                           # target not reached within the projection period
		else         yrT=as.numeric(names(z)[z][1])    # first year when target reached
	}
	N=yrT - yr0                                       # number of years to reach target
	if (plotit) {
		par(mar=c(4,5,0.5,0.5))
		#ylim=c(0, max(p.hi,ratio))
		ylim=c(min(p.hi,0.5),1)
		plot(yr0:yrN,p.hi,type="n",ylim=ylim,ylab="",mgp=c(2.0,0.75,0))
		lines(yr0:yrN,p.hi,col="grey")
		points(yr0:yrN,p.hi,pch=20,col="orange",cex=1.2)
		mtext(text=expression(p~~frac(B[t],B[Target]) ), side=2, line=1.5, cex=1.5)
		abline(h=conf,v=yrT,col="dodgerblue")
	}
	eval(parse(text=paste("return(",retVal,")",sep=""))) 
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~findTarget


#importProjRec--------------------------2013-02-13
# Imports the projected recruitments
#  (actually what's saved is the N(0,1) random numbers,
#   which for a particular MCMC sample are the same for all the catch strategies),
#   which 'importProj' does not do. Need this for YMR. 4th July 2011.
# Also importing vulnerable biomass, for Jaclyn's requested projected
#  exploitation rate decision tables for POP 3CD and 5DE, 2013.
# importProjRecAndy.r - extending importProjRec to include VB, to
#  then calculate projected exploitation rates. 13th Feb 2013
#----------------------------------------------AME
importProjRec=function (dir, info="", coda=FALSE, ngear=1, quiet=TRUE) 
{
	get.Policies <- function() {
		if (!quiet) cat("Policies  ")
		Policies <- read.table(paste(dir, "strategy.out", sep="/"), skip=1)
		if (!quiet) cat("file...")
		Policies <- unique(as.vector(as.matrix(Policies)))
		if (!quiet) cat("unique...OK\n")
		return(Policies)
	}
	get.Years <- function() {
		if (!quiet) cat("Years...")
		Years <- read.table(paste(dir, "strategy.out", sep="/"), nrows=1)
		if (!quiet) cat("file...")
		Years <- unlist(strsplit(as.matrix(Years), "_"))
		if (!quiet) cat("labels...")
		Years <- unique(matrix(Years, nrow=3)[2, ])
		if (!quiet) cat("unique...OK\n")
		return(Years)
	}
	get.B <- function(Policies, Years) {
		if (!quiet) cat("Biomass   ")
		B <- read.table(paste(dir, "projspbm.out", sep="/"), header=TRUE)[, -c(1, 2)]
		if (!quiet) cat("file...")
		Blist <- list()
		for (p in 1:length(Policies)) {
			from <- (p - 1) * length(Years) + 1
			to <- p * length(Years)
			Blist[[p]] <- B[, from:to]
			names(Blist[[p]]) <- Years
		}
		names(Blist) <- Policies
		B <- Blist
		if (!quiet) cat("list...OK\n")
		return(B)
	}
	get.Y <- function(Policies, Years) {
		if (!quiet) cat("Landings...")
		Y <- read.table(paste(dir, "procatch.out", sep="/"), header=TRUE)
		if (!quiet) cat("file...")
		Ylist <- list()
		for (p in 1:length(Policies)) {
			from <- (p - 1) * length(Years) + 1
			to <- p * length(Years)
			Ylist[[p]] <- Y[, from:to]
			names(Ylist[[p]]) <- Years
		}
		names(Ylist) <- Policies
		Y <- Ylist
		if (!quiet) cat("list...OK\n")
		return(Y)
	}
	# AME adding to load in projected recruitment residuals.
	#  Bit different to the others as the file is saved
	#  differently. NOTE that this returns the 'tempdev' values
	#  from Awatea, which are just N(0,1) values, without
	#  multiplying by the sigma_R. Call this epsilon_t
	#  and maybe multiply here by sigmaR, which is
	#	obj$extra$residuals$p_log_RecDev[6]
	get.eps <- function(Policies, Years) {
		if (!quiet) cat("Recruitment...")
		eps <- read.table(paste(dir, "RecRes.out", sep="/"), header=FALSE, skip=1)
		if (!quiet) cat("file...")
		nRow=dim(eps)[1] / length(Policies)
		epslist <- list()
		for (p in 1:length(Policies)) {
			rows <- (0:(nRow-1)) * length(Policies) + p 
			epslist[[p]] <- sigmaR * eps[rows, ]
			names(epslist[[p]]) <- Years
		}
		names(epslist) <- Policies
		eps <- epslist
		if (!quiet) cat("list...OK\n")
		return(eps)
	}
	# AME adding to load in Vulnerable Biomass, editing get.B
	get.VB <- function(Policies, Years) {
		if (!quiet) cat("Vulnerable Biomass   ")
		VB <- read.table(paste(dir, "Projbiom.out", sep="/"), header=TRUE)[, -c(1:(2*ngear))]  # RH added `ngear' to deal with multiple Virgin VBs
		if (!quiet) cat("file...")
		VBlist <- list()
		for (p in 1:length(Policies)) {
			from <- (p - 1) * length(Years) + 1
			to <- p * length(Years)
			VBlist[[p]] <- VB[, from:to]
			names(VBlist[[p]]) <- Years
		}
		names(VBlist) <- Policies
		VB <- VBlist
		if (!quiet) cat("list...OK\n")
		return(VB)
	}
	files <- paste(dir, c("strategy.out", "projspbm.out", "procatch.out", "recres.out", "projbiom.out"), sep="/")
	sapply(files, function(f) if (!file.exists(f)) 
		stop("File ", f, " does not exist. Please check the 'dir' argument.", call.=FALSE))
	if (!quiet) cat("\nParsing files in directory ", dir, ":\n\n", sep="")
	Policies <- get.Policies()
	Years <- get.Years()
	B <- get.B(Policies, Years)
	Y <- get.Y(Policies, Years)
	eps <- get.eps(Policies, Years)
      	VB <- get.VB(Policies, Years)
	if (!quiet) cat("\n")
	output <- list(B=B, Y=Y, eps=eps, VB=VB)
	## Deprecate the use of coda's mcmc function
	#if (coda) {
	#	#eval(parse(text="require(coda, quietly=TRUE, warn.conflicts=FALSE)"))
	#	output <- lapply(output, function(x) lapply(x, mcmc))
	#}
	attr(output, "call") <- match.call()
	attr(output, "info") <- info
	return(output)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~importProjRec


#srFun----------------------------------2011-08-31
# Stock-recruitment function. From ProjRecCalcs.r
#  To input a vector of spawners in year t-1 and calculate recruits in year t.
#  Output for recruits is vector, each element corresponds to spawners the
#  the year before, so will usually want to shift the output by 1 so that 
#  recruits in year t are based on spawners in year t-1.
#  Can also have each input as a vector (used when calculating a single 
#  year but multiple MCMCs, as in first year of projections is based on 
#  penultimate year of MCMC calcualtions.
#----------------------------------------------AME
srFun=function(spawners, h=h.mpd, R0=R0.mpd, B0=B0.mpd) {
# to input a vector of spawners in year t-1 and calculate recruits in year t 
	4 * h * R0 * spawners / ( ( 1 - h) * B0 + (5 * h - 1) * spawners)
}


