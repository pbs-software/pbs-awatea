##==============================================================================
## PBSawatea utility functions:
##  combGear         : combine catches by year from multiple gear types.
##  cquantile        : cumulative quantile, from cumuplot
##  cquantile.vec    : get one probability at a time
##  dfoAxis          : create a DFO axis? whatever (ad hoc)
##  dropLast         : drop the last value in a vector? pulease (ad hoc)
##  findTarget       : derive decision tables for MSY-based ref.pts. and for moving windows
##  getNpan          : get panel number when inside a multi-panel plot.
##  importCor        : import Awatea parameter correlations.
##  importEva        : import Awatea Hessian eigenvalues.
##  importLik        : import Awatea likelihoods.
##  importPar        : import Awatea parameters (all).
##  importRes        : import Awatea results.
##  importStd        : import Awatea output parameter standard deviations.
##  MAfun            : mean age function (Chris Francis, 2011, weighting assumption T3.4, p.1137)
##  makeCmat         : make a 1-column matrix
##  makeRmat         : make a -row matrix
##  makeErrMat       : make simple ageing error matrix for Awatea.
##  redo.currentMCMC : re-run the R components of 'run-masterMCMC.Snw' to get MCMC binary objects
##  repeatMPD        : repeat MPDs for likelihood profiles.
##  splitGear        : split catches, U, and VB by year from multiple gear types.
##  tabSAR           : generate comma-del., 2-D tables from reference point objects.
##  tex.that vec     : convert a vector to a phrase 'x, y, and z' (see texThatVec in PBStools for more advanced version)
##==============================================================================


## combGear ----------------------------2019-05-01
## Combine catches by year from multiple gear types.
## Specifically to collapse VB output from Awatea.
## ---------------------------------------------RH
combGear = function(dat, fn=function(x){sum(x,na.rm=TRUE)})
{
	ayrs = as.numeric(substring(colnames(dat),1,4))
	#aval = cut(ayrs,breaks=(min(ayrs)-1):max(ayrs),labels=FALSE)
	yrs  = unique(ayrs); #names(yrs) = unique(aval)
	if (length(ayrs)==length(yrs)) {
		out = dat
		dimnames(out) = list(rownames(dat), yrs)
	} else {
		out = array(0, dim=c(nrow(dat),length(yrs)), dimnames=list(rownames(dat), yrs))
		for (i in yrs) {
			ii = as.character(i)
			jj = grep(ii, colnames(dat))
			out[,ii] = apply(dat[,jj,drop=FALSE],1,fn)
		}
	}
	return(out)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~combGear


#cquantile------------------------------2010-10-20
# cumulative quantile, from cumuplot
#----------------------------------------------AME
cquantile <- function(z, probs)  
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
      width = 1, digits = 7), "%", sep = "")
  return(cquant)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~cquantile


#cquantile.vec--------------------------2010-10-20
# AME doing this, just do one prob at a time 
# (so it returns a vector not a matrix)
#----------------------------------------------AME
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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~cquantile.vec


## dfoAxis------------------------------2022-01-17
##  Create a DFO axis? Whatever.
##  Ad hoc code originally in PBStools.
## ---------------------------------------------RH
dfoAxis = function(side, vec, prng=0.10)
{
	vec.len = length(vec)
	new.vec = vec
	names(new.vec) = gsub("^\\s+|$\\s+", "", format(new.vec, scientific=FALSE, big.mark=options()$big.mark))
	end.pos   = c(ifelse(side %in% c(1,3), 1, 3), ifelse(side %in% c(1,3), 2, 4))
	dif.axs  = diff(par()$usr[end.pos])
	for (i in 1:length(end.pos)) {
		ii = end.pos[i]
		iii = switch (i, 1, vec.len)
		targ = par()$usr[ii]
		#zone = targ + c(-1,1) * targ * prng
		zone = targ + c(-1,1) * dif.axs * prng
		val  = vec[iii]  ## start or end value
		do.drop = val >= zone[1] && val <=zone[2]
		if (!do.drop || round(val,5)==0) next
		names(new.vec)[iii] = ""
	}
	return(new.vec)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~dfoAxis


## dropLast-----------------------------2022-01-14
##  Drop the last value in a vector? Pulease.
##  Ad hoc code originally in PBStools.
## ---------------------------------------------RH
dropLast = function(vec, target, prng=0.10)
{
	end.zone  = target + c(-1,1) * target * prng
	fin.val   = rev(vec)[1]
	do.drop   = fin.val >= end.zone[1] && fin.val <=end.zone[2]
	if (do.drop)
		return( rev(rev(vec)[-1]) )
	else
		return(vec)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~dropLast


## findTarget --------------------------2019-12-04
##  To derive decision tables for moving windows and find
##  the times to achieve recovery with given confidence.
##  Note: only seems to be used in RPA situations (RH 181015)
##   Vmat   = matrix of projected B-values (MCMC projections x Year)
##   yrP    = user-specified projection years
##   yrG    = number of years for moving target window (e.g, 90y=3 YMR generations). Might not work for all possibilities.
##   ratio  = recovery target ratio
##   target = recovery target values (e.g., B0, Bmsy).
##          = B0.MCMC for ratios of B0
##          = Bmsy.MCMC for ratios of Bmsy
##          = Bt.MCMC for moving window
##   conf   = confidence level required
##   plotit = logical to plot the probability of Bt/target
##   retVal = character name of object to return
##          = "N", look for the global object "Ttab" (number of years to acheive target)
#           = "p.hi" gives global object "Ptab", a list of decision tables where row is the catch option and column is the year.
##  Values are probabilities of acheiving target.
##  (2012-02-20) 'xhi=x>=r' change to 'xhi=x>r'
## ---------------------------------------------RH
findTarget=function(Vmat, yrU=as.numeric(dimnames(Vmat)[[2]]), yrG=90, 
   ratio=0.5, target=B0.MCMC, conf=0.95, plotit=FALSE, retVal="N", op=">")
{
	## oldpar=par(no.readonly=TRUE);  on.exit(par(oldpar))
	yrA   =as.numeric(dimnames(Vmat)[[2]])   ## years available
	yrP   =sort(intersect(yrA,yrU))          ## years for proj
	yr0   =yrP[1]; yrN=rev(yrP)[1]

	vmat=Vmat[,is.element(dimnames(Vmat)[[2]],as.character(yrP))]             ## include only yrP years
#browser();return()

	## Check if COSEWIC reference criterion
	if (is.data.frame(target) || is.matrix(target)) {
		yrM  = yrP - yrG                                                       ## moving target years
		yrM1 = intersect(as.numeric(dimnames(target)[[2]]),yrM)                ## available target years from MCMC
		if (length(yrM1)==0) {                                                 ## projection not long enough for any overlap with 3 generations
			if (retVal=="N") return(NA)
			else {p.hi=rep(NA,length(yrP)); names(p.hi)=yrP }; return(p.hi) }
		yrMr =range(yrM1)                                                      ## range of years to use from MCMC
		targM=target[,as.character(yrM1)]                                      ## target data from MCMC
		yrM2 =setdiff(yrM,yrM1)                                                ## missing target years (can occur before and after the MCMC years)

		if (length(yrM2)>0) {
			nrow=dim(target)[1]
			if (any(yrM2<yrMr[1])) {
				yrMo =yrM2[yrM2<yrMr[1]]                                         ## years of data older than MCMCs
				ncol =length(yrMo)
				targ0=matrix(rep(target[,as.character(yrM1[1])],ncol),
					nrow=nrow, ncol=ncol, dimnames=list(1:nrow,yrMo))             ## repeat B0 (first column)
				targM=cbind(as.data.frame(targ0),targM)                          ## moving target
			}
			if (any(yrM2>yrMr[2])) {
				yrMn =yrM2[yrM2>yrMr[2]]                                         ## years of data newer than MCMCs
				ncol =length(yrMn)
				targN=vmat[,as.character(yrMn)]                                  ## start using projections
				targM=cbind(targM,targN)                                         ## moving target
			}
		}
		rats=vmat/targM                                                        ## matrix of ratios Bt/ moving target
	}
	else    ## if it's a vector, so no moving window
		rats=apply(vmat,2,function(x,targ){x/targ},targ=target)                ## matrix of ratios Bt/ target (B0 or Bmsy)

	#p.hi=apply(rats,2,function(x,r){xhi=x>r; sum(xhi)/length(xhi)},r=ratio)  ## vector of probabilities Bt/B0 > target ratio for each year.

	## vector of probabilities Bt/B0 op (>|<) target ratio for each year.
	p.hi=apply(rats,2,function(x,r){xhi= eval(call(op,x,r)); sum(xhi)/length(xhi)}, r=ratio)

	## p.hi can become each row of a decision table (AME checked)
	##  the numbers for 0.4 Bmsy match my existing
	##  independent calculations). Need to save this for moving window.

	z.hi = p.hi >= conf                               ## logical: is p.hi >= confidence limit specified

#browser();return()
	if (all(z.hi) ||                                  ## all p.hi exceed the confidence level, also check:
		(all(rats[,1]==1) && all(z.hi[-1]))) yrT=yr0   ## if values in first year of projection = values in target (e.g., Bcurr) before proceeding
	else if (!any(z.hi)) yrT=yrN                      ## no  p.hi exceed the confidence level
	else {
		pdif = round(diff(p.hi),5)                     ## one-year change in trend
		z1 = pdif >= 0                                 ## logical: trend increasing? -- sometimes it just remains flat
		## z2 = c(pdif[-1],rev(pdif)[1]) >= 0          ## logical: trend one period later increasing? (or flat) -- PJS does not like this requirement
		z3 = z.hi[-1]                                  ## does the probability of equalling or exceeding the target ratio exceed the confidence level?
		## z  = z1 & z2 & z3                           ## logical: potential years when target reached
		z  = z1 & z3                                   ## logical: potential years when target reached
		if (!any(z)) yrT=yrN                           ## target not reached within the projection period
		else {
			yrT=as.numeric(names(z)[z][1])              ## first year when target reached
			
		}
	}
	N=yrT - yr0                                       ## number of years to reach target
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
#browser();return()
	eval(parse(text=paste("return(",retVal,")",sep=""))) 
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~findTarget


## getNpan------------------------------2019-05-10
##  Get panel number when inside a multi-panel plot.
## ---------------------------------------------RH
getNpan = function()
{
	mfg=par()$mfg
	mfg[2]+(mfg[1]-1)*mfg[4]
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~getNpan


#importCor------------------------------2012-07-30
# Import Awatea parameter correlations.
#-----------------------------------------------RH
importCor = function(cor.file) {
	cfile = readLines(cor.file)
	out = list(par=cfile)
	header = cfile[1]
	cfile = cfile[-1]
	hessian=as.numeric(strsplit(header,split=" = ")[[1]][2])
	cfile = sub("std dev","std.dev",cfile)
	cfile = gsub("]","",gsub("\\[",".",cfile))
	nI = length(cfile)-1 # number of indices
	for (i in 1:nI) {
		ii = i + 1
		cfile[ii] = paste(c(cfile[ii],rep("NA",nI-i)),collapse=" ")
	}
	writeLines(cfile,"cor.tmp")
	cor = read.table("cor.tmp",header=TRUE,sep="")
	iii = paste("i",pad0(1:nI,n=ceiling(log10(nI))),sep="")
	names(cor)[grep("X",names(cor))] = iii
	row.names(cor) = iii
	zi = upper.tri(cor[iii,iii])
	cor[iii,iii][zi] = t(cor[iii,iii])[zi] # populate upper right triangle with lower left values
#browser();return()
	cor.mat = as.matrix(cor[iii,iii])
	cor.name = cor$name
	cor.value = cor$value
	cor.std.dev = cor$std.dev
	out = list(cfile=cfile,cor=cor,cor.mat=cor.mat,index=iii,cor.name=cor.name,cor.value=cor.value,cor.std.dev=cor.std.dev,hessian_log_determinant=hessian)
	file.remove("cor.tmp")
	return(out) }
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^importCor
#test=importCor("something.cor")


#importEva------------------------------2012-08-10
# Import Awatea Hessian eigenvlaues.
#-----------------------------------------------RH
importEva = function(eva.file)
{
	efile=as.numeric(read.table(eva.file,header=FALSE,sep=""))
	out = list(eva=efile)
	return(out)
}


#importLik------------------------------2012-08-08
# Import Awatea likelihoods.
#-----------------------------------------------RH
importLik = function(lik.file)
{
	lfile = readLines(lik.file)
	out = list(lik=lfile)
	liks = lfile[!is.element(lfile,c("","**Likelihoods**"))]
	liks = gsub("@","",gsub("  "," ",gsub("   "," ",liks)))
	liksplit = strsplit(liks,split=" ")
	liknams = sapply(liksplit,function(x){x[1]})
	likvals = sapply(liksplit,function(x,nf){
		mess=paste("assign(\"",x[1],"\",c(",paste(x[-1],collapse=","),"),envir=sys.frame(which=nf))",sep="")
		eval(parse(text=mess))}, nf=sys.nframe())
	for (i in liknams) out[[i]] = get(i)
	return(out)
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^importLik
#test=importLik("likelihood.dat")


#importPar------------------------------2012-07-27
# Import all Awatea parameters.
#-----------------------------------------------RH
importPar = function(par.file) 
{
	pfile = readLines(par.file)
	out = list(par=pfile)
	header = pfile[1]
	pfile = pfile[-1]
	mess = gsub(" ","_",gsub("  ",";",gsub(" = ","=",sub("# ","",header))))
	mess = sub("Maximum_gradient_component","maxgrad",mess)
	mess = sub("Number_of_parameters","npars",mess)
	mess = sub("Objective_function_value","fval",mess)
	eval(parse(text=mess))
	out = c(out, list(npars=npars,fval=fval,maxgrad=maxgrad))
	vnampos = grep("#",pfile)
	vvalbeg = vnampos + 1
	vvalend = vnampos+c(diff(vnampos)-1,length(pfile)-rev(vnampos)[1])
	pfile = gsub("]","",gsub("\\[",".",pfile))
	pfile = gsub(":","",gsub("# ","",pfile))
	pfile = gsub(" ",",",gsub("^ ","",pfile))
	vals  = sapply(1:length(vnampos),function(x,f,n,v1,v2){
		#mess = paste("out[\"",f[n[x]],"\"]= c(", paste(f[v1[x]:v2[x]],collapse=","),")",sep="")
		mess = paste(f[n[x]]," = c(", paste(f[v1[x]:v2[x]],collapse=","),")",sep="")
		eval(parse(text=mess))
	},f=pfile,n=vnampos,v1=vvalbeg,v2=vvalend)
	names(vals)=pfile[vnampos]
	out=c(out,vals)
#browser();return()
	return(out)
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^importPar
#test=importPar("something.par")


## importRes----------------------------2021-03-03
## Awatea res file has following structure (some elements may be
## missing dependent on model configuration and importCol details.
##
## N predicted numbers at age
## B predicted biomass, recruitment, and observed landings
## Sel predicted selectivity and observed maturity (age things)
## Dev predicted recruitment deviates from the stock-recruitment curve
## CPUE, Survey commercial and survey abundance index and fit
## CAc, CAs commercial and survey C@A (catch at age) and fit
## CLc, CLs commercial and survey C@L (catch at length) and fit
## LA observed L@A and fit
## extra = bits and bobs requested by Andy
##----------------------------------------------RH

importRes <- function (res.file, info="", Dev=FALSE, CPUE=FALSE, 
   Survey=FALSE, CAc=FALSE, CAs=FALSE, CLc=FALSE, CLs=FALSE, 
   LA=FALSE, quiet=TRUE, extra=TRUE, sep=" ")
{
	#---SUBFUNCTIONS-------------------------------
	readVector <- function(keyword, same.line = TRUE, file = res.file,
		vector = res.vector) {
		line <- match(keyword, substring(vector, 1, nchar(keyword)))
		v <- if (same.line)
			as.numeric(scan(file, what = "", skip = line - 1, nlines = 1, quiet = TRUE)[-1])
		else as.numeric(scan(file, what = "", skip = line, nlines = 1, quiet = TRUE))
		if (!quiet) .flash.cat("vector...")
		return(v)
	}
	readMatrix <- function(keyword, nrow, header = FALSE, stripe = c("no",
		"left", "right", "upper", "lower"), file = res.file, vector = res.vector) {
		stripe <- match.arg(stripe)
		line <- match(keyword, substring(vector, 1, nchar(keyword))) + as.numeric(header)
#if (keyword=="methodyearsamplesizesex1a1sex1a2sex1a3" && nrow==29) {browser(); return()}
		m <- scan(file, skip = line, nlines = nrow, quiet = TRUE)
		m <- matrix(m, byrow = TRUE, nrow = nrow)
		m <- switch(stripe, left = m[, seq(1, ncol(m)/2)], right = m[,
			seq(ncol(m)/2 + 1, ncol(m))], upper = m[seq(1, nrow(m) -
			1, by = 2), ], lower = m[seq(2, nrow(m), by = 2),], m)
		if (!quiet) .flash.cat("matrix...")
		return(m)
	}
	getN <- function(sexes, years, ages) {
		if (!quiet) .flash.cat("N         ")
		nsexes <- length(sexes)
		nyears <- length(years)
		nages <- length(ages)
		if (nsexes == 1) {
			Nu <- readMatrix("Numbers_at_age_by_Year,sex_and_age", nrow = nyears * nsexes)
			N <- data.frame(Sex=rep(sexes, nyears * nages), Year=rep(years, each = nages), 
				Age=rep(ages,nyears), N=as.vector(t(Nu)))
		}
		if (nsexes == 2) {
			Nf <- readMatrix("Numbers_at_age_by_Year,sex_and_age", nrow = nyears * nsexes, stripe = "upper")
			Nm <- readMatrix("Numbers_at_age_by_Year,sex_and_age", nrow = nyears * nsexes, stripe = "lower")
			N <- data.frame(Sex = rep(sexes, each = nyears * nages), Year = rep(rep(years, each = nages),2), 
				Age = rep(ages, 2 * nyears), N = as.vector(t(rbind(Nf,Nm))))
		}
		if (!quiet) .flash.cat("OK\n")
		return(N)
	}
	getB <- function(years, gears) {
		ngears <- length(gears)
		if (!quiet) .flash.cat("B         ")
		vb <- readMatrix("Vulnerable_Biomass_by_Method_and_Year", nrow = ngears)
		sb <- readVector("Spawning_Biomass_by_Year", same.line = FALSE)
		## *** ADD in the exploitation rate, last year is missing and need
		##     to pad the matrix with missing values to match other "B" columns.
		U <- readMatrix("Exploitation_Rate_by_Method_and_Year", nrow=ngears )
		U <- cbind( U, rep(NA,nrow(U)) )
		#y <- c(readVector("Total_Catch_by_Method_and_Year", same.line = FALSE), NA)
		## BUG FIX: Appears should call readMatrix to accommodate multiple gear series.
		##          Then, sum over the gears to get total catch.
		y <- readMatrix( "Total_Catch_by_Method_and_Year", nrow=ngears )
		#y <- apply( y,2,sum,na.rm=TRUE )
		#y <- c( y,NA )
		y <- cbind( y, rep(NA,nrow(y)) )

		B <- as.data.frame( cbind(years, t(vb), sb, t(y), t(U)) )
		names(B) <- if (ngears == 1) c("Year", "VB", "SB", "Y", "U")
		else c("Year", paste0("VB.",gears), "SB", paste0("Y.",gears), paste0("U.",gears) )
		if (!quiet) .flash.cat("OK\n")
		return(B)
	}
	getSel <- function(gears, surveys, years, sexes, ages) {
		if (!quiet)
			.flash.cat("Sel       ")
		ngears <- length(gears)
		nsurveys <- length(surveys)
		nyears <- length(years)
		nsexes <- length(sexes)
		nages <- length(ages)
		com <- readMatrix("Commercial_age-specific_selectivity_by_method,Year,sex_and_age", nrow=ngears*nyears*nsexes)
		com <- com[seq(1, to = ngears * nyears * nsexes, by = nyears),]
		srv <- readMatrix("Survey_age-specific_selectivity_by_survey,Year,sex_and_age", nrow = nsurveys * nsexes)
		fecundity <- readVector("Fecundity_by_year_and_age", same.line = FALSE)
		weight <- readVector("Weight_by_year,sex_and_age", same.line = FALSE)
		mat <- rep(ifelse(weight > 0, fecundity/weight, 0), nsexes)
		if (is.numeric(gears))
			gears <- paste("Gear", gears)
		if (is.numeric(surveys))
			surveys <- paste("Survey", surveys)
		Sel <- data.frame(Series = c(rep(gears, each = nsexes * nages), 
			rep(surveys, each = nsexes * nages), rep("Maturity", nsexes * nages)), 
			Sex = rep(rep(sexes, each = nages), ngears + nsurveys + 1), 
			Age = rep(ages, (ngears + nsurveys + 1) * nsexes), P = c(t(com), t(srv), mat))
		if (!quiet) .flash.cat("OK\n")
		return(Sel)
	}
	getDev <- function(ages, years) {
		if (!quiet) .flash.cat("Dev       ")
		Dev <- list()
		Dev$Initial <- readVector("log_InitialDev", same.line = TRUE)
		names(Dev$Initial) <- ages[-c(1, length(ages))]
		Dev$Annual <- readVector("log_RecDev", same.line = TRUE)
		names(Dev$Annual) <- years[-length(years)]
		if (!quiet) .flash.cat("OK\n")
		return(Dev)
	}
	getCPUE <- function(gears, years) {
		if (!quiet) .flash.cat("CPUE      ")
		nseries <- readVector("NCPUEindex")
		ngears <- length(gears)
		nyears <- length(years)
		obs <- readMatrix("indexmethodyearvaluecv", nrow = readVector("Number_of_CPUE_data", same.line = FALSE))
		obs <- data.frame(Series = obs[, 1], Gear = obs[, 2], Year = obs[, 3], Obs = obs[, 4], CV = obs[, 5])
		fit <- readMatrix("CPUE_Index_Trajectories", nrow = nseries)
		fit <- data.frame(Series = rep(1:nseries, each = nyears), Year = rep(years, nseries), Fit = as.vector(t(fit)))
		CPUE <- merge(obs[, names(obs) != "Gear"], fit, all = TRUE)
		sgkey <- unique(obs[, c("Series", "Gear")])
		CPUE <- merge(sgkey, CPUE)
		CPUE <- data.frame(Series = paste("Series ", CPUE$Series, "-", CPUE$Gear, sep = ""), 
			Year = as.integer(CPUE$Year), Obs = CPUE$Obs, CV = CPUE$CV, Fit = CPUE$Fit)
		if (!quiet) .flash.cat("OK\n")
		return(CPUE)
	}
	getSurvey <- function(years) {
		if (!quiet) .flash.cat("Survey    ")
		nyears <- length(years)
		nseries <- readVector("Nsurveyindex")
		obs <- readMatrix("indexyearvaluecv", nrow = readVector("Number_of_survey_data", same.line = FALSE))
		obs <- data.frame(Series = obs[, 1], Year = obs[, 2], Obs = obs[, 3], CV = obs[, 4])
		fit <- readMatrix("Survey_Index_Trajectories", nrow = nseries)
		fit <- data.frame(Series = rep(1:nseries, each = nyears), Year = rep(years, nseries), Fit = as.vector(t(fit)))
		Survey <- merge(obs, fit, all = TRUE)
		Survey$Series <- as.integer(Survey$Series)
		Survey$Year <- as.integer(Survey$Year)
		if (!quiet) .flash.cat("OK\n")
		return(Survey)
	}
	getCAc <- function(sexes, ages) {
		if (!quiet) .flash.cat("CAc       ")
		nsexes <- length(sexes)
		nages <- length(ages)
		nobs <- readVector("Number_of_Commercial_C@A", same.line=FALSE)
		obs <- readMatrix("methodyearsamplesizesex1a1sex1a2sex1a3", nrow=nobs)                     # "Observed_C@A"  not unique
		fit <- readMatrix("methodyearsamplesizesex1a1sex1a2sex1a3", nrow=nobs, header=2*(nobs+1))  # "Predicted_C@A" not unique
		CAc <- data.frame(Series=rep(obs[,1],each=nsexes*nages), Year=rep(obs[,2],each=nsexes*nages),
				SS=rep(obs[,3],each=nsexes*nages), Sex=rep(rep(sexes,each=nages),nobs),
				startL=rep(obs[,4],each=nsexes*nages),endL=rep(obs[,5],each=nsexes*nages),
				Age=rep(ages,nsexes*nobs), Obs=as.vector(t(as.matrix(obs[,-(1:5)]))),
				Fit=as.vector(t(as.matrix(fit))))
		# loads in okay with next line commented. Gears not
		#  included in getCAs(), and we don't need for POP. AME
		#   CAc$Gear <- as.integer(CAc$Gear)
		CAc$Year <- as.integer(CAc$Year)
		CAc$Age <- as.integer(CAc$Age)
		if (!quiet) .flash.cat("OK\n")
		return(CAc)
	}
	getCAs <- function(sexes, ages) {
		if (!quiet) .flash.cat("CAs       ")
		nsexes <- length(sexes)
		nages <- length(ages)
		nobs <- readVector("Number_of_survey_C@A", same.line = FALSE)
		obs <- readMatrix("surveyyearsamplesizesex1a1sex1a2sex1a3", nrow = nobs)
		fit <- readMatrix("surveyyearsamplesizesex1a1sex1a2sex1a3", nrow = nobs, header = 2 * (nobs + 1))
		CAs <- data.frame(Series=rep(obs[,1],each=nsexes*nages), Year=rep(obs[,2],each=nsexes*nages),
				SS=rep(obs[,3],each=nsexes*nages), Sex=rep(rep(sexes,each=nages),nobs),
				startL=rep(obs[,4],each=nsexes*nages),endL=rep(obs[,5],each=nsexes*nages),
				Age=rep(ages,nsexes*nobs), Obs=as.vector(t(as.matrix(obs[,-(1:5)]))),
				Fit=as.vector(t(as.matrix(fit))))
		CAs$Series <- as.integer(CAs$Series)
		CAs$Year <- as.integer(CAs$Year)
		CAs$Age <- as.integer(CAs$Age)
		if (!quiet) .flash.cat("OK\n")
		return(CAs)
	}
	getCLc <- function(sexes, lengths) {
		if (!quiet) .flash.cat("CLc       ")
		nsexes <- length(sexes)
		nlengths <- length(lengths)
		nobs <- readVector("Number_of_Commercial_C@L", same.line=FALSE)
		obs <- readMatrix("methodyearsamplesizesex1l1sex1l2sex1l3", nrow=nobs)                 # "Observed_C@L"  not unique
		fit <- readMatrix("methodyearsamplesizesex1l1sex1l2sex1l3", nrow=nobs, header=nobs+1)  # "Predicted_C@L" not unique
		CLc <- data.frame(Series=rep(obs[,1],each=nsexes*nlengths), Year=rep(obs[,2],each=nsexes*nlengths),
				SS=rep(obs[,3],each=nsexes*nlengths), Sex=rep(rep(sexes,each=nlengths),nobs),
				startL=rep(obs[,4],each=nsexes*nlengths),endL=rep(obs[,5],each=nsexes*nlengths),
				Length=rep(lengths,nsexes*nobs), Obs=as.vector(t(as.matrix(obs[,-(1:5)]))),
				Fit=as.vector(t(as.matrix(fit))))
		CLc$Series <- as.integer(CLc$Series)
		CLc$Year <- as.integer(CLc$Year)
		CLc$Length <- as.integer(CLc$Length)
		if (!quiet) .flash.cat("OK\n")
		return(CLc)
	}
	getCLs <- function(sexes, lengths) {
		if (!quiet) .flash.cat("CLs       ")
		nsexes <- length(sexes)
		nlengths <- length(lengths)
		nobs <- readVector("Number_of_surveyC@L",same.line=FALSE)
		obs <- readMatrix("surveyyearsamplesizesex1l1sex1l2sex1l3", nrow=nobs)                     # "Observed_C@L"  not unique
		fit <- readMatrix("surveyyearsamplesizesex1l1sex1l2sex1l3", nrow=nobs, header=2*(nobs+1))  # "Predicted_C@L" not unique
		CLs <- data.frame(Series=rep(obs[,1],each=nsexes*nlengths), Year=rep(obs[,2],each=nsexes*nlengths),
				SS=rep(obs[,3],each=nsexes*nlengths), Sex=rep(rep(sexes,each=nlengths),nobs),
				startL=rep(obs[,4],each=nsexes*nlengths),endL=rep(obs[,5],each=nsexes*nlengths),
				Length=rep(lengths,nsexes*nobs), Obs=as.vector(t(as.matrix(obs[,-(1:5)]))),
				Fit=as.vector(t(as.matrix(fit))))
		CLs$Series <- as.integer(CLs$Series)
		CLs$Year <- as.integer(CLs$Year)
		CLs$Length <- as.integer(CLs$Length)
		if (!quiet) .flash.cat("OK\n")
		return(CLs)
	}
	getLA <- function(sexes, ages) {
		if (!quiet) .flash.cat("LA        ")
		nsexes <- length(sexes)
		nages <- length(ages)
		nobs <- readVector("#femalesmales", same.line = FALSE, file = latage.file, vector = latage.vector)
		obs <- readMatrix("VonBertalanfy--Lenght-at-agefit--Likelihood", nrow = sum(nobs), header = 8)
		obs <- data.frame(Sex = rep(sexes, nobs), Age = obs[,1], Obs = obs[,2])
		owarn <- options(warn = -1)
		Linf <- readVector("VonBeratalanfy:Linf")[-(1:3)]
		K <- readVector("VonBeratalanfy:k")[-(1:3)]
		t0 <- readVector("VonBeratalanfy:to")[-(1:3)]
		CV1 <- readVector("cvoftheFitbysex")[-(1:5)]
		CVratio <- readVector("ratioofcv(L_an)/cv(L_a1)oftheFitbysex")[-(1:7)]
		options(owarn)
		sigmaLA <- readVector("#LinearrelationshipofsigmaL@A(1=age;2=length---ignoreifW@Aissupplied)",
			same.line = FALSE, file = txt.file, vector = txt.vector)[1]
		max.age <- c(max(obs$Age[obs$Sex == sexes[1]]), max(obs$Age[obs$Sex==sexes[2]]))
		fit <- data.frame(Sex = rep(sexes, max.age), Age = c(1:max.age[1],1:max.age[2]))
		fit$Fit[fit$Sex == sexes[1]] <- Linf[1] * (1 - exp(-K[1] * (fit$Age[fit$Sex == sexes[1]] - t0[1])))
		fit$Fit[fit$Sex == sexes[2]] <- Linf[2] * (1 - exp(-K[2] * (fit$Age[fit$Sex == sexes[2]] - t0[2])))
		if (sigmaLA == 1) {
			A <- rep(max(ages), 2)
			a <- cbind(fit$Age[fit$Sex == sexes[1]], fit$Age[fit$Sex==sexes[2]])
			fit$CV[fit$Sex == sexes[1]] <- CV1[1] + CV1[1] * (CVratio[1] - 1)/(A[1] - 1) * (a[, 1] - 1)
			fit$CV[fit$Sex == sexes[2]] <- CV1[2] + CV1[2] * (CVratio[2] - 1)/(A[2] - 1) * (a[, 2] - 1)
		}
		if (sigmaLA == 2) {
			L1 <- Linf * (1 - exp(-K * (1 - t0)))
			Ln <- Linf * (1 - exp(-K * (max(ages) - t0)))
			fit$CV[fit$Sex == sexes[1]] <- CV1[1] + CV1[1] * (CVratio[1] - 1)/(Ln[1] - L1[1]) * (fit$Fit[fit$Sex==sexes[1]] - L1[1])
			fit$CV[fit$Sex == sexes[2]] <- CV1[2] + CV1[2] * (CVratio[2] - 1)/(Ln[2] - L1[2]) * (fit$Fit[fit$Sex==sexes[2]] - L1[2])
		}
		LA <- merge(obs, fit, by = c("Sex", "Age"), all = TRUE)
		LA$Age <- as.integer(LA$Age)
		LA$Fit <- LA$Fit
		LA$CV <- LA$CV
		if (!quiet) .flash.cat("OK\n")
		return(LA)
	}
	getExtra = function(resvec,sep=" ") {
		# only for variables with data on same line
		extra = list(
			likelihoods = c(
			"CPUE","Survey_Index","C@A_Commercial","C@A_survey","Prior_penalties"),# Likelihoods
			parameters = c(
			"R0","avgR0","h",                                                      # Parameters (P)
			"M1","M2",                                                             #  (P) Sex_specific
			"Sfullest","SfullDelta","log_varLest","log_varRest",                   #  (P) Method_specific
			"log_qCPUE","log_BetaCPUE",                                            #  (P) CPUE_index_specific
			"log_qsurvey","surveySfull","survey_SfullDelta",
			"log_surveyvarL","log_surveyvarR",                                     #  (P)  Survey_index_specific
			"log_RecDev"),                                                         #  (P) Recruitment_residuals
			priors = c(
			"R0_prior","h_prior",                                                  # Priors (I)
			"M1_prior","M2_prior","Rinit_prior","uinit_prior","p_plusscale",       #  (I) Sex_specific
			"p_Sfullest","p_Sfulldelta","log_varLest_prior","log_varRest_prior",   #  (I) Method_specific
			"errSfull_prior","errvarL_prior","errvarR_prior",                      #  (I) Method_specific_and_annual
			"log_qCPUE_prior","log_BetaCPUE_prior",                                #  (I) CPUE_index_specific
			"qCPUEerr_prior",                                                      #  (I) CPUE_index_specific_and_annual
			"log_qsurvey_prior","surveySfull_prior","p_surveySfulldelta",
			"log_surveyvarL_prior","log_surveyvarR_prior"),                        #  (I) Survey_index_specific
			residuals = c(
			"p_log_InitialDev","p_log_RecDev")                                     # Recruitment_residuals
		)
		Nsexes       = as.numeric(rev(strsplit(resvec[grep("^Nsexes",resvec)],split=sep)[[1]])[1])
		Nages        = as.numeric(rev(strsplit(resvec[grep("^Nages",resvec)],split=sep)[[1]])[1])
		Nmethods     = as.numeric(rev(strsplit(resvec[grep("^Nmethods",resvec)],split=sep)[[1]])[1])
		NCPUEindex   = as.numeric(rev(strsplit(resvec[grep("^NCPUEindex",resvec)],split=sep)[[1]])[1])
		Nsurveyindex = as.numeric(rev(strsplit(resvec[grep("^Nsurveyindex",resvec)],split=sep)[[1]])[1])
		StartYear    = as.numeric(rev(strsplit(resvec[grep("^StartYear",resvec)],split=sep)[[1]])[1])
		EndYear      = as.numeric(rev(strsplit(resvec[grep("^EndYear",resvec)],split=sep)[[1]])[1])
#browser();return()
		glist = list(general=list(Nsexes=Nsexes,Nages=Nages,Nmethods=Nmethods,
			NCPUEindex=NCPUEindex,Nsurveyindex=Nsurveyindex,StartYear=StartYear,EndYear=EndYear))
		elist = sapply(extra,function(x,resvec) { 
			ex = as.list(x); names(ex) = x
#if (x=="Virgin_Vulnerable_Biomass") {browser();return() }
			for (i in ex) {
				expr=paste("index = grep(\"^",i,sep,"\",resvec)",sep="")
				eval(parse(text=expr))
				if (length(index)==1) {
					if (i %in% c("M1_prior","M2_prior","uinit_prior","p_plusscale")) index = index + (1:Nsexes) - 1
					if (i %in% c("p_Sfullest","p_Sfulldelta","log_varLest_prior","log_varRest_prior",
						"errSfull_prior","errvarL_prior","errvarR_prior")) index = index + (1:Nmethods) - 1
					if (i %in% c("log_qCPUE_prior","log_BetaCPUE_prior","qCPUEerr_prior")) index = index + (1:NCPUEindex) - 1
					if (i %in% c("log_qsurvey_prior","surveySfull_prior","p_surveySfulldelta",
						"log_surveyvarL_prior","log_surveyvarR_prior")) index = index + (1:Nsurveyindex) - 1
					exvec = strsplit(resvec[index],split=sep)
					exmat = t(sapply(exvec,function(x){as.numeric(x[!is.element(x,c(i,""))])}))
#if (i=="log_qCPUE_prior") {browser();return() }
					if (!sapply(extra,function(e,i){is.element(i,e)},i=i)["priors"] && nrow(exmat)==1 )
						exres=as.vector(exmat)
						else exres=exmat
					ex[[i]] = exres
				}
				else ex[[i]] = "not found"
			}
			return(ex) }, resvec=resvec, simplify=FALSE)
		# for variables with names on single line followed by chunks of values
		chunk = list(
			biomass = c(
			"Virgin_Vulnerable_Biomass","Virgin_Spawning_Biomass")
		)
		clist = sapply(chunk,function(x,resvec) { 
			ex = as.list(x); names(ex) = x
			for (i in ex) {
				expr=paste("index = grep(\"^",i,"\",resvec)",sep="")
				eval(parse(text=expr))
				if (length(index)==1) {
					if (i %in% c("Virgin_Vulnerable_Biomass","Virgin_Spawning_Biomass")) index = index + 1
					exvec = strsplit(resvec[index],split=sep)
					exmat = t(sapply(exvec,function(x){as.numeric(x[!is.element(x,c(i,""))])}))
					if (nrow(exmat)==1) exres=as.vector(exmat)
					else                exres=exmat
					ex[[i]] = exres
				}
				else ex[[i]] = "not found"
			}
			return(ex) }, resvec=resvec, simplify=FALSE)
	return(c(glist,elist,clist))
	}
	#---END SUBFUNCTIONS---------------------------

	if (!file.exists(res.file))
		stop("File ", res.file, " not found. Use / or \\\\ separators.")
	file.copy(from=res.file, to=paste(res.file,".temp",sep=""), overwrite=TRUE)  ## save temporary version
	file.rename(from=res.file, to=paste(res.file,".original",sep=""))            ## rename original file
	file.copy(from=paste(res.file,".temp",sep=""), to=res.file, overwrite=TRUE)  ## save temporary version
	exitfun = function(resfile) {
		file.rename(from=paste(resfile,".original",sep=""), to=resfile)
		file.remove(paste(resfile,".temp",sep="")) }
	on.exit( exitfun(res.file) )                                                 ## restore original file and remove the temporary copy
	res.vector <- readLines(res.file)
	res.vector <- gsub("[?]","",res.vector,perl=TRUE)
	writeLines(res.vector,con=res.file)

	res.vector <- gsub("\t", " ", gsub(sep, " ", res.vector))
	if (extra)
		extra = getExtra(res.vector,sep=sep)
	else extra = NULL

	res.vector <- gsub("\"", "", gsub("\t", "", gsub(" ", "", res.vector)))
	if (!quiet)
		.flash.cat("\nParsing text file ", res.file, ":\n\nPreamble  ", sep = "")
	sexes   <- if (readVector("Nsexes") == 1) "Unisex"
		else c("Female", "Male")
	gears   <- seq(1, length.out = readVector("Nmethods"))
	surveys <- seq(1, length.out = readVector("Nsurveyindex"))
	years   <- seq(from = readVector("StartYear"), to = readVector("EndYear") + 1)
	ages    <- seq(from = 1, to = readVector("Nages"))
	lengths <- seq(from = readVector("First_length"), by = readVector("Length_class_increment"),
		length.out = readVector("Number_of_length_classes"))
	if (!quiet) .flash.cat("OK\n")
	model   <- list()
	model$N <- getN(sexes, years, ages)
	model$B <- getB(years, gears)
	rec     <- model$N[model$N$Age == 1, ]
	rec     <- tapply(rec$N, rec$Year, sum)
	#model$B$R <- c(rec[-1], NA) # RH: Lagged to represent age 1?
	model$B$R <- rec             # RH: Just use original number as AME reverts to this in `run-Master.Snw`
	model$Sel <- getSel(gears, surveys, years, sexes, ages)
	if (Dev)
		model$Dev <- getDev(ages, years)
	if (CPUE)
		model$CPUE <- getCPUE(gears, years)
	if (Survey)
		model$Survey <- getSurvey(years)
	if (CAc)
		model$CAc <- getCAc(sexes, ages)
	if (CAs)
		model$CAs <- getCAs(sexes, ages)
	if (CLc)
		model$CLc <- getCLc(sexes, lengths)
	if (CLs)
		model$CLs <- getCLs(sexes, lengths)
	if (LA) {
		latage.file <- paste(dirname(res.file), "l_at_age.dat", sep = "/")
		if (!file.exists(latage.file))
			stop("File ", latage.file, " not found. Use / or \\\\ separators.")
		latage.vector <- readLines(latage.file)
		latage.vector <- gsub("\"", "", gsub("\t", "", gsub(" ","", latage.vector)))
		txt.file <- gsub("\\.res", "\\.txt", res.file)
		if (!file.exists(txt.file))
			stop("File ", txt.file, " not found. Use / or \\\\ separators.")
		txt.vector <- readLines(txt.file)
		txt.vector <- gsub("\"", "", gsub("\t", "", gsub(" ", "", txt.vector)))
		model$LA <- getLA(sexes, ages)
	}
	model$extra = extra
	if (!quiet) .flash.cat("\n")
	attr(model, "call") <- match.call()
	attr(model, "scape.version") <- installed.packages()["scape", "Version"]
	attr(model, "info") <- info
	class(model) <- "scape"
	return(model)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~importRes


#importStd------------------------------2012-07-27
# Import Awatea table of estimated parameters.
#-----------------------------------------------RH
importStd = function(std.file, vnam="name") 
{
	sfile = readLines(std.file)
	sfile = sub("std dev","std.dev",sfile)
	sfile = gsub("]","",gsub("\\[",".",sfile))
	writeLines(sfile,"std.tmp")
	std = read.table("std.tmp",header=TRUE,sep="")
	out = list(std=std)
	onam =setdiff(names(std),vnam)
	vars = unique(std[,vnam])
	for (v in vars) {
		mess = paste(v,"=std[is.element(std[,vnam],\"",v,"\"),onam,drop=FALSE]",sep="")
		mess = c(mess, paste("out[\"",v,"\"] = list(",v,"=",v,")",sep=""))
		eval(parse(text=paste(mess,collapse=";")))
	}
	file.remove("std.tmp")
	return(out)
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^importStd
#test=importStd("something.std")


## MAfun--------------------------------2021-02-25
##  Mean age function (Chris Francis, 2011, weighting assumption T3.4, p.1137)
##  MAfun is used in both `runADMB` and in `runSweave`.
##  Francis: j=composition dataset, y=years, b=bins (ages)
## ---------------------------------------------RH
MAfun = function(padata,brks=NULL)
{
	padata=padata[padata$Age>=padata$startL & padata$Age<=padata$endL,]# Choose only ages older than `startL`
	# S = series, y = year, a = age bin, O = observed proportions, P = Predicted (fitted) proportions, SS = N = Sample Size
	S=padata$Series; y=padata$Year; a=padata$Age; O=padata$Obs; E=padata$Fit; SS=padata$SS   # note: SD and NR not used
	if (is.null(brks)) {
		b = paste(S,y,sep="-"); J = unique(S) } ## b = bins = no. age classes * no. sexes
	else {
		B = cut(y, breaks=brks, include.lowest=TRUE, labels=FALSE)
		b = paste(S,B,y,sep="-"); J = unique(paste(S,B,sep="-")) }
	# make sure that input age proportions are standardised (especially if females only) (they should be)
	O    = as.vector(sapply(split(O,b),function(x){x/sum(x)})) # standardise O (obs props)
	#E    = as.vector(sapply(split(E,b),function(x){x/sum(x)})) #standardise E (fitted props)  -- causes instability
	## Note: sum(O) sould equal sum(E) and both should be the number of years

	## From Francis (2011, p.1137) -- Methods allowing for correlations (bottom left of page)
	Oay  = a * O; Oay2 = a^2 * O                 ## Francis: xb*Ojby, xb*Ejby, xb^2 * Ebjy
	Eay  = a * E; Eay2 = a^2 * E                 ## Francis: xb*Ojby, xb*Ejby, xb^2 * Ebjy
	mOy  = sapply(split(Oay,b),sum,na.rm=TRUE)   ## mean observed age by year (*** used in TA1.8 w calc ***)
	mOy2 = sapply(split(Oay2,b),sum,na.rm=TRUE)  ## component used in variance calculation
	Vobs = mOy2 - mOy^2                          ## variance of observed age distribution
	mEy  = sapply(split(Eay,b),sum,na.rm=TRUE)   ## component used in variance calculation (*** used in TA1.8 w calc ***)
	mEy2 = sapply(split(Eay2,b),sum,na.rm=TRUE)  ## weird inflated mean age by year for variance calc
	Vexp = mEy2 - mEy^2                          ## variance of expected age distribution (*** used in TA1.8 w calc ***) v_jy = Sum_b(xb^2 * Ebjy) - mean(Ejy)^2
	N    = sapply(split(SS,b),mean,na.rm=TRUE)   ## sample size (*** used in TA1.8 w calc ***)
	Yr   = as.numeric(substring(names(mOy),nchar(names(mOy))-3))
	m    = sapply(split(SS,b),length)
	CI   = 1.96 * sqrt(Vobs)/sqrt(m)
	#CLlo = mOy-CI; CLhi=mOy+CI
	return(list(MAobs=mOy, MAexp=mEy, Vobs=Vobs, Vexp=Vexp, Yr=Yr, N=N, CI=CI, J=J)) # observed and expected mean ages, variance of expected ages
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MAfun


makeCmat =function(x,colname="Y") {
	matrix(x,ncol=1,dimnames=list(names(x),colname)) }
makeRmat =function(x,rowname="Y") {
	matrix(x,nrow=1,dimnames=list(rowname,names(x))) }

#makeErrMat-----------------------------2011-05-05
# Make a simple ageing error matrix for Awatea.
#-----------------------------------------------RH
makeErrMat = function(N=60, ondiag=0.8, offdiag=0.1, corner=0.9) 
{
	errMat = diag(ondiag,N,N)
	for (i in 1:(N-1))
		errMat[i,i+1] = offdiag
	for (j in 1:(N-1))
		errMat[j+1,j] = offdiag
	errMat[1,1] = errMat[N,N] = corner
	write.table(errMat,file="errmat.dat",sep="\t",row.names=FALSE,col.names=FALSE)
	return(errMat)
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^makeErrMat


## redo.currentMCMC---------------------2023-01-24
## Re-run the R components of 'run-masterMCMC.Snw' to get MCMC binary objects:
##  currentRes, currentMCMC, currentProj, currentMSY, Bmcmc
## ---------------------------------------------RH
redo.currentMCMC = function(strSpp, assYr, stock, mpdir, mcdir, mcsub=201:1200) 
{
	run.dir = dirname(mpdir)
	run.rwt =  sub("^MPD\\.","",basename(mpdir))
	mpd.res = list.files(run.dir, pattern=paste0(".+",run.rwt,"\\.res"))
	model.name = paste0(strSpp, "-", stock, ".", run.rwt)
	spp.name = strSpp

	## Function availability
	d.awatea = "C:/Users/haighr/Files/Projects/R/Develop/PBSawatea/Authors/Rcode/develop/"
	r.awatea = c("importRes","importProjRec","msyCalc","findTarget","refPointsHist","calc.refProbsHist")
	for (i in r.awatea) source(paste0(d.awatea,i,".r"))

	## Object availability
	catch.type = ""
	area.name  = stock
	
	## Awatea Current Results
	## ============================================
	currentRes = importRes(paste0(run.dir,"/",mpd.res), Dev=TRUE, CPUE=TRUE, Survey=TRUE, CLc=TRUE, CLs=TRUE, CAs=TRUE, CAc=TRUE, extra=TRUE)
	#assign("currentRes", currentRes, pos=1)   ## maybe need the assign (maybe because you are jumping from code chunk to code chunk)
	save("currentRes",file=paste0(mcdir,"/currentRes.rda"))  ## useful to have when compiling the final results appendix
	## ============================================

	## Awatea Current MCMC
	## ============================================
	load(paste0(mcdir,"/currentRes.rda"))  ## just for debugging
	currentMCMC <- scape::importMCMC( dir=mcdir, quiet=FALSE )
	#currentMCMC.orig = currentMCMC
	
	## importMCMC (a scape function) seems to get years wrong on the recruitment, see popScape2.r for details, here is the fix
	colnames(currentMCMC$R) = as.integer(colnames(currentMCMC$R)) + 1  ## currentRes$B$R seems one off
	
	Lcpue = currentRes$extra$likelihoods$CPUE
	Ncpue = if (length(Lcpue)==1 && Lcpue==0) 0 else length(Lcpue)
	Ncomm = Ngear = currentRes$extra$general$Nmethods
	Nsurv = currentRes$extra$general$Nsurveyindex
	Sser  = unique(currentRes$Survey$Series)

	if ( Ncpue==0 ){
		useP = !is.element(names(currentMCMC$P),findPat(c("log_qCPUE","LogBetaCPUE"),names(currentMCMC$P)))
		currentMCMC$P = currentMCMC$P[,useP]
	}
	Pnames = names(currentMCMC$P)
	## can take out paste once get updated PBSawatea.
	new.Pnames = gsub("M1_2",paste("M_",2,sep=""),sub("M1_1",paste("M_",1,sep=""),sub("R0","R_0",Pnames)))

	if (Ncpue>0) {
		for (i in 1:Ncpue) {
			ii = i + Nsurv
			new.Pnames = sub(paste("log_qCPUE_",i,sep=""),paste("log q_",ii,sep=""),new.Pnames)
		}
		new.Pnames = sub("log_BetaCPUE","log beta",new.Pnames)
	}
	for (i in 1:Ncomm) {
		ii = i + Nsurv
		new.Pnames = sub(paste("Sfullest_",i,sep=""),paste("mu_",ii,sep=""),new.Pnames)
		new.Pnames = sub(paste("log_varLest_",i,sep=""),paste("log v_",ii,"L",sep=""),new.Pnames)
		#@rmdome  new.Pnames = sub(paste("log_varRest_",i,sep=""),paste("log v_",ii,"R",sep=""),new.Pnames)
		new.Pnames = sub(paste("Sfulldelta_",i,sep=""),paste("Delta_",ii,sep=""),new.Pnames)
	}
	for (i in 1:Nsurv) {
		ii = Sser[i]
		new.Pnames = sub(paste("surveySfull_",i,sep=""),paste("mu_",ii,sep=""),new.Pnames)
		new.Pnames = sub(paste("log_surveyvarL_",i,sep=""),paste("log v_",ii,"L",sep=""),new.Pnames)
		#@rmdome  new.Pnames = sub(paste("log_surveyvarR_",i,sep=""),paste("log v_",ii,"R",sep=""),new.Pnames)
		new.Pnames = sub(paste("surveySfulldeltaest_",i,sep=""),paste("Delta_",ii,sep=""),new.Pnames)
		new.Pnames = sub(paste("log_qsurvey_",i,sep=""),paste("log q_",ii,sep=""),new.Pnames)
	}
	tab.Pnames = c("R_0","M_1","M_2","h", paste("log q_",1:(Nsurv+Ncpue),sep=""))
	if (Ncpue>0) tab.Pnames = c(tab.Pnames, "log beta")
	tab.Pnames = c(tab.Pnames,
	  paste("mu_",1:(Nsurv+Ncomm),sep=""),
	  paste("Delta_",1:(Nsurv+Ncomm),sep=""),
	  paste("log v_",1:(Nsurv+Ncomm),"L",sep="")
		#@rmdome  , paste("log v_",1:(Nsurv+Ncomm),"R",sep="")
	)

	## Put them into the standard order:
	use.Pnames = tab.Pnames[is.element(tab.Pnames,new.Pnames)]
	
	if(length(use.Pnames) != length(new.Pnames))
		{stop("Check tab.Pnames in run-masterMCMC.Snw")}
	
	## Assign new names
	names(currentMCMC$P) = new.Pnames    
	# Re-order to match Paul's results table 3:
	currentMCMC$P = currentMCMC$P[,use.Pnames]
	
	## First value of MCMC = MPD so collect before burn-in removal
	## Can double-check with values in currentRes$extra$parameters
	mpd.P = currentMCMC$P[1,]
	mpd.B = currentMCMC$B[1,]
	mpd.R = currentMCMC$R[1,]

	## Subset MCMC to remove burn-in (e.g., 201:1200)
	#currentMCMC <- sapply(currentMCMC,function(x,s){x[s,]},s=mcsub,simplify=FALSE)
	currentMCMC <- lapply(currentMCMC,function(x,s){x[s,]},s=mcsub)
	#assign( "currentMCMC", currentMCMC, pos=1 ) # you do need the assign (maybe because you are jumping from code chunk to code chunk)
	num.MCMC = dim(currentMCMC$P)[1]   ## number of MCMC samples

	##--- Vulnerable biomass ----------------------------------------------
	## Also have to import vulnerable biomass from vulnBiom.pst, as it's
	##  not done in importMCMC. Can just do as a data table. Has columns
	##  representing years, and each of 1000 rows is an MCMC sample. Same
	##  size as currentMCMC$B. It is
	##  calculated as denominator of (D.11) in QCS POP model appendix.
	
	vbMCMC = read.table(paste0(mcdir,"/vulnBiom.pst"), header=TRUE)
	#names(vbMCMC) = names(currentMCMC$B)      # to make them simply years (only works if Ngear=1)
	names(vbMCMC) = paste(rep(names(currentMCMC$B),Ngear),rep(1:Ngear,each=dim(currentMCMC$B)[2]),sep="_") # standardise names to "year_gear"
	currentMCMC[["VB"]] = vbMCMC[mcsub,]

	## Now going to use q_1, q_2, q_3, not log q_1 etc. So rename and then change values:
	qnames = use.Pnames[grep("log q_",use.Pnames)]
	for (i in qnames) {
		currentMCMC$P[,i] = exp(currentMCMC$P[,i])
		mpd.P[i] = exp(mpd.P[i])  # exponentiates the MCMC values
	}
	names(currentMCMC$P)[is.element(names(currentMCMC$P),qnames)] = substring(qnames,5)
	names(mpd.P)[is.element(names(mpd.P),qnames)] = substring(qnames,5)
	
	## Save the MPD values for use in Appendix F
	currentMCMC$mpd = list(mpd.P=mpd.P, mpd.B=mpd.B, mpd.R=mpd.R)

	#assign( "currentMCMC", currentMCMC, pos=1 ) ## you do need the assign (maybe because you are jumping from code chunk to code chunk)
	## currentMCMC is updated with U on L253, UoverUmsy on L261 and BoverBmsy on L263, so save again later on L270
	save("currentMCMC",file=paste0(mcdir,"/currentMCMC.rda"))  ## useful to have when compiling the final results appendix
	## ============================================

	## Awatea Current Projections
	## ============================================
	prj.dir = paste0(mcdir,"/PRJ.", run.rwt)
	years = currentRes$B[,"Year"]
	sigmaR = currentRes$extra$residuals$p_log_RecDev[6]

	## Assume constant catch (CC) policy:
	if (file.exists(paste0(sub("/$","",prj.dir),"/CC"))) {
		currentProj <- importProjRec( dir=paste0(sub("/$","",prj.dir),"/CC"), ngear=Ngear, sigmaR=sigmaR, quiet=FALSE )
	} else {
		currentProj <- importProjRec( dir=prj.dir, ngear=Ngear, sigmaR=sigmaR, quiet=FALSE )
	}
	## Check for Harvest Rate (HR) policy:
	if (file.exists(paste0(sub("/$","",prj.dir),"/HR"))) {
		currentProj2 <- importProjRec( dir=paste0(sub("/$","",prj.dir),"/HR"), ngear=Ngear, sigmaR=sigmaR, quiet=FALSE )
	} else {
		currentProj2 <- NULL
	}
	## RH added `ngear' to deal with multiple Virgin VBs
	## importProjRec includes epsilons for recruitments
	## Then below we calcualte the actual recruitments.
	## RH: Note that only one VB series is projected because input file specifies:
	## Gear used in projections (0=read in; 1=use same gear proportions as last year of data)

	currentProj <- lapply(currentProj,function(X){lapply(X,function(x,s){x[s,]},s=mcsub)})
	if (!is.null(currentProj2))
		currentProj2 <- lapply(currentProj2,function(X){lapply(X,function(x,s){x[s,]},s=mcsub)})

	## **** COMMENTING THIS OUT TO GET POP 5ABC WORKING
	## Rowan, can we put an 'if' statement in to ignore this only for POP 5ABC.
	## Take off final year of projection
	if (!(strSpp=="POP" && assYr==2012)) { 
		currentProj$B = lapply(currentProj$B, function(x) {  x[ ,1:(dim(x)[[2]]-1)]  })
		currentProj$Y = lapply(currentProj$Y, function(x) { x[ ,1:(dim(x)[[2]]-1)]  })
		currentProj$eps = lapply(currentProj$eps, function(x) {  x[ ,1:(dim(x)[[2]]-1)] })
		currentProj$VB = lapply(currentProj$VB, function(x) {  x[ ,1:(dim(x)[[2]]-1)]  })
		if (!is.null(currentProj2)) {
			currentProj2$B = lapply(currentProj2$B, function(x) {  x[ ,1:(dim(x)[[2]]-1)]  })
			currentProj2$Y = lapply(currentProj2$Y, function(x) { x[ ,1:(dim(x)[[2]]-1)]  })
			currentProj2$eps = lapply(currentProj2$eps, function(x) {  x[ ,1:(dim(x)[[2]]-1)] })
			currentProj2$VB = lapply(currentProj2$VB, function(x) {  x[ ,1:(dim(x)[[2]]-1)]  })
		}
	}
	## Calculate projected exploitation rates.
	currentProj$U = currentProj$VB    ## Want the same size
	if (!is.null(currentProj2)) {
		currentProj2$U = currentProj2$VB    ## Want the same size
	}
	catchProj = names(currentProj$VB)
	for(i in catchProj)
	{
		currentProj$U[[i]] = currentProj$Y[[i]] / currentProj$VB[[i]]
	}
	if (!is.null(currentProj2)) {
		catchProj2 = names(currentProj2$VB)
		for(i in catchProj2)
		{
			currentProj2$U[[i]] = currentProj2$Y[[i]] / currentProj2$VB[[i]]
		}
	}
	save("currentProj",file=paste0(mcdir,"/currentProj.rda"))  ## useful to have when compiling the final results appendix
	if (!is.null(currentProj2))
		save("currentProj2",file=paste0(mcdir,"/currentProj2.rda"))  ## useful to have when compiling the final results appendix
	## ============================================

	## Awatea Current Maximum Sustainable Yield
	## ============================================
	msy.dir = paste0(mcdir,"/MSY.", run.rwt)

	currentMSY = msyCalc(dir=msy.dir, error.rep=0, despike=TRUE)  ## RH 200427: Added argument 'despike'
	currentMSY = lapply(currentMSY, function(x,s){x[s]}, s=mcsub) ## mcsub subsets the MCMC samples
	#assign( "currentMSY", currentMSY, pos=1 ) ## Forces it global (or something)
	save("currentMSY",file=paste0(mcdir,"/currentMSY.rda"))  ## useful to have when compiling the final results appendix
	## ============================================

	## Do these for ease of showing statistics in tables. 
	## Each should be a vector with value for each MCMC draw

	if(currentRes$extra$priors$Rinit_prior[1] > 0 & currentRes$extra$priors$Rinit_prior[7] != 1)
		stop("Not starting from unfished equilibrium, so need to fix B0 values")

	## AME changing RH's Year to currYearChar
	currYearChar  = rev(dimnames(currentMCMC$B)[[2]])[1]      
	                # character current year (start for projections)
	## AME changing RH's Year0 to currYear
	currYear      = as.numeric(currYearChar)          ## numeric current year
	prevYear      = currYear - 1   ## previous year for final exploitation
	prevYearChar  = as.character(prevYear)
	startYearChar = dimnames(currentMCMC$B)[[2]][1]   ## start year
	startYear     = as.numeric(startYearChar)

	##-----Gather B0, Bcurr, VB0, VBcurr--------------
	B0.MCMC     = currentMCMC$B[,1,drop=FALSE]
	Bcurr.MCMC  = currentMCMC$B[,currYearChar,drop=FALSE]
	VB0.MCMC    = currentMCMC$VB[,grep(startYearChar,names(currentMCMC$VB)),drop=FALSE]
	VBcurr.MCMC = currentMCMC$VB[,grep(currYearChar,names(currentMCMC$VB)),drop=FALSE]

	## MSY procedure appears to only calculate one VB and U per MCMC sample
	Bmsy.MCMC   = currentMSY$B
	VBmsy.MCMC  = currentMSY$VB # only one VB is calculated in projections
	msy.MCMC    = currentMSY$yield
	umsy.MCMC   = currentMSY$u

	##-----Gather Ucurr-------------------------------
	## Need to calculate exploitation rates over time for MCMC (MPD's are included in `currentRes', but nothing in `currentMCMC'.
	## Going to add currentMCMC$U to currentMCMC. 
	## After doing this realised it was sort of done in popScapeRuns2.r, but only internally for plotting figures.

	catch = currentRes$B[,-1][,grep("Y",names(currentRes$B[,-1])),drop=FALSE]
	if ( !all(is.na(catch[nrow(catch),])) )
		stop("Check catch =   and  currentMCMC$U =     in run-masterMCMC.Snw") 
	catch = catch[-nrow(catch),,drop=FALSE]
	row.names(catch) = years[-length(years)]
	CATCH = apply(catch,1,sum)
	
	currentMCMC$U = currentMCMC$VB[,grep(currYearChar,names(currentMCMC$VB),invert=TRUE)] ## RH need to be tricky when Ngear>1
	names.cmu     = names(currentMCMC$U)
	VBcatch       = unlist(catch);    names(VBcatch)=names.cmu
	currentMCMC$U = as.data.frame(t(apply(currentMCMC$U, 1, function(x,y){ y/x }, y=VBcatch))) ## Need transpose to get right way round again

	## So currentMCMC$U is now exploitation rate for MCMC output.
	Ucurr.MCMC    = currentMCMC$U[,grep(prevYearChar,names(currentMCMC$U)),drop=FALSE]  ## same as AME's upenult.MCMC
	Umax.MCMC     = apply(currentMCMC$U, 1, max)
	## Calculate geomentric mean across gears (for now) as target in Decision Tables (RH 190425)
	UcurrGM.MCMC  = apply(Ucurr.MCMC,1,calcGM)

	currentMCMC$UoverUmsy = as.data.frame(apply(currentMCMC$U, 2, function(x,y){ x/y },  y=umsy.MCMC))    # No transpose
	UoverUmsy.med = apply(currentMCMC$UoverUmsy, 2, median)
	currentMCMC$BoverBmsy = as.data.frame(apply(currentMCMC$B, 2, function(x,y){ x/y },  y=Bmsy.MCMC))    # No transpose, also agrees with table below
	BoverBmsy.med = apply(currentMCMC$BoverBmsy, 2, median)

	##-----Gather MCMC values suggested by Trevor Davies (Dalhousie) for comparison
	trevorMCMC = currentMCMC$P[,c("h","M_1","M_2")[is.element(c("h","M_1","M_2"),names(currentMCMC$P))],drop=FALSE]
	trevorMCMC = data.frame(trevorMCMC, B0=B0.MCMC[,1], MSY=msy.MCMC, Bmsy=Bmsy.MCMC, umsy=umsy.MCMC, f=currentMCMC$L$f)
	currentMCMC[["trevor"]] = trevorMCMC
	save("currentMCMC",file=paste0(mcdir,"/currentMCMC.rda"))  ## useful to have when compiling the final results appendix
	## ============================================

	## To calculate the actual projected recruitments (from the eps)
	currentProj$R = list()
	if (!is.null(currentProj2)) {
		currentProj2$R = list()
	}
	NN = matrix(1:num.MCMC, nrow=1)      # to use to populate each data.frame
	projYearsNames = names(currentProj$B[[1]])
	projYearsNum   = length(projYearsNames)
	
	## First calculate recruits for first projection year, which is based
	##  on penultimate MCMC year's spawners (which is 2010, final year
	##  of that is 2011, first proj year is 2011, and I checked that
	##  values for last year of MCMC equal those for first yr of proj:
	##  > range(currentMCMC$B[, 72] - currentProj$B$'250'[, 1]) 0 0.
	
	## Do this here as does not depend on projections (and so is same for all catch strategies).
	Bpen.MCMC = currentMCMC$B[, rev(names(currentMCMC$B))[2]]
	  ## spawning biomass for penultimate year of MCMC
	
	## Need a vector of h for projections, so if h not estimated
	##  make hForProj just replicate the mpd value:
	
	if (!is.element("h",use.Pnames)) {
		hForProj = rep(currentRes$extra$parameters$h, num.MCMC)
	} else {
		hForProj = currentMCMC$P$h
	}
	RfirstProj = srFun(Bpen.MCMC, R0 = currentMCMC$P$R_0, h=hForProj, B0=B0.MCMC)
	
	currPros = sapply(paste0("currentProj",c("","2")), function(x){eval(parse(text=paste0("!is.null(",x,")")))})
	
	for(i in names(currPros)[currPros]) {
		cP = get(i)
		## Stochastic multiplier, will be same for all strategies as random
		##  numbers currentProj$eps[[j]] are the same for each strategy j
		stochMult = exp(cP$eps$'0' - sigmaR^2/2)
	
		for(j in 1:length(cP$eps) )    ## loop over policies
		{
			junk = apply(NN, 2, function(i, B, R0, h, B0) {
				srFun(as.numeric(B[i,]), h = h[i], R0=R0[i], B0=B0[1] )
			},
			B  = cP$B[[j]],
			R0 = currentMCMC$P$R_0, h = hForProj, B0=B0.MCMC[,1] )  ## RH: I changed B0.MCMC to be a one-column matrix above to be comparable with VB0.MCMC
			## Rowan's trick for using apply on each row.
			junk = t(junk)
			junk = as.data.frame(junk)
			## This gives data frame,
			## rows are MCMC samples, columns are recruits for the next year.
			## Need to insert RfirstProj as first column, and remove final column
			## (which corresponds to recruits for the year after final projection year).
			detR = cbind(RfirstProj, junk)
			detR = detR[, -dim(detR)[2] ]         ## take off final column
			names(detR) = projYearsNames          ## detR is deterministic values
			stochR = detR * stochMult
			cP$R[[j]] = stochR
		}
		names(cP$R) = names(cP$eps)
		assign( i, cP, pos=1 )
		mess = paste0("save(\"", i, "\", file=\"", mcdir, "/", i, ".rda\")")
		eval(parse(text=mess))  ## useful to have when compiling the final results appendix
	}

#browser();return()

	## Refernce Points
	##=============================================
	
	## --------------------------------------------------------------
	## 'BRPprobsList' contains probabilities for the main
	##   decision tables (usually only use MSY-based RefPts).
	## Will include Bmsy and B0 RefPts already calculated above,
	##   e.g., BRPprobsList$'0.4Bmsy' - refProbs$LRP = matrix of 0's
	## --------------------------------------------------------------
	BRPprobsList = BRPprobsList2 = URPprobsList = URPprobsList2 = RPprobTabs = list()
	
	## List of targets:
	Tlst = list(
	  list(ratio=0.4,target=Bmsy.MCMC),
	  list(ratio=0.8,target=Bmsy.MCMC),
	  list(ratio=1.0,target=Bmsy.MCMC),
	  list(ratio=1.0,target=Bcurr.MCMC[,1]),
	  list(ratio=0.2,target=B0.MCMC[,1]),
	  list(ratio=0.4,target=B0.MCMC[,1]),
	  ##---COSEWIC---
	  list(ratio=0.5,target=currentMCMC$B),
	  list(ratio=0.7,target=currentMCMC$B),
	  list(ratio=0.5,target=B0.MCMC[,1]),
	  list(ratio=0.7,target=B0.MCMC[,1])
	)
	names(Tlst)=c("0.4Bmsy","0.8Bmsy","Bmsy","Bcurr","0.2B0","0.4B0","0.5Gen3","0.7Gen3","0.5B0","0.7B0")
	
	## And have to do separately for u_t < u_MSY:
	Ulst = list(
	  list(ratio=1, target=umsy.MCMC),
	  list(ratio=1, target=UcurrGM.MCMC)
	)
	names(Ulst) = c("umsy","ucurr")
	
	## Generation times (y): YMR=30, BOR=20
	## ------------------------------------
	gen1 = ifelse(spp.name %in% c("BOR"), 20, ifelse(spp.name %in% c("YMR"), 30, 20))
	Ngen = 3
	
	for(i in names(currPros)[currPros]) {
		cP = get(i)
		ii = sub("^currentProj","",i)
		#print(ii)
	
		BRPpList = sapply(Tlst,function(x){sapply(cP$B, findTarget, ratio=x$ratio, target=x$target, retVal="p.hi", op=">", yrG=Ngen*gen1)}, simplify=FALSE)
		BRPpList = sapply(BRPpList,t,simplify=FALSE)  ## transpose matrices
		assign( paste0("BRPprobsList",ii), BRPpList )
		RPprobTabs[[paste0("BRPprobsList",ii)]] = BRPpList
	
		## create table for u < umsy
		## Project u must combine both commercial gears using 'tcall(Usplit)'
		## ------------------------------------------------------------------
		URPpList = sapply(Ulst, function(x){sapply(cP$U, findTarget, ratio=x$ratio, target=x$target, retVal="p.hi", op="<", yrG=Ngen*gen1)}, simplify=FALSE)
		URPpList = sapply(URPpList,t,simplify=FALSE)  ## transpose matrices
		assign( paste0("URPprobsList",ii), URPpList )
		RPprobTabs[[paste0("URPprobsList",ii)]] = URPpList
	
		## Also calculate number of years to reach reference target points with a specified confidence
		Ttab.temp = Utab.temp = list()
		for (j in c(0.50,0.65,0.80,0.95)) {
			jj = formatC(j,digits=2,format="f")
			Ttab.temp[[jj]]  = pmin(sapply(Tlst, function(x){sapply(cP$B, findTarget, ratio=x$ratio, target=x$target, conf=j,  retVal="N", op=">", yrG=Ngen*gen1)}), Ngen*gen1)
			Utab.temp[[jj]]  = pmin(sapply(Ulst, function(x){sapply(cP$U, findTarget, ratio=x$ratio, target=x$target, conf=j,  retVal="N", op="<", yrG=Ngen*gen1)}), Ngen*gen1)
		}
		assign( paste0("Ttab.conf",ii), Ttab.temp )
		RPprobTabs[[paste0("Ttab.conf",ii)]] = Ttab.temp
		assign( paste0("Utab.conf",ii), Utab.temp )
		RPprobTabs[[paste0("Utab.conf",ii)]] = Utab.temp
	}
#browser();return()

	## Catch for last num.recentCatchYears years of data:
	num.recentCatchYears = 5
	
	#if( diff(range(rev(catch)[1:2])) > 0.1 ) 
	#   { stop("lastFiveCatch in run-masterMCMC.Snw assumes final year 
	#       equals penultimate, but this isn't the case here so fix it") }
	
	if( diff(range(rev(CATCH)[1:2])) < 0.7 )  # for 5DE the 2012, was set to 2011 rounded 
	{
		recentCatch = rev(rev(CATCH)[2:(2+num.recentCatchYears-1)])
		# The last num.recentCatchYears without the final year as have assumed that's not real data
	} else {
	  recentCatch = rev(rev(CATCH)[1:num.recentCatchYears])
	}
	
	recentCatchMean = mean(recentCatch)
	lab.recentCatchYears = paste(names(recentCatch)[c(1,num.recentCatchYears)],collapse="-")
	refCatSentence = paste0("For reference, the average ", catch.type, " catch over the last ", num.recentCatchYears, " years (", lab.recentCatchYears, ") was ", round(recentCatchMean, digits=0), "~t. ")
	
	maxCatch = max(CATCH)
	maxCatchYear = as.numeric(names(CATCH)[grep(maxCatch,CATCH)[1]])
	if (any(spp.name=="ROL")){
		maxCatSentence = paste("The maximum historical female catch estimate in Area ", area.name, "  was ", round(maxCatch), "~t in ", maxCatchYear,". ",sep="")
	} else { maxCatSentence = "" }
	
	## u.MCMC.med = apply(currentMCMC$U, 2, median)  # median for each year
	## For snail plots:
	
	A = max(currentRes$Sel[,"Age"])
	#T = diff(range(currentRes$B[,"Year"]))+1       # =72
	Ct = currentRes$B$Y[-length(currentRes$B$Y)]    # Takes off catch in final year
	# years = currentRes$B[,"Year"]    # moving earlier
	ages = sort(unique(currentRes$CAc$Age))
	
	selgeqComm = currentRes$Sel[currentRes$Sel$"Series" == "Gear 1",] 
	## comm sel, selgeq4 for POP, presumably woudl be 6 for YMR as 5 survey yrs, so just write Comm
	## mat = currentRes$Sel[currentRes$Sel$"Series" == "Maturity" & currentRes$Sel$"Sex" == "Female",]      # Female maturity
	mat = currentRes$Sel[is.element(currentRes$Sel$"Series","Maturity") & is.element(currentRes$Sel$"Sex",c("Female","Unisex")),]  # Female maturity
	# M1 = currentMCMC$P[1,"M_1"]                 # MPD is first line of MCMC
	# M2 = currentMCMC$P[1,"M_2"]                 # MPD is first line of MCMC
	Rt = currentRes$B$R[-length(currentRes$B$R)]  # Remove final NA for the last year,
	Rt.mpd = Rt       # don't think Rt gets used elsewhere, but leave it valid just in case.
	## For confirmation:
	## R0.mpd = currentMCMC$P[1,"R_0"] #Matches numbers from Ro_So_VB.pst, but wasn't going to use that before?
	# h.mpd = currentMCMC$P[1,"h"]
	Nats.mpd = currentRes$N
	#ut.mpd = currentRes$B$U[-length(currentRes$B$U)]  # remove final NA for the last year
	ut.mpd = currentRes$B[1:(nrow(currentRes$B)-1),grep("U",names(currentRes$B)),drop=FALSE]  # RH: need in case Ngear > 1
	Bt.mpd = currentRes$B[,"SB",drop=FALSE]
	B0.mpd = Bt.mpd[1,1]       #****CHANGE*** if change init cdts.
	Vt.mpd = currentRes$B[,grep("V",names(currentRes$B)),drop=FALSE]  # RH: need in case Ngear > 1
	logRecDev.mpd = currentRes$Dev$Annual
	## Also equals currentRes$extra$parameters$log_RecDev from Rowan's 'extra' sublist  in currentRes.
	
	## AME deleted lots of commented out code. 20th August 2012.

	##-----------------------------------------------------
	## RH adding historical refrerence points for Rock Sole
	##-----------------------------------------------------
	if (is.element(spp.name,c("ROL","rol","621"))) {
		if (area.name=="5CD") HRP.YRS = list(blimYrs=1966:2005, btarYrs=1971:1980, ulimYrs=NULL, utarYrs=1966:2005)
		else                  HRP.YRS = list(blimYrs=1966:2005, btarYrs=1977:1985, ulimYrs=NULL, utarYrs=1966:2005)
	} else {
		HRP.YRS = list(blimYrs=1966:2005, btarYrs=1977:1985, ulimYrs=NULL, utarYrs=1966:2005)
	}
	HRPlist = refPointsHist(mcmcObj=currentMCMC, HRP.YRS=HRP.YRS)  # pass in list of obscure years for calculating Rock Sole historical ref points (B & u MCMCs)
	
	if (is.element(spp.name,c("ROL","rol","621"))) {
		trevorMCMC = data.frame(trevorMCMC, 
			Blim = HRPlist$blimHRP,
			Btar = HRPlist$btarHRP,
			utar = HRPlist$utarHRP ) }
	
	hrp.B = calc.refProbsHist(projObj=currentProj$B, refPlist=HRPlist[c("blimHRP","btarHRP")], op=">", verbose=F) ## default values (function in `PBSscape.r')
	hrp.U = calc.refProbsHist(projObj=currentProj$U, refPlist=HRPlist[c("ulimHRP","utarHRP")], op="<", verbose=F) ## function in `PBSscape.r'
	HRPprobsList = c(hrp.B,hrp.U)
	
	## That's a list, each element is the full table for P> LRP, USR or Bmsy, rows are const catch scenarios and columns are years.
	
	## Adding in tenYear projections for POP 3CD and 5DE 2012.
	## scenarioSubset = c("0", "500", "1000", "1500", "2000" "2250"....
	fiveYears   = as.character(currYear+seq(0,5,1))
	tenYears    = as.character(currYear+seq(0,10,1))
	twentyYears = as.character(currYear+seq(0,20,5))
	ninetyYears = as.character(currYear+seq(0,90,15))
	
	HRPprobs5 = list()    # Probs for 5 years
	HRPprobs10 = list()   # Probs for 10 years
	HRPprobs20 = list()   # Probs for 20 years
	HRPprobs90 = list()   # Probs for 90 years
	for(i in 1:length(HRPprobsList)) {
		if (is.null(HRPprobsList[[i]])) next
		if (projYearsNum>=5) {
			HRPprobs5[[i]] = HRPprobsList[[i]][,fiveYears]
			names(HRPprobs5)[i] = names(HRPprobsList)[i]
		}
		if (projYearsNum>=10) {
			HRPprobs10[[i]] = HRPprobsList[[i]][,tenYears]
			names(HRPprobs10)[i] = names(HRPprobsList)[i]
		}
		if (projYearsNum>=20) {
			HRPprobs20[[i]] = HRPprobsList[[i]][,twentyYears]
			names(HRPprobs20)[i] = names(HRPprobsList)[i]
		}
		if (projYearsNum>=90) {
			HRPprobs90[[i]] = HRPprobsList[[i]][,ninetyYears]
			colnames(HRPprobs90[[i]]) =  as.numeric(colnames(HRPprobs90[[i]])) - currYear
			names(HRPprobs90)[i] = names(HRPprobsList)[i]
		}
	}

	include = FALSE  ## include code has not been checked for function or object availability. (RH 230124)
	if (include) {
	## Moved from above, as now have currentMCMC$U calcs
	## plt.idx( currentRes$Survey,main="Survey Indices") # wasn't called 
	##  in plt.mpdGraphs. Doing it here spits out SD of standardised 
	##  residuals also. May be useful for iterative reweighting?
	
	#meanPolicy = round(mean(recentCatch),-floor(log10(mean(recentCatch))))
	#meanCatch  = mean(recentCatch) # already calculated as `recentCatchMean` on line 903
	projPolicy  = as.numeric(names(currentProj$Y))
	nPolicy     = length(projPolicy)
	## https://stat.ethz.ch/pipermail/r-help/2008-July/167216.html :
	onePos      = which(abs(projPolicy-recentCatchMean)==min(abs(projPolicy-recentCatchMean)))
	nPanel      = 6    ## assume 6 panels in plot
	if (nPolicy <= nPanel) {
		prjPos   = 1:nPolicy
	} else {
		incPos   = floor(nPolicy/nPanel)
		prjPos   = onePos #starting point
		while(length(prjPos)!=(nPanel)) {
			if (min(prjPos)-incPos > 0) prjPos = c(min(prjPos)-incPos,prjPos)
			if (max(prjPos)+incPos <= nPolicy) prjPos = c(prjPos,max(prjPos)+incPos)
			if (length(prjPos)>(nPanel)) prjPos=prjPos[1:(nPanel)]
			.flash.cat(prjPos,"\n")
		}
	}
	## always have 0-catch policy
	if (prjPos[1]==onePos) {
		prjPos = c(1,prjPos[-nPanel])
	} else {
		prjPos[1]  = 1
	}
	if (spp.name=="RSR" && area.name=="5DE")  ## in 2018, recentCatchMean = 109 t so use seq(0,500,100)
		prjPos = 1:nPanel
	onePolicy  = projPolicy[onePos]
	plotPolicy = projPolicy[prjPos]
	
	if(redo.Graphs)
		plt.mcmcGraphs( mcmcObj=currentMCMC, projObj=currentProj, mpdObj=currentRes, save=TRUE, 
			ptypes=tcall(PBSawatea)$ptype, pngres=400, ngear=Ngear,
			plotPolicies = as.character(plotPolicy), #plotPolicies = names(currentProj$Y[1:6]),
			onePolicy    = as.character(onePolicy),  #onePolicy = names(currentProj$Y[2]), 
			mpd=list(mpd.P=mpd.P, mpd.B=mpd.B, mpd.R=mpd.R),
			trevObj = trevorMCMC, lang=lang
		)
	## Change policy options if want other catch policies shown.
	## Set up for length(plotPolicies)=6
	## See help for other options.
	## close.allWin()
	
	## Priors tabulation for priors in MPD table.
	## Must read in a vector of length, and outputs it in the format for the table.
	ptab = function(xx) {
		print(paste0(c(xx[1], " & [", xx[2], ",", xx[3], "] & ", xx[4], " & [", xx[5], ",", xx[6], "] & ", xx[7]), collapse="")) 
	}
	## Quantile tabulation summary using decimal places
	qtab = function(xx.MCMC, dig=0) {  ## dig is number of dec places
		print(paste0( c( prettyNum(round(quantile(xx.MCMC, quants3[1]), digits=dig), big.mark=options()$big.mark),
			" & ", prettyNum(round(quantile(xx.MCMC, quants3[2]), digits=dig), big.mark=options()$big.mark),
			" & ", prettyNum(round(quantile(xx.MCMC, quants3[3]), digits=dig), big.mark=options()$big.mark)), collapse=""))
	}
	## Quantile tabulation summary using significant digits
	stab = function(xx.MCMC, dig=3) {  ## dig is number sig digits
		print(paste0( c( prettyNum(signif(quantile(xx.MCMC, quants3[1]), digits=dig), big.mark=options()$big.mark), 
			" & ", prettyNum(signif(quantile(xx.MCMC, quants3[2]), digits=dig), big.mark=options()$big.mark),
			" & ", prettyNum(signif(quantile(xx.MCMC, quants3[3]), digits=dig), big.mark=options()$big.mark)), collapse=""))
	}
	## to give median (5%-95%) to put in text.
	med5.95 = function(xx.MCMC, dig=0){  ## dig is number of dec places
		print(paste0( c( prettyNum(round(quantile(xx.MCMC, quants3[2]), digits=dig), big.mark=options()$big.mark),
			"~(", prettyNum(round(quantile(xx.MCMC, quants3[1]), digits=dig), big.mark=options()$big.mark), 
			"-", prettyNum(round(quantile(xx.MCMC, quants3[3]), digits=dig), big.mark=options()$big.mark), ")"), collapse=""))
	}
	
	## For q_i in table, values can vary between runs, so set to number
	##  of significant digits. NOT USED, just doing 4 decimal places for all. 
	## sapply(signif(x,3), sprintf, fmt="%#.3g")   # try that for 
	##  1.3001  to be 1.300.   Haven't played with yet.
	
	## not saving for YMR for now (don't have all these variables, though
	## just MPDs so don't need MCMC output).  [These were from MCMC's]
	## save(A, T, Ct, selgeqComm, mat, M1, M2, Rt, R0.mpd, h.mpd, Nats.mpd, ut.mpd,  Bt.mpd, B0.mpd, Vt.mpd, logRecDev.mpd, file="run23values.RData")
	## save.image(file="run23all.RData")
	
	## See popScape2.r for pairs plots, from:
	## Copy and run this for pairs plots        to
	## text(currentRes$B$SB, currentRes$B$U, 1:72)
	
	## Want to report the mean of the median recruitments for past and projections.
	recMed = apply(currentMCMC$R, 2, median)
	meanRecMed = mean(recMed)
	
	muMed = apply(currentMCMC$P[grep("mu_",names(currentMCMC$P))],2,median)
	
	## Do for one policy
	## recProjMed1500 = apply(currentProj$R$'1500', 2, median)
	## meanRecProjMed1500 = mean(recProjMed1500)
	
	## Want to report the year of first age data:
	CAfirstYear = min(min(currentRes$CAc$Year), min(currentRes$CAs$Year))
	
	if (resdoc) {
		if (is.element(spp.name,c("SGR","POP","YMR"))) {
			mpdfigs = c("survIndSer2","meanage","commAgeResids","stockRecruit","recruits","selectivity","exploit")
			mpdfigs = paste(mpdfigs,".eps",sep="")
			mpdfigs = c(mpdfigs,list.files(mpd.dir,pattern=c("ageComm")))
			mpdfigs = c(mpdfigs,list.files(mpd.dir,pattern=c("ageSurv")))
			mpdfigs = c(mpdfigs,list.files(mpd.dir,pattern=c("survRes")))
			mpdfigs = c(mpdfigs,list.files(mpd.dir,pattern=c("survAgeResSer")))
			mpdfigs = c(mpdfigs,list.files(mpd.dir,pattern=c("recDev")))
		}
		if (is.element(spp.name,c("ROL"))) {
			mpdfigs = c("survIndSer2","CPUEser","commAgeResids","stockRecruit")
			mpdfigs = paste(mpdfigs,".eps",sep="")
			mpdfigs = c(mpdfigs,list.files(mpd.dir,pattern=c("ageComm")))
			mpdfigs = c(mpdfigs,list.files(mpd.dir,pattern=c("ageSurv")))
			mpdfigs = c(mpdfigs,list.files(mpd.dir,pattern=c("survRes")))
			mpdfigs = c(mpdfigs,list.files(mpd.dir,pattern=c("survAgeResSer")))
			mpdfigs = c(mpdfigs,list.files(mpd.dir,pattern=c("recDev")))
		}
		## There's no longer a need to copy these figures as the Model Results Appendix
		##   can navigate to various directories using \graphicspath. For example:
		## \graphicspath{{C:/Users/haighr/Files/GFish/PSARC17/POP/Data/Awatea/5ABC/POPrun08/MPD.08.03/}}
		## \input{"POPrun08-3"}
		## file.copy(from=paste(mpd.dir,mpdfigs,sep="/"),to=fig.dir,overwrite=TRUE)
	}

	if (length(muMed)>0) {
		muCI = findPat((Nsurv+1):(Nsurv+Ncomm),names(muMed))
		names(muCI) = sapply(strsplit(muCI,"_"),function(x){x[2]})
		muSI=findPat(1:Nsurv,names(muMed))
		names(muSI) = sapply(strsplit(muSI,"_"),function(x){x[2]})
	} else {
		muCI = currentRes$extra$priors$p_Sfullest[,7]
		muSI = currentRes$extra$priors$surveySfull_prior[,7]
	}

	if (is.null(assYrs)) {
		assSentence = ""
	} else {
		NassY = length(assYrs)
		assSentence = paste0("The filled gold circle", ifelse(NassY>1,"s",""), " indicate", ifelse(NassY>1,"","s"), " the status in ", texThatVec(assYrs), ", which coincide", ifelse(NassY>1,"","s"), " with", ifelse(NassY>1,""," a"), " previous assessment", ifelse(NassY>1,"s",""), " for this species. ")
	}
	
	}

	## savingObjects-------------------------------
	processObj = function(x) {
		if (is.null(dim(x))) return(x)
		else {
			dimx = dim(x); nr=dimx[1]; nc=dimx[2]
			return(x[1:nr,1:nc])
		}
	}
	Bmcmc <- list()
	SA=strsplit(model.name,split="-")  # Species & Area
	Bmcmc[[toupper(SA[[1]][1])]][[SA[[1]][2]]] <- list(
		B0.MCMC        = processObj(B0.MCMC),
		Bt.MCMC        = processObj(Bcurr.MCMC),
		Bmsy.MCMC      = processObj(Bmsy.MCMC),
		Blim.MCMC      = processObj(HRPlist$blimHRP),
		Btar.MCMC      = processObj(HRPlist$btarHRP),
		Utar.MCMC      = processObj(HRPlist$utarHRP),
		P.MCMC         = processObj(currentMCMC$P),
		B.MCMC         = processObj(currentMCMC$B),
		R.MCMC         = processObj(currentMCMC$R),
		VB.MCMC        = processObj(currentMCMC$VB),
		U.MCMC         = processObj(currentMCMC$U),
		UoverUmsy.MCMC = processObj(currentMCMC$UoverUmsy),
		BoverBmsy.MCMC = processObj(currentMCMC$BoverBmsy),
		MSY.MCMC       = processObj(msy.MCMC),
		Umsy.MCMC      = processObj(umsy.MCMC),
		RPprobTabs     = processObj(RPprobTabs)
		)
#browser();return()
	save("Bmcmc",file=paste0(mcdir,"/Bmcmc-",model.name,".rda"))

	invisible(Bmcmc)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~redo.currentMCMC


## repeatMPD----------------------------2021-05-10
##  Repeat MPDs for likelihood profiles.
## ---------------------------------------------RH
repeatMPD = function(M=seq(0.07,0.08,0.01), A=c(40,45,50), R0=NULL,
	prefix="WWR-MMM-", clean=FALSE, argsMPD="", dpY)
{
	## Start subfunctions
	## Determine number of decimal places with non-trailing zeroes
	## See user 'darocsig' @ https://stackoverflow.com/questions/5173692/how-to-return-number-of-decimal-places-in-r
	decimalplaces <- function(x) {
		dc = as.character(); dp = as.numeric()
		for (xx in x) {
			if (abs(xx - round(xx)) > .Machine$double.eps^0.5) {
				dc = c(dc, strsplit(sub('0+$', '', as.character(xx)), ".", fixed = TRUE)[[1]][2] )
				dp = c(dp, nchar(strsplit(sub('0+$', '', as.character(xx)), ".", fixed = TRUE)[[1]][2]) )
			} else {
				dc = c(dc, NA)
				dp = c(dp, 0)
			}
		}
		## Automatically determine minimum number of secimal places to yield unique values
		udp = max(dp) + 1
		for (i in max(dp):1){
			#print(length(unique(round(x,i)))==length(x))
			if (length(unique(round(x,i))) == length(x))
				udp = udp - 1
		}
#browser();return()
		return(list(dc=dc, dp=dp, mindp=udp))
	}
	cleanup = function(){
		junkpat = c("^[Aa]watea","^admodel","\\.log$","\\.pst$","\\.out$","\\.rpt$","\\.tmp$","^variance$","^results.dat$","^likelihood.dat$")
		junkit  = sapply(junkpat,function(x){list.files(pattern=x)})
		junkit  = sapply(junkit,setdiff,"Awatea.exe")
		junk    = sapply(junkit,function(x){ if (length(x)>0) for (i in x) if (file.exists(i)) file.remove(i)})
	}
	## End subfunctions

	if (clean) cleanup()
	ncA   = max(nchar(A))
	if (!is.null(M)) {
		Y = Ydec = Yval = M; Ypref="MM-"; Yvar="MMM"
	} else if (!is.null(R0)) {  ## R0 expressed as log(R0)
		Y = R0; Ydec=Y/100; Yval=exp(Y); Ypref="RR-"; Yvar="RRR"
	}
	dcY = decimalplaces(Ydec)
	if (missing(dpY))
		dpY = dcY$mindp
	fileX = c(paste0(c("results","likelihood"),".dat"), paste("Awatea",c("par","std","cor","eva"),sep=".")) ## MPD files to save

	for (a in A){
		afile = paste0(prefix,"A",pad0(a,ncA),".txt")
		aline = readLines(afile)
		for (i in 1:length(Y)){
			yval = Yval[i]
			ychr = round(Ydec[i], dpY)
			ipref = paste0(sub(Ypref,paste0(sub("0\\.","",show0(ychr,dpY)),"-"),prefix),"A",pad0(a,ncA))
			.flush.cat(paste0("Processing run '", ipref, "' ..."), "\n")
			ifile = paste0(ipref,".txt")
#browser();return()

			## for some reason, certain combos of M & A do not converge unless M is nudged (WWR 2019)
			if (yval==0.03 && A==45) smudge = 0.00001
			else smudge = 0 
			iline = gsub(Yvar, yval+smudge, aline)
			writeLines(iline, con=ifile)
#browser();return()
			expr=paste("mess = shell(cmd=\"awatea -ind ",ifile,argsMPD,"\", wait=TRUE, intern=TRUE)",sep="")
			.flush.cat("   ", expr, "\n")
			eval(parse(text=expr))

			for (jfile in fileX){
				if (file.exists(jfile)){
					if (substring(jfile,1,6)=="Awatea")
						suffix  = sapply(strsplit(jfile,"\\."),tail,1)
					else if (jfile=="likelihood.dat") suffix = "lik"
					else if (jfile=="results.dat") suffix = "res"
					else suffix = "tmp"
					file.copy(jfile, paste0(ipref,".",suffix), overwrite=TRUE)
				}
			}
			if (clean) cleanup()
			rubbish = gc(verbose=FALSE)
		} ## end m loop (natural mortality)
	} ## end a loop (maximum age for plus class)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~repeatMPD


## splitGear ---------------------------2019-12-02
## Split catches by year from multiple gear types.
## Specifically to separate VB output from Awatea.
## ---------------------------------------------RH
splitGear = function(dat, fn=function(x){sum(x,na.rm=TRUE)})
{
	out  = list()
	ayrs = as.numeric(substring(colnames(dat),1,4))
	#aval = cut(ayrs,breaks=(min(ayrs)-1):max(ayrs),labels=FALSE)
	yrs  = unique(ayrs); #names(yrs) = unique(aval)
	if (length(ayrs)==length(yrs)) {
		adat = dat
		dimnames(adat) = list(rownames(dat), yrs)
		out[[1]] = adat
	} else {
		gear = as.numeric(substring(colnames(dat),6))
		for (i in gear){
			idat = dat[,is.element(gear,i)]
			colnames(idat) = yrs
			out[[i]] = idat
		}
	}
#browser();return()
	return(out)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~splitGear


#tabSAR---------------------------------2011-07-18
# Generate comma-delimited, two-dimensional output tables from reference point objects.
#  models - names of binary system files that store the decision tables.
#  pnam   - name of list object containing matrices of reference probabilities.
#  tnam   - names of matrices reporting times to reach reference points/criteria.
#  cats   - catch strategies (subset) to report in output tables.
#  digits - number of digits to retain after the decimal.
#-----------------------------------------------RH
tabSAR = function(models=paste("input-ymr",pad0(c(29,30),2),pad0(1,2),sep="."),
    #prefix="input-ymr", run=c(29,30), rwt=1,
    pnam = "refProbs3Gen90", tnam=c("Ttab0.5", "Ttab0.8", "Ttab0.95"),
    cats = seq(0,2500,500), digits=2 )
{
	#models = paste(prefix,pad0(run,2),pad0(rwt,2),sep=".")
	files  = paste(models,"Tables.RData",sep="")
	nfiles = length(files)
	for (i in 1:nfiles) {
		ifile = files[i]
		load(ifile)

		pcsv = gsub("Tables\\.RData","_prob.csv",ifile)  # output CSV name for probabilities
		cat(models[i],"\n",file=pcsv)
		cat("Annual catch,,,,Projection Year,,,\n",file=pcsv,append=TRUE)
		probs = get(pnam)
		cat(paste(c("strategy",dimnames(probs[[1]])[[2]]),collapse=","),"\n",file=pcsv,append=TRUE)
		for (j in names(probs)) {
			cat(paste("P(Bt > ",j,")",sep=""),"\n",file=pcsv,append=TRUE)
			ptab = probs[[j]]
			ptab = ptab[as.character(cats),]
			mess = paste(paste(dimnames(ptab)[[1]],apply(ptab,1,function(x){
				paste(show0(round(x,digits),digits,add2int=TRUE),collapse=",")}),sep=","),collapse="\n") # flatten table
			cat(mess,"\n",file=pcsv,append=TRUE)
		}

		tcsv = gsub("Tables\\.RData","_targ.csv",ifile)  # output CSV name for years to target
		cat(models[i],"\n",file=tcsv)
		cat("Annual catch,,,,Target Reference,,\n",file=tcsv,append=TRUE)
		for (k in tnam) {
			ttab = get(k)
			if (k==tnam[1])
				cat(paste(c("strategy",dimnames(ttab)[[2]]),collapse=","),"\n",file=tcsv,append=TRUE)
#browser();return()
			cat(paste(as.numeric(substring(k,5))*100,"% confidence",sep=""),"\n",file=tcsv,append=TRUE)
			ttab = ttab[as.character(cats),]
			mess = paste(paste(dimnames(ttab)[[1]],apply(ttab,1,function(x){
				paste(show0(round(x,digits),digits,add2int=TRUE),collapse=",")}),sep=","),collapse="\n") # flatten table
			cat(mess,"\n",file=tcsv,append=TRUE)
		}
		
	}
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^tabSAR


## tex.that.vec-------------------------2019-01-09
##  Convert a vector to a phrase 'x, y, and z'
##  Generally use 'texThatVec' in PBStools if avaialable,
##   but this package does not rely on PBStools.
##  The code below is not maintained as it is in PBStools.
## ---------------------------------------------RH
tex.that.vec = function(vec, simplify=TRUE)
{
	if (exists("texThatVec")) {
		texvec = texThatVec(vec, simplify=simplify)
	} else {
		if (length(vec)==1) return(paste0(vec))
		if (simplify) {
			if (is.character(vec))  vec = as.numeric(vec)
			uvec = sort(unique(vec))
			## User: A5C1D2H2I1M1N2O1R2T1 (140719)
			## https://stackoverflow.com/questions/24837401/find-consecutive-values-in-vector-in-r
			lvec = split(uvec,cumsum(c(1, diff(uvec) != 1)))
			cvec = sapply(lvec, function(x) {
				if (length(x)==1) return(paste0(x))
				else return(paste0(x[1],"-",rev(x)[1]))
			})
		} else cvec = vec
		texvec = paste0(paste0(cvec[1:(length(cvec)-1)], collapse=", "), ifelse(length(cvec)>2,",",""), " and ", rev(cvec)[1])
	}
	return(texvec)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~tex.that.vec
