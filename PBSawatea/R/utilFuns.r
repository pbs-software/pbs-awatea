#===============================================================================
# PBSawatea utility functions:
#  importCor        : import Awatea parameter correlations.
#  importEva        : import Awatea Hessian eigenvalues.
#  importLik        : import Awatea likelihoods.
#  importPar        : import Awatea parameters (all).
#  importRes        : import Awatea results.
#  importStd        : import Awatea output parameter standard deviations.
#  makeErrMat       : mMake simple ageing error matrix for Awatea.
#  tabSAR           : generate comma-del., 2-D tables from reference point objects.
#===============================================================================

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
importEva = function(eva.file){
	efile=as.numeric(read.table(eva.file,header=FALSE,sep=""))
	out = list(eva=efile)
	return(out) }

#importLik------------------------------2012-08-08
# Import Awatea likelihoods.
#-----------------------------------------------RH
importLik = function(lik.file) {
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
importPar = function(par.file) {
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
	return(out) }
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^importPar
#test=importPar("something.par")

#importRes--------------------------------------------------2013-09-04
# Awatea res file has following structure (some elements may be      #
# missing dependent on model configuration and importCol details.    #
#                                                                    #
# N predicted numbers at age                                         #
# B predicted biomass, recruitment, and observed landings            #
# Sel predicted selectivity and observed maturity (age things)       #
# Dev predicted recruitment deviates from the stock-recruitment curve#
# CPUE, Survey commercial and survey abundance index and fit         #
# CAc, CAs commercial and survey C@A (catch at age) and fit          #
# CLc, CLs commercial and survey C@L (catch at length) and fit       #
# LA observed L@A and fit                                            #
# extra = bits and bobs requested by Andy                            #
#-------------------------------------------------------------------RH

importRes <- function (res.file, info="", Dev=FALSE, CPUE=FALSE, 
     Survey=FALSE, CAc=FALSE, CAs=FALSE, CLc=FALSE, CLs=FALSE, 
     LA=FALSE, quiet=TRUE, extra=TRUE, sep=" ") {

	#---SUBFUNCTIONS-------------------------------
	# Flush the cat down the console
	.flush.cat = function(...) { cat(...); flush.console() }

	readVector <- function(keyword, same.line = TRUE, file = res.file,
		vector = res.vector) {
		line <- match(keyword, substring(vector, 1, nchar(keyword)))
		v <- if (same.line)
			as.numeric(scan(file, what = "", skip = line - 1, nlines = 1, quiet = TRUE)[-1])
		else as.numeric(scan(file, what = "", skip = line, nlines = 1, quiet = TRUE))
		if (!quiet) .flush.cat("vector...")
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
		if (!quiet) .flush.cat("matrix...")
		return(m)
	}
	getN <- function(sexes, years, ages) {
		if (!quiet) .flush.cat("N         ")
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
		if (!quiet) .flush.cat("OK\n")
		return(N)
	}
	getB <- function(years, gears) {
		ngears <- length(gears)
		if (!quiet) .flush.cat("B         ")
		vb <- readMatrix("Vulnerable_Biomass_by_Method_and_Year", nrow = ngears)
		sb <- readVector("Spawning_Biomass_by_Year", same.line = FALSE)
		# *** ADD in the exploitation rate, last year is missing and need
		#     to pad the matrix with missing values to match other "B" columns.
		U <- readMatrix("Exploitation_Rate_by_Method_and_Year", nrow=ngears )
		U <- cbind( U, rep(NA,nrow(U)) )
#		y <- c(readVector("Total_Catch_by_Method_and_Year", same.line = FALSE), NA)
		# BUG FIX: Appears should call readMatrix to accommodate multiple gear series.
		#          Then, sum over the gears to get total catch.
		y <- readMatrix( "Total_Catch_by_Method_and_Year", nrow=ngears )
		y <- apply( y,2,sum,na.rm=T )
		y <- c( y,NA )

		B <- as.data.frame( cbind(years, t(vb), sb, y, t(U)) )
		names(B) <- if (ngears == 1) c("Year", "VB", "SB", "Y", "U")
		else c("Year", paste("VB", gears, sep = "."), "SB", "Y", paste("U",gears,sep=".") )
		if (!quiet) .flush.cat("OK\n")
		return(B)
	}
	getSel <- function(gears, surveys, years, sexes, ages) {
		if (!quiet)
			.flush.cat("Sel       ")
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
		if (!quiet) .flush.cat("OK\n")
		return(Sel)
	}
	getDev <- function(ages, years) {
		if (!quiet) .flush.cat("Dev       ")
		Dev <- list()
		Dev$Initial <- readVector("log_InitialDev", same.line = TRUE)
		names(Dev$Initial) <- ages[-c(1, length(ages))]
		Dev$Annual <- readVector("log_RecDev", same.line = TRUE)
		names(Dev$Annual) <- years[-length(years)]
		if (!quiet) .flush.cat("OK\n")
		return(Dev)
	}
	getCPUE <- function(gears, years) {
		if (!quiet) .flush.cat("CPUE      ")
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
		if (!quiet) .flush.cat("OK\n")
		return(CPUE)
	}
	getSurvey <- function(years) {
		if (!quiet) .flush.cat("Survey    ")
		nyears <- length(years)
		nseries <- readVector("Nsurveyindex")
		obs <- readMatrix("indexyearvaluecv", nrow = readVector("Number_of_survey_data", same.line = FALSE))
		obs <- data.frame(Series = obs[, 1], Year = obs[, 2], Obs = obs[, 3], CV = obs[, 4])
		fit <- readMatrix("Survey_Index_Trajectories", nrow = nseries)
		fit <- data.frame(Series = rep(1:nseries, each = nyears), Year = rep(years, nseries), Fit = as.vector(t(fit)))
		Survey <- merge(obs, fit, all = TRUE)
		Survey$Series <- as.integer(Survey$Series)
		Survey$Year <- as.integer(Survey$Year)
		if (!quiet) .flush.cat("OK\n")
		return(Survey)
	}
	getCAc <- function(sexes, ages) {
		if (!quiet) .flush.cat("CAc       ")
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
		if (!quiet) .flush.cat("OK\n")
		return(CAc)
	}
	getCAs <- function(sexes, ages) {
		if (!quiet) .flush.cat("CAs       ")
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
		if (!quiet) .flush.cat("OK\n")
		return(CAs)
	}
	getCLc <- function(sexes, lengths) {
		if (!quiet) .flush.cat("CLc       ")
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
		if (!quiet) .flush.cat("OK\n")
		return(CLc)
	}
	getCLs <- function(sexes, lengths) {
		if (!quiet) .flush.cat("CLs       ")
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
		if (!quiet) .flush.cat("OK\n")
		return(CLs)
	}
	getLA <- function(sexes, ages) {
		if (!quiet) .flush.cat("LA        ")
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
		LA <- merge(obs, fit, by = c("Sex", "Age"), all = T)
		LA$Age <- as.integer(LA$Age)
		LA$Fit <- LA$Fit
		LA$CV <- LA$CV
		if (!quiet) .flush.cat("OK\n")
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
		glist = list(general=list(Nsexes=Nsexes,Nages=Nages,Nmethods=Nmethods,NCPUEindex=NCPUEindex,Nsurveyindex=Nsurveyindex))
		elist = sapply(extra,function(x,resvec) { 
			ex = as.list(x); names(ex) = x
#if (x=="Virgin_Vulnerable_Biomass") {browser();return() }
			for (i in ex) {
				expr=paste("index = grep(\"^",i,sep,"\",resvec)",sep="")
				eval(parse(text=expr))
				if (length(index)==1) {
					if (i %in% c("M1_prior","M2_prior","Rinit_prior","uinit_prior","p_plusscale")) index = index + (1:Nsexes) - 1
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
#browser();return()
	return(c(glist,elist,clist))
	}
	#---END SUBFUNCTIONS---------------------------

	if (!file.exists(res.file))
		stop("File ", res.file, " not found. Use / or \\\\ separators.")
	file.copy(from=res.file, to=paste(res.file,".temp",sep=""), overwrite=TRUE)            # save temporary version
	file.rename(from=res.file, to=paste(res.file,".original",sep=""))                      # rename original file
	file.copy(from=paste(res.file,".temp",sep=""), to=res.file, overwrite=TRUE)            # save temporary version
	exitfun = function(resfile) {
		file.rename(from=paste(resfile,".original",sep=""), to=resfile)
		file.remove(paste(resfile,".temp",sep="")) }
	on.exit( exitfun(res.file) )                                                           # restore original file and remove the temporary copy
	res.vector <- readLines(res.file)
	res.vector <- gsub("[?]","",res.vector,perl=T)
	writeLines(res.vector,con=res.file)

	res.vector <- gsub("\t", " ", gsub(sep, " ", res.vector))
#browser();return()
	if (extra)
		extra = getExtra(res.vector,sep=sep)
	else extra = NULL

	res.vector <- gsub("\"", "", gsub("\t", "", gsub(" ", "", res.vector)))
	if (!quiet)
		.flush.cat("\nParsing text file ", res.file, ":\n\nPreamble  ", sep = "")
	sexes   <- if (readVector("Nsexes") == 1) "Unisex"
		else c("Female", "Male")
	gears   <- seq(1, length.out = readVector("Nmethods"))
	surveys <- seq(1, length.out = readVector("Nsurveyindex"))
	years   <- seq(from = readVector("StartYear"), to = readVector("EndYear") + 1)
	ages    <- seq(from = 1, to = readVector("Nages"))
	lengths <- seq(from = readVector("First_length"), by = readVector("Length_class_increment"),
		length.out = readVector("Number_of_length_classes"))
	if (!quiet) .flush.cat("OK\n")
	model   <- list()
	model$N <- getN(sexes, years, ages)
	model$B <- getB(years, gears)
	rec     <- model$N[model$N$Age == 1, ]
	rec     <- tapply(rec$N, rec$Year, sum)
	model$B$R <- c(rec[-1], NA)
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
#browser();return()
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
	if (!quiet) .flush.cat("\n")
	attr(model, "call") <- match.call()
	attr(model, "scape.version") <- installed.packages()["scape", "Version"]
	attr(model, "info") <- info
	class(model) <- "scape"
	return(model)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~importRes


#importStd------------------------------2012-07-27
# Import Awatea table of estimated parameters.
#-----------------------------------------------RH
importStd = function(std.file, vnam="name") {
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
	return(out) }
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^importStd
#test=importStd("something.std")

#makeErrMat-----------------------------2011-05-05
# Make a simple ageing error matrix for Awatea.
#-----------------------------------------------RH
makeErrMat = function(N=60, ondiag=0.8, offdiag=0.1, corner=0.9) {
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
    cats = seq(0,2500,500), digits=2 ) {

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


#MAfun----------------------------------2013-09-04
# Mean age function (Chris Francis, 2011, weighting assumption T3.4, p.1137)
# MAfun is used in both `runADMB` and in `runSweave`.
#-----------------------------------------------RH
MAfun = function(padata,brks=NULL)
{
	padata=padata[padata$Age>=padata$startL & padata$Age<=padata$endL,]# Choose only ages older than `startL`
	# S = series, y = year, a = age bin, O = observed proportions, P = Predicted (fitted) proportions, N=sample size
	S=padata$Series; y=padata$Year; a=padata$Age; O=padata$Obs; E=padata$Fit; SS=padata$SS   # note: SD and NR not used
	if (is.null(brks)) {
		f = paste(S,y,sep="-"); J = unique(S) }
	else {
		B = cut(y, breaks=brks, include.lowest=TRUE, labels=FALSE)
		f = paste(S,B,y,sep="-"); J = unique(paste(S,B,sep="-")) }
	# make sure that input age proportions are standardised (especially if females only)
	O    = as.vector(sapply(split(O,f),function(x){x/sum(x)})) #standardise O (obs props)
	#E    = as.vector(sapply(split(E,f),function(x){x/sum(x)})) #standardise E (fitted props)  -- causes instability
	Oay  = a * O; Eay = a * E; Eay2 = a^2 * E
	mOy  = sapply(split(Oay,f),sum,na.rm=TRUE)   # mean observed age
	mEy  = sapply(split(Eay,f),sum,na.rm=TRUE)   # expected mean age
	mEy2 = sapply(split(Eay2,f),sum,na.rm=TRUE)
	Vexp = mEy2-mEy^2                            # variance of expected mean age
	N    = sapply(split(SS,f),mean,na.rm=TRUE)
	Yr   = as.numeric(substring(names(mOy),nchar(names(mOy))-3))
	return(list(MAobs=mOy, MAexp=mEy, Vexp=Vexp, N=N, J=J, Yr=Yr)) # observed and expected mean ages, variance of expected ages
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MAfun

