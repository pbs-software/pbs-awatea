#--------------------------------------------------------------------#
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
# extra = bits and bobs requested by Andy                            #
#--------------------------------------------------------------------#

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
		else as.numeric(scan(file, what = "", skip = line, nlines = 1,
			quiet = TRUE))
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
		Nsurveyindex = as.numeric(rev(strsplit(resvec[grep("^Nsurveyindex",resvec)],split=sep)[[1]])[1])
		NCPUEindex   = as.numeric(rev(strsplit(resvec[grep("^NCPUEindex",resvec)],split=sep)[[1]])[1])
		elist = sapply(extra,function(x,resvec) { 
			ex = as.list(x); names(ex) = x
			for (i in ex) {
				expr=paste("index = grep(\"^",i,sep,"\",resvec)",sep="")
				eval(parse(text=expr))
				if (length(index)==1) {
					if (i %in% c("M1_prior","M2_prior","uinit_prior","p_plusscale")) index = index + (1:Nsexes) - 1
					if (i %in% c("log_qsurvey_prior","surveySfull_prior","p_surveySfulldelta",
						"log_surveyvarL_prior","log_surveyvarR_prior")) index = index + (1:Nsurveyindex) - 1
					if (i %in% c("p_Sfullest","p_Sfulldelta","log_varLest_prior","log_varRest_prior",
						"errSfull_prior","errvarL_prior","errvarR_prior",
						"log_qCPUE_prior","log_BetaCPUE_prior","qCPUEerr_prior")) index = index + (1:NCPUEindex) - 1
					exvec = strsplit(resvec[index],split=sep)
#browser();return()
					exmat = t(sapply(exvec,function(x){as.numeric(x[!is.element(x,c(i,""))])}))
					if (nrow(exmat)==1) exres=as.vector(exmat)
					else                exres=exmat
#if (i=="log_qsurvey_prior") {browser();return() }
					ex[[i]] = exres
				}
				else ex[[i]] = "not found"
			}
			return(ex) }, resvec=resvec, simplify=FALSE)

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
#test=importRes(res.file="input5-ymr.001.res", Dev=T, CPUE=T, Survey=T, CLc=T, CLs=T, CAs=T, CAc=T, extra=T)
#test=importRes("s3age-estmh00.res", CAc=T, extra=T)

