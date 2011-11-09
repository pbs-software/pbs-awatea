setClass ("AWATEAdata", 
     representation( txtnam="character", input="character", vlst="list", 
     dnam="character", nvars="numeric", vdesc="character", vars="list", 
     gcomm="character", vcomm="character", output="list" , reweight="list") )

#require(PBSmodelling)
#source("importRes.r",local=FALSE)

# Flush the cat down the console
.flush.cat = function(...) { cat(...); flush.console() }

#runADMB--------------------------------2011-05-31
# Run AD Model Builder code for Awatea
#-----------------------------------------------RH
runADMB = function(filename.ext, wd=getwd(), strSpp="YMR", runNo=25, rwtNo=0,
     doMPD=FALSE, N.reweight=0, cvpro=FALSE, mean.age=TRUE, 
     doMCMC=FALSE, mcmc=1e6, mcsave=1e3, ADargs=NULL, verbose=FALSE, 
     doMSY=FALSE, msyMaxIter=15000., msyTolConv=0.01, endStrat=0.301, stepStrat=0.001, ...) {

	ciao = function(){setwd(cwd); gc(verbose=FALSE)} # exit function
	cwd = getwd(); on.exit(ciao())
	runNoStr = pad0(runNo,2)
	runname  = paste(strSpp,"run",runNoStr,sep="")
	rundir   = paste(wd,runname,sep="/")
	if (file.exists("results.dat")) file.remove("results.dat")
	ext = sapply(strsplit(filename.ext,"\\."),tail,1)
	prefix = substring(filename.ext,1,nchar(filename.ext)-nchar(ext)-1)
	if (substring(prefix,6,7) != runNoStr) showAlert("input file appears mismatched with Run")
	prefix = gsub(runNoStr,"",prefix) # get rid of superfluous run number in name
	argsMPD = argsMCMC = NULL
	if (!is.null(ADargs)) {
		if(!is.list(ADargs)) stop("'ADargs' must be either a list or NULL")
		for (i in ADargs) {
			if (!grepl("mc",i))     argsMPD = c(argsMPD,list(i)) 
			if (!grepl("nohess",i) & !grepl("mcmc",i) & !grepl("mcsave",i)) argsMCMC = c(argsMCMC,list(i)) 
		}
	}
	argsMPD  = paste(paste(" ",unlist(argsMPD),sep=""),collapse="")
	argsMCMC = paste(paste(" ",unlist(argsMCMC),sep=""),collapse="")

	if (doMPD) {
		.flush.cat("Reweighting surveys and proportions-at-age...\n")
		if (!file.exists(filename.ext)) stop("Specified input file does not exist")
		sdnrfile = paste(prefix,runNoStr,"sdnr",sep=".")
		cat("#SDNR (Surveys, CPUE, CAc, CAs)\n",file=sdnrfile)
		#file0 = gsub(paste("\\.",ext,sep=""),paste(".000.",ext,sep=""),filename.ext)
		file0 = paste(prefix,runNoStr,"00.txt",sep=".")
		for (i in 0:N.reweight) {
			ii = pad0(i,2)
			if (i==0) {
				fileN = file0
				file.copy(from=filename.ext,to=file0,overwrite=TRUE) }
			else {
				desc = Robj@vdesc; dat = Robj@vars; newdat=Robj@reweight
				rvar1 = names(findPat("Survey abundance indices",desc))
				rvar2 = names(findPat("CPUE \\(Index",desc))
				rvar3 = names(findPat("Commercial catch at age data",desc))
				rvar4 = names(findPat("Survey C@A data",desc))

				survey = dat[[rvar1]]; #dimnames(index) = list(1:nrow(index),c("series","year","obs","CV"))
				survey[,4] = newdat$survey$CVnew
				cpue = dat[[rvar2]]
				if (!is.null(newdat$cpue))
					cpue[,5] = newdat$cpue$CVnew
				CAc = dat[[rvar3]]
				CAc[,3] = newdat$wNcpa
				CAs = dat[[rvar4]]
				CAs[,3] = newdat$wNspa

				Robjnew = fix(Robj,c(rvar1,rvar2,rvar3,rvar4),list(survey,cpue,CAc,CAs))
				#hh=pad0(i-1,3); fileN = gsub(hh,ii,fileN)
				fileN = paste(prefix,runNoStr,ii,"txt",sep=".")
				write(Robjnew,fileN)
			}
			expr=paste("mess = shell(cmd=\"awatea -ind ",fileN,argsMPD,"\", wait=TRUE, intern=TRUE)",sep=""); eval(parse(text=expr))
			if (verbose)  .flush.cat(mess, sep="\n")
			if (length(mess)<10) stop("Abnormal program termination")
			file.copy("results.dat",gsub(paste("\\.",ext,sep=""),".res",fileN),overwrite=TRUE)
			file.copy("Awatea.par",gsub(paste("\\.",ext,sep=""),".par",fileN),overwrite=TRUE)
			eval(parse(text=paste("Robj = readAD(\"",fileN,"\")",sep="")))
			Robj@reweight = list(nrwt=i)
			Robj = reweight(Robj, cvpro=cvpro, mean.age=mean.age, sfile=sdnrfile, fileN=fileN)
		}
		if (!file.exists(rundir)) dir.create(rundir)
#browser();return()
		filesN = paste(prefix,runNoStr,rep(pad0(0:N.reweight,2),each=3),rep(c(ext,"par","res"),N.reweight+1),sep=".")
		file.copy(paste(wd,c(filename.ext,filesN,sdnrfile),sep="/"),rundir,overwrite=TRUE)
		file.remove(paste(wd,c(filesN,sdnrfile),sep="/"))
	}
	if (doMCMC | doMSY) {
		rwtNoStr = pad0(rwtNo,2)
		fileN    = paste(prefix,runNoStr,rwtNoStr,ext,sep=".")
		mcname   = paste("MCMC",runNoStr,rwtNoStr,sep=".")
	}
	if (doMCMC) {
		.flush.cat(paste("Running",mcmc,"MCMC iterations...\n"))
		mcdir    = paste(wd,runname,mcname,sep="/")
		if (!file.exists(mcdir)) dir.create(mcdir)
		filesN = paste(prefix,runNoStr,rwtNoStr,c(ext,"res","par"),sep=".")
		file.copy(paste(rundir,filesN,sep="/"),mcdir,overwrite=TRUE); setwd(mcdir)
#browser();return()
		#if (!doMPD) {
		#	eval(parse(text=paste("Robj = readAD(\"",fileN,"\")",sep="")))
		#	Robj = reweight(Robj, cvpro=cvpro, mean.age=mean.age) 
		#}
		#expr=paste("mess = shell(cmd=\"awatea -ind ",fileN," -mcmc ",format(mcmc,scientific=FALSE)," -mcsave ",
		#	format(mcsave,scientific=FALSE),argsMCMC,"\", wait=TRUE, intern=TRUE)",sep="")
		#.flush.cat(expr, sep="\n")
		#eval(parse(text=expr))
		#if (verbose)  .flush.cat(mess, sep="\n")
		expr=paste("shell(cmd=\"awatea -ind ",fileN," -mceval\" , wait=TRUE, intern=FALSE)",sep="")
		.flush.cat(expr, sep="\n")
		eval(parse(text=expr))
		
		Robj="dummy4now"
	}
	if (doMSY) {
		.flush.cat(paste("Running MSY yield calculations...\n"))
		msyname  = paste("MSY",runNoStr,rwtNoStr,sep=".")
		msydir   = paste(wd,runname,mcname,msyname,sep="/")
		if (!file.exists(msydir)) dir.create(msydir)
		filesN = c(paste(wd,runname,fileN,sep="/"),paste(wd,runname,mcname,"Awatea.psv",sep="/"))
		file.copy(filesN,msydir,overwrite=TRUE); setwd(msydir)
		ctlfile  = paste(c("#MSY control file",
			"#Maximum number of iterations",format(msyMaxIter,scientific=FALSE),
			"#Tolerance for convergence",format(msyTolConv,scientific=FALSE) ), collapse="\n")
		cat(ctlfile,"\n",sep="",file="Yields.ctl")

		infile = readAD(fileN)
		strategy = view(infile,pat=c("Strategy","End year"),see=FALSE)
		# Reset the strategy and express in terms of U
		strategy[[grep("Strategy Type",names(strategy))]] = 2
		strategy[[grep("End year of projections",names(strategy))]]  = -99
		strategy[[grep("End Strategy",names(strategy))]]  = endStrat
		strategy[[grep("Step Strategy",names(strategy))]] = stepStrat
		vnam = substring(names(strategy),1,4)
		infile = fix(infile,vnam,strategy) # replace the contents of infile with updated strategy
		write(infile,fileN)                # overwrite the input file for MSY calculations
		expr=paste("shell(cmd=\"awatea -ind ",fileN," -mceval\" , wait=TRUE, intern=FALSE)",sep="")
		.flush.cat(expr, sep="\n")
		eval(parse(text=expr))
		
		Robj=list(ctlfile,strategy)
			
	}
	return(Robj) }
#------------------------------------------runADMB

#readAD---------------------------------2011-05-31
# Read the ADMB input file and create an AWATEA class object.
#-----------------------------------------------RH
readAD = function(txt) {
	#require(PBSmodelling)
	txtnam = as.character(substitute(txt))
	otxt = readLines(txt) # original text
	xlst = strsplit(otxt,"")
	xlst = xlst[sapply(xlst,length)>0]
	ntxt = sapply(xlst,function(x){paste(clipVector(x,clip="\t",end=2),collapse="")})
	ntxt = gsub("\\\t\\\t","\t",ntxt)   # internal cleanup
	vlst = list(); vcom=NULL; acom=NULL
	for (i in 1:length(ntxt)) {
		if (substring(ntxt[i],1,1)=="#") {
			expr = paste("vlst[[\"L",pad0(i,3),": Comment\"]] = ",deparse(ntxt[i],width=500),sep="")
			comm = substring(ntxt[i],2,) ; names(comm)=i
			acom = c(acom,comm) }
		else {
			vcom = c(vcom,comm)
			expr = paste("vlst[[\"L",pad0(i,3),": Data {",comm,"}\"]] = c(",gsub("\\\t",",",ntxt[i]),")",sep="") }
		eval(parse(text=expr))
	}
	
	# description of variables with inputs
	vdesc = unique(vcom) 
	gcomm = acom[!is.element(acom,vdesc)] # general comments
	vcomm = acom[is.element(acom,vdesc)]  # variable comments (may be duplicated)
	nvars = length(vdesc)
	names(vdesc) = paste("v",pad0(1:nvars,3),sep="")
	vars = as.list(vdesc)

	dnam = names(vlst); names(dnam)=1:length(dnam)
	dnam = dnam[grep("Data",dnam)]
	dnam = substring(dnam,13,nchar(dnam)-1)
	
	for (i in names(vdesc)) {
		ivar = vlst[as.numeric(names(dnam[is.element(dnam,vdesc[i])]))]
		ilen = length(ivar)
		if (ilen==1) vars[[i]] = ivar[[1]]
		else if (ilen>=2) {
			jvar=NULL
			for (j in 1:length(ivar))
				jvar=rbind(jvar,ivar[[j]])
			vars[[i]] = jvar
		}
		else {browser();return()}
	}
	#writeList(vars,gsub("\\.txt",".pbs.txt",txtnam),format="P") # write to a PBS formatted text file
	resnam = gsub("\\.txt$",".res",txtnam)
	if (!file.exists(resnam))
		output = list()
	else {
		output = importRes( res.file=resnam, Dev=TRUE, CPUE=TRUE, Survey=TRUE, CLc=FALSE, CLs=FALSE, CAs=TRUE, CAc=TRUE, extra=TRUE)
		attr(output,"class")="list" }
	Data=new("AWATEAdata",txtnam=txtnam, input=ntxt, vlst=vlst, dnam=dnam, 
		nvars=nvars, vdesc=vdesc, vars=vars, gcomm=gcomm, vcomm=vcomm, output=output,reweight=list(nrwt=0))
	return(Data) }
#-------------------------------------------readAD

#setMethod.view-------------------------2011-11-09
# Set the method for 'view' when using an AWATEA class.
#-----------------------------------------------RH
setMethod("view", signature(obj = "AWATEAdata"),
     function (obj, n = 1:5, last = FALSE, random = FALSE, print.console=TRUE, see=TRUE, ...) {
	dat = obj@vars; desc= obj@vdesc; nvars=obj@nvars
	dots = list(...)
	if (!is.null(dots$pat)) {
		patty = findPat(dots$pat,desc)
		if (length(patty)==0) stop("Specified pattern not found")
		else {
			seedat = dat[names(patty)]; names(seedat) = paste(names(patty),patty,sep=": ")
		}
	}
	else if (!last & !random) {
		seedat = dat[n]; names(seedat) = paste(names(desc[n]),desc[n],sep=": ") }
	else if (last) {
		seedat = dat[rev(nvars-n+1)]; names(seedat) = paste(names(desc[rev(nvars-n+1)]),desc[rev(nvars-n+1)],sep=": ") }
	else seedat="Random option not implemented for AWATEAdata class"
	if (see) print(seedat)
	else invisible(seedat)
	} )
#-----------------------------------setMethod.view

#setMethod.fix--------------------------2011-05-31
# Set the method for 'fix' when using an AWATEA class.
#-----------------------------------------------RH
setMethod("fix", signature(x = "AWATEAdata"),
    function (x, varN, xnew, ...) {
	dat = x@vars; desc=x@vdesc
	dots = list(...)
	nvar = length(varN)
	if (nvar>1 & !is.list(xnew) && length(xnew)!=nvar)
		stop("'xnew' must be a list of objects to match multiple 'varN'")
	for ( i in 1:nvar) {
		if (nvar==1 & !is.list(xnew)) oo = xnew
		else oo = xnew[[i]]
		dat[[varN[i]]] = oo
	}
	x@vars = dat
	return(x) } )
#------------------------------------setMethod.fix

#setMethod.write------------------------2011-05-31
# Set the method for 'write' when using an AWATEA class.
#-----------------------------------------------RH
setMethod("write", signature(x = "AWATEAdata"),
    function(x, file="something.txt", ncolumns = if(is.character(x)) 1 else 5,
    append = FALSE, sep = "\t") {
	vlst=x@vlst; dnam=x@dnam; vdesc=x@vdesc; vars=x@vars; gcomm=x@gcomm; vcomm=x@vcomm
	olst = list()
	for (i in 1:length(vlst)) {
		ii = as.character(i); iii=as.character(i+1)
		if (is.element(ii,names(gcomm)))
			olst[[ii]] = paste("#",gcomm[ii],sep="")
		else if (is.element(ii,names(vcomm))) {
			olst[[ii]] = paste("#",vcomm[ii],sep="")
#if(i==50) {browser();return()}
			if (is.na(dnam[iii])) next  # if variable header is duplicated and one is commented out
			olst[[iii]] = vars[[names(vdesc[is.element(vdesc,dnam[iii])])]]
		}
		else next
	}
	nlin = sum(sapply(olst,function(x){if(is.matrix(x)) nrow(x) else 1}))
	nvec = rep("",nlin)
	ii   = 0
	for (i in 1:length(olst)) {
		oo = olst[[i]]
		ii = ii + 1
		if (is.matrix(oo)) {
			for (j in 1:nrow(oo)) {
				ii = ii + ifelse(j==1,0,1)
				nvec[ii] = paste(oo[j,],collapse=sep)
			}
		}
		else if (is.numeric(oo)) 
			nvec[ii] = paste(oo,collapse=sep)
		else
			nvec[ii] = oo
	}
	writeLines(nvec,con=file)
	invisible(nvec) } )
#----------------------------------setMethod.write

#setMethod.reweight---------------------2011-05-31
# Set the method for 'reweight' when using an AWATEA class.
#-----------------------------------------------RH
reweight <- function(obj, cvpro=FALSE, mean.age=TRUE, ...) return(obj)
setMethod("reweight", signature(obj = "AWATEAdata"),
    function (obj, cvpro=FALSE,  mean.age=TRUE, ...) {
	nrwt = obj@reweight$nrwt; dat = obj@vars; desc=obj@vdesc; res=obj@output
	dots = list(...)
#	NRfun  = function(O,P,dO)  { (O-P)/dO }        # Normal residual for an observation (F.26)
	NRfun  = function(O,P,CV,maxR=0,logN=TRUE)  {  # Normal residual for an observation (Coleraine Excel)
		if (logN) {
			num = log(O) - log(P) + 0.5*log(1+CV^2)^2
			den = sqrt(log(1+CV^2))
			NR  = num/den }
		else
			NR = (O-P)/CV
		if (round(maxR,6)!=0) {
			maxR = c(-abs(maxR),abs(maxR))
			NR[NR < maxR[1]] = maxR[1]
			NR[NR > maxR[2]] = maxR[2] }
		return(NR)
	}
	SDfun = function(p,N,A) {                      # Standard deviaton of observed proportion-at-ages (F.27)
		# p = observed proportions, N = sample size, A = no. sexes times number of age bins
		SD = sqrt((p * (1-p) + 0.1/A) / N) 
		SD[!is.finite(SD)] = NA
		return(SD)
	}
	eNfun = function(S,Y,O,P,Nmax=200) {           #Re-weighted (effective) N for porportions-at-age (Coleraine Excel)
		# S = Series, Y = Year, O = Observed proportions, P = Predicted (fitted) proportions, Nmax = maximum sample size allowed
		f   = paste(S,Y,sep="-")
		num = sapply(split(P,f),function(x) { sum( x * (1 - x) ) } )
		den = sapply(split(O-P,f),function(x){ sum( x^2 ) } )
		N  = pmin(num/den, Nmax)
		return( N )
	}
	MRfun = function(y,a,O,P)  {                   # Mean age residuals (Chris Francis) (not used)
		Oay = a * O; Pay = a * P
		mOy = sapply(split(Oay,y),sum,na.rm=TRUE)
		mPy = sapply(split(Pay,y),sum,na.rm=TRUE)
		ry  = mOy - mPy
		return(ry) 
	}
	CFfun = function(y,N,a,O,P) {                  # Correction factor based on mean age residuals (not used)
		Oay = a * O; OAy = a^2 * O
		mOy = sapply(split(Oay,y),sum,na.rm=TRUE)
		MOy = sapply(split(OAy,y),sum,na.rm=TRUE)
		Vy  = MOy - mOy^2  # variance of the age frequency in year y
		ry  = MRfun(y,a,O,P)
		CF  = 1 / var(ry*(N/Vy)^.5)
		return(CF) 
	}
	MAfun = function(padata,brks=NULL)  {                   # Mean age function (Chris Francis, 2011, submitted to CJFAS))
		# S = series, y = year, a = age bin, O = observed proportions, P = Predicted (fitted) proportions, N=sample size
		S=padata$Series; y=padata$Year; a=padata$Age; O=padata$Obs; E=padata$Fit; N=padata$SS
		if (is.null(brks)) {
			f = paste(S,y,sep="-"); J = unique(S) }
		else {
			B = cut(y, breaks=brks, include.lowest=TRUE, labels=FALSE)
			f = paste(S,B,y,sep="-"); J = unique(paste(S,B,sep="-")) }
		Oay  = a * O; Eay = a * E; Eay2 = a^2 * E
		mOy  = sapply(split(Oay,f),sum,na.rm=TRUE)
		mEy  = sapply(split(Eay,f),sum,na.rm=TRUE)
		mEy2 = sapply(split(Eay2,f),sum,na.rm=TRUE)
		N    = sapply(split(N,f),mean,na.rm=TRUE)
		return(list(MAobs=mOy, MAexp=mEy, Vexp=mEy2-mEy^2, N=N, J=J)) # observed and expected mean ages, variance of expected ages
	}
	wfun = function (MAlist) {                 # weighting function TA1.8 (Francis 2011, CJFAS)
		# y=year, a=age bin, O=observed proportions, P=Predicted (fitted) proportions, N=sample size, J=series
		unpackList(MAlist,scope="L")
		w    = rep(NA,length(J)); names(w)=J
		jvec = substring(names(N),1,nchar(names(N))-5)
		wN   = N
		for (j in J) {
			z = is.element(jvec,j)
			w[j] = 1 / var((MAobs[z]-MAexp[z])/((Vexp[z]/N[z])^0.5) )  
			wN[z] = N[z] * w[j]
			}
		return(list(w=w,wN=wN))
	}

	#---Start reweight-------------------
	survey =  res$Survey[!is.na(res$Survey[,"Obs"]),]
	Sseries = sort(unique(survey[,"Series"]))
	SDNR   = rep(0,length(Sseries)); names(SDNR) = Sseries
	survey$NR = NRfun(survey$Obs,survey$Fit,survey$CV)
	survey$CVnew = rep(0,nrow(survey))
	for (s in Sseries) {
		ss = as.character(s)
		zs = is.element(survey$Series,s)
		sser = survey[zs,]
		SDNR[ss] = sd(sser$NR)
		if (cvpro && nrwt==0)
			survey$CVnew[zs] = sqrt(survey$CV[zs]^2 + cvpro^2)
		else if (cvpro && nrwt>0)
			survey$CVnew[zs] = survey$CV[zs]
		else
			survey$CVnew[zs] = sser$CV * SDNR[ss]
	}

	if (view(obj,pat="CPUE likelihood",see=FALSE)==0) 
		cpue = NULL
	else {
		cpue =  res$CPUE[!is.na(res$CPUE[,"Obs"]),]
		Utmp = sort(unique(cpue[,"Series"])); Useries = gsub("Series","cpue",Utmp); names(Useries) = Utmp
		sdnr = rep(0,length(Useries)); names(sdnr) = Useries; SDNR = c(SDNR,sdnr)
		cpue$NR = NRfun(cpue$Obs,cpue$Fit,cpue$CV)
		cpue$CVnew = rep(0,nrow(cpue))
		for (u in names(Useries)) {
			uu=as.character(u)
			zu = is.element(cpue$Series,u)
			user = cpue[zu,]
			SDNR[Useries[uu]] = sd(user[,"NR"])
			if (cvpro && nrwt==0)
				cpue$CVnew[zu] = sqrt(user$CV^2 + cvpro^2)
			else if (cvpro && nrwt>0)
				cpue$CVnew[zu] = user$CV
			else
				cpue$CVnew[zu] = user$CV * SDNR[Useries[uu]]
		}
	}

	sdnr = rep(0,2); names(sdnr) = c("cpa","spa"); SDNR = c(SDNR,sdnr)
	wj   = NULL
	CAc  = res$CAc  # commercial proportions at age
	cpa  = CAc[!is.na(CAc$SS),]
	cpa$SD = SDfun(cpa$Obs, cpa$SS, A=length(unique(cpa$Sex))*length(unique(cpa$Age)) )
	cpa$NR = NRfun(cpa$Obs, cpa$Fit, cpa$SD, maxR=3, logN=FALSE)
	if (mean.age) {
		#MAc   = MAfun(cpa,brks=c(1979,1985,1997,2004,2009)) # commercial mean ages
		MAc   = MAfun(cpa) # commercial mean ages
		Wc    = wfun(MAc)
		wNcpa = Wc$wN
		#SDNR["cpa"] = NA  # No formulae appropriate for composition-data likelihoods due to correlations (Francis 2011, Appendix B, CJFAS)
		wtemp = Wc$w; names(wtemp)=paste("cpa-",names(wtemp),sep="")
		wj = c(wj,wtemp)
	}
	else
		wNcpa  = eNfun(cpa$Series,cpa$Year,cpa$Obs,cpa$Fit)
	SDNR["cpa"] = sd(cpa$NR,na.rm=TRUE)

	CAs = res$CAs  # survey proportions at age
	spa = CAs[!is.na(CAs$SS),]
	spa$SD = SDfun(spa$Obs, spa$SS, A=length(unique(spa$Sex))*length(unique(spa$Age)) )
	spa$NR = NRfun(spa$Obs, spa$Fit, spa$SD, maxR=3, logN=FALSE)
	if (mean.age) {
		MAs  = MAfun(spa) # survey mean ages
		Ws    = wfun(MAs)
		wNspa = Ws$wN
		#SDNR["spa"] = NA  # No formulae appropriate for composition-data likelihoods due to correlations (Francis 2011, Appendix B, CJFAS)
		wtemp = Ws$w; names(wtemp)=paste("spa-",names(wtemp),sep="")
		wj = c(wj,wtemp)
	}
	else
		wNspa  = eNfun(spa$Series,spa$Year,spa$Obs,spa$Fit)
	SDNR["spa"] = sd(spa$NR,na.rm=TRUE)

#browser();return()

	if (!is.null(dots$sfile)) {
		sfile=dots$sfile
		.flush.cat(paste(dots$fileN," (",paste(names(SDNR),collapse=", "),")\n",sep=""))
		.flush.cat(paste(round(SDNR,5),collapse="\t"),"\n\n",sep="")  # cat to console
		cat("\n",dots$fileN,"\n",file=sfile,append=TRUE,sep="")
		cat("SDNR: ",paste(round(SDNR,5),collapse="\t"),"\n",file=sfile,append=TRUE,sep="")
		if (!is.null(wj)) {
			.flush.cat(paste("wj (",paste(names(wj),collapse=", "),")\n",sep=""))
			.flush.cat(paste(round(wj,5),collapse="\t"),"\n\n",sep="")  # cat to console
			cat("wj:   ",paste(round(wj,5),collapse="\t"),"\n",file=sfile,append=TRUE,sep="")
		}
	} 
	#sMAR = MRfun(spa$Year,spa$Age,spa$Obs,spa$Fit) # srvey mean age residuals
	#w = c((1/SDNR)[as.character(c(Sseries,Useries))],(1/SDNR^2)[c("cpa","spa")]) # lognormal and multinomial (Francis)
	#f = c( cpa=CFfun(cpa$Year,cpa$SS,cpa$Age,cpa$Obs,cpa$Fit), spa=CFfun(spa$Year,spa$SS,spa$Age,spa$Obs,spa$Fit) )
	obj@reweight = list(nrwt=nrwt+1, survey=survey,cpue=cpue,wNcpa=wNcpa,wNspa=wNspa,SDNR=SDNR,wj=wj)
	return(obj)
	} )
#-------------------------------setMethod.reweight

#=================================================
#out = readAD("input.txt")
#nvec = write(out)
#popin = readAD("s3age-estmh00.txt")
#popin = reweight(popin)

#=== POP ===
#out=runADMB("s3age-estmh02.txt",doMPD=F,doMCMC=T,mcmc=1000,mcsave=100)
#out=runADMB("s3age-estmh.txt",doMPD=T,N.reweight=2,ADargs=list("-nohess"))
#out=runADMB("s3age-estmh.002.txt",doMPD=F,doMCMC=T,mcmc=1000,mcsave=100)
#out=runADMB("s3age-estmh.txt",doMPD=T,N.reweight=1,doMCMC=T,mcmc=1000,mcsave=100,ADargs=list("-nohess"))

#=== YMR ===
#out=runADMB("input36-ymr.txt",runNo=36,doMPD=TRUE,N.reweight=6,ADargs=list("-nohess"),mean.age=TRUE,cvpro=0.2)
#out=runADMB("input25-ymr.txt",runNo=25,rwtNo=1,doMCMC=TRUE,mcmc=1e6,mcsave=1e3,mean.age=TRUE,cvpro=0.2)
#out=runADMB("input24-ymr.txt",runNo=24,rwtNo=1,doMSY=TRUE,msyMaxIter=15000,msyTolConv=0.01,endStrat=0.301,stepStrat=0.001)

#for (i in 29:30)
#	out=runADMB(paste("input",pad0(i,2),"-ymr.txt",sep=""),runNo=i,rwtNo=1,doMSY=TRUE,msyMaxIter=15000,msyTolConv=0.01,endStrat=0.301,stepStrat=0.001)

#out=runADMB("input26-ymr.txt",runNo=26,rwtNo=1,doMSY=TRUE)
#out=runADMB("input27-ymr.txt",runNo=27,rwtNo=1,doMSY=TRUE)
#out=runADMB("input28-ymr.txt",runNo=28,rwtNo=1,doMSY=TRUE)



