#runSweave------------------------------2013-09-06
# Create and run customised Sweave files for Awatea MCMC runs.
# Updated 'runSweave.r' to parallel 'runADMB.r'  5/10/11
# Updated 'runSweaveMCMC.r' to parallel 'runADMB.r'  5/10/11
#-----------------------------------------------RH
runSweaveMCMC = function(wd=getwd(), strSpp="XYZ",
		filename="spp-area-00.txt",           # Name of Awatea .txt file in 'run.dir' to run
		runNo = 1,
		rwtNo = 0,
		running.awatea=0,   # running.awatea=0 : load previous '.rep'; =1 : rerun Awatea
		Nsex    = 2,                              # if 1 then Unisex, if 2 Males & Females
		Ncpue   = 0,
		Nsurvey = 3,
		estM  = TRUE,
		mcsub = 1:1000,
		delim = "-",
		locode = FALSE,                          # if source as local code (for debugging)
		awateaPath = "C:/Users/haighr/Files/Projects/ADMB/Coleraine",
		codePath = "C:/Users/haighr/Files/Projects/R/Develop/PBSawatea/Authors/Rcode/develop",
		sexlab  = c("Females","Males")
	) {
	on.exit(setwd(wd))
	remove(list=setdiff(ls(1,all.names=TRUE),c("runMCMC","runSweaveMCMC","awateaCode","toolsCode")),pos=1)
	if (locode) { 
		getFile(gfcode,path=system.file("data",package="PBSawatea"))
		require(PBSmodelling, quietly=TRUE)
		require(gplots, quietly=TRUE)
		require(xtable, quietly=TRUE) 
		require(lattice, quietly=TRUE)
		require(scape, quietly=TRUE)     # Arni Magnusson's support functions for Awatea.
		require(scapeMCMC, quietly=TRUE) # Arni Magnusson's support functions for Awatea MCMC.
		require(gdata, quietly=TRUE)     # Data manipulation functions from CRAN.
		source(paste(codePath,"PBSscape.r",sep="/"),local=FALSE)
		source(paste(codePath,"runADMB.r",sep="/"),local=FALSE)
		source(paste(codePath,"runSweave.r",sep="/"),local=FALSE)
		source(paste(codePath,"plotFuns.r",sep="/"),local=FALSE)
		source(paste(codePath,"utilFuns.r",sep="/"),local=FALSE)
		source(paste(codePath,"menuFuns.r",sep="/"),local=FALSE)
		#assign("importCol2",importRes,envir=.GlobalEnv) # RH: removed importCol2 (2013-09-13)
	}
	cpue     = Ncpue > 0
	runNoStr = pad0(runNo,2)
	rwtNoStr = pad0(rwtNo,2)
	run.name = paste(strSpp,"run",runNoStr,sep="")
	run.dir  = paste(wd,run.name,sep="/")
	ext      = sapply(strsplit(filename,"\\."),tail,1)
	prefix   = substring(filename,1,nchar(filename)-nchar(ext)-1)
	#prefix   = gsub(runNoStr,"",prefix)      # get rid of superfluous run number in name
	prefix   = gsub(paste(delim,runNoStr,sep=""),"",prefix)        # get rid of superfluous run number in name
	model.name = paste(prefix,runNoStr,rwtNoStr,sep=".")
	mcname   = paste("MCMC",runNoStr,rwtNoStr,sep=".")
	mc.dir   = paste(run.dir,mcname,sep="/")  # Directory where all the postscript crap happens
	msyname  = paste("MSY",runNoStr,rwtNoStr,sep=".")
	msy.dir  = paste(mc.dir,msyname,sep="/")
	prjname  = paste("PRJ",runNoStr,rwtNoStr,sep=".")
	prj.dir  = paste(mc.dir,prjname,sep="/")
	if (file.exists(mc.dir)) 
		setwd(mc.dir)
	else {
		stop(paste("MCMC directory << ",mc.dir," >> does not exist.\n",sep="")) }
		#dir.create(mc.dir); setwd(mc.dir) }
	if (!file.exists("run-masterMCMC.Snw"))
		file.copy(paste(system.file(package="PBSawatea"),"/snw/run-masterMCMC.Snw",sep=""),wd)
	masterSweave = readLines(paste(ifelse(locode,codePath,wd),"run-masterMCMC.Snw",sep="/"))
	tfile = masterSweave
#browser();return()

	# First, get rid of those annoying comments and disabled code
	notcode = union(grep("^%",tfile),grep("^#",tfile))
	tfile = tfile[setdiff(1:length(tfile),notcode)]

	tfile = gsub("@cwd",wd,tfile)
	tfile = gsub("@model.name",model.name,tfile)
	tfile = gsub("@run.dir",run.dir,tfile)
	tfile = gsub("@fig.dir",mc.dir,tfile)
	tfile = gsub("@msy.dir",msy.dir,tfile)
	tfile = gsub("@prj.dir",prj.dir,tfile)
	tfile = gsub("@running.awatea",running.awatea,tfile)
	tfile = gsub("@mcsub",deparse(mcsub),tfile)
	tfile = gsub("@sppcode",strSpp,tfile)
	data(gfcode,package="PBSawatea")
	tfile = gsub("@sppname", gfcode[is.element(gfcode$code3,strSpp),"name"],tfile)
	mcmc  = read.table("params.pst",header=TRUE)
	ncol  = dim(mcmc)[2];  Nfigs = ceiling(ncol/6)

	if (Nsex==1) {
		z0    = grep("@rmsex",tfile)
		tfile = tfile[setdiff(1:length(tfile),z0)]
		z1    = intersect(grep("Female",tfile),grep("onefig",tfile))
		z2    = intersect(grep("Female",tfile),grep("twofig",tfile))
		tfile[z2] = sapply(tfile[z2],function(x){xx=strsplit(x,split="\\{"); paste(xx[[1]][c(1,2,4)],collapse="{")})
		z12 = c(z1,z2)
		tfile[z12] = gsub("twofig","onefig",gsub("Female","Unisex",tfile[z12]))
	} else {
		tfile = gsub("@rmsex ","",tfile) # assumes space after @rmsex for readability in `run-Master.Snw`
	}

	if (!cpue) {
		z0    = grep("@rmcpue",tfile)
		tfile = tfile[setdiff(1:length(tfile),z0)]
	} else {
		tfile = gsub("@rmcpue ","",tfile) # assumes space after @rmcpue for readability in `run-Master.Snw`
	}

	# Start expanding lines using bites
	SpriorBites = c("log_qsurvey_prior\\[1,]","surveySfull_prior\\[1,]","p_surveySfulldelta\\[1,]","log_surveyvarL_prior\\[1,]")
	CpriorBites = c("log_qCPUE_prior\\[1,]","p_Sfullest\\[1,]","p_Sfulldelta\\[1,]","log_varLest_prior\\[1,]")
	SpostBites  = c("currentMCMC\\$P\\$q_1","currentMCMC\\$P\\$mu_1","currentMCMC\\$P\\$Delta_1","currentMCMC\\$P\\$\"log v_1L\"")
	CpostBites  = c("currentMCMC\\$P\\$q_999","currentMCMC\\$P\\$mu_999","currentMCMC\\$P\\$Delta_999","currentMCMC\\$P\\$\"log v_999L\"")
	#MpriorBites = c("p_Sfullest\\[1,]","p_Sfulldelta\\[1,]","log_varLest_prior\\[1,]") # would bite with Nmethods (if need be)
	figBites =c("pairs1")

SApos=c(TRUE,TRUE)
	biteMe = function(infile, bites, N) {
		if (N==0) return(infile)
		for (b in bites) {
			Nline = grep(b,infile)
			if (length(Nline)==0) next
			aline = infile[ Nline ]
			alines=NULL
			NApos = 0
			for ( i in 1:N) {
				iline = aline
				NApos = NApos + as.numeric(SApos[i])
				#if ((grepl("[Aa]ge",b) || grepl("CAs",b)) && !SApos[i]) next
				#if ((grepl("[Aa]ge",b) || grepl("CAs",b) || grepl("muvec",b)) && !SApos[i]) next
				#if (N==1 && !any(b==figBites)) alines=c(alines,gsub("\\[1,]","",aline))
				if (any(b==c(CpriorBites,CpostBites)) && grepl("999",iline)) iline = gsub("999",Nsurvey+i,iline)
				if (any(b==figBites)) {
					if (i==2)      iline=gsub("\\{st}","{nd}",iline)
					else if (i==3) iline=gsub("\\{st}","{rd}",iline)
					else if (i>=4) iline=gsub("\\{st}","{th}",iline)
				}
				if (grepl("CAs",b) && SApos[i]){
					iline = gsub("1",i,iline)
					alines = c(alines, gsub(paste("\\[",i,"]",sep=""),paste("\\[",NApos,"]",sep=""),iline))
				}
				else if (N==1) alines=c(alines,gsub("\\[1,]","",iline))
				else alines=c(alines,gsub("1",i,iline))
			}
			infile=c(infile[1:(Nline-1)],alines,infile[(Nline+1):length(infile)])
		}
		return(infile)
	}
	tfile = biteMe(tfile,SpriorBites,Nsurvey)
	tfile = biteMe(tfile,CpriorBites,Ncpue)
	tfile = biteMe(tfile,figBites,Nfigs)
	tfile = biteMe(tfile,SpostBites,Nsurvey)
	tfile = biteMe(tfile,CpostBites,Ncpue)
#browser();return()
	#if (any(SApos)) tfile = biteMe(tfile,"CAs 1",Nsurvey)
	#else tfile = tfile[-grep("^CAs 1",tfile)]

	localHistory = paste(mc.dir,"/runHistory.tex",sep="")
	if(file.exists(paste(wd,"/runHistory.tex",sep=""))) {
		is.history=file.copy(paste(wd,"/runHistory.tex",sep=""),localHistory,overwrite=TRUE)
	} else
		tfile = gsub("\\\\input","%\\\\input",tfile)

	if (length(grep("CUT HERE",tfile))>0)
		tfile = tfile[1:grep("CUT HERE",tfile)[1]]

	localName   = paste(mcname,".Snw",sep="")
	localSweave = paste(mc.dir,"/",localName,sep="")

	writeLines(tfile,con=localSweave)
#browser();return()
	Sweave(localSweave)
	shell(cmd=paste("latex ",mcname,".tex",sep=""),wait=TRUE)
	shell(cmd=paste("latex ",mcname,".tex",sep=""),wait=TRUE)
	shell(cmd=paste("dvips ",mcname,".dvi",sep=""),wait=TRUE)
	shell(cmd=paste("ps2pdf ",mcname,".ps",sep=""),wait=TRUE)
	invisible(tfile) }
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~runSweaveMCMC

#runMCMC------------------------------ -2012-08-23
# Wrapper to function 'runSweaveMCMC' for MCMCs.
#-----------------------------------------------RH
runMCMC = function(strSpp="XYZ", prefix=c("spp","area"), runs=7, rewts=0:6, cpue=FALSE, estM=TRUE, delim="-", mcsub=1:1000) {
	for (i in runs) {
		for (j in rewts) {
			#runSweaveMCMC(filename=paste("input",pad0(i,2),"-ymr.txt",sep=""), runNo=i,rwtNo=j, cpue=cpue, estM=estM)
			runSweaveMCMC(strSpp=strSpp,filename=paste(paste(c(prefix,pad0(i,2)),collapse=delim),".txt",sep=""),runNo=i,rwtNo=j,cpue=cpue,estM=estM,mcsub=mcsub)
}	}	}

#runMCMC(strSpp="POP",prefix=c("pop","wcvi"),runs=29,rewts=1,cpue=FALSE,estM=TRUE)


#runMCMC(c(24:26),1,cpue=FALSE,estM=TRUE)
# runMCMC(c(27:28),1,cpue=FALSE,estM=FALSE)
# runMCMC(26,1,cpue=FALSE,estM=TRUE)
#runMCMC(29,1,cpue=FALSE,estM=TRUE) # The two main runs for
#runMCMC(30,1,cpue=FALSE,estM=FALSE) # YMR11 submission

# To check MSY convergence:
# test = msyCalc( dir=msy.dir, error.rep = 1 )

# If need to use a higher U, do
# sum(currentMSY$u > 0.2999)      # to see for how many (automate into function when time)
