#runSweave------------------------------2014-08-26
# Create and run customised Sweave files for Awatea MCMC runs.
# Updated 'runSweave.r' to parallel 'runADMB.r'  5/10/11
# Updated 'runSweaveMCMC.r' to parallel 'runADMB.r'  5/10/11
#-----------------------------------------------RH
runSweaveMCMC = function(wd=getwd(), strSpp="XYZ",
   filename="spp-area-00.txt",  # name of Awatea .txt file in 'run.dir' to run
   runNo   = 1,
   rwtNo   = 0,
   running.awatea=0,            # running.awatea=0 : load previous '.rep'; =1 : rerun Awatea
   Nsex    = 2,                 # if 1 then Unisex, if 2 Males & Females
   Ncpue   = 0,
   Nsurvey = 3,
   SApos   = rep(TRUE,Nsurvey), # surveys with age composition data
   mcsub   = 1:1000,
   delim   = "-",
   locode  = FALSE,             # source this function as local code (for development)
   awateaPath = "C:/Users/haighr/Files/Projects/ADMB/Coleraine",
   codePath   = "C:/Users/haighr/Files/Projects/R/Develop/PBSawatea/Authors/Rcode/develop",
   histRP  = FALSE,             # historical reference points
   wpaper  = FALSE,             # working paper
   resdoc  = FALSE,             # research document
   redo.Graphs = TRUE,          # recreate all the figures (.eps, .wmf, .png)
   skip.last.year = TRUE        # remove last year of projections (set to FALSE for POP 5ABC in 2010)
	) {
	ciao = function(wd){setwd(wd);gc(verbose=FALSE)}
	on.exit(ciao(wd))
	remove(list=setdiff(ls(1,all.names=TRUE),c("runMCMC","runSweaveMCMC","awateaCode","toolsCode")),pos=1)
	if (locode) { 
		getFile(gfcode,path=system.file("data",package="PBSawatea"))
		mess = c(
		"require(PBSmodelling, quietly=TRUE, warn.conflicts=FALSE)",
		"require(gplots, quietly=TRUE)",
		"require(xtable, quietly=TRUE)",
		"require(lattice, quietly=TRUE)",
		"require(scape, quietly=TRUE)",     # Arni Magnusson's support functions for Awatea.
		"require(scapeMCMC, quietly=TRUE)", # Arni Magnusson's support functions for Awatea MCMC.
		"require(gdata, quietly=TRUE)"     # Data manipulation functions from CRAN.
		)
		eval(parse(text=mess))
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
	mpdname  = paste("MPD",runNoStr,rwtNoStr,sep=".")
	mpd.dir  = paste(run.dir,mpdname,sep="/") # MPD directory with figures
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

	# First, get rid excess lines, annoying comments, and disabled code
	if (length(grep("CUT HERE",tfile))>0)
		tfile = tfile[1:grep("CUT HERE",tfile)[1]]
	notcode = union(grep("^%",tfile),grep("^#",tfile))
	tfile = tfile[setdiff(1:length(tfile),notcode)]
#browser();return()

	# Bring in the results files (if they exist)
	rfiles=c("resultsMCMC","resultsMPDfigs","resultsMPDtabs","resultsMPD")
	for (r in rfiles) {
		if (any(grepl(paste("@",r,sep=""),tfile))) {
			rescode = grep(paste("@",r,sep=""),tfile)[1]
			resfile = paste(wd,"/",r,"-run",runNoStr,".tex",sep="")
			if (file.exists(resfile)) {
				rfile = readLines(resfile)
				tfile = c(tfile[1:(rescode-1)],rfile,tfile[(rescode+1):length(tfile)])
			}
			else 
				tfile = tfile[setdiff(1:length(tfile),rescode)]
		}
	}
#browser();return()
	tfile = gsub("@cwd",wd,tfile)
	tfile = gsub("@model.name",model.name,tfile)
	tfile = gsub("@run.dir",run.dir,tfile)
	tfile = gsub("@fig.dir",mc.dir,tfile)
	tfile = gsub("@mpd.dir",mpd.dir,tfile)
	tfile = gsub("@msy.dir",msy.dir,tfile)
	tfile = gsub("@prj.dir",prj.dir,tfile)
	tfile = gsub("@running.awatea",running.awatea,tfile)
	tfile = gsub("@redo.Graphs",redo.Graphs,tfile)
	tfile = gsub("@skip.last.year",skip.last.year,tfile)
	tfile = gsub("@mcsub",deparse(mcsub),tfile)
	tfile = gsub("@nsex",Nsex,tfile)
	tfile = gsub("@sppcode",strSpp,tfile)
	if (!locode) data(gfcode,package="PBSawatea")
	sppname = gfcode[is.element(gfcode$code3,strSpp),"name"]
	sppname = sapply(sapply(strsplit(sppname," "),function(x){paste(toupper(substring(x,1,1)),tolower(substring(x,2)),sep="")},simplify=FALSE),paste,collapse=" ")
	spplatin = gfcode[is.element(gfcode$code3,strSpp),"latin"]
	tfile = gsub("@sppname", sppname ,tfile)
	mcmc  = read.table("params.pst",header=TRUE)
	ncol  = dim(mcmc)[2];  Nfigs = ceiling(ncol/6)

	if (Nsex==1) {
		z0    = grep("@rmsex",tfile)
		tfile = tfile[setdiff(1:length(tfile),z0)]
		z1    = intersect(grep("Female",tfile),grep("onefig",tfile))
		z2    = intersect(grep("Female",tfile),grep("twofig",tfile))
		if (length(z2) > 0) {
			tfile[z2] = sapply(tfile[z2],function(x){xx=strsplit(x,split="\\{"); paste(xx[[1]][c(1,2,4)],collapse="{")})
			z12 = c(z1,z2)
			if (length(z1z2) > 0)
				tfile[z12] = gsub("twofig","onefig",gsub("Female","Unisex",tfile[z12]))
		}
	} else {
		tfile = gsub("@rmsex ","",tfile) # assumes space after @rmsex for readability in `run-MasterMCMC.Snw`
	}
	if (!cpue) {
		z0    = grep("@rmcpue",tfile)
		if (length(z0) > 0)
			tfile = tfile[setdiff(1:length(tfile),z0)]
	} else {
		tfile = gsub("@rmcpue ","",tfile) # assumes space after @rmcpue for readability in `run-MasterMCMC.Snw`
	}
	if (!histRP) {
		z0    = grep("@rmhrp",tfile)
		if (length(z0) > 0)
			tfile = tfile[setdiff(1:length(tfile),z0)]
	} else {
		tfile = gsub("@rmhrp ","",tfile) # assumes space after @rmhrp for readability in `run-MasterMCMC.Snw`
	}
	if (wpaper||resdoc) {
		resdoc = TRUE
		tfile = gsub("@resdoc",TRUE,tfile)
		if (any(grepl("@rmresdoc",tfile))) {
			z0 = grep("@rmresdoc",tfile)
			tfile = tfile[setdiff(1:length(tfile),z0)] }
		if (strSpp=="ROL" && any(grepl("@rmROL",tfile))) {   # Kendra-specific removals
			z0 = grep("@rmROL",tfile)
			tfile = tfile[setdiff(1:length(tfile),z0)] }
	}
	else {
		resdoc = FALSE
		tfile = gsub("@resdoc",FALSE,tfile)
		tfile = gsub("@rmresdoc ","",tfile) # assumes space after @rmresdoc for readability in `run-MasterMCMC.Snw`
		tfile = gsub("@rmROL ","",tfile) # assumes space after @rmresdoc for readability in `run-MasterMCMC.Snw`
	}

	# Start expanding lines using bites
	SpriorBites = c("log_qsurvey_prior\\[1,]","surveySfull_prior\\[1,]","p_surveySfulldelta\\[1,]","log_surveyvarL_prior\\[1,]")
	CpriorBites = c("p_Sfullest\\[1,]","p_Sfulldelta\\[1,]","log_varLest_prior\\[1,]")
	cpueBites   = c("log q_999")
	figBites    = c("onefig\\{pairs1\\}")

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
#if (length(iline)>1) {browser();return()}
				if (grepl("999",iline)){ iline = gsub("999",Nsurvey+i,iline) }#; browser()}
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
	tfile = biteMe(tfile,CpriorBites,1)
	if (cpue)
		tfile = biteMe(tfile,cpueBites,Ncpue)
#browser();return()
	tfile = biteMe(tfile,figBites,Nfigs)
	#tfile = biteMe(tfile,SpostBites,Nsurvey)
	#tfile = biteMe(tfile,CpostBites,max(Ncpue,1))
	#if (any(SApos)) tfile = biteMe(tfile,"CAs 1",Nsurvey)
	#else tfile = tfile[-grep("^CAs 1",tfile)]

#browser();return()
	localHistory = paste(mc.dir,"/runHistory.tex",sep="")
	if(file.exists(paste(wd,"/runHistory.tex",sep=""))) {
		is.history=file.copy(paste(wd,"/runHistory.tex",sep=""),localHistory,overwrite=TRUE)
	} else
		tfile = gsub("\\\\input","%\\\\input",tfile)

	# Final clean-up of empty lines ( cannot because paragraphs delineation is squashed)
	#tfile = tfile[setdiff(1:length(tfile),grep("^$",tfile))]

	mcname      = paste(mcname,ifelse(wpaper,"-Wpaper",ifelse(resdoc,"-ResDoc","")),sep="")
	localName   = paste(mcname,".Snw",sep="")
	localSweave = paste(mc.dir,"/",localName,sep="")

	writeLines(tfile,con=localSweave)
	Sweave(localSweave)
	shell(cmd=paste("latex ",mcname,".tex",sep=""),wait=TRUE)
	shell(cmd=paste("latex ",mcname,".tex",sep=""),wait=TRUE)
	shell(cmd=paste("dvips ",mcname,".dvi",sep=""),wait=TRUE)
	shell(cmd=paste("ps2pdf ",mcname,".ps",sep=""),wait=TRUE)
	invisible(tfile) }
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~runSweaveMCMC

#runMCMC------------------------------ -2013-11-12
# Wrapper to function 'runSweaveMCMC' for MCMCs. (not tested recently)
#-----------------------------------------------RH
runMCMC = function(prefix=c("spp","area"), runs=1, rwts=0, ...) {
	# (...) pass in arguments specific to runSweaveMCMC if different from the defaults:
	# If prefix=NULL, filename will be taken from (...) or set to the default.
	dots = list(...)
	if (is.null(dots$delim)) delim="-" else delim=dots$delim
	for (i in runs) {
		if (!is.null(prefix)) filename=paste(paste(c(prefix,pad0(i,2)),collapse=delim),".txt",sep="")
		for (j in rwts) {
			runSweaveMCMC(filename=filename, runNo=i, rwtNo=j, ...)
}	}	}

#runMCMC = function(strSpp="XYZ", prefix=c("spp","area"), runs=7, rewts=0:6, Nsex=2, Ncpue=0, Nsurvey=3, SApos=rep(TRUE,Nsurvey), delim="-", mcsub=1:1000) {
#	for (i in runs) {
#		for (j in rewts) {
#			runSweaveMCMC(strSpp=strSpp,filename=paste(paste(c(prefix,pad0(i,2)),collapse=delim),".txt",sep=""),runNo=i,rwtNo=j,Nsex=Nsex,Ncpue=Ncpue,Nsurvey=Nsurvey,SApos=SApos,mcsub=mcsub)
#}	}	}

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
