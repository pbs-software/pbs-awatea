#runSweave------------------------------2013-09-11
# Create and run customised Sweave files for Awatea runs.
# Updated 'runSweave.r' to parallel 'runADMB.r'  5/10/11
#-----------------------------------------------RH
runSweave = function( wd = getwd(), strSpp="XYZ",
		filename = "spp-area-00.txt",             # Name of Awatea .txt file in 'run.dir' to run
		runNo   = 1,
		rwtNo   = 0,
		running.awatea =0,                        # 0 if just loading previous '.rep'; 1 if rerunning Awatea
		Nsex    = 2,                              # if 1 then Unisex, if 2 Males & Females
		Ncpue   = 0,
		Nsurvey = 3,
		Snames  = paste("Ser",1:Nsurvey,sep=""),  # survey names (w/out spaces)
		SApos   = rep(TRUE,Nsurvey),              # surveys with age composition data
		delim   = "-",
		debug   = FALSE,
		locode  = FALSE,                          # if source as local code (for debugging)
		awateaPath="C:/Users/haighr/Files/Projects/ADMB/Coleraine",
		codePath="C:/Users/haighr/Files/Projects/R/Develop/PBSawatea/Authors/Rcode/develop"
	) {
	on.exit(setwd(wd))
	remove(list=setdiff(ls(1,all.names=TRUE),c("runMPD","runSweave","awateaCode","toolsCode")),pos=1)
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
		source(paste(codePath,"runSweaveMCMC.r",sep="/"),local=FALSE)
		source(paste(codePath,"plotFuns.r",sep="/"),local=FALSE)
		source(paste(codePath,"utilFuns.r",sep="/"),local=FALSE)
		source(paste(codePath,"menuFuns.r",sep="/"),local=FALSE)
		assign("importCol2",importRes,envir=.GlobalEnv)
	}
	cpue     = Ncpue > 0
	runNoStr = pad0(runNo,2)
	rwtNoStr = pad0(rwtNo,2)
	run.name = paste(strSpp,"run",runNoStr,sep="")
	run.dir  = paste(wd,run.name,sep="/")
	ext      = sapply(strsplit(filename,"\\."),tail,1)
	prefix   = substring(filename,1,nchar(filename)-nchar(ext)-1)
	prefix   = gsub(paste(delim,runNoStr,sep=""),"",prefix)        # get rid of superfluous run number in name
	model.name = paste(prefix,runNoStr,rwtNoStr,sep=".")
	mpdname  = paste("MPD",runNoStr,rwtNoStr,sep=".")
	mpd.dir  = paste(run.dir,mpdname,sep="/")  # directory where all the postscript crap happens
	if (file.exists(mpd.dir)) 
		setwd(mpd.dir)
	else {
		dir.create(mpd.dir); setwd(mpd.dir) }
	if (!file.exists("run-master.Snw"))
		file.copy(paste(system.file(package="PBSawatea"),"/snw/run-master.Snw",sep=""),wd)
	masterSweave = readLines(paste(ifelse(locode,codePath,wd),"run-master.Snw",sep="/"))
	tfile = masterSweave

	# First, get rid of those annoying comments and disabled code
	notcode = union(grep("^%",tfile),grep("^#",tfile))
	tfile = tfile[setdiff(1:length(tfile),notcode)]

	tfile = gsub("@cwd",wd,tfile)
	tfile = gsub("@model.name",model.name,tfile)
	tfile = gsub("@run.dir",run.dir,tfile)
	tfile = gsub("@fig.dir",mpd.dir,tfile)
	tfile = gsub("@running.awatea",running.awatea,tfile)
	tfile = gsub("@sppcode",strSpp,tfile)
	if (!locode) data(gfcode)
	tfile = gsub("@sppname", gfcode[is.element(gfcode$code3,strSpp),"name"],tfile)
#browser();return()
	packList(stuff=c("Snames"), target="PBSawatea")
	#if (exists("tput")) tput(Snames)
	snames = rep(Snames,Nsurvey)[1:Nsurvey] # enforce same number of names as surveys
	snames = gsub(" ","",snames)
	tfile = gsub("@surveys",paste(snames,collapse="\",\""),tfile)

	if (Nsex==1) {
		z0    = grep("@rmsex",tfile)
		tfile = tfile[setdiff(1:length(tfile),z0)]
		z1    = intersect(grep("Female",tfile),grep("sppfig",tfile))
		z2    = intersect(grep("Female",tfile),grep("twofig",tfile))
		tfile[z2] = sapply(tfile[z2],function(x){xx=strsplit(x,split="\\{"); paste(xx[[1]][c(1,2,4)],collapse="{")})
		z12 = c(z1,z2)
		tfile[z12] = gsub("twofig","sppfig",gsub("Female","Unisex",tfile[z12]))
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
	priorBites = c("logqvec\\.prior\\[1,]","muvec\\.prior\\[1,]","logvvec\\.prior\\[1,]","deltavec\\.prior\\[1,]")
	figBites   = c("survIndSer4-1","ageSurv","survRes","survAgeRes")
	cpueBites  = c("logqCPUE\\.prior\\[1,]","CPUE 1")
	survBites  = c("Survey 1") 
	# note assume only one method, otherwise need to expand "CAc 1"

	biteMe = function(infile, bites, N) {
		for (b in bites) {
			Nline = grep(b,infile)
			if (length(Nline)==0) next
			aline = infile[ Nline ]
			alines=NULL
			NApos = 0
			for ( i in 1:N) {
				NApos = NApos + as.numeric(SApos[i])
				if ((grepl("[Aa]ge",b) || grepl("CAs",b)) && !SApos[i]) next
				#if ((grepl("[Aa]ge",b) || grepl("CAs",b) || grepl("muvec",b)) && !SApos[i]) next
				#if (N==1 && !any(b==figBites)) alines=c(alines,gsub("\\[1,]","",aline))
				if (grepl("CAs",b) && SApos[i]){
					bline = gsub("1",i,aline)
					alines = c(alines, gsub(paste("\\[",i,"]",sep=""),paste("\\[",NApos,"]",sep=""),bline))
				}
				else if (N==1) alines=c(alines,gsub("\\[1,]","",aline))
				else alines=c(alines,gsub("1",i,aline))
			}
			infile=c(infile[1:(Nline-1)],alines,infile[(Nline+1):length(infile)])
		}
		return(infile)
	}
	tfile = biteMe(tfile,priorBites,Nsurvey)
	tfile = biteMe(tfile,figBites,Nsurvey)
	tfile = biteMe(tfile,cpueBites,Ncpue)
	tfile = biteMe(tfile,survBites,Nsurvey)
	if (any(SApos)) tfile = biteMe(tfile,"CAs 1",Nsurvey)
	else tfile = tfile[-grep("^CAs 1",tfile)]
#browser();return()

	for (i in 1:Nsurvey)
		tfile = gsub(paste("@survey",i,sep=""),snames[i],tfile)

	localHistory = paste(mpd.dir,"/runHistory.tex",sep="")
	if(file.exists(paste(wd,"/runHistory.tex",sep=""))) {
		is.history=file.copy(paste(wd,"/runHistory.tex",sep=""),localHistory,overwrite=TRUE)
	} else
		tfile = gsub("\\\\input","%\\\\input",tfile)

	if (length(grep("CUT HERE",tfile))>0)
		tfile = tfile[1:grep("CUT HERE",tfile)[1]]

	localName   = paste(run.name,"-",rwtNo,sep="")
	localSweave = paste(mpd.dir,"/",localName,".Snw",sep="")

	writeLines(tfile,con=localSweave)
	if (debug) { browser();return() }
	
	Sweave(localSweave)
	#mkpath="C:\\Miktex\\miktex\\bin\\"
	shell(cmd=paste("latex ",localName,".tex",sep=""),wait=TRUE)
	shell(cmd=paste("latex ",localName,".tex",sep=""),wait=TRUE)
	shell(cmd=paste("dvips ",localName,".dvi",sep=""),wait=TRUE)
	shell(cmd=paste("ps2pdf ",localName,".ps",sep=""),wait=TRUE)
	# if `pstopdf` not working, try:
	#system(cmd=paste("mgs -sDEVICE=pdfwrite -o ",localName,".pdf ",localName,".ps",sep=""),wait=TRUE)
	#http://tex.stackexchange.com/questions/49682/how-to-configure-ps2pdf-in-miktex-portable
	invisible(tfile) }
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~runSweave

#runMPD---------------------------------2012-07-25
# Wrapper to function 'runSweave' for MPDs.
#-----------------------------------------------RH
#runMPD = function(strSpp="XYZ",prefix=c("spp","area"), runs=1, rewts=0:6, cpue=FALSE, delim="-") {
runMPD = function(prefix=c("spp","area"), runs=1, rwts=0, ...) {
	# (...) pass in arguments specific to runSweave if different from the defaults:
	# Args: c(wd, cpue, strSpp, filename, runNo, rwtNo, running.awatea, Nsurvey, Snames, delim)
	# If prefix=NULL, filename will be taken from (...) or set to the default.
	dots = list(...)
	if (is.null(dots$delim)) delim="-" else delim=dots$delim
	for (i in runs) {
		if (!is.null(prefix)) filename=paste(paste(c(prefix,pad0(i,2)),collapse=delim),".txt",sep="")
		for (j in rwts) {
			runSweave(filename=filename, runNo=i, rwtNo=j, ...)
}	}	}

#runMPD(strSpp="POP",prefix=c("pop","3CD"),runs=3,rwts=1,cpue=TRUE,Nsurvey=2,Snames=c("NMFS Triennial","WCVI Synoptic"),SApos=c(F,T),delim="-",debug=T)
#runMPD(strSpp="POP",prefix=c("ymr","cst"),runs=29, rwts=0:1, cpue=FALSE, Nsurvey=5, Snames=c("GIG","QCSsyn","QCSshr","WCHGsyn","WCVIsyn"))
# runMPD(c(24, 26:28), 1,cpue=FALSE)     # No switch for EstM. AME doesn't have run25 in the new format
                                         # (but not using it anyway, MCMCs weren't good).
# runMPD(29, 0:6,cpue=FALSE)             # No switch for EstM. AME doesn't
# runMPD(36, 0:6, cpue=FALSE, Nsurvey=6) # No switch for EstM. AME doesn't

#=== ROL 5DE 2013 ===
#runMPD(strSpp="ROL",prefix=c("ROL","5CD"),runs=1,rwts=3,Nsex=1,Ncpue=2,Nsurvey=2,Snames=c("HS Assemblage","HS Synoptic"),SApos=c(T,T),locode=T)

#=== ROL 5AB 2013 ===
#runMPD(strSpp="ROL",prefix=c("ROL","5AB"),runs=7,rwts=3,Nsex=1,Ncpue=2,Nsurvey=1,Snames=c("QCS Synoptic"),SApos=c(TRUE),locode=TRUE)
#runMPD(strSpp="ROL",prefix=c("ROL","5AB"),runs=8,rwts=3,Nsex=1,Ncpue=2,Nsurvey=1,Snames=c("QCS Synoptic"),SApos=c(TRUE),locode=TRUE)

#=== ROL 5ABCD 2013 ===
#runMPD(strSpp="ROL",prefix=c("ROL","5ABCD"),runs=1,rwts=3,Nsex=1,Nsurvey=3,Ncpue=2,Snames=c("HS Assemblage","HS Synoptic","QCS Synoptic"),SApos=c(TRUE,TRUE,TRUE),locode=TRUE)

#=== SGR CST 2013 ===
#runMPD(strSpp="SGR",prefix=c("SGR","CST"),runs=1,rwts=3,Nsex=2,Ncpue=0,Nsurvey=6,Snames=c("Historic GB Reed","WCHG Synoptic","HS Synoptic","QC Sound Synoptic","US Triennial","WCVI Synoptic"),SApos=c(FALSE,TRUE,TRUE,TRUE,FALSE,TRUE),locode=TRUE)
#runMPD(strSpp="SGR",prefix=c("SGR","CST"),runs=2,rwts=3,Nsex=2,Ncpue=0,Nsurvey=6,Snames=c("Historic GB Reed","WCHG Synoptic","HS Synoptic","QC Sound Synoptic","US Triennial","WCVI Synoptic"),SApos=c(FALSE,TRUE,TRUE,TRUE,FALSE,TRUE),locode=TRUE)
