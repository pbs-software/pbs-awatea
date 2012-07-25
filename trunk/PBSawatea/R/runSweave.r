#runSweave------------------------------2012-07-19
# Create and run customised Sweave files for Awatea runs.
# Updated 'runSweave.r' to parallel 'runADMB.r'  5/10/11
#-----------------------------------------------RH
runSweave = function( wd = getwd(), cpue=FALSE, strSpp="XYZ",
		filename = "spp-area-00.txt",    # Name of Awatea .txt file in 'run.dir' to run
		runNo = 1,
		rwtNo = 0,
		running.awatea =0,                # 0 if just loading previous '.rep'; 1 if rerunning Awatea
		Nsurvey = 5,
		Snames = paste("Ser",1:Nsurvey,sep=""),
		delim="-"
		) {
	on.exit(setwd(wd))
	remove(list=setdiff(ls(1,all.names=TRUE),c("runMPD","runSweave")),pos=1)
	#require(PBSmodelling, quietly=TRUE)
	#require(xtable, quietly=TRUE) 
	#require(lattice, quietly=TRUE)
	#require(scape, quietly=TRUE)     # Arni Magnusson's support functions for Awatea.
	#require(scapeMCMC, quietly=TRUE) # Arni Magnusson's support functions for Awatea MCMC.
	#require(gdata, quietly=TRUE)     # Data manipulation functions from CRAN.

	#source("ymrScape.r",local=FALSE)
	#source("utilFuns.r",local=FALSE)
	#source("plotFuns.r",local=FALSE)
	#assign("importCol2",importRes,envir=.GlobalEnv)
	
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
	masterSweave = readLines(paste(wd,"run-master.Snw",sep="/"))
	tfile = gsub("@cwd",wd,masterSweave)
	tfile = gsub("@model.name",model.name,tfile)
	tfile = gsub("@run.dir",run.dir,tfile)
	tfile = gsub("@fig.dir",mpd.dir,tfile)
	tfile = gsub("@running.awatea",running.awatea,tfile)
	tfile = gsub("@sppcode",strSpp,tfile)
	data(gfcode)
	tfile = gsub("@sppname", gfcode[is.element(gfcode$code3,strSpp),"name"],tfile)
	snames = rep(Snames,Nsurvey)[1:Nsurvey] # enforce same number of names as surveys
	tfile = gsub("@surveys",paste(snames,collapse="\",\""),tfile)

#browser();return()
	
	if (!cpue) {
		rmcpue = sapply(c("CPUEser","CPUEfit"),function(x,y){
			x=paste("ymrfig\\{",x,sep="")
			z=grep(x,y)
			if (length(z)==0) 0 else z },y=tfile,simplify=TRUE)
		tfile=tfile[!is.element(1:length(tfile),rmcpue)] }
	priorBites = c("logqvec\\.prior\\[1,]","muvec\\.prior\\[1,]",
		"logvvec\\.prior\\[1,]","deltavec\\.prior\\[1,]")
	figBites =c("survIndSer4-1","ageSurv","survResSer1","survAgeResSer1")
	for (b in c(priorBites,figBites)) {
		Nline = grep(b,tfile)
		aline = tfile[ Nline ]
		alines=NULL
		for ( i in 1:Nsurvey) {
			if (Nsurvey==1 && !any(b==figBites)) alines=c(alines,gsub("\\[1,]","",aline))
			else alines=c(alines,gsub("1",i,aline))
		}
		tfile=c(tfile[1:(Nline-1)],alines,tfile[(Nline+1):length(tfile)])
	}
	for (i in 1:Nsurvey)
		tfile = gsub(paste("@survey",i,sep=""),snames[i],tfile)
	localName   = paste(run.name,"-",rwtNo,".Snw",sep="")
	localSweave = paste(mpd.dir,"/",localName,sep="")
	if (length(grep("CUT HERE",tfile))>0)
		tfile = tfile[1:grep("CUT HERE",tfile)[1]]
	writeLines(tfile,con=localSweave)
#browser();return()
	Sweave(localSweave)
	shell(cmd=paste("latex -interaction=nonstopmode ",gsub("\\.Snw$",".tex",localName),sep=""),wait=TRUE)
	shell(cmd=paste("dvips -q ",gsub("\\.Snw$",".dvi",localName),sep=""),wait=TRUE)
	shell(cmd=paste("ps2pdf ",gsub("\\.Snw$",".ps",localName),sep=""),wait=TRUE)
	invisible() }
#----------------------------------------runSweave

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

#runMPD(strSpp="POP",prefix=c("ymr","cst"),runs=29, rwts=0:1, cpue=FALSE, Nsurvey=5, Snames=c("GIG","QCSsyn","QCSshr","WCHGsyn","WCVIsyn"))

# runMPD(c(24, 26:28), 1,cpue=FALSE)     # No switch for EstM. AME doesn't have run25 in the new format
                                         # (but not using it anyway, MCMCs weren't good).
# runMPD(29, 0:6,cpue=FALSE)             # No switch for EstM. AME doesn't
# runMPD(36, 0:6, cpue=FALSE, Nsurvey=6) # No switch for EstM. AME doesn't


