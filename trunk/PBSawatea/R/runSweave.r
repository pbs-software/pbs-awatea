#runSweave------------------------------2011-05-31
# Create and run customised Sweave files for Awatea runs.
# Updated 'runSweave.r' to parallel 'runADMB.r'  5/10/11
#-----------------------------------------------RH
runSweave = function( wd = getwd(), cpue=FALSE, strSpp="YMR",
		filename = "input25-ymr.txt",    # Name of Awatea .txt file in 'run.dir' to run
		runNo = 25,
		rwtNo = 0,
		running.awatea =0,                # 0 if just loading previous '.rep'; 1 if rerunning Awatea
		Nsurvey = 5
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
	#source("importRes.r",local=FALSE)
	#assign("importCol2",importRes,envir=.GlobalEnv)
	#source("plotDensPOP.r",local=FALSE)      # AME's version that modifies defaults to give better recruitment and Bt plots
	#source("plotDensPOPpars.r",local=FALSE)  # AME's version to add MPD for params
	#source("plotTracePOP.r",local=FALSE)     # AME's version that adds on MPD
	#source("plotBVBnorm.r",local=FALSE)     # B/B0 and V/V0, as lattice so can use for multiple runs. NOT as lattice now.
	
	runNoStr = pad0(runNo,2)
	rwtNoStr = pad0(rwtNo,2)
	run.name = paste(strSpp,"run",runNoStr,sep="")
	run.dir  = paste(wd,run.name,sep="/")
	ext      = sapply(strsplit(filename,"\\."),tail,1)
	prefix   = substring(filename,1,nchar(filename)-nchar(ext)-1)
	prefix   = gsub(runNoStr,"",prefix)        # get rid of superfluous run number in name
	model.name = paste(prefix,runNoStr,rwtNoStr,sep=".")
	mpdname  = paste("MPD",runNoStr,rwtNoStr,sep=".")
	mpd.dir  = paste(run.dir,mpdname,sep="/")  # directory where all the postscript crap happens
	if (file.exists(mpd.dir)) 
		setwd(mpd.dir)
	else {
		dir.create(mpd.dir); setwd(mpd.dir) }
	masterSweave = readLines(paste(wd,"ymrrun-master.Snw",sep="/"))
	tfile = gsub("@cwd",wd,masterSweave)
	tfile = gsub("@model.name",model.name,tfile)
	tfile = gsub("@run.dir",run.dir,tfile)
	tfile = gsub("@fig.dir",mpd.dir,tfile)
	tfile = gsub("@running.awatea",running.awatea,tfile)
	if (!cpue) {
		rmcpue = sapply(c("CPUEser","CPUEfit"),function(x,y){
			x=paste("ymrfig\\{",x,sep="")
			z=grep(x,y)
			if (length(z)==0) 0 else z },y=tfile,simplify=TRUE)
		tfile=tfile[!is.element(1:length(tfile),rmcpue)] }
	survBites = c("logqvec\\.prior\\[1,]","muvec\\.prior\\[1,]",
		"logvvec\\.prior\\[1,]","deltavec\\.prior\\[1,]",
		"survIndSer4-1")
	for (b in survBites) {
		Nline = grep(b,tfile)
		aline = tfile[ Nline ]
		alines=NULL
		for ( i in 1:Nsurvey) {
			alines=c(alines,gsub("1",i,aline))
		}
		tfile=c(tfile[1:(Nline-1)],alines,tfile[(Nline+1):length(tfile)])
	}
	localName   = paste(run.name,"-",rwtNo,".Snw",sep="")
	localSweave = paste(mpd.dir,"/",localName,sep="")
	if (length(grep("CUT HERE",tfile))>0)
		tfile = tfile[1:grep("CUT HERE",tfile)[1]]
	writeLines(tfile,con=localSweave)
	Sweave(localSweave)
	shell(cmd=paste("latex -interaction=nonstopmode ",gsub("\\.Snw$",".tex",localName),sep=""),wait=TRUE)
	shell(cmd=paste("dvips -q ",gsub("\\.Snw$",".dvi",localName),sep=""),wait=TRUE)
	shell(cmd=paste("ps2pdf ",gsub("\\.Snw$",".ps",localName),sep=""),wait=TRUE)
	invisible() }
#----------------------------------------runSweave

#runMPD----------------------------------2011-05-31
# Wrapper to function 'runSweave' for MPDs.
#-----------------------------------------------RH
runMPD = function(runs=1, rewts=0:6, cpue=FALSE) {
	for (i in runs) {
		for (j in rewts) {
			runSweave(filename=paste("input",pad0(i,2),"-ymr.txt",sep=""), runNo=i,rwtNo=j, cpue=cpue)
}	}	}

# runMPD(24:28, 0:6,cpue=FALSE)          # No switch for EstM. 
# runMPD(c(24, 26:28), 1,cpue=FALSE)     # No switch for EstM. AME doesn't have run25 in the new format
                                         # (but not using it anyway, MCMCs weren't good).
# runMPD(29, 0:6,cpue=FALSE)             # No switch for EstM. AME doesn't
# runMPD(36, 0:6, cpue=FALSE, Nsurvey=6) # No switch for EstM. AME doesn't


