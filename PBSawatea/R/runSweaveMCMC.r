#runSweave------------------------------2011-05-31
# Create and run customised Sweave files for Awatea MCMC runs.
# Updated 'runSweave.r' to parallel 'runADMB.r'  5/10/11
# Updated 'runSweaveMCMC.r' to parallel 'runADMB.r'  5/10/11
#-----------------------------------------------RH
runSweaveMCMC = function(wd=getwd(), cpue=FALSE, estM=TRUE, strSpp="YMR",
    filename="input25-ymr.txt",           # Name of Awatea .txt file in 'run.dir' to run
    runNo=25, rwtNo=0, running.awatea=0   # running.awatea=0 : load previous '.rep'; =1 : rerun Awatea
	) {
	on.exit(setwd(wd))
	remove(list=setdiff(ls(1,all=TRUE),c("runMCMC","runSweaveMCMC")),pos=1)
	#require(PBSmodelling, quietly=TRUE)
	#require(xtable, quietly=TRUE) 
	#require(lattice, quietly=TRUE)
	#require(scape, quietly=TRUE)             # Arni Magnusson's support functions for Awatea.
	#require(scapeMCMC, quietly=TRUE)         # Arni Magnusson's support functions for Awatea MCMC.
	#require(gdata, quietly=TRUE)             # Data manipulation functions from CRAN.

	#source("ymrScape.r",local=FALSE)
	#source("importRes.r",local=FALSE)
	#assign("importCol2",importRes,envir=.GlobalEnv)
	#source("plotDensPOP.r",local=FALSE)      # AME's version that modifies defaults to give better
                                             # recruitment and Bt plots
	#source("plotDensPOPpars.r",local=FALSE)  # AME's version to add MPD for params
	#source("plotTracePOP.r",local=FALSE)     # AME's version that adds on MPD
	#source("plotBVBnorm.r",local=FALSE)      # B/B0 and V/V0, as lattice so can use for
                                             # multiple runs. NOT as lattice now.

	runNoStr = pad0(runNo,2)
	rwtNoStr = pad0(rwtNo,2)
	run.name = paste(strSpp,"run",runNoStr,sep="")
	run.dir  = paste(wd,run.name,sep="/")
	ext      = sapply(strsplit(filename,"\\."),tail,1)
	prefix   = substring(filename,1,nchar(filename)-nchar(ext)-1)
	prefix   = gsub(runNoStr,"",prefix)      # get rid of superfluous run number in name
	model.name = paste(prefix,runNoStr,rwtNoStr,sep=".")
	mcname   = paste("MCMC",runNoStr,rwtNoStr,sep=".")
	mc.dir   = paste(run.dir,mcname,sep="/")  # Directory where all the postscript crap happens
	msyname  = paste("MSY",runNoStr,rwtNoStr,sep=".")
	msy.dir  = paste(mc.dir,msyname,sep="/")
	if (file.exists(mc.dir)) 
		setwd(mc.dir)
	else {
		dir.create(mc.dir); setwd(mc.dir) }
	masterSweave = readLines(paste(wd,"ymrrun-masterMCMC.Snw",sep="/"))
	tfile = gsub("@cwd",wd,masterSweave)
	tfile = gsub("@model.name",model.name,tfile)
	tfile = gsub("@run.dir",run.dir,tfile)
	tfile = gsub("@fig.dir",mc.dir,tfile)
	tfile = gsub("@msy.dir",msy.dir,tfile)
	tfile = gsub("@running.awatea",running.awatea,tfile)
	if (!estM){
		tfile = gsub("\"M_1\",","",tfile)
		tfile = gsub("\"M_2\",","",tfile) }
	if (!cpue) {
		rmcpue = sapply(c("CPUEser","CPUEfit"),function(x,y){
			x=paste("ymrfig\\{",x,sep="")
			z=grep(x,y)
			if (length(z)==0) 0 else z },y=tfile,simplify=TRUE)
		tfile=tfile[!is.element(1:length(tfile),rmcpue)] }
	localName   = paste(mcname,".Snw",sep="")
	localSweave = paste(mc.dir,"/",localName,sep="")
	if (length(grep("CUT HERE",tfile))>0)
		tfile = tfile[1:grep("CUT HERE",tfile)[1]]
	writeLines(tfile,con=localSweave)
	Sweave(localSweave)
	shell(cmd=paste("latex -interaction=nonstopmode ",gsub("\\.Snw$",".tex",localName),sep=""),wait=TRUE)
	shell(cmd=paste("latex -interaction=nonstopmode ",gsub("\\.Snw$",".tex",localName),sep=""),wait=TRUE)      # latex twice to get labels correct
	shell(cmd=paste("dvips -q ",gsub("\\.Snw$",".dvi",localName),sep=""),wait=TRUE)
	shell(cmd=paste("ps2pdf ",gsub("\\.Snw$",".ps",localName),sep=""),wait=TRUE)
	invisible() }
#------------------------------------runSweaveMCMC

#runMCMC--------------------------------2011-05-31
# Wrapper to function 'runSweaveMCMC' for MCMCs.
#-----------------------------------------------RH
runMCMC = function(runs=7, rewts=0:6, cpue=FALSE, estM=TRUE) {
	for (i in runs) {
		for (j in rewts) {
			runSweaveMCMC(filename=paste("input",pad0(i,2),"-ymr.txt",sep=""), runNo=i,rwtNo=j, cpue=cpue, estM=estM)
}	}	}

#runMCMC(c(24:26),1,cpue=FALSE,estM=TRUE)
# runMCMC(c(27:28),1,cpue=FALSE,estM=FALSE)
# runMCMC(26,1,cpue=FALSE,estM=TRUE)
#runMCMC(29,1,cpue=FALSE,estM=TRUE) # The two main runs for
#runMCMC(30,1,cpue=FALSE,estM=FALSE) # YMR11 submission

# To check MSY convergence:
# test = msyCalc( dir=msy.dir, error.rep = 1 )

# If need to use a higher U, do
# sum(currentMSY$u > 0.2999)      # to see for how many (automate into function when time)
