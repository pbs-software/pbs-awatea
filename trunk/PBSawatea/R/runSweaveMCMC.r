#runSweave------------------------------2012-07-18
# Create and run customised Sweave files for Awatea MCMC runs.
# Updated 'runSweave.r' to parallel 'runADMB.r'  5/10/11
# Updated 'runSweaveMCMC.r' to parallel 'runADMB.r'  5/10/11
#-----------------------------------------------RH
runSweaveMCMC = function(wd=getwd(), cpue=FALSE, estM=TRUE, strSpp="XYZ",
    filename="spp-area-00.txt",           # Name of Awatea .txt file in 'run.dir' to run
    runNo=1, rwtNo=0, running.awatea=0,   # running.awatea=0 : load previous '.rep'; =1 : rerun Awatea
    delim="-"
	) {
	on.exit(setwd(wd))
	remove(list=setdiff(ls(1,all.names=TRUE),c("runMCMC","runSweaveMCMC")),pos=1)
	#require(PBSmodelling, quietly=TRUE)
	#require(xtable, quietly=TRUE) 
	#require(lattice, quietly=TRUE)
	#require(scape, quietly=TRUE)             # Arni Magnusson's support functions for Awatea.
	#require(scapeMCMC, quietly=TRUE)         # Arni Magnusson's support functions for Awatea MCMC.
	#require(gdata, quietly=TRUE)             # Data manipulation functions from CRAN.

	#source("PBSscape.r",local=FALSE)
	#source("utilFuns.r",local=FALSE)
	#source("plotFuns.r",local=FALSE)
	#assign("importCol2",importRes,envir=.GlobalEnv)

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
		dir.create(mc.dir); setwd(mc.dir) }
	if (!file.exists("run-masterMCMC.Snw"))
		file.copy(paste(system.file(package="PBSawatea"),"/snw/run-masterMCMC.Snw",sep=""),wd)
	masterSweave = readLines(paste(wd,"run-masterMCMC.Snw",sep="/"))
	tfile = gsub("@cwd",wd,masterSweave)
	tfile = gsub("@model.name",model.name,tfile)
	tfile = gsub("@run.dir",run.dir,tfile)
	tfile = gsub("@fig.dir",mc.dir,tfile)
	tfile = gsub("@msy.dir",msy.dir,tfile)
	tfile = gsub("@prj.dir",prj.dir,tfile)
	tfile = gsub("@running.awatea",running.awatea,tfile)
	tfile = gsub("@sppcode",strSpp,tfile)
	data(gfcode,package="PBSawatea")
	tfile = gsub("@sppname", gfcode[is.element(gfcode$code3,strSpp),"name"],tfile)

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
#browser();return()
	Sweave(localSweave)
	shell(cmd=paste("latex -interaction=nonstopmode ",gsub("\\.Snw$",".tex",localName),sep=""),wait=TRUE)
	shell(cmd=paste("latex -interaction=nonstopmode ",gsub("\\.Snw$",".tex",localName),sep=""),wait=TRUE)      # latex twice to get labels correct
	shell(cmd=paste("dvips -q ",gsub("\\.Snw$",".dvi",localName),sep=""),wait=TRUE)
	shell(cmd=paste("ps2pdf ",gsub("\\.Snw$",".ps",localName),sep=""),wait=TRUE)
	invisible() }
#------------------------------------runSweaveMCMC

#runMCMC------------------------------ -2012-07-18
# Wrapper to function 'runSweaveMCMC' for MCMCs.
#-----------------------------------------------RH
runMCMC = function(strSpp="XYZ", prefix=c("spp","area"), runs=7, rewts=0:6, cpue=FALSE, estM=TRUE, delim="-") {
	for (i in runs) {
		for (j in rewts) {
			#runSweaveMCMC(filename=paste("input",pad0(i,2),"-ymr.txt",sep=""), runNo=i,rwtNo=j, cpue=cpue, estM=estM)
			runSweaveMCMC(strSpp=strSpp,filename=paste(paste(c(prefix,pad0(i,2)),collapse=delim),".txt",sep=""),runNo=i,rwtNo=j,cpue=cpue,estM=estM)
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
