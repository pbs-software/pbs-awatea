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
		Ncpue = 0,
		estM  = TRUE,
		mcsub = 1:1000,
		delim = "-",
		locode = FALSE                          # if source as local code (for debugging)
	) {
	on.exit(setwd(wd))
	remove(list=setdiff(ls(1,all.names=TRUE),c("runMCMC","runSweaveMCMC")),pos=1)
	if (locode) { 
		getFile(gfcode,path=system.file("data",package="PBSawatea"))
		require(PBSmodelling, quietly=TRUE)
		require(gplots, quietly=TRUE)
		require(xtable, quietly=TRUE) 
		require(lattice, quietly=TRUE)
		require(scape, quietly=TRUE)     # Arni Magnusson's support functions for Awatea.
		require(scapeMCMC, quietly=TRUE) # Arni Magnusson's support functions for Awatea MCMC.
		require(gdata, quietly=TRUE)     # Data manipulation functions from CRAN.
		source("PBSscape.r",local=FALSE)
		source("utilFuns.r",local=FALSE)
		source("plotFuns.r",local=FALSE)
		assign("importCol2",importRes,envir=.GlobalEnv)
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
	masterSweave = readLines(paste(wd,"run-masterMCMC.Snw",sep="/"))
	tfile = gsub("@cwd",wd,masterSweave)
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
	figBites =c("pairs1")
	for (b in c(figBites)) {
		Nline = grep(b,tfile)
		if (length(Nline)==0) next
		aline = tfile[ Nline ]
		alines=NULL
		for ( i in 1:Nfigs) {
			nline = aline
			if (i==2)      nline=gsub("\\{st}","{nd}",nline)
			else if (i==3) nline=gsub("\\{st}","{rd}",nline)
			else if (i>=4) nline=gsub("\\{st}","{th}",nline)
			alines=c(alines,gsub("1",i,nline))
		}
		tfile=c(tfile[1:(Nline-1)],alines,tfile[(Nline+1):length(tfile)])
	}
	

# Smarter Sweave will now deal with non-estimated parameters
#	if (!estM){
#		tfile = gsub("\"M_1\",","",tfile)
#		tfile = gsub("\"M_2\",","",tfile) }
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
	# Finally, get rid of those annoying comments and disabled code
	notcode = union(grep("^%",tfile),grep("^#",tfile))
	tfile = tfile[setdiff(1:length(tfile),notcode)]

	writeLines(tfile,con=localSweave)
#browser();return()
	Sweave(localSweave)
	shell(cmd=paste("latex ",localName,".tex",sep=""),wait=TRUE)
	shell(cmd=paste("latex ",localName,".tex",sep=""),wait=TRUE)
	shell(cmd=paste("dvips ",localName,".dvi",sep=""),wait=TRUE)
#browser();return()
	shell(cmd=paste("ps2pdf ",localName,".ps",sep=""),wait=TRUE)
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
