#runSweave------------------------------2018-04-04
# Create and run customised Sweave files for Awatea runs.
# Updated 'runSweave.r' to parallel 'runADMB.r'  5/10/11
#-----------------------------------------------RH
runSweave = function(
   wd = getwd(), 
   strSpp="XYZ",
   filename = "spp-area-00.txt",      ## Name of Awatea .txt file in 'run.dir' to run
   runNo   = 1,
   rwtNo   = 0,
   running.awatea =0,                 ## 0 if just loading previous '.rep'; 1 if rerunning Awatea
   Nsex    = 2,                       ## if 1 then Unisex, if 2 then Males & Females
   Ncpue   = 0,
   Nsurvey = 3,
   Ngear   = 1,                       ## number of commercial gear types
   NCAset  = 2,                       ## number of commercial catch-age-age plot sets (>1 when #CA years > 20)
   Snames  = paste0("Ser",1:Nsurvey), ## survey names (w/out spaces)
   SApos   = rep(TRUE,Nsurvey),       ## surveys with age composition data
   Cnames  = paste0("Gear",1:Ngear),  ## survey names (w/out spaces)
   CApos   = rep(TRUE,Ngear),         ## commercial gears with age composition
   delim   = "-",
   debug   = FALSE,
   locode  = FALSE,                   ## source this function as local code (for development)
   codePath = "C:/Users/haighr/Files/Projects/R/Develop/PBSawatea/Authors/Rcode/develop",
   sexlab  = c("Females","Males"),
   resdoc  = FALSE,                   ## is this build for a research document?
   redo.Graphs = TRUE,                ## recreate all the figures (.eps, .png)
   ptype = "png"                      ## plot type --  either "eps" or "png"
) {
	on.exit(setwd(wd))
	remove(list=setdiff(ls(1,all.names=TRUE),c("runMPD","runSweave","Rcode","Scode","qu","so",".First")),pos=1)
	if (locode) { 
		mess = c(
		"require(PBSmodelling, quietly=TRUE, warn.conflicts=FALSE)",
		"require(gplots, quietly=TRUE)",
		"require(xtable, quietly=TRUE)",
		"require(lattice, quietly=TRUE)",
		"require(scape, quietly=TRUE)",     # Arni Magnusson's support functions for Awatea.
		"require(plotMCMC, quietly=TRUE)", # Arni Magnusson's plot functions for Awatea MCMC.
		"require(gdata, quietly=TRUE)"     # Data manipulation functions from CRAN.
		)
		eval(parse(text=mess))
#browser();return()
		source(paste(codePath,"PBSscape.r",sep="/"),local=FALSE)
		source(paste(codePath,"runADMB.r",sep="/"),local=FALSE)
		source(paste(codePath,"runSweaveMCMC.r",sep="/"),local=FALSE)
		source(paste(codePath,"plotFuns.r",sep="/"),local=FALSE)
		source(paste(codePath,"utilFuns.r",sep="/"),local=FALSE)
		source(paste(codePath,"menuFuns.r",sep="/"),local=FALSE)
		#assign("importCol2",importRes,envir=.GlobalEnv) # RH: removed importCol2 (2013-09-13)
	}
	cpue     = Ncpue > 0
	sexlab   = rep(sexlab,Nsex)[1:Nsex]
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
#browser();return()
	if (file.exists(mpd.dir)) 
		setwd(mpd.dir)
	else {
		dir.create(mpd.dir); setwd(mpd.dir) }
	#if (!file.exists("run-master.Snw")) ## it will almost never exist in the MPD run directory
	## This way, someone only using the package has a fighting chance of modifying the Snw:
	if (!file.exists(paste(wd,"run-master.Snw",sep="/")))
		file.copy(paste(system.file(package="PBSawatea"),"/snw/run-master.Snw",sep=""),wd)
	masterSweave = readLines(paste(ifelse(locode,codePath,wd),"run-master.Snw",sep="/"))
	tfile = masterSweave

	# First, get rid excess lines, annoying comments, and disabled code
	if (length(grep("CUT HERE",tfile))>0)
		tfile = tfile[1:grep("CUT HERE",tfile)[1]]
	notcode = union(grep("^%",tfile),grep("^#",tfile))
	tfile = tfile[setdiff(1:length(tfile),notcode)]

	tfile = gsub("@cwd",wd,tfile)
	tfile = gsub("@model.name",model.name,tfile)
	tfile = gsub("@run.dir",run.dir,tfile)
	tfile = gsub("@fig.dir",mpd.dir,tfile)
	tfile = gsub("@running.awatea",running.awatea,tfile)
	tfile = gsub("@redo.Graphs",redo.Graphs,tfile)
	tfile = gsub("@sexlab",deparse(sexlab),tfile)
	tfile = gsub("@sppcode",strSpp,tfile)
	tfile = gsub("@ptype",ptype,tfile)
	if (locode) {
		if (any(strSpp==c("POP","pop","396"))) sppname = "Pacific Ocean Perch"
		else if (any(strSpp==c("YMR","ymr","440"))) sppname = "Yellowmouth Rockfish"
		else if (any(strSpp==c("ROL","rol","621"))) sppname = "Rock Sole"
		else if (any(strSpp==c("SGR","sgr","405"))) sppname = "Silvergray Rockfish"
		else if (any(strSpp==c("YTR","ytr","418"))) sppname = "Yellowtail Rockfish"
		else if (any(strSpp==c("RBR","rbr","401"))) sppname = "Redbanded Rockfish"
		else if (any(strSpp==c("ARF","arf","602"))) sppname = "Arrowtooth Flounder"
		else if (any(strSpp==c("WAP","wap","228"))) sppname = "Walleye Pollock"
		else if (any(strSpp==c("RSR","rsr","439"))) sppname = "Redstripe Rockfish"
		else sppname="Unspecified species"
	} else {
		data(gfcode,package="PBSawatea")
		sppname = gfcode[is.element(gfcode$code3,strSpp),"name"]
	}
	tfile = gsub("@sppname", sppname, tfile)
	tfile = gsub("@strSpp", strSpp, tfile)

	packList(stuff=c("Snames","SApos","Cnames","CApos","ptype"), target="PBSawatea")
	#if (exists("tput")) tput(Snames)
	snames = rep(Snames,Nsurvey)[1:Nsurvey] # enforce same number of names as surveys
	snames = gsub(" ","",snames)
	tfile  = gsub("@surveys",paste(snames,collapse="\",\""),tfile)
	if (Ncpue>0)
		Cnames = rep(Cnames,Ncpue)[1:Ncpue] # enforce same number of names as surveys
	cnames = gsub(" ","",Cnames)
	tfile  = gsub("@cpues",paste(cnames,collapse="\",\""),tfile)

#browser();return()
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
	
	# Remove catch-at-age lines if no catch-at-age data
	if (sum(CApos)==0) {
		z0    = grep("@rmCA",tfile)
		tfile = tfile[setdiff(1:length(tfile),z0)]
	} else {
		tfile = gsub("@rmCA ","",tfile) # assumes space after @rmCA for readability in `run-master.Snw`
	}
	if (sum(SApos)==0) {
		z0    = grep("@rmSA",tfile)
		tfile = tfile[setdiff(1:length(tfile),z0)]
	} else {
		tfile = gsub("@rmSA ","",tfile) # assumes space after @rmSA for readability in `run-master.Snw`
	}
	if (sum(CApos)==0 && sum(SApos)==0) {
		z0    = grep("@rmCSA",tfile)
		tfile = tfile[setdiff(1:length(tfile),z0)]
	} else {
		tfile = gsub("@rmCSA ","",tfile) # assumes space after @rmCA for readability in `run-master.Snw`
	}
	if (resdoc) {
		#tfile = gsub("@resdoc",TRUE,tfile)
		if (any(grepl("@rmresdoc",tfile))) {
			z0 = grep("@rmresdoc",tfile)
			tfile = tfile[setdiff(1:length(tfile),z0)] }
		if (strSpp=="ROL" && any(grepl("@rmROL",tfile))) {   # Kendra-specific removals
			z0 = grep("@rmROL",tfile)
			tfile = tfile[setdiff(1:length(tfile),z0)] }
	}
	else {
		#tfile = gsub("@resdoc",FALSE,tfile)
		tfile = gsub("@rmresdoc ","",tfile) # assumes space after @rmresdoc for readability in `run-Master.Snw`
		tfile = gsub("@rmROL ","",tfile)    # assumes space after @rmROL for readability in `run-Master.Snw`
	}

	## Deal with CA figures that have been split by selecting @rmCA1, @rmCA2, etc. that is appropriate based on NCAset argument
	tfile = tfile[-grep(paste0("@rmCA[",paste0(setdiff(0:9,NCAset),collapse=""),"]"),tfile)] ## get rid of @rmCA's without the NCAset suffix
	tfile = gsub(paste0("@rmCA",NCAset," "),"",tfile) ## assumes space after @rmCAN for readability in `run-master.Snw`
#browser();return()

	# Start expanding lines using bites
	# IMPORTANT: each element string below must be a unique match to a place in `run-Master.Smw'
	SpriorBites = c("logqvec\\.prior\\[1,]","muvec\\.prior\\[1,]","logvvec\\.prior\\[1,]","deltavec\\.prior\\[1,]")
	CpriorBites = c("Vy.mpd\\[1]","muC.prior\\[1,]","logvC.prior\\[1,]","deltaC.prior\\[1,]")
	figBites   = c("survIndSer4-1","onefig\\{survRes","twofig\\{ageSurv","onefig\\{survAgeResSer1}",
		"onefig\\{survAgeResSer1Female}", "onefig\\{survAgeResSer1Male}")
	cpueBites  = c("logqCPUE\\.prior\\[1,]","CPUE 1")
	survBites  = c("Survey 1") 
	#gearBits   = c("onefig\\{commAgeResSer1}") #change only first `1'
	gearBites  = c(#"Vy.mpd\\[1]","muC.prior\\[1,]","logvC.prior\\[1,]","deltaC.prior\\[1,]",
		"onefig\\{cpueRes",
		"onefig\\{ageCommFemaleSer1}","onefig\\{ageCommFemaleSer1A}","onefig\\{ageCommFemaleSer1B}",
		"onefig\\{ageCommMaleSer1}","onefig\\{ageCommMaleSer1A}","onefig\\{ageCommMaleSer1B}",
		"onefig\\{commAgeResSer1}","onefig\\{commAgeResSer1Female}","onefig\\{commAgeResSer1Male}")

	biteMe = function(infile, bites, N, CSpos=SApos, allsub=TRUE) { # bug fix: need to supply SApos or CApos
		if (N==0) return(infile)
		if (allsub) subfun =gsub else subfun=sub
		for (b in bites) {
			Nline = grep(b,infile)
			if (length(Nline)==0) next
			aline = infile[ Nline ]
			alines=NULL
			NApos = 0; 
			for ( i in 1:N) {
				NApos = NApos + as.numeric(CSpos[i])
				if ((grepl("[Aa]ge",b) || grepl("CAs",b) || grepl("CAc",b)) && !CSpos[i]) next
				if (grepl("CAs",b) && CSpos[i]){
					bline = subfun("1",i,aline)
					alines = c(alines, subfun(paste("\\[",i,"]",sep=""),paste("\\[",NApos,"]",sep=""),bline))
				}
				else if (N==1) alines=c(alines,subfun("\\[1,]","",aline))
				else alines=c(alines,subfun("1",i,aline))
			}
			infile=c(infile[1:(Nline-1)],alines,infile[(Nline+1):length(infile)])
#browser();return()
		}
		return(infile)
	}
	tfile = biteMe(tfile,SpriorBites,Nsurvey)
	tfile = biteMe(tfile,CpriorBites,Ngear)
	tfile = biteMe(tfile,figBites,Nsurvey)
	tfile = biteMe(tfile,cpueBites,Ncpue)
	tfile = biteMe(tfile,survBites,Nsurvey)
#browser();return()

	if (any(SApos)) tfile = biteMe(tfile,"CAs 1",Nsurvey)
	#else            tfile = tfile[-grep("^CAs 1",tfile)]  # alreday been removed above
	if (any(CApos)){
		tfile = biteMe(tfile,gearBites,Ngear,CSpos=CApos)
		tfile = biteMe(tfile,"CAc 1",Ngear,CSpos=CApos)
	}
	#else tfile = tfile[-grep("^CAc 1",tfile)]  # alreday been removed above

	tfile = gsub("@one","1",tfile)  # to restore true values of `1' in expanded lines

	for (i in 1:Nsurvey)
		tfile = gsub(paste("@survey",i,sep=""),snames[i],tfile)
	if (Ncpue>0) {
		for (i in 1:Ncpue)
			tfile = gsub(paste("@cpue",i,sep=""),cnames[i],tfile)
	}

	localHistory = paste(mpd.dir,"/runHistory.tex",sep="")
	if(file.exists(paste(wd,"/runHistory.tex",sep=""))) {
		is.history=file.copy(paste(wd,"/runHistory.tex",sep=""),localHistory,overwrite=TRUE)
	} else
		tfile = gsub("\\\\input","%\\\\input",tfile)


	localName   = paste(run.name,"-",rwtNo,sep="")
	localSweave = paste(mpd.dir,"/",localName,".Snw",sep="")
#browser();return()

	writeLines(tfile,con=localSweave)
	if (debug) { browser();return() }
	
#browser();return()
	Sweave(localSweave)
	#mkpath="C:\\Miktex\\miktex\\bin\\"
	if (ptype=="eps") {
		shell(cmd=paste("latex ",localName,".tex",sep=""),wait=TRUE)
		shell(cmd=paste("latex ",localName,".tex",sep=""),wait=TRUE)
		shell(cmd=paste("dvips ",localName,".dvi",sep=""),wait=TRUE)
		shell(cmd=paste("ps2pdf ",localName,".ps",sep=""),wait=TRUE)
	} else {
		shell(cmd=paste0("texify --pdf --synctex=1 --clean ",localName,".tex"),wait=TRUE)
	}
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

#=== YTR CST 2014 ===
#outdat=runMPD(strSpp="YTR",prefix=c("YTR","CST"),runs=2,rwts=1,Nsex=2,Ncpue=0,Nsurvey=7,Ngear=2,Snames=c("HS Synoptic","QC Sound Synoptic","WCVI Synoptic","Historic GB Reed","WCHG Synoptic","US Triennial","WCVI Shrimp"),SApos=c(T,T,T,F,F,F,F), Cnames=c("Bottom Trawl","Midwater Trawl"),locode=T)
#outdat=runMPD(strSpp="YTR",prefix=c("YTR","CST2F"),runs=5,rwts=2,Nsex=2,Ncpue=0,Nsurvey=6,Ngear=2,Snames=c("HS Synoptic","QC Sound Synoptic","WCVI Synoptic","Historic GB Reed","WCHG Synoptic","US Triennial"),SApos=c(T,T,T,F,F,F), Cnames=c("Bottom Trawl","Midwater Trawl"),locode=T)
#outMPD = runSweave(strSpp="YTR",filename="YTR-CST1F-05.txt",runNo=5,rwtNo=2,Nsex=2,Ncpue=0,Nsurvey=6,Ngear=1, Snames=c("HS Synoptic","QC Sound Synoptic","WCVI Synoptic","Historic GB Reed","WCHG Synoptic","US Triennial"), SApos=c(T,T,T,F,F,F), Cnames=c("Trawl"),locode=T)

#=== RBR CST 2014 ===
#outMPD=runSweave(strSpp="RBR",filename="RBR-CST2F-01.txt",runNo=1,rwtNo=1,Nsex=2,Ncpue=0,Nsurvey=8,Ngear=2,Snames=c("QC Sound Synoptic","WCVI Synoptic","QC Sound Shrimp","WCHG Synoptic","HS Synoptic","US Triennial","Historic GB Reed","IPHC Longline"),SApos=c(T,T,T,F,F,F,F,F),Cnames=c("Bottom Trawl","Longline"),locode=T)

#=== RBR CST 2014 ===
#outMPD=runMPD(strSpp="POP",prefix=c("POP","5ABC"),runs=1,rwts=1,Nsex=2,Ncpue=0,Nsurvey=3,Ngear=1,Snames=c("GIG Historical","QC Sound Synoptic","QC Sound Shrimp"),Cnames=c("Bottom Trawl"),SApos=c(T,T,F),locode=T)
