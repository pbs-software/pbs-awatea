##==============================================================================
## PBSawatea plot functions:
##  compB0               : compare reference points relative to B0
##  mochaLatte           : alternative to lattice plots (mockLattice)
##  panelBoxes           : Plot quantile plots using 'nchains' to delimit separate boxes.
##  panelChains..........: Plot cumulative fequency of 'nchains' by partitioning one trace
##  panelTraces..........: Plot sequential trace of MCMC samples with running median and (0.05, 0.95) quantiles
##  plotACFs             : plot ACFs for the estimated parameters
##  plotAges             : plot the MPD model fits to age data
##  plotB2               : [deprecated?] only gets called in menuFuns.r
##  plotBars             : barplots of specific year age proportions
##  plotBVBnorm          : plot spawnining & vulnerable biomass relative to unfished equilibrium
##  plotCI               : plot points with confidence intervals
##  plotCPUE             : plot CPUE and fit with error bars
##  plotDensPOP          : edited plotMCMC::plotDens function
##  plotDensPOPpars      : edited plotMCMC::plotDens for parameters
##  plotDensPOPparsPrior : adding the prior automatically
##  plotIndexNotLattice  : taking some of plt.idx, but doing plot.Index NOT as lattice
##  plotMeanAge          : plot obs. & exp. mean ages from comm. & survey C@A data
##  plotRmcmcPOP         : plot recruitment posteriors quantiles as one graph over time
##  plotSnail            : plot snail-trails for MCMC stock status
##  plotTracePOP         : plot traces with running median
##  plotTraj             : show all median trajectories (base+sens) in one figure
##  plotVBcatch          : plot vulnerable biomass trajectory and catch history
##  plt.biomass          : plot small biomass figures
##  plt.bubbles          : bubble plots of observed and fitted ages
##  plt.catch            : plot small catch figures
##  plt.cpue             : plot crude CPUE figure
##  plt.initagedev       : initial age deviations figure
##  plt.quantBio         : plot quantiles of reconstructed and projected biomass|recruits
##  plt.recdev           : log recruitment deviations figure 
##  plt.recdevacf        : auto-correlation function of the log recruitment deviations
##  quantBox             : redefine boxplot to show quantiles
##  plotBox              : modified boxplot with quantile whiskers
##  splineCPUE           : Fit spline curves through CPUE data to determine CV process error.

##==============================================================================

#compB0---------------------------------2018-07-09
# Compare reference points and criteria relative to B0.
#-----------------------------------------------RH
compB0=function(B, Mnams=NULL, ratios=c(0.4,0.8), 
   include=list(A1=TRUE, A2=TRUE, SSPM=TRUE, Bmsy=TRUE, Bt=TRUE),
   t.yr=2011, boxwidth=0.6, figgy=FALSE, width=12, height=9, 
   pngres=400, lang=c("e","f"),...) {

	oldpar = par(no.readonly=TRUE); oldpso = grDevices::ps.options()
	ciao = function() {
		par(oldpar)
		mess = paste("grDevices::ps.options(",paste(paste(names(oldpso),sapply(oldpso,deparse),sep="="),collapse=","),")",sep="")
		eval(parse(text=mess))
		gc(verbose=FALSE) }
	on.exit(ciao())

	nmods = length(B)       # number of model runs
	mnams = names(B)        # default model names
	if (is.null(mnams)) { mnams = paste("mod",1:nmods,sep=""); names(B) = mnams }
	nBmsy = length(ratios)  # number of Bmsy ratios
	uBmsy = rep(unlist(include["Bmsy"]),nBmsy); names(uBmsy) = paste(ratios,names(uBmsy),sep="") # use Bmsy?
	include[["Bmsy"]] <- uBmsy
	nBars = unlist(sapply(include[1:3],function(x){x[x]}))
	nBoxs = unlist(sapply(include[4:5],function(x){x[x]}))
	names(nBoxs) = gsub("^Bmsy\\.","",names(nBoxs))
	mBoxs = rep(nBoxs,nmods); names(mBoxs) = paste(rep(mnams,each=length(mBoxs)/nmods),names(mBoxs),sep="~")
	nBarBox = c(nBars,mBoxs); 

	BarBox = as.list(nBarBox)                                                      # list object for bar and box data
	namBB = names(nBarBox)                                                         # names of bars and boxes that will be plotted
	namBar = namBB[is.element(namBB,c("A1","A2","SSPM"))]                          # namse of bars that will be plotted
	namBox = namBB[!is.element(namBB,namBar)]                                      # names of boxplots
	Bars = list(A1=c(0.3,0.5),A2=c(0.5,0.7),SSPM=c(0.2,0.4,0.5))                   # COSEWIC criteria A1 and A2, Schaefer surplus proction model
	BarBox[namBar] = Bars[namBar]                                                  # populate master list with bar information
	BmsyB0 = sapply(B,function(x){x[["Bmsy.MCMC"]]/x[["B0.MCMC"]]},simplify=FALSE) # sample ratios Bmsy/B0
	BtB0   = sapply(B,function(x){x[["Bt.MCMC"]]/x[["B0.MCMC"]]},simplify=FALSE)   # sample ratios Bt/B0
	for (i in namBox) {
		ii = strsplit(i,split="~")  # (model name, ratio of B0)
		imod = ii[[1]][1]; irat = ii[[1]][2]
		if (irat == "Bt")
			BarBox[[i]] = BtB0[[imod]]                                               # populate master list with Bt/B0
		if (grepl("Bmsy$",irat)) {
			rat = as.numeric(strsplit(irat,split="B")[[1]][1])
			BarBox[[i]] = rat * BmsyB0[[imod]]                                       # populate master list with ratios of Bmsy/B0
		}
	}
	xBox = BarBox; xBox[namBar] = NA  # for boxplots only

	if (is.null(Mnams)) {
		Mnams = rep("",nmods); names(Mnams) = mnams; ratsep="" }
	else {
		if (length(Mnams)!=nmods) 
			stop("if 'Mnams' is specified, it must be a string vector equal to the number of models")
		names(Mnams) = mnams; ratsep = "\n" }
	for (i in mnams)
		names(xBox) <- gsub(paste(i,"~",sep=""),paste(Mnams[i],ratsep,sep=""),names(xBox))
	names(xBox) <- gsub("Bt",paste("B",t.yr,sep=""),names(xBox))

	xlim=c(0.25,length(xBox)+0.75); ylim=c(0,1)
	medcol=whiskcol=staplecol=c("red","blue","black"); medcol[3]="green4"
	nClr =  sapply(include,sum); CLRS=clrs=character(0)
	if (nClr["Bmsy"]>0) {
		CLRS = c(CLRS,c("red","blue","black")[1:nClr["Bmsy"]])
		clrs = c(clrs,c("pink","lightblue1","gainsboro")[1:nClr["Bmsy"]]) }
	if (nClr["Bt"]>0) {
		CLRS = c(CLRS,c("black")[1:nClr["Bt"]])
		clrs = c(clrs,c("moccasin")[1:nClr["Bt"]]) }
	CLRS = c(rep("white",length(nBars)),rep(CLRS,nmods))
	clrs = c(rep("white",length(nBars)),rep(clrs,nmods))
	medcol=whiskcol=staplecol=CLRS; #medcol[3]="green4"
	boxfill = clrs
	dots = list(...); unpackList(dots)
	spp  = attributes(B)$spp
	if (!exists("spp")) spp="UNK"
	else if (is.null(spp) || spp=="") spp="UNK"

	#--------------PLOT ROUTINE----------------
	if (figgy) figout = c("eps","pdf","png","wmf","win") else figout="win"
	fout = fout.e = paste("CompB0-",spp,"-(",paste(names(xBox),collapse=","),")",sep="")

	for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
		changeLangOpts(L=l)
		fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
		for (f in figout) {
			if (f=="png") {
				png(paste(fout,".png",sep=""), units="in", res=pngres, width=width, height=height)
			} else if (f=="eps") {
				grDevices::ps.options(horizontal = TRUE)
				postscript(file=paste(fout,".eps",sep=""),width=width,height=height,fonts="mono")
			} else if (f=="pdf") {
				grDevices::ps.options(horizontal = TRUE)
				pdf(file=paste(fout,".pdf",sep=""),width=width,height=height,fonts="mono")
			} else if (f=="wmf" && .Platform$OS.type=="windows") {
				do.call("win.metafile",list(filename=paste0(fout,".wmf"), width=width, height=height))
			}
			par(mar=c(3.5,5,0.5,1),cex=ifelse(f%in%c("png","eps"),1,1.2),xaxs="i")
			quantbox(xBox,xlim=xlim,ylim=ylim,yaxs="i",las=1,xaxt="n",yaxt="n",xlab="",ylab="",
				pars=list(boxwex=boxwidth,medlwd=2,whisklty=1,medcol=medcol,staplecol=staplecol,whiskcol=whiskcol,boxfill=boxfill,...)) 
			ypos = par()$usr[4]-.025*diff(par()$usr[3:4])
			cex.txt = ifelse(f%in%"win",1.0,1.2)
			xaxislab = names(xBox)
			if (ratsep!="\n") { # if Mnams is supplied, they cannot be easily incorporated into an expression
				xaxislab = paste("expression(",paste(xaxislab,collapse=","),")",sep="")
				xaxislab = gsub("1\\*Bmsy","Bmsy",gsub("Bmsy","*Bmsy",xaxislab))
				xaxislab = gsub("2011","[2011]",gsub("msy","[MSY]",xaxislab))
				xaxislab = gsub("B","italic(B)",xaxislab) 
			}
			mess = paste(c("axis(1,at=1:length(xBox),labels=",deparse(xaxislab),",tick=FALSE,padj=0.5,mgp=c(2,0.75,0),cex.axis=1.2)"),collapse="")
			if (ratsep!="\n")  mess = gsub("\\\"","",mess)
			eval(parse(text=mess))  # plaster the mess along the x-axis
			axis(2,at=seq(0,1,0.1),mgp=c(2,0.75,0),las=1,cex.axis=1.2)
			if (any(nBars)) {
				for (i in 1:length(nBars)) {
					ii = namBar[i]
					ylev = BarBox[[ii]]; nlev = length(ylev)
					bxw = 0.5 * boxwidth
					xi = rep(c(i - bxw, i + bxw, NA),nlev)
					yi = as.vector(sapply(as.list(ylev),function(x){c(rep(x,2),NA)}))
					if (ii=="SSPM") ymax=median(ylev) else ymax=max(ylev)
					lines(c(c(i,i)- bxw,NA,c(i,i)+ bxw),c(0,ymax,NA,0,ymax),lwd=2,col="grey") # grey sides of bar
					if (ii=="SSPM") {
						for (j in 1:2) lines(c(i-bxw,i+bxw),rep(ylev[j],2),lwd=3,col=switch(j,"red","blue","green4")) 
						text(i,0.1,ifelse(l=="f","critique","Critical"),srt=90,cex=1.4,col="red")
						text(i,0.3,ifelse(l=="f","pr\u{00E9}cautionneux","Cautious"),srt=90,cex=1.4,col="blue") 
						text(i,ypos,ifelse(l=="f","Schaefer\nmod\u{00E8}le de\nproduction\nexc\u{00E9}dentaire","Schaefer\nsurplus\nproduction\nmodel"),cex=cex.txt,adj=c(.5,1))
						abline(v = i + 0.5,lty=2,col="grey",lwd=2)
					}
					else {
						lines(xi,yi,lwd=3) # vertical zone delimiters
						text(i,0.15,ifelse(l=="f","endanger","Endangered"),srt=90,cex=1.4,col="chocolate3")
						text(i,ifelse(ii=="A1",0.4,0.6),ifelse(l=="f","menac\u{00E9}es","Threatened"),srt=90,cex=1.4,col="purple")
						nCOS = nBars[intersect(names(nBars),c("A1","A2"))]
						text(sqrt(sum(nCOS)),ypos,"COSEWIC\ncriteria",cex=cex.txt,adj=c(.5,1))
						if (names(rev(nCOS))[1]==ii) abline(v = i + 0.5,lty=2,col="grey",lwd=2)
					}
				}
			}
			if (l=="f")
				mess =paste(c("mtext(expression(paste(\"Crit\u{00E8}res de r\u{00E9}f\u{00E9}rence et points relatifs \u{00E0} 	",ifelse(f%in%c("win","wmf"),"  ",""),
					"\",italic(B)[0],sep=\"\")),side=2,line=3.25,cex=",ifelse(f%in%c("win","wmf"),1.75,1.75),")"),collapse="")
			else
				mess =paste(c("mtext(expression(paste(\"Reference criteria and points relative to ",ifelse(f%in%c("win","wmf"),"  ",""),
					"\",italic(B)[0],sep=\"\")),side=2,line=3.25,cex=",ifelse(f%in%c("win","wmf"),1.75,1.75),")"),collapse="")
			eval(parse(text=mess)) # plaster the mess along the y-axis
			
			abline(v = length(nBars)+(1:(nmods-1))*length(nBoxs)+0.5,lty=2,col="grey",lwd=2)
			nB = sum(nClr[c("Bmsy","Bt")])
			modlab = length(nBars)+median(1:nB)+nB*((1:nmods)-1)
			for (i in 1:nmods)
				if (l=="f")
					text(modlab[i],ypos,paste(switch(i,"estimer","fixer"),"la mortalit\u{00E9} naturelle"),cex=cex.txt,adj=c(.5,1))
				else
					text(modlab[i],ypos,paste(switch(i,"Estimating","Fixing"),"natural mortality"),cex=cex.txt,adj=c(.5,1))
			box()
			if (f!="win") dev.off()
		} ## end f (figout) loop
	}; eop()
	invisible(list(BarBox=BarBox,xBox=xBox)) 
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~compB0


## mochaLatte---------------------------2019-05-21
##  An alternative to lattice plots (mockLattice)
##----------------------------------------------RH
mochaLatte = function(dat, xfld, yfld, ffld, panel,
   strip=list(col=lucent("black",0.5), bg=lucent("moccasin",0.5), height=0.1, cex=1.4), 
   fn.ylim=range, ...)
{
	opar = par(no.readonly=TRUE)
	on.exit(par(opar))
	options(scipen=5)

	getlab = function(lim, p=0.1) {
		tck = pretty(lim, n=6)
		if (p<=0) return(as.character(tck))
		lim = lim + (diff(lim) * p * c(1,-1))
		pos = tck>=lim[1] & tck<=lim[2]
		sho = rep("",length(tck))
		sho[pos] = as.character(tck[pos])
		return(list(tck=tck, sho=sho))
	}
	facs = as.character(unique(dat[,ffld]))
	nfac = length(facs)
	rc   = .findSquare(nfac)
	#if (rc[1]>rc[2]) rc = rev(rc)  ## more visually pleasing when more columns than rows
	if (rc[1]<rc[2]) rc = rev(rc)  ## more visually pleasing when more columns than rows
	strip$cex = strip$cex * (rc[1]^(-0.2))
	dots = list(...)
	if (is.null(dots$mar)) mar=c(3,3,0.5,0.5) else mar = dots$mar
	if (is.null(dots$oma)) oma=c(4,4,0.5,1)   else oma = dots$oma
	if (is.null(dots$mgp)) mgp=c(2,0.5,0)     else mgp = dots$mgp
#browser();return()
	par(mfrow=rc, mar=mar, oma=oma, mgp=mgp)
	for (i in facs) {
		idat = dat[is.element(dat[,ffld],i),]
		if (is.null(dots$xlim)) xlim = range(idat[,xfld], na.rm=TRUE) else xlim = dots$xlim
		if (is.null(dots$ylim)){
			yval = idat[,yfld]
			names(yval) = idat[,xfld]
			ylim = round(fn.ylim(yval),5)
		} else {
			ylim = dots$ylim
		}
#browser(); return()
		yticks    = getlab(ylim, p=0.05)
		ylim[2]   = ylim[2] + diff(ylim)*strip$height ## add space for latte foam
		hzero     = mar[2]==0  ## plots joined horizontally
		vzero     = mar[1]==0  ## plots joined vertically
		#do.call("plot", c(list(x=0, y=0, type="n", xlim=xlim, ylim=ylim, xaxt=ifelse(vzero,"n","n"), yaxt=ifelse(hzero,"n","n"), xlab="", ylab=""), dots[setdiff(names(dots), c("xlim","ylim","xlab","ylab","xfac","yfac"))]))
		exclude = c("xlim","ylim","xlab","ylab","xfac","yfac","outline")
		evalCall(plot, c(list(x=0, y=0, type="n", xlim=xlim, ylim=ylim, xaxt=ifelse(vzero,"n","n"), yaxt=ifelse(hzero,"n","n"), xlab="", ylab=""), dots[setdiff(names(dots), exclude)]), checkdef=T, checkpar=T)
		#dots[setdiff(names(dots), c("xlim","ylim","xlab","ylab","xfac","yfac"))]), checkdef=T, checkpar=T)
#browser(); return()
		if ((!vzero&&!hzero) || vzero || (hzero && par()$mfg[2]==1)){
			exclude = c("xlim","ylim","xfac","yfac","outline")
			#do.call("axis", c(list(side=2, at=yticks[["tck"]], labels=yticks[["sho"]]), dots[setdiff(names(dots),c("xlim","ylim","xfac","yfac"))]))
			evalCall(axis, c(list(side=2, at=yticks[["tck"]], labels=yticks[["sho"]]), dots[setdiff(names(dots),exclude)]), checkdef=T, checkpar=T)
		}
		if ((!vzero&&!hzero) || (hzero && !vzero) || i%in%rev(facs)[1:rc[2]] ) {
			if (!is.null(dots$xfac)) {
				xticks = unique(dat[,xfld])
				exclude = c("xlim","ylim","xfac","yfac", "outline")
				#do.call("axis", c(list(side=1, at=xticks, labels=dots$xfac), dots[setdiff(names(dots),c("xlim","ylim","xfac","yfac"))]))
				evalCall(axis, c(list(side=1, at=xticks, labels=dots$xfac), dots[setdiff(names(dots),exclude)]), checkdef=T, checkpar=T)
			} else {
				xticks = getlab(xlim, p=0.025)
#browser();return()
				exclude = c("xlim","ylim","xfac","yfac","exclude")
				#do.call("axis", c(list(side=1, at=xticks[["tck"]], labels=xticks[["sho"]]), dots[setdiff(names(dots),c("xlim","ylim","xfac","yfac"))]))
				evalCall(axis, c(list(side=1, at=xticks[["tck"]], labels=xticks[["sho"]]), dots[setdiff(names(dots),exclude)]), checkdef=T, checkpar=T)
			}
		}
		panel(x=idat[,xfld], y=idat[,yfld], dots[setdiff(names(dots),c("xfac","yfac"))])
		#legend("topleft", legend=i, x.intersp=0, box.col=strip$col, bg=strip$bg)
		strip$xbox = par()$usr[c(1,1,2,2)]
		strip$ybox = rep(par()$usr[4],4)
		strip$yoff = c(-diff(par()$usr[3:4])*strip$height,0,0,-diff(par()$usr[3:4])*strip$height)
		strip$ybox = strip$ybox + strip$yoff
		polygon(strip$xbox,strip$ybox,col=strip$bg, border=strip$col)
		text(x=mean(strip$xbox), y=mean(strip$ybox), labels=i, cex=strip$cex)
#if (i=="M_1") {browser();return()}
		box()
	}
	if (!is.null(dots$xlab))
		mtext(text=dots$xlab, side=1, line=par()$oma[1]*0.6, outer=TRUE, cex=ifelse(is.null(dots$cex.lab),1.5,dots$cex.lab))
	if (!is.null(dots$ylab))
		mtext(text=dots$ylab, side=2, line=par()$oma[2]*0.4, outer=TRUE, cex=ifelse(is.null(dots$cex.lab),1.5,dots$cex.lab))
#browser();return()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~mochaLatte


## panelBoxes---------------------------2019-11-24
##  Plot quantile plots using 'nchains' to delimit separate boxes.
##  mcmc=data.frame e.g, 'currentMCMC$P' from object created by 'importMCMC'.
##  Very difficult to manipulate trellis plots (RH)
## -----------------------------------------------
panelBoxes=function (mcmc, nchains=9, pdisc=0, 
   axes=FALSE, same.limits=FALSE, between=list(x=axes,y=axes), div=1,
   log=FALSE, base=10, main=NULL, xlab=NULL, ylab=NULL, xlim=c(0.25,9.75), 
   boxfill=paste0(rep(c("cyan","green","coral"),each=3), rep(1:3,3)),
   cex.main=1.2, cex.lab=1.2, cex.strip=0.9, cex.axis=0.9, las=0, 
   tck=0.4, tick.number=5, xfac=paste0("B",1:nchains), outline=TRUE,
	lang="e", ...)
{
	panel.box <- function(x, y, ...) {
		dots = list(...)[[1]]
		unpackList(dots)
		if (is.null(dots$outline)) outline = TRUE
		## xlim and ylim determined by 'mochaLatte'
		#if (is.null(dots$xlim)) xlim = range(x,na.rm=TRUE)
		#if (is.null(dots$ylim)) ylim = if (outline) range(y,na.rm=TRUE) else quantile(y, quants5[c(1,5)])
		#if (is.null(dots$xfac)) xfac = unique(x)
		basecol = "slategray"
		#boxfill = paste0(rep(c("cyan","green","coral"),each=3)
		boxpars = list(boxwex=0.5, boxfill=boxfill, boxcol=basecol, outpch=3, outcex=0.3, outcol=lucent(basecol,0.25), medlwd=1, whisklty=1, whiskcol=basecol)
		chainlink = rep(1:nchains,ff)
		chainbox  = split(y,x)
		#quantbox(chainbox, add=T, xaxt="n", yaxt="n", pars=boxpars, ...)
		mess =  sapply(dots,deparse)
		mess = paste0(paste0(names(mess),"=",mess),collapse=",") 
		messy = paste0("list(add=TRUE, xaxt=\"n\", yaxt=\"n\", pars=boxpars, ", mess,")")
		argos = eval(parse(text=messy))
		#do.call("quantbox", args=list(x=chainbox, add=TRUE, xaxt="n", yaxt="n", pars=boxpars, deparse(mess)) )
		do.call("quantbox", args=c(x=list(chainbox),argos) )
		#evalCall(quantbox, args=argos, checkdef=T, checkpar=T )
#browser(); return()
	}

	if (pdisc>0 && pdisc<1)
		mcmc = mcmc[(round(pdisc*nrow(mcmc))+1):nrow(mcmc),]  # get rid of the first 'pdisc' (e.g., 10%)
	relation <- if (same.limits) "same" else "free"
	if (is.null(dim(mcmc))) {
		mcmc.name <- rev(as.character(substitute(mcmc)))[1]
		mcmc <- matrix(mcmc, dimnames=list(NULL, mcmc.name))
	}
	mcmc <- if (log) 
		log(mcmc/div, base=base)
	else mcmc/div
	mcmc <- as.data.frame(mcmc)
	n <- nrow(mcmc)
	ff = rep(round(n/nchains),nchains-1)
	ff = c(ff,n-sum(ff))
	p <- ncol(mcmc)
	dat <- data.frame(Factor=ordered(rep(names(mcmc), each=n), 
		names(mcmc)), Draw=rep(1:n, p), Chain=rep(rep(1:nchains,ff),p), Value=as.vector(as.matrix(mcmc)))
	dat$Index = paste(dat$Factor,dat$Chain,sep="-")

	fn.ylim =
		if (outline) function(x){range(x, na.rm=TRUE)} 
		else         function(x){extendrange(sapply(split(x,names(x)), quantile, tcall(quants5)[c(1,5)], na.rm=TRUE))}
#browser();return()
	mochaLatte(dat, xfld="Chain", yfld="Value", ffld="Factor", panel=panel.box, xlim=xlim, mar=c(0,3.8,0,0), oma=c(4,3,0.5,1), tcl=-0.3, las=1, cex.axis=cex.axis, cex.lab=cex.lab, xlab=linguaFranca(xlab,lang), ylab=linguaFranca(ylab,lang), xfac=xfac, fn.ylim=fn.ylim , outline=outline)
	gc(verbose=FALSE)
	invisible(dat)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~panelBoxes


## panelChains---------------------------2019-05-09
##  Plots cumulative fequency of 'nchains' by partitioning one trace.
##  Revised from 'plotTracePOP'
##  mcmc=data.frame e.g, 'currentMCMC$P' from object created by 'importMCMC'.
##  Very difficult to manipulate trellis plots (RH)
## -----------------------------------------------
panelChains = function (mcmc, nchains=3, pdisc=0.1, 
   axes=FALSE, same.limits=FALSE, between=list(x=axes,y=axes), div=1, span=1/4,
   log=FALSE, base=10, main=NULL, xlab=NULL, ylab=NULL, 
   cex.main=1.2, cex.lab=1, cex.strip=0.8, cex.axis=0.8, las=0, 
   tck=0.4, tick.number=5, lty.trace=1, lwd.trace=1, col.trace="grey", 
   lty.median=1, lwd.median=1, col.median="black", lty.quant=2, lwd.quant=1, 
   col.quant="black", plot=TRUE, probs=tcall(quants3), lang="e", ...)
{
	panel.chain <- function(x, y, ...) {
		dots = list(...)
		unpackList(dots)
		if (is.null(dots$xlim)) xlim = range(x,na.rm=TRUE)
		if (is.null(dots$ylim)) ylim = range(y,na.rm=TRUE)
		abline (h=0.5, lty=3, lwd=1, col="grey")
		chainlink = rep(1:nchains,ff)
		for (i in 1:nchains) {
			z = is.element(chainlink,i)
			lines(x[z], y[z], lty=rep(lty.trace,nchains)[i], lwd=2, col=rep(col.trace,nchains)[i])
			#lines(x[z], y[z], lty=1, lwd=2, col=c("red","green4","blue"))
		}
	}

	if (pdisc>0 && pdisc<1)
		mcmc = mcmc[(round(pdisc*nrow(mcmc))+1):nrow(mcmc),]  # get rid of the first 'pdisc' (e.g., 10%)
	relation <- if (same.limits) "same" else "free"
	if (is.null(dim(mcmc))) {
		mcmc.name <- rev(as.character(substitute(mcmc)))[1]
		mcmc <- matrix(mcmc, dimnames=list(NULL, mcmc.name))
	}
	mcmc <- if (log) 
		log(mcmc/div, base=base)
	else mcmc/div
	mcmc <- as.data.frame(mcmc)
	n <- nrow(mcmc)
	ff = rep(round(n/nchains),nchains-1)
	ff = c(ff,n-sum(ff))
	p <- ncol(mcmc)
	dat <- data.frame(Factor=ordered(rep(names(mcmc), each=n), 
		names(mcmc)), Draw=rep(1:n, p), Chain=rep(rep(1:nchains,ff),p), Value=as.vector(as.matrix(mcmc)))

	dat$Index = paste(dat$Factor,dat$Chain,sep="-")
	vList     = split(dat$Value,dat$Index)
	qList     = sapply(vList,function(x){
		xsort  = sort(x)
		xscal  = xsort - min(xsort)
		ycumu  = cumsum(xscal)/sum(xscal)
		out    = cbind(x=xsort,y=ycumu)
		return(out)
	}, simplify = FALSE )
	dat$CumFreq = dat$ValueSort = NA
	for (i in names(qList)) {
		z = is.element(dat$Index,i)
		dat$ValueSort[z] = qList[[i]][,"x"]
		dat$CumFreq[z]   = qList[[i]][,"y"]
	}
	#mar=c(0,3,0,0), oma=c(4,3,0.5,1)
	mochaLatte(dat,xfld="ValueSort",yfld="CumFreq",ffld="Factor", panel=panel.chain, ylim=c(0,1), mar=c(2,2,0,0), oma=c(4,4,0.5,1), tcl=-0.3, las=1, cex.axis=cex.axis, cex.lab=cex.lab, xlab=linguaFranca(xlab,lang), ylab=linguaFranca(ylab,lang))
	#mochaLatte(dat,xfld="ValueSort",yfld="CumFreq",ffld="Factor", panel=panel.chain, ylim=c(0,1), mar=c(2,0,0,0), oma=c(1.5,4.5,0.5,1), tcl=-0.3, las=1, cex.axis=cex.axis, cex.lab=cex.lab, xlab=linguaFranca(xlab,lang), ylab=linguaFranca(ylab,lang))
	invisible(dat)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~panelChains


## panelTraces-------------------------2019-05-09
##  Plots the sequential  trace of MCMC samples 
##  with running median and (0.05, 0.95) quantiles.
## ---------------------------------------------RH
panelTraces=function (mcmc, mpd=mcmc[1,], nchains=1, pdisc=0, 
   axes=FALSE, same.limits=FALSE, between=list(x=axes,y=axes), div=1,
   log=FALSE, base=10, main=NULL, xlab=NULL, ylab=NULL, 
   cex.main=1.2, cex.lab=1.2, cex.strip=0.9, cex.axis=0.9, las=0, 
   tck=0.4, tick.number=5, xfac=NULL, 
	lang="e", ...)
{
	panel.trace <- function(x, y, ...) {
		dots = list(...)
		unpackList(dots)
		if (is.null(dots$xlim)) xlim = range(x,na.rm=TRUE)
		if (is.null(dots$ylim)) ylim = range(y,na.rm=TRUE)
		lty.trace  = 1;  lwd.trace  = 1;  col.trace  = "grey"
		lty.quant  = 2;  lwd.quant  = 1;  col.quant  = "black"
		lty.median = 1;  lwd.median = 1.5;  col.median = "blue"
		lines(x, y, lty=lty.trace, lwd=lwd.trace, col=col.trace)
		if (any(is.finite(y)) && var(y) > 0) {
			lines(x, cquantile.vec(y, prob=tcall(quants3)[1]), lty=lty.quant,  lwd=lwd.quant,  col=col.quant)
			lines(x, cquantile.vec(y, prob=tcall(quants3)[2]), lty=lty.median, lwd=lwd.median, col=col.median)
			lines(x, cquantile.vec(y, prob=tcall(quants3)[3]), lty=lty.quant,  lwd=lwd.quant,  col=col.quant)
			points(x[1], mpd[getNpan()], pch=21, col="black", bg="red", cex=2)
		}
	}

	if (pdisc>0 && pdisc<1)
		mcmc = mcmc[(round(pdisc*nrow(mcmc))+1):nrow(mcmc),]  # get rid of the first 'pdisc' (e.g., 10%)
	if (is.null(dim(mcmc))) {
		mcmc.name <- rev(as.character(substitute(mcmc)))[1]
		mcmc <- matrix(mcmc, dimnames=list(NULL, mcmc.name))
	}
	mcmc <- if (log) 
		log(mcmc/div, base=base)
	else mcmc/div
	mcmc <- as.data.frame(mcmc)
	ylim = NULL
	if (same.limits)
		ylim = range(mcmc,mpd,na.rm=TRUE)
	n <- nrow(mcmc)
	ff = rep(round(n/nchains),nchains-1)
	ff = c(ff,n-sum(ff))
	p <- ncol(mcmc)
	dat <- data.frame(Factor=ordered(rep(names(mcmc), each=n), 
		names(mcmc)), Draw=rep(1:n, p), Chain=rep(rep(1:nchains,ff),p), Value=as.vector(as.matrix(mcmc)))
	#dat$Index = paste(dat$Factor,dat$Chain,sep="-")

#browser();return()
	mochaLatte(dat, xfld="Draw", yfld="Value", ffld="Factor", panel=panel.trace, xlim=c(0,nrow(mcmc)), ylim=ylim, mar=c(0,3,0,0), oma=c(4,4,0.5,1), tcl=-0.3, las=1, cex.axis=cex.axis, cex.lab=cex.lab, xlab=linguaFranca(xlab,lang), ylab=linguaFranca(ylab,lang) )
	invisible(dat)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~panelTraces


## plotACFs-----------------------------2018-07-09
##  Plot ACFs for the estimated parameters.
##  Control eps and png from PBScape.r in plt.mcmcGraphs
##----------------------------------------------RH
plotACFs =function(mcmc, lag.max=60, lang="e") #, ptypes=tcall(PBSawatea)$ptype, pngres=400)
{
	#if (!is.null(dev.list())) on.exit(expandGraph(mfrow=c(1,1)))
	acfs  = apply(mcmc, 2, function(x){acf(x,plot=FALSE)$acf})
	ylim  = range(acfs[round(acfs,5)>-1 & round(acfs,5)<1])
	idx   = apply(mcmc, 2, allEqual)
	mcmcP = mcmc[,!idx,drop=FALSE]
	rc    = .findSquare(ncol(mcmc[,!idx]))
	#for (p in ptypes) {
		#if (p=="eps") postscript("paramAcf.eps", width=8, height=8, horizontal=FALSE,  paper="special")
		#else if (p=="png") png("paramAcf.png", width=8, height=8, units="in", res=pngres)
		expandGraph(mfrow=rc, mar=c(1,1,0,0), oma=c(3,3,0.5,0.5))
		sapply(1:ncol(mcmcP), function(i){
			ii  = colnames(mcmcP)[i]
			mcP = mcmcP[,i]
			acf(mcP, lag.max=lag.max, ylim=ylim, xaxt="n", yaxt="n", xlab="", ylab="")
			axis(1,labels=ifelse(i > (ncol(mcmcP)-rc[2]),TRUE,FALSE),cex.axis=1.2,las=1)
			axis(2,labels=ifelse(par()$mfg[2]==1,TRUE,FALSE),cex.axis=1.2,las=1)
			addLabel(0.95,0.95,ii,cex=1.5,adj=c(1,1), col="blue")
			box(lwd=1.5)
		})
		mtext(linguaFranca("Lag",lang),side=1,outer=TRUE,line=1.5,cex=1.5)
		mtext(ifelse(lang=="f","FAC","ACF"),side=2,outer=TRUE,line=1.25,cex=1.5)
		#dev.off()
	#}
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotACFs


## plotAges-----------------------------2019-08-12
##  Plot the MPD model fits to age data
##  (commercial or survey) using the awkward 
##  scape function `plotCA'.
##------------------------------------------AME/RH
plotAges = function(obj, what="c", maxcol=5, sexlab=c("Females","Males"),
   ptypes=tcall(PBSawatea)$ptype, pngres=400, lang=c("e","f"), ...)
{
	seriesType = paste0("CA",what)
	seriesList = sort( unique( obj[[seriesType]][["Series"]]) )
	if (what=="c") seriesName = tcall(PBSawatea)$Cnames[tcall(PBSawatea)$CApos]
	else           seriesName = tcall(PBSawatea)$Snames[tcall(PBSawatea)$SApos]
	#seriesName.f = gsub("Historical","historique", gsub("Triennial","triennal", gsub("Synoptic","synoptique", seriesName)))
	#sexlab.f = gsub("Males",eval(parse(text=deparse("m\u{00E2}les"))),gsub("Females","femelles",sexlab)) ## already defined globally in run-master.Snw
	CA.yrs  = lapply(split(obj[[seriesType]][["Year"]], obj[[seriesType]][["Series"]]), unique)
	CA.nyrs = sapply(CA.yrs,length)
	CA.fit  = sapply(split(obj[[seriesType]][["Fit"]], obj[[seriesType]][["Series"]]), function(x){!all(x==0)})

	for ( i in 1:length(seriesList) )  {
		ii  = seriesList[i]    ## sometimes a survey age series is missing
		if (!CA.fit[ii]) next  ## don't plot age fits if they are not fitted
		yrs = CA.yrs[i]; nyrs = CA.nyrs[i]
		if (!is.null(list(...)$years)) {
			yrs = intersect(CA.yrs[[i]],list(...)$years)
			nyrs = length(yrs)
		}
		ncols   = min(maxcol,max(nyrs,1))
		nrows   = ceiling(nyrs/ncols)
		age.layout = rev(c(nrows,ncols)) # backwards in stupid lattice
		## page width & page height
		pwidth  = ifelse(ncols>2,8,ifelse(ncols>1,6,4))
		pheight = ifelse(nrows>2,8,ifelse(nrows>1,6,4))
		#pheight = nrows/ncols*pwidth  ## this does not work consistently
#browser();return()
		CA.sex = unique(obj[[seriesType]][["Sex"]])
		for(plot.sex in CA.sex) {
			j = grep(plot.sex,CA.sex)
			fout = fout.e = paste0(ifelse(what=="c","ageComm","ageSurv"), plot.sex,"Ser",ii) ## need the right series for plot
			for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
				changeLangOpts(L=l)
				fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
				## legend key:
				CA.key = list(text = list(lab= c(ifelse(l=="f", "obs", "Obs"), ifelse(l=="f", "pr\u{00E9}d", "Pred"))), lines=list(col=c("black",ifelse(plot.sex=="Male","blue","red")), cex = c(1,NA)), type=c("p","l"), x=ifelse(ncols>2,0.85,ifelse(ncols>1,0.8,0.75)), y=ifelse(nrows>2,-0.04,ifelse(nrows>1,-0.06,-0.12)), pch=c(20,NA), lwd=2, between=0.8)
				for (p in ptypes) {
					#pname = paste0(ifelse(what=="c","ageComm","ageSurv"), plot.sex,"Ser",ii) ## need the right series for plot
					set  = if (!is.null(list(...)$set)) list(...)$set else ""
					pnames = paste0(fout,set)
					if (p=="eps") postscript(paste0(pnames,".eps"), width=pwidth, height=pheight, horizontal=FALSE,  paper="special", onefile=FALSE)
					else if (p=="png") png(paste0(pnames,".png"), units="in", res=pngres, width=pwidth, height=pheight)
					plotCA( obj, what=what, ylab=linguaFranca("Proportion",l), xlab=linguaFranca("Age class",l), sex=plot.sex, layout=age.layout, key=CA.key, main=paste0(linguaFranca(seriesName[i],l), " - ", linguaFranca(sexlab[j],l)), pch=16, col.lines=ifelse(plot.sex=="Male","dodgerblue","red"), lwd.lines=2 , series=ii, ...)  ## need to pass the right series ii
					if (p %in% c("eps","png")) dev.off()
				} ## end of plot type loop
			}; eop()
		} ## end of plot.sex loop
	} ## end of seriesList loop
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotAges


#plotB2---------------------------------2011-08-31
# This is an alteration of Arni Magnussons "plotB" function to accommodate
# PJS request not to show biomass prior to fishery and survey indices period.
#----------------------------------------------AME
plotB2 <- function (model, what="d", series=NULL, years=NULL, axes=TRUE,
    div=1, legend="bottom", main="", xlab="", ylab="",
    cex.main=1.2, cex.legend=1, cex.lab=1, cex.axis=0.8,
    las=1, tck=c(1, what == "d")/2, tick.number=5, lty.grid=3,
    col.grid="white", pch=16, cex.points=0.8, col.points="black",
    lty.lines=1:3, lwd.lines=2, col.lines="black", ratio.bars=3,
    col.bars="grey", plot=TRUE, ...)
{
    panel.linebar <- function(x, y, bars, ...) {
        panel.abline(h=pretty(y, tick.number), lty=lty.grid,
            col=col.grid)
        panel.superpose(x, y, ...)
        panel.barchart(bars$Year, bars$Value, horizontal=FALSE,
            box.ratio=ratio.bars, col=col.bars)
    }
    panel.bar <- function(x, y, ...) {
        panel.abline(h=pretty(y, tick.number), lty=lty.grid,
            col=col.grid)
        panel.barchart(x, y, horizontal=FALSE, box.ratio=ratio.bars,
            col=col.bars)
    }
    if (class(model) != "scape")
        stop("The 'model' argument should be a scape object, not ",
            chartr(".", " ", class(model)), ".")
    what <- match.arg(what, c("d", "s", "l"))
    las <- as.numeric(las)
    x <- model$B
    x <- data.frame(Year=rep(x$Year, ncol(x) - 1), Series=rep(names(x)[-1],
        each=nrow(x)), Value=as.vector(as.matrix(x[, -1])))
    x$Value <- x$Value
    if (is.null(series))
        series <- unique(as.character(x$Series))
    if (is.null(years))
        years <- unique(x$Year)
    ok.series <- x$Series %in% series
    if (!any(ok.series))
        stop("Please check if the 'series' argument is correct.")
    ok.years <- x$Year %in% years
    if (!any(ok.years))
        stop("Please check if the 'years' argument is correct.")
    x <- x[ok.series & ok.years, ]


    Bframe <- x[x$Series %in% grep("B", series, value=TRUE),]
    Bframe$Series <- factor(Bframe$Series)

	## Find the first year where there are fishery CPUE data.
    cpueYear1 <- min( model$CPUE$Year[ !is.na( model$CPUE$Obs ) ] )

	## Set all SB and VB values to NA for years less than cpueYear1.
    Bframe$Value[ Bframe$Year < cpueYear1 ] <- NA

    Rframe <- x[x$Series == "R", ]
    Yframe <- x[x$Series == "Y", ]

    Bframe$Value <- Bframe$Value/div[1]
    Rframe$Value <- Rframe$Value/rep(div, length.out=2)[2]
    Yframe$Value <- Yframe$Value/div[1]

    #mess = c(
    #"require(grid, quietly=TRUE, warn.conflicts=FALSE)",
    #"require(lattice, quietly=TRUE, warn.conflicts=FALSE)"
    #)
    #eval(parse(text=mess))
    if (trellis.par.get()$background$col == "#909090") {
        for (d in dev.list()) dev.off()
        trellis.device(color=FALSE)
    }
    main <- rep(main, length.out=2)
    xlab <- rep(xlab, length.out=2)
    ylab <- rep(ylab, length.out=2)
    las <- rep(las, length.out=2)
    mymain <- list(label=main[1], cex=cex.main)
    myxlab <- list(label=xlab[1], cex=cex.lab)
    myylab <- list(label=ylab[1], cex=cex.lab)
    myrot <- switch(as.character(las[1]), "0"=list(x=list(rot=0),
        y=list(rot=90)), "1"=list(x=list(rot=0), y=list(rot=0)),
        "2"=list(x=list(rot=90), y=list(rot=0)), "3"=list(x=list(rot=90),
            y=list(rot=90)))
    myscales <- c(list(draw=axes, cex=cex.axis, tck=tck,
        tick.number=tick.number), myrot)
    lty.lines <- rep(lty.lines, length.out=nlevels(Bframe$Series))
    lwd.lines <- rep(lwd.lines, length.out=nlevels(Bframe$Series))
    col.lines <- rep(col.lines, length.out=nlevels(Bframe$Series))
    mykey <- list(space=legend, text=list(lab=levels(Bframe$Series),
        cex=cex.legend), lines=list(lty=lty.lines, lwd=lwd.lines,
        col=col.lines))
    if (what == "s") {
        graph <- xyplot(Rframe$Value ~ Bframe$Value[Bframe$Series ==
            "SB"], main=mymain, xlab=myxlab, ylab=myylab,
            scales=myscales, pch=pch, cex=cex.points, col=col.points,
            ...)
        graph$x.limits[1] <- 0
    }
    else if (what == "d" && nrow(Bframe) > 0) {
        graph <- xyplot(Value ~ Year, groups=Series, data=Bframe,
            panel=panel.linebar, type="l", bars=Yframe,
            main=mymain, xlab=myxlab, ylab=myylab, scales=myscales,
            key=mykey, lty=lty.lines, lwd=lwd.lines, col=col.lines,
            ...)
    }
    else {
        graph <- xyplot(Value ~ Year, data=Yframe, panel=panel.bar,
            main=mymain, xlab=myxlab, ylab=myylab, scales=myscales,
            ...)
    }
    graph$y.limits[1] <- 0
    if (plot) {
        print(graph)
        invisible(x)
    }
    else {
        invisible(graph)
    }
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotB2


## plotBars-----------------------------2018-07-30
## Plot barplots of specific year age proportions.
## ---------------------------------------------RH
plotBars = function(res, type="N", prop=TRUE, year=min(res[[type]][["Year"]]), 
    sex=c(2,1), # sex 2 =females (gfbio) 1 = males, 3 = unisex (awatea)
    age=NULL, fill=c("orange","cyan","green"), 
    eps=FALSE, png=FALSE, win=TRUE, pngres=400, lang="e", ...)
{
	if (!any(type==names(res))) stop("Choose another object in the list")
	nyear = length(year)
	SEX = c("Male","Female","Unisex")
	Sex = SEX[sex]; nsex=length(Sex)
	M1  = res$extra$parameters$M1
	M   = if (length(M1)==1) M1 else rev(M1)[sex]
	names(M) = Sex
	fill = rep(fill,nsex)[1:nsex]; names(fill) = Sex
	dat = res[[type]]
	dat = dat[is.element(dat$Year,year),]
	dat = dat[is.element(dat$Sex,Sex),]
	if (is.null(age)) age = sort(unique(dat$Age))
	AGE = sort(unique(dat$Age))
	nage = length(age); Nage=length(AGE)
	mat = array(0,dim=c(Nage,nsex,nyear), dimnames=list(age=AGE,sex=Sex,year=year))
	for (i in year) {
		ii = as.character(i)
		idat = dat[is.element(dat$Year,i),]
		for (s in Sex) {
			sdat = idat[is.element(idat$Sex,s),]
			x = sdat[,type]; names(x) = sdat[,"Age"]
			if (prop) x = x/sum(x)
			mat[names(x),s,ii] = x
			
		}
	}
	#ncol=floor(sqrt(nyear)); nrow=ceiling(nyear/ncol)
	#nrow=min(nyear,4); ncol=ceiling(nyear/nrow)
	rc =  .findSquare(nyear)
	nrow=rc[1]; ncol=rc[2]
	fout = list(...)$fout
	if (is.null(fout)) fout = "ageBars"
	fnam = paste0(fout, paste(Sex,collapse=""))
	figs = c(eps=eps,png=png,win=win)
#browser();return()
	for (k in names(figs)){
		if (!figs[k]) next
		if (k=="eps") postscript(paste0(fnam,".eps"), horizontal=FALSE, paper="special", width=6.5, height=2.5*nrow)
		else if (k=="png") png(paste0(fnam,".png"), units="in", res=pngres, width=6.75, height=2.5*nrow)
		else resetGraph()
		par(mfcol=c(nrow,ncol),mgp=c(2,0.5,0),las=1,xaxs="i",
			mar = if(nyear==1) c(4,6,1,1) else c(4,3,1,1),
			oma = if(nyear==1) c(0,0,0,0) else c(0,3,0,0))
		for (i in year) {
			ii = as.character(i); aa=as.character(age)
			xy = barplot(t(mat[aa,,ii]), beside=TRUE, space=c(0,0), xaxt=ifelse(nage>15,"n","s"),
				col=fill, xlab=linguaFranca("Age",lang), ylab="", cex.lab=1.5, cex.axis=1.25, cex=1.25)
			ylab = linguaFranca(ifelse(prop,"Proportions",type), lang)
			if (nyear==1) mtext(ylab,side=2,outer=FALSE,line=4,cex=2,las=0)
			else          mtext(ylab,side=2,outer=TRUE,line=0.5,cex=1.75,las=0)
			xpos  = apply(xy,2,mean)
			xdiff = diff(xpos)[1]
			for (s in Sex) {
				m = M[s]; f = fill[s]
				a = 0:max(age)
				b = exp(-m*a)
				za=is.element(a,age)
				ypos = (b[za]/max(b))* max(mat[,,ii])
				lines(xpos,ypos,lwd=3,col="black")
				#lines(xpos,ypos,lwd=1,col=f)
			}
			if (nage>15) {
				AGE = seq(5,200,5)
				agelab = intersect(AGE,age)
				agepos = xpos[is.element(age,agelab)]
				axis(1,tick=FALSE,at=agepos,labels=agelab,mgp=c(2,0.1,0),cex.axis=1.25)
			}
			addLabel(0.95, 0.95, linguaFranca(ii,lang), adj=c(1,1), cex=1.75, font=2)
			if (i==year[1])
				addLegend(0.975, 0.875, legend=paste0(linguaFranca(Sex,lang)," M = ",M), fill=fill, bty="n", yjust=1, xjust=1)
		}
		if (any(k==c("eps","png"))) dev.off()
	}
	invisible(list(dat=dat,mat=mat,xpos=xpos))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotBars


## plotBVBnorm--------------------------2018-07-12
## AME doing, tried in separate file, but then changed that to
##  lattice and wouldn't be good format for Arni's boxplots.
##  Based on plotVBcatch (tweaking some)
##  currentMCMC$B.  currentRes1 is local currentRes.
##  xLab - x position for label, etc.
## -----------------------------------------AME/RH
plotBVBnorm=function(mcmcObj,
   p = tcall(quants5),
   xyType="quantBox",
   lineType=c(3,2,1,2,3),
   refLines=NULL, xLim=NULL, yLim=NULL,
   userPrompt=FALSE, save=TRUE, xLeg=0.05, yLeg=0.2,  # legend in relative (0,1) space
   yaxis.by=0.05, tcl.val=-0.2,
   B.col="black", VB.col="black", ngear=1, lang="e", ...)
   # xLab - x position for label, etc.
{
	BVBlist = as.list(0:ngear); names(BVBlist)=c("Spawning Biomass",Cnames[1:ngear])
	BVBlist[[1]] = mcmcObj$B
	for (g in 1:ngear) {
		gfile = mcmcObj$VB[,grep(paste0("_",g),names(mcmcObj$VB))]
		names(gfile) = substring(names(gfile),1,4)
		BVBlist[[g+1]] = gfile
	}
	## Calculate medians to be plotted
	BVB0list = sapply(BVBlist,function(x){x/x[,1]},simplify=FALSE)  #B/B0 and VB/VB0 each chain
	BVB0med  = sapply(BVB0list,function(x){apply(x,2,median)},simplify=FALSE)  # median each year

	## Plot quantiles of biomass using the posterior densities.
	if ( is.null(yLim) )
		yLim = c(0,max(sapply(BVB0med,max)))
	yrs = sapply(BVB0med,function(x){as.numeric(names(x))},simplify=FALSE)
	if ( is.null(xLim) )
		xLim = range(sapply(yrs,range))
	VB.col = rep(VB.col,ngear)[1:ngear]

	plot(0, 0, xlim=xLim, ylim=yLim, type="n", xlab=linguaFranca("Year",lang), ylab=linguaFranca("Biomass relative to unfished equilibrium",lang), cex.lab=1.25)
	axis(1, at=intersect(seq(1900,3000,5),xLim[1]:xLim[2]), tcl=tcl.val, labels=FALSE)
	axis(2, at=seq(0, yLim[2], by=yaxis.by), tcl=tcl.val, labels=FALSE)
	points(yrs[[1]], BVB0med[[1]], type="p", col=B.col) 
	for (i in 1:ngear)
		points(yrs[[i+1]], BVB0med[[i+1]], type="l", col=VB.col[i])
	legtxt = bquote(italic(B[t])/italic(B)[0] ~ .(linguaFranca("    Spawning",lang)))
	for (i in 1:ngear)
		legtxt = c(legtxt, bquote(italic(V[t])/italic(V)[0] ~ .(linguaFranca("    Vulnerable - ",lang)) ~ .(linguaFranca(tolower(Cnames[i]),lang)) ) )
	addLegend(xLeg,yLeg,legend=as.expression(legtxt),bty="n", lty=c(0,rep(1,ngear)), pch=c(1, rep(NA,ngear)), col=c(B.col, VB.col))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotBVBnorm


##.plotChains.lattice-------------------2018-04-10
## Candidate for deprecation !!!
## Plots cumulative fequency of 'nchains' by partitioning one trace.
## Revised from 'plotTracePOP'
## mcmc=data.frame e.g, 'currentMCMC$P' from object created by 'importMCMC'.
## -----------------------------------------------
.plotChains.lattice=function (mcmc, nchains=3, pdisc=0.1, 
   axes=FALSE, same.limits=FALSE, between=list(x=axes,y=axes), div=1, span=1/4,
   log=FALSE, base=10, main=NULL, xlab=NULL, ylab=NULL, 
   cex.main=1.2, cex.lab=1, cex.strip=0.8, cex.axis=0.8, las=0, 
   tck=0.5, tick.number=5, lty.trace=1, lwd.trace=1, col.trace="grey", 
   lty.median=1, lwd.median=1, col.median="black", lty.quant=2, lwd.quant=1, 
   col.quant="black", plot=TRUE, probs=tcall(quants3), ...)  # AME probs
{
	panel.trace <- function(x, y, ...) {
		panel.xyplot(x, y, type="n")
		chainlink = rep(1:nchains,f)
		for (i in 1:nchains) {
			z = is.element(chainlink,i)
			panel.xyplot(x[z], y[z], type="l", lty=lty.trace, lwd=2, col=rep(col.trace,nchains)[i])
		}
		#panel.xyplot(x, y, type="l", lty=lty.trace, lwd=lwd.trace, col=col.trace)
	}

	if (pdisc>0 && pdisc<1)
		mcmc = mcmc[(round(pdisc*nrow(mcmc))+1):nrow(mcmc),]  # get rid of the first 'pdisc' (e.g., 10%)
	relation <- if (same.limits) "same" else "free"
	if (is.null(dim(mcmc))) {
		mcmc.name <- rev(as.character(substitute(mcmc)))[1]
		mcmc <- matrix(mcmc, dimnames=list(NULL, mcmc.name))
	}
	mcmc <- if (log) 
		log(mcmc/div, base=base)
	else mcmc/div
	mcmc <- as.data.frame(mcmc)
	n <- nrow(mcmc)
	f=rep(round(n/nchains),nchains-1); f=c(f,n-sum(f))
	p <- ncol(mcmc)
	dat <- data.frame(Factor=ordered(rep(names(mcmc), each=n), 
		names(mcmc)), Draw=rep(1:n, p), Chain=rep(rep(1:nchains,f),p), Value=as.vector(as.matrix(mcmc)))

	if (trellis.par.get()$background$col == "#909090") {
		for (d in dev.list()) dev.off()
		trellis.device(color=FALSE)
	}
	mymain <- list(label=main, cex=cex.main)
	myxlab <- list(label=xlab, cex=cex.lab)
	myylab <- list(label=ylab, cex=cex.lab)
	myrot  <- switch(as.character(las), `0`=0, `1`=0, `2`=0, `3`=90)     # AME changed '0'=90 to 0
	myscales <- list(x=list(draw=axes, relation=relation, cex=cex.axis, tck=tck, tick.number=tick.number, rot=myrot), 
		y=list(draw=axes, relation=relation, cex=cex.axis, tck=tck, tick.number=tick.number, rot=myrot))
	mystrip  <- list(cex=cex.strip)

	dat$Index = paste(dat$Factor,dat$Chain,sep="-")
	vList     = split(dat$Value,dat$Index)
	qList     = sapply(vList,function(x){
		xsort  = sort(x)
		xscal  = xsort - min(xsort)
		ycumu  = cumsum(xscal)/sum(xscal)
		out    = cbind(x=xsort,y=ycumu)
		return(out)
	}, simplify = FALSE )
	dat$CumFreq = dat$ValueSort=NA
	for (i in names(qList)) {
		z = is.element(dat$Index,i)
		dat$ValueSort[z] = qList[[i]][,"x"]
		dat$CumFreq[z]   = qList[[i]][,"y"]
	}
	## RH fiddling (to test how trellis works)
	layout.widths = trellis.par.get("layout.widths")
	layout.widths$right.padding = 0  ## default = 1 = changes right outer margin
	trellis.par.set("layout.widths",layout.widths)
	layout.heights = trellis.par.get("layout.heights")
	layout.heights$top.padding =0  ## default=1
	trellis.par.set("layout.heights",layout.heights)

	graph <- xyplot(CumFreq ~ ValueSort  | Factor, panel=panel.trace, 
		data=dat, as.table=TRUE, between=between, main=mymain, 
		xlab=myxlab, ylab=myylab, par.strip.text=mystrip, 
		scales=myscales, ylim=c(0,1), ...)
	if (plot) {
		print(graph)
		invisible(dat)
	}
	else {
		invisible(graph)
	}
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.plotChains.lattice


##plotCPUE------------------------------2019-05-07
## Plotting CPUE and fit with error bars
##  (copying plotIndexNotLattice).
## obj=currentRes$CPUE
##----------------------------------------------RH
plotCPUE <- function(obj, main="", save=NULL, bar=1.96, yLim=NULL,
   ptypes=tcall(PBSawatea)$ptype, pngres=400, lang=c("e","f"), ...)
{
	seriesList = sort( unique( obj$Series ) )  ## sort is risky if not always in same order
	nseries  = length(seriesList)
	surveyHeadName = c("CPUE")
	cvpro    = tcall(PBSawatea)$cvpro
	if (is.null(cvpro))
		cvpro = rep(0, Nsurv + nseries)
	unpackList(tcall(PBSawatea)[c("runNo","rwtNo")])

	pwidth=6.0;  pheight=switch(nseries,5,8,9)

	fout = fout.e = "CPUEser"
	for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
		changeLangOpts(L=l)
		fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
		for (p in ptypes) {
			#pname = "CPUEser"
			if (p=="eps") postscript(paste0(fout,".eps"), width=pwidth, height=pheight, horizontal=FALSE, paper="special", onefile=FALSE)
			else if (p=="png") png(paste0(fout,".png"), units="in", res=pngres, width=pwidth, height=pheight)
			par(mfrow=c(nseries,1),mai=c(0.75,0.75,0,0.1),omi=c(0,0,0.25,0),mgp=c(2,.75,0))
			# par(mai=c(0.25, 0.5, 0.3, 0.2)) # JAE changed  for each figure was for POP 0.45, 0.5, 0.3, 0.2
			# par(omi=c(0.45,0.1,0,0))  # Outer margins of whole thing, inch
			yrTicks=as.numeric( obj$Year)

			for ( i in 1:nseries ) {
				ii = Nsurv + i  # to index the CPUE cvpro
				idx <- seriesList[i]==obj$Series
				seriesVals=obj[idx,]
				# seriesvals$Obs=seriesvals$Obs   # /q[i] - set to 1 anyway
				seriesVals$Hi <- seriesVals$Obs * exp(bar * seriesVals$CV)
				seriesVals$Lo <- seriesVals$Obs/exp(bar * seriesVals$CV)
				# browser(); return()
				yearsnotNA=seriesVals[ !is.na(seriesVals$Obs), ]$Year
				# yearsPlot=min(yearsnotNA):max(yearsnotNA)
				xLim=range(yearsnotNA)
				if(i==1)
					xLimAll=xLim    # range to use in next plot
				xLimAll=range(xLim, xLimAll)    # range to use in next plot
				#if(is.null(yLim))     
				yLim=c(0, max(seriesVals$Hi, na.rm=TRUE))
				plotCI(seriesVals$Year, seriesVals$Obs, ui=seriesVals$Hi,
					li=seriesVals$Lo, xlim=xLim, ylim=yLim, xlab=linguaFranca("Year",l),
					ylab=paste0(linguaFranca("CPUE index",l), ": ", linguaFranca(seriesList[i],l)), 
					gap=0, pch=21, col="blue", bg="cyan", lwd=2)
				lines(seriesVals$Year, seriesVals$Fit, lwd=2)
				axis( side=1, at=yrTicks, tcl=-0.2, labels=FALSE )
#browser();return()
				if (Ncpue==0)
					addLabel(0.95, 0.95, linguaFranca("CPUE not used",l), adj=c(1,0), cex=0.8, col="slategrey")
				else if (is.numeric(cvpro[ii]) && round(cvpro[ii],5)!=0 && rwtNo>0)
					addLabel(0.95, 0.95, paste0(linguaFranca("+ CV process error = ",l), cvpro[ii]), adj=c(1,0), cex=0.8, col="slategrey")
				# mtext( side=3, line=0.25, cex=0.8, outer=FALSE, surveyHeadName[i]) #  outer=TRUE
				# if(i==3)  mtext( side=2, line=-0.5, cex=1, outer=TRUE,"Relative biomass")
				# if(i==5)  mtext(side=1, line=0, cex=1, outer=TRUE, "Year")
			}  # cex was 0.8 for POP
		if (p %in% c("eps","png")) dev.off()
		} ## end p (ptypes) loop
	}; eop()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotCPUE


#plotCI---------------------------------2018-04-06
# Lifted and modified from gplots::plotCI
#-----------------------------------------------RH
plotCI = function (x, y=NULL, ui, li, uiw=0.05, liw=uiw, clipNA=TRUE,
   gap=1, col=par("col"), barcol=col, lwd=par("lwd"), lty=par("lty"), ...)
{
	xobjnam = deparse(substitute(x))
	yobjnam = deparse(substitute(y))
	if (is.list(x)) {
		y <- x$y
		x <- x$x
	}
	if (is.null(y)) {
		if (is.null(x)) 
			stop("both x and y NULL")
		y <- as.numeric(x)
		x <- seq(along = x)
	}
	if (clipNA) {
		y  = clipVector(y, NA)
		i  = as.numeric(names(y))
		x  = x[i]
		ui = ui[i]
		li = li[i]
	}
	dots = list(...)
	need = c("xlim","ylim","xlab","ylab","pch","bg","cex.axis","cex.lab","cex.main")
	narg = as.list(rep(NA,length(need))); names(narg) = need
	for (a in need)
		eval(parse(text=paste0(a," = dots$",a)))
	if (is.null(xlim))      xlim      = range(x, na.rm = TRUE)
	if (is.null(ylim))      ylim      = range(c(y, ui, li), na.rm=TRUE)
	if (is.null(xlab))      xlab      = xobjnam
	if (is.null(ylab))      ylab      = yobjnam
	if (is.null(pch))       pch       = 21
	if (is.null(bg))        bg        = "white"
	if (is.null(cex.axis))  cex.axis  = 1.2
	if (is.null(cex.lab))   cex.lab   = 1.5
	if (is.null(cex.main))  cex.main  = 1.5
	for (a in need)
		narg[[a]] = get(a)
	marg = dots[setdiff(names(dots),need)]

	do.call(plot,args=c(list(x=x,y=y),narg,marg))
#if (pch==20) {browser();return()}
	#plot(x, y, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, type="n", cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main)
	arrows(x, li, x, pmax(y-gap,li), col=barcol, lwd=lwd, lty=lty, angle=90, length=liw, code=1)
	arrows(x, ui, x, pmin(y+gap,ui), col=barcol, lwd=lwd, lty=lty, angle=90, length=uiw, code=1)
	points(x, y, pch=pch, col=col, bg=bg, cex=1.2)
	invisible(list(x = x, y = y))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotCI


## plotDensPOP--------------------------2018-07-12
##  editing plotMCMC::plotDens function to have
##  less whitesapce, not repeat x axis labels, and make y axes
##  the same scales. Can't just do through options. For Recruits
##  and Biomass. See plotDensPOPpar.r for parameters.
##  Tried y axes the same scales, but 1973-1975 are so narrow that
##  it makes all the others really small: same.limits=TRUE,
##  ylim=c(0,0.0005).
##  Andrew Edwards. Edited lines indicated by AME. 19 October 2010
## -----------------------------------------AME/RH
plotDensPOP = function (mcmc, probs=tcall(quants3)[c(1,3)], points = FALSE, axes = TRUE, 
   same.limits = FALSE, between = list(x = axes, y = axes), 
   div = 1, log = FALSE, base = 10, main = NULL, xlab = NULL, 
   ylab = NULL, cex.main = 1.2, cex.lab = 1, cex.strip = 0.8, 
   cex.axis = 0.7, las = 0, tck = 0.5, tick.number = 5, lty.density = 1, 
   lwd.density = 3, col.density = "black", lty.median = 2, lwd.median = 1, 
   col.median = "darkgrey", lty.outer = 3, lwd.outer = 1, col.outer = "darkgrey", 
   pch = "|", cex.points = 1, col.points = "black", plot = TRUE,
   MPD.height = 0.04, mpd=mcmc[1,], lang="e", ...)     #MPD.height, how far up to put MPD
{
	panel.dens <- function(x, ...) { ## x here seems to a vector
		if (any(is.finite(x)) && var(x) > 0) ## for each panel
			panel.densityplot(x, lty = lty.density, lwd = lwd.density, 
				col.line = col.density, plot.points = points, 
				pch = pch, cex = cex.points, col = col.points, ...)
		else panel.densityplot(x, type = "n", ...)

		panel.abline(v = quantile(x, probs=probs), lty=lty.outer, lwd=lwd.outer, col=col.outer)
		panel.abline(v = median(x), lty=lty.median, lwd=lwd.median, col=col.median)
		# scan(); print(current.panel.limits()$ylim[2])  - max of y
		# print(graph$y.limits) is list of all panels
		panel.xyplot(mpd[panel.number()], current.panel.limits()$ylim[2]*MPD.height, pch=19, col="red") # AME, MPD. 0.04 of way up; assumes ..ylim[1]=0
		panel.xyplot(mpd[panel.number()], current.panel.limits()$ylim[2]*MPD.height, pch=1, col="black") #AME
		# scan(); print(summary(x))		# Yes, here x is just vector
	}
	relation <- if (same.limits) "same" else "free"
	if (is.null(dim(mcmc))) {
		mcmc.name <- rev(as.character(substitute(mcmc)))[1]
		mcmc <- matrix(mcmc, dimnames = list(NULL, mcmc.name))
	}
	mcmc <- if (log) 
		log(mcmc/div, base = base)
	else mcmc/div
	mcmc <- as.data.frame(mcmc)
	n <- nrow(mcmc)
	p <- ncol(mcmc)
	x <- data.frame(Factor = ordered(rep(names(mcmc), each = n), 
		names(mcmc)), Draw = rep(1:n, p), Value = as.vector(as.matrix(mcmc)))
	# scan(); print(summary(x))
	#mess = c(
	#"require(grid, quietly = TRUE, warn.conflicts = FALSE)",
	#"require(lattice, quietly = TRUE, warn.conflicts = FALSE)"
	#)
	#eval(parse(text=mess))
	if (trellis.par.get()$background$col == "#909090") {
		for (d in dev.list()) dev.off()
			trellis.device(color = FALSE)
	}
	mymain <- list(label = linguaFranca(main,lang), cex = cex.main)
	myxlab <- list(label = linguaFranca(xlab,lang), cex = cex.lab)
	myylab <- list(label = linguaFranca(ylab,lang), cex = cex.lab)
	myrot <- switch(as.character(las), `0` = 0, `1` = 0, `2` = 90, `3` = 90)
	myscales <- list(y = list(draw = FALSE, relation = "free"), 
		x = list(draw = axes, relation = "same", cex = cex.axis, tck = tck,  rot = myrot, alternating = TRUE))

	# AME: for y, relation = "same" -> relation = "free"
	# AME: for x, draw = axes -> draw = FALSE, but then no
	# marks, so back to axes (which =TRUE)
	# alternating = TRUE, relation="same"
	#   at=c(0, 50000, 100000, 150000, 200000))/1000
	#   took out tick.number = tick.number,  

	mystrip <- list(cex = cex.strip)
	graph <- densityplot(~Value | Factor, panel = panel.dens, 
		data = x, as.table = TRUE, between = between, main = mymain, 
		xlab = myxlab, ylab = myylab, par.strip.text = mystrip, scales = myscales, ...)
	if (!log) {
		if (is.list(graph$y.limits)) 
			graph$y.limits <- lapply(graph$y.limits, function(y) {
				y[1] <- 0
				return(y)
			})
		else graph$y.limits[1] <- 0
	}
	if (plot) {
		print(graph) #panel.height=10, panel.width=10)
		invisible(x)
	}
	else {
		invisible(graph)
	}
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotDensPOP
#mcmcObj=currentMCMC
#plotDensPOP(mcmcObj$B[,getYrIdx(names(mcmcObj$B))]/1000, xlab="Female spawning biomass, Bt (1000 t)", 
#	between=list(x=0.2, y=0.2), ylab="Density", lwd.density=2, #panel.height=list(x=rep(1,5),unit="inches"), #*****Needs resolving
#	same.limits=TRUE, lty.outer=2, mpd=mpd.B[getYrIdx(names(mcmcObj$B))]/1000) #, layout=c(4,5)) 


#plotDensPOPpars------------------------2010-10-26
# editing plotMCMC::plotDens for parameters,
#  to put MPDs on. AME. 26th Oct 2010.
#----------------------------------------------AME
plotDensPOPpars =
    function (mcmc, probs=tcall(quants3)[c(1,3)], points = FALSE, axes = TRUE, 
    same.limits = FALSE, between = list(x = axes, y = axes), 
    div = 1, log = FALSE, base = 10, main = NULL, xlab = NULL, 
    ylab = NULL, cex.main = 1.2, cex.lab = 1, cex.strip = 0.8, 
    cex.axis = 0.7, las = 0, tck = 0.5, tick.number = 5, lty.density = 1, 
    lwd.density = 3, col.density = "black", lty.median = 2, lwd.median = 1, 
    col.median = "darkgrey", lty.outer = 3, lwd.outer = 1, col.outer = "darkgrey", 
    pch = "|", cex.points = 1, col.points = "black", plot = TRUE,
    MPD.height = 0.04, mpd=mcmc[1,], lang="e", ...)    # MPD.height, how far up to put MPD
{
	panel.dens <- function(x, ...) {
		if (any(is.finite(x)) && var(x) > 0) 
			panel.densityplot(x, lty = lty.density, lwd = lwd.density, 
				col.line = col.density, plot.points = points, 
				pch = pch, cex = cex.points, col = col.points, ...)
		else panel.densityplot(x, type = "n", ...)

		panel.abline(v = quantile(x, probs = probs), lty = lty.outer, lwd = lwd.outer, col = col.outer)
		panel.abline(v = median(x), lty = lty.median, lwd = lwd.median, col = col.median)
		panel.xyplot(mpd[panel.number()], current.panel.limits()$ylim[2]*MPD.height,  pch=19, col="red") # AME, MPD. 0.04 of way up; assumes ..ylim[1]=0
		panel.xyplot(mpd[panel.number()], current.panel.limits()$ylim[2]*MPD.height, pch=1, col="black") #AME
	}
	relation <- if (same.limits) "same" else "free"
	if (is.null(dim(mcmc))) {
		mcmc.name <- rev(as.character(substitute(mcmc)))[1]
		mcmc <- matrix(mcmc, dimnames = list(NULL, mcmc.name))
	}
	mcmc <- if (log) 
		log(mcmc/div, base = base)
	else mcmc/div
	mcmc <- as.data.frame(mcmc)
	n <- nrow(mcmc)
	p <- ncol(mcmc)
	x <- data.frame(Factor = ordered(rep(names(mcmc), each = n), 
		names(mcmc)), Draw = rep(1:n, p), Value = as.vector(as.matrix(mcmc)))
	#mess = c(
	#"require(grid, quietly = TRUE, warn.conflicts = FALSE)",
	#"require(lattice, quietly = TRUE, warn.conflicts = FALSE)"
	#)
	#eval(parse(text=mess))
	if (trellis.par.get()$background$col == "#909090") {
		for (d in dev.list()) dev.off()
			trellis.device(color = FALSE)
	}
	mymain <- list(label = linguaFranca(main,lang), cex = cex.main)
	myxlab <- list(label = linguaFranca(xlab,lang), cex = cex.lab)
	myylab <- list(label = linguaFranca(ylab,lang), cex = cex.lab)
	myrot <- switch(as.character(las), `0` = 0, `1` = 0, `2` = 90, `3` = 90)
	myscales <- list(y = list(draw = FALSE, relation = "free"), 
		x = list(draw = axes, relation = relation, cex = cex.axis, 
			tck = tck, tick.number = tick.number, rot = myrot))
	mystrip <- list(cex = cex.strip)
	graph <- densityplot(~Value | Factor, panel = panel.dens, 
		data = x, as.table = TRUE, between = between, main = mymain, 
		xlab = myxlab, ylab = myylab, par.strip.text = mystrip, 
		scales = myscales, ...)
	if (!log) {
		if (is.list(graph$y.limits)) 
			graph$y.limits <- lapply(graph$y.limits, function(y) {
				y[1] <- 0
				return(y)
			})
		else graph$y.limits[1] <- 0
	}
	if (plot) {
		print(graph)
		invisible(x)
	}
	else {
		invisible(graph)
	}
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotDensPOPpars


## plotDensPOPparsPrior-----------------2018-07-12
## Adding the prior automatically.
## -----------------------------------------AME/RH
plotDensPOPparsPrior <-
    function (mcmc, probs=tcall(quants3)[c(1,3)], points=FALSE, axes=TRUE, 
    same.limits=FALSE, between=list(x=axes, y=axes), 
    div=1, log=FALSE, base=10, main=NULL, xlab=NULL, 
    ylab=NULL, cex.main=1.2, cex.lab=1, cex.strip=0.8, 
    cex.axis=0.7, las=0, tck=0.5, tick.number=5, lty.density=1, 
    lwd.density=3, col.density="black", lty.median=2, lwd.median=1, 
    col.median="darkgrey", lty.outer=3, lwd.outer=1, col.outer="darkgrey", 
    pch="|", cex.points=1, col.points="black", plot=TRUE,
    MPD.height=0.04, mpd=mcmc[1,], lang="e", ...)    # MPD.height, how far up to put MPD
{
	panel.dens <- function(x, ...) {
		if (any(is.finite(x)) && var(x) > 0) 
			panel.densityplot(x, lty=lty.density, lwd=lwd.density, 
				col.line=col.density, plot.points=points, 
				pch=pch, cex=cex.points, col=col.points, ...)
		else panel.densityplot(x, type="n", ...)

		panel.abline(v=quantile(x, probs=probs), lty=lty.outer, lwd=lwd.outer, col=col.outer)
		panel.abline(v=median(x), lty=lty.median, lwd=lwd.median, col=col.median)
		panel.xyplot(mpd[panel.number()], current.panel.limits()$ylim[2]*MPD.height, pch=19, col="red") # AME, MPD. 0.04 of way up; assumes ..ylim[1]=0
		panel.xyplot(mpd[panel.number()], current.panel.limits()$ylim[2]*MPD.height, pch=1, col="black") #AME
		# panel.curve(priorDistList[[panel.number()]], min=-1, current.panel.limits()$xlim[1], current.panel.limits()$xlim[2], col="blue")
		# panel.curve(priorDistList[[1]], col="blue")
		# panel.curve(priorDistList[[panel.number()]](x), from=max(priorBoundsList[[panel.number()]][1],
		#		  current.panel.limits()$xlim[1]), to=min(priorBoundsList[[panel.number()]][2], 
		#		  current.panel.limits()$xlim[2]), col="blue") # need the bounds, from max of lower bound and panel xlim[1], to min of upper bound and panel xlim[2]
		panel.curve(priorDistList[[panel.number()]]
			(x, priorInput[panel.number(), ] ),
			from = max(priorInput[panel.number(), 2] , current.panel.limits()$xlim[1]),
			to = min(priorInput[panel.number(), 3], current.panel.limits()$xlim[2]), col="blue")
	}
	relation <- if (same.limits) "same" else "free"
	if (is.null(dim(mcmc))) {
		mcmc.name <- rev(as.character(substitute(mcmc)))[1]
		mcmc <- matrix(mcmc, dimnames=list(NULL, mcmc.name))
	}
	mcmc <- if (log) 
		log(mcmc/div, base=base)
	else mcmc/div
	mcmc <- as.data.frame(mcmc)
	n <- nrow(mcmc)
	p <- ncol(mcmc)
	x <- data.frame(Factor=ordered(rep(names(mcmc), each=n), 
		names(mcmc)), Draw=rep(1:n, p), Value=as.vector(as.matrix(mcmc)))
	#mess = c(
	#"require(grid, quietly=TRUE, warn.conflicts=FALSE)",
	#"require(lattice, quietly=TRUE, warn.conflicts=FALSE)"
	#)
	#eval(parse(text=mess))
	if (trellis.par.get()$background$col == "#909090") {
		for (d in dev.list()) dev.off()
			trellis.device(color=FALSE)
	}
	mymain <- list(label = linguaFranca(main,lang), cex = cex.main)
	myxlab <- list(label = linguaFranca(xlab,lang), cex = cex.lab)
	myylab <- list(label = linguaFranca(ylab,lang), cex = cex.lab)
	myrot <- switch(as.character(las), `0`=0, `1`=0, `2`=90, `3`=90)
	myscales <- list(y=list(draw=FALSE, relation="free"), 
		x=list(draw=axes, relation=relation, cex=cex.axis, 
		tck=tck, tick.number=tick.number, rot=myrot))
	mystrip <- list(cex=cex.strip)
	graph <- densityplot(~Value | Factor, panel=panel.dens, 
		data=x, as.table=TRUE, between=between, main=mymain, 
		xlab=myxlab, ylab=myylab, par.strip.text=mystrip, scales=myscales, ...)
	if (!log) {
		if (is.list(graph$y.limits)) 
			graph$y.limits <- lapply(graph$y.limits, function(y) {
				y[1] <- 0
				return(y)
			})
		else graph$y.limits[1] <- 0
	}
	if (plot) {
		print(graph)
		invisible(x)
	}
	else {
		invisible(graph)
	}
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotDensPOPparsPrior


## plotIndexNotLattice------------------2019-05-07
## Taking some of plt.idx, but doing plot.Index NOT as lattice
## obj = currentRes
## plotCI is now custom function in PBSawatea (not gplots)
## -----------------------------------------AME/RH
plotIndexNotLattice <- function(obj, main="", save=NULL,
   bar=1.96, ssnames=paste("Ser",1:9,sep=""),
   ptypes=tcall(PBSawatea)$ptype, pngres=400, lang=c("e","f"), ...)
{
	cvcol="slategrey"
	objSurv        = obj$Survey
	objCPUE        = obj$CPUE
	seriesList     = sort( unique( objSurv$Series ) ) ## sort is risky if not always in same order
	nseries        = length(seriesList)
	surveyHeadName = if (!exists(".PBSmodEnv")) PBSawatea$Snames else tcall(PBSawatea)$Snames
	surveyHeadName = gsub("QC Sound", "QCS", surveyHeadName)
	#surveyHeadName.f = gsub("Historical","historique", gsub("Triennial","triennal", gsub("Synoptic","synoptique", surveyHeadName)))
	cvpro = tcall(PBSawatea)$cvpro
	#if (is.null(cvpro) || all(cvpro==FALSE)) cvpro=0 ## cvpro can no longer be logical
	if (is.null(cvpro))
		cvpro = rep(0, nseries)
	unpackList(tcall(PBSawatea)[c("runNo","rwtNo")])

	# (1) Plot the survey indices 
	yrTicks=min(objSurv$Year):max(objSurv$Year)
	rc = .findSquare(nseries)
	fout = fout.e = "survIndSer"
	for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
		changeLangOpts(L=l)
		fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
		for (p in ptypes) {
			if (p=="eps") postscript(paste0(fout,".eps"), width=6.5, height=8.5, horizontal=FALSE,  paper="special")
			else if (p=="png") png(paste0(fout,".png"), units="in", res=pngres, width=6.5, height=8.5)
			par(mfrow=c(rc[1],rc[2]), mar=c(2,2,1.25,0.5), oma=c(1.5,1.5,1,0.5), mgp=c(1.75,0.5,0))
			for ( i in 1:nseries ) {
				idx <- seriesList[i]==objSurv$Series
				seriesVals=objSurv[idx,]
				# seriesvals$Obs=seriesvals$Obs   # /q[i] - set to 1 anyway
				seriesVals$Hi <- seriesVals$Obs * exp(bar * seriesVals$CV)
				seriesVals$Lo <- seriesVals$Obs/exp(bar * seriesVals$CV)
				yearsnotNA=seriesVals[ !is.na(seriesVals$Obs), ]$Year
				# yearsPlot=min(yearsnotNA):max(yearsnotNA)
				xLim=range(yearsnotNA)
				if(i==1)
					xLimAll=xLim    # range to use in next plot
				xLimAll=range(xLim, xLimAll)    # range to use in next plot
				yLim=c(0, max(seriesVals$Hi, na.rm=TRUE))
				plotCI(seriesVals$Year, seriesVals$Obs, ui=seriesVals$Hi, li=seriesVals$Lo,
					xlim=xLim, ylim=yLim, xlab="", ylab="", gap=0, pch=19) # restrict years for plot, does error bars
				lines(seriesVals$Year, seriesVals$Fit, lwd=2)
				axis( side=1, at=yrTicks, tcl=-0.2, labels=FALSE )
				mtext( side=3, line=0.25, cex=ifelse(nseries>2,1.2,1.5), outer=FALSE, linguaFranca(surveyHeadName[i],l)) #  outer=TRUE
				if (is.numeric(cvpro[i]) && round(cvpro[i],5)!=0 && rwtNo>0)
					addLabel(0.95, 0.95, paste0(linguaFranca("+ CV process error = ",l), cvpro[i]), adj=c(1,1),cex=0.8,col=cvcol)
				if(i==nseries) {
					mtext(side=2, line=0, cex=1.5, outer=TRUE, text=linguaFranca("Relative biomass",l))
					mtext(side=1, line=0, cex=1.5, outer=TRUE, text=linguaFranca("Year",l))
				}
			}
			if (p %in% c("eps","png")) dev.off()
		}  # cex was 0.8 for POP
	}; eop()
#browser();return()

	# (2) And again, but with the same year axis for each. Think will be instructive to see
	ymaxsurvIndSer3=0           # For ymaxsurvIndSer3.eps
	xLimmaxsurvIndSer3=NA
	XLIM=numeric()
	fout = fout.e = "survIndSer2"
	for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
		changeLangOpts(L=l)
		fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
		for (p in ptypes) {
			if (p=="eps") postscript(paste0(fout,".eps"), width=6.5, height=8.5, horizontal=FALSE,  paper="special")
			else if (p=="png") png(paste0(fout,".png"), units="in", res=pngres, width=6.5, height=8.5)
			par(mfrow=c(nseries,1), mar=c(0,2,0,1), oma=c(3.2,3,0.5,0), mgp=c(1.75,0.5,0))
			for ( i in 1:length(seriesList) ) {
				idx <- seriesList[i]==objSurv$Series
				yrTicks=as.numeric( objSurv$Year[idx])        ## all years
				#YrTicks=intersect(seq(1900,2100,10),yrTicks) ## decadal ticks (doesn't work for short survey in odd years)
				XLIM = range(c(XLIM,seriesVals$Year[!is.na(seriesVals$Obs)]))
				YrTicks=pretty(XLIM)
				seriesVals=objSurv[idx,]
				# seriesvals$Obs=seriesvals$Obs   # /q[i] - set to 1 anyway
				seriesVals$Hi <- seriesVals$Obs * exp(bar * seriesVals$CV)
				seriesVals$Lo <- seriesVals$Obs/exp(bar * seriesVals$CV)
				yearsnotNA=seriesVals[ !is.na(seriesVals$Obs), ]$Year
				# yearsPlot=min(yearsnotNA):max(yearsnotNA)
				# xLim=range(yearsnotNA)
				yLim=c(0, max(seriesVals$Hi, na.rm=TRUE))
				# For axis for survIndSer3.eps:
				ymaxsurvIndSer3=max(ymaxsurvIndSer3, max(seriesVals$Obs, na.rm=TRUE)/ mean(seriesVals$Obs, na.rm=TRUE) )
				xLimmaxsurvIndSer3=range(xLimmaxsurvIndSer3, yearsnotNA, na.rm=TRUE)     # setting xLimmaxsurvIndSer3=NA above
				#gplots::plotCI(seriesVals$Year, seriesVals$Obs, ui=seriesVals$Hi,
				plotCI(seriesVals$Year, seriesVals$Obs, ui=seriesVals$Hi, li=seriesVals$Lo, xaxt="n", xlim=xLimAll, ylim=yLim, xlab="", ylab="", gap=0, pch=19)
				# restrict years for plot, does error bars
				lines(seriesVals$Year, seriesVals$Fit, lwd=2, col="blue")
				axis( side=1, at=yrTicks, tcl=-0.2, labels=FALSE)
				axis( side=1, at=YrTicks, tcl=-0.4, labels=ifelse(i==nseries,TRUE,FALSE),cex.axis=1.2)
				#mtext(side=4, line=1.5, cex=0.8, outer=FALSE, paste0(strwrap(surveyHeadName[i],10),collapse="\n"))
				mtext(side=2, line=1.75, cex=ifelse(nseries>2,1,1.5), outer=FALSE, linguaFranca(surveyHeadName[i],l), col="blue")
				if (is.numeric(cvpro[i]) && round(cvpro[i],5)!=0 && rwtNo>0)
					addLabel(0.95, 0.95, paste0(linguaFranca("+ CV process error = ",l), cvpro[i]), adj=c(1,1), cex=.8+(.05*(nseries-1)), col=cvcol)
				if(i==nseries) {
					#axis( side=1, at=yrTicks, tcl=-0.2, labels=FALSE)
					#axis( side=1, at=YrTicks, tcl=-0.4, labels=TRUE)
					mtext(side=1, line=2,   cex=1.5, outer=TRUE, text=linguaFranca("Year",l))
					mtext(side=2, line=1.2, cex=1.5, outer=TRUE, text=linguaFranca("Relative biomass",l))
				}
			} ## cex was 0.8 for POP
			if (p %in% c("eps","png")) dev.off()
		} ## end p (ptypes) loop
	}; eop()
#browser();return()

	# (3) And again, but all series on same plot, normalised to their means to see trends. Maybe add CPUE also.
	# Calculate max of normalised surveys and CPUEs
	objsurv=objSurv[!is.na(objSurv$Obs),]
	norsurv=sapply(split(objsurv$Obs,objsurv$Series),function(x){xx=x[!is.na(x)]; if (length(xx)==0) 0 else xx/mean(xx)},simplify=FALSE)
	maxsurv=max(sapply(norsurv,max))
	yrssurv=range(objsurv$Year); yrsspan=yrssurv[1]:yrssurv[2]
	objcpue=objCPUE[is.element(objCPUE$Year,yrsspan),]
	norcpue=sapply(split(objcpue$Obs,objcpue$Series),function(x){xx=x[!is.na(x)]; if (length(xx)==0) 0 else xx/mean(xx)},simplify=FALSE)
	maxcpue=max(sapply(norcpue,max))

	fout = fout.e = "survIndSer3"
	for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
		changeLangOpts(L=l)
		fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
		for (p in ptypes) {
			if (p=="eps") postscript(paste0(fout,".eps"), width=6.5, height=6.5, horizontal=FALSE,  paper="special")
			else if (p=="png") png(paste0(fout,".png"), units="in", res=pngres, width=6.5, height=6.5)
			par(mfrow=c(1,1), mar=c(3,3,0.5,1), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
			#postscript("survIndSer3.eps", height=6.0, width=6.0, horizontal=FALSE,  paper="special")   # height was 6 for POP
			yrTicks=yrsspan
			#YrTicks=intersect(seq(1900,2100,10),yrTicks) ## decadal ticks (doesn't work for short survey in odd years)
			YrTicks=pretty(xLim)
			NSL=1:length(seriesList)
			CLRS=c("black","blue","red","green4","orange","purple","navy","salmon")
			clrs=rep(CLRS,max(NSL))[NSL]
			# Set up plot
			plot(NA, xlim=yrssurv, ylim=c(0, max(maxsurv,maxcpue)), xlab=linguaFranca("Years",l), ylab=ifelse(l=="f", "Les indices d'enqu\u{00EA}te normalis\u{00E9}s par des moyens", "Survey indices normalised by means"), xaxt="n")
			abline(h=1, col="grey")
			axis(1,at=yrTicks,tck=-0.01,labels=FALSE)
			axis(1,at=YrTicks,tck=-0.02,labels=TRUE)
			for ( i in NSL ) {
				idx <- is.element(objsurv$Series,seriesList[i])
				seriesVals=objsurv[idx,]
				yearsnotNA=seriesVals[ !is.na(seriesVals$Obs), ]$Year
				points(yearsnotNA, seriesVals$Obs[ !is.na(seriesVals$Obs)] / mean(seriesVals$Obs, na.rm=TRUE), pch=i, col=clrs[i], type="o")
			}
			legtxt = if (exists(".PBSmodEnv")) tcall(PBSawatea)$Snames else PBSawatea$Snames
			legtxt = gsub("QC Sound","QCS",legtxt)
			# Now draw on CPUE series also:
			nseries = length(seriesList)  #  Need nseries from surveys for col
			seriesListCPUE <- sort( unique( objcpue$Series ) )
			ncpue=length(seriesListCPUE)
			cpue = as.logical(obj$extra$likelihoods$CPUE)
			# sort risky if not always in same order.
			lty=rep(1,length(NSL))
			if(cpue) {
				lty=c(lty,rep(2,ncpue))
				NSL=c(NSL, (max(NSL)+1):(max(NSL)+ncpue))
				clrs=c(clrs,rep(CLRS,ncpue)[1:ncpue])
				legtxt=c(legtxt, gsub("Series","CPUE",seriesListCPUE))
				for ( i in 1:ncpue ) {
					idx <- is.element(objcpue$Series,seriesListCPUE[i])
					seriesVals=objcpue[idx,]
					yearsnotNA=seriesVals[ !is.na(seriesVals$Obs), ]$Year
					if (length(yearsnotNA)>0)
						points(yearsnotNA, seriesVals$Obs[ !is.na(seriesVals$Obs)] / mean(seriesVals$Obs, 
							na.rm=TRUE), pch=i+nseries, col=clrs[i+nseries], type="o", lty=2)
				}
			}
			addLegend(0.075, 0.975, bty="n", col=clrs, pch=NSL, legend=linguaFranca(legtxt,l), cex=0.8)
			if (p %in% c("eps","png")) dev.off()
		} ## end p (ptypes) loop
	}; eop()

	# (4) And again, but big figures for ease of viewing. 
	for ( i in 1:length(seriesList) ) {
		idx <- seriesList[i]==objSurv$Series
		seriesVals=objSurv[idx,]
		# seriesvals$Obs=seriesvals$Obs   # /q[i] - set to 1 anyway
		seriesVals$Hi <- seriesVals$Obs * exp(bar * seriesVals$CV)
		seriesVals$Lo <- seriesVals$Obs/exp(bar * seriesVals$CV)
		yearsnotNA=seriesVals[ !is.na(seriesVals$Obs), ]$Year
		# yearsPlot=min(yearsnotNA):max(yearsnotNA)
		xLim=range(yearsnotNA)
		yLim=c(0, max(seriesVals$Hi, na.rm=TRUE))

		fout = fout.e = "survIndSer4-"
		for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
			changeLangOpts(L=l)
			fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
			for (p in ptypes) {
				if (p=="eps") postscript(paste0(fout, i, ".eps"), width=6.5, height=6.5, horizontal=FALSE,  paper="special")
				else if (p=="png") png(paste0(fout, i, ".png"), units="in", res=pngres, width=6.5, height=6.5)
				par(mfrow=c(1,1), mar=c(3,3.25,2,1), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
				plotCI(seriesVals$Year, seriesVals$Obs, ui=seriesVals$Hi, li=seriesVals$Lo, 
					xlim=xLim, ylim=yLim, xlab=linguaFranca("Year",l), ylab=linguaFranca("Relative biomass",l), gap=0, pch=19, cex.lab=1.5)
				lines(seriesVals$Year, seriesVals$Fit, lwd=2)
				axis( side=1, at=yrTicks, tcl=-0.2, labels=FALSE )
				mtext( side=3, line=0.25, cex=1.5, outer=FALSE, linguaFranca(surveyHeadName[i],l)) #  outer=TRUE
				if (is.numeric(cvpro[i]) && round(cvpro[i],5)!=0 && rwtNo>0)
					addLabel(0.95, 0.95, paste0(linguaFranca("+ CV process error = ",l), cvpro[i]), adj=c(1,0), cex=0.8, col=cvcol)
				if (p %in% c("eps","png")) dev.off()
			} ## end p (ptypes) loop
		}; eop()
	}
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotIndexNotLattice


## plotMeanAge--------------------------2018-04-16
##  Plot observed and expected mean ages from 
##  commercial and survey C@A data.
##----------------------------------------------RH
plotMeanAge =function(obj, useCA=TRUE, useSA=TRUE, CAnames, lang="e")
{
	## Here plot the mean age for catch and surveys (`MAfun` in `utilsFun.r`)
#browser();return()
	MAc = MAfun(obj$CAc)       ## catch mean age
	MAs = MAfun(obj$CAs)       ## surveys mean age
	nseries = 0
	MAp = character()
	if (useCA) {
		nseries = nseries + length(sort(unique(obj$CAc$Series)))
		MAp = c(MAp, "MAc") }
	if (useSA) {
		nseries = nseries + length(sort(unique(obj$CAs$Series)))
		MAp = c(MAp, "MAs") }
	
	if (nseries>0) {
		MA.pjs = list()
		par(mfrow=c(nseries,1), mar=c(1.75,3,2,1), oma=c(2.5,2.5,0,0), mgp=c(2,0.75,0))
		for (m in MAp) {
			MA   = get(m)
			last = regexpr("-",names(MA$MAobs))-1
			for ( i in 1:length(MA$J)) {   #1:length(unique(MAsSurvNum)) ) 
				ii  = MA$J[i]
				iii = substring(names(MA$MAobs),1,last)
				z   = is.element(iii,ii)
				for (k in setdiff(names(MA),"J"))
					MA.pjs[[paste0(m,ii)]][[k]] = MA[[k]][z]
#if (m=="MAs" && ii==2) {browser();return()}
				ylim = extendrange(c(MA$MAobs[z]+MA$CI[z],MA$MAobs[z]-MA$CI[z],MA$MAexp[z]),f=0.1)
				plot(MA$Yr[z],MA$MAobs[z],type="n",xlab="",ylab="", cex.axis=1.4, ylim=ylim, las=1)
				lines(MA$Yr[z], MA$MAexp[z], col="blue", lwd=2)
				points(MA$Yr[z], MA$MAexp[z], pch=22, col="blue", bg="cyan", cex=1.75, lwd=1)
				CLlo = MA$MAobs[z]-MA$CI[z]
				CLhi = MA$MAobs[z]+MA$CI[z]
				xCI  = as.vector(rbind(MA$Yr[z],MA$Yr[z],rep(NA,length(MA$Yr[z]))))
				yCI  = as.vector(rbind(CLlo,CLhi,rep(NA,length(CLlo))))
				lines(xCI,yCI,col="green4",lwd=2)
				points(MA$Yr[z],MA$MAobs[z], pch=21, col="green4", bg="green", cex=1.5)
				if (m=="MAc") {
					if (missing(CAnames)) CAnames=tcall(PBSawatea)$Cnames[MAc$J]
						mtext(linguaFranca(CAnames[i],lang), side=3, line=0.25, cex=1, outer=FALSE)
					} else {
						surveyHeadName = if (!exists("tcall")) ssnames[MAs$J] else tcall(PBSawatea)$Snames[MAs$J]
						mtext(linguaFranca(surveyHeadName[i],lang), side=3, line=0.25, cex=1, outer=FALSE)
					}
			}
			mtext(linguaFranca("Mean Age (years)",lang), side=2, line=0.5, cex=1.5, outer=TRUE)
			if (par()$mfg[1]==par()$mfg[3]) mtext(linguaFranca("Year",lang), side=1, line=1.0, cex=1.5, outer=TRUE)
		} ## end MAp
#browser();return()
		do.call("assign", args=list(x="MA.pjs", value=MA.pjs, envir=.GlobalEnv))
		dump("MA.pjs",file="MA.pjs.r")  ## for Paul
		save("MA.pjs",file="MA.pjs.rda")
	} ## end nseries
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotMeanAge


#plotBmcmcPOP [defunct]-----------------2011-08-31
# AME writing plotBmcmcPOP(), to plot spawning biomass and vulnerable biomass
#  from posterior, so do as boxplot, and also the catch on same graph. So
#  combining some of plt.quantBio and plotB2. Don't need lattice, just one
#  figure, no panels..     Vulnerable biomass has no posterior saved, which
#  must be why it's not been done before. Hmmm.... still worth seeing spawning
#  though?
# Taking what's needed from plt.quantBio, this basically works:
#        plt.quantBio( currentMCMC$B, xyType=rpType ), though does 2x3 plots
# obj should be the specficic MCMC posterior by year (so just a data.frame),
#  e.g. currentMCMC$B.  currentRes1 is local currentRes.
#----------------------------------------------AME

## plotRmcmcPOP-------------------------2018-07-12
## AME adding, plotting recruitment posteriors quantiles as one graph over time.
##  Already have the full posterior densities done.
##  Using plotBmcmcPOP as template, but will be simpler as no extra stuff. Prob
##   not simplifying down as much as could due to time constraints.
## Adding yLab and then using for exploitation plot also
## -----------------------------------------AME/RH
plotRmcmcPOP=function(obj, 
   p = tcall(quants5),
   xyType="quantBox",
   lineType=c(3,2,1,2,3),
   refLines=NULL, xLim=NULL, yLim=NULL,
   userPrompt=FALSE, save=TRUE, tcl.val=-0.2,
   yaxis.by=10000, yLab="Recruitment", lang="e", ...)
{
	# See plt.quantBio if want other xyTypes, as took out here:
	plt.qB <- function( obj, xyType="lines", new=TRUE, xLim, yLim,... ) 
		{
		if ( new )
			plot( xLim, yLim, type="n", xlab=linguaFranca("Year",lang), ylab=linguaFranca(yLab,lang) )

		yrs <- as.numeric(dimnames(obj)[[2]])

		# Quantile boxplots - assumes five quantiles.
		if ( xyType=="quantBox" )
		{
			delta <- 0.25 ## width of half-box
			# Draw the outer whiskers.
			segments( yrs,obj[1,], yrs,obj[5,], lty=1,col=1 )
			# Overlay the box.
			for ( i in 1:length(yrs) )
				#rect( yrs[i]-delta,obj[2,i], yrs[i]+delta, obj[4,i],... ) ## AME
				polygon(x=c(rep(yrs[i]-delta,2),rep(yrs[i]+delta,2)), y=obj[c(2,4,4,2),i], border="black", col="gainsboro", ...) ## RH (190620)
			# Add the median.
			segments( yrs-delta,obj[3,],yrs+delta,obj[3,],lty=1,col=1 )
		}
	}
	# Plot quantiles of biomass using the posterior densities.
	yrs1 = yrs2 = result1 = result2 = NULL

	# Calculate the quantiles of the reconstructed biomass.
	result1 <- apply( obj,2,quantile,probs=p )
	yrs1 <- as.numeric(dimnames(result1)[[2]])

	if ( is.null(yLim) )
		{
			yLim <- range(result1)
		}
		if ( is.null(xLim) )
		{
			xLim=range(yrs1)
		}
	plt.qB( result1, xLim=xLim, yLim=yLim, xyType=xyType )
	axis(1, at=intersect(seq(1900,3000,5), xLim[1]:xLim[2]), tcl=tcl.val, labels=FALSE)
	axis(2, at=seq(0, yLim[2], by=yaxis.by), tcl=tcl.val, labels=FALSE)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotRmcmcPOP


## plotSnail----------------------------2019-11-22
## Plot snail-trail plots (aka 'Kobe plots') for MCMC analysis.
##  AME: replacing "2010" with as.character(currYear - 1)
##  RH: added assYrs = years past with estimated Bcurr
##      from previous assessment(s), e.g., 5ABC QCS c(2011, 2017)
## -----------------------------------------AME/RH
plotSnail=function (BoverBmsy, UoverUmsy, p=c(0.05,0.95), xLim=NULL, yLim=NULL, 
	Lwd=1.5, ngear=1, currYear=2020, assYrs=NULL, outs=FALSE, Cnames, lang="e") ## outs = outliers
{
	if (missing(Cnames))
		Cnames  = tcall("PBSawatea")$Cnames        ## names of commercial gear
	## BU -- B = spawning biomass, U = harvest rate (or exploitation rate)
	BUlist = as.list(0:ngear); names(BUlist)=c("Spawning Biomass",Cnames[1:ngear])
	#BUlist[[1]] = BoverBmsy[,-length(BoverBmsy)]
	BUlist[[1]] = BoverBmsy[,-1]  ## conversation with PJS: we both agree that B2017/Bmsy should be paired with U2016/Umsy
	for (g in 1:ngear) {
		if (any(grepl("_",names(UoverUmsy)))) {
			gfile = UoverUmsy[,grep(paste0("_",g),names(UoverUmsy))]
			names(gfile) = substring(names(gfile),1,4)
		} else if (is.list(UoverUmsy)) {
			gfile = UoverUmsy[[g]]
		} else {
			gfile = UoverUmsy
		}
		BUlist[[g+1]] = gfile
	}
	# Calculate medians to be plotted
	BUmed    = sapply(BUlist,function(x){apply(x,2,median)},simplify=FALSE)  # median each year
	colPal   = colorRampPalette(c("grey95", "grey30"))
	colSlime = rep(c("grey","slategray2"),ngear)[1:ngear]
	colStart = rep(c("green","yellowgreen"),ngear)[1:ngear]
	colStop  = rep(c("cyan","thistle"),ngear)[1:ngear]  #,"cyan"
	colLim   = rep(c("blue2","purple"),ngear)[1:ngear]
	colAss   = rep(c("gold","orange"),ngear)[1:ngear]

	nB = length(BUmed[[1]])
	if (is.null(xLim))
		xLim=c(0, max(c(BUmed[[1]], rev(apply(BUlist[[1]],2,quantile,ifelse(outs,1,p[2])))[1], 1)))
	if (is.null(yLim))
		yLim=c(0, max(c(sapply(BUmed[(1:ngear)+1],max), rev(sapply(BUlist[(1:ngear)+1],function(x,p){apply(x,2,quantile,p)},p=ifelse(outs,1,p[2])))[1], 1)))
	P = list(p=p)
	if (outs)
		P = list (o=c(0,1), p=p)

	plot(0,0, xlim=xLim, ylim=yLim, type="n", 
		xlab = switch(lang, 'e'=expression(paste(italic(B[t])/italic(B)[MSY])),   'f'=expression(paste(italic(B[t])/italic(B)[RMS])) ),
		ylab = switch(lang, 'e'=expression(paste(italic(u[t-1])/italic(u)[MSY])), 'f'=expression(paste(italic(u[t-1])/italic(u)[RMS])) ),
		cex.lab=1.25, cex.axis=1.0, las=1)
	abline(h=1, col=c("grey20"), lwd=2, lty=3)
	abline(v=c(0.4,0.8), col=c("red","green4"), lwd=2, lty=2)
	for (i in ngear:1) {
		lines(BUmed[[1]], BUmed[[i+1]], col=colSlime[i], lwd=Lwd)
		points(BUmed[[1]], BUmed[[i+1]], type="p", pch=19, col=colPal(nB))
		points(BUmed[[1]][1], BUmed[[i+1]][1], pch=21, col=1, bg=colStart[i], cex=1.2)
		xend = rev(BUmed[[1]])[1]
		yend = rev(BUmed[[i+1]])[1]
		if (!is.null(assYrs))
			points(BUmed[[1]][as.character(assYrs)], BUmed[[i+1]][as.character(assYrs)], pch=21, col=1, bg=colAss[i], cex=1.2)
		for (j in 1:length(P)) {
			q = P[[j]]
			lty = ifelse(j==1 && outs, 3, 1)
			lwd = ifelse(j==1 && outs, 1, 2)
			## Plot horizontal (BtB0) quantile range
			xqlo = quantile(BUlist[[1]][,as.character(currYear)],q[1])
			xqhi = quantile(BUlist[[1]][,as.character(currYear)], q[2])
			segments(xqlo, yend, xqhi, yend, col=lucent(colLim[i],0.5), lty=lty, lwd=lwd)
			## Plot vertical (UtU0) quantile range
			yqlo = quantile(BUlist[[i+1]][, as.character(currYear-1)], q[1])
			yqhi = quantile(BUlist[[i+1]][, as.character(currYear-1)], q[2])
			segments(xend, yqlo, xend, yqhi, col=lucent(colLim[i],0.5), lty=lty, lwd=lwd)
		}
		points(xend, yend, pch=21, cex=1.2, col=colLim[i], bg=colStop[i])
	}
	if (ngear>1)  addLegend(0.95, 0.80, legend=linguaFranca(Cnames,lang), lty=1, lwd=Lwd, col=colSlime, seg.len=4, xjust=1, bty="n", cex=0.8)
	box()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotSnail


## plotTracePOP-------------------------2018-07-12
## Now adding running median, and taking off overall
## median and lowess line. Using cquantile from cumuplot.
## Trying to add in the MPD as a big circle for
## trace plots. 20th Oct 2010 (20/10/2010!)
## -----------------------------------------AME/RH
plotTracePOP = function (mcmc, axes = FALSE, same.limits = FALSE, between = list(x = axes, 
   y = axes), div = 1, span = 1/4, log = FALSE, base = 10, main = NULL, 
   xlab = NULL, ylab = NULL, cex.main = 1.2, cex.lab = 1, cex.strip = 0.8, 
   cex.axis = 0.8, las = 0, tck = 0.5, tick.number = 5, lty.trace = 1, 
   lwd.trace = 1, col.trace = "grey", lty.median = 1, lwd.median = 1, 
   col.median = "black", lty.quant = 2, lwd.quant = 1, col.quant = "black", 
   plot = TRUE, probs=tcall(quants3), mpd=mcmc[1,], lang="e", ...)  # AME probs
{
	panel.trace <- function(x, y, ...) {
		panel.xyplot(x, y, type = "l", lty = lty.trace, lwd = lwd.trace, col = col.trace)
		if (any(is.finite(y)) && var(y) > 0) {
			# print(x)  # gives 1 2 3 ... 1000 for each parameter/yr
			# panel.xyplot(range(x), rep(median(y), 2), type = "l", 
			#  lty = lty.median, lwd = lwd.median, col = col.median)
			panel.xyplot(x, cquantile.vec(y, prob=tcall(quants3)[1]),
				type = "l", lty = lty.quant, lwd = lwd.quant, col = col.quant)
			panel.xyplot(x, cquantile.vec(y, prob=tcall(quants3)[2]),
				type = "l", lty = lty.median, lwd = lwd.median, col = col.median)
			panel.xyplot(x, cquantile.vec(y, prob=tcall(quants3)[3]),
				type = "l", lty = lty.quant, lwd = lwd.quant, col = col.quant)
			panel.xyplot(x[1], mpd[panel.number()], pch=19, col="red") # AME
			panel.xyplot(x[1], mpd[panel.number()], pch=1, col="black") 
			# AME, based on plt.trace, assume x[1]=1
			# suppressWarnings(panel.loess(x, y, span = span,
			#  lty = lty.loess, lwd = lwd.loess, col =col.loess,...))
		}
	}
	relation <- if (same.limits) "same" else "free"
	if (is.null(dim(mcmc))) {
		mcmc.name <- rev(as.character(substitute(mcmc)))[1]
		mcmc <- matrix(mcmc, dimnames = list(NULL, mcmc.name))
	}
	mcmc <- if (log) 
		log(mcmc/div, base = base)
	else mcmc/div
	mcmc <- as.data.frame(mcmc)
	n <- nrow(mcmc)
	p <- ncol(mcmc)
	x <- data.frame(Factor = ordered(rep(names(mcmc), each = n), 
		names(mcmc)), Draw = rep(1:n, p), Value = as.vector(as.matrix(mcmc)))
	#mess = c(
	#"require(grid, quietly = TRUE, warn.conflicts = FALSE)",
	#"require(lattice, quietly = TRUE, warn.conflicts = FALSE)"
	#)
	#eval(parse(text=mess))
	if (trellis.par.get()$background$col == "#909090") {
		for (d in dev.list()) dev.off()
			trellis.device(color = FALSE)
	}
	mymain <- list(label = linguaFranca(main,lang), cex = cex.main)
	myxlab <- list(label = linguaFranca(xlab,lang), cex = cex.lab)
	myylab <- list(label = linguaFranca(ylab,lang), cex = cex.lab)
	myrot <- switch(as.character(las), `0` = 90, `1` = 0, `2` = 0, `3` = 90)
	myscales <- list(x = list(draw = FALSE), y = list(draw = axes, 
		relation = relation, cex = cex.axis, tck = tck, tick.number = tick.number, rot = myrot))
	mystrip <- list(cex = cex.strip)

	graph <- xyplot(Value ~ Draw | Factor, panel = panel.trace, 
		data = x, as.table = TRUE, between = between, main = mymain, 
		xlab = myxlab, ylab = myylab, par.strip.text = mystrip, scales = myscales, ...)
	if (plot) {
		print(graph)
		invisible(x)
	}
	else {
		invisible(graph)
	}
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotTracePOP


## plotTraj-----------------------------2019-12-06
## Show all median trajectories (base+sens) in one figure.
## ---------------------------------------------RH
plotTraj = function(sdat, index, traj="B", bdat=NULL, y0=FALSE, lab.stock, 
   startYear=1935, currYear=2020, sen.lab,
   col = c("black","green4","blue","red","purple","orange"), 
   lty = c(2:6), logR=FALSE, 
   png=FALSE, pngres=400, PIN=c(8,8), lang=c("e","f"), ...)
{
	opar = par(no.readonly=TRUE); on.exit(par(opar))
	Ntraj     = length(traj)
	trajYears = startYear:currYear
	Nyears    = length(trajYears)
	if (is.null(tcall(run.sens)))
		run.sens  = names(sdat$currentMCMC.sens)[index]
	else
		run.sens = tcall(run.sens)[index]
#browser();return()
	do.call("assign", args=list(x="run.sens", value=run.sens, envir=.GlobalEnv))
	Nruns     = length(run.sens)
	trajMat   = array(NA, dim=c(Nyears,Nruns,Ntraj), dimnames=list(year=trajYears, run=run.sens, traj=traj))
	#trajMat   = array(NA, dim=c(Nyears,Nruns,Ntraj+ifelse("BtB0"%in%traj,1,0)), dimnames=list(year=trajYears, run=run.sens, traj=if("BtB0"%in%traj) c(traj,"BmsyB0") else traj))
	Bmsy.q5   = list()  ## not always needed
	if (!is.null(bdat))
		col = setdiff(col,"black") ## save black for base case
	## ylabs needs to be a list if mixing character elements with expressions
	ylabs     = list("Spawning Biomass","Vulnerable Biomass","Recruitment","Exploitation Rate",expression(italic(B[t])/italic(B)[0]),"Unknown")
	names(ylabs) = c("B","VB","R","U","BtB0","NA")
	for (i in 1:Nruns) {
		ii = run.sens[i]
		for (j in 1:Ntraj) {
			jj = traj[j]
			if (jj=="BtB0") {
				jdat = sdat[["currentMCMC.sens"]][[ii]][["B"]]
				sdat[["currentMCMC.sens"]][[ii]][[jj]] = sweep(jdat,1,jdat[,1],"/")
				## Add Bmsy stuff (thanks PJS for the complication); and then he changed his mind...
				sdat[["currentMCMC.sens"]][[ii]][["BmsyB0"]] = sdat[["currentMSY.sens"]][[ii]][["B"]]/jdat[,1]
				Bmsy.q5[[ii]] = quantile(0.8*sdat[["currentMCMC.sens"]][[ii]][["BmsyB0"]],quants5)
			}
			jtmp = sdat[["currentMCMC.sens"]][[ii]][[jj]]
			if (jj %in% c("U","VB")) jtmp = splitGear(jtmp)[[1]]
			jval = apply(jtmp,2,median)
			if (jj=="R" && logR) jval = log10(jval)
			trajMat[substring(names(jval),1,4),ii,jj] = jval
		}
	}
	if (missing(lab.stock)) lab.stock = istock
	
	createFdir(lang)
	fout = fout.e = paste0(lab.stock,".sens.traj.",paste0(traj,collapse="+"))
	for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
		changeLangOpts(L=l)
		fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )

		## tcol and legtxt need to be reset inside language loop because they are altered when bdat is supplied
		tcol    = rep(col,Nruns)[1:Nruns]
		legtxt  = sen.lab[index]  ## defined in global data object 'stock'
#browser();return()
		tlty    = rep(lty,Nruns)[1:Nruns]

		if (png) png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
		par(mfrow=.findSquare(Ntraj), mar=c(3,3.5,0.5,0), oma=c(0,0,0,1), mgp=c(2,0.5,0))
		x = trajYears; xlim = range(x)
		for (j in 1:Ntraj) {
			jj   = traj[j]
			jmat = trajMat[,,jj]
			ylim = range(jmat,na.rm=TRUE)
			#if (jj %in% c("BtB0")) ylim[1]=0
			if (y0) ylim[1] = 0
			if (jj=="R" && logR) ylim[1]=1
			plot(0,0,xlim=xlim,ylim=ylim,type="n",xlab="",ylab="",log=ifelse(jj=="RRR","y",""))
			if (jj %in% c("BtB0")){
				#abline(h=seq(0.2,1,0.2),col="gainsboro")
				abline(h=c(0.2,0.4,1), lty=5, col="grey20") #col=c("salmon","darkorchid","navy"))
			}
			x  = as.numeric(rownames(jmat))
			for (k in 1:ncol(jmat)){
				kk = run.sens[k]
				y  = jmat[,k]
				lines(x,y, col=tcol[k], lwd=ifelse(i==1,2,2), lty=tlty[k]) ## no longer include base case in first position
			}
			if (!is.null(bdat)){  ## Assume for now that this is the Base Case (could be other runs)
				bline = bdat[[jj]]
#browser();return()
				lines(x=as.numeric(names(bline)), y=if (jj=="R" && logR) log10(bline) else bline, col="black", lty=1, lwd=3)
			}
			mtext(linguaFranca("Year",l), side=1, line=1.75, cex=ifelse(Ntraj==1,1.5,1))
			#ylab = ifelse(jj=="B","Spawning Biomass",ifelse(jj=="VB","Vulnerable Biomass",ifelse(jj=="R","Recruitment",ifelse(jj=="U","Exploitation Rate","Unknown"))))
			mtext(linguaFranca(ylabs[[jj]],l), side=2, line=1.8, cex=ifelse(Ntraj==1,1.5,1.2))
			if (j==1){
				legtxt = gsub("_"," ",legtxt)
				if (!is.null(bdat)){
					legtxt = c("Central run", legtxt[1:Nruns])
					tcol   = c("black",tcol)
					tlty    = c("solid", tlty)
				}
#browser():return()
				addLegend(ifelse(jj%in%c("R"),0.025,0.05), ifelse(jj%in%c("BtB0"),0.025,0.975), col=tcol, seg.len=5, legend=linguaFranca(legtxt,l), bty="o", box.col="grey", bg="white", xjust=ifelse(jj%in%c("R"),0,0), yjust=ifelse(jj%in%c("BtB0"),0,1), lwd=2, lty=tlty, ...)
			}
			if (Ntraj==1) {
				axis(1,at=intersect(seq(1900,2500,5),x),labels=FALSE,tcl=-0.2)
				axis(1,at=intersect(seq(1900,2500,10),x),labels=FALSE)
			}
			box()
		}
		if (png) dev.off()
	}; eop()
	do.call("assign", args=list(x="trajSens", value=trajMat, envir=.GlobalEnv))
	invisible()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotTraj


## plotVBcatch--------------------------2018-07-12
## AME adding, based on plotBmcmcPOP (tweaking some)
##  currentMCMC$B.  currentRes1 is local currentRes.
## -----------------------------------------AME/RH
plotVBcatch=function(obj, currentRes1=currentRes,
   p = tcall(quants5),
   xyType="quantBox",
   lineType=c(3,2,1,2,3),
   refLines=NULL, xLim=NULL, yLim=NULL,
   userPrompt=FALSE, save=TRUE, 
   xLab="Year",
   yLab="Catch and vulnerable biomass (t)",
   textLab=c("catch", "vulnerable"),
   yaxis.by=10000, tcl.val=-0.2,
   gear=1, lang="e", ...)
   # xLab - x position for label, etc.
{
	# See plt.quantBio if want other xyTypes, as took out here:
	plt.qB <- function( obj, xyType="lines", new=TRUE, xLim, yLim, ... ) 
	{
		if ( new )
			plot( xLim,yLim, type="n", xlab=linguaFranca(xLab,lang), ylab=linguaFranca(yLab,lang), ... )
		yrs <- as.numeric(dimnames(obj)[[2]])

		## Quantile boxplots - assumes five quantiles.
		if ( xyType=="quantBox" )
		{
			delta <- 0.25
			## Draw the outer whiskers.
			segments( yrs,obj[1,], yrs,obj[5,], lty=1,col=1 )
			## Overlay the box.
			for ( i in 1:length(yrs) )
				#rect( yrs[i]-delta,obj[2,i], yrs[i]+delta, obj[4,i], ... ) ## AME
				polygon(x=c(rep(yrs[i]-delta,2),rep(yrs[i]+delta,2)), y=obj[c(2,4,4,2),i], border="black", col="gainsboro", ...) ## RH (190620)
			## Add the median.
			segments( yrs-delta,obj[3,],yrs+delta,obj[3,],lty=1,col=1 )
		}
	}

	## Plot quantiles of biomass using the posterior densities.
	yrs1 <- yrs2 <- result1 <- result2 <- NULL

	## Calculate the quantiles of the reconstructed biomass.
	result1 <- apply( obj,2,quantile,probs=p )
	yrs1 <- as.numeric(dimnames(result1)[[2]])

	if ( is.null(yLim) )
		{
			yLim <- c(0, max(c(max(result1), max(currentRes1$B$VB)))) #range(result1)
		}
	if ( is.null(xLim) )
		{
			xLim=range(yrs1)
		}

	# xLegPos=xLeg*diff(xLim)+xLim[1]  ## position of xLeg
	# yLegPos=yLeg*diff(yLim)+yLim[1]

	plt.qB( result1,xLim=xLim,yLim=yLim, xyType=xyType, yaxt="n", ...)
	#points(currentRes1$B$Year, currentRes1$B$Y, type="h", lwd=3)	 # catch -- RH: won't work when Ngear > 1
	Cgears = currentRes1$B[,-1][,grep("Y",names(currentRes1$B[,-1])),drop=FALSE]
	Ctotal = apply(Cgears,1,sum)
	zpos = Cgears[,gear]>0
	points(currentRes1$B$Year[zpos], Cgears[,gear][zpos], type="h", lwd=3)	 ## gear-specific catch -- RH: use in case Ngear > 1
	# points(obj$B$Year, currentRes1$B$VB, type="p") ## was vuln biom MPD
	# text( xLab, yLab, linguaFranca(textLab,lang), pos=4, offset=0)   # Taking out as not really needed if give a decent caption
	axis(1, at=intersect(seq(1900,3000,5),xLim[1]:xLim[2]), tcl=tcl.val, labels=FALSE)
	axis(2, at=seq(0, yLim[2], by=yaxis.by), tcl=tcl.val, labels=FALSE)
	axis(2, at=pretty(yLim), labels=format(pretty(yLim), scientific=FALSE, big.mark=options()$big.mark))

	# legend(xLegPos, yLegPos, linguaFranca(c("Vulnerable","Spawning","Catch"),lang), bty="n")
	# points(xLegPos-2, yLegPos, type="p")
	# mtext( side=1, line=2, cex=1.0, linguaFranca("Year",lang))
	# mtext( side=2, line=2, cex=1.0, linguaFranca("Biomass",lang))
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotVBcatch
# RH -- following 3 lines for debugging only
#g=1; ngear=1
#test = currentMCMC$VB[,grep(paste0("_",g),names(currentMCMC$VB))]; names(test) = substring(names(test),1,4)
#plotVBcatch(test, currentRes, gear=g, yLab=ifelse(ngear==1,"Catch and vulnerable biomass (t)",Cnames[g]), yLim=c(0,max(sapply(test,quantile,tcall(quants5)[5]))),cex.lab=1.25)


#--------------------------------------------------------------------#
#                         Plotting Functions (plt)                   #
#--------------------------------------------------------------------#

# AME introducing for YMR. Bubble plots of data:
# plt.bubbles=function( obj,
# For now just do in Sweave as quicker.

## QUANTILE BOXES of AGE-FIT RESIDUALS by AGE, YEAR, COHORT ==========

## plt.ageResidsPOP---------------------2019-07-18
## AME changing for POP, just plotting age class resids
##  here, not qq-plot. Moving that to new function (so can do 4 on
##  page). See popScape.r for original.
## -----------------------------------------AME|RH
plt.ageResidsPOP <- function( obj, ages=NULL, main=NULL, lang="e")
{
	## Input is the output from stdRes.CA
	## par( oma=c(2,1,1,1), mar=c(2,2,2,1), mfrow=c(2,1) )
	## Subset to required ages.
	if (!is.null(ages))
		obj <- obj[ (obj$Age >= ages[1]) & (obj$Age <=ages[2]), ]
	else
		ages = range(obj$Age)
	if( max(diff(sort(unique(obj$Age)))) > 1) {
		allAges=min(obj$Age):max(obj$Age)
		nodataAges=allAges[ !(allAges %in% obj$Age)]
		xx=split(c(obj$stdRes, rep(NA, length(nodataAges))), c(obj$Age, nodataAges))
		#xpos <- boxplot( xx, whisklty=1, xlab="", ylab="", outline=FALSE )     #AME outline=FALSE removes outliers
		xpos = quantbox(xx, xaxt="n", xlab="", ylab="", pars=tcall(boxpars))
		xage = setdiff(allAges, nodataAges)
		if (length(xage)>20)
			xage = sort(unique(ceiling(xage/5)*5))
		axis(1, at=match(xage,as.numeric(xpos$names)), labels=xage, cex=1, mgp=c(2,0.5,0))
	} else {
		#xpos <- boxplot( split( obj$stdRes, obj$Age ), whisklty=1, xlab="", ylab="", outline=FALSE )     #AME outline=FALSE removes outliers
		xpos = quantbox(split(obj$stdRes,obj$Age), xaxt="n", xlab="", ylab="", pars=tcall(boxpars))
		xage = sort(unique(ceiling(obj$Age/5)*5))
		axis(1, at=match(xage,as.numeric(xpos$names)), labels=xage, cex=1, mgp=c(2,0.5,0))
	}
	abline( h=0, lty=2, col="red" )
	mtext( side=1, line=2, cex=0.8, text=linguaFranca("Age class",lang) )
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.ageResidsPOP


## plt.yearResidsPOP--------------------2019-07-18
## AME adding to plot age residuals by year. Is called for comm and survs.
##  fill.in=TRUE is to add the missing years for boxplot
##  ..POP does not do qq plot. See popScape.r for previous.
## -----------------------------------------AME|RH
plt.yearResidsPOP <- function(obj, ages=NULL, main=NULL, fill.in=TRUE, lang="e", ...)
{
	# Subset to required ages - still do as don't want age 1.
	if (!is.null(ages))
		obj <- obj[ (obj$Age >= ages[1]) & (obj$Age <=ages[2]), ]
	if(fill.in) {
		allYears = min(obj$Year):max(obj$Year)
		nodataYears = allYears[ !(allYears %in% obj$Year)]
		xx = split(c(obj$stdRes, rep(NA, length(nodataYears))), c(obj$Year, nodataYears))
		#xpos <- boxplot( xx, whisklty=1, xlab="", ylab="", outline=FALSE, ... )     #AME outline=FALSE removes outliers
		xpos = quantbox(xx, xaxt="n", xlab="", ylab="", pars=tcall(boxpars), ...)
		xyrs = setdiff(allYears, nodataYears)
		axis(1, at=match(xyrs,as.numeric(xpos$names)), labels=xyrs, cex=1, mgp=c(2,0.5,0))
	} else {
		#xpos <- boxplot( split( obj$stdRes, obj$Year ), whisklty=1, xlab="", ylab="", outline=FALSE, ... ) #AME outline=FALSE removes outliers
		xpos = quantbox( split(obj$stdRes,obj$Year), xlab="", ylab="", pars=tcall(boxpars), ...)
		xyrs = sort(unique(ceiling(obj$Year/5)*5))
		axis(1, at=match(xyrs,as.numeric(xpos$names)), labels=xyrs, cex=1, mgp=c(2,0.5,0))
	}
	abline( h=0, lty=2, col="red" )
	mtext( side=1, line=2, cex=0.8, text=linguaFranca("Year",lang) )
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.yearResidsPOP


## plt.cohortResids---------------------2019-07-18
## Plot age residuals by cohort.
## -----------------------------------------AME|RH
plt.cohortResids <- function( obj, ages=NULL, main=NULL, lang="e" )
{
	## Input is the CAc object from a Awatea res file. Ages to 59 as
	##  plus-age class will mess up year-of-birth calculation. Not automated.
	## par( oma=c(2,1,1,1), mar=c(2,2,2,1), mfrow=c(2,1) )
	## Subset to required ages - still do as don't want age 1 or 60 for cohorts
	if (!is.null(ages))
		obj <- obj[ (obj$Age >= ages[1]) & (obj$Age <=ages[2]), ]
	## obj$stdRes has residuals for each age, year and both sexes. Need
	##  to assign a year of birth for each age as an extra column, then
	##  presumably just do the boxplot split using that.
	obj$birthyr = obj$Year - obj$Age
	if( max(diff(sort(unique(obj$birthyr)))) > 1) {
		allYears = min(obj$birthyr):max(obj$birthyr)
		nodataYears = allYears[ !(allYears %in% obj$birthyr)]
		xx = split(c(obj$stdRes, rep(NA, length(nodataYears))), c(obj$birthyr, nodataYears))
		#xpos <- boxplot( xx, whisklty=1, xlab="", ylab="", outline=FALSE )     #AME outline=FALSE removes outliers
		xpos = quantbox(xx, xaxt="n", xlab="", ylab="", pars=tcall(boxpars))
		xyrs = setdiff(allYears, nodataYears)
		if (length(xyrs)>20)
			xyrs = sort(unique(ceiling(xyrs/5)*5))
		axis(1, at=match(xyrs,as.numeric(xpos$names)), labels=xyrs, cex=1, mgp=c(2,0.5,0))
	} else {
		#xpos=boxplot( split( obj$stdRes, obj$birthyr ), whisklty=1, xlab="", ylab="", outline=FALSE )     #AME outline=FALSE removes outliers
		xpos = quantbox(split(obj$stdRes,obj$birthyr), xaxt="n", xlab="", ylab="", pars=tcall(boxpars))
		xyrs = sort(unique(ceiling(obj$birthyr/5)*5))
		axis(1, at=match(xyrs,as.numeric(xpos$names)), labels=xyrs, cex=1, mgp=c(2,0.5,0))
	}
	abline( h=0, lty=2, col="red" )
	mtext( side=1, line=2, cex=0.8, text=linguaFranca("Year of birth",lang) )
	if ( !is.null(main) )
		mtext( side=3, line=-0.5, cex=1.0, outer=TRUE, text=linguaFranca(main,lang) )
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.cohortResids

##=================================Quantile boxes of age-fit residuals


## plt.ageResidsqqPOP-------------------2018-07-11
## Plotting qq plot for age class resids.
## -----------------------------------------AME/RH
plt.ageResidsqqPOP <- function( obj, ages=c(2,60),
   pct=c(5,25,50,75,95),  main=NULL, lang="e" )
{
  # Subset to required ages.
  obj <- obj[ (obj$Age >= ages[1]) & (obj$Age <=ages[2]), ]

  # Plot the q-normal plot of the standardised residuals 
  qqnorm( obj$stdRes, pch=20, col=lucent("black",0.5), xlab="", ylab="", main="" )
  abline( a=0, b=1 )
  abline( h=quantile(obj$stdRes,p=pct/100,na.rm=TRUE),lty=c(3,2,1,2,3) )
  mtext( side=1, line=2, cex=0.8, linguaFranca("Theoretical quantiles",lang) )

  #mtext( side=2, line=-1, cex=0.8, outer=TRUE, linguaFranca("Standardised Residuals",lang) )
  if ( !is.null(main) )
    mtext( side=3, line=-0.5, cex=1.0, outer=TRUE, text=linguaFranca(main,lang) )
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.ageResidsqqPOP


#plt.allTraces--------------------------2011-08-31
# Plot all MCMC traces
#----------------------------------------------AME
plt.allTraces <- function( obj, bioYrList=NULL, recYrList=NULL, save=TRUE )
{
  # Input a Awatea MCMC object.

  plt.trace <- function( obj )
  {
    # Input "obj" is a vector of MCMC samples.
    # Produces one panel trace plot.

    nSample <- length( obj )
    plot( c(1:nSample), obj, type="n", axes=FALSE, xlab="", ylab="" )
    points( c(1:nSample),obj, cex=0.25, pch=16, col="darkgray" )

    # Plot MPD point (1st row).
    points( 1,obj[1], cex=2.0, pch=16, col="green" )
    points( 1,obj[1], cex=2.0, pch=1 )

    lines( lowess( c(1:nSample),obj,f=1/4), lty=1, lwd=1 )
    abline( h=mean(obj), lty=2 )
    axis( side=2 )
    box()
  }

  # Find the active parameters.
  #iPars <- apply( obj$P,2,allEqual )
  #ACH: changed allEqual because you want to find the params that are not all the same (the ones that are estimated)
  iPars <- !apply( obj$P,2,allEqual )
  nPars <- sum( iPars )

  # (1) Plot biomass traces, every 5th year by default (getYrIdx).
  if ( is.null(bioYrList) )
    bioYrList <- getYrIdx( names(obj$B) )
  nYrs <- length( bioYrList )

  graphics( view="landscape" )
  par( oma=c(2,2,1,1), mar=c(1,2,1,1), mfrow=c(3,max(round(nYrs/3),4)) )
  for ( i in 1:nYrs )
  {
    biomass <- obj$B[ ,bioYrList[i] ]
    plt.trace( biomass )
    panLab( 0.5,0.95, cex=1.0, bioYrList[i] )

    mtext( side=1, line=0.5, cex=1.0, outer=TRUE, "Sample" )
    mtext( side=2, line=0.5, cex=1.0, outer=TRUE, "Biomass" )
  }
  if ( save )
    savePlot( paste("biomassTracePJS"), type="jpg" )


  # (2) Plot recruitment traces, every 5th year by default (getYrIdx).

  if ( is.null(recYrList) )
    recYrList <- getYrIdx( names(obj$R) )
  nYrs <- length( recYrList )

  graphics( view="landscape" )
  par( oma=c(2,2,1,1), mar=c(1,2,1,1), mfrow=c(3,max(round(nYrs/3),4)) )
  for ( i in 1:nYrs )
  {
    recruits <- obj$R[ ,recYrList[i] ]
    plt.trace( recruits )
    panLab( 0.5,0.95, cex=1.0, recYrList[i] )

    mtext( side=1, line=0.5, cex=1.0, outer=TRUE, "Sample" )
    mtext( side=2, line=0.5, cex=1.0, outer=TRUE, "Recruitment" )
  }
  if ( save )
    savePlot( paste("recruitTracePJS"), type="jpg" )

  # Plot the active parameters.
  graphics( view="landscape" )
  par( oma=c(2,2,1,1), mar=c(1,2,1,1), mfrow=c(3,min(round(nPars/3),3)) )

  # Find the active parameters.
  iPars <- apply( obj$P,2,allEqual )
  activePars <- obj$P[,!iPars]
  parNames <- names( activePars )
  nPars <- ncol( activePars )

  for ( i in 1:nPars )
  {
    parVec <- activePars[ ,i]
    plt.trace( parVec )
    panLab( 0.5,0.95, cex=1.0, parNames[i] )

    mtext( side=1, line=0.5, cex=1.0, outer=TRUE, "Sample" )
  }
  if ( save )
    savePlot( paste("parTracePJS"), type="jpg" )

  par( mfrow=c(1,1) )
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.allTraces


## plt.biomass--------------------------2018-07-11
## Small biomass figures
## transferred from Sweave `run-master.Snw'
## <<Btplot, results=hide, echo=FALSE>>=
## <<BtB0plot, results=hide, echo=FALSE>>=
## -----------------------------------------AME/RH
plt.biomass = function(years, Bt, xint=5, yint=2500,
   ptypes=tcall(PBSawatea)$ptype, pngres=400, pname="Bt",
   xlab="Year", ylab="Spawning biomass (t), Bt", lang=c("e","f"))
{
	#pname = gsub("\\.mpd","",as.character(substitute(Bt)))
	x = years; xlim = range(x); xsmall = intersect(seq(1900,2100,xint),x)
	if (is.null(dim(Bt))) y = matrix(Bt,ncol=1) else y=Bt
	ylim = c(0,max(y)); ysmall = seq(yint,ylim[2],yint)
	ngear=ncol(y)
	pchGear = seq(21,20+ngear,1)
	colGear = rep(c("black","blue"),ngear)[1:ngear]
	fout = fout.e = pname
	for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
		changeLangOpts(L=l)
		fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
		for (p in ptypes) {
			if (p=="eps") postscript(paste0(fout,".eps"), width=6, height=5, horizontal=FALSE,  paper="special")
			else if (p=="png") png(paste0(fout,".png"), units="in", res=pngres, width=6, height=5)
			par(mfrow=c(1,1), mar=c(3.2,3.2,1,1), oma=c(0,0,0,0), mgp=c(2,0.75,0))
			plot(0,0, xlim=xlim, ylim=ylim, type="n", xlab=linguaFranca(xlab,l), ylab=linguaFranca(ylab,l))
			sapply(1:ngear, function(g,x,y){
				points(x, y[,g], pch=pchGear[g], col=colGear[g], bg="white", cex=0.8) }, x=x,y=y)
			tcl.val = -0.2
			axis(1, at=xsmall, tcl=tcl.val, labels=FALSE)
			axis(2, at=ysmall, tcl=tcl.val, labels=FALSE)
			if (p %in% c("eps","png")) dev.off()
		} ## end p (ptypes) loop
	}; eop()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.biomass


## plt.bubbles--------------------------2018-07-11
## Bubble plots of observed and fitted ages
## transferred from Sweave `run-master.Snw'
## <<bubbleplots, results=hide, echo=FALSE>>=
## -----------------------------------------AME/RH
plt.bubbles = function(mpdObj, nsex=2,
   ptypes=tcall(PBSawatea)$ptype, pngres=400, redo.Graphs=TRUE, lang=c("e","f"))
{
	blow.bubbles = function(obj, cac, mod, sex, surnames=NULL) {
		series = sort(unique(obj[[cac]][["Series"]]))
		ages   = sort(unique(obj[[cac]][["Age"]]))
		CAlist = as.list(series); names(CAlist)=series
		for (i in names(CAlist)) {
			iCA = obj[[cac]][is.element(obj[[cac]][["Series"]],i),]
			CAlist[[i]] = matrix(iCA[[mod]][is.element(iCA[["Sex"]], sex)], nrow=length(ages))
			yrCA = sort(unique(iCA[["Year"]]))
			dimnames(CAlist[[i]])[[1]] = ages
			dimnames(CAlist[[i]])[[2]] = yrCA
		}
		if (!is.null(surnames) && length(surnames)==length(CAlist))
			names(CAlist) = surnames
		return(CAlist)
	}
	# These next two lines are't strictly necessary as they are now globally available
	#SAnames   = tcall("PBSawatea")$Snames[tcall("PBSawatea")$SApos] # names of surveys with ages
	#CAnames   = tcall("PBSawatea")$Cnames[tcall("PBSawatea")$CApos] # names of commercial gear with ages
	
#browser();return()
	CAcObsFem = blow.bubbles(mpdObj,"CAc","Obs",c("Female","Unisex"),surnames=CAnames)
	CAcFitFem = blow.bubbles(mpdObj,"CAc","Fit",c("Female","Unisex"),surnames=CAnames)
	CAsObsFem = blow.bubbles(mpdObj,"CAs","Obs",c("Female","Unisex"),surnames=SAnames)
	CAsFitFem = blow.bubbles(mpdObj,"CAs","Fit",c("Female","Unisex"),surnames=SAnames)
	if (nsex >1) {
		CAcObsMale = blow.bubbles(mpdObj,"CAc","Obs","Male",surnames=CAnames)
		CAcFitMale = blow.bubbles(mpdObj,"CAc","Fit","Male",surnames=CAnames)
		CAsObsMale = blow.bubbles(mpdObj,"CAs","Obs","Male",surnames=SAnames)
		CAsFitMale = blow.bubbles(mpdObj,"CAs","Fit","Male",surnames=SAnames)
	}

	CAplot = character()
	if (useCA) CAplot = c(CAplot,"CAc")
	if (useSA) CAplot = c(CAplot,"CAs")
	for (i in CAplot) {
		if (i=="CAc") { nr = length(CAnames); inames = CAnames }
		else          { nr = length(SAnames); inames = SAnames }
		for (j in c("Obs","Fit")) {
			if (nsex>1) kk = c("Fem","Male") else kk = "Fem"
			for (k in kk) {
				ijk = paste0(i,j,k); ijk.list = get(ijk)
				do.call("assign", args=list(x=ijk, value=ijk.list, envir=.GlobalEnv))
				if (redo.Graphs) {
					fout = fout.e = ijk
					for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
						changeLangOpts(L=l)
						fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
						for (p in ptypes) {
							if (p=="eps") postscript(paste0(fout,".eps"), width=6, height=ifelse(nr==1,6,8), horizontal=FALSE,  paper="special")
							else if (p=="png") png(paste0(fout,".png"), units="in", res=pngres, width=6, height=ifelse(nr==1,6,8))
							par(mfrow=c(nr,1), mar=c(2,3.5,2,1), oma=c(0,0,0,0), mgp=c(2,0.75,0))
							junk=sapply(1:nr,function(s,x,n){ # nr = no. rows = ngear or nsurv
								plotBubbles(x[[s]], dnam=TRUE, size=0.10, hide0=TRUE, main=linguaFranca(n[s],l), prettyaxis=TRUE, las=1)
								mtext(linguaFranca("Age",l), side=2, line=2, cex=1.2)
							}, x=ijk.list, n=inames)
							if (p %in% c("eps","png")) dev.off()
						} ## end p (ptypes) loop
					}; eop()
				}
			}
		}
	}
	if (useCA) {
		do.call("assign", args=list(x="residsCAcFem", value=sapply(CAcFitFem,function(x){prod(dim(x)+c(-1,0))}), envir=.GlobalEnv))
		if (nsex>1)
			do.call("assign", args=list(x="residsCAcMale", value=sapply(CAcFitMale,function(x){prod(dim(x)+c(-1,0))}), envir=.GlobalEnv))
	}
	if (useSA) {
		do.call("assign", args=list(x="residsCAsFem", value=sapply(CAsFitFem,function(x){prod(dim(x)+c(-1,0))}), envir=.GlobalEnv))
		if (nsex>1)
			do.call("assign", args=list(x="residsCAsMale", value=sapply(CAsFitMale,function(x){prod(dim(x)+c(-1,0))}), envir=.GlobalEnv))
	}
	#browser();return()
	invisible()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.bubbles


## plt.catch----------------------------2018-07-11
## Small catch figures
## transferred from Sweave `run-master.Snw'
## <<catch, results=hide, echo=FALSE>>=
## -----------------------------------------AME/RH
plt.catch = function(years, Ct, xint=5, yint=250,
   ptypes=tcall(PBSawatea)$ptype, pngres=400, lang=c("e","f"))
{
	x = years[-length(years)]; xlim = range(x)
	xsmall = intersect(seq(1900,2100,xint),x)
	if (is.null(dim(Ct))) y = matrix(Ct,ncol=1) else y=Ct
	ylim = c(0,max(y)); 
	ysmall = seq(yint,ylim[2],yint)
	ngear=ncol(y)
	pchGear=seq(21,20+ngear,1)
	colGear=rep(c("black","blue"),ngear)[1:ngear]
	fout = fout.e = "catch"
	for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
		changeLangOpts(L=l)
		fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
		for (p in ptypes) {
			if (p=="eps") postscript(paste0(fout,".eps"), width=6.5, height=4.5, horizontal=FALSE,  paper="special")
			else if (p=="png") png(paste0(fout,".png"), units="in", res=pngres, width=6.5, height=4.5)
			par(mfrow=c(1,1), mar=c(3.2,3.2,1,1), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
			#plot(x, y, type="h", xlim=xlim, ylim=ylim, xlab="Year", ylab="Catch (t)")
			xy = barplot(t(y), space=0.5, beside=FALSE, col=colGear, border="gainsboro", xlab=linguaFranca("Year",l), ylab=linguaFranca("Catch (t)",l), yaxs="i", names.arg=rep("",length(x)))
			lines(c(xy[1],rev(xy)[1]),c(0,0))
			axis(1, at=xy[match(xsmall,x)], tcl=-0.2, labels=xsmall, pos=0)
			axis(2, at=ysmall, tcl=-0.2, labels=FALSE)
			if (ngear>1) addLegend(0.05, 0.80, fill=colGear, linguaFranca(Cnames[1:ngear],l), yjust=0, bty="n")
#browser();return()
			if (p %in% c("eps","png")) dev.off()
		} ## end p (ptypes) loop
	}; eop()
	fout = fout.e = "catchSmall"
	for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
		changeLangOpts(L=l)
		fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
		for (p in ptypes) {
			if (p=="eps") postscript(paste0(fout,".eps"), width=6, height=3, horizontal=FALSE,  paper="special")
			else if (p=="png") png(paste0(fout,".png"), units="in", res=pngres, width=6, height=3)
			par(mfrow=c(1,1), mar=c(3.2,3.2,1,1), oma=c(0,0,0,0), mgp=c(2,0.75,0))
			plot(x, apply(y,1,sum), type="h", xlab=linguaFranca("Year",l), ylab=linguaFranca("Catch (t)",l))
			if (p %in% c("eps","png")) dev.off()
		} ## end p (ptypes) loop
	}; eop()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.catch


## plt.cpue-----------------------------2018-07-11
## Crude CPUE figure 
## transferred from Sweave `run-master.Snw'
## <<CPUEfig, results=hide, echo=FALSE>>= 
## -----------------------------------------AME/RH
plt.cpue = function(cpueObj, #xint=5, yint=2.5,
   ptypes=tcall(PBSawatea)$ptype, pngres=400, lang=c("e","f"))
{
	zobs = !is.na(cpueObj$Obs)
	xlim = range(cpueObj$Year[zobs])
	ylim = range(c(cpueObj$Obs[zobs],cpueObj$Fit[zobs]))
	fout = fout.e = "CPUEfit"
	for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
		changeLangOpts(L=l)
		fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
		for (p in ptypes) {
			if (p=="eps") postscript(paste0(fout,".eps"), width=6, height=6, horizontal=FALSE,  paper="special")
			else if (p=="png") png(paste0(fout,".png"), units="in", res=pngres, width=6, height=5)
			par(mfrow=c(1,1), mar=c(3.2,3.2,1,1), oma=c(0,0,0,0), mgp=c(2,0.75,0))
			plot(cpueObj$Year, cpueObj$Obs, xlim=xlim, ylim=ylim, type="n", xlab=linguaFranca("Year",l), ylab=linguaFranca("CPUE: Observed & Fitted",l))
			series = unique(cpueObj$Series)
			nseries = length(series)
			for (i in 1:nseries) {
				ii = series[i]; z = is.element(cpueObj$Series,ii)
				points(cpueObj$Year[z], cpueObj$Obs[z], pch=21, bg=i+1, cex=1.2)
				lines(cpueObj$Year[z], cpueObj$Fit[z], col=i+1, lwd=2)
			}
			legend("topright", bty="n", col=(1:nseries)+1, lwd=2, legend=linguaFranca(gsub("Series","CPUE",series),l), cex=0.8)
			if (p %in% c("eps","png")) dev.off()
		} ## end p (ptypes) loop
	}; eop()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.cpue


#plt.expRate----------------------------2011-08-31
# Input an object from "load.allResFiles".
# Plot exploitation rate against year.
#----------------------------------------------AME
plt.expRate <- function( obj, yLim=c(0,0.5), xLim=c(1954,2005) )
{
  nPanels <- length(obj)

  par( oma=c(2,2,1,1), mar=c(2,2,1,1), mfrow=c(2,round(nPanels/2)) )

  for ( i in 1:nPanels )
  {
    res <- obj[[i]]
    year <- res$B$Year

    plot( year, res$B$U, type="n",
      xlab="", xlim=xLim, ylab="", ylim=yLim )
    lines( year, res$B$U, lwd=2, lty=1 )

    yrTicks <- as.numeric( year[ year %% 10==0 ] )

    axis( side=1, at=yrTicks )
    axis( side=2 )
    axis( side=4, labels=FALSE )
    box()

    panLab( 0.5,0.95, names(obj)[i] )

    mfg <- par( "mfg" )
    if ( mfg[1]==mfg[3] & mfg[2]==mfg[4] | i==nPanels )
    {
      mtext( side=1, line=0.5, cex=1.0, outer=TRUE, "Year" )
      mtext( side=2, line=0.5, cex=1.0, outer=TRUE, "Exploitation Rate" )
    }
  }

  par( mfrow=c(1,1) )
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.expRate


## plt.idx------------------------------2019-09-26
## AME doing postscript for POP.
## RH adapting to five surveys for YMR.
## -----------------------------------------AME|RH
plt.idx <- function(obj, main="Residuals", save=NULL, ssnames=paste("Ser",1:9,sep=""),
   ptypes=tcall(PBSawatea)$ptype, pngres=400, lang=c("e","f"), ...)
{
	sType = substring(rev(as.character(substitute(obj)))[1],1,1)
	seriesList = sort( unique( obj$Series ) )
	nseries    = length(seriesList)
	surveyFigName = paste0(ifelse(sType=="S","surv",ifelse(sType=="C","cpue","unkn")),"Res",ssnames)
	surveyFigName = gsub(" ","",surveyFigName)
	surveyHeadName=if (!exists("tcall") || is.null(tcall(PBSawatea)[[paste0(sType,"names")]])) ssnames else tcall(PBSawatea)[[paste0(sType,"names")]]

	for ( i in 1:nseries )
	{
		idx <- seriesList[i]==obj$Series
		result <- stdRes.index( obj[idx,], label=paste(main,"Series",i) )
		#pname = surveyFigName[i]
		fout = fout.e = surveyFigName[i]
		for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
			changeLangOpts(L=l)
			fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
			for (p in ptypes) {
				if (p=="eps") postscript(paste0(fout,".eps"), height=6, width=5, horizontal=FALSE,  paper="special")
				else if (p=="png") png(paste0(fout,".png"), units="in", res=pngres, height=6, width=5)
				par(mfrow=c(nseries,1), mar=c(1.75,2,2,1), oma=c(2,2,0,0), mgp=c(2,0.75,0))
				plt.stdResids( result, xLim=range(result[!is.na(result$Obs), ]$Year), lang=l) # restrict years for plot
				mtext( side=3, line=0, cex=1.0, outer=TRUE, text=linguaFranca(surveyHeadName[i],l))
				if (p %in% c("eps","png")) dev.off()
			} ## end p (ptypes) loop
		}; eop()
	}
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.idx


## plt.initagedev-----------------------2018-07-11
## Initial age deviations figure
## transferred from Sweave `run-master.Snw'
## <<initagedevplot, results=hide, echo=FALSE>>=
## -----------------------------------------AME/RH
plt.initagedev = function(logInitAgeDev, 
   ptypes=tcall(PBSawatea)$ptype, pngres=400, lang=c("e","f"))
{
	fout = fout.e = "initAgeDev"
	for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
		changeLangOpts(L=l)
		fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
		for (p in ptypes) {
			if (p=="eps") postscript(paste0(fout,".eps"), width=6.5, height=4, horizontal=FALSE,  paper="special")
			else if (p=="png") png(paste0(fout,".png"), units="in", res=pngres, height=5, width=6)
			par(mfrow=c(1,1), mar=c(3.2,3.2,1,1), oma=c(0,0,0,0), mgp=c(2,0.75,0))
			plot(names(logInitAgeDev), logInitAgeDev, xlab=linguaFranca("Age",l), ylab=linguaFranca("Log initial age deviations",l))
			abline(h=0, col="grey")
			if (p %in% c("eps","png")) dev.off()
		} ## end p (ptypes) loop
	}; eop()
}
##----------------------------------plt.initagedev


#plt.numR-------------------------------2011-08-31
# Input an object from "load.allResFiles".
# Plot numbers at age at equilibrium.
# Plot recruitment age 1's.
#----------------------------------------------AME
plt.numR <- function( obj, minYr=NULL )
{
  nPanels <- length(obj)
  par( oma=c(2,2,1,1), mar=c(2,2,1,1), mfcol=c(2,2) )
  for ( i in 1:nPanels )
  {
    res <- obj[[i]]
    # Plot numbers at age for initial conditions.
    x <- res$N
    x <- x[ x$Year==min(x$Year), ]
    plot( x$Age, x$N, type="n", xlab="Age", ylab="", ylim=c(0,max(x$N)) )
    delta <- 0.4
    rect( x$Age-delta, 0, x$Age+delta, x$N, density=-1, col="gray" )
    axis( side=4, labels=FALSE )
    box()
    panLab( 0.5,0.95, names(obj)[i] )
    mfg <- par( "mfg" )
    if ( mfg[1]==1 & mfg[2]==1 )
      mtext( side=2, line=2.5, cex=1.0, "Initial numbers (000s)" )
    if ( is.null( minYr ) )
      minYr <- min(x$Year)
    # Age-1 recruits against year.
    x <- res$N
    x <- x[ x$Year<max(x$Year),]
    x <- x[ x$Age==min(x$Age), ]
    plot( x$Year, x$N, type="n", xlab="", ylab="", ylim=c(0,max(x$N)) )
    abline( v=seq(1955,2005,5 ), lty=3 )
    abline( v=seq(1950,2010,10), lty=2 )
    x <- x[ x$Year >= minYr, ]
    delta <- 0.35
    rect( x$Year-delta, 0, x$Year+delta, x$N, density=-1, col="gray" )
    mfg <- par( "mfg" )
    if ( mfg[1]==2 & mfg[2]==1 )
      mtext( side=2, line=2.5, cex=1.0, "Age 1's (000s)" )
    if ( mfg[1]==mfg[3] & mfg[2]==mfg[4] )
    {
      mtext( side=1, line=0, cex=1.0, outer=TRUE, "Year" )
    }
  }
  par( mfrow=c(1,1) )
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.numR


## plt.quantBio-------------------------2018-07-12
## Plot quantiles of reconstructed and projected biomass|recruits.
## ----------------------
## From popScapeRuns2.r -- AME now replaced yLim to force 0.
## This prints out tables (if run from command line), so be good to
##   use as template for decisions tables once we have MSY.
## --------------------------------------------AME
plt.quantBio <- function( obj, projObj=NULL, policy=NULL,
   p=tcall(quants5), xyType="lines", lineType=c(3,2,1,2,3),
   refLines=NULL, xLim=NULL, yLim=NULL,
   userPrompt=FALSE, save=TRUE, yaxis.lab="Spawning biomass", lang="e")
{
	opar = par(no.readonly=TRUE); on.exit(par(opar))
	## Subfunction
	plt.qB <- function( obj, xyType="lines", new=TRUE, xLim, yLim, line.col="black", ...)  #AME line.col="black", col is for filling rect
	{
		yrs <- as.numeric(dimnames(obj)[[2]])
		if ( new ) {
			plot( xLim,yLim, type="n", xlab="",ylab="", cex.axis=1.2)
			abline( v=yrs2[1]-0.5, lty=2, col=lucent("blue", 0.5), lwd=1 )
		}

		## Connect the quantiles with lines.
		if ( xyType=="lines" ) {
			for ( i in 1:nrow(obj) ) {
				# Plot reconstructed biomass.
				lines( yrs,obj[i,], lty=lineType[i],... )
			}
		}
		# ARK vertical line-dot plot.
		## Assumes that five quantiles requested, with median as one.
		if ( xyType=="lineDot" )
		{
			points( yrs,obj[2,], pch=3, cex=0.5,... )
			points( yrs,obj[3,], pch=1,... )
			points( yrs,obj[4,], pch=3, cex=0.5,... )
			segments( yrs,obj[1,], yrs,obj[5,], lty=1,... )
		}
		## Quantile boxplots - assumes five quantiles.
		if ( xyType=="quantBox" ) {
			delta <- 0.25
			# Draw the outer whiskers.
			segments( yrs,obj[1,], yrs,obj[5,], lty=1, col=line.col) #, ... )
			## AME col=1 removed, col=line.col ... so can have red for projs
			## Overlay the box.
			for ( i in 1:length(yrs) )
				rect( yrs[i]-delta,obj[2,i], yrs[i]+delta, obj[4,i], border=line.col, col="white")#AME border,col=NA (empty)
			## Add the median.
			segments( yrs-delta,obj[3,],yrs+delta,obj[3,],lty=1,col="black") #line.col )	 #AME black
		}
	} ## end subfunction `plt.qB'

	## Plot quantiles of biomass using the posterior densities.
	## If proj!=NULL then add the projections for all policies.
	## If policy!=NULL than plot only the specified policy.

	## Plotting ranges for reconstruction (1) and projection (2).
	nCol <- min( length(policy),3 )
	if ( !is.null(policy) )
		nRow <- length( policy )/nCol
	nRow <- 3
	nCol <- 2
	yrs1 = yrs2 = result1 = result2 = NULL

	par( oma=c(1.75,2.5,0.5,0), mar=c(2,1,1.5,1), mfrow=c(nRow,nCol), userPrompt, mgp=c(2,0.5,0) ) # AME mfcol -> mfrow to fill left to right

	# Calculate the quantiles of the reconstructed biomass.
	result1 <- apply( obj,2,quantile,probs=p )
	yrs1 <- as.numeric(dimnames(result1)[[2]])
	if ( is.null(yLim) )
		# yLim <- range( result1 )
		yLim=c(0, max( result1 ))

	## Reconstructed biomass.
	if ( is.null(projObj) ) {
		plt.qB( result1,xLim=range(yrs1),yLim=yLim, xyType=xyType )
		mtext( side=1, line=2, cex=1.0, text=linguaFranca("Year",lang) )
		mtext( side=2, line=2, cex=1.0, text=linguaFranca("Biomass",lang) )
	}

	## Reconstructed biomass and projections.
	if ( !is.null(projObj) ) {
		## Get the policies to be plotted.
		if ( is.null(policy) )
			policyList <- names( projObj )
		else
			policyList <- policy

		## Loop over the policies.
		result2 <- as.list(policyList)
		names( result2 ) <- policyList

		iPage <- 1
		nPolicies <- length(policyList)
		for ( j in 1:nPolicies ) {
			## Calculate quantiles of the projected biomass for policy.
			pol <- policyList[j]
			result2[[j]] <- apply( projObj[[pol]],2,quantile,probs=p )
			# cat( "\n\nQuantiles of projection for policy=",
			# policyList[j],"\n" )
			# print( result2[[j]] )
			yrs2 <- as.numeric( dimnames(result2[[j]])[[2]] )
			if ( is.null(xLim) )
				xLim <- range( yrs1,yrs2 )
			# yLim <- range( result1,result2[[j]] )	# AME to get same axes
			# yLim <- c(0,yLim[2])

			## Plot the quantiles of the biomass.
			if ( xyType=="quantBox" )
				plt.qB( result1, xyType=xyType, new=TRUE, xLim,yLim, col=NA, line.col="black", med.col="red" )
			else
				plt.qB( result1, xyType=xyType, new=TRUE, xLim,yLim, col="black" )
			if ( !is.null(refLines) )
				abline( v=refLines,lty=4 )

			## Plot the quantiles of the projected biomass.
			if ( xyType=="quantBox" )
				plt.qB( result2[[j]], xyType=xyType, new=FALSE, xLim,yLim, line.col="red")	# AME: col fills in box, I want to change line cols 
			else
				plt.qB( result2[[j]], xyType=xyType, new=FALSE, xLim,yLim, col="red" )
			#for ( i in 1:nrow(result2[[j]]) )
			#	lines( yrs2,result2[[j]][i,],lty=lineType[i],lwd=2,col=2 )
			#abline( v=yrs2[1]-0.5, lty=2, col=lucent("blue", 0.5) )

			mfg <- par( "mfg" )
			if ( mfg[1]==1 & mfg[2]==1 ) {
				mtext( side=1, line=0.5, cex=1.2, outer=TRUE, text=linguaFranca("Year",lang) )   ## AME line=0 changed
				mtext( side=2, line=1.0, cex=1.2, outer=TRUE, text=linguaFranca(yaxis.lab,lang)) ## " and Spawning
			}
			mtext( side=3, line=0.25, cex=0.9, text=paste0( linguaFranca("Catch strategy: ",lang), policyList[j]) )

			## Last panel on page or last policy.
			if ( mfg[1]==mfg[3] & mfg[2]==mfg[4] | j==nPolicies ) {
				if ( save )
					savePlot( paste0( "policyProj", iPage), type="png" )
				iPage <- iPage + 1
				if ( j < nPolicies ) {
					do.call("windows", list(record=TRUE))
					par( oma=c(2,2,1,1), mar=c(2,2,1.5,1), mfcol=c(nRow,nCol), userPrompt )
				}
			}
		}
	}
	val <- list( recon=result1, proj=result2 )
	return(val)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.quantBio


## plt.quantBioBB0----------------------2018-07-12
## This is for B/B0 for each run, with just one projection. Doing one
##  plot here instead of multiple, so taking some from plotBVBnorm
##  which did one for each run. Don't think used for biomasses.
##  Now using for single recruitment projection plot
## -----------------------------------------AME/RH
plt.quantBioBB0 <- function( obj, projObj=NULL, policy=NULL,
   p = tcall(quants5),
   xyType="lines",
   lineType=c(3,2,1,2,3),
   refLines=NULL, xLim=NULL, yLim=NULL,
   userPrompt=FALSE, save=TRUE, main="", cex.main="",
   tcl.val=-0.2, xaxis.by=1, yaxis.by=10000,
   xaxis.lab="Year", yaxis.lab= "Spawning biomass", lang="e" )
{
	plt.qB <- function( obj, xyType="lines", new=TRUE, xLim, yLim, line.col="black", med.col="black", ... )   #AME line.col="black", col is for filling rect med.col is for median
	{
		if ( new )
			plot( xLim,yLim, type="n", xlab=linguaFranca(xaxis.lab,lang), ylab=linguaFranca(yaxis.lab,lang) )
		yrs <- as.numeric(dimnames(obj)[[2]])
		## Connect the quantiles with lines.
		if ( xyType=="lines" )
		{
			for ( i in 1:nrow(obj) )
			{
				# Plot reconstructed biomass.
				lines( yrs,obj[i,], lty=lineType[i],... )
			}
		}
		## ARK vertical line-dot plot.
		## Assumes that five quantiles requested, with median as one.
		if ( xyType=="lineDot" )
		{
			points( yrs,obj[2,], pch=3, cex=0.5,... )
			points( yrs,obj[3,], pch=1,... )
			points( yrs,obj[4,], pch=3, cex=0.5,... )
			segments( yrs,obj[1,], yrs,obj[5,], lty=1,... )
		}
		## Quantile boxplots - assumes five quantiles.
		if ( xyType=="quantBox" )
		{
			delta <- 0.25
			# Draw the outer whiskers.
			segments( yrs,obj[1,], yrs,obj[5,], lty=1, col=line.col) #, ... )
			## AME col=1 removed, col=line.col ... so can have red for projs
			## Overlay the box.
			for ( i in 1:length(yrs) )
				rect( yrs[i]-delta,obj[2,i], yrs[i]+delta, obj[4,i], border=line.col, col="white") ## AME border,col=NA (empty)
			## Add the median.
			# segments( yrs-delta,obj[3,],yrs+delta,obj[3,],lty=1,col=line.col )   #AME black
			segments( yrs-delta, obj[3,], yrs+delta, obj[3,], lty=1, col=med.col )   #AME black
		}
	}
	## Plot quantiles of biomass using the posterior densities.
	## If proj!=NULL then add the projections for all policies.
	## If policy!=NULL than plot only the specified policy.
	## Just one plot now. Deleted some.
	yrs1 <- yrs2 <- result1 <- result2 <- NULL

	## Calculate the quantiles of the reconstructed biomass.
	result1 <- apply( obj,2,quantile,probs=p )
	yrs1 <- as.numeric(dimnames(result1)[[2]])
	if ( is.null(yLim) )
		yLim <- range( result1 )

	## Reconstructed biomass.
	if ( is.null(projObj) )
	{
		plt.qB( result1, xLim=range(yrs1), yLim=yLim, xyType=xyType )
		# mtext( side=1, line=2, cex=1.0, "Year" )
		# mtext( side=2, line=2, cex=1.0, "Biomass" )
	}
	## Reconstructed biomass and projections.
	if ( !is.null(projObj) )
	{
		## Get the policies to be plotted.
		if ( is.null(policy) )
			policyList <- names( projObj )
		else
			policyList <- policy

		## Loop over the policies.
		result2 <- as.list(policyList)
		names( result2 ) <- policyList
		iPage <- 1
		nPolicies <- length(policyList)
		for ( j in 1:nPolicies )
		{
			## Calculate quantiles of the projected biomass for policy.
			pol <- policyList[j]
			result2[[j]] <- apply( projObj[[pol]],2,quantile,probs=p )
			# cat( "\n\nQuantiles of projection for policy=",
			# policyList[j],"\n" )
			# print( result2[[j]] )
			yrs2 <- as.numeric( dimnames(result2[[j]])[[2]] )
			if ( is.null(xLim) )
				xLim <- range( yrs1,yrs2 )
			# yLim <- range( result1,result2[[j]] )  # AME to get same axes
			# yLim <- c(0,yLim[2])

			## Plot the quantiles of the biomass.
			if ( xyType=="quantBox" )
				plt.qB( result1, xyType=xyType, new=TRUE, xLim,yLim, col="black", line.col="black", med.col="red" )
			else
				plt.qB( result1, xyType=xyType, new=TRUE, xLim,yLim, col="red") # line.col="red")
			if ( !is.null(refLines) )
				abline( v=refLines,lty=4 )

			## Plot the quantiles of the projected biomass.
			if ( xyType=="quantBox" )
				plt.qB( result2[[j]], xyType=xyType, new=FALSE, xLim,yLim, line.col="red", med.col="black")  # AME: col fills in box, I want to change line cols 
			else
				plt.qB( result2[[j]], xyType=xyType, new=FALSE, xLim,yLim, col="red" )
			#for ( i in 1:nrow(result2[[j]]) )
			#  lines( yrs2,result2[[j]][i,],lty=lineType[i],lwd=2,col=2 )

			abline( v=yrs2[1]-0.5, lty=2 )
			axis(1, at=seq(xLim[1], xLim[2], by=xaxis.by), tcl=tcl.val, labels=FALSE)
			axis(2, at=seq(0, yLim[2], by=yaxis.by), tcl=tcl.val, labels=FALSE)

			# mfg <- par( "mfg" )
			# if ( mfg[1]==1 & mfg[2]==1 )
			# {
			#  mtext( side=1, line=0.8, cex=1.0, outer=TRUE, "Year" )   #AME line=0 changed
			#  mtext( side=2, line=0.8, cex=1.0, outer=TRUE,
			#			  "Spawning biomass" ) # "   and Spawning
			#}
			#mtext( side=3, line=0.25, cex=0.8,
			#  paste( "Policy:",policyList[j]) )

			# Last panel on page or last policy.
			#if ( mfg[1]==mfg[3] & mfg[2]==mfg[4] | j==nPolicies )
			#{
			#  if ( save )
			#		savePlot( paste( "policyProj",iPage,sep=""),type="png" )
			#  iPage <- iPage + 1
			#
			#  if ( j < nPolicies )
			#  {
			#		windows()
			#		par( oma=c(2,2,1,1), mar=c(2,2,1.5,1),
			#			   mfcol=c(nRow,nCol), userPrompt )
			#  }
			# }
		}
	}
	# par( mfrow=c(1,1) )
	val <- list( recon=result1, proj=result2 )
	return(val)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.quantBioBB0


## plt.recdev---------------------------2018-07-11
## Log recruitment deviations figure 
## transferred from Sweave `run-master.Snw'
## <<recdevplot, results=hide, echo=FALSE>>=
##------------------------------------------AME/RH
plt.recdev = function(logRecDev, xint=5, #yint=0.1,
   ptypes=tcall(PBSawatea)$ptype, pngres=400, lang=c("e","f"))
{
	x = as.numeric(names(logRecDev)); xlim = range(x); xsmall = intersect(seq(1900,2100,xint),x)
	fout = fout.e = "recDev"
	for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
		changeLangOpts(L=l)
		fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
		for (p in ptypes) {
			if (p=="eps") postscript(paste0(fout,".eps"), width=6.5, height=4, horizontal=FALSE,  paper="special")
			else if (p=="png") png(paste0(fout,".png"), units="in", res=pngres, width=6, height=4)
			par(mfrow=c(1,1), mar=c(3,3.5,1,1), oma=c(0,0,0,0), mgp=c(2,0.5,0))
			plot(x, logRecDev, type="n", xlab=linguaFranca("Year",l), ylab=bquote(.(linguaFranca("Log recruitment deviations, ",l)) ~ italic(epsilon[t])), las=1)
			abline(h=0, lwd=1, lty=3, col="grey")
			lines(x, logRecDev, col="gainsboro")
			zpos = logRecDev > 0
			points(x[zpos], logRecDev[zpos],   pch=21, cex=0.8, col="blue", bg="cyan")
			points(x[!zpos], logRecDev[!zpos], pch=21, cex=0.8, col="red", bg="pink")
			tcl.val = -0.2
			axis(1, at=xsmall, tcl=tcl.val, labels=FALSE)
			if (p %in% c("eps","png")) dev.off()
		} ## end p (ptypes) loop
	}; eop()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.recdev


## plt.recdevacf------------------------2018-07-11
## Auto-correlation function of the log recruitment deviations
## transferred from Sweave `run-master.Snw'
## <<recdevacf, results=hide, echo=FALSE>>=
## -----------------------------------------AME/RH
plt.recdevacf = function(logRecDev, muC, logvC, A, years, yr1, 
   ptypes=tcall(PBSawatea)$ptype, pngres=400, redo.Graphs=TRUE, lang=c("e","f"))
{
	ageHalfFemSelComm = min(round(muC - sqrt(exp(logvC) * log(2) ))) # take the min when Ngears>1
	do.call("assign", args=list(x="ageHalfFemSelComm", value=ageHalfFemSelComm, envir=.GlobalEnv))
	# Take selectivity equation and find a which satisifies s_{ag1} = 0.5. 
	# Working out saved in POP12 folder.
	# This may be 1 year off, given I've now fixed the recruitment indexing issue.
	yearsForACF = max((yr1 - A + ageHalfFemSelComm), years[1]) : (as.numeric(max(names(logRecDev))) - ageHalfFemSelComm)
	do.call("assign", args=list(x="yearsForACF", value=yearsForACF, envir=.GlobalEnv))
	# max() to start no earlier than first year of model
	logRecDevForACF = logRecDev[as.character(yearsForACF)]
	if (redo.Graphs) {
		fout = fout.e = "recDevAcf"
		for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
			changeLangOpts(L=l)
			fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
			for (p in ptypes) {
				if (p=="eps") postscript(paste0(fout,".eps"), width=6.5, height=4, horizontal=FALSE,  paper="special")
				else if (p=="png") png(paste0(fout,".png"), width=6, height=5, units="in", res=pngres)
				par(mfrow=c(1,1), mar=c(3.25,3.5,1,1), oma=c(0,0,0,0), mgp=c(2,0.75,0))
				if (all(logRecDevForACF==0)) {
					plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
					text(0,0,"No ACF plot.\nAll recruitment deviations = 0",cex=1.5,col="red")
				} else {
					acf(logRecDevForACF, lag.max=30, main="", xlab=linguaFranca("Lag",l), ylab=bquote(.(linguaFranca("Auto-correlation function of ",l)) ~ italic(epsilon[t])), na.action=na.pass)
				}
				if (p %in% c("eps","png")) dev.off()
			} ## end p (ptypes) loop
		}; eop()
	}
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.recdevacf


## plt.selectivity----------------------2019-07-18
## Transferred selectivity code from PBSscape.r
## into plt function (RH 190718)
## -----------------------------------------AME|RH
plt.selectivity <- function( obj, mainTitle="Rockfish", 
   ptypes=tcall(PBSawatea)$ptype, pngres=400, lang=c("e","f"))
{
	## Plot the selectivity.
	ageP = obj$Sel$P; names(ageP)=obj$Sel$Age
	selP = split(ageP,paste(obj$Sel$Series,obj$Sel$Sex,sep="."))
	xmax = max(as.numeric(sapply(selP,function(x){names(x[is.element(x,1)])[1]})),na.rm=TRUE) #maximum minimum age when P first hits 1
	if (is.na(xmax)) xmax = obj$extra$general$Nages
	#if (any(round(unlist(obj$extra$parameters[c("log_varRest","log_surveyvarR")]),5)<100)) ## temporary fix for dome-selectivity
	## Check for dome-selectivity estimation
	if ( any(unlist(sapply(currentRes$extra$priors[c("log_varRest_prior","log_surveyvarR_prior")],function(x){x[,1]}))>0) )
		xmax = obj$extra$general$Nages
	obj$Sel = obj$Sel[obj$Sel$Age <= xmax,]
	sel.f = obj$Sel

	## linguafranca takes too long on big vectors so shorten to unique
	user.f  = unique(sel.f$Series)
	## scape::plotSel looks for records labelled "Maturity" so cannot convert to french
	not.mat = !is.element(user.f,"Maturity")
	user.f  = linguaFranca(user.f[not.mat],"f")
	not.maturity = !is.element(sel.f$Series,"Maturity")
	sel.f$Series[not.maturity] = user.f[sel.f$Series[not.maturity]]

	usex.f  = unique(sel.f$Sex)
	usex.f  = linguaFranca(usex.f,"f")
	sel.f$Sex = usex.f[sel.f$Sex]

	obj.f = obj
	obj.f$Sel = sel.f
	fout = fout.e = "selectivity"
	for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
		changeLangOpts(L=l)
		fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
		if (is.null(ptypes)) ptypes="win"
		for (p in ptypes) {
			if (p=="eps")      postscript(paste0(fout,".eps"), width=6.5, height=4.5, horizontal=FALSE,  paper="special")
			else if (p=="png") png(paste0(fout,".png"), units="in", res=pngres, width=6.5, height=4.5)
			## regular par settings are ignored in lattice::xyplot -- need to fiddle with lattice settings
			#par(mfrow=c(1,1), mar=c(3.2,3.2,0.5,1), oma=c(0,0,0,0), mgp=c(2,0.75,0))
			par.sel = list(layout.widths = list(left.padding=-0.5, right.padding=-0.5),
				layout.heights=list(top.padding=0.5, main.key.padding=-0.75, bottom.padding=-2))
			if (l=="f")
				plotSel(obj.f, main=paste0("s\u{00E9}lectivit\u{00E9} du ", linguaFranca(mainTitle,l)), xlim=c(0,xmax), par.settings=par.sel)
			else
				plotSel(obj, main=paste0(mainTitle, " Selectivity"), xlim=c(0,xmax), par.settings=par.sel)
			if (p %in% c("eps","png")) dev.off()
		} ## end p (ptypes) loop
	}; eop()

}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.selectivity


#plt.ssbVbCatch-------------------------2011-08-31
# Plot the spawning stock and vulnerable biomass.
#----------------------------------------------AME
plt.ssbVbCatch <- function( obj, x1=1966, xLim=c(1954,2005), yLim=c(0,25000) )
{
  # Input an object from "load.allResFiles".
  # Plot spawning biomass, vulnerable biomass, catch against year.

  nPanels <- length(obj)

  par( oma=c(2,2,1,1), mar=c(2,2,1,1), mfrow=c(2,round(nPanels/2)) )

  for ( i in 1:nPanels )
  {
    res <- obj[[i]]
    year <- res$B$Year
    plot( year, res$B$SB, type="n", axes=FALSE,
      xlab="", xlim=xLim, ylab="", ylim=yLim )

    idx <- year >= x1
    lines( year[idx], res$B$SB[idx], lwd=2, lty=1 )
    lines( year[idx], res$B$VB[idx], lwd=2, lty=5 )

    yrTicks <- as.numeric( year[ year %% 10==0 ] )

    axis( side=1, at=yrTicks )
    axis( side=2 )
    axis( side=4, labels=FALSE )
    box()

    delta <- 0.35
    rect( year-delta, 0, year+delta, res$B$Y, density=-1, col="gray" )

    panLab( 0.5,0.95, names(obj)[i] )
    panLegend( 0.6, 0.95, legTxt=c("SB","VB"),
               ncol=1, lty=c(1,5), lwd=c(2,2), bty="n" )

    mfg <- par( "mfg" )
    if ( mfg[1]==mfg[3] & mfg[2]==mfg[4] | i==nPanels )
    {
      mtext( side=1, line=0.5, cex=1.0, outer=TRUE, "Year" )
      mtext( side=2, line=0.5, cex=1.0, outer=TRUE, "Metric tonnes" )
    }
  }
  par( mfrow=c(1,1) )
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.ssbVbCatch


## stdRes.index-------------------------2011-08-31
## Compute the standardised residuals.
## --------------------------------------------AME
stdRes.index <- function( obj, label=NULL, prt=TRUE )
{
  res <- log(obj$Obs) - log(obj$Fit)
  stdRes <- res / obj$CV
  result <- cbind( obj,stdRes )

  if ( prt )
  {
    sdRes <- sqrt( var( stdRes,na.rm=TRUE ) )
    cat( "\n",label,"\n" )
    cat( "\n     Std. Dev. of standardised Residuals=",sdRes,"\n" )
  }
  result
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~stdRes.index



## plt.stdResids------------------------2019-09-26
## Plot standardised residuals (AME adding xlim's for POP).
## Panels 2 & 3 show points labelled as years (RH 190926)
## -----------------------------------------AME|RH
plt.stdResids <- function( obj, pct=c(5,25,50,75,95),
   main=NULL, yLim=NULL, xLim=xLim, lang="e")
{
  ## Input is a data.frame with columns "Year", "stdRes", "Fit".
  ## cex.val=1.2       # size of text, if want to shrink down
  par( oma=c(0.5,1,2,0.5), mar=c(3,2,1,0.5), mfrow=c(3,1), mgp=c(1.5, 0.5, 0) )

  if ( is.null(yLim) )
    yLim <- range( obj$stdRes, na.rm=TRUE )
	ptcex=1; ptpch=19

  # Plot the standardised residuals against time.
  plot( obj$Year, obj$stdRes, type="n", xlab="", ylab="", ylim=yLim, xlim=xLim)
  abline( h=0, lty=2 )
  points( obj$Year, obj$stdRes, cex=ptcex, pch=ptpch) #, bg="orange" )
  mtext( side=1, line=2, cex=0.8, text=linguaFranca("Year",lang) )

  ## Plot the standardised residuals against predicted values.
  zin = !is.na(obj$Fit) & !is.na(obj$stdRes)
  plot(log(obj$Fit[zin]), obj$stdRes[zin], type="n", xlab="", ylab="", ylim=yLim )
  abline( h=0, lty=2 )
  #points( log(obj$Fit[zin]), obj$stdRes[zin], cex=ptcex, pch=ptpch) #, bg="orange" )
  #points( log(obj$Fit[zin]), obj$stdRes[zin], cex=3, pch=21, col="blue", bg="gainsboro") #, bg="orange" )
  text( log(obj$Fit[zin]), obj$stdRes[zin], labels=substring(obj$Year[zin],3,4), col="blue", font=2) #, bg="orange" )
  mtext( side=1, line=2, cex=0.8, text=linguaFranca("Predicted",lang) )

  ## Plot the q-normal plot of the standardised residuals.
  #wiggle = qqnorm(obj$stdRes, xlab="", ylab="", main="", pch=ptpch, cex=ptcex)
  stdres = obj$stdRes[zin]; names(stdres) = obj$Year[zin]
  wiggle = qqnorm(stdres, xlab="", ylab="", main="", pch=ptpch, cex=ptcex, type="n")
  abline(a=0, b=1)
  abline(h=quantile(stdres, p=pct/100, na.rm=TRUE), lty=c(3,2,1,2,3), lwd=0.75 )
#browser();return()
  #points(wiggle, pch=ptpch, cex=ptcex) #, bg="orange")
  text(wiggle, labels=substring(names(wiggle$y),3,4), col="blue", font=2) #, bg="orange" )
  mtext( side=1, line=2, cex=0.8, text=linguaFranca("Theoretical quantiles",lang) )

  mtext( side=2, line=-0.5, cex=0.8, outer=TRUE, text=linguaFranca("Standardised residuals",lang) )
  if ( !is.null(main) )
    mtext( side=3, line=0, cex=1.0, outer=TRUE, text=linguaFranca(main,lang) )
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.stdResids


## quantBox-----------------------------2018-04-03
##  Redefine boxplot to show quantiles (RH 150910)
##  http://r.789695.n4.nabble.com/Box-plot-with-5th-and-95th-percentiles-instead-of-1-5-IQR-problems-implementing-an-existing-solution-td3456123.html
##  Use PBStools solution without requiring PBStools.
##  This needs to be in `zzz.r' for package compilation but will repeat
##  in `plotFuns.r' for running code locally without loading package.
## ---------------------------------------------RH
local(envir=.PBSmodEnv,expr={
	myboxplot.stats <- function (x, coef=NULL, do.conf=TRUE, do.out=TRUE)
	{
		nna <- !is.na(x)
		n <- sum(nna)
		if (!exists("quants5"))
			quants5 = c(0.05,0.25,0.50,0.75,0.95)
		stats <- quantile(x, quants5, na.rm=TRUE) ## one day figure out how to make this dynamic
		iqr <- diff(stats[c(2, 4)])
		out <- x < stats[1] | x > stats[5]
		conf <- if (do.conf)
			stats[3] + c(-1.58, 1.58) * diff(stats[c(2, 4)])/sqrt(n)
		list(stats = stats, n = n, conf = conf, out = x[out & nna])
	}
	boxcode = deparse(boxplot.default)
	boxcode = gsub("boxplot\\.stats","tcall(myboxplot.stats)",boxcode)
	eval(parse(text=c("qboxplot=",boxcode)))
})

##quantBox------------------------------2016-03-24
## Redefine boxplot to show quantiles (RH 150910)
## http://r.789695.n4.nabble.com/Box-plot-with-5th-and-95th-percentiles-instead-of-1-5-IQR-problems-implementing-an-existing-solution-td3456123.html
## Use PBStools' 'quantBox' but call it 'quantbox'
##----------------------------------------------RH
quantbox = function (x, use.cols = TRUE, ...) ## taken from boxplot.matrix
{
	tget(qboxplot)
	if (rev(class(x))[1]=="matrix") {
		groups <- if (use.cols) 
			split(x, rep.int(1L:ncol(x), rep.int(nrow(x), ncol(x))))
		else split(x, seq(nrow(x)))
		if (length(nam <- dimnames(x)[[1 + use.cols]])) 
		names(groups) <- nam
		qboxplot(groups, ...)
	}
	else qboxplot(x, ...)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotbox


## splineCPUE---------------------------2019-01-03
## Fit spline curves through CPUE data to determine 
## optimal degrees of freedom (balance between rigorously
## fitting indices while not removing the majority of the signal)
## and to calculate CV process error from RSS at optimal DF.
## ---------------------------------------------RH
splineCPUE = function(dat, ndf=50,
   png=FALSE, pngres=400, PIN=c(8,8))
{
	## Collect residual sum of squares by degrees of freedom
	DF  = seq(2,nrow(dat),length.out=ndf)
	#DF  = seq(0,1,length.out=100)
	RSS = rep(NA,length(DF))#; names(RSS) = DF
	for (i in 1:length(DF)){
		ii = DF[i]
		RSS[i] = smooth.spline(dat[,"year"], dat[,"cpue"], df=ii, all.knots=TRUE)$pen.crit
	}
	dRSS = c(0,diff(RSS))
	df.opt = DF[findPV(min(dRSS),dRSS)]

	if (png) png("CPUEres-CVpro.png", units="in", res=pngres, width=PIN[1], height=PIN[2])
	expandGraph(mfrow=c(2,2), mar=c(3,3,0.5,0.5), cex=1)

	plot(DF, RSS, type="n")
	addLabel(0.5, 0.95, "Residual Sum of Squares", adj=0.5)
	abline(v=df.opt, col="green4", lty=3)
	lines(DF, RSS, col="red", lwd=2)

	plot(DF, dRSS, type="n")
	addLabel(0.5, 0.95, "Change in RSS (~slope)", adj=0.5)
	abline(v=df.opt, col="green4", lty=3)
	lines(DF, dRSS, col="blue", lwd=2)
	if (ndf<=50)
		points(DF,dRSS, pch=21, cex=.8, col="blue", bg="yellow")

	CVpro = sqrt(RSS[findPV(min(dRSS),dRSS)]/(nrow(dat)-2))/mean(dat$cpue)

	##  Plot the data and the 'optimal' fit
	plot(cpue ~ year, data = dat, pch=21, col="green4", bg="green", cex=1.1) #, main = "data(cpue) & smoothing splines")
	cpue.spl <- with(dat, smooth.spline(year, cpue, all.knots=TRUE))
	lines(cpue.spl, col="blue", lty=2, lwd=2)

	cpue.df <- smooth.spline(dat[,"year"], dat[,"cpue"], df=df.opt, all.knots=TRUE)
	lines(cpue.df, col="red", lty=1, lwd=2)
	addLegend(0.025, 0.97, legend=c(paste0("default df = ", round(cpue.spl$df,4)), paste0("with df = ", round(df.opt,4),", CVpro = ", round(CVpro,4))), col = c("blue","red"), lty = 2:1, bty="n", bg="transparent", adj=0)

	## Residual (Tukey Anscombe) plot:
	plot(residuals(cpue.df) ~ fitted(cpue.df), ylim=c(-0.35,0.45), pch=21, col="red", bg="pink", cex=1.1)
	abline(h = 0, col = "red")
	## consistency check:
	stopifnot(all.equal(dat$cpue, fitted(cpue.df) + residuals(cpue.df)))
	## The chosen inner knots in original x-scale :
	with(cpue.df$fit, min + range * knot[-c(1:3, nk+1 +1:3)]) # == unique(cpue$year)

	if (png) dev.off()
	return(list(DF=df.opt, RSS=RSS[findPV(min(dRSS),dRSS)], CVpro=CVpro))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~splineCPUE


#==============H I D D E N========================

#.compB0.get----------------------------2011-07-13
# This function (now in PBSawatea) loads the MCMCs from specified models.
# The MCMCs are typically available in the Sweave routine.
#-----------------------------------------------RH
.compB0.get =function(run.iter, path) { 
	#require(PBSawatea)
	B = list() #as.list(rep(NA,length(run.iter))); names(B)=run.iter
	for (i in run.iter) {
		run = strsplit(i,split="\\.")[[1]][1]
		mcmc.dir    <- paste(path,run,"/MCMC.",i,sep="")
		msy.dir     <- paste(mcmc.dir,"/MSY.",i,sep="")
		currentMCMC <- importMCMC( dir=mcmc.dir, quiet=FALSE )
		B[[i]][["B0.MCMC"]]     <- currentMCMC$B[,1]
		B[[i]][["Bt.MCMC"]]     <- currentMCMC$B[,dim(currentMCMC$B)[[2]]]
		currentMSY  <- msyCalc( dir=msy.dir, error.rep = 0 )
		B[[i]][["Bmsy.MCMC"]]   <- currentMSY$B
	}
	return(B) }

## Just for testing whether namespaces really work (RH)
.flushPlot=function(){
   plotBubbles(round(genMatrix(40,20),0),clrs=c("green","grey","red"))
   plotCI(1:10,runif(10,3,7),3,2)
   plotTrace(data.frame(this=rnorm(1000),that=rnorm(1000)))
}

