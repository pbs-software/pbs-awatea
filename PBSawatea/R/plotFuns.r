##==============================================================================
## PBSawatea plot functions:
##  compB0          : compare reference points relative to B0
##  cquantile       : cumulative quantile, from cumuplot
##  cquantile.vec   : get one probability at a time
##  mochaLatte      : alternative to lattice plots (mockLattice)
##  plotACFs        : plot ACFs for the estimated parameters
##  plotAges        : plot the MPD model fits to age data
##  plotBars        : barplots of specific year age proportions
##  plotBox         : modified boxplot with quantile whiskers
##  plotChains      : plot cumulative fequency of 'nchains' by partitioning one trace
##  plotCI          : plot points with confidence intervals
##  plotCPUE        : plot CPUE and fit with error bars
##  plotDensPOP     : edited plotMCMC::plotDens function
##  plotDensPOPpars : edited plotMCMC::plotDens for parameters
##  plotMeanAge     : plot obs. & exp. mean ages from comm. & survey C@A data
##  plotSnail       : plot snail-trails for MCMC stock status
##  plotTracePOP    : plot traces with running median
##  plotTraj        : show all median trajectories (base+sens) in one figure
##  plt.biomass     : plot small biomass figures
##  plt.bubbles     : bubble plots of observed and fitted ages
##  plt.catch       : plot small catch figures
##  plt.cpue        : plot crude CPUE figure
##  plt.initagedev  : initial age deviations figure
##  plt.quantBio    : plot quantiles of reconstructed and projected biomass|recruits
##  plt.recdev      : log recruitment deviations figure 
##  plt.recdevacf   : auto-correlation function of the log recruitment deviations

##==============================================================================

#compB0---------------------------------2011-12-15
# Compare reference points and criteria relative to B0.
#-----------------------------------------------RH
compB0=function(B, Mnams=NULL, ratios=c(0.4,0.8), 
   include=list(A1=TRUE, A2=TRUE, SSPM=TRUE, Bmsy=TRUE, Bt=TRUE),
   t.yr=2011, boxwidth=0.6, figgy=FALSE, width=12, height=9, pngres=400, ...) {

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
	fout = paste("CompB0-",spp,"-(",paste(names(xBox),collapse=","),")",sep="")

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
		par(mar=c(3.5,5,0.5,0.5),cex=ifelse(f%in%c("png","eps"),1,1.2),xaxs="i")
		plotBox(xBox,xlim=xlim,ylim=ylim,yaxs="i",las=1,xaxt="n",yaxt="n",xlab="",ylab="",
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
					text(i,0.1,"Critical",srt=90,cex=1.4,col="red")
					text(i,0.3,"Cautious",srt=90,cex=1.4,col="blue") 
					text(i,ypos,"Schaefer\nsurplus\nproduction\nmodel",cex=cex.txt,adj=c(.5,1))
					abline(v = i + 0.5,lty=2,col="grey",lwd=2)
				}
				else {
					lines(xi,yi,lwd=3) # vertical zone delimiters
					text(i,0.15,"Endangered",srt=90,cex=1.4,col="chocolate3")
					text(i,ifelse(ii=="A1",0.4,0.6),"Threatened",srt=90,cex=1.4,col="purple")
					nCOS = nBars[intersect(names(nBars),c("A1","A2"))]
					text(sqrt(sum(nCOS)),ypos,"COSEWIC\ncriteria",cex=cex.txt,adj=c(.5,1))
					if (names(rev(nCOS))[1]==ii) abline(v = i + 0.5,lty=2,col="grey",lwd=2)
				}
			}
		}
		mess =paste(c("mtext(expression(paste(\"Reference criteria and points relative to ",ifelse(f%in%c("win","wmf"),"  ",""),
			"\",italic(B)[0],sep=\"\")),side=2,line=3.25,cex=",ifelse(f%in%c("win","wmf"),1.75,1.75),")"),collapse="")
		eval(parse(text=mess)) # plaster the mess along the y-axis
		
		abline(v = length(nBars)+(1:(nmods-1))*length(nBoxs)+0.5,lty=2,col="grey",lwd=2)
		nB = sum(nClr[c("Bmsy","Bt")])
		modlab = length(nBars)+median(1:nB)+nB*((1:nmods)-1)
		for (i in 1:nmods)
			text(modlab[i],ypos,paste(switch(i,"Estimating","Fixing"),"natural mortality"),cex=cex.txt,adj=c(.5,1))
		box()
		if (f!="win") dev.off()
	}
	invisible(list(BarBox=BarBox,xBox=xBox)) }
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~compB0


#cquantile------------------------------2010-10-20
# cumulative quantile, from cumuplot
#----------------------------------------------AME
cquantile <- function(z, probs)  
  {
  cquant <- matrix(0, nrow = length(z), length(probs))
  for (i in seq(along = z)) if (is.R())
    {
    cquant[i, ] <- quantile(z[1:i], probs = probs, names = FALSE)
    }
  else {
        cquant[i, ] <- quantile(z[1:i], probs = probs)
       }
  cquant <- as.data.frame(cquant)
  names(cquant) <- paste(formatC(100 * probs, format = "fg", 
      width = 1, digits = 7), "%", sep = "")
  return(cquant)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~cquantile


#cquantile.vec--------------------------2010-10-20
# AME doing this, just do one prob at a time 
# (so it returns a vector not a matrix)
#----------------------------------------------AME
cquantile.vec <- function(z, prob)  # cumulative quantile of vector
  {                                 #  prob is a single number
  cquant <- rep(NA, length(z))
  if(length(prob) != 1) stop("length prob should be 1")
  for (i in 1:length(z))
    {
    cquant[i] <- quantile(z[1:i], probs = prob, names = FALSE)
    }
  return(cquant)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~cquantile.vec


## Function: mochaLatte-----------------2018-04-24
##  An alternative to lattice plots (mockLattice)
##----------------------------------------------RH
mochaLatte = function(dat, xfld, yfld, ffld, panel, 
   strip=list(col=lucent("black",0.5), bg=lucent("moccasin",0.5), height=0.1, cex=1.4), ...)
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
	dots = list(...)
	if (is.null(dots$mar)) mar=c(3,3,0.5,0.5) else mar = dots$mar
	if (is.null(dots$oma)) oma=c(0,0,0,0)     else oma = dots$oma
	if (is.null(dots$mgp)) mgp=c(2,0.5,0)     else mgp = dots$mgp
	par(mfrow=rc, mar=mar, oma=oma, mgp=mgp)
	for (i in facs) {
		idat = dat[is.element(dat[,ffld],i),]
		if (is.null(dots$xlim)) xlim = range(idat[,xfld], na.rm=TRUE) else xlim = dots$xlim
		if (is.null(dots$ylim)) ylim = range(idat[,yfld], na.rm=TRUE) else ylim = dots$ylim
		ylim[2] = ylim[2] + diff(ylim)*strip$height ## add space for latte foam
		hzero = mar[2]==0  ## plots joined horizontally
		vzero = mar[1]==0  ## plots joined vertically
		do.call("plot", c(list(x=0, y=0, type="n", xlim=xlim, ylim=ylim, xaxt=ifelse(hzero,"n","s"), yaxt=ifelse(hzero,"n","s"), xlab="", ylab=""),
			dots[setdiff(names(dots),c("xlim","ylim","xlab","ylab"))]))
		if (hzero && par()$mfg[2]==1)
			do.call("axis", c(list(side=2), dots[setdiff(names(dots),c("xlim","ylim"))]))
		if (hzero && !vzero) {
			xticks = getlab(xlim, p=0.025)
			do.call("axis", c(list(side=1, at=xticks[["tck"]], labels=xticks[["sho"]]), dots[setdiff(names(dots),c("xlim","ylim"))]))
		}
		panel(x=idat[,xfld], y=idat[,yfld], ...)
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
		mtext(text=dots$xlab, side=1, line=0.25, outer=TRUE, cex=ifelse(is.null(dots$cex.lab),1.5,dots$cex.lab))
	if (!is.null(dots$ylab))
		mtext(text=dots$ylab, side=2, line=2.75, outer=TRUE, cex=ifelse(is.null(dots$cex.lab),1.5,dots$cex.lab))
#browser();return()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~mochaLatte


## plotACFs-----------------------------2017-05-29
##  Plot ACFs for the estimated parameters.
##  Control eps and png from PBScape.r in plt.mcmcGraphs
##----------------------------------------------RH
plotACFs =function(mcmcObj, lag.max=60) #, ptypes=tcall(PBSawatea)$ptype, pngres=400)
{
	#if (!is.null(dev.list())) on.exit(expandGraph(mfrow=c(1,1)))
	acfs  = apply(mcmcObj$P, 2, function(x){acf(x,plot=FALSE)$acf})
	ylim  = range(acfs[round(acfs,5)>-1 & round(acfs,5)<1])
	idx   = apply(mcmcObj$P, 2, allEqual)
	mcmcP = mcmcObj$P[,!idx,drop=FALSE]
	rc    = .findSquare(ncol(mcmcObj$P[,!idx]))
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
		mtext("Lag",side=1,outer=TRUE,line=1.5,cex=1.5)
		mtext("ACF",side=2,outer=TRUE,line=1.25,cex=1.5)
		#dev.off()
	#}
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotACFs


## plotAges-----------------------------2018-04-17
##  Plot the MPD model fits to age data
##  (commercial or survey) using the awkward 
##  scape function `plotCA'.
##------------------------------------------AME/RH
plotAges = function(obj, what="c", maxcol=5, sexlab=c("Females","Males"),
   ptypes=tcall(PBSawatea)$ptype, pngres=400, ...)
{
	seriesType = paste0("CA",what)
	seriesList = sort( unique( obj[[seriesType]][["Series"]]) )
	if (what=="c") seriesName = tcall(PBSawatea)$Cnames[tcall(PBSawatea)$CApos]
	else           seriesName = tcall(PBSawatea)$Snames[tcall(PBSawatea)$SApos]
	CA.yrs  = sapply(split(obj[[seriesType]][["Year"]], obj[[seriesType]][["Series"]]), unique, simplify=FALSE)
	CA.nyrs = sapply(CA.yrs,length)

	for ( i in 1:length(seriesList) )  {
		ii  = seriesList[i] ## sometimes a survey age series is missing
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
		CA.sex = unique(obj[[seriesType]][["Sex"]])
		for(plot.sex in CA.sex) {
			j = grep(plot.sex,CA.sex)
			## legend key:
			CA.key = list(text=list(lab=c("Obs","Pred")), lines=list(col=c("black",ifelse(plot.sex=="Male","blue","red")),
				cex = c(1,NA)), type=c("p","l"), x=ifelse(ncols>2,0.85,ifelse(ncols>1,0.8,0.75)), y=ifelse(nrows>2,-0.04,ifelse(nrows>1,-0.06,-0.12)), pch=c(20,NA), lwd=2, between=0.8)
			## pch[2] doesn't get used, as type[2]="l". Have to match up if change options in plotCA(currentRes, ....)
			for (p in ptypes) {
				pname = paste0(ifelse(what=="c","ageComm","ageSurv"), plot.sex,"Ser",ii) ## need the right series for plot
				set = if (!is.null(list(...)$set)) list(...)$set else ""
				pname = paste0(pname,set)
#browser();return()
				if (p=="eps") postscript(paste0(pname,".eps"), width=pwidth, height=pheight, horizontal=FALSE,  paper="special", onefile=FALSE)
				else if (p=="png") png(paste0(pname,".png"), units="in", res=pngres, width=pwidth, height=pheight)
				plotCA( obj, what=what, ylab="Proportion", xlab="Age class", sex=plot.sex, layout= age.layout, key=CA.key, 
					main=paste0(seriesName[i]," - ",sexlab[j]), pch=16, col.lines=ifelse(plot.sex=="Male","dodgerblue","red"), lwd.lines=2 , series=ii, ...)  ## need to pass the right series ii
				if (p %in% c("eps","png")) dev.off()
			} ## end of plot type loop
		} ## end of plot.sex loop
	} ## end of seriesList loop
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotAges


#plotBars-------------------------------2013-09-11
# Plot barplots of specific year age proportions.
#-----------------------------------------------RH
plotBars = function(res, type="N", prop=TRUE, year=min(res[[type]][["Year"]]), 
    sex=c(2,1), # sex 2 =females (gfbio) 1 = males, 3 = unisex (awatea)
    age=NULL, fill=c("orange","cyan","green"), eps=FALSE, png=FALSE, win=TRUE, pngres=400, ...) {

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
	fnam = paste("ageBars",paste(Sex,collapse=""),sep="")
	figs = c(eps=eps,png=png,win=win)
#browser();return()
	for (k in names(figs)){
		if (!figs[k]) next
		if (k=="eps") postscript(paste(fnam,"eps",sep="."), horizontal=FALSE, paper="special", width=6.5, height=2.5*nrow)
		else if (k=="png") png(paste(fnam,"png",sep="."), units="in", res=pngres, width=6.75, height=2.5*nrow)
		else resetGraph()
		par(mfcol=c(nrow,ncol),mgp=c(2,0.5,0),las=1,xaxs="i",
			mar = if(nyear==1) c(4,6,1,1) else c(4,3,1,1),
			oma = if(nyear==1) c(0,0,0,0) else c(0,3,0,0))
		for (i in year) {
			ii = as.character(i); aa=as.character(age)
			xy = barplot(t(mat[aa,,ii]),beside=TRUE,space=c(0,0), xaxt=ifelse(nage>15,"n","s"),
				col=fill, xlab="Age",ylab="",cex.lab=1.5,cex.axis=1.25,cex=1.25)
			ylab = ifelse(prop,"Proportions",type)
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
				lines(xpos,ypos,lwd=1,col=f)
			}
			if (nage>15) {
				AGE = seq(5,200,5)
				agelab = intersect(AGE,age)
				agepos = xpos[is.element(age,agelab)]
				axis(1,tick=FALSE,at=agepos,labels=agelab,mgp=c(2,0.1,0),cex.axis=1.25)
			}
			addLabel(0.95,0.95,ii,adj=c(1,1),cex=1.75,font=2)
			if (i==year[1])
				addLegend(0.975,0.875,legend=paste(Sex," M = ",M,sep=""),fill=fill,bty="n",yjust=1,xjust=1)
		}
		if (any(k==c("eps","png"))) dev.off()
	}
	invisible(list(dat=dat,mat=mat,xpos=xpos)) }
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotBars


#plotBox--------------------------------2011-12-15
# Modified boxplot with quantile whiskers.
#-----------------------------------------------RH
plotBox = function (x, ..., range=1.5, width=NULL, varwidth=FALSE, 
    notch=FALSE, outline=TRUE, names, plot=TRUE, 
    border=par("fg"), col=NULL, log="", 
    pars=list(boxwex=0.8, staplewex=0.5, outwex=0.5, whisklty=1), 
    horizontal=FALSE, add=FALSE, at=NULL,
    quants=get("quants5"), outliers=FALSE) 
{
	# RH tweaks for non-list inputs (e.g., vectors)
	if (!is.list(x)) {
		xnam = gsub(" ","",deparse(substitute(x)))
		x = list(x); attr(x,"names") <- xnam }
	#----------------------------------------------
    args <- list(x, ...)
    namedargs <- if (!is.null(attributes(args)$names)) 
        attributes(args)$names != ""
    else rep(FALSE, length.out = length(args))
    groups <- if (is.list(x)) 
        x
    else args[!namedargs]
    if (0L == (n <- length(groups))) 
        stop("invalid first argument")
    if (length(class(groups))) 
        groups <- unclass(groups)
    if (!missing(names)) 
        attr(groups, "names") <- names
    else {
        if (is.null(attr(groups, "names"))) 
            attr(groups, "names") <- 1L:n
        names <- attr(groups, "names")
    }
    cls <- sapply(groups, function(x) class(x)[1L])
    cl <- if (all(cls == cls[1L])) 
        cls[1L]
    else NULL
    for (i in 1L:n) groups[i] <- list(boxplot.stats(unclass(groups[[i]]), 
        range))
    stats <- matrix(0, nrow = 5L, ncol = n)
    conf <- matrix(0, nrow = 2L, ncol = n)
    ng <- out <- group <- numeric(0L)
    ct <- 1
    for (i in groups) {
        stats[, ct] <- i$stats
        conf[, ct] <- i$conf
        ng <- c(ng, i$n)
        if ((lo <- length(i$out))) {
            out <- c(out, i$out)
            group <- c(group, rep.int(ct, lo))
        }
        ct <- ct + 1
    }
    #----RH tweaks for qantile plots----
    stats = sapply(x,quantile,quants,na.rm=TRUE)
    if (!outliers) {
    	out = NULL; group = NULL }
    #-----------------------------------
    if (length(cl) && cl != "numeric") 
        oldClass(stats) <- cl
    z <- list(stats = stats, n = ng, conf = conf, out = out, 
        group = group, names = names)
    if (plot) {
        if (is.null(pars$boxfill) && is.null(args$boxfill)) 
            pars$boxfill <- col
        do.call("bxp", c(list(z, notch = notch, width = width, 
            varwidth = varwidth, log = log, border = border, 
            pars = pars, outline = outline, horizontal = horizontal, 
            add = add, at = at), args[namedargs]))
        invisible(z)
    }
    else z 
    invisible(x)}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotBox


##plotChains----------------------------2018-04-10
## Plots cumulative fequency of 'nchains' by partitioning one trace.
## Revised from 'plotTracePOP'
## mcmc=data.frame e.g, 'currentMCMC$P' from object created by 'importMCMC'.
## Very difficult to manipulate trellis plots (RH)
##------------------------------------------------
plotChains=function (mcmc, nchains=3, pdisc=0.1, 
   axes=FALSE, same.limits=FALSE, between=list(x=axes,y=axes), div=1, span=1/4,
   log=FALSE, base=10, main=NULL, xlab=NULL, ylab=NULL, 
   cex.main=1.2, cex.lab=1, cex.strip=0.8, cex.axis=0.8, las=0, 
   tck=0.4, tick.number=5, lty.trace=1, lwd.trace=1, col.trace="grey", 
   lty.median=1, lwd.median=1, col.median="black", lty.quant=2, lwd.quant=1, 
   col.quant="black", plot=TRUE, probs=get("quants3"), ...)
{
	panel.trace <- function(x, y, ...) {
		dots = list(...)
		unpackList(dots)
		if (is.null(dots$xlim)) xlim = range(x,na.rm=TRUE)
		if (is.null(dots$ylim)) ylim = range(y,na.rm=TRUE)
		abline (h=0.5, lty=3, lwd=1, col="gainsboro")
		chainlink = rep(1:nchains,ff)
		for (i in 1:nchains) {
			z = is.element(chainlink,i)
			lines(x[z], y[z], lty=lty.trace, lwd=2, col=rep(col.trace,nchains)[i])
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
	mochaLatte(dat,xfld="ValueSort",yfld="CumFreq",ffld="Factor", panel=panel.trace, ylim=c(0,1), mar=c(2,0,0,0), oma=c(1.5,4.5,0.5,1), tcl=-0.3, las=1, cex.axis=cex.axis, cex.lab=cex.lab, xlab=xlab, ylab=ylab)
	invisible(dat)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotChains


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
   col.quant="black", plot=TRUE, probs=get("quants3"), ...)  # AME probs
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


##plotCPUE------------------------------2018-04-13
## Plotting CPUE and fit with error bars
##  (copying plotIndexNotLattice).
## obj=currentRes$CPUE
##----------------------------------------------RH
plotCPUE <- function(obj, main="", save=NULL, bar=1.96, yLim=NULL,
   ptypes=tcall(PBSawatea)$ptype, pngres=400, ...)
{
	seriesList <- sort( unique( obj$Series ) )   # sort is risky if not always in same order
	nseries=length(seriesList)
	# surveyFigName =c("survIndGIG.eps", "survIndQCSsyn.eps", "survIndQCSshr.eps")
	surveyHeadName=c("CPUE")
	cvpro = tcall(PBSawatea)$cvpro
	if (is.null(cvpro) || all(cvpro==FALSE)) cvpro=0
	pwidth=6.0;  pheight=switch(nseries,5,8,9)
	for (p in ptypes) {
		pname = "CPUEser"
		if (p=="eps") postscript(paste0(pname,".eps"), width=pwidth, height=pheight, horizontal=FALSE,  paper="special", onefile=FALSE)
		else if (p=="png") png(paste0(pname,".png"), units="in", res=pngres, width=pwidth, height=pheight)
		par(mfrow=c(nseries,1),mai=c(0.75,0.75,0,0.1),omi=c(0,0,0.25,0),mgp=c(2,.75,0))
		# par(mai=c(0.25, 0.5, 0.3, 0.2)) # JAE changed  for each figure was for POP 0.45, 0.5, 0.3, 0.2
		# par(omi=c(0.45,0.1,0,0))  # Outer margins of whole thing, inch
		yrTicks=as.numeric( obj$Year)

		for ( i in 1:nseries )
		{
			ii=Nsurv + i  # to index the CPUE cvpro
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
				li=seriesVals$Lo, xlim=xLim, ylim=yLim, xlab="Year",
				ylab=paste("CPUE index:",seriesList[i]), gap=0, pch=21, col="blue", bg="cyan", lwd=2)
			lines(seriesVals$Year, seriesVals$Fit, lwd=2)
			axis( side=1, at=yrTicks, tcl=-0.2, labels=FALSE )
#browser();return()
			if (Ncpue==0)
				addLabel(0.95,0.95,"CPUE not used",adj=c(1,0),cex=0.8,col="slategrey")
			else if (is.numeric(cvpro[ii]) && round(cvpro[ii],5)!=0)
				addLabel(0.95,0.95,paste("+ CV process error ",cvpro[ii],sep=""),adj=c(1,0),cex=0.8,col="slategrey")
			# mtext( side=3, line=0.25, cex=0.8, outer=FALSE, surveyHeadName[i]) #  outer=TRUE
			# if(i==3)  mtext( side=2, line=-0.5, cex=1, outer=TRUE,"Relative biomass")
			# if(i==5)  mtext(side=1, line=0, cex=1, outer=TRUE, "Year")
		}  # cex was 0.8 for POP
	if (p %in% c("eps","png")) dev.off()
	}
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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotCI


#plotDensPOP----------------------------2010-10-19
#  editing plotMCMC::plotDens function to have
#  less whitesapce, not repeat x axis labels, and make y axes
#  the same scales. Can't just do through options. For Recruits
#  and Biomass. See plotDensPOPpar.r for parameters.
#  Tried y axes the same scales, but 1973-1975 are so narrow that
#  it makes all the others really small: same.limits=TRUE,
#  ylim=c(0,0.0005).
#  Andrew Edwards. Edited lines indicated by AME. 19 October 2010
#----------------------------------------------AME
plotDensPOP = function (mcmc, probs=get("quants3")[c(1,3)], points = FALSE, axes = TRUE, 
    same.limits = FALSE, between = list(x = axes, y = axes), 
    div = 1, log = FALSE, base = 10, main = NULL, xlab = NULL, 
    ylab = NULL, cex.main = 1.2, cex.lab = 1, cex.strip = 0.8, 
    cex.axis = 0.7, las = 0, tck = 0.5, tick.number = 5, lty.density = 1, 
    lwd.density = 3, col.density = "black", lty.median = 2, lwd.median = 1, 
    col.median = "darkgrey", lty.outer = 3, lwd.outer = 1, col.outer = "darkgrey", 
    pch = "|", cex.points = 1, col.points = "black", plot = TRUE,
    MPD.height = 0.04, mpd=mcmc[1,], ...)     #MPD.height, how far up to put MPD
{
    panel.dens <- function(x, ...) {     # x here seems to a vector
        if (any(is.finite(x)) && var(x) > 0)        # for each panel
            panel.densityplot(x, lty = lty.density, lwd = lwd.density, 
                col.line = col.density, plot.points = points, 
                pch = pch, cex = cex.points, col = col.points, 
                ...)
        else panel.densityplot(x, type = "n", ...)
        
        panel.abline(v = quantile(x, probs = probs), lty = lty.outer, 
            lwd = lwd.outer, col = col.outer)
        panel.abline(v = median(x), lty = lty.median, lwd = lwd.median, 
            col = col.median)
         # scan(); print(current.panel.limits()$ylim[2])  - max of y
         # print(graph$y.limits) is list of all panels
        panel.xyplot(mpd[panel.number()], current.panel.limits()$ylim[2]*MPD.height,
                     pch=19, col="red") # AME, MPD. 0.04 of way up
                                        #  assumes ..ylim[1]=0
        panel.xyplot(mpd[panel.number()], current.panel.limits()$ylim[2]*MPD.height,
                     pch=1, col="black") #AME
        # scan(); print(summary(x))    # Yes, here x is just vector
        
    }
    relation <- if (same.limits) 
        "same"
    else "free"
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
    mymain <- list(label = main, cex = cex.main)
    myxlab <- list(label = xlab, cex = cex.lab)
    myylab <- list(label = ylab, cex = cex.lab)
    myrot <- switch(as.character(las), `0` = 0, `1` = 0, `2` = 90, 
        `3` = 90)
    myscales <- list(y = list(draw = FALSE, relation = "free"), 
        x = list(draw = axes, relation = "same", cex = cex.axis, 
            tck = tck,  rot = myrot,
            alternating = TRUE))
            # AME: for y, relation = "same" -> relation = "free"
            # AME: for x, draw = axes -> draw = FALSE, but then no
            #       marks, so back to axes (which =TRUE)
            #       alternating = TRUE, relation="same"
            #      at=c(0, 50000, 100000, 150000, 200000))/1000
            #      took out tick.number = tick.number,  
    
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
#browser()
    if (plot) {
        print(graph) #panel.height=10, panel.width=10)
        invisible(x)
    }
    else {
        invisible(graph)
    }
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotDensPOP
#mcmcObj=currentMCMC
#plotDensPOP(mcmcObj$B[,getYrIdx(names(mcmcObj$B))]/1000, xlab="Female spawning biomass, Bt (1000 t)", 
#	between=list(x=0.2, y=0.2), ylab="Density", lwd.density=2, #panel.height=list(x=rep(1,5),unit="inches"), #*****Needs resolving
#	same.limits=TRUE, lty.outer=2, mpd=mpd.B[getYrIdx(names(mcmcObj$B))]/1000) #, layout=c(4,5)) 


#plotDensPOPpars------------------------2010-10-26
# editing plotMCMC::plotDens for parameters,
#  to put MPDs on. AME. 26th Oct 2010.
#----------------------------------------------AME
plotDensPOPpars =
    function (mcmc, probs=get("quants3")[c(1,3)], points = FALSE, axes = TRUE, 
    same.limits = FALSE, between = list(x = axes, y = axes), 
    div = 1, log = FALSE, base = 10, main = NULL, xlab = NULL, 
    ylab = NULL, cex.main = 1.2, cex.lab = 1, cex.strip = 0.8, 
    cex.axis = 0.7, las = 0, tck = 0.5, tick.number = 5, lty.density = 1, 
    lwd.density = 3, col.density = "black", lty.median = 2, lwd.median = 1, 
    col.median = "darkgrey", lty.outer = 3, lwd.outer = 1, col.outer = "darkgrey", 
    pch = "|", cex.points = 1, col.points = "black", plot = TRUE,
    MPD.height = 0.04, mpd=mcmc[1,],  ...)    # MPD.height, how far up to put MPD
{
    panel.dens <- function(x, ...) {
        if (any(is.finite(x)) && var(x) > 0) 
            panel.densityplot(x, lty = lty.density, lwd = lwd.density, 
                col.line = col.density, plot.points = points, 
                pch = pch, cex = cex.points, col = col.points, 
                ...)
        else panel.densityplot(x, type = "n", ...)

        panel.abline(v = quantile(x, probs = probs), lty = lty.outer, 
            lwd = lwd.outer, col = col.outer)
        panel.abline(v = median(x), lty = lty.median, lwd = lwd.median, 
            col = col.median)
                panel.xyplot(mpd[panel.number()], current.panel.limits()$ylim[2]*MPD.height,
                     pch=19, col="red") # AME, MPD. 0.04 of way up
                                        #  assumes ..ylim[1]=0
        panel.xyplot(mpd[panel.number()], current.panel.limits()$ylim[2]*MPD.height,
                     pch=1, col="black") #AME
    }
    relation <- if (same.limits) 
        "same"
    else "free"
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
    mymain <- list(label = main, cex = cex.main)
    myxlab <- list(label = xlab, cex = cex.lab)
    myylab <- list(label = ylab, cex = cex.lab)
    myrot <- switch(as.character(las), `0` = 0, `1` = 0, `2` = 90, 
        `3` = 90)
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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotDensPOPpars


## plotMeanAge--------------------------2018-04-16
##  Plot observed and expected mean ages from 
##  commercial and survey C@A data.
##----------------------------------------------RH
plotMeanAge =function(obj, useCA=TRUE, useSA=TRUE, CAnames)
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
		par(mfrow=c(nseries,1), mar=c(1.75,3,2,0.5), oma=c(2.5,2.5,0,0), mgp=c(2,0.75,0))
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
						mtext(CAnames[i], side=3, line=0.25, cex=1, outer=FALSE)
					} else {
						surveyHeadName = if (!exists("tcall")) ssnames[MAs$J] else tcall(PBSawatea)$Snames[MAs$J]
						mtext(surveyHeadName[i], side=3, line=0.25, cex=1, outer=FALSE)
					}
			}
			mtext("Mean Age (y)",side=2,line=0.5,cex=1.5,outer=TRUE)
			if (par()$mfg[1]==par()$mfg[3]) mtext("Year",side=1,line=1.0,cex=1.5,outer=TRUE)
		} ## end MAp
#browser();return()
		assign("MA.pjs",MA.pjs,envir=.GlobalEnv)
		dump("MA.pjs",file="MA.pjs.r")  ## for Paul
		save("MA.pjs",file="MA.pjs.rda")
	} ## end nseries
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotMeanAge


## plotSnail----------------------------2018-05-04
## Plot snail-trail plots for MCMC analysis.
##  AME: replacing "2010" with as.character(currYear - 1)
##  RH: added assYrs = years past with estimated Bcurr
##      from previous assessment(s), e.g., 5ABC QCS c(2011, 2017)
## -----------------------------------------AME/RH
plotSnail=function (BoverBmsy, UoverUmsy, p=c(0.1,0.9), xLim=NULL, yLim=NULL, 
	Lwd=1.5, ngear=1, assYrs=2011, outs=FALSE) ## outs = outliers
{
	## BU -- B = spawning biomass, U = harvest rate (or exploitation rate)
	BUlist = as.list(0:ngear); names(BUlist)=c("Spawning Biomass",Cnames[1:ngear])
	#BUlist[[1]] = BoverBmsy[,-length(BoverBmsy)]
	BUlist[[1]] = BoverBmsy[,-1]  ## conversation with PJS: we both agree that B2017/Bmsy should be paired with U2016/Umsy
	for (g in 1:ngear) {
		gfile = UoverUmsy[,grep(paste0("_",g),names(UoverUmsy))]
		names(gfile) = substring(names(gfile),1,4)
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
		xlab = expression(paste(italic(B[t])/italic(B)[MSY])), 
		ylab = expression(paste(italic(u[t-1])/italic(u)[MSY])),
		cex.lab=1.25,cex.axis=1.0,las=1)
	abline(h=1, col=c("grey20"), lwd=2, lty=3)
	abline(v=c(0.4,0.8), col=c("red","green4"), lwd=2, lty=2)
	for (i in ngear:1) {
		lines(BUmed[[1]], BUmed[[i+1]], col=colSlime[i], lwd=Lwd)
		points(BUmed[[1]], BUmed[[i+1]], type="p", pch=19, col=colPal(nB))
		points(BUmed[[1]][1], BUmed[[i+1]][1], pch=21, col=1, bg=colStart[i], cex=1.2)
		xend = rev(BUmed[[1]])[1]
		yend = rev(BUmed[[i+1]])[1]
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
#browser();return()
	}
	if (ngear>1)  addLegend(0.95,0.80,legend=Cnames,lty=1,lwd=Lwd,col=colSlime,seg.len=4,xjust=1,bty="n",cex=0.8)
	box()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotSnail


#plotTracePOP---------------------------2012-08-23
# Now adding running median, and taking off overall
# median and lowess line. Using cquantile from cumuplot.
# Trying to add in the MPD as a big circle for
# trace plots. 20th Oct 2010 (20/10/2010!)
#----------------------------------------------AME
plotTracePOP = function (mcmc, axes = FALSE, same.limits = FALSE, between = list(x = axes, 
    y = axes), div = 1, span = 1/4, log = FALSE, base = 10, main = NULL, 
    xlab = NULL, ylab = NULL, cex.main = 1.2, cex.lab = 1, cex.strip = 0.8, 
    cex.axis = 0.8, las = 0, tck = 0.5, tick.number = 5, lty.trace = 1, 
    lwd.trace = 1, col.trace = "grey", lty.median = 1, lwd.median = 1, 
    col.median = "black", lty.quant = 2, lwd.quant = 1, col.quant = "black", 
    plot = TRUE, probs=get("quants3"), mpd=mcmc[1,], ...)  # AME probs
{
    panel.trace <- function(x, y, ...) {
        panel.xyplot(x, y, type = "l", lty = lty.trace, lwd = lwd.trace, 
            col = col.trace)
        if (any(is.finite(y)) && var(y) > 0) {
            # print(x)  # gives 1 2 3 ... 1000 for each parameter/yr
            # panel.xyplot(range(x), rep(median(y), 2), type = "l", 
            #  lty = lty.median, lwd = lwd.median, col = col.median)
            panel.xyplot(x, cquantile.vec(y, prob=get("quants3")[1]),
              type = "l", lty = lty.quant, lwd = lwd.quant,
              col = col.quant)
            panel.xyplot(x, cquantile.vec(y, prob=get("quants3")[2]),
              type = "l", lty = lty.median, lwd = lwd.median,
              col = col.median)
            panel.xyplot(x, cquantile.vec(y, prob=get("quants3")[3]),
              type = "l", lty = lty.quant, lwd = lwd.quant,
              col = col.quant)
            panel.xyplot(x[1], mpd[panel.number()], pch=19, col="red") # AME
            panel.xyplot(x[1], mpd[panel.number()], pch=1, col="black") 
                     # AME, based on plt.trace, assume x[1]=1
            # suppressWarnings(panel.loess(x, y, span = span,
            #  lty = lty.loess, lwd = lwd.loess, col =col.loess,...))
        }
    }
    relation <- if (same.limits) 
        "same"
    else "free"
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
    mymain <- list(label = main, cex = cex.main)
    myxlab <- list(label = xlab, cex = cex.lab)
    myylab <- list(label = ylab, cex = cex.lab)
    myrot <- switch(as.character(las), `0` = 90, `1` = 0, `2` = 0, 
        `3` = 90)
    myscales <- list(x = list(draw = FALSE), y = list(draw = axes, 
        relation = relation, cex = cex.axis, tck = tck, tick.number = tick.number, 
        rot = myrot))
    mystrip <- list(cex = cex.strip)

    graph <- xyplot(Value ~ Draw | Factor, panel = panel.trace, 
        data = x, as.table = TRUE, between = between, main = mymain, 
        xlab = myxlab, ylab = myylab, par.strip.text = mystrip, 
        scales = myscales, ...)
    if (plot) {
        print(graph)
        invisible(x)
    }
    else {
        invisible(graph)
    }
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotTracePOP


## plotTraj-----------------------------2018-05-16
## Show all median trajectories (base+sens) in one figure.
## ---------------------------------------------RH
plotTraj = function(dat, index, traj="B", y0=FALSE, lab.stock, 
   col=c("black","green4","blue","red","purple","orange"), 
   png=FALSE, pngres=400, PIN=c(8,8), ...)
{
	opar = par(no.readonly=TRUE); on.exit(par(opar))
	Ntraj     = length(traj)
	trajYears = startYear:currYear
	Nyears    = length(trajYears)
	run.rwts  = names(dat$currentMCMC.sens)[index]
	assign("run.rwts",run.rwts,envir=.GlobalEnv)
	Nruns     = length(run.rwts)
	trajMat   = array(NA, dim=c(Nyears,Nruns,Ntraj), dimnames=list(year=trajYears, run=run.rwts, traj=traj))
#browser();return()
	#trajMat   = array(NA, dim=c(Nyears,Nruns,Ntraj+ifelse("BtB0"%in%traj,1,0)), dimnames=list(year=trajYears, run=run.rwts, traj=if("BtB0"%in%traj) c(traj,"BmsyB0") else traj))
	Bmsy.q5   = list()  ## not always needed
	tcol      = rep(col,Nruns)[1:Nruns]
	legtxt    = sen.lab  ## defined in global data object 'stock'
	ylabs     = c("Spawning Biomass","Vulnerable Biomass","Recruitment","Exploitation Rate",expression(italic(B[t])/italic(B)[0]),"Unknown")
	names(ylabs) = c("B","VB","R","U","BtB0","NA")
	for (i in 1:Nruns) {
		ii = run.rwts[i]
		for (j in 1:Ntraj) {
			jj = traj[j]
			if (jj=="BtB0") {
				jdat = dat[["currentMCMC.sens"]][[ii]][["B"]]
				dat[["currentMCMC.sens"]][[ii]][[jj]] = sweep(jdat,1,jdat[,1],"/")
				## Add Bmsy stuff (thanks PJS for the complication); and then he changed his mind...
				dat[["currentMCMC.sens"]][[ii]][["BmsyB0"]] = dat[["currentMSY.sens"]][[ii]][["B"]]/jdat[,1]
				Bmsy.q5[[ii]] = quantile(0.8*dat[["currentMCMC.sens"]][[ii]][["BmsyB0"]],quants5)  ## this won't work if more than one commercial gear
			}
			jval = apply(dat[["currentMCMC.sens"]][[ii]][[jj]],2,median)  ## this won't work if more than one commercial gear
			#if (jj=="R") jval = log10(jval)
			trajMat[substring(names(jval),1,4),ii,jj] = jval
		}
	}
	if (missing(lab.stock)) lab.stock = istock
	if (png) png(paste0(lab.stock,".Sens.Traj.",paste0(traj,collapse="+"),".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
	par(mfrow=.findSquare(Ntraj), mar=c(3,3.5,0.5,0), oma=c(0,0,0,1), mgp=c(2,0.5,0))
	x = trajYears; xlim = range(x)
	for (j in 1:Ntraj) {
		jj   = traj[j]
		jmat = trajMat[,,jj]
		ylim = range(jmat,na.rm=TRUE)
		#if (jj %in% c("BtB0")) ylim[1]=0
		if (y0) ylim[1] = 0
		plot(0,0,xlim=xlim,ylim=ylim,type="n",xlab="",ylab="",log=ifelse(jj=="R","y",""))
		if (jj %in% c("BtB0"))
			abline(h=seq(0.2,1,0.2),col="gainsboro")
		sapply(ncol(jmat):1,function(i){
			ii  = run.rwts[i]
			lty = list(...)$lty
			## Add median BmsyB0 if available (PJS says it's too cluttered)
			#if (jj=="BtB0" && length(Bmsy.q5)>0)
			#	abline(h=Bmsy.q5[[ii]][3], col=tcol[i], lwd=1, lty=lty[i])
			y = jmat[,ii]
			if (is.null(lty)) lty=rep(1,Nruns)
			lines(x,y, col=tcol[i], lwd=ifelse(i==1,1.5,1), lty=lty[i])
		})
		mtext("Year", side=1, line=1.75, cex=1)
		#ylab = ifelse(jj=="B","Spawning Biomass",ifelse(jj=="VB","Vulneerable Biomass",ifelse(jj=="R","Recruitment",ifelse(jj=="U","Exploitation Rate","Unknown"))))
		mtext(ylabs[jj], side=2, line=1.8, cex=1.2)
		if (j==1) addLegend(0.02, 0, col=tcol, seg.len=3, legend=legtxt[1:Nruns], bty="n", xjust=0, yjust=0, lwd=1, ...)
		if (Ntraj==1) {
			axis(1,at=intersect(seq(1900,2500,5),x),labels=FALSE,tcl=-0.2)
			axis(1,at=intersect(seq(1900,2500,10),x),labels=FALSE)	
		}
		box()
	}
#browser();return()
	if (png) dev.off()
	assign("trajSens",trajMat,envir=.GlobalEnv)
	invisible()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotTraj


#plt.biomass----------------------------2014-09-15
# Small biomass figures
# transferred from Sweave `run-master.Snw'
#<<Btplot, results=hide, echo=FALSE>>=
#<<BtB0plot, results=hide, echo=FALSE>>=
#-------------------------------------------AME/RH
plt.biomass = function(years, Bt, xint=5, yint=2500,
   ptypes=tcall(PBSawatea)$ptype, pngres=400, 
   pname="Bt", xlab="Year", ylab="Spawning biomass (t), Bt")
{
	#pname = gsub("\\.mpd","",as.character(substitute(Bt)))
	x = years; xlim = range(x); xsmall = intersect(seq(1900,2100,xint),x)
	if (is.null(dim(Bt))) y = matrix(Bt,ncol=1) else y=Bt
	ylim = c(0,max(y)); ysmall = seq(yint,ylim[2],yint)
	ngear=ncol(y)
	pchGear=seq(21,20+ngear,1)
	colGear=rep(c("black","blue"),ngear)[1:ngear]
	for (p in ptypes) {
		if (p=="eps") postscript(paste0(pname,".eps"), width=6, height=5, horizontal=FALSE,  paper="special")
		else if (p=="png") png(paste0(pname,".png"), units="in", res=pngres, width=6, height=5)
		par(mfrow=c(1,1), mar=c(3.2,3.2,1,1), oma=c(0,0,0,0), mgp=c(2,0.75,0))
		plot(0,0, xlim=xlim, ylim=ylim, type="n", xlab=xlab, ylab=ylab)
		sapply(1:ngear, function(g,x,y){
			points(x, y[,g], pch=pchGear[g], col=colGear[g], bg="white", cex=0.8) }, x=x,y=y)
		tcl.val = -0.2
		axis(1, at=xsmall, tcl=tcl.val, labels=FALSE)
		axis(2, at=ysmall, tcl=tcl.val, labels=FALSE)
		if (p %in% c("eps","png")) dev.off()
	}
}


#plt.bubbles----------------------------2014-09-18
# Bubble plots of observed and fitted ages
# transferred from Sweave `run-master.Snw'
#<<bubbleplots, results=hide, echo=FALSE>>=
#-------------------------------------------AME/RH
plt.bubbles = function(mpdObj, nsex=2,
   ptypes=tcall(PBSawatea)$ptype, pngres=400, redo.Graphs=TRUE)
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
				assign(ijk, ijk.list, pos=1)
				if (redo.Graphs) {
					for (p in ptypes) {
						if (p=="eps") postscript(paste0(ijk,".eps"), width=6, height=ifelse(nr==1,6,8), horizontal=FALSE,  paper="special")
						else if (p=="png") png(paste0(ijk,".png"), units="in", res=pngres, width=6, height=ifelse(nr==1,6,8))
						par(mfrow=c(nr,1), mar=c(2,3.5,2,0.5), oma=c(0,0,0,0), mgp=c(2,0.75,0))
						junk=sapply(1:nr,function(s,x,n){ # nr = no. rows = ngear or nsurv
							plotBubbles(x[[s]], dnam=TRUE, size=0.10, hide0=TRUE, main=n[s], prettyaxis=TRUE, las=1)
							mtext("Age",side=2,line=2,cex=1.2)},
							x=ijk.list, n=inames)
						if (p %in% c("eps","png")) dev.off()
					}
				}
			}
		}
	}
	if (useCA) {
		assign("residsCAcFem",sapply(CAcFitFem,function(x){prod(dim(x)+c(-1,0))}), pos=1)
		if (nsex>1)
			assign("residsCAcMale",sapply(CAcFitMale,function(x){prod(dim(x)+c(-1,0))}), pos=1)
	}
	if (useSA) {
		assign("residsCAsFem",sapply(CAsFitFem,function(x){prod(dim(x)+c(-1,0))}), pos=1)
		if (nsex>1)
			assign("residsCAsMale",sapply(CAsFitMale,function(x){prod(dim(x)+c(-1,0))}), pos=1)
	}
	#browser();return()
	invisible()
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.bubbles


#plt.catch------------------------------2014-09-25
# Small catch figures
# transferred from Sweave `run-master.Snw'
# <<catch, results=hide, echo=FALSE>>=
#-------------------------------------------AME/RH
plt.catch = function(years, Ct, xint=5, yint=250,
   ptypes=tcall(PBSawatea)$ptype, pngres=400)
{
	x = years[-length(years)]; xlim = range(x)
	xsmall = intersect(seq(1900,2100,xint),x)
	if (is.null(dim(Ct))) y = matrix(Ct,ncol=1) else y=Ct
	ylim = c(0,max(y)); 
	ysmall = seq(yint,ylim[2],yint)
	ngear=ncol(y)
	pchGear=seq(21,20+ngear,1)
	colGear=rep(c("black","blue"),ngear)[1:ngear]
	for (p in ptypes) {
		if (p=="eps") postscript("catch.eps", width=6.5, height=4.5, horizontal=FALSE,  paper="special")
		else if (p=="png") png("catch.png", units="in", res=pngres, width=6.5, height=4.5)
		par(mfrow=c(1,1), mar=c(3.2,3.2,1,1), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
		#plot(x, y, type="h", xlim=xlim, ylim=ylim, xlab="Year", ylab="Catch (t)")
		xy = barplot(t(y),space=0.5,beside=FALSE,col=colGear,border="gainsboro",xlab="Year", ylab="Catch (t)",yaxs="i",names.arg=rep("",length(x)))
		lines(c(xy[1],rev(xy)[1]),c(0,0))
		axis(1, at=xy[match(xsmall,x)], tcl=-0.2, labels=xsmall,pos=0)
		axis(2, at=ysmall, tcl=-0.2, labels=FALSE)
		if (ngear>1) addLegend(0.05,0.80,fill=colGear,Cnames[1:ngear],yjust=0,bty="n")
#browser();return()
		if (p %in% c("eps","png")) dev.off()
	}
	for (p in ptypes) {
		if (p=="eps") postscript("catchSmall.eps", width=6, height=3, horizontal=FALSE,  paper="special")
		else if (p=="png") png("catchSmall.png", res=pngres, width=6*pngres, height=3*pngres)
		par(mfrow=c(1,1), mar=c(3.2,3.2,1,1), oma=c(0,0,0,0), mgp=c(2,0.75,0))
		plot(x, apply(y,1,sum), type="h", xlab="Year", ylab="Catch (t)")
		if (p %in% c("eps","png")) dev.off()
	}
}


#plt.cpue-------------------------------2014-09-16
# Crude CPUE figure 
# transferred from Sweave `run-master.Snw'
#<<CPUEfig, results=hide, echo=FALSE>>=    # Crude for now
#-------------------------------------------AME/RH
plt.cpue = function(cpueObj, #xint=5, yint=2.5,
   ptypes=tcall(PBSawatea)$ptype, pngres=400)
{
	zobs = !is.na(cpueObj$Obs)
	xlim = range(cpueObj$Year[zobs])
	ylim = range(c(cpueObj$Obs[zobs],cpueObj$Fit[zobs]))
	for (p in ptypes) {
		if (p=="eps") postscript("CPUEfit.eps", width=6, height=6, horizontal=FALSE,  paper="special")
		else if (p=="png") png("CPUEfit.png", units="in", res=pngres, width=6, height=5)
		par(mfrow=c(1,1), mar=c(3.2,3.2,1,1), oma=c(0,0,0,0), mgp=c(2,0.75,0))
		plot(cpueObj$Year, cpueObj$Obs, xlim=xlim, ylim=ylim, type="n", xlab="Year",ylab="CPUE: Observed & Fit")
		series = unique(cpueObj$Series)
		nseries = length(series)
		for (i in 1:nseries) {
			ii = series[i]; z = is.element(cpueObj$Series,ii)
			points(cpueObj$Year[z], cpueObj$Obs[z], pch=21, bg=i+1, cex=1.2)
			lines(cpueObj$Year[z], cpueObj$Fit[z], col=i+1, lwd=2)
		}
		legend("topright",bty="n",col=(1:nseries)+1,lwd=2,legend=gsub("Series","CPUE",series),cex=0.8)
		if (p %in% c("eps","png")) dev.off()
	}
}


## plt.initagedev-----------------------2014-09-16
## Initial age deviations figure
## transferred from Sweave `run-master.Snw'
## <<initagedevplot, results=hide, echo=FALSE>>=
## -----------------------------------------AME/RH
plt.initagedev = function(logInitAgeDev, 
   ptypes=tcall(PBSawatea)$ptype, pngres=400)
{
	for (p in ptypes) {
		if (p=="eps") postscript("initAgeDev.eps", width=6.5, height=4, horizontal=FALSE,  paper="special")
		else if (p=="png") png("initAgeDev.png", units="in", res=pngres, height=5, width=6)
		par(mfrow=c(1,1), mar=c(3.2,3.2,1,1), oma=c(0,0,0,0), mgp=c(2,0.75,0))
		plot(names(logInitAgeDev), logInitAgeDev, xlab="Age", ylab="Log initial age deviations")
		abline(h=0, col="grey")
		if (p %in% c("eps","png")) dev.off()
	}
}


## plt.quantBio-------------------------2018-05-02
## Plot quantiles of reconstructed and projected biomass|recruits.
## ----------------------
## From popScapeRuns2.r -- AME now replaced yLim to force 0.
## This prints out tables (if run from command line), so be good to
##   use as template for decisions tables once we have MSY.
## --------------------------------------------AME
plt.quantBio <- function( obj, projObj=NULL, policy=NULL,
   p=get("quants5"), xyType="lines", lineType=c(3,2,1,2,3),
   refLines=NULL, xLim=NULL, yLim=NULL,
   userPrompt=FALSE, save=TRUE, yaxis.lab="Spawning biomass")
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
		mtext( side=1, line=2, cex=1.0, "Year" )
		mtext( side=2, line=2, cex=1.0, "Biomass" )
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
				mtext( side=1, line=0.5, cex=1.2, outer=TRUE, "Year" )   ## AME line=0 changed
				mtext( side=2, line=1.0, cex=1.2, outer=TRUE, yaxis.lab) ## " and Spawning
			}
			mtext( side=3, line=0.25, cex=0.9, paste0( "Catch strategy: ",policyList[j]) )

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


## plt.recdev---------------------------2018-04-24
## Log recruitment deviations figure 
## transferred from Sweave `run-master.Snw'
## <<recdevplot, results=hide, echo=FALSE>>=
##------------------------------------------AME/RH
plt.recdev = function(logRecDev, xint=5, #yint=0.1,
   ptypes=tcall(PBSawatea)$ptype, pngres=400)
{
	x = as.numeric(names(logRecDev)); xlim = range(x); xsmall = intersect(seq(1900,2100,xint),x)
	for (p in ptypes) {
		if (p=="eps") postscript("recDev.eps", width=6.5, height=4, horizontal=FALSE,  paper="special")
		else if (p=="png") png("recDev.png", units="in", res=pngres, width=6, height=4)
		par(mfrow=c(1,1), mar=c(3,3.5,1,1), oma=c(0,0,0,0), mgp=c(2,0.5,0))
		plot(x, logRecDev, type="n", xlab="Year", ylab=expression(paste("Log recruitment deviations, ", italic(epsilon[t]))), las=1)
		abline(h=0, lwd=1, lty=3, col="grey")
		lines(x, logRecDev, col="gainsboro")
		zpos = logRecDev > 0
		points(x[zpos], logRecDev[zpos],   pch=21, cex=0.8, col="blue", bg="cyan")
		points(x[!zpos], logRecDev[!zpos], pch=21, cex=0.8, col="red", bg="pink")
		tcl.val = -0.2
		axis(1, at=xsmall, tcl=tcl.val, labels=FALSE)
		if (p %in% c("eps","png")) dev.off()
	}
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.recdev

#plt.recdevacf--------------------------2014-09-16
# Auto-correlation function of the log recruitment deviations
# transferred from Sweave `run-master.Snw'
#<<recdevacf, results=hide, echo=FALSE>>=
#-------------------------------------------AME/RH
plt.recdevacf = function(logRecDev, muC, logvC, A, years, yr1, 
   ptypes=tcall(PBSawatea)$ptype, pngres=400, redo.Graphs=TRUE)
{
	ageHalfFemSelComm = min(round(muC - sqrt(exp(logvC) * log(2) ))) # take the min when Ngears>1
	assign("ageHalfFemSelComm", ageHalfFemSelComm, pos=1)
	# Take selectivity equation and find a which satisifies s_{ag1} = 0.5. 
	# Working out saved in POP12 folder.
	# This may be 1 year off, given I've now fixed the recruitment indexing issue.
	yearsForACF = max((yr1 - A + ageHalfFemSelComm), years[1]) : (as.numeric(max(names(logRecDev))) - ageHalfFemSelComm)
	assign("yearsForACF", yearsForACF, pos=1)
	# max() to start no earlier than first year of model
	logRecDevForACF = logRecDev[as.character(yearsForACF)]
	if (redo.Graphs) {
		for (p in ptypes) {
			if (p=="eps") postscript("recDevAcf.eps", width=6.5, height=4, horizontal=FALSE,  paper="special")
			else if (p=="png") png("recDevAcf.png", width=6, height=5, units="in", res=pngres)
			par(mfrow=c(1,1), mar=c(3.25,3.5,1,1), oma=c(0,0,0,0), mgp=c(2,0.75,0))
			if (all(logRecDevForACF==0)) {
				plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
				text(0,0,"No ACF plot.\nAll recruitment deviations = 0",cex=1.5,col="red")
			} else {
				acf(logRecDevForACF, lag.max=30, main="", ylab="Auto-correlation function of epsilon_t", na.action=na.pass)
			}
			if (p %in% c("eps","png")) dev.off()
		}
	}
}


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


