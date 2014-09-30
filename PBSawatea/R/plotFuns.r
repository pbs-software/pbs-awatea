#===============================================================================
# PBSawatea plot functions:
#  compB0           : compare reference points relative to B0.
#  compBmsy         : compare biomass posteriors relative to Bmsy.
#  cquantile        : cumulative quantile, from cumuplot.
#  cquantile.vec    : get one probability at a time.
#  plotBars         : barplots of specific year age proportions
#  plotBox          : modified boxplot with quantile whiskers.
#  plotDensPOP      : edited plotMCMC::plotDens function.
#  plotDensPOPpars  : edited plotMCMC::plotDens for parameters.
#  plotTracePOP     : plot traces with running median.
#===============================================================================

#compB0---------------------------------2011-12-15
# Compare reference points and criteria relative to B0.
#-----------------------------------------------RH
compB0=function(B, Mnams=NULL, ratios=c(0.4,0.8), 
   include=list(A1=TRUE, A2=TRUE, SSPM=TRUE, Bmsy=TRUE, Bt=TRUE),
   t.yr=2011, boxwidth=0.6, figgy=FALSE, width=12, height=9, ...) {

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
	if (figgy) figout = c("eps","pdf","pix","wmf","win") else figout="win"
	fout = paste("CompB0-",spp,"-(",paste(names(xBox),collapse=","),")",sep="")

	for (f in figout) {
		if (f=="eps"){    grDevices::ps.options(horizontal = TRUE)
		                  postscript(file=paste(fout,".eps",sep=""),width=width,height=height,fonts="mono") }
		if (f=="pdf"){    grDevices::ps.options(horizontal = TRUE)
		                  pdf(file=paste(fout,".pdf",sep=""),width=width,height=height,fonts="mono") }
		else if (f=="pix") png(paste(fout,".png",sep=""), units="in", res=300, width=width, height=height)
		else if (f=="wmf") win.metafile(paste(fout,".wmf",sep=""), width=width, height=height)
		par(mar=c(3.5,5,0.5,0.5),cex=ifelse(f%in%c("pix","eps"),1,1.2),xaxs="i")
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


#compBmsy-------------------------------2013-10-31
# Compare biomass posteriors relative to Bmsy.
#-----------------------------------------------RH
compBmsy = function(Bspp, spp="POP", Mnams=c("Est M","Fix M"),
   ratios=c(0.4,0.8), t.yr=2013, figgy=FALSE, width=12, height=9, 
   spplabs=TRUE, ...)
{
	oldpar = par(no.readonly=TRUE); oldpso = grDevices::ps.options()
	ciao = function() {
		par(oldpar)
		mess = paste("grDevices::ps.options(",paste(paste(names(oldpso),sapply(oldpso,deparse),sep="="),collapse=","),")",sep="")
		eval(parse(text=mess))
		gc(verbose=FALSE) }
	on.exit(ciao())

	spp = unique(spp)
	spp = spp[is.element(spp,names(Bspp))]
	if (length(spp)==0) 
		stop (paste("Input list 'Bspp' contains no species labelled: (\"",paste(spp,collapse="\",\""),"\")",sep=""))
	Nmods = sapply(Bspp,length); SPP = names(Nmods); nmods=sum(Nmods[spp])
	nrats = length(ratios)
	Bmsy=list()
	for (i in spp) {
		iBspp = Bspp[[i]]
		if (is.null(Mnams)) Mnams = names(iBspp)
		names(iBspp) = if (spplabs) paste(i,"\n",Mnams,sep="") else Mnams
		iBmsy = sapply(iBspp,function(x) { x[["Bt.MCMC"]] / x[["Bmsy.MCMC"]] }, simplify=FALSE) 
		for (j in names(iBmsy))
			Bmsy[[j]] = iBmsy[[j]]
	}
	Bmsy = Bmsy[rev(names(Bmsy))]

	ylim=c(0,max(sapply(Bmsy,quantile,0.96)))
	dots=list(...)
	unpackList(dots)
	if (!is.null(dots$medcol)) medcol=rev(medcol)
	if (!is.null(dots$boxfill)) boxfill=rev(boxfill)

	if (figgy) figout = c("eps","pdf","pix","wmf","win") else figout="win"
	fout = paste("CompBmsy-",paste(spp,collapse="+"),"-(",paste(gsub(" ","",Mnams),collapse=","),")",sep="")

	for (f in figout) {
		if (f=="eps"){    grDevices::ps.options(horizontal = FALSE)
		                  postscript(file=paste(fout,".eps",sep=""),width=width*1.25,height=height,fonts="mono",paper="special") }
		if (f=="pdf"){    grDevices::ps.options(horizontal = TRUE)
		                  pdf(file=paste(fout,".pdf",sep=""),width=width*1.25,height=height*1.25,fonts="mono") }
		else if (f=="pix") png(paste(fout,".png",sep=""), units="in", res=300, width=width, height=height)
		else if (f=="wmf") win.metafile(paste(fout,".wmf",sep=""), width=width*1.25, height=height*1.25)
		yaxis.space = (max(nchar(names(Bmsy)))-ifelse(spplabs,nchar(spp),0))^0.9
		par(mar=c(4,yaxis.space,0.5,0.5),cex=ifelse(f%in%c("pix","eps"),1,1.2),mgp=c(1.6,0.6,0))
		plotBox(Bmsy, horizontal=TRUE, las=1, xlim=c(0.5,nmods+1),ylim=ylim, cex.axis=1.2, yaxs="i", 
			pars=list(boxwex=boxwidth,medlwd=2,whisklty=1),quants=c(0.05,0.25,0.5,0.75,0.95))
		abline(v=ratios,col=rep(c("red","green4","blue"),nrats)[1:nrats],lty=2,lwd=2)
		plotBox(Bmsy, horizontal=TRUE, las=1, xlim=c(0.5,nmods+1),ylim=ylim, cex.axis=1.2, yaxs="i", 
			pars=list(boxwex=boxwidth,medlwd=2,whisklty=1,medcol=medcol,boxfill=boxfill,...),add=TRUE,quants=c(0.05,0.25,0.5,0.75,0.95))
		y2 = par()$usr[4] - 0.2*diff(par()$usr[3:4])
		text(c(0.2,0.6,1.2),rep(y2,3),c("Critical","Cautious","Healthy"),col=c("red","darkorange","green4"),font=2,cex=1.1,srt=90,adj=c(0,0.5))
		#axis(1,at=ratios,labels=paste(ratios,"",sep=""),tick=FALSE,cex.axis=1.2,line=2)
		text(ratios,par()$usr[3],labels=ratios,adj=c(1.1,-.5),col=c("red","green4"))
		#mess = paste("mtext(expression(italic(B)[",t.yr,"]/italic(B)[MSY]),side=1,line=2.5,cex=2)",sep="")
		mess = paste("mtext(expression(italic(B)[italic(t)]/italic(B)[MSY]),side=1,line=2.5,cex=1.5)",sep="")
#browser();return()
		eval(parse(text=mess))
		if (f!="win") dev.off()
	}
	invisible(Bmsy)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~compBmsy


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


#plotBars-------------------------------2013-09-11
# Plot barplots of specific year age proportions.
#-----------------------------------------------RH
plotBars = function(res, type="N", prop=TRUE, year=min(res[[type]][["Year"]]), 
    sex=c(2,1), # sex 2 =females (gfbio) 1 = males, 3 = unisex (awatea)
    age=NULL, fill=c("orange","cyan","green"), eps=FALSE, pix=FALSE, win=TRUE, ...) {

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
	figs = c(eps=eps,pix=pix,win=win)
#browser();return()
	for (k in names(figs)){
		if (!figs[k]) next
		if (k=="eps") postscript(paste(fnam,"eps",sep="."), horizontal=FALSE, paper="special", width=6.5, height=2.5*nrow)
		else if (k=="pix") png(paste(fnam,"png",sep="."), res=100, width=6.75*100, height=2.5*100*nrow, pointsize=12)
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
		if (any(k==c("eps","pix"))) dev.off()
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
    quants=c(0.025,0.25,0.5,0.75,0.975), outliers=FALSE) 
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
plotDensPOP = function (mcmc, probs = c(0.025, 0.975), points = FALSE, axes = TRUE, 
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
    mess = c(
    "require(grid, quietly = TRUE, warn.conflicts = FALSE)",
    "require(lattice, quietly = TRUE, warn.conflicts = FALSE)"
    )
    eval(parse(text=mess))
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
    function (mcmc, probs = c(0.025, 0.975), points = FALSE, axes = TRUE, 
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
    mess = c(
    "require(grid, quietly = TRUE, warn.conflicts = FALSE)",
    "require(lattice, quietly = TRUE, warn.conflicts = FALSE)"
    )
    eval(parse(text=mess))
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
    plot = TRUE, probs=c(0.025, 0.5, 0.975) , mpd=mcmc[1,], ...)  # AME probs
{
    panel.trace <- function(x, y, ...) {
        panel.xyplot(x, y, type = "l", lty = lty.trace, lwd = lwd.trace, 
            col = col.trace)
        if (any(is.finite(y)) && var(y) > 0) {
            # print(x)  # gives 1 2 3 ... 1000 for each parameter/yr
            # panel.xyplot(range(x), rep(median(y), 2), type = "l", 
            #  lty = lty.median, lwd = lwd.median, col = col.median)
            panel.xyplot(x, cquantile.vec(y, prob=0.025),
              type = "l", lty = lty.quant, lwd = lwd.quant,
              col = col.quant)
            panel.xyplot(x, cquantile.vec(y, prob=0.5),
              type = "l", lty = lty.median, lwd = lwd.median,
              col = col.median)
            panel.xyplot(x, cquantile.vec(y, prob=0.975),
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
    mess = c(
    "require(grid, quietly = TRUE, warn.conflicts = FALSE)",
    "require(lattice, quietly = TRUE, warn.conflicts = FALSE)"
    )
    eval(parse(text=mess))
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


#plt.catch------------------------------2014-09-25
# Small catch figures
# transferred from Sweave `run-master.Snw'
# <<catch, results=hide, echo=FALSE>>=
#-------------------------------------------AME/RH
plt.catch = function(years, Ct, xint=5, yint=250,
   ptypes=c("eps","png"), pngres=150)
{
	x = years[-length(years)]; xlim = range(x); xsmall = intersect(seq(1900,2100,xint),x)
	if (is.null(dim(Ct))) y = matrix(Ct,ncol=1) else y=Ct
	ylim = c(0,max(y)); ysmall = seq(yint,ylim[2],yint)
	ngear=ncol(y)
	pchGear=seq(21,20+ngear,1)
	colGear=rep(c("black","blue"),ngear)[1:ngear]
	for (p in ptypes) {
		if (p=="eps") postscript("catch.eps", width=6.5, height=4.5, horizontal=FALSE,  paper="special")
		else if (p=="png") png("catch.png", res=pngres, width=6.5*pngres, height=4.5*pngres)
		par(mfrow=c(1,1), mar=c(3.2,3.2,1,1), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
		#plot(x, y, type="h", xlim=xlim, ylim=ylim, xlab="Year", ylab="Catch (t)")
		xy = barplot(t(y),space=0.5,beside=FALSE,col=colGear,border="gainsboro",xlab="Year", ylab="Catch (t)",yaxs="i",names.arg=rep("",length(x)))
		lines(c(xy[1],rev(xy)[1]),c(0,0))
		axis(1, at=xy[match(xsmall,x)], tcl=-0.2, labels=xsmall,pos=0)
		axis(2, at=ysmall, tcl=-0.2, labels=FALSE)
		if (ngear>1) addLegend(0.05,0.80,fill=colGear,Cnames[1:ngear],yjust=0,bty="n")
#browser();return()
		dev.off()
	}
	for (p in ptypes) {
		if (p=="eps") postscript("catchSmall.eps", width=6, height=3, horizontal=FALSE,  paper="special")
		else if (p=="png") png("catchSmall.png", res=pngres, width=6*pngres, height=3*pngres)
		par(mfrow=c(1,1), mar=c(3.2,3.2,1,1), oma=c(0,0,0,0), mgp=c(2,0.75,0))
		plot(x, apply(y,1,sum), type="h", xlab="Year", ylab="Catch (t)")
		dev.off()
	}
}

#plt.biomass----------------------------2014-09-15
# Small biomass figures
# transferred from Sweave `run-master.Snw'
#<<Btplot, results=hide, echo=FALSE>>=
#<<BtB0plot, results=hide, echo=FALSE>>=
#-------------------------------------------AME/RH
plt.biomass = function(years, Bt, xint=5, yint=2500,
   ptypes=c("eps","png"), pngres=150, 
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
		else if (p=="png") png(paste0(pname,".png"), res=pngres, width=6*pngres, height=5*pngres)
		par(mfrow=c(1,1), mar=c(3.2,3.2,1,1), oma=c(0,0,0,0), mgp=c(2,0.75,0))
		plot(0,0, xlim=xlim, ylim=ylim, type="n", xlab=xlab, ylab=ylab)
		sapply(1:ngear, function(g,x,y){
			points(x, y[,g], pch=pchGear[g], col=colGear[g], bg="white", cex=0.8) }, x=x,y=y)
		tcl.val = -0.2
		axis(1, at=xsmall, tcl=tcl.val, labels=FALSE)
		axis(2, at=ysmall, tcl=tcl.val, labels=FALSE)
		dev.off()
	}
}

#plt.cpue-------------------------------2014-09-16
# Crude CPUE figure 
# transferred from Sweave `run-master.Snw'
#<<CPUEfig, results=hide, echo=FALSE>>=    # Crude for now
#-------------------------------------------AME/RH
plt.cpue = function(cpueObj, #xint=5, yint=2.5,
   ptypes=c("eps","png"), pngres=150 )
{
	zobs = !is.na(cpueObj$Obs)
	xlim = range(cpueObj$Year[zobs])
	ylim = range(c(cpueObj$Obs[zobs],cpueObj$Fit[zobs]))
	for (p in ptypes) {
		if (p=="eps") postscript("CPUEfit.eps", width=6, height=6, horizontal=FALSE,  paper="special")
		else if (p=="png") png("CPUEfit.png", res=pngres, width=6*pngres, height=5*pngres)
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
		dev.off()
	}
}

#plt.recdev-----------------------------2014-09-16
# Log recruitment deviations figure 
# transferred from Sweave `run-master.Snw'
#<<recdevplot, results=hide, echo=FALSE>>=
#-------------------------------------------AME/RH
plt.recdev = function(logRecDev, xint=5, #yint=0.1,
   ptypes=c("eps","png"), pngres=150 )
{
	x = as.numeric(names(logRecDev)); xlim = range(x); xsmall = intersect(seq(1900,2100,xint),x)
	for (p in ptypes) {
		if (p=="eps") postscript("recDev.eps", width=6.5, height=4, horizontal=FALSE,  paper="special")
		else if (p=="png") png("recDev.png", res=pngres, width=6*pngres, height=4*pngres)
		par(mfrow=c(1,1), mar=c(3.25,3.5,1,1), oma=c(0,0,0,0), mgp=c(2,0.75,0))
		plot(x, logRecDev, xlab="Year", ylab="Log recruitment deviations, epsilon_t")
		abline(h=0, col="grey")
		tcl.val = -0.2
		axis(1, at=xsmall, tcl=tcl.val, labels=FALSE)
		dev.off()
	}
}

#plt.recdevacf--------------------------2014-09-16
# Auto-correlation function of the log recruitment deviations
# transferred from Sweave `run-master.Snw'
#<<recdevacf, results=hide, echo=FALSE>>=
#-------------------------------------------AME/RH
plt.recdevacf = function(logRecDev, muC, logvC, A, years, yr1, 
   ptypes=c("eps","png"), pngres=150 )
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
	for (p in ptypes) {
		if (p=="eps") postscript("recDevAcf.eps", width=6.5, height=4, horizontal=FALSE,  paper="special")
		else if (p=="png") png("recDevAcf.png", res=pngres, width=6*pngres, height=5*pngres)
		par(mfrow=c(1,1), mar=c(3.25,3.5,1,1), oma=c(0,0,0,0), mgp=c(2,0.75,0))
		if (all(logRecDevForACF==0)) {
			plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
			text(0,0,"No ACF plot.\nAll recruitment deviations = 0",cex=1.5,col="red")
		} else {
			acf(logRecDevForACF, lag.max=30, main="", ylab="Auto-correlation function of epsilon_t", na.action=na.pass)
		}
		dev.off()
	}
}

#plt.initagedev-------------------------2014-09-16
# Initial age deviations figure
# transferred from Sweave `run-master.Snw'
#<<initagedevplot, results=hide, echo=FALSE>>=
#-------------------------------------------AME/RH
plt.initagedev = function(logInitAgeDev, 
   ptypes=c("eps","png"), pngres=150 )
{
	for (p in ptypes) {
		if (p=="eps") postscript("initAgeDev.eps", width=6.5, height=4, horizontal=FALSE,  paper="special")
		else if (p=="png") png("initAgeDev.png", res=pngres, height=5*pngres, width=6*pngres)
		par(mfrow=c(1,1), mar=c(3.2,3.2,1,1), oma=c(0,0,0,0), mgp=c(2,0.75,0))
		plot(names(logInitAgeDev), logInitAgeDev, xlab="Age", ylab="Log initial age deviations")
		abline(h=0, col="grey")
		dev.off()
	}
}

#plt.bubbles----------------------------2014-09-18
# Bubble plots of observed and fitted ages
# transferred from Sweave `run-master.Snw'
#<<bubbleplots, results=hide, echo=FALSE>>=
#-------------------------------------------AME/RH
plt.bubbles = function(mpdObj, nsex=2,
   ptypes=c("eps","png"), pngres=150 )
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
	SAnames   = tcall("PBSawatea")$Snames[tcall("PBSawatea")$SApos] # names of surveys with ages
	CAnames   = tcall("PBSawatea")$Cnames[tcall("PBSawatea")$CApos] # names of commercial gear with ages
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

	for (i in c("CAc","CAs")) {
		if (i=="CAc") { nr = length(CAnames); inames = CAnames }
		else          { nr = length(SAnames); inames = SAnames }
		for (j in c("Obs","Fit")) {
			if (nsex>1) kk = c("Fem","Male") else kk = "Fem"
			for (k in kk) {
				ijk = paste0(i,j,k); ijk.list = get(ijk)
				assign(ijk, ijk.list, pos=1)
				for (p in ptypes) {
					if (p=="eps") postscript(paste0(ijk,".eps"), width=6, height=ifelse(nr==1,6,8), horizontal=FALSE,  paper="special")
					else if (p=="png") png(paste0(ijk,".png"), res=pngres, width=6*pngres, height=ifelse(nr==1,6,8)*pngres)
					par(mfrow=c(nr,1), mar=c(2,3.5,2,0.5), oma=c(0,0,0,0), mgp=c(2,0.75,0))
					junk=sapply(1:nr,function(s,x,n){ # nr = no. rows = ngear or nsurv
						plotBubbles(x[[s]], dnam=TRUE, size=0.10, hide0=TRUE, main=n[s], prettyaxis=TRUE, las=1)
						mtext("Age",side=2,line=2,cex=1.2)},
						x=ijk.list, n=inames)
					dev.off()
				}
			}
		}
	}
	assign("residsCAcFem",sapply(CAcFitFem,function(x){prod(dim(x)+c(-1,0))}), pos=1)
	assign("residsCAsFem",sapply(CAsFitFem,function(x){prod(dim(x)+c(-1,0))}), pos=1)
	if (nsex>1) {
		assign("residsCAcMale",sapply(CAcFitMale,function(x){prod(dim(x)+c(-1,0))}), pos=1)
		assign("residsCAsMale",sapply(CAsFitMale,function(x){prod(dim(x)+c(-1,0))}), pos=1)
	}
	#browser();return()
	invisible()
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.bubbles

#plotAges-------------------------------2014-09-30
#  Plot the MPD model fits to age data (commercial
#  or survey) using the scape function `plotCA'.
#-------------------------------------------AME/RH
plotAges = function(obj, what="c", maxcol=4, sexlab=c("Females","Males"),
   ptypes = c("eps","png"), pngres=150)
{
	seriesType = paste0("CA",what)
	seriesList = sort( unique( obj[[seriesType]][["Series"]]) )
	if (what=="c") seriesName = tcall(PBSawatea)$Cnames[tcall(PBSawatea)$CApos]
	else           seriesName = tcall(PBSawatea)$Snames[tcall(PBSawatea)$SApos]
	CA.yrs  = sapply(split(obj[[seriesType]][["Year"]], obj[[seriesType]][["Series"]]), unique, simplify=FALSE)
	CA.nyrs = sapply(CA.yrs,length)

	for ( i in 1:length(seriesList) )  {
		ncols   = min(maxcol,max(CA.nyrs[i],1))
		nrows   = ceiling(CA.nyrs[i]/ncols)
		age.layout = rev(c(nrows,ncols)) # backwards in stupid lattice
		pwidth  = switch(ncols,4,6,6.5,6.5,6.5,6.5,6.5,6.5,6.5,6.5,6.5,6.5)
		pheight = switch(nrows,4,6,9,9,9,9,9,9,9,9,9,9)
		CA.sex = unique(obj[[seriesType]][["Sex"]])
		for(plot.sex in CA.sex) {
			j = grep(plot.sex,CA.sex)
			# legend key:
			CA.key = list(text=list(lab=c("Obs","Pred")), lines=list(col=c("black",ifelse(plot.sex=="Male","blue","red")),
				cex = c(0.5,0.5)), type=c("p","l"), x=0.78, y=ifelse(nrows==1,-0.10,-0.03), pch=c(20,20), lwd=1, between=0.3)
			# pch[2] doesn't get used, as type[2]="l". Have to match up if change options in plotCA(currentRes, ....)
			for (p in ptypes) {
				pname = paste0(ifelse(what=="c","ageComm","ageSurv"), plot.sex, i)
				if (p=="eps") postscript(paste0(pname,".eps"), width=pwidth, height=pheight, horizontal=FALSE,  paper="special", onefile=FALSE)
				else if (p=="png") png(paste0(pname,".png"), res=pngres, width=pwidth*pngres, height=pheight*pngres)
				par(mar=c(5,3.5,1,1), oma=c(0,0,0,0), mgp=c(2,0.5,0)) # this line may have no effect on `plotCA'
				plotCA( obj, what=what, ylab="Proportion", xlab="Age class", sex=plot.sex, layout= age.layout, key=CA.key, 
					main=paste0(seriesName[i]," - ",sexlab[j]), pch=20, cex.points=0.5, col.lines=ifelse(plot.sex=="Male","blue","red"), lwd.lines=2 ,series=i)
				if (p %in% c("eps","png")) dev.off()
			} # end of plot type loop
		} # end of plot.sex loop
	} # end of seriesList loop
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotAges


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
