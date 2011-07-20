#.compB0.get----------------------------2011-07-13
# This function (now in PBSawatea) loads the MCMCs from specified models.
# The MCMCs are typically available in the Sweave routine.
#-----------------------------------------------RH
.compB0.get =function(run.iter, path) { 
	require(PBSawatea)
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
	return(B)
}
#run.iter=c("29.01","30.01")
#B = .compB0.get(run.iter, path="E:/GFish/PSARC11/YMR/Awatea/YMRrun")

#compB0---------------------------------2011-07-20
# Compare reference points and criteria relative to B0.
#-----------------------------------------------RH
compB0=function(B, Mnam=c("Est M","Fix M"), ratios=c(0.4,0.8), 
   incl.Bt=TRUE, boxwidth=0.6, figgy=FALSE, width=12, height=9, ...) {
	oldpar = par(no.readonly=TRUE); oldpso = grDevices::ps.options()
	ciao = function() {
		par(oldpar)
		paste("grDevices::ps.options(",paste(paste(names(oldpso),sapply(oldpso,deparse),sep="="),collapse=","),")",sep="")
		eval(parse(text=mess))
		gc(verbose=FALSE) }
	on.exit(ciao())
	x = list(A1=c(0.3,0.5),A2=c(0.5,0.7),SSPM=c(0.2,0.4,0.5)) # COSEWIC criteria A1 and A2, Schaefer surplus proction model
	M = sapply(B,function(x){x[["Bmsy.MCMC"]]/x[["B0.MCMC"]]},simplify=FALSE)
	S = sapply(B,function(x){x[["Bt.MCMC"]]/x[["B0.MCMC"]]},simplify=FALSE)
	M1 = sapply(M,function(x,r){xlst=list();for (j in 1:length(r)) xlst[[j]]=r[j]*x; names(xlst)=r; return(xlst)},r=ratios,simplify=FALSE)
	M2 = list()
	for (i in 1:length(B)) {
		M2 = c(M2,M1[[i]])
		if (incl.Bt) M2 = c(M2,list(Bt=S[[i]]))
	}
	if (is.null(Mnam)) {
		Mnam = rep("",length(B)); ratsep="" }
	else {
		Mnam = sapply(Mnam,function(x){paste("'",x,":'",sep="")}); ratsep = "~" }
	if (incl.Bt) ratlab=c(ratios,99) else ratlab=ratios
	Mnams=paste(rep(Mnam,each=length(ratlab)),ratsep,rep(ratlab,length(Mnam)),"Bmsy",sep="")
	Mnams = gsub("99Bmsy","B2011",gsub("1Bmsy","Bmsy",Mnams))
	names(M2)= Mnams 
	xM = xBox = c(x,M2)
	xBox[1:3]=NA
	xlim=c(0.25,length(xBox)+0.75); ylim=c(0,1)

	if (figgy) figout = c("eps","pdf","pix","wmf","win") else figout="win"
	fout = "CompB0"
	for (f in figout) {
		if (f=="eps"){    grDevices:::ps.options(horizontal = TRUE)
		                  postscript(file=paste(fout,".eps",sep=""),width=width,height=height,fonts="mono") }
		if (f=="pdf"){    grDevices:::ps.options(horizontal = TRUE)
		                  pdf(file=paste(fout,".pdf",sep=""),width=width,height=height,fonts="mono") }
		else if (f=="pix") png(paste(fout,".png",sep=""), units="in", res=300, width=width, height=height)
		else if (f=="wmf") win.metafile(paste(fout,".wmf",sep=""), width=width, height=height)
		par(mar=c(3.5,5,0.5,0.5),cex=ifelse(f%in%c("pix","eps"),1,1.2),xaxs="i")
		plotBox(xBox,xlim=xlim,ylim=ylim,yaxs="i",las=1,xaxt="n",yaxt="n",xlab="",ylab="",
			pars=list(boxwex=boxwidth,medlwd=2,whisklty=1,...)) 
		xaxislab = names(xBox)
		if (ratsep!="\n") {
			xaxislab = paste("expression(",paste(xaxislab,collapse=","),")",sep="")
			xaxislab = gsub("8B","8*B",gsub("4B","4*B",xaxislab))
			xaxislab = gsub("2011","[2011]",gsub("msy","[MSY]",xaxislab))
			xaxislab = gsub("B","italic(B)",xaxislab) 
		}
		mess = paste(c("axis(1,at=1:length(xBox),labels=",deparse(xaxislab),",tick=FALSE,padj=0.5,mgp=c(2,0.75,0),cex.axis=1.2)"),collapse="")
		if (ratsep!="\n")  mess = gsub("\\\"","",mess)
		eval(parse(text=mess))
		axis(2,at=seq(0,1,0.1),mgp=c(2,0.75,0),las=1,cex.axis=1.2)
		for (i in 1:3) {
			ylev = x[[i]]; ni = length(ylev)
			bxw=0.5*boxwidth
			xi=rep(c(i - bxw, i + bxw, NA),ni)
			yi = as.vector(sapply(as.list(ylev),function(x){c(rep(x,2),NA)}))
			if (i==3) ymax=median(ylev) else ymax=max(ylev)
			lines(c(c(i,i)- bxw,NA,c(i,i)+ bxw),c(0,ymax,NA,0,ymax),lwd=2,col="grey") # grey sides of bar
			if (i==3) 
				for (j in 1:2) lines(c(i-bxw,i+bxw),rep(ylev[j],2),lwd=3,col=switch(j,"red","blue","green4")) 
			else
				lines(xi,yi,lwd=3) # vertical zone delimiters
		}
		text(1:2,0.15,"Endangered",srt=90,cex=1.4,col="chocolate3")
		text(1:2,c(0.4,0.6),"Threatened",srt=90,cex=1.4,col="purple")
		text(3,0.1,"Critical",srt=90,cex=1.4,col="red")
		text(3,0.3,"Cautious",srt=90,cex=1.4,col="blue")
		mess =paste(c("mtext(expression(paste(\"Reference criteria and points relative to ",ifelse(f%in%c("win","wmf"),"  ",""),
			"\",italic(B)[0],sep=\"\")),side=2,line=3.25,cex=",ifelse(f%in%c("win","wmf"),1.75,1.75),")"),collapse="")
		eval(parse(text=mess))
		abline(v=c(2.5,3.5,6.5),lty=2,col="grey",lwd=2)
		ypos = par()$usr[4]-.025*diff(par()$usr[3:4])
		cex.txt = ifelse(f%in%"win",1.0,1.2)
		text(1.5,ypos,"COSEWIC criteria  ",cex=cex.txt,adj=c(.5,1))
		text(3,ypos,"Schaefer\nsurplus\nproduction\nmodel",cex=cex.txt,adj=c(.5,1))
		text(5,ypos,"Estimating natural mortality",cex=cex.txt,adj=c(.5,1))
		text(8,ypos,"Fixing natural mortality",cex=cex.txt,adj=c(.5,1))
		#text(5,ypos,expression(paste("Run 'Estimate ",italic(M),"'",sep="")),cex=cex.txt,adj=c(.5,1))
		#text(8,ypos,expression(paste("Run 'Fix ",italic(M),"'",sep="")),cex=cex.txt,adj=c(.5,1))
		box()
		if (f!="win") dev.off()
	}
	invisible(xBox)
}
#-------------------------------------------compB0
#medcol=whiskcol=staplecol=c("red","blue","black"); medcol[3]="green4"
#out=compB0(B,incl.Bt=T,boxwidth=0.6,figgy=F,medcol=medcol,whiskcol=whiskcol,staplecol=staplecol, boxfill=c("pink","lightblue1","honeydew"),whisklwd=2,staplelwd=2,Mnam=NULL) #c("M1","M2"))