.compB0.getStuff =function(run.iter) { # this is done in the Sweave file
	require(PBSawatea)
	#source("ymrscape.r")
	B = list() #as.list(rep(NA,length(run.iter))); names(B)=run.iter
	for (i in run.iter) {
		run = strsplit(i,split="\\.")[[1]][1]
		mcmc.dir    <- paste("E:/GFish/PSARC11/YMR/Awatea/YMRrun",run,"/MCMC.",i,sep="")
		msy.dir     <- paste(mcmc.dir,"/MSY.",i,sep="")
		#proj.dir    <<- paste(mcmc.dir,"Proj.90yrs",sep="/")
		currentMCMC <- importMCMC( dir=mcmc.dir, quiet=FALSE )
		#Bt.MCMC     <<- currentMCMC$B
		B[[i]][["B0.MCMC"]]     <- currentMCMC$B[,1]
		B[[i]][["Bt.MCMC"]]     <- currentMCMC$B[,dim(currentMCMC$B)[[2]]]
#browser();return()
		currentMSY  <- msyCalc( dir=msy.dir, error.rep = 0 )
		B[[i]][["Bmsy.MCMC"]]   <- currentMSY$B
		#currentProj <<- importProj( dir=proj.dir, quiet=FALSE )
	}
	return(B)
}
#run.iter=c("29.01","30.01"); 
#B = .compB0.getStuff(run.iter)

#compB0---------------------------------2011-07-13
# Compare reference points and criteria in B0 space.
#-----------------------------------------------RH
compB0=function(B,Mnam=c("Est M","Fix M"), ratios=c(0.4,0.8), 
   incl.Bt=TRUE, boxwidth=0.6, figgy=FALSE, width=12, height=9) {
	source("plotBox.r")
	oldpar=par(no.readonly=TRUE); on.exit(par(oldpar))
	x = list(A1=c(0.3,0.5),A2=c(0.5,0.7),SSPM=c(0.2,0.4,0.5)) # COSEWIC criteria A1 and A2, Schaeffer surplus proction model
	M = sapply(B,function(x){x[["Bmsy.MCMC"]]/x[["B0.MCMC"]]},simplify=FALSE)
	S = sapply(B,function(x){x[["Bt.MCMC"]]/x[["B0.MCMC"]]},simplify=FALSE)
	M1 = sapply(M,function(x,r){xlst=list();for (j in 1:length(r)) xlst[[j]]=r[j]*x; names(xlst)=r; return(xlst)},r=ratios,simplify=FALSE)
	M2 = list()
	for (i in 1:length(B)) {
		M2 = c(M2,M1[[i]])
		if (incl.Bt) M2 = c(M2,list(Bt=S[[i]]))
	}
	if (!is.null(Mnam)) {
		if (incl.Bt) ratlab=c(ratios,99) else ratlab=ratios
		Mnams=paste(rep(Mnam,each=length(ratlab)),"\n",rep(ratlab,length(Mnam)),"Bmsy",sep="")
		Mnams = gsub("99Bmsy","B2011",gsub("1Bmsy","Bmsy",Mnams))
		names(M2)= Mnams }
	xM = xBox = c(x,M2)
	xBox[1:3]=NA
	#ylim = c(0,max(sapply(xM,max)))
	xlim=c(0.25,length(xBox)+0.75); ylim=c(0,1)

	if (figgy) figout = c("eps","pix","wmf","win") else figout="win"
	fout = "CompB0"
	for (f in figout) {
		if (f=="eps"){    grDevices:::ps.options(horizontal = TRUE)
		                  postscript(file=paste(fout,".eps",sep=""),width=width,height=height,fonts="mono") }
		else if (f=="pix") png(paste(fout,".png",sep=""), units="in", res=300, width=width, height=height)
		else if (f=="wmf") win.metafile(paste(fout,".wmf",sep=""), width=width, height=height)
		par(mar=c(4,5,0.5,0.5),cex=ifelse(f%in%c("pix","eps"),1,1.2),xaxs="i")
		plotBox(xBox,xlim=xlim,ylim=ylim,yaxs="i",las=1,xaxt="n",yaxt="n",xlab="",ylab="",
			pars=list(boxwex=boxwidth,medlwd=2,whisklty=1,medcol="navyblue", boxfill=c("moccasin","lightcyan","#d2ffd2")))
#browser();return()
		axis(1,at=1:length(xBox),labels=names(xBox),tick=F,padj=0.5,mgp=c(2,1,0),cex.axis=1.1)
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
#if (i==3) {browser();return()}
		}
		text(1:2,0.15,"Endangered",srt=90,cex=1.4,col="red")
		text(1:2,c(0.4,0.6),"Threatened",srt=90,cex=1.4,col="blue")
		text(3,0.1,"Critical",srt=90,cex=1.4,col="red")
		text(3,0.3,"Cautious",srt=90,cex=1.4,col="blue")
		mtext("Reference Criteria / Points in Terms of B0",side=2,line=3.25,cex=ifelse(f%in%c("win","wmf"),1.75,1.75))
#browser();return()
#abline(h=.859466)
		abline(v=2.5,col="grey10",lwd=2)
		abline(v=c(3.5,6.5),lty=2,col="grey",lwd=2)
		box()
		if (f!="win") dev.off()
	}
	invisible(xBox)
}
#-------------------------------------------compB0
#out=compB0(B,incl.Bt=TRUE,boxwidth=0.6,figgy=F)