#-----------------------------------------------------------2013-09-13
#                           Menu Functions                           #
#--------------------------------------------------------------------#

get.resFile <- function( resFile=NULL )
{
  resList <- dir( path=".", pattern=".res$" )

  choice <- 1
  while ( choice > 0 )
  {
    cat( "\nCurrent res file:",resFile," Select new file (0 to exit):\n" )
    choice <- menu( resList )
    if ( choice > 0 )
    {
      resFile <- resList[choice]
      choice <- 0
    }
  }

   #  newRes <- importCol2( res.file=resFile, Dev=TRUE, CPUE=TRUE, Survey=TRUE, CAc=TRUE)
   #  AME changed to the following, for which the options agree
   #  with those in load.allResFiles, which works. Now it loads
   #  in okay, previously had trouble when overwriting.
   #  load.allResFiles does load in them all at the beginning, but
   #  the reloading here was probably done so an updated .res file
   #  can be loaded in. Also added CAc=TRUE, CAs=TRUE
	#  RH: removed importCol2 (2013-09-13) because importRes is preferred and maintained.
	#  RH: removed `plotMCMC::' from `plotTrace' and `plotDens'; regardless,
	#    PBSawatea uses the functions from plotMCMC over PBSmodelling and issues WARNINGS just because.

	newRes <- importRes( res.file=resFile, Dev=TRUE, CPUE=TRUE,
		Survey=TRUE, CLc=TRUE, CLs=TRUE, CAs=TRUE, CAc=TRUE)
	assign( "currentRes", newRes, pos=1 )
	assign( "resFile",resFile )
	cat( "\nLoaded Awatea res file: ",resFile,"\n\n" )
	#print( ll( currentRes ) ) ## RH: this is the only `gdata' function used and importing gdata causes problems

	resFile
}


mainMenu <- function()
{
  menuItems <- c( "Import files",
                  "MPD plots",
                  "Plot all MPD graphs",
                  "Save all MPD plots to PNG",
                  "MCMC plots",
                  "Plot all MCMC plots",
                  "Save all MCMC plots to PNG",
                  "Close all graphics windows",
                  "Help & Utilities" )

  choice <- 1
  while ( choice > 0 )
  {
    # cat( "\nCanary rockfish assessment 2007\n" )
    cat( "\nPacific ocean perch assessment 2010\n" )  # AME
    cat( "Select menu item by number (0 to exit):\n" )
    choice <- menu( menuItems )
    switch( choice,
      {
         loadMenu()
#        resFile <- get.resFile( resFile )
#        assign( "resFile",resFile,pos=1 )
#        currentRes <- importCol( res.file=resFile, Dev=TRUE, CPUE=TRUE, Survey=TRUE, CAc=TRUE )
#        assign( "currentRes", currentRes, pos=1 )
#        cat( "\nLoaded Awatea res file: \n\n" )
#        print( ll( currentRes ) )
      },
      mpdMenu(),
      {
        plt.mpdGraphs( currentRes, save=FALSE )
      },
      plt.mpdGraphs( currentRes, save=TRUE ),
      mcmcMenu(),
      plt.mcmcGraphs( currentMCMC, currentProj, save=FALSE ),
      plt.mcmcGraphs( currentMCMC, currentProj, save=TRUE ),
      closeAllWin(),
      utilMenu()
    )
  }
  choice
}


loadMenu <- function()
{
  menuItems <- c( "Get Awatea res file",
                  "Get Awatea MCMC file",
                  "Get Awatea projection file",
                  "Load all res files in working directory",
                  "Get PJS Delay Difference MCMC+Projection" )

  choice <- 1
  while ( choice > 0 )
  {
    cat( "\nSelect file for import\n" )
    cat( "Select menu item by number (0 to exit):\n" )
    choice <- menu( menuItems )
    switch( choice,
      {
        # Awatea res file.
        resFile <- get.resFile( resFile )
        assign( "resFile",resFile,pos=1 )
      },
      {
        # Awatea MCMC.
        currentMCMC <- importMCMC( dir=".", quiet=FALSE )
        assign( "currentMCMC", currentMCMC, pos=1 )
      },
      {
        # Awatea projection.
        currentProj <- importProj( dir=".", quiet=FALSE )
        assign( "currentProj", currentProj, pos=1 )
      },
      {
        tmp <- load.allResFiles()
        assign( "resFileList",tmp,pos=1 )
        cat( "\nAll Awatea res files in working directory loaded...\n" )
      },
      {
        # Delay difference.
        currentMCMC <- importMCMC.ddiff()
        assign( "currentMCMC", currentMCMC, pos=1 )
        currentProj <- importProj.ddiff()
        assign( "currentProj", currentProj, pos=1 )
      }
    )
  }
  choice
}


mpdMenu <- function()
{
  menuItems <- c( "Plot biomass, recruitment, catch",
                  "Plot numbers at age",
                  "Plot selectivity and maturity",
                  "Plot commercial catch-at-age results",
                  "Plot survey catch-at-age results",
                  #"Plot survey catch-at-length results",
                  "Plot abundance index",
                  "All residual plots",
                  "Plot multi-panel biomass, recruitment, catch",
                  "Plot multi-panel exploitation rate",
                  "Plot alternative numbers at age" )


  choice <- 1
  while ( choice > 0 )
  {
    cat( "\nModel fit plots (0 to exit):\n" )
    # cat( "\n*** WARNING: NOT AVAILABLE FOR PJS DDIFF *** \n" )#AME
    choice <- menu( menuItems )
    switch( choice,
      {
        if ( exists( "currentRes" ) )
          plotB2( currentRes,main=mainTitle )   # Doesn't actually agree with
                                                #  menu
      },
      {
        if ( exists( "currentRes" ) )
          plotN( currentRes, main=mainTitle, ages=c(2:60) )
          #plotN( currentRes, main=mainTitle)
      },
      {
        if ( exists( "currentRes" ) )
           currentResRed = currentRes
            # Reduced currentRes, just plotting Selectivity to 20
           currentResRed$Sel = currentResRed$Sel[currentResRed$Sel$Age < 21,]
          plotSel( currentResRed, main=mainTitle, xlim=c(0,20) )
      },
      {
        if ( exists( "currentRes" ) )
        {
          # NOTE: There is a bug in plotCA that prevents plotting multiple
          #       series given a list of character vectors in series. AME: ?
          seriesList <- sort( unique( obj$CAc$Series) )
          for ( i in 1:length(seriesList) )   # for catch data. AME - POP=1
          {
            # AME dividing years into 4 groups - note no 1985, 86,
            #  or 88 data
            windows()
            # plotCA( currentRes, series=i, main=paste("Comm",mainTitle,"Series",i), what="c" )
            plotCA( currentRes, series=i, main=paste("Comm",mainTitle,"Series",i), what="c", years=1978:1984 )
            windows()
            plotCA( currentRes, series=i, main=paste("Comm",mainTitle,"Series",i, "(no 1985, 86 or 88)"), what="c", years=1985:1994 )  # no 85, 86 or 88
            windows()
            plotCA( currentRes, series=i, main=paste("Comm",mainTitle,"Series",i), what="c", years=1995:2001 )
            windows()
            plotCA( currentRes, series=i, main=paste("Comm",mainTitle,"Series",i), what="c", years=2002:2009 )
          }
        }

      },
      # AME adding - plotting the CA survey data for all three series
      {
        if ( exists( "currentRes" ) )
        {
          seriesList <- sort( unique( obj$CAs$Series) )
          for ( i in 1:length(seriesList) )
          {
            windows()
            plotCA( currentRes, series=i, main=paste("Survey",mainTitle,"Series",i), what="s" )
          }
        }
      },
      #{
      #  if ( exists( "currentRes" ) )
      #  {
          # NOTE: There is a bug in plotCA that prevents plotting multiple
          #       series given a list of character vectors in series.

      #    seriesList <- sort( unique( obj$CLs$Series) )
      #    for ( i in 1:length(seriesList) )
      #    {
      #      windows()
      #      plotCL( currentRes, series=i, main=paste("Surv",mainTitle,"Series",i), what="c" )
      #    }
      #  }

      #},         
      {
        if ( exists( "currentRes" ) )
        {
        #ACH: I couldn't get plotIndex2 function to work in the limited time I had. AME changed these two to plotIndex
          plotIndex( currentRes, main=mainTitle, what="c", bar=1.96 )
          # Think that plots the CPUE, which we're not using (?)
          windows()
          plotIndex( currentRes, main=mainTitle, what="s", bar=1.96 )
        }
      },
      {
        if ( exists( "currentRes" ) )
        {
          plt.ageResids( stdRes.CA( obj$CAc ),main=paste(mainTitle,"Comm Ages"))
      # AME adding - plotting the CA residuals for all three surveys:
          seriesList <- sort( unique( obj$CAs$Series) )
          for ( i in 1:length(seriesList) )
          {
            windows()
            plt.ageResids( stdRes.CA( obj$CAs[
                obj$CAs$Series == i,] ),
                main=paste(mainTitle,"Survey",i))
          }
          # windows()
          #plt.lengthResids( stdRes.CL( obj$CLs ),main=paste(mainTitle,"Survey lengths"))
          #plt.idx( obj$CPUE,  main="Commercial Fishery CPUE")
               # AME - not for POP
          plt.idx( obj$Survey,main="Survey Indices")
        }
      },
      {
        plt.ssbVbCatch( resFileList )
      },
      {
        plt.expRate( resFileList )
      },
      {
        plt.numR( resFileList[1:2],minYr=minCpueYr )
      }
    )
  }
  choice
}


mcmcMenu <- function()
{
  menuItems <- c( "Plot biomass and projections by policy",
                  "Probability of projection biomass > reference",
                  "Expectation of projection biomass / reference",
                  "Plot biomass posterior densities (plotDens)",
                  "Plot recruitment posterior densities (plotDens)",
                  "Plot parameter posterior densities (plotDens)",
                  "Plot cumulative quantiles (plotCumu)",
                  "Plot traces (plotTrace)",
                  "Plot PJS traces (plt.allTraces)" )
#plt.ssbVbCatch( junk )
#windows()
#plt.expRate( junk )
#plt.numR( junk[1:2] )

  choice <- 1
  while ( choice > 0 )
  {
    cat( "\nMCMC Results (0 to exit):\n" )
    choice <- menu( menuItems )
    switch( choice,
#      {
#        currentMCMC <- importMCMC( dir=".", quiet=FALSE )
#        assign( "currentMCMC", currentMCMC, pos=1 )
#      },
#      {
#        currentProj <- importProj( dir=".", quiet=FALSE )
#        assign( "currentProj", currentProj, pos=1 )
#      },
      {
        plt.quantBio( currentMCMC$B,currentProj$B, policy=policy,
                      userPrompt=FALSE, xyType=rpType )
      },
      {
        currentProbs <- calc.projProbs2( currentMCMC$B, currentProj$B, refs )
        assign( "currentProbs",currentProbs, pos=1 )
        out.pmTables( currentProbs,fileName="pmProbs" )
        # AME: This works:
        # test = calc.projProbs( currentMCMC$B, currentProj$B,
        #        refYrs=c(1972,1985) )
      },
      {
#        calc.projExpect( currentMCMC$B, currentProj$B,
#                         refYrs=c(1972,1985) )

        currentExpt <- calc.projExpect2( currentMCMC$B, currentProj$B, refs )
        assign( "currentExpt",currentExpt, pos=1 )
        out.pmTables( currentExpt,fileName="pmExpt",dec=2 )
        # AME: after runSweaveMCMC, this works:   ...
        # currentExpt <- calc.projExpect( currentMCMC$B, currentProj$B, c(1940, 1990) )
      },
      {
        plotDens( currentMCMC$B[,getYrIdx(
          names(currentMCMC$B))],
          tick.number=5,
          xlab="Biomass", xlim=c(range(currentMCMC$B)),
          ylab="Posterior density" )
      },
      {
       #  plotDens( currentMCMC$R[,getYrIdx(
       #    names(currentMCMC$R))],
       #    xlab="Recruitment", xlim=c(range(currentMCMC$R)) )
       # AME replacing with:
       plotDensPOP( currentMCMC$R[,getYrIdx(names(currentMCMC$R))] ,          xlab="Recruitment", xlim=c(0, 150000),
         between = list(x=0.2, y=0.2), ylab="Density",
         lwd.density=2, same.limits=TRUE)
       # AME: between is gaps between panels,
       #  same.limits makes axes same for all panels
      },
      {
        idx <- apply( currentMCMC$P,2,allEqual )
        plotDens( currentMCMC$P[,!idx] )
      }, 
      {
        if ( length( getYrIdx(names(currentMCMC$B)) ) <= 9 )
          mfRow <- c(3,3)
        else
          #ACH 9/20/07: I'm changing the mfRow becasue the canary model has 68 years
          mfRow <- c(5,3)
        windows(record=TRUE)
        par( oma=c(2,2,1,1), mar=c(2,2,2,1), mfrow=mfRow )
        plotCumu( currentMCMC$B[,getYrIdx(names(currentMCMC$B))],
          auto.layout=FALSE, xlab="", ylab="", main="" )
        mtext( side=1, line=0.5, cex=1.0, outer=TRUE, "Iteration" )
        mtext( side=2, line=1, cex=1.0, outer=TRUE, "Biomass" )

        windows()
        par(mfrow=c(3,4), oma=c(2,2,1,1), mar=c(2,2,2,1) )
        idx <- apply( currentMCMC$P,2,allEqual )
        #ACH: Added the next few lines as a quick fix. Awatea now only outputs estimated parameters, thus the allEqual stuff is not necessary
        #plotCumu( (currentMCMC$P[,!idx])[,1:12] , auto.layout=FALSE)
        npars <- length(idx)
        plotCumu( (currentMCMC$P[,!idx])[,1:(min(npars,12))] , auto.layout=FALSE)
        mtext( side=1, line=0.5, cex=1.0, outer=TRUE, "Iteration" )
        mtext( side=2, line=1, cex=1.0, outer=TRUE, "Value" )
        if(npars > 12) {
            windows()
            par(mfrow=c(3,4), oma=c(2,2,1,1), mar=c(2,2,2,1) )
            plotCumu( (currentMCMC$P[,!idx])[,13:22], auto.layout=FALSE )
            mtext( side=1, line=0.5, cex=1.0, outer=TRUE, "Iteration" )
            mtext( side=2, line=1, cex=1.0, outer=TRUE, "Value" )
        }
      },
      {
        plotTrace( currentMCMC$R[,getYrIdx(
          names(currentMCMC$R))],axes=TRUE,
          xlab="Recruitment" )
        windows()
        par( oma=c(2,2,1,1), mar=c(2,2,2,1), mfrow=mfRow )
        plotTrace( currentMCMC$B[,getYrIdx(
          names(currentMCMC$B))],axes=TRUE,
          xlab="Biomass" )
        windows()
        idx <- apply( currentMCMC$P,2,allEqual )
        plotTrace( currentMCMC$P[,!idx], axes=TRUE )
      },
      {
        plt.allTraces( currentMCMC )
      }
    )
  }
  choice
}


utilMenu <- function()
{
  menuItems <- c( "scape Help",
                  "plotMCMC Help",
                  "Portrait graphsheet",
                  "Landscape graphsheet" )

  choice <- 1
  while ( choice > 0 )
  {
    cat( "\nUtilities:\n" )
    choice <- menu( menuItems )
    switch( choice,
      print( help( "scape", help_type="html" ) ),
      print( help( "plotMCMC", help_type="html" ) ),
      graphics( view="portrait" ),
      graphics( view="landscape" )
    )
  }
  choice
}
