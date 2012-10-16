#PBSscape-------------------------------2012-08-21
#  Modified functions from Arni Magnussen's 
#  packages 'scape' and 'scapeMCMC'
#-------------------------------------------AME/RH

#---History---------------------------------------
# ymrScape.r - for ymr, just the function definitions here. Calls to
#  them will go into .Snw. Call this from .Snw. 23rd February 2011
# SEE popScapeRuns2.r for figures to do for both runs at once, as
#  went into the POP SAR. Wait until we have final results.

# popScape2.r - going to do all postscript files as better than .png.
#  Getting rid of menus and automatically doing all required figures
#   each time I run it, then automatically into latex to check, then
#   when all good put into the Word file with captions etc. Will be
#   way more efficient for me to work, rather than all the time it
#   takes opening and closing windows.
# **Jump to end to load in MCMC and plot pairs plots and acf.
#
# popScape.r notes:
# Always choose one of the "Save all ** plots to PNG" options, as I'm not
#  repeating the new figures. 
# Go through figures, see what else to add. Some will just need tweaking - see
#  notes on Appendix D printout.
# Checking that first value for MCMC chain is what's in currentRes as
#  the "MPD". Do get:
# > currentRes$B[30, ]
#   Year     VB    SB       Y         U        R
#30 1969 145424 88107 10382.4 0.0713941 17866.36
#> currentMCMC$B[1, "1969"]
#[1] 88107    # so this is correct (they call it SB then B)
# But recruitment isn't - damn, maybe it's a year off again, like what I
#  figured out for MCMC:
# > currentMCMC$R[1, "1969"]
# [1] 12773.1        # so disagrees
# > currentMCMC$R[1, "1970"]
# [1] 17866.4
#                    # so it is a year off.
# But recruits aren't actually plotted from MPD stuff, so think is okay.

# Plots like biomass.png use the MPD, but shouldn't they really use the median
#  of posterior (or mean), and maybe show credible intervals?? Check difference
#  between median of Bt and MPD.  Yes, they should, so doing that now


# popScape.r - developing further in assess7/, from version in
#  assess6/.  Had to use postscript for the recruitment marginal
#  posterior densities due to resizing. They will look fine when
#  the .doc is converted to .pdf.
#  Using the lattice approaches, not the three ....PJS.jpg sections
#   which are not lattice, and more fiddly to change. 18-Oct-10
# popScape.r - renaming carScapeAndy.r, and continuing to develop. Haven't
#  kept track of the many changes in the list below - could always do later
#  using CSDiff if really necessary. 1-Sep-10
# carScapeAndy.r - AME editing so that menu works for POP.
#  Search for AME for other changes. 12-Aug-10
# Added some plots: residuals by age w/o outliers, and residuals by age w/o outliers. May only work for save all MPD plots as doing quickly.   27th August 2010
# AME: all three types of plots of residuals of catch-at-age data,
#  with x-axis being year, age or year of birth, fill in NA's to give
#  white space for years/ages in which there are no data. 30th August
#  2010
# goto HERE - ask Michael for Trellis book.
# AME changing all wmf -> png
# For Paul: I fixed the problem with the loading in of .res files,
#  that menu now works. So first load in the .res file from the
#  main menu (I think this
#  could be simplified as it does automatically load in all .res
#  files in the directory, but for now we'll stick with this).
# For the MPD plots (option 2 on main menu): options 1-3 work.
#  For option 1 (Plot biomass, recruitment, catch),
#  we need to discuss whether to use the plotB or plotB2 function.
#  Currently plotB2 is used, for which below it says:
#  "function to accommodate PJS request not to show biomass
#  prior to fishery and survey indices period." Consequently the 
#  biomass estimates are only show for the later years.
#  option 2 just plots initial age structure and recruitments,
#  but by changing the options it can show more ages, so we can
#  maybe play with that.
#  option 4 (Plot commercial catch-at-age results) does not work, 
#  but the code mentions a bug that needs sorting.

#  option 5 now Plot survey catch-at-age results, following have
#   thus been incremented by 1.
#  option 6 (Plot abundance index) does not work, but the similar
#  command:
#    plotIndex(currentRes, what="s", xlim=c(1960, 2012))
#  does work for the survey indices, so I could just change the
#  the options to select something similar if you like. 
#  option 7 (All residual plots) doesn't give an error now.
#  options 8-10 plot results from each .res file in the directory
#  so that they can be easily compared - looks useful.
# 
#  
#--------------------------------------------------------------------#
# rsScape.r : Graphical analyses for rocksole Coleraine ouput.       #
# Developers: A.R. Kronlund, P.J. Starr                              #
# Modified by Allan C. Hicks to use output of Awatea forSPO7         #
# Required libraries: gdata, scape, scapeMCMC, PBSmodelling          #
#                                                                    #
# Date Revised:                                                      #
# 08-Nov-05  Initial implementation.                                 #
# 10-Nov-05  Added generic residual plot plus stdRes calc functions. #
# 11-Nov-05  Added plot and save all functions.                      #
# 11-Nov-05  Added standardised age residuals a'la PJS doc.          #
# 12-Nov-05  Added res file selection function.                      #
# 12-Nov-05  Added MCMC load and selected plotting.                  #
# 13-Nov-05  Added biomass and projection plots. Probability tables. #
# 14-Nov-05  PJS adds a very small number of minor edits.            #
# 15-Nov-05  Addressed PJS improvements list.                        #
# 15-Nov-05  Added policy list at end of file so you don't have to   #
#            go hunting for it in the function code.                 #
# 15-Nov-05  Added console prompt to projection plots to step thru   #
#            multiple pages of plots depending on number of policies.#
# 17-Nov-05  Fixed age residuals, fixed multple fishery CA plots.    #
#            Fixed importCol to accommodate multiple gears.          #
#            Replace "importCol" with "importCol2" below.            #
# 21-Nov-05  Fixed bug in saving projection plots, Arni's scales to  #
#            trace plots.                                            #
# 21-Nov-05  Added quantile box plots for reconstruction-projections.#
# 22-Nov-05  Added observation-based reference points funciton.      #
# 23-Nov-05  Added import funs for PJS Delay Difference model output.#
# 24-Nov-05  Fixed bug in saved policyProjection plots.              #
# 25-Nov-05  Added new biomass trajectory plot to plot biomass only  #
#            over the period of tuning data.                         #
# 09-Dec-05  Added multi-panel plots for SSB/VB.                     #
# 09-Dec-05  Add PJS trace+mpd plots to rsScape.r.                   #
#            Correct bad recruits histogram and plot age 1's.        #
# 13-Dec-05  Finished PJS trace plots.                               #
# 16-Dec-05  Revised plots to not include the pre-1966 period.       #
#                                                                    #
# 16-Apr-06  Modified importCol2 to take Awatea output (ACH)         #
#            Added 0.5 to Year when plotting indices                 #
# 05-Sep-07  Minor modifications to plot age instead of length by    #                                    PJS                                       #
# 12-Aug-10  Commenting out mainMenu() to run from a script. AME.    #
# **-Aug-10  Many further changes specific to POP assessment. AME    #
# **-Sep-10   "                                                      #
# NOTE:                                                              #
# (1) For most plots there must be a "res" file loaded.              #
# (2) For MCMC plots you must have MCMC output loaded (*.pst files). #
# (3) For Projection plots and tables you must have MCMC and also    #
#     Projection output loaded.  See the MCMC menu.                  #
#                                                                    #
# To do:                                                             #
#                                                                    #
# (1) Check with Paul regarding calculations for age residuals.      #
#--------------------------------------------------------------------#
#                                                                    #
# Awatea res file has following structure (some elements may be      #
# missing dependent on model configuration and importCol details).   #
#                                                                    #
# N predicted numbers at age                                         #
# B predicted biomass, recruitment, and observed landings            #
# Sel predicted selectivity and observed maturity (age things)       #
# Dev predicted recruitment deviates from the stock-recruitment curve#
# CPUE, Survey commercial and survey abundance index and fit         #
# CAc, CAs commercial and survey C@A (catch at age) and fit          #
# CLc, CLs commercial and survey C@L (catch at length) and fit       #
# LA observed L@A and fit                                            #
#                                                                    #
# MCMC                                                               #
#                                                                    #
# The importProj function loads a list with elements "B" and "Y" for #
# biomass by catch policy and year, and catch by harvest policy and  #
# year.  The "B" element is itself a list of matrices with a matrix  #
# for each level of the catch.  This matrix has rows equal to the    #
# length of the chain and columns corresponding to projection years. #
# There are no plotting routines for these data.                     #
#--------------------------------------------------------------------#

#--------------------------------------------------------------------#
#                    Awatea Related Functions                     #
#--------------------------------------------------------------------#

# BUG FIX: Original importCol appears not to recognize multiple
#          fishery series.

#rm(list=ls())
#require("PBSmodelling")  # Doing it here before scapeMCMC, as both
                         #  have plotDens(). Though if already loaded
                         #  scapeMCMC that won't always work, so
                         #  replacing plotDens with
                         #  scapeMCMC::plotDens throughout code
                         # scapeMCMC::plotTrace() also
#source("plotDensPOP.r")  # AME's version that modifies defaults
                         #  to give better recruitment and Bt plots
#source("plotDensPOPpars.r")  # AME's version to add MPD for params
#source("plotTracePOP.r") # AME's version that adds on MPD
# source("plotBVBnorm.r") # B/B0 and V/V0, as lattice so can use for
                        #  multiple runs. NOT as lattice now.

#importCol2---------------------------- DEPRECATED
# Modified 'importCol' with RH's 'extra' section.
# Use function 'importRes' which is maintained.
#----------------------------------------------AME
importCol2 <- function (res.file, info = "", Dev = FALSE, CPUE = FALSE, Survey = FALSE,
    CAc = FALSE, CAs = FALSE, CLc = FALSE, CLs = FALSE, LA = FALSE, quiet = TRUE, extra = TRUE)
{
    readVector <- function(keyword, same.line = TRUE, file = res.file,
        vector = res.vector) {
        line <- match(keyword, substring(vector, 1, nchar(keyword)))
        v <- if (same.line)
            as.numeric(scan(file, what = "", skip = line - 1,
                nlines = 1, quiet = TRUE)[-1])
        else as.numeric(scan(file, what = "", skip = line, nlines = 1,
            quiet = TRUE))
        if (!quiet)
            cat("vector...")
        return(v)
    }
    readMatrix <- function(keyword, nrow, header = FALSE, stripe = c("no",
        "left", "right", "upper", "lower"), file = res.file,
        vector = res.vector) {
        stripe <- match.arg(stripe)
        line <- match(keyword, substring(vector, 1, nchar(keyword))) +
            as.numeric(header)
        m <- scan(file, skip = line, nlines = nrow, quiet = TRUE)
        m <- matrix(m, byrow = TRUE, nrow = nrow)
        m <- switch(stripe, left = m[, seq(1, ncol(m)/2)], right = m[,
            seq(ncol(m)/2 + 1, ncol(m))], upper = m[seq(1, nrow(m) -
            1, by = 2), ], lower = m[seq(2, nrow(m), by = 2),
            ], m)
        if (!quiet)
            cat("matrix...")
        return(m)
    }
    getN <- function(sexes, years, ages) {
        if (!quiet)
            cat("N         ")
        nsexes <- length(sexes)
        nyears <- length(years)
        nages <- length(ages)
        if (nsexes == 1) {
            Nu <- readMatrix("Numbers_at_age_by_Year,sex_and_age",
                nrow = nyears * nsexes)
            N <- data.frame(Sex = rep(sexes, nyears * nages),
                Year = rep(years, each = nages), Age = rep(ages,
                  nyears), N = as.vector(t(Nu)))
        }
        if (nsexes == 2) {
            Nf <- readMatrix("Numbers_at_age_by_Year,sex_and_age",
                nrow = nyears * nsexes, stripe = "upper")
            Nm <- readMatrix("Numbers_at_age_by_Year,sex_and_age",
                nrow = nyears * nsexes, stripe = "lower")
            N <- data.frame(Sex = rep(sexes, each = nyears *
                nages), Year = rep(rep(years, each = nages),
                2), Age = rep(ages, 2 * nyears), N = as.vector(t(rbind(Nf,
                Nm))))
        }
        if (!quiet)
            cat("OK\n")
        return(N)
    }
    getB <- function(years, gears) {
        ngears <- length(gears)
        if (!quiet)
            cat("B         ")
        vb <- readMatrix("Vulnerable_Biomass_by_Method_and_Year",
            nrow = ngears)
        sb <- readVector("Spawning_Biomass_by_Year", same.line = FALSE)

        # *** ADD in the exploitation rate, last year is missing and need
        #     to pad the matrix with missing values to match other "B" columns.

        U <- readMatrix("Exploitation_Rate_by_Method_and_Year", nrow=ngears )
        U <- cbind( U, rep(NA,nrow(U)) )

#        y <- c(readVector("Total_Catch_by_Method_and_Year", same.line = FALSE),
#            NA)

        # BUG FIX: Appears should call readMatrix to accommodate multiple gear series.
        #          Then, sum over the gears to get total catch.
        y <- readMatrix( "Total_Catch_by_Method_and_Year", nrow=ngears )
        y <- apply( y,2,sum,na.rm=TRUE )
        y <- c( y,NA )

        B <- as.data.frame( cbind(years, t(vb), sb, y, t(U)) )
        names(B) <- if (ngears == 1)
            c("Year", "VB", "SB", "Y", "U")
        else c("Year", paste("VB", gears, sep = "."), "SB", "Y",
                       paste("U",gears,sep=".") )
        if (!quiet)
            cat("OK\n")
        return(B)
    }
    getSel <- function(gears, surveys, years, sexes, ages) {
        if (!quiet)
            cat("Sel       ")
        ngears <- length(gears)
        nsurveys <- length(surveys)
        nyears <- length(years)
        nsexes <- length(sexes)
        nages <- length(ages)
        com <- readMatrix("Commercial_age-specific_selectivity_by_method,Year,sex_and_age",
            nrow = ngears * nyears * nsexes)
        com <- com[seq(1, to = ngears * nyears * nsexes, by = nyears),
            ]
        srv <- readMatrix("Survey_age-specific_selectivity_by_survey,Year,sex_and_age",
            nrow = nsurveys * nsexes)
        fecundity <- readVector("Fecundity_by_year_and_age",
            same.line = FALSE)
        weight <- readVector("Weight_by_year,sex_and_age", same.line = FALSE)
        mat <- rep(ifelse(weight > 0, fecundity/weight, 0), nsexes)
        if (is.numeric(gears))
            gears <- paste("Gear", gears)
        if (is.numeric(surveys))
            surveys <- paste("Survey", surveys)
        Sel <- data.frame(Series = c(rep(gears, each = nsexes *
            nages), rep(surveys, each = nsexes * nages), rep("Maturity",
            nsexes * nages)), Sex = rep(rep(sexes, each = nages),
            ngears + nsurveys + 1), Age = rep(ages, (ngears +
            nsurveys + 1) * nsexes), P = c(t(com), t(srv), mat))
        if (!quiet)
            cat("OK\n")
        return(Sel)
    }
    getDev <- function(ages, years) {
        if (!quiet)
            cat("Dev       ")
        Dev <- list()
        Dev$Initial <- readVector("log_InitialDev", same.line = TRUE)
        names(Dev$Initial) <- ages[-c(1, length(ages))]
        Dev$Annual <- readVector("log_RecDev", same.line = TRUE)
        names(Dev$Annual) <- years[-length(years)]
        if (!quiet)
            cat("OK\n")
        return(Dev)
    }
    getCPUE <- function(gears, years) {
        if (!quiet)
            cat("CPUE      ")
        nseries <- readVector("NCPUEindex")
        ngears <- length(gears)
        nyears <- length(years)
        obs <- readMatrix("indexmethodyearvaluecv", nrow = readVector("Number_of_CPUE_data",
            same.line = FALSE))
        obs <- data.frame(Series = obs[, 1], Gear = obs[, 2],
            Year = obs[, 3], Obs = obs[, 4], CV = obs[, 5])
        fit <- readMatrix("CPUE_Index_Trajectories", nrow = nseries)
        fit <- data.frame(Series = rep(1:nseries, each = nyears),
            Year = rep(years, nseries), Fit = as.vector(t(fit)))
        CPUE <- merge(obs[, names(obs) != "Gear"], fit, all = TRUE)
        sgkey <- unique(obs[, c("Series", "Gear")])
        CPUE <- merge(sgkey, CPUE)
        CPUE <- data.frame(Series = paste("Series ", CPUE$Series,
            "-", CPUE$Gear, sep = ""), Year = as.integer(CPUE$Year),
            Obs = CPUE$Obs, CV = CPUE$CV, Fit = CPUE$Fit)
        if (!quiet)
            cat("OK\n")
        return(CPUE)
    }
    getSurvey <- function(years) {
        if (!quiet)
            cat("Survey    ")
        nyears <- length(years)
        nseries <- readVector("Nsurveyindex")
        obs <- readMatrix("indexyearvaluecv", nrow = readVector("Number_of_survey_data",
            same.line = FALSE))
        obs <- data.frame(Series = obs[, 1], Year = obs[, 2],
            Obs = obs[, 3], CV = obs[, 4])
        fit <- readMatrix("Survey_Index_Trajectories", nrow = nseries)
        fit <- data.frame(Series = rep(1:nseries, each = nyears),
            Year = rep(years, nseries), Fit = as.vector(t(fit)))
        Survey <- merge(obs, fit, all = TRUE)
        Survey$Series <- as.integer(Survey$Series)
        Survey$Year <- as.integer(Survey$Year)
        if (!quiet)
            cat("OK\n")
        return(Survey)
    }
    getCAc <- function(sexes, ages) {
        if (!quiet)
            cat("CAc       ")
        nsexes <- length(sexes)
        nages <- length(ages)
        nobs <- readVector("Number_of_Commercial_C@A", same.line=FALSE)
        obs <- readMatrix("methodyearsamplesizesex1a1sex1a2sex1a3", nrow=nobs)                     # "Observed_C@A"  not unique
        fit <- readMatrix("methodyearsamplesizesex1a1sex1a2sex1a3", nrow=nobs, header=2*(nobs+1))  # "Predicted_C@A" not unique
        CAc <- data.frame(Series=rep(obs[,1],each=nsexes*nages), Year=rep(obs[,2],each=nsexes*nages),
                      SS=rep(obs[,3],each=nsexes*nages), Sex=rep(rep(sexes,each=nages),nobs),
                      startL=rep(obs[,4],each=nsexes*nages),endL=rep(obs[,5],each=nsexes*nages),
                      Age=rep(ages,nsexes*nobs), Obs=as.vector(t(as.matrix(obs[,-(1:5)]))),
                      Fit=as.vector(t(as.matrix(fit))))
         # loads in okay with next line commented. Gears not
         #  included in getCAs(), and we don't need for POP. AME
         #   CAc$Gear <- as.integer(CAc$Gear)
        CAc$Year <- as.integer(CAc$Year)
        CAc$Age <- as.integer(CAc$Age)
        if (!quiet)
            cat("OK\n")
        return(CAc)
    }
    getCAs <- function(sexes, ages) {
        if (!quiet)
            cat("CAs       ")
        nsexes <- length(sexes)
        nages <- length(ages)
        nobs <- readVector("Number_of_survey_C@A", same.line = FALSE)
        obs <- readMatrix("surveyyearsamplesizesex1a1sex1a2sex1a3",
            nrow = nobs)
        fit <- readMatrix("surveyyearsamplesizesex1a1sex1a2sex1a3",
            nrow = nobs, header = 2 * (nobs + 1))
        CAs <- data.frame(Series=rep(obs[,1],each=nsexes*nages), Year=rep(obs[,2],each=nsexes*nages),
                      SS=rep(obs[,3],each=nsexes*nages), Sex=rep(rep(sexes,each=nages),nobs),
                      startL=rep(obs[,4],each=nsexes*nages),endL=rep(obs[,5],each=nsexes*nages),
                      Age=rep(ages,nsexes*nobs), Obs=as.vector(t(as.matrix(obs[,-(1:5)]))),
                      Fit=as.vector(t(as.matrix(fit))))
        CAs$Series <- as.integer(CAs$Series)
        CAs$Year <- as.integer(CAs$Year)
        CAs$Age <- as.integer(CAs$Age)
        if (!quiet)
            cat("OK\n")
        return(CAs)
    }
    getCLc <- function(sexes, lengths) {
        if (!quiet)
            cat("CLc       ")
        nsexes <- length(sexes)
        nlengths <- length(lengths)
        nobs <- readVector("Number_of_Commercial_C@L", same.line=FALSE)
        obs <- readMatrix("methodyearsamplesizesex1l1sex1l2sex1l3", nrow=nobs)                 # "Observed_C@L"  not unique
        fit <- readMatrix("methodyearsamplesizesex1l1sex1l2sex1l3", nrow=nobs, header=nobs+1)  # "Predicted_C@L" not unique
    
        CLc <- data.frame(Series=rep(obs[,1],each=nsexes*nlengths), Year=rep(obs[,2],each=nsexes*nlengths),
                      SS=rep(obs[,3],each=nsexes*nlengths), Sex=rep(rep(sexes,each=nlengths),nobs),
                      startL=rep(obs[,4],each=nsexes*nlengths),endL=rep(obs[,5],each=nsexes*nlengths),
                      Length=rep(lengths,nsexes*nobs), Obs=as.vector(t(as.matrix(obs[,-(1:5)]))),
                      Fit=as.vector(t(as.matrix(fit))))
        CLc$Series <- as.integer(CLc$Series)
        CLc$Year <- as.integer(CLc$Year)
        CLc$Length <- as.integer(CLc$Length)
        if (!quiet)
            cat("OK\n")
        return(CLc)
    }
    getCLs <- function(sexes, lengths) {
        if (!quiet)
            cat("CLs       ")
        nsexes <- length(sexes)
        nlengths <- length(lengths)
        nobs <- readVector("Number_of_surveyC@L",same.line=FALSE)
        obs <- readMatrix("surveyyearsamplesizesex1l1sex1l2sex1l3", nrow=nobs)                     # "Observed_C@L"  not unique
        fit <- readMatrix("surveyyearsamplesizesex1l1sex1l2sex1l3", nrow=nobs, header=2*(nobs+1))  # "Predicted_C@L" not unique
        CLs <- data.frame(Series=rep(obs[,1],each=nsexes*nlengths), Year=rep(obs[,2],each=nsexes*nlengths),
                      SS=rep(obs[,3],each=nsexes*nlengths), Sex=rep(rep(sexes,each=nlengths),nobs),
                      startL=rep(obs[,4],each=nsexes*nlengths),endL=rep(obs[,5],each=nsexes*nlengths),
                      Length=rep(lengths,nsexes*nobs), Obs=as.vector(t(as.matrix(obs[,-(1:5)]))),
                      Fit=as.vector(t(as.matrix(fit))))
        CLs$Series <- as.integer(CLs$Series)
        CLs$Year <- as.integer(CLs$Year)
        CLs$Length <- as.integer(CLs$Length)
        if (!quiet)
            cat("OK\n")
        return(CLs)
    }
    getLA <- function(sexes, ages) {
        if (!quiet)
            cat("LA        ")
        nsexes <- length(sexes)
        nages <- length(ages)
        nobs <- readVector("#femalesmales", same.line = FALSE,
            file = latage.file, vector = latage.vector)
        obs <- readMatrix("VonBertalanfy--Lenght-at-agefit--Likelihood",
            nrow = sum(nobs), header = 8)
        obs <- data.frame(Sex = rep(sexes, nobs), Age = obs[,
            1], Obs = obs[, 2])
        owarn <- options(warn.conflicts = -1)
        Linf <- readVector("VonBeratalanfy:Linf")[-(1:3)]
        K <- readVector("VonBeratalanfy:k")[-(1:3)]
        t0 <- readVector("VonBeratalanfy:to")[-(1:3)]
        CV1 <- readVector("cvoftheFitbysex")[-(1:5)]
        CVratio <- readVector("ratioofcv(L_an)/cv(L_a1)oftheFitbysex")[-(1:7)]
        options(owarn)
        sigmaLA <- readVector("#LinearrelationshipofsigmaL@A(1=age;2=length---ignoreifW@Aissupplied)",
            same.line = FALSE, file = txt.file, vector = txt.vector)[1]
        max.age <- c(max(obs$Age[obs$Sex == sexes[1]]), max(obs$Age[obs$Sex ==
            sexes[2]]))
        fit <- data.frame(Sex = rep(sexes, max.age), Age = c(1:max.age[1],
            1:max.age[2]))
        fit$Fit[fit$Sex == sexes[1]] <- Linf[1] * (1 - exp(-K[1] *
            (fit$Age[fit$Sex == sexes[1]] - t0[1])))
        fit$Fit[fit$Sex == sexes[2]] <- Linf[2] * (1 - exp(-K[2] *
            (fit$Age[fit$Sex == sexes[2]] - t0[2])))
        if (sigmaLA == 1) {
            A <- rep(max(ages), 2)
            a <- cbind(fit$Age[fit$Sex == sexes[1]], fit$Age[fit$Sex ==
                sexes[2]])
            fit$CV[fit$Sex == sexes[1]] <- CV1[1] + CV1[1] *
                (CVratio[1] - 1)/(A[1] - 1) * (a[, 1] - 1)
            fit$CV[fit$Sex == sexes[2]] <- CV1[2] + CV1[2] *
                (CVratio[2] - 1)/(A[2] - 1) * (a[, 2] - 1)
        }
        if (sigmaLA == 2) {
            L1 <- Linf * (1 - exp(-K * (1 - t0)))
            Ln <- Linf * (1 - exp(-K * (max(ages) - t0)))
            fit$CV[fit$Sex == sexes[1]] <- CV1[1] + CV1[1] *
                (CVratio[1] - 1)/(Ln[1] - L1[1]) * (fit$Fit[fit$Sex ==
                sexes[1]] - L1[1])
            fit$CV[fit$Sex == sexes[2]] <- CV1[2] + CV1[2] *
                (CVratio[2] - 1)/(Ln[2] - L1[2]) * (fit$Fit[fit$Sex ==
                sexes[2]] - L1[2])
        }
        LA <- merge(obs, fit, by = c("Sex", "Age"), all = TRUE)
        LA$Age <- as.integer(LA$Age)
        LA$Fit <- LA$Fit
        LA$CV <- LA$CV
        if (!quiet)
            cat("OK\n")
        return(LA)
    }
    	getExtra = function(resvec) {
		extra = list(
			likelihoods = c(
			"CPUE","Survey_Index","C@A_Commercial","C@A_survey","Prior_penalties"),  # Likelihoods
			parameters = c(
			"R0","avgR0","h",                                                      # Parameters (P)
			"M1","M2",                                                             #  (P) Sex_specific
			"Sfullest","SfullDelta","log_varLest","log_varRest",                   #  (P) Method_specific
			"log_qCPUE","log_BetaCPUE",                                            #  (P) CPUE_index_specific
			"log_qsurvey","surveySfull","survey_SfullDelta",
			"log_surveyvarL","log_surveyvarR",                                     # (P)  Survey_index_specific
			"log_RecDev"),                                                         #  (P) Recruitment_residuals
			priors = c(
			"R0_prior","h_prior",                                                  # Priors (I)
			"M1_prior","M2_prior","Rinit_prior","uinit_prior","p_plusscale",       #  (I) Sex_specific
			"p_Sfullest","p_Sfulldelta","log_varLest_prior","log_varRest_prior",   #  (I) Method_specific
			"errSfull_prior","errvarL_prior","errvarR_prior",                      #  (I) Method_specific_and_annual
			"log_qCPUE_prior","log_BetaCPUE_prior",                                #  (I) CPUE_index_specific
			"qCPUEerr_prior",                                                      #  (I) CPUE_index_specific_and_annual
			"log_qsurvey_prior","surveySfull_prior","p_surveySfulldelta",
			"log_surveyvarL_prior","log_surveyvarR_prior"),                        #  (I) Survey_index_specific
			residuals = c(
			"p_log_InitialDev","p_log_RecDev")                                     # Recruitment_residuals
		)
		Nsexes = as.numeric(rev(strsplit(resvec[grep("^Nsexes",resvec)],split=" ")[[1]])[1])
		Nsurveyindex = as.numeric(rev(strsplit(resvec[grep("^Nsurveyindex",resvec)],split=" ")[[1]])[1])
		elist = sapply(extra,function(x,resvec) { 
			ex = as.list(x); names(ex) = x
			for (i in ex) {
				expr=paste("index = grep(\"^",i," \",resvec)",sep="")
				eval(parse(text=expr))
				if (length(index)==1) {
					if (i %in% c("M1_prior","M2_prior","uinit_prior","p_plusscale")) index = index + (1:Nsexes) - 1
					if (i %in% c("log_qsurvey_prior","surveySfull_prior","p_surveySfulldelta",
						"log_surveyvarL_prior","log_surveyvarR_prior")) index = index + (1:Nsurveyindex) - 1
					exvec = strsplit(resvec[index],split=" ")
					exmat = t(sapply(exvec,function(x){as.numeric(x[!is.element(x,c(i,""))])}))
					if (nrow(exmat)==1) exres=as.vector(exmat)
					else                exres=exmat
#if (i=="log_qsurvey_prior") {browser();return() }
					ex[[i]] = exres
				}
				else ex[[i]] = "not found"
			}
			return(ex) }, resvec=resvec, simplify=FALSE)

	}
	#---END SUBFUNCTIONS---------------------------
    if (!file.exists(res.file))
        stop("File ", res.file, " not found. Use / or \\\\ separators.")
    res.vector <- readLines(res.file)
    if (extra)
	extra = getExtra(res.vector)
       else extra = NULL

    res.vector <- gsub("\"", "", gsub("\t", "", gsub(" ", "",
        res.vector)))
    if (!quiet)
        cat("\nParsing text file ", res.file, ":\n\nPreamble  ",
            sep = "")
    sexes <- if (readVector("Nsexes") == 1)
        "Unisex"
    else c("Female", "Male")
    gears <- seq(1, length.out = readVector("Nmethods"))
    surveys <- seq(1, length.out = readVector("Nsurveyindex"))
    years <- seq(from = readVector("StartYear"), to = readVector("EndYear") +
        1)
    ages <- seq(from = 1, to = readVector("Nages"))
    lengths <- seq(from = readVector("First_length"), by = readVector("Length_class_increment"),
        length.out = readVector("Number_of_length_classes"))
    if (!quiet)
        cat("OK\n")
    model <- list()
    model$N <- getN(sexes, years, ages)
    model$B <- getB(years, gears)
    rec <- model$N[model$N$Age == 1, ]
    rec <- tapply(rec$N, rec$Year, sum)
    model$B$R <- c(rec[-1], NA)
    model$Sel <- getSel(gears, surveys, years, sexes, ages)
    if (Dev)
        model$Dev <- getDev(ages, years)
    if (CPUE)
        model$CPUE <- getCPUE(gears, years)
    if (Survey)
        model$Survey <- getSurvey(years)
    if (CAc)
        model$CAc <- getCAc(sexes, ages)
    if (CAs)
        model$CAs <- getCAs(sexes, ages)
    if (CLc)
        model$CLc <- getCLc(sexes, lengths)
    if (CLs)
        model$CLs <- getCLs(sexes, lengths)
    if (LA) {
        latage.file <- paste(dirname(res.file), "l_at_age.dat",
            sep = "/")
        if (!file.exists(latage.file))
            stop("File ", latage.file, " not found. Use / or \\\\ separators.")
        latage.vector <- readLines(latage.file)
        latage.vector <- gsub("\"", "", gsub("\t", "", gsub(" ",
            "", latage.vector)))
        txt.file <- gsub("\\.res", "\\.txt", res.file)
        if (!file.exists(txt.file))
            stop("File ", txt.file, " not found. Use / or \\\\ separators.")
        txt.vector <- readLines(txt.file)
        txt.vector <- gsub("\"", "", gsub("\t", "", gsub(" ",
            "", txt.vector)))
        model$LA <- getLA(sexes, ages)
    }
    model$extra = extra
    if (!quiet)
        cat("\n")
    attr(model, "call") <- match.call()
    attr(model, "scape.version") <- installed.packages()["scape",
        "Version"]
    attr(model, "info") <- info
    class(model) <- "scape"
    return(model) }
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^importCol2


load.allResFiles <- function( resList=NULL )  
{
  # Loads all Awatea "res" files in working directory into a list.
  if ( is.null(resList) )
  resList <- dir( path="../.", pattern="results.dat$" )
         #YMR AME path was "."
         # YMR pattern was ".res$", but need to go through Excel
         #  to get .res, and .dat has same numbers (see .Snw)
         # NOT called now by ymrrun1dos.Snw,just import results.dat

  result <- as.list( c(1:length(resList)) )
  names( result) <- resList

  for ( i in 1:length(result) )
    result[[i]] <- importCol2( res.file=resList[i], Dev=TRUE, CPUE=TRUE,
                               Survey=TRUE, CLc=TRUE, CLs=TRUE, CAs=TRUE, CAc=TRUE)
  result                             # AME added CAc=TRUE, CAs=TRUE.
}


#plotB2---------------------------------2011-08-31
# This is an alteration of Arni Magnussons "plotB" function to accommodate
# PJS request not to show biomass prior to fishery and survey indices period.
#----------------------------------------------AME
plotB2 <- function (model, what = "d", series = NULL, years = NULL, axes = TRUE,
    div = 1, legend = "bottom", main = "", xlab = "", ylab = "",
    cex.main = 1.2, cex.legend = 1, cex.lab = 1, cex.axis = 0.8,
    las = 1, tck = c(1, what == "d")/2, tick.number = 5, lty.grid = 3,
    col.grid = "white", pch = 16, cex.points = 0.8, col.points = "black",
    lty.lines = 1:3, lwd.lines = 2, col.lines = "black", ratio.bars = 3,
    col.bars = "grey", plot = TRUE, ...)
{
    panel.linebar <- function(x, y, bars, ...) {
        panel.abline(h = pretty(y, tick.number), lty = lty.grid,
            col = col.grid)
        panel.superpose(x, y, ...)
        panel.barchart(bars$Year, bars$Value, horizontal = FALSE,
            box.ratio = ratio.bars, col = col.bars)
    }
    panel.bar <- function(x, y, ...) {
        panel.abline(h = pretty(y, tick.number), lty = lty.grid,
            col = col.grid)
        panel.barchart(x, y, horizontal = FALSE, box.ratio = ratio.bars,
            col = col.bars)
    }
    if (class(model) != "scape")
        stop("The 'model' argument should be a scape object, not ",
            chartr(".", " ", class(model)), ".")
    what <- match.arg(what, c("d", "s", "l"))
    las <- as.numeric(las)
    x <- model$B
    x <- data.frame(Year = rep(x$Year, ncol(x) - 1), Series = rep(names(x)[-1],
        each = nrow(x)), Value = as.vector(as.matrix(x[, -1])))
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


    Bframe <- x[x$Series %in% grep("B", series, value = TRUE),]
    Bframe$Series <- factor(Bframe$Series)

    # Find the first year where there are fishery CPUE data.
    cpueYear1 <- min( model$CPUE$Year[ !is.na( model$CPUE$Obs ) ] )

    # Set all SB and VB values to NA for years less than cpueYear1.
    Bframe$Value[ Bframe$Year < cpueYear1 ] <- NA

    Rframe <- x[x$Series == "R", ]
    Yframe <- x[x$Series == "Y", ]

    Bframe$Value <- Bframe$Value/div[1]
    Rframe$Value <- Rframe$Value/rep(div, length.out = 2)[2]
    Yframe$Value <- Yframe$Value/div[1]

    require(grid, quietly = TRUE, warn.conflicts = FALSE)
    require(lattice, quietly = TRUE, warn.conflicts = FALSE)
    if (trellis.par.get()$background$col == "#909090") {
        for (d in dev.list()) dev.off()
        trellis.device(color = FALSE)
    }
    main <- rep(main, length.out = 2)
    xlab <- rep(xlab, length.out = 2)
    ylab <- rep(ylab, length.out = 2)
    las <- rep(las, length.out = 2)
    mymain <- list(label = main[1], cex = cex.main)
    myxlab <- list(label = xlab[1], cex = cex.lab)
    myylab <- list(label = ylab[1], cex = cex.lab)
    myrot <- switch(as.character(las[1]), "0" = list(x = list(rot = 0),
        y = list(rot = 90)), "1" = list(x = list(rot = 0), y = list(rot = 0)),
        "2" = list(x = list(rot = 90), y = list(rot = 0)), "3" = list(x = list(rot = 90),
            y = list(rot = 90)))
    myscales <- c(list(draw = axes, cex = cex.axis, tck = tck,
        tick.number = tick.number), myrot)
    lty.lines <- rep(lty.lines, length.out = nlevels(Bframe$Series))
    lwd.lines <- rep(lwd.lines, length.out = nlevels(Bframe$Series))
    col.lines <- rep(col.lines, length.out = nlevels(Bframe$Series))
    mykey <- list(space = legend, text = list(lab = levels(Bframe$Series),
        cex = cex.legend), lines = list(lty = lty.lines, lwd = lwd.lines,
        col = col.lines))
    if (what == "s") {
        graph <- xyplot(Rframe$Value ~ Bframe$Value[Bframe$Series ==
            "SB"], main = mymain, xlab = myxlab, ylab = myylab,
            scales = myscales, pch = pch, cex = cex.points, col = col.points,
            ...)
        graph$x.limits[1] <- 0
    }
    else if (what == "d" && nrow(Bframe) > 0) {
        graph <- xyplot(Value ~ Year, groups = Series, data = Bframe,
            panel = panel.linebar, type = "l", bars = Yframe,
            main = mymain, xlab = myxlab, ylab = myylab, scales = myscales,
            key = mykey, lty = lty.lines, lwd = lwd.lines, col = col.lines,
            ...)
    }
    else {
        graph <- xyplot(Value ~ Year, data = Yframe, panel = panel.bar,
            main = mymain, xlab = myxlab, ylab = myylab, scales = myscales,
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
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^plotB2


#plotBmcmcPOP---------------------------2011-08-31
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
plotBmcmcPOP = function(obj, currentRes1 = currentRes,
                  p=c(0.025,0.25,0.5,0.75,0.975),
                  xyType="quantBox",
                  lineType=c(3,2,1,2,3),
                  refLines=NULL, xLim=NULL, yLim=NULL,
                  userPrompt=FALSE, save=TRUE, xLab = c(1939, 1939, 1939),
                  yLab = c(10000, 70000, 170000),
                  textLab = c("catch", "spawning", "vulnerable"),
                  yaxis.by=10000, tcl.val=-0.2, ...)
                               # xLab - x position for label, etc.
  {
                  # See plt.quantBio if want other xyTypes, as took out here:
  plt.qB <- function( obj, xyType="lines", new=TRUE, xLim, yLim,... ) 
    {
    if ( new )
     plot( xLim,yLim, type="n", xlab="Year",ylab="Biomass or catch (t)" )

    yrs <- as.numeric(dimnames(obj)[[2]])

    # Quantile boxplots - assumes five quantiles.
    if ( xyType=="quantBox" )
    {
      delta <- 0.25
      # Draw the outer whiskers.
      segments( yrs,obj[1,], yrs,obj[5,], lty=1,col=1 )
      # Overlay the box.
      for ( i in 1:length(yrs) )
        rect( yrs[i]-delta,obj[2,i], yrs[i]+delta, obj[4,i],... )
      # Add the median.
      segments( yrs-delta,obj[3,],yrs+delta,obj[3,],lty=1,col=1 )
    }
  }

  # Plot quantiles of biomass using the posterior densities.

  yrs1 <- NULL
  yrs2 <- NULL
  result1 <- NULL
  result2 <- NULL

  # Calculate the quantiles of the reconstructed biomass.
  result1 <- apply( obj,2,quantile,probs=p )
  yrs1 <- as.numeric(dimnames(result1)[[2]])

  if ( is.null(yLim) )
    {
      yLim <- c(0, max(c(max(result1), max(currentRes1$B$VB)))) #range(result1)
      xLim=range(yrs1)
    }
  # xLegPos = xLeg*diff(xLim)+xLim[1]    # position of xLeg
  # yLegPos = yLeg*diff(yLim)+yLim[1]

  plt.qB( result1,xLim=xLim,yLim=yLim, xyType=xyType )
  points(obj$B$Year, currentRes1$B$Y, type="h", lwd=3)   # catch
  points(obj$B$Year, currentRes1$B$VB, type="p")         # vuln biom
  text( xLab, yLab, textLab, pos=4, offset=0)
  axis(1, at=xLim[1]:xLim[2], tcl=tcl.val, labels=FALSE)
  axis(2, at = seq(0, yLim[2], by=yaxis.by), tcl=tcl.val, labels=FALSE)

  # legend(xLegPos, yLegPos, c("Vulnerable", "Spawning", "Catch"), bty="n")
  # points(xLegPos-2, yLegPos, type="p")
    # mtext( side=1, line=2, cex=1.0, "Year" )
    # mtext( side=2, line=2, cex=1.0, "Biomass" )
  # }
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^plotBmcmcPOP


#plotVBcatch----------------------------2011-08-31
# AME adding, based on plotBmcmcPOP (tweaking some)
#  currentMCMC$B.  currentRes1 is local currentRes.
#----------------------------------------------AME
plotVBcatch = function(obj, currentRes1 = currentRes,
                  p=c(0.025,0.25,0.5,0.75,0.975),
                  xyType="quantBox",
                  lineType=c(3,2,1,2,3),
                  refLines=NULL, xLim=NULL, yLim=NULL,
                  userPrompt=FALSE, save=TRUE, xLab = c(1939, 1939),
                  yLab = c(10000, 220000),
                  textLab = c("catch", "vulnerable"),
                  yaxis.by=10000, tcl.val=-0.2, ...)
                               # xLab - x position for label, etc.
  {
                  # See plt.quantBio if want other xyTypes, as took out here:
  plt.qB <- function( obj, xyType="lines", new=TRUE, xLim, yLim,... ) 
    {
    if ( new )
     plot( xLim,yLim, type="n", xlab="Year",ylab="Catch and vulnerable biomass (t)" )
    
    yrs <- as.numeric(dimnames(obj)[[2]])

    # Quantile boxplots - assumes five quantiles.
    if ( xyType=="quantBox" )
    {
      delta <- 0.25
      # Draw the outer whiskers.
      segments( yrs,obj[1,], yrs,obj[5,], lty=1,col=1 )
      # Overlay the box.
      for ( i in 1:length(yrs) )
        rect( yrs[i]-delta,obj[2,i], yrs[i]+delta, obj[4,i],... )
      # Add the median.
      segments( yrs-delta,obj[3,],yrs+delta,obj[3,],lty=1,col=1 )
    }
  }

  # Plot quantiles of biomass using the posterior densities.

  yrs1 <- NULL
  yrs2 <- NULL
  result1 <- NULL
  result2 <- NULL

  # Calculate the quantiles of the reconstructed biomass.
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

  # xLegPos = xLeg*diff(xLim)+xLim[1]    # position of xLeg
  # yLegPos = yLeg*diff(yLim)+yLim[1]

  plt.qB( result1,xLim=xLim,yLim=yLim, xyType=xyType )
  points(currentRes1$B$Year, currentRes1$B$Y, type="h", lwd=3)   # catch

  # points(obj$B$Year, currentRes1$B$VB, type="p")
                          # was vuln biom MPD
  # text( xLab, yLab, textLab, pos=4, offset=0)   # Taking out
  #   as not really needed if give a decent caption
  axis(1, at=xLim[1]:xLim[2], tcl=tcl.val, labels=FALSE)
  axis(2, at = seq(0, yLim[2], by=yaxis.by), tcl=tcl.val, labels=FALSE)

  # legend(xLegPos, yLegPos, c("Vulnerable", "Spawning", "Catch"), bty="n")
  # points(xLegPos-2, yLegPos, type="p")
    # mtext( side=1, line=2, cex=1.0, "Year" )
    # mtext( side=2, line=2, cex=1.0, "Biomass" )
  # }
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^plotVBcatch


#plotBVBnorm----------------------------2011-08-31
# AME doing, tried in separate file, but then changed that to
#  lattice and wouldn't be good format for Arni's boxplots.
#  Based on plotVBcatch (tweaking some)
#  currentMCMC$B.  currentRes1 is local currentRes.
#  xLab - x position for label, etc.
#----------------------------------------------AME
plotBVBnorm = function(mcmcObj,
                   p=c(0.025,0.25,0.5,0.75,0.975),
                   xyType="quantBox",
                   lineType=c(3,2,1,2,3),
                   refLines=NULL, xLim=NULL, yLim=NULL,
                   userPrompt=FALSE, save=TRUE, xLeg = 0.7, yLeg=0.9,
                   yaxis.by=0.02, tcl.val=-0.2,
                   B.col="black", VB.col="black", ...)
                                # xLab - x position for label, etc.
   {
   # Calculate medians to be plotted
   BoverB0 = currentMCMC$B / currentMCMC$B[,1]     # B/B0  each chain
   VBoverVB0 = currentMCMC$VB / currentMCMC$VB[,1] #VB/VB0 each chain
   
   BoverB0med = apply(BoverB0, 2, median)         # median each year
   VBoverVB0med = apply(VBoverVB0, 2, median)     # median each year
 
 
   # Plot quantiles of biomass using the posterior densities.
 
   yrs1 <- NULL
   yrs2 <- NULL
   result1 <- NULL
   result2 <- NULL
 
   result1 = BoverB0med
   result2 = VBoverVB0med
   yrs1 <- as.numeric(names(result1))    # dimnames(result1)[[2]])
 
   if ( is.null(yLim) )
     {
       yLim <- c(0, max(c(max(result1), max(result2))))
     }
   if ( is.null(xLim) )
     {
       xLim=range(yrs1)
     }
 
   # xLegPos = xLeg*diff(xLim)+xLim[1]    # position of xLeg
   # yLegPos = yLeg*diff(yLim)+yLim[1]
 
   plot( xLim,yLim, type="n", xlab="Year",ylab="Biomasses relative to virgin")
   points(yrs1, result1, type="p", col=B.col) 
   points(yrs1, result2, type="l", col=VB.col)
   # text( xLab, yLab, textLab, pos=4, offset=0)
   axis(1, at=xLim[1]:xLim[2], tcl=tcl.val, labels=FALSE)
   axis(2, at = seq(0, yLim[2], by=yaxis.by), tcl=tcl.val, labels=FALSE)
 
   xLegPos = xLeg*diff(xLim)+xLim[1]    # position of xLeg
   yLegPos = yLeg*diff(yLim)+yLim[1]
 
   legend(xLegPos, yLegPos, c(expression(B[t]/B[0]), expression(V[t]/V[0])), bty="n", lty=c(0,1), pch=c(1, NA), col=c(B.col, VB.col))   # pch=c(1, NA) okay, but adds o
   # points(xLegPos-2, yLegPos, type="p")
     # mtext( side=1, line=2, cex=1.0, "Year" )
     # mtext( side=2, line=2, cex=1.0, "Biomass" )
   # }
 }
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^plotBVBnorm


#plotRmcmcPOP---------------------------2011-08-31
# AME adding, plotting recruitment posteriors quantiles as one graph over time.
#  Already have the full posterior densities done.
#  Using plotBmcmcPOP as template, but will be simpler as no extra stuff. Prob
#   not simplifying down as much as could due to time constraints.
# Adding yLab and then using for exploitation plot also
#----------------------------------------------AME
plotRmcmcPOP = function(obj, 
                  p=c(0.025,0.25,0.5,0.75,0.975),
                  xyType="quantBox",
                  lineType=c(3,2,1,2,3),
                  refLines=NULL, xLim=NULL, yLim=NULL,
                  userPrompt=FALSE, save=TRUE, tcl.val=-0.2,
                  yaxis.by=10000, yLab="Recruitment", ...)
  {
                  # See plt.quantBio if want other xyTypes, as took out here:
  plt.qB <- function( obj, xyType="lines", new=TRUE, xLim, yLim,... ) 
    {
    if ( new )
     plot( xLim,yLim, type="n", xlab="Year",ylab=yLab )

    yrs <- as.numeric(dimnames(obj)[[2]])

    # Quantile boxplots - assumes five quantiles.
    if ( xyType=="quantBox" )
    {
      delta <- 0.25
      # Draw the outer whiskers.
      segments( yrs,obj[1,], yrs,obj[5,], lty=1,col=1 )
      # Overlay the box.
      for ( i in 1:length(yrs) )
        rect( yrs[i]-delta,obj[2,i], yrs[i]+delta, obj[4,i],... )
      # Add the median.
      segments( yrs-delta,obj[3,],yrs+delta,obj[3,],lty=1,col=1 )
    }
  }

  # Plot quantiles of biomass using the posterior densities.

  yrs1 <- NULL
  yrs2 <- NULL
  result1 <- NULL
  result2 <- NULL

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
  plt.qB( result1,xLim=xLim,yLim=yLim, xyType=xyType )
  axis(1, at=xLim[1]:xLim[2], tcl=tcl.val, labels=FALSE)
  axis(2, at = seq(0, yLim[2], by=yaxis.by), tcl=tcl.val, labels=FALSE)
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^plotRmcmcPOP


#plotIndex2-----------------------------2011-08-31
# Revised version of Arni's function to confine plotting to data region.
#----------------------------------------------AME
plotIndex2 <- function (model, what = "c", series = NULL, axes = TRUE, same.limits = FALSE,
    between = list(x = axes, y = axes), ylim = NULL, q = 1, bar = 1,
    log = FALSE, base = 10, main = "", xlab = "", ylab = "",
    cex.main = 1.2, cex.lab = 1, cex.strip = 0.8, cex.axis = 0.8,
    las = 1, tck = c(1, 0)/2, tick.number = 5, lty.grid = 3,
    col.grid = "white", pch = 16, cex.points = 1.2, col.points = "black",
    lty.lines = 1, lwd.lines = 4, col.lines = "dimgrey", lty.bar = 1,
    plot = TRUE, ...)
{
    panel.index <- function(x, y, subscripts, yobs, yfit, col.points,
        col.lines, ...) {
        y.range <- range(c(rep(0, !log), attr(yobs, "other"),
            y), na.rm = TRUE)
        panel.abline(v = pretty(x, tick.number), h = pretty(y.range,
            tick.number), lty = lty.grid, col = col.grid)
        panel.xyplot(x, y, type = "n", ...)
        panel.xyplot(x, yfit[subscripts], type = "l", lty = lty.lines,
            lwd = lwd.lines, col = col.lines[subscripts], ...)
        ok.Y <- !is.na(yobs[subscripts])
        if (lty.bar == 0)
            panel.xyplot(x, yobs[subscripts], col = col.points,
                ...)
        else panel.xYplot(x[ok.Y], yobs[subscripts][ok.Y], subscripts = subscripts[ok.Y],
            col = col.points[subscripts], lty.bar = lty.bar,
            ...)
    }
    if (class(model) != "scape")
        stop("The 'model' argument should be a scape object, not ",
            chartr(".", " ", class(model)), ".")
    what <- match.arg(what, c("c", "s"))
    relation <- if (same.limits)
        "same"
    else "free"
    if (what == "c") {
        if (any(names(model) == "CPUE"))
        {
        # Truncate fit to show only data years.
          minYr <- min( model$CPUE$Year[ !is.na(model$CPUE$Obs) ] )
          x <- model$CPUE[ model$CPUE$Year >= minYr, ]
        }
        else {
            what <- "s"
            cat("Commercial CPUE data (", substitute(model),
                "$CPUE) not found. Assuming user intended what=\"s\".\n",
                sep = "")
        }
    }
    if (what == "s") {
        if (any(names(model) == "Survey"))
        {
        # Truncate fit to show only data years.
            x <- model$Survey
            minYr <- min( model$Survey$Year[ !is.na(model$Survey$Obs) ] )
            x <- model$Survey[ model$Survey$Year >= minYr, ]
        }
        else stop("Found neither commercial CPUE data (", substitute(model),
            "$CPUE) nor survey abundance data (", substitute(model),
            "$Survey).\nPlease verify that ", substitute(model),
            " is a 'scape' model that contains abundance index data.")
    }
    if (is.null(series))
        series <- unique(x$Series)
    ok.series <- x$Series %in% series
    if (!any(ok.series))
        stop("Please check if the 'series' argument is correct.")
    x <- x[ok.series, ]
    if (is.numeric(x$Series))
        x$Series <- factor(paste("Series", x$Series))
    x$Obs <- x$Obs/q
    x$Fit <- x$Fit/q
    x$Hi <- x$Obs * exp(bar * x$CV)
    x$Lo <- x$Obs/exp(bar * x$CV)
    if (log) {
        x$Obs <- log(x$Obs, base = base)
        x$Fit <- log(x$Fit, base = base)
        x$Hi <- log(x$Hi, base = base)
        x$Lo <- log(x$Lo, base = base)
    }
    require(grid, quietly = TRUE, warn.conflicts = FALSE)
    require(Hmisc, quietly = TRUE, warn.conflicts = FALSE)
    require(lattice, quietly = TRUE, warn.conflicts = FALSE)
    if (trellis.par.get()$background$col == "#909090") {
        for (d in dev.list()) dev.off()
        trellis.device(color = FALSE)
    }
    col.points <- rep(col.points, length.out = nlevels(x$Series))
    col.lines <- rep(col.lines, length.out = nlevels(x$Series))
    mymain <- list(label = main, cex = cex.main)
    myxlab <- list(label = xlab, cex = cex.lab)
    myylab <- list(label = ylab, cex = cex.lab)
    myrot <- switch(as.character(las), "0" = list(x = list(rot = 0),
        y = list(rot = 90)), "1" = list(x = list(rot = 0), y = list(rot = 0)),
        "2" = list(x = list(rot = 90), y = list(rot = 0)), "3" = list(x = list(rot = 90),
            y = list(rot = 90)))
    myscales <- c(list(draw = axes, relation = relation, cex = cex.axis,
        tck = tck, tick.number = tick.number), myrot)
    mystrip <- list(cex = cex.strip)
    graph <- xyplot(Fit ~ Year+0.5 | Series, data = x, panel = panel.index,
        yobs = Cbind(x$Obs, x$Hi, x$Lo), yfit = x$Fit, as.table = TRUE,
        between = between, main = mymain, xlab = myxlab, ylab = myylab,
        par.strip.text = mystrip, scales = myscales, pch = pch,
        cex = cex.points, col.points = col.points[x$Series],
        col.lines = col.lines[x$Series], ...)
    if (is.list(ylim))
        graph$y.limits <- rep(ylim, length.out = length(series))
    else if (is.numeric(ylim))
        graph$y.limits <- rep(list(ylim), length(series))
    else {
        extremes <- as.data.frame(lapply(split(x[, c("Fit", "Hi",
            "Lo")], x$Series), range, na.rm = TRUE))
        if (!log)
            graph$y.limits <- lapply(extremes, function(x) c(0,
                1.04 * x[2]))
        else graph$y.limits <- lapply(extremes, function(x) range(x) +
            c(-0.04, 0.04) * diff(range(x)))
    }
    if (same.limits) 
        graph$y.limits <- range(unlist(graph$y.limits))
    if (plot) {
        print(graph)
        invisible(x)
    }
    else {
        invisible(graph)
    }
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^plotIndex2


#--------------------------------------------------------------------#
#                        Utility Functions                           #
#--------------------------------------------------------------------#
allEqual <- function(x)
{
  result <- all( x==x[1] )
  result
}

closeAllWin <- function()
{
  winList <- dev.list()
  if ( !is.null(winList) )
    for ( i in 1:length(winList) )
      dev.off()
}

graphics <- function( view="portrait" )
{
  if ( view=="portrait" )
    windows( width=8.5, height=11 )
  if ( view=="landscape" )
    windows( width=11, height=8.5 )
}

panLab <- function( x, y, txt, ... )
{
#  orig.par <- par()
  usr <- par( "usr" )
  par( usr=c(0,1,0,1) )
  text( x, y, txt, ... )
  par( usr=usr )
#  par( orig.par )
  return( NULL )
}

panLegend <- function( x, y, legTxt, ... )
{
#  orig.par <- par()
  usr <- par( "usr" )
  par( usr=c(0,1,0,1) )
  legend( x, y, legend=legTxt, ... )
  par( usr=usr )
#  par( orig.par )
  return( NULL )
}

#--------------------------------------------------------------------#
#                        Calculation Functions                       #
#--------------------------------------------------------------------#

#calc.projExpect------------------------2011-08-31
# Calculate the expectation of projection to reference.
# Compare refYears to projection years.
#----------------------------------------------AME
calc.projExpect <- function( obj, projObj, refYrs )
{
  policyList <- names(projObj)
  projYrs <- dimnames( projObj[[1]] )[[2]]
  refYears <- as.character(refYrs)

  nPolicies <- length(policyList)
  nProjYrs <- length(projYrs)
  nRefYrs <- length(refYears)

  # Final results.
  result <- as.list( c(1:nRefYrs) )
  names( result ) <- refYears

  # Loop over the reference years.
  for ( i in 1:nRefYrs )
  {
    # Create a results matrix for nPolicies (j), nProjYr (k).
    val <- matrix( NA, nrow=nPolicies,ncol=nProjYrs )
    dimnames( val ) <- list( policyList,projYrs )

    # Loop over the policies.
    for ( j in 1:nPolicies )
    {
      # This are the projection results for policy j.
      proj <- projObj[[j]]

      # Loop over the projection years.
      for ( k in 1:ncol(proj) )
      {
        tmp <- proj[,k] / obj[,refYears[i]]
        val[j,k] <- mean( tmp )
      }
    }
    result[[i]] <- val
  }
  cat( "\n\nExpectation of projection biomass > reference year:\n\n" )
  print( result )
  result
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^calc.projExpect


#calc.projExpect2-----------------------2011-08-31
# Calculate expectation (projection biomass / reference biomass).
#----------------------------------------------AME
calc.projExpect2 <- function( obj, projObj, refList )
{
  policyList <- names(projObj)
  nPolicies <- length(policyList)

  projYrs <- dimnames( projObj[[1]] )[[2]]
  nProjYrs <- length(projYrs)

  refYears <- refList$refYrs
  nRefYrs <- nrow( refYears )
  funs <- refList$funVec

  # Final results, each reference period stored as list element.
  result <- as.list( c(1:nRefYrs) )
  names( result ) <- dimnames(refYears)[[1]]

  # Loop over the reference years.
  for ( i in 1:nRefYrs )
  {
    # Create a results matrix for nPolicies (j), nProjYr (k).
    val <- matrix( NA, nrow=nPolicies,ncol=nProjYrs )
    dimnames( val ) <- list( policyList,projYrs )

    # Build reference years and coerce to character for indexing.
    period <- as.character( c( refYears[i,1]:refYears[i,2] ) )

    # Calculate the reference value for the performance measure.
    refVal <- calc.refVal( obj,period,funs[[i]] )

    # Loop over the catch level policies.
    for ( j in 1:nPolicies )
    {
      # These are the projection results for policy j.
      proj <- projObj[[j]]

      # Loop over the projection years.
      for ( k in 1:ncol(proj) )
      {
        tmp <- proj[,k] / refVal
        val[j,k] <- mean( tmp )
      }
    }
    result[[i]] <- val
  }
  cat( "\n\nExpectation of (projection biomass / reference period):\n\n" )
  print( data.frame( refYears,funs ) )
  cat( "\n" )
  print( result )
  result$refs <- refList
  result
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^calc.projExpect2


#calc.projProbs-------------------------2011-08-31
# Calculate the probability of being greater than refYears.
# Compare refYears to projection years.
#----------------------------------------------AME
calc.projProbs <- function( obj, projObj, refYrs )
{
  policyList <- names(projObj)
  projYrs <- dimnames( projObj[[1]] )[[2]]
  refYears <- as.character(refYrs)

  nPolicies <- length(policyList)
  nProjYrs <- length(projYrs)
  nRefYrs <- length(refYears)

  # Final results.
  result <- as.list( c(1:nRefYrs) )
  names( result ) <- refYears

  # Loop over the reference years.
  for ( i in 1:nRefYrs )
  {
    # Create a results matrix for nPolicies (j), nProjYr (k).
    val <- matrix( NA, nrow=nPolicies,ncol=nProjYrs )
    dimnames( val ) <- list( policyList,projYrs )

    # Loop over the policies.
    for ( j in 1:nPolicies )
    {
      # This are the projection results for policy j.
      proj <- projObj[[j]]

      # Loop over the projection years.
      for ( k in 1:ncol(proj) )
      {
        tmp <- proj[,k] > obj[,refYears[i]]
        val[j,k] <- sum( tmp ) / length(tmp)
      }
    }
    result[[i]] <- val
  }
  cat( "\n\nProbability of projection biomass > reference year:\n\n" )
  print( result )
  result
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^calc.projProbs


#calc.projProbs2------------------------2011-08-31
# Calculate the probability of being greater than refYears.
# Compare refYears to projection years.
#----------------------------------------------AME
calc.projProbs2 <- function( obj, projObj, refList )
{
  policyList <- names(projObj)
  nPolicies <- length(policyList)

  projYrs <- dimnames( projObj[[1]] )[[2]]
  nProjYrs <- length(projYrs)

  refYears <- refList$refYrs
  nRefYrs <- nrow( refYears )
  funs <- refList$funVec

  # Final results, each reference period stored as list element.
  result <- as.list( c(1:nRefYrs) )
  names( result ) <- dimnames(refYears)[[1]]

  # Loop over the reference years.
  for ( i in 1:nRefYrs )
  {
    # Create a results matrix for nPolicies (j), nProjYr (k).
    val <- matrix( NA, nrow=nPolicies,ncol=nProjYrs )
    dimnames( val ) <- list( policyList,projYrs )

    # Build reference years and coerce to character for indexing.
    period <- as.character( c( refYears[i,1]:refYears[i,2] ) )

    # Calculate the reference value for the performance measure.
    refVal <- calc.refVal( obj,period,funs[[i]] )

    # Loop over the catch level policies.
    for ( j in 1:nPolicies )
    {
      # These are the projection results for policy j.
      proj <- projObj[[j]]

      # Loop over the projection years.
      for ( k in 1:ncol(proj) )
      {
         tmp <- proj[,k] > refVal
         val[j,k] <- sum( tmp ) / length(tmp)
      }
    }
    result[[i]] <- val
  }
  cat( "\n\nProbability of projection biomass > reference period:\n\n" )
  print( data.frame( refYears,funs ) )
  cat( "\n" )
  print( result )
  result$refs <- refList
  result
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^calc.projProbs2


#calc.refProbs--------------------------2011-08-31
# Calculate the reference probabilities (basing on calc.projProbs2)
#----------------------------------------------AME
calc.refProbs <- function( projObj=currentProj$B, refPlist=refPointsList )  
{
  # refPlist will be a list, LRP, URP and Bmsy are defaults with values for each draw
  policyList <- names(projObj)
  nPolicies <- length(policyList)

  projYrs <- dimnames( projObj[[1]] )[[2]]
  nProjYrs <- length(projYrs)
  
  nRefPoints = length(names(refPlist))  # refYears <- refList$refYrs
    # nRefYrs <- nrow( refYears )
    # funs <- refList$funVec
  # browser(); return()
  # Final results, each reference point stored as list element.
  result <- as.list(1:nRefPoints)   # to get right number of list elements (1,2,3, but
                                    #  they'll get overwritten)
  names( result ) <- names(refPlist)

  # Loop over the reference points.
  for ( i in 1:nRefPoints )
  {
    # Create a results matrix for nPolicies (j), nProjYr (k).
    val <- matrix( NA, nrow=nPolicies,ncol=nProjYrs )
    dimnames( val ) <- list( policyList,projYrs )

    # Build reference years and coerce to character for indexing.
    # period <- as.character( c( refYears[i,1]:refYears[i,2] ) )

    # Calculate the reference value for the performance measure.
    # refVal <- calc.refVal( obj,period,funs[[i]] )
    refVal = refPlist[[i]]
    
    # Loop over the catch level policies.
    for ( j in 1:nPolicies )
    {
      # These are the projection results for policy j.
      proj <- projObj[[j]]

      # Loop over the projection years.
      for ( k in 1:ncol(proj) )
      {
         tmp <- proj[,k] > refVal
         val[j,k] <- sum( tmp ) / length(tmp)
      }
    }
    result[[i]] <- val
  }
  cat( "\n\nProbability of projection biomass > reference point:\n\n" )
  # print( data.frame( refPointsrefYears,funs ) )
  # cat( "\n" )
  print( result )
  # result$refs <- refPlist
  result
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^calc.refProbs


#calc.refVal----------------------------2011-08-31
# Calculates the reference value for performance measures.
# Returns a vector of length nrow(obj) reference values.
#----------------------------------------------AME
calc.refVal <- function( obj, refYrs, fun=mean )
{
  # Input:
  # obj: scape Biomass matrix with n rows and nYr columns.
  # refYrs: numeric years in reference period.
  # fun: The function to apply to reference period i.

  # Coerce obj to a matrix.
  obj <- as.matrix( obj )

  # Change reference years to character to allow named indexing.
  refYrs <- as.character( refYrs )

  # Extract relevant columns and apply function to rows.
  # Coerce to a matrix to accommodate single year reference.

  tmp <- matrix( obj[,refYrs],ncol=length(refYrs) )
  result <- apply( tmp,1,fun )
  result
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^calc.refVal


#getYrIdx-------------------------------2011-08-31
# Purpose is to return selected years for plotting.
# Default is to select 5 year increments.
#----------------------------------------------AME
getYrIdx <- function( yrNames,mod=5 )
{
  # Coerce to numerice and select the years modulo "mod".
  yrVals <- as.numeric( yrNames )
  idx <- yrVals %% mod==0

  # Select years from character vector yrNames.
  result <- yrNames[ idx ]
  result
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^getYrIdx


#out.pmTables---------------------------2011-08-31
# Write the decision tables to a comma delimited file.
#----------------------------------------------AME
out.pmTables <- function( obj, fileName="pm", dec=3 )
{
  nPM <- length( obj )
  for ( i in 1:(nPM-1) )
  {
    tmp <- obj[[i]]
    tmp <- format( round(tmp,digits=dec) )

    rowNames <- dimnames(tmp)[[1]]
    colNames <- dimnames(tmp)[[2]]
    catch <- as.numeric( rowNames )
    dimnames( tmp ) <- list( NULL,NULL )
    tmp <- cbind( catch, tmp )
    dimnames( tmp ) <- list( NULL, c("Catch",colNames) )
    write.table( tmp, quote=FALSE, row.names=FALSE,
      file=paste( fileName,i,".csv",sep=""),sep="," )
  }
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^out.pmTables


#stdRes.CA------------------------------2011-08-31
# This implements standardised residuals for the Awatea
# implementation of the Fournier robustified normal likelihood
# for proportions at length. Based on PJS summary of CASAL doc and ACH change to length.
#----------------------------------------------AME
stdRes.CA <- function( obj, trunc=3, myLab="Age Residuals", prt=TRUE )
{
  # Make a column for the standardised residuals.
  result <- cbind( obj,stdRes=rep(NA,nrow(obj)) )

  # Raw residuals.
  res <- obj$Obs - obj$Fit

  # Number of age bins.
  # QUESTION: Should this be from FIRST age to plus group?
  n <- length( unique(obj$Age) )

  # Kludgy, but loop through years.
  # Could reformat into matrix and vectorize.

  yrList <- sort( unique( obj$Year ) )
  for ( i in 1:length( yrList ) )
  {
    idx <- yrList[i]==obj$Year
    Z <- obj$Obs[idx]*(1.0-obj$Obs[idx]) + 0.1/n
    N <- min( obj$SS[idx],1000)
    SD <- sqrt( Z/N )
    result$stdRes[idx] <- res[idx]/SD
  }

  if ( prt )
  {
    sdRes <- sqrt( var( result$stdRes,na.rm=TRUE ) )
    sdTrunc <- ifelse( result$stdRes > trunc, trunc, result$stdRes )
    sdTrunc <- ifelse( result$stdRes < -trunc, -trunc, result$stdRes )
    sdResTrunc <- sqrt( var( sdTrunc,na.rm=TRUE ) )
    cat( "\n",myLab,"\n" )
    cat( "\n     Std. Dev. of standardised Residuals = ",sdRes,"\n" )
    cat( "     Std. Dev. of Truncated standardised Residuals = ",sdResTrunc,"\n" )
  }
  result
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^stdRes.CA


#stdRes.index---------------------------2011-08-31
# Compute the standardised residuals.
#----------------------------------------------AME
stdRes.index <- function( obj, label=NULL, prt=TRUE )
{
  res <- log(obj$Obs) - log(obj$Fit)
  stdRes <- res / obj$CV
  result <- cbind( obj,stdRes )

  if ( prt )
  {
    sdRes <- sqrt( var( stdRes,na.rm=TRUE ) )
    cat( "\n",label,"\n" )
    cat( "\n     Std. Dev. of standardised Residuals = ",sdRes,"\n" )
  }
  result
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^stdRes.index


#MAfun2---------------------------------2011-05-04
# to work out mean ages, from Rowan's MAfun from runADMB.r.
#  4th May 2011. Don't think much changed from MAfun.
# Call will be MAfun2(obj$CAc); f is regime, don't worry for here.
#-------------------------------------------RH/AME
MAfun2 = function(padata,brks=NULL)  {
       # Mean age function (Chris Francis, 2011, submitted to CJFAS))
       # S = series, y = year, a = age bin, O = observed proportions,
       # P = Predicted (fitted) proportions, N=sample size
       S=padata$Series; y=padata$Year; a=padata$Age;
       O=padata$Obs; E=padata$Fit; N=padata$SS
       if (is.null(brks)) {
	 f = paste(S,y,sep="-"); J = unique(S) } else
       {
  	 B = cut(y, breaks=brks, include.lowest=TRUE, labels=FALSE)
	 f = paste(S,B,y,sep="-"); J = unique(paste(S,B,sep="-"))
       }
       Oay  = a * O; Eay = a * E; Eay2 = a^2 * E
       mOy  = sapply(split(Oay,f),sum,na.rm=TRUE)
       mEy  = sapply(split(Eay,f),sum,na.rm=TRUE)
       mEy2 = sapply(split(Eay2,f),sum,na.rm=TRUE)
       N    = sapply(split(N,f),mean,na.rm=TRUE)
       return(list(MAobs=mOy, MAexp=mEy, Vexp=mEy2-mEy^2, N=N, J=J))
        # observed and expected mean ages, variance of expected ages
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^MAfun2


#--------------------------------------------------------------------#
#                         Plotting Functions                         #
#--------------------------------------------------------------------#

# AME introducing for YMR. Bubble plots of data:
# plt.bubbles = function( obj,
# For now just do in Sweave as quicker.

#plt.ageResidsPOP-----------------------2011-08-31
# AME changing for POP, just plotting age class resids
#  here, not qq-plot. Moving that to new function (so can do 4 on
#  page). See popScape.r for original.
#----------------------------------------------AME
plt.ageResidsPOP <- function( obj, ages=c(2,60), pct=c(5,25,50,75,95)
                             ,  main=NULL )
{
  # Input is the output from stdRes.CA
  # par( oma=c(2,1,1,1), mar=c(2,2,2,1), mfrow=c(2,1) )
  # Subset to required ages.
  obj <- obj[ (obj$Age >= ages[1]) & (obj$Age <=ages[2]), ]
  if( max(diff(sort(unique(obj$Age)))) > 1)
    {
      allAges = min(obj$Age):max(obj$Age)
      nodataAges = allAges[ !(allAges %in% obj$Age)]
      xx = split(c(obj$stdRes, rep(NA, length(nodataAges))), c(obj$Age, nodataAges))
      xpos <- boxplot( xx, whisklty=1, xlab="", ylab="", outline=FALSE )     #AME outline = FALSE removes outliers
    } else
    {            
    xpos <- boxplot( split( obj$stdRes, obj$Age ), whisklty=1, xlab="", ylab="", outline=FALSE )     #AME outline = FALSE removes outliers
    }
  abline( h=0, lty=2, col="red" )
  mtext( side=1, line=2, cex=0.8, "Age class" )
  # These wouldn't show up:
  # legend(x="topleft", "(a)", bty="n")
  # text( par("usr")[1] + 0.95*diff(par("usr")[1:2]),
  #     par("usr")[3] + 0.05*diff(par("usr")[3:4]), "(a)") 
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^plt.ageResidsPOP


#plt.ageResidsqqPOP---------------------2011-08-31
# Plotting qq plot for age class resids.
#----------------------------------------------AME
plt.ageResidsqqPOP <- function( obj, ages=c(2,60),
                pct=c(5,25,50,75,95),  main=NULL )
{
  # Subset to required ages.
  obj <- obj[ (obj$Age >= ages[1]) & (obj$Age <=ages[2]), ]

  # Plot the q-normal plot of the standardised residuals 
  qqnorm( obj$stdRes,xlab="",ylab="",main="" )
  abline( a=0,b=1 )
  abline( h=quantile(obj$stdRes,p=pct/100,na.rm=TRUE),lty=2 )
  mtext( side=1, line=2, cex=0.8, "Theoretical quantiles" )

  mtext( side=2, line=-1, cex=0.8, outer=TRUE, "Standardised Residuals" )
  if ( !is.null(main) )
    mtext( side=3, line=-0.5, cex=1.0, outer=TRUE, main )
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^plt.ageResidsqqPOP


#plt.yearResidsPOP----------------------2011-08-31
# AME adding to plot age residuals by year. Is called for comm and survs.
#  fill.in = TRUE is to add the missing years for boxplot
#  ..POP does not do qq plot. See popScape.r for previous.
#----------------------------------------------AME
plt.yearResidsPOP <- function( obj, ages=c(2,60), pct=c(5,25,50,75,95),
                           main=NULL, fill.in = TRUE, ... )
{
  # Subset to required ages - still do as don't want age 1.
  obj <- obj[ (obj$Age >= ages[1]) & (obj$Age <=ages[2]), ]
  if(fill.in)
    {
      allYears = min(obj$Year):max(obj$Year)
      nodataYears = allYears[ !(allYears %in% obj$Year)]
      xx = split(c(obj$stdRes, rep(NA, length(nodataYears))), c(obj$Year, nodataYears))
      xpos <- boxplot( xx, whisklty=1, xlab="", ylab="",
          outline=FALSE, ... )     #AME outline = FALSE removes outliers
      # browser()
    } else
    {  
      xpos <- boxplot( split( obj$stdRes, obj$Year ), whisklty=1, xlab="", ylab="", outline=FALSE, ... )     #AME outline = FALSE removes outliers
    }
  abline( h=0, lty=2, col="red" )
  mtext( side=1, line=2, cex=0.8, "Year" )
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^plt.yearResidsPOP


#plt.cohortResids-----------------------2011-08-31
# Plot age residuals by cohort.
#----------------------------------------------AME
plt.cohortResids <- function( obj, ages=c(2,59), pct=c(5,25,50,75,95),                     main=NULL )
{
  # Input is the CAc object from a Awatea res file. Ages to 59 as
  #  plus-age class will mess up year-of-birth calculation. Not automated.

  # par( oma=c(2,1,1,1), mar=c(2,2,2,1), mfrow=c(2,1) )

  # Subset to required ages - still do as don't want age 1 or 60
  #  for cohorts
  obj <- obj[ (obj$Age >= ages[1]) & (obj$Age <=ages[2]), ]
  # obj$stdRes has residuals for each age, year and both sexes. Need
  #  to assign a year of birth for each age as an extra column, then
  #  presumably just do the boxplot split using that.
  obj$birthyr = obj$Year - obj$Age
  if( max(diff(sort(unique(obj$birthyr)))) > 1)
    {
      allYears = min(obj$birthyr):max(obj$birthyr)
      nodataYears = allYears[ !(allYears %in% obj$birthyr)]
      xx = split(c(obj$stdRes, rep(NA, length(nodataYears))), c(obj$birthyr, nodataYears))
      xpos <- boxplot( xx, whisklty=1, xlab="", ylab="", outline=FALSE )     #AME outline = FALSE removes outliers
    } else
    {            
    xpos = boxplot( split( obj$stdRes, obj$birthyr ), whisklty=1, xlab="", ylab="", outline=FALSE )     #AME outline = FALSE removes outliers
    }
  abline( h=0, lty=2, col="red" )
  mtext( side=1, line=2, cex=0.8, "Year of birth" )

  # Plot the q-normal plot of the standardised residuals.
  # qqnorm( obj$stdRes,xlab="",ylab="",main="" )
  # abline( a=0,b=1 )
  # abline( h=quantile(obj$stdRes,p=pct/100,na.rm=TRUE),lty=2 )
  # mtext( side=1, line=2, cex=0.8, "Theoretical quantiles" )

  # mtext( side=2, line=0, cex=1.0, outer=TRUE, "Standardised Residuals" )
  # if ( !is.null(main) )
    mtext( side=3, line=-0.5, cex=1.0, outer=TRUE, main )

  # par( mfrow=c(1,1) )
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^plt.cohortResids


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
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^plt.allTraces


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
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^plt.expRate


#plt.idx--------------------------------2011-08-31
# AME doing postscript for POP. Adapting to five surveys for YMR
#----------------------------------------------AME
plt.idx <- function( obj,main="Residuals",save=NULL,ssnames=paste("Ser",1:9,sep=""),...)
{
  seriesList <- sort( unique( obj$Series ) )
  nseries = length(seriesList)
  #surveyFigName =c("survGIG.eps", "survQCSsyn.eps", "survQCSshr.eps", "survWCHG.eps", "survWCVI.eps", "survNFMS.eps")
  #surveyHeadName = c("GIG historical", "QCS synoptic", "QCS shrimp", "WCHG synoptic", "WCVI synoptic", "NMFS Triennial")
  surveyFigName =paste("survRes",ssnames,".eps",sep="")
  surveyHeadName = ssnames
  for ( i in 1:nseries )
  {
    idx <- seriesList[i]==obj$Series
    result <- stdRes.index( obj[idx,],
                label=paste(main,"Series",i) )
    # windows()
     postscript(surveyFigName[i],  height = 6, width = 5,
         horizontal=FALSE,  paper="special")  # was 7.5, 6.0
    # png(surveyFigName[i],  height = 4, width = 3.2,
    #     horizontal=FALSE,  paper="special")
    plt.stdResids( result,
               xLim = range(result[!is.na(result$Obs), ]$Year))
               # restrict years for plot
      mtext( side=3, line=0, cex=1.0, outer=TRUE,
           surveyHeadName[i])
    dev.off()
    # if ( !is.null(save) )
    #   savePlot( paste(save,i,sep=""), type="png" )
  }
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^plt.idx


#plotIndexNotLattice--------------------2012-00-06
# Taking some of plt.idx, but doing plot.Index NOT as lattice
# YMR - have five survey series, so adapting as required.
#----------------------------------------------AME
plotIndexNotLattice <- function( obj,objCPUE,main="",save=NULL,bar=1.96,
    ssnames=paste("Ser",1:9,sep=""), ...)
{
  # obj=obj$Survey, objCPUE=obj$CPUE
  seriesList <- sort( unique( obj$Series ) )   # sort is risky if not always in same order
  nseries = length(seriesList)
  # surveyFigName =c("survIndGIG.eps", "survIndQCSsyn.eps",
  #   "survIndQCSshr.eps")
  #surveyHeadName = c("GIG historical", "QCS synoptic", "QCS shrimp", "WCHG synoptic", "WCVI synoptic", "NMFS Triennial") # Appears again below
  surveyHeadName = ssnames
  postscript("survIndSer.eps",
    height = 8.0, width = 6.0,
    horizontal=FALSE,  paper="special")   # height was 6 for POP
  par(mfrow=c(nseries,1),mgp=c(2,0.75,0)) # Do c(3,2) for 6 surveys
  par(mai = c(0.3,0.25, 0.4,0.1))         # RH changed space & below; JAE changed for each figure # was for POP 0.45, 0.5, 0.3, 0.2
  par(omi = c(0.25,0.25,0,0.1))           # Outer margins of whole thing, inch
  yrTicks = as.numeric( obj$Year)
   
  for ( i in 1:nseries )
  {
    idx <- seriesList[i]==obj$Series
    seriesVals = obj[idx,]
    # seriesvals$Obs = seriesvals$Obs   # /q[i] - set to 1 anyway
    seriesVals$Hi <- seriesVals$Obs * exp(bar * seriesVals$CV)
    seriesVals$Lo <- seriesVals$Obs/exp(bar * seriesVals$CV)

    yearsnotNA = seriesVals[ !is.na(seriesVals$Obs), ]$Year
    # yearsPlot = min(yearsnotNA):max(yearsnotNA)
    xLim = range(yearsnotNA)
    if(i==1)
      xLimAll = xLim    # range to use in next plot
    xLimAll = range(xLim, xLimAll)    # range to use in next plot
    yLim = c(0, max(seriesVals$Hi, na.rm=TRUE))
    # postscript(surveyFigName[i],  height = 4, width = 3.2,
    #     horizontal=FALSE,  paper="special")
    #gplots::plotCI(seriesVals$Year, seriesVals$Obs, ui=seriesVals$Hi,
    plotCI(seriesVals$Year, seriesVals$Obs, ui=seriesVals$Hi,
           li=seriesVals$Lo, xlim = xLim, ylim=yLim, xlab="",
           ylab="", gap=0, pch=19)
         # restrict years for plot, does error bars
    lines(seriesVals$Year, seriesVals$Fit, lwd=2)
    axis( side=1, at=yrTicks, tcl=-0.2, labels=FALSE )
    mtext( side=3, line=0.25, cex=0.8, outer=FALSE,
          surveyHeadName[i]) #  outer=TRUE
    if(i==3)
       mtext( side=2, line=0.5, cex=1, outer=TRUE,"Relative biomass")
    if(i==5)
       mtext(side=1, line=0.5, cex=1, outer=TRUE, "Year")      
  }                         # cex was 0.8 for POP
  dev.off()

  # And again, but with the same year axis for each. Think will
  #  be instructive to see
  postscript("survIndSer2.eps",
    height = 8.0, width = 6.0,
    horizontal=FALSE,  paper="special")   # height was 6 for POP
  par(mfrow=c(nseries,1),mgp=c(2,0.75,0)) #  do c(3,2) for 6 series, and change axes labelling
  par(mai = c(0.3,0.25, 0.4,0.1))         # RH changed space & below; JAE changed for each figure # was for POP 0.45, 0.5, 0.3, 0.2
  par(omi = c(0.25,0.25,0,0.1))           # Outer margins of whole thing, inch
  yrTicks = as.numeric( obj$Year)
  ymaxsurvIndSer3 = 0           # For ymaxsurvIndSer3.eps
  xLimmaxsurvIndSer3 = NA
  for ( i in 1:length(seriesList) )
  {
    idx <- seriesList[i]==obj$Series
    seriesVals = obj[idx,]
    # seriesvals$Obs = seriesvals$Obs   # /q[i] - set to 1 anyway
    seriesVals$Hi <- seriesVals$Obs * exp(bar * seriesVals$CV)
    seriesVals$Lo <- seriesVals$Obs/exp(bar * seriesVals$CV)

    yearsnotNA = seriesVals[ !is.na(seriesVals$Obs), ]$Year
    # yearsPlot = min(yearsnotNA):max(yearsnotNA)
    # xLim = range(yearsnotNA)
    yLim = c(0, max(seriesVals$Hi, na.rm=TRUE))
    # For axis for survIndSer3.eps:
    ymaxsurvIndSer3 = max(ymaxsurvIndSer3,
         max(seriesVals$Obs, na.rm=TRUE)/
         mean(seriesVals$Obs, na.rm=TRUE) )
    xLimmaxsurvIndSer3 = range(xLimmaxsurvIndSer3, yearsnotNA,
      na.rm=TRUE)     # setting xLimmaxsurvIndSer3 = NA above

    #gplots::plotCI(seriesVals$Year, seriesVals$Obs, ui=seriesVals$Hi,
    plotCI(seriesVals$Year, seriesVals$Obs, ui=seriesVals$Hi,
           li=seriesVals$Lo, xlim = xLimAll, ylim=yLim, xlab="",
           ylab="", gap=0, pch=19)
         # restrict years for plot, does error bars
    lines(seriesVals$Year, seriesVals$Fit, lwd=2)
    axis( side=1, at=yrTicks, tcl=-0.2, labels=FALSE )
    mtext( side=3, line=0.25, cex=0.8, outer=FALSE,
          surveyHeadName[i]) #  outer=TRUE
    if(i==3)
       mtext( side=2, line=0, cex=1, outer=TRUE,"Relative biomass")
    if(i==5)
       mtext(side=1, line=0, cex=1, outer=TRUE, "Year")      
  }                         # cex was 0.8 for POP
  dev.off()
  
  # And again, but all series on same plot, normalised to their means
  #  to see trends. Maybe add CPUE also.
  postscript("survIndSer3.eps",
    height = 6.0, width = 6.0,
    horizontal=FALSE,  paper="special")   # height was 6 for POP
    yrTicks = as.numeric( obj$Year)
    # Set up plot
  plot(NA, xlim = xLimmaxsurvIndSer3, ylim = c(0, ymaxsurvIndSer3),
       xlab="Years", ylab="Survey indices normalised by means")
  for ( i in 1:length(seriesList) )
    {
      idx <- seriesList[i]==obj$Series
      seriesVals = obj[idx,]
      # seriesVals$Hi <- seriesVals$Obs * exp(bar * seriesVals$CV)
      # seriesVals$Lo <- seriesVals$Obs/exp(bar * seriesVals$CV)

      yearsnotNA = seriesVals[ !is.na(seriesVals$Obs), ]$Year
      points(yearsnotNA, seriesVals$Obs[ !is.na(seriesVals$Obs)] /
          mean(seriesVals$Obs, na.rm=TRUE), pch=i, col=i, type="o")
    }
  # Now draw on CPUE series also:
  nseries = length(seriesList)  #  Need nseries from surveys for col
  # obj=objCPUE
  seriesListCPUE <- sort( unique( objCPUE$Series ) )
    # sort risky if not always in same order.
  for ( i in 1:length(seriesListCPUE) )
    {
      idx <- seriesListCPUE[i]==objCPUE$Series
      seriesVals = objCPUE[idx,]
      yearsnotNA = seriesVals[ !is.na(seriesVals$Obs), ]$Year
      points(yearsnotNA, seriesVals$Obs[ !is.na(seriesVals$Obs)] /
          mean(seriesVals$Obs, na.rm=TRUE), pch=i+nseries,
          col=i+nseries, type="o", lty=2)
    }
  dev.off()

  # And again, but just first series and CPUE to see conflicting
  #  signal. Axes not automated.
  postscript("survIndSer4.eps",
    height = 5.0, width = 5.0,
    horizontal=FALSE,  paper="special")   # height was 6 for POP
    yrTicks = as.numeric( obj$Year)
    # Set up plot
  plot(NA, xlim = c(1965, 1995), ylim = c(0, 2.5),
       xlab="Years", ylab="Normalised GIG survey index and CPUE")
  abline(h=1, col="grey")
  for ( i in 1)   # Just doing first series
    {
      idx <- seriesList[i]==obj$Series
      seriesVals = obj[idx,]
      # seriesVals$Hi <- seriesVals$Obs * exp(bar * seriesVals$CV)
      # seriesVals$Lo <- seriesVals$Obs/exp(bar * seriesVals$CV)

      yearsnotNA = seriesVals[ !is.na(seriesVals$Obs), ]$Year
      points(yearsnotNA, seriesVals$Obs[ !is.na(seriesVals$Obs)] /
          mean(seriesVals$Obs, na.rm=TRUE), pch=19, type="o")
    }
  # Now draw on CPUE series also:
  nseries = length(seriesList)  #  Need nseries from surveys for col
  # obj=objCPUE
  seriesListCPUE <- sort( unique( objCPUE$Series ) )
    # sort risky if not always in same order.
  for ( i in 1:length(seriesListCPUE) )
    {
      idx <- seriesListCPUE[i]==objCPUE$Series
      seriesVals = objCPUE[idx,]
      yearsnotNA = seriesVals[ !is.na(seriesVals$Obs), ]$Year
      points(yearsnotNA, seriesVals$Obs[ !is.na(seriesVals$Obs)] /
          mean(seriesVals$Obs, na.rm=TRUE), pch=5,
          col="blue", type="o", lty=2)
    }
  dev.off()

  
  # And again, but big figures for ease of viewing. 
  #  postscript("survIndSer4-%d.eps",
  #  height = 8.0, width = 6.0,
  #  horizontal=FALSE,  paper="special")   # height was 6 for POP
  #par(mfrow=c(2,1),mgp=c(2,0.75,0))
  #par(mai = c(0.25,0.25, 0.25, 0.1))      # JAE changed  for each figure # was for POP 0.45, 0.5, 0.3, 0.2
  #par(omi = c(0.25,0.25,0,0.1))           # Outer margins of whole thing, inch
  yrTicks = as.numeric( obj$Year)
   
  for ( i in 1:length(seriesList) )
  {
    idx <- seriesList[i]==obj$Series
    seriesVals = obj[idx,]
    # seriesvals$Obs = seriesvals$Obs   # /q[i] - set to 1 anyway
    seriesVals$Hi <- seriesVals$Obs * exp(bar * seriesVals$CV)
    seriesVals$Lo <- seriesVals$Obs/exp(bar * seriesVals$CV)

    yearsnotNA = seriesVals[ !is.na(seriesVals$Obs), ]$Year
    # yearsPlot = min(yearsnotNA):max(yearsnotNA)
    xLim = range(yearsnotNA)
    yLim = c(0, max(seriesVals$Hi, na.rm=TRUE))
    # postscript(surveyFigName[i],  height = 4, width = 3.2,
    #     horizontal=FALSE,  paper="special")
    postscript(paste("survIndSer4-", i, ".eps", sep=""),
    height = 6.0, width = 6.0,
    horizontal=FALSE,  paper="special")   # height was 6 for POP
    #gplots::plotCI(seriesVals$Year, seriesVals$Obs, ui=seriesVals$Hi,
    plotCI(seriesVals$Year, seriesVals$Obs, ui=seriesVals$Hi,
           li=seriesVals$Lo, xlim = xLim, ylim=yLim, xlab="Year",
           ylab="Relative biomass", gap=0, pch=19)
         # restrict years for plot, does error bars
    lines(seriesVals$Year, seriesVals$Fit, lwd=2)
    axis( side=1, at=yrTicks, tcl=-0.2, labels=FALSE )
    mtext( side=3, line=0.25, cex=0.8, outer=FALSE,
          surveyHeadName[i]) #  outer=TRUE
    #if(i==3)
    #   mtext( side=2, line=0.5, cex=1, outer=TRUE,"Relative biomass")
    #if(i==5)
    #   mtext(side=1, line=0.5, cex=1, outer=TRUE, "Year")
    dev.off()
  }                         # cex was 0.8 for POP
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^plotIndexNotLattice


#plotChains-----------------------------2011-05-11
# Plots cumulative fequency of 'nchains' by partitioning one trace.
# Revised from 'plotTracePOP'
# mcmc = data.frame e.g, 'currentMCMC$P' from object created by 'importMCMC'.
#-------------------------------------------------
plotChains = function (mcmc, nchains=3, pdisc=0.1, 
     axes=FALSE, same.limits=FALSE, between=list(x=axes,y=axes), div=1, span=1/4,
     log=FALSE, base=10, main=NULL, xlab=NULL, ylab=NULL, 
     cex.main=1.2, cex.lab=1, cex.strip=0.8, cex.axis=0.8, las=0, 
     tck=0.5, tick.number=5, lty.trace=1, lwd.trace=1, col.trace="grey", 
     lty.median=1, lwd.median=1, col.median="black", lty.quant=2, lwd.quant=1, 
     col.quant="black", plot=TRUE, probs=c(0.025, 0.5, 0.975), ...)  # AME probs
{
	panel.trace <- function(x, y, ...) {

		panel.xyplot(x, y, type="n")
		chainlink = rep(1:nchains,f)
		for (i in 1:nchains) {
			z = is.element(chainlink,i)
			panel.xyplot(x[z], y[z], type="l", lty=lty.trace, lwd=2, col=rep(col.trace,nchains)[i])
		}
		#panel.xyplot(x, y, type = "l", lty = lty.trace, lwd = lwd.trace, col = col.trace)
	}

	mcmc=mcmc[(round(pdisc*nrow(mcmc))+1):nrow(mcmc),]  # get rid of the first 'pdisc' (e.g., 10%)
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
	f = rep(round(n/nchains),nchains-1); f = c(f,n-sum(f))
	p <- ncol(mcmc)
	dat <- data.frame(Factor = ordered(rep(names(mcmc), each = n), 
		names(mcmc)), Draw = rep(1:n, p), Chain = rep(rep(1:nchains,f),p), Value = as.vector(as.matrix(mcmc)))
	require(grid, quietly = TRUE, warn.conflicts = FALSE)
	require(lattice, quietly = TRUE, warn.conflicts = FALSE)
	if (trellis.par.get()$background$col == "#909090") {
		for (d in dev.list()) dev.off()
		trellis.device(color = FALSE)
	}
	mymain <- list(label = main, cex = cex.main)
	myxlab <- list(label = xlab, cex = cex.lab)
	myylab <- list(label = ylab, cex = cex.lab)
	myrot <- switch(as.character(las), `0` = 0, `1` = 0, `2` = 0, `3` = 90)     # AME changed '0' = 90 to 0
	myscales <- list(x = list(draw=axes, relation=relation, cex=cex.axis, tck=tck, tick.number=tick.number, rot=myrot), 
		y= list(draw=axes, relation=relation, cex=cex.axis, tck=tck, tick.number=tick.number, rot=myrot))
	mystrip <- list(cex = cex.strip)

	dat$Index = paste(dat$Factor,dat$Chain,sep="-")
	vList = split(dat$Value,dat$Index)
	qList = sapply(vList,function(x){
		xsort = sort(x)
		xscal = xsort - min(xsort)
		ycumu = cumsum(xscal)/sum(xscal)
		out = cbind(x=xsort,y=ycumu)
		return(out) }, simplify=FALSE )
	dat$CumFreq = dat$ValueSort = NA
	for (i in names(qList)) {
		z = is.element(dat$Index,i)
		dat$ValueSort[z] = qList[[i]][,"x"]
		dat$CumFreq[z]   = qList[[i]][,"y"]
	}
	graph <- xyplot(CumFreq ~ ValueSort  | Factor, panel = panel.trace, 
		data = dat, as.table = TRUE, between = between, main = mymain, 
		xlab = myxlab, ylab = myylab, par.strip.text = mystrip, 
		scales = myscales, ylim=c(0,1), ...)
#browser();return()
	if (plot) {
		print(graph)
		invisible(dat)
	}
	else {
		invisible(graph)
	}
}
#P=currentMCMC$P
#plotChains(mcmc=P,axes=TRUE,between=list(x=0.15,y=0.2),col.trace=c("green","red","blue"),xlab="Sample",ylab="Cumulative Frequency")


# Plotting CPUE and fit with error bars, copying plotIndexNotLattice
# from above.
plotCPUE <- function( obj,main="",save=NULL,bar=1.96, yLim=NULL, ...)
{                            # obj=currentRes$CPUE
  seriesList <- sort( unique( obj$Series ) )   # sort is risky if not always in same order
  nseries = length(seriesList)
  # surveyFigName =c("survIndGIG.eps", "survIndQCSsyn.eps",
  #   "survIndQCSshr.eps")
  surveyHeadName = c("CPUE")
  postscript("CPUEser.eps",
        height = switch(nseries,5,8,9), width = 6.0,
        horizontal=FALSE,  paper="special")   # height was 6 for POP
  par(mfrow=c(nseries,1),mai=c(0.75,0.75,0,0.1),omi=c(0,0,0.25,0),mgp=c(2,.75,0))
  # par(mai = c(0.25, 0.5, 0.3, 0.2)) # JAE changed  for each figure was for POP 0.45, 0.5, 0.3, 0.2
  # par(omi = c(0.45,0.1,0,0))  # Outer margins of whole thing, inch
  yrTicks = as.numeric( obj$Year)
   
  for ( i in 1:nseries )
  {
    idx <- seriesList[i]==obj$Series
    seriesVals = obj[idx,]
    # seriesvals$Obs = seriesvals$Obs   # /q[i] - set to 1 anyway
    seriesVals$Hi <- seriesVals$Obs * exp(bar * seriesVals$CV)
    seriesVals$Lo <- seriesVals$Obs/exp(bar * seriesVals$CV)
    # browser(); return()
    yearsnotNA = seriesVals[ !is.na(seriesVals$Obs), ]$Year
    # yearsPlot = min(yearsnotNA):max(yearsnotNA)
    xLim = range(yearsnotNA)
    if(i==1)
      xLimAll = xLim    # range to use in next plot
    xLimAll = range(xLim, xLimAll)    # range to use in next plot
    #if(is.null(yLim))     
    yLim = c(0, max(seriesVals$Hi, na.rm=TRUE))
    # postscript(surveyFigName[i],  height = 4, width = 3.2,
    #     horizontal=FALSE,  paper="special")
    #gplots::plotCI(seriesVals$Year, seriesVals$Obs, ui=seriesVals$Hi,
    plotCI(seriesVals$Year, seriesVals$Obs, ui=seriesVals$Hi,
           li=seriesVals$Lo, xlim = xLim, ylim=yLim, xlab="Year",
           ylab=paste("CPUE index:",seriesList[i]), gap=0, pch=19)
         # restrict years for plot, does error bars
    lines(seriesVals$Year, seriesVals$Fit, lwd=2)
    axis( side=1, at=yrTicks, tcl=-0.2, labels=FALSE )
    # mtext( side=3, line=0.25, cex=0.8, outer=FALSE, surveyHeadName[i]) #  outer=TRUE
    # if(i==3)  mtext( side=2, line=-0.5, cex=1, outer=TRUE,"Relative biomass")
    # if(i==5)  mtext(side=1, line=0, cex=1, outer=TRUE, "Year")
  }                         # cex was 0.8 for POP
  dev.off()
}


#plt.mcmcGraphs-------------------------2012-10-16
#  AME editing (with *AME*) to give a policy option 
#  to be specified in run-master.Snw. 16th August 2012.
#-------------------------------------------AME/AM
plt.mcmcGraphs <-
function (mcmcObj, projObj, save = FALSE, 
          ylim.recruitsMCMC = NULL, ylim.exploitMCMC = NULL,
          ylim.VBcatch = NULL, ylim.BVBnorm = NULL,
          xlim.snail = NULL, ylim.snail = NULL,
          plotPolicies = names(projObj$Y[1:6]), # *AME*
          onePolicy=names(projObj$Y[2]),
          mpd=list() )        # *AME*

          # plotPolicies is 6 policies projections to plot *AME*
          # onePolicy is one to use for some figures *AME*
          #*AME*xlim.pdfrec was =c(0, 200000). Put options for others
          #  that will be useful if want to scale two model runs
          #  to the same ylim. If NULL then fits ylim automatically.
  {     
	panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
	{
		usr <- par("usr"); on.exit(par(usr))
		par(usr = c(0, 1, 0, 1))
		r <- abs(cor(x, y))
		txt <- format(c(r, 0.123456789), digits=digits)[1]
		txt <- paste(prefix, txt, sep="")
		text(0.5, 0.5, txt, cex = 1.4)
	}
    postscript("recruitsMCMC.eps", height = 5, width = 6.2, horizontal = FALSE, 
        paper = "special")
    plotRmcmcPOP(currentMCMC$R, yLim = ylim.recruitsMCMC) # *AME*
    dev.off()
    postscript("exploitMCMC.eps", height = 5, width = 6.2, horizontal = FALSE, 
        paper = "special")
    plotRmcmcPOP(currentMCMC$U, yLab = "Exploitation rate",
                 yLim = ylim.exploitMCMC, yaxis.by = 0.01)
    dev.off()
    postscript("pdfBiomass%d.eps", height = 7, width = 6.2, horizontal = FALSE, 
        paper = "special", onefile = FALSE)
    plotDensPOP(currentMCMC$B/1000, xlab = "Female spawning biomass, Bt (1000 t)", 
        between = list(x = 0.2, y = 0.2), ylab = "Density", lwd.density = 2, 
        same.limits = TRUE, layout = c(4, 6), lty.outer = 2, mpd=mpd[["mpd.B"]]/1000)
    dev.off()
    postscript("pdfRecruitment%d.eps", height = 7, width = 6.2, 
        horizontal = FALSE, paper = "special", onefile = FALSE)
    plotDensPOP(currentMCMC$R/1000, xlab = "Recruitment, Rt (1000s)", 
        between = list(x = 0.2, y = 0.2), ylab = "Density", lwd.density = 2, 
        same.limits = TRUE, layout = c(4, 6), lty.median = 2, 
        lty.outer = 2, mpd=mpd[["mpd.R"]]/1000)
    dev.off()
    postscript("pdfParameters.eps", horizontal = FALSE, paper = "special", 
        height = 7, width = 6.2)
    plotDensPOPparsPrior(mcmcObj$P, lty.outer = 2, between = list(x = 0.3, 
        y = 0.2),mpd=mpd[["mpd.P"]])
    dev.off()
    postscript("selectivityMCMC.eps", height = 5, width = 6.2, 
        horizontal = FALSE, paper = "special")
    plot(1:10, main = "This will be selectivities from MCMC")
    dev.off()
    postscript("traceRecruits.eps", horizontal = FALSE, paper = "special", 
        height = 7, width = 6.2)
    plotTracePOP(mcmcObj$R[, getYrIdx(names(mcmcObj$R))]/1000, 
        axes = TRUE, between = list(x = 0.2, y = 0.2), xlab = "Sample", 
        ylab = "Recruitment, Rt (1000s)", mpd=mpd[["mpd.R"]][getYrIdx(names(mcmcObj$R))]/1000)
    dev.off()
    postscript("traceBiomass.eps", horizontal = FALSE, paper = "special", 
        height = 7, width = 6.2)
    plotTracePOP(mcmcObj$B[, getYrIdx(names(mcmcObj$B))]/1000, 
        axes = TRUE, between = list(x = 0.2, y = 0.2), xlab = "Sample", 
        ylab = "Female spawning biomass, Bt (1000 t)", mpd=mpd[["mpd.B"]][getYrIdx(names(mcmcObj$B))]/1000)
    dev.off()
    postscript("traceParams.eps", horizontal = FALSE, paper = "special", 
        height = 7, width = 6.2)
    idx <- apply(mcmcObj$P, 2, allEqual)
    plotTracePOP(mcmcObj$P[, !idx], axes = TRUE, between = list(x = 0.2, 
        y = 0.2), xlab = "Sample", ylab = "Parameter estimate", mpd=mpd[["mpd.P"]][!idx])
    dev.off()
    postscript("splitChain.eps", horizontal = FALSE, paper = "special", 
        height = 7, width = 6.2)
    plotChains(mcmc = currentMCMC$P, axes = TRUE, between = list(x = 0.15, 
        y = 0.2), col.trace = c("green", "red", "blue"), xlab = "Sample", 
        ylab = "Cumulative Frequency", pdisc = 0.001)
    dev.off()
    postscript("VBcatch.eps", horizontal = FALSE, paper = "special", 
        height = 5, width = 6.2)
    plotVBcatch(currentMCMC$VB, currentRes, yLim = ylim.VBcatch)
    dev.off()
    postscript("BVBnorm.eps", horizontal = FALSE, paper = "special", 
        height = 5, width = 6.2)
    plotBVBnorm(currentMCMC, xLeg = 0.02, yLeg = 0.2,
                yLim = ylim.BVBnorm)
    dev.off()
    options(scipen = 10)
    postscript("Bproj.eps", horizontal = FALSE, paper = "special", 
        height = 7, width = 6.2)
    plt.quantBio(currentMCMC$B, currentProj$B, xyType = "quantBox", 
        policy = plotPolicies,    # *AME*
        save = FALSE)
    dev.off()
    postscript("Rproj.eps", horizontal = FALSE, paper = "special", 
        height = 7, width = 6.2)
    plt.quantBio(currentMCMC$R, currentProj$R, xyType = "quantBox", 
        policy = plotPolicies,    # *AME*
        save = FALSE, yaxis.lab = "Recruitment (1000s)")
    dev.off()
    postscript("RprojOnePolicy.eps", horizontal = FALSE, paper = "special",                    # *AME* (filename)
        height = 5, width = 6.2)
    plt.quantBioBB0(currentMCMC$R, currentProj$R, xyType = "quantBox", 
        policy = onePolicy, save = FALSE, xaxis.by = 10, yaxis.lab = "Recruitment (1000s)")    # *AME* (onePolicy)
    dev.off()
    postscript("snail.eps", height = 5, width = 6.2, horizontal = FALSE, 
        paper = "special")
    plotSnail(currentMCMC$BoverBmsy, currentMCMC$UoverUmsy, p = c(0.1, 
        0.9), xLim = xlim.snail, yLim = ylim.snail)
    dev.off()
    # Doing 6 on a page:
    npr = 6
    nuP = length(use.Pnames)
    npp = ceiling(nuP/npr) # number of pairs plots
    for (i in 1:npp) {
      if (i<npp) ii = (1:npr)+(i-1)*npr
      else ii = (nuP-npr+1):nuP 
      postscript(paste("pairs",i,".eps",sep=""),height=7,width=7,horizontal=FALSE,paper="special")
      #pairs(currentMCMC$P[, ii], pch=20, cex=0.2, gap=0)
      pairs(currentMCMC$P[, ii], pch = 20, cex = 0.2, gap = 0, lower.panel=panel.cor)
      dev.off()
    }
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^plt.mcmcGraphs


#plt.mpdGraphs--------------------------2012-10-16
# Plot the MPD graphs to encapsulated postscript files.
#----------------------------------------------AME
plt.mpdGraphs <- function( obj, save=FALSE, ssnames=paste("Ser",1:9,sep=""))
{
	#AME some actually MCMC. # Doing as postscript now. # Taking some out for ymr.
	# Does all MPD graphs below.  If save=TRUE then PNG file saved.
	closeAllWin()

  # Plot the biomass and catch.
  # plotB2( obj,main=mainTitle )              
  # plotBmcmcPOP(currentMCMC$B, currentRes)    # AME new one
  # if ( save )                           # though not using now as
  #  savePlot( "biomass",type="png" )    # spawning and catch on same
                                        # axis isn't really sensible
                                        # Doing plotVBcatch in
                                        #  plt.mcmcGraphs above

  # AME adding, plot exploitation rate, not writing new function:
  postscript("exploit.eps", height = 5, width = 6.2, horizontal=FALSE,  paper="special")
  B = obj$B
  xlim = range(B$Year,na.rm=TRUE)
  uflds = findPat("U",names(B))
  B$U = apply(B[,uflds,drop=FALSE],1,sum)
  ylim  = range(B$U,na.rm=TRUE)
  par(mfrow=c(1,1),mgp=c(2,.75,0))
  plot(0,0,xlim=xlim,ylim=ylim,type="n",xlab="Year",ylab="Exploitation rate")
  lines(B$Year, B$U, col="black")
  points(B$Year, B$U, cex=0.9, pch=21, col="black", bg="white") 
  dev.off()
  # if ( save )
  #  savePlot( "exploit",type="png" )
  
  # AME had added recruits.eps for POP, but that's MCMC, so moving
  #  to plt.mcmcGraphs, changing that filename to recruitsMCMC.eps
  #  and just adding here to do MPD for recruits.
  postscript("recruits.eps", height = 5, width = 6.2,
              horizontal=FALSE,  paper="special")
  plot(obj$B$Year, obj$B$R, type="o", xlab="Year",
       ylab="Recruitment, Rt (1000s)", ylim=c(0, max(obj$B$R, na.rm=TRUE)))
  dev.off()

  
  # Plot the numbers at age and recruits. Taking out for ymr
  # windows()
	  # plotN( obj, main=mainTitle, ages=c(2:60) )
  # #plotN( obj, main=mainTitle )
  # if ( save )
  #   savePlot( "numAtAge", type="png" )

  # Plot the selectivity.
  # windows()
  # objRed = obj      # Reduced object, just plotting Selectivity to 20
  # objRed$Sel = objRed$Sel[objRed$Sel$Age < 21,]
  # plotSel( objRed, main=paste(mainTitle,"Selectivity"), xlim=c(0,20))
  # if ( save )
  #   savePlot( "selectivity", type="png" )

  postscript("selectivity.eps", height = 5, width = 6.2,
              horizontal=FALSE,  paper="special")
  objRed = obj      # Reduced object, just plotting Selectivity to 20
  objRed$Sel = objRed$Sel[objRed$Sel$Age < 21,]
  plotSel( objRed, main=paste(mainTitle,"Selectivity"), xlim=c(0,20))
  dev.off()
  
  # Plot the catch at age
  # NOTE: There is a bug in plotCA that prevents plotting multiple
  #       series given a list of character vectors in series.
  #         ACH: I'm not sure if this applies to CL

  seriesList <- sort( unique( obj$CAc$Series) )
  maxcol = 4
  CAc.yrs = sapply(split(obj$CAc$Year,obj$CAc$Series),unique,simplify=FALSE)
  CAc.nyrs = sapply(CAc.yrs,length)

  for ( i in 1:length(seriesList) )
  {
            # AME dividing years into 4 groups - note no 1985, 86,
            #  or 88 data - see below, doing more automatically saving as .eps
    # windows()
    # plotCA( obj, series=i, main=paste("Comm",mainTitle,"Series",i), what="c" )
    # if ( save )
    #   savePlot( paste("catchAgeComm",i,sep=""), type="png" )
    # windows()
    # AME trying to do two on one, from Paul Murrel's chap4:
    # HERE
    # plot1 <-  plotCA( obj, series=i, main=paste("Comm",mainTitle,"Series",i), what="c", years=1978:1984 )
    # plot2 <-  plotCA( obj, series=i, main=paste("Comm",mainTitle,"Series",i), what="c", years=1978:1984 )
    # print(plot1, position=c(0,0.2,1,1), more=TRUE)  #didn't work.
    # ----
#    par(mfrow=c(2,1))
#    plotCA( obj, series=i, main=paste("Comm",mainTitle,"Series",i), what="c", years=1978:1984 )
   # if(0==1) {   # these were four files total, now doing less per page below:
   # windows()
   # plotCA( obj, series=i, main=paste("Comm",mainTitle,"Series",i), what="c", years=1978:1984 )
   # if ( save )
      # savePlot( paste("catchAgeComm",i,sep=""), type="png" )
   #   savePlot( paste("catchAgeComm1"), type="png" ) # AME took out i
   # windows()
   # plotCA( obj, series=i, main=paste("Comm",mainTitle,"Series",i, "(no 1985, 86 or 88)"), what="c", years=1985:1994 )  # no 85, 86 or 88
   # if ( save )
   #   savePlot( paste("catchAgeComm2"), type="png" ) # AME took out i
   # windows()
   # plotCA( obj, series=i, main=paste("Comm",mainTitle,"Series",i), what="c", years=1995:2001 )
   # if ( save )
   #   savePlot( paste("catchAgeComm3"), type="png" ) # AME took out i
   # 
   # windows()
   # plotCA( obj, series=i, main=paste("Comm",mainTitle,"Series",i), what="c", years=2002:2009 )
   # if ( save )
   #   savePlot( paste("catchAgeComm4"), type="png" ) # AME took out i
   #} # end if(0==1)
   # unique(currentRes$CAc$Year) is, do four a page:
   # 1978 1979 1980 1981 1982 1983 1984 1987 1989 1990 1991 1992 1993 1994 1995
   #  1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009
  # uniqueCAcyears = unique(currentRes$CAc$Year)  # length=29, so do groups of 4
  # # numpanels = 4
  # for(jj in 0:7) { 
  #    windows()
  #    startyr = uniqueCAcyears[jj*4 + 1]
  #    # print(startyr)
  #    endyr = uniqueCAcyears[(jj+1)*4]
  #      if(jj == 7) {endyr = uniqueCAcyears[length(uniqueCAcyears)]}
  #    plotCA( obj, series=i, main=paste("Comm",mainTitle,"Series",i), what="c", years=startyr:endyr )
      # browser()
  #  if ( save )
      # savePlot( paste("catchAgeComm",i,sep=""), type="png" )
  #    savePlot( paste("catchAgeComm", jj, sep=""), type="png" ) # AME took out i
  # } # end for jj in 0:7
  # END of AME manually doing multiple plots. This is much easier:
  # AME doing again, but more automatically, saving as .eps as for plotDens
  #  as multiple pages. Starting writing plotCAPOP, but could do all
  #  by changing arguments to plotCA(...). Haven't referred to i here, as for
  #  POP there's only one commercial set of data.
  CA.key = list(text = list(lab= c("Observed", "Predicted")),
    lines = list(col= c("black", "red"),  cex= c(0.3, 0.3)),
    type=c("p", "l"), x = 0.78, y = -0.04, pch=c(20,20), lwd=1, between=0.3)
                    # pch[2] doesn't get used, as type[2]="l". Have to match
                    #  up if change options in plotCA(currentRes, ....)
  CAc.sex = unique(obj$CAc$Sex)
  for(plot.sex in CAc.sex)
    {                       # For YMR, changing height from 4.5, and
                            #  layout from c(4,4), to put all on one
      age.layout = rev(c(ceiling(CAc.nyrs[i]/maxcol),maxcol)) # backwards in stupid lattice
      postscript(paste("ageComm", plot.sex,i,".eps", sep=""),
        height = switch(age.layout[2],4,6,9,9,9,9,9,9,9), width = 6.5,
        horizontal=FALSE,  paper="special", onefile=FALSE)
      plotCA( obj, what="c", ylab="Proportion", xlab="Age class",
        sex=plot.sex, layout= age.layout, key=CA.key, main=plot.sex,
        pch=20, cex.points=0.5, col.lines=c("red", "red"), lwd.lines=2 ,series=i)
                 # col.lines otherwise boys are blue
                 # Tried using rbind to add dummy data for 1985, 1986
                 #  and 1988, but didn't work:
                 # add = c(1, 1985, 0, "Female", 1, 60, 60, 0, 0)
                 #   xxx$CAc = rbind(xxx$CAc, add)
                 #  
    dev.off()
    }    # end of plot.sex loop
  }

  # AME adding - plotting the CA survey data for all three series.
  #  doing .eps below
    # {
      if ( exists( "currentRes" ) )
      {
        seriesList <- sort( unique( obj$CAs$Series) )
        # for ( i in 1:length(seriesList) )
        # {
        #   windows()
        #   plotCA( currentRes, series=i, main=paste("Survey",mainTitle,"Series",i), what="s" )
        #   if(save)
        #     savePlot( paste("catchAgeSurvey", i, sep=""), type="png" )
        #  }
      }
      # }

  # AME adding - CA survey data for all three (no fitting for 3)
  #   postscript files.
  #      seriesList <- sort( unique( obj$CAs$Series) )
  #      for ( i in 1:length(seriesList) )
  #      {
  #        windows()
  #        plotCA( currentRes, series=i, main=paste("Survey",mainTitle,"Series",i), what="s" )
  #        if(save)
  #          savePlot( paste("catchAgeSurvey", i, sep=""), type="png" )
  #        }

  # Survey Age Fits (modified by RH 2012-08-01)
  # Shifting the key here down slightly
  CAs.key = list(text = list(lab= c("Obs", "Pred")),
    lines = list(col= c("black", "red"),  cex= c(0.3, 0.3)),
    type=c("p", "l"), x = 0.78, y = -0.16, pch=c(20,20), lwd=1, between=0.3)
        # AME: pch[2] doesn't get used, as type[2]="l". Have to match up if change options in plotCA(currentRes, ....)

  CAs.sex = unique(obj$CAs$Sex)
  # GIG, survey 1      # should do a loop,         only first two
  #ageSurveyFigName =c("ageSurvGIG", "ageSurvQCSsyn", "ageSurvQCSshr")
  age.height = c(3.0, 2.8, 3)    # For .eps figs:
  age.width = c(3.5, 7, 3)    # For .eps figs:
  seriesList <- sort( unique( obj$CAs$Series) )
  ageSurveyFigName = paste("ageSurv",ssnames[seriesList],sep="") # friggin nightmare
  age.layout = list()
  age.layout[[1]] = c(2,1)
  age.layout[[2]] = c(4,1)
  #age.layout[[3]] = c(1,1)
    for ( i in 1:length(seriesList) ) {
      ii = seriesList[i]; zi=is.element(obj$CAs$Series,ii)
      if (!any(zi)) next
      iyr = unique(obj$CAs$Year[zi]); nyr = length(iyr)
      ncol = min(nyr,4); nrow=ceiling(nyr/ncol)
      for(plot.sex in CAs.sex) {
        postscript(paste(ageSurveyFigName[i], plot.sex,ii,".eps", sep=""),
          #sep=""),  height = age.height[i], width = age.width[i], # RH disabled
          height=3.5, width=2*ncol, horizontal=FALSE, paper="special", onefile=FALSE)
        plotCA( obj, what="s", series = ii, ylab="Proportion",
          #xlab="Age class", sex=plot.sex, layout=age.layout[[i]], key=CAs.key, main=plot.sex, # RH disabled
          xlab="Age class", sex=plot.sex, layout=c(ncol,nrow), key=CAs.key, main=plot.sex, # RH: Stupid trellis appears to take (columns,rows) for the layout.
          pch=20, cex.points=0.5, col.lines=c("red", "red"), lwd.lines=2 )
              # AME: col.lines otherwise boys are blue
              # AME: Tried using rbind to add dummy data for 1985, 1986, and 1988, but didn't work:
              #      add = c(1, 1985, 0, "Female", 1, 60, 60, 0, 0)
              #      xxx$CAc = rbind(xxx$CAc, add)
    dev.off()
    }    # end of plot.sex loop
  }      # end of seriesList loop
  
  
  # Plot the fishery index (CPUE data I think)
  # windows()
  # plotIndex(obj, main=paste(mainTitle,"CPUE"), what="c", bar=1.96 )
  # #plotIndex2( obj, main=mainTitle, what="c", bar=1.96 )
  # if ( save )
  #   savePlot( "fishIndex", type="png" )
  
  # Plot the survey indices.
  # windows()
  # plotIndex( obj, main=paste(mainTitle,"Survey Indices"), what="s", bar=1.96 )
  # #plotIndex2( obj, main=mainTitle, what="s", bar=1.96 )
  # if ( save )
  #   savePlot( "surveyIndex", type="png" )
  
  # Now do on one plot as postscript:
  plotIndexNotLattice(obj$Survey, obj$CPUE, ssnames=ssnames)

  # Single plot of CPUE:
  plotCPUE(obj$CPUE, yLim=c(0, 200))

  # Plot standardised residuals. Now doing four plots on one page.
  #  Commercial.
  postscript("commAgeResids.eps",
        height = 8.5, width = 6.8,
        horizontal=FALSE,  paper="special")
  par(mai = c(0.45, 0.5, 0.1, 0.2)) # JAE changed  for each figure
  par(omi = c(0.45,0.1,0,0))      # Outer margins of whole thing, inch

  stdRes.CA.CAc = stdRes.CA( obj$CAc )
  #  Outliers don't get plotted, except for qq plot
  par(mfrow=c(4,1))
  plt.ageResidsPOP( stdRes.CA.CAc, main="" ) 
  plt.yearResidsPOP(stdRes.CA.CAc)
  plt.cohortResids(stdRes.CA.CAc)   # cohort resid, by year of birth
  plt.ageResidsqqPOP(stdRes.CA.CAc)
  dev.off()

  # And now for surveys.
  # AME adding - plotting the CA residuals for the two surveys:
  seriesList <- sort( unique( obj$CAs$Series) )  
  nseries = length(seriesList)
  for ( i in 1:nseries ) {   # POP no fits for survey 3
		ii = seriesList[i]
    postscript(paste("survAgeResSer", ii, ".eps", sep=""),
        height = 8.5, width = 6.8,
        horizontal=FALSE,  paper="special")
    par(mai = c(0.45, 0.5, 0.1, 0.2)) # JAE changed  for each figure
    par(omi = c(0.45,0.1,0,0))      # Outer margins of whole thing, inch
    stdRes.CA.CAs = stdRes.CA( obj$CAs[
          obj$CAs$Series == ii,] )
    #  Outliers don't get plotted, except for qq plot
    par(mfrow=c(4,1))
    plt.ageResidsPOP(stdRes.CA.CAs, main="" ) 
    plt.yearResidsPOP(stdRes.CA.CAs)
    plt.cohortResids(stdRes.CA.CAs)   # cohort resid, by year of birth
    plt.ageResidsqqPOP(stdRes.CA.CAs)
    dev.off()
    }
  # Here plot the mean age for catch and surveys
  postscript("meanAge.eps", height = 6, width = 6.2,
              horizontal=FALSE,  paper="special")
  par(mfrow=c(nseries+1,1), mai=c(0.75,0.75,0.5,0.1))
  MAc = MAfun2(obj$CAc)       # catch mean age
  MAs = MAfun2(obj$CAs)       # surveys mean age
  plot(as.numeric(substring(names(MAc$MAobs), 3)), MAc$MAobs,
       xlab="Year", ylab="Mean age",
       ylim=c(0, max(MAc$MAobs, MAc$MAexp)),mgp=c(2,0.75,0))
      # substring takes off the 1- at the start. This is commercial.
  points(as.numeric(substring(names(MAc$MAexp), 3)), MAc$MAexp,
         pch=20, type="o")
  mtext( side=3, line=0.25, cex=1.2, outer=FALSE,
          "Commercial") #  outer=TRUE

  MAsSurvNum = as.numeric(substring(names(MAs$MAobs), 1, last=1))
         # 1 or 2 for YMR
  #surveyHeadName = c("GIG historical", "QCS synoptic", "QCS shrimp", "WCHG synoptic", "WCVI synoptic")
  surveyHeadName = ssnames[MAs$J]

  for ( i in 1:length(unique(MAsSurvNum) ))
    {
      plot(as.numeric(substring(names(MAs$MAobs[MAsSurvNum ==
          unique(MAsSurvNum)[i]]), 3)), MAs$MAobs[MAsSurvNum==
          unique(MAsSurvNum)[i]], xlab="Year", ylab="Mean age",
          ylim=c(0, max(MAs$MAobs[MAsSurvNum==unique(MAsSurvNum)[i]],
          MAs$MAexp[MAsSurvNum== unique(MAsSurvNum)[i]] )),mgp=c(2,0.75,0))
      # ylim=range(MAs$MAobs[MAsSurvNum==unique(MAsSurvNum)[i]],
      #     MAs$MAexp[MAsSurvNum== unique(MAsSurvNum)[i]]  ))
      points(as.numeric(substring(names(MAs$MAobs[MAsSurvNum ==
          unique(MAsSurvNum)[i]]), 3)), MAs$MAexp[MAsSurvNum==
          unique(MAsSurvNum)[i]],
          pch=20, type="o")
      mtext( side=3, line=0.25, cex=1.2, outer=FALSE,
          surveyHeadName[i]) #  outer=TRUE
    }
  dev.off()

  # Plot stock-recruitment function (based on MPD's)
  # xLimSR and yLimSR fixed here for YMR to have Run 26 and 27 figs
  #  on same scales. Use these first two values to scale to data:
  # xLimSR =c(0, max(obj$B$SB))
  # yLimSR = c(0, max(obj$B$R, na.rm=TRUE))
  #xLimSR = c(0, max(c(max(obj$B$SB),45000)))   # so it draw bigger if necessary
  #yLimSR = c(0, max(c(max(obj$B$R, na.rm=TRUE),55000)))
  xLimSR = c(0, 1.5*max(obj$B$SB,na.rm=TRUE))   # so it draw bigger if necessary
  xxx = (seq(0, xLimSR[2], length.out=100))
  yyy = srFun(xxx)
  yLimSR = c(0, 1.1*max(c(yyy,obj$B$R),na.rm=TRUE))
  postscript("stockRecruit.eps", height = 5, width = 6.2,
              horizontal=FALSE,  paper="special")
  plot(xxx, yyy, lwd=2, xlim=xLimSR,
       ylim=yLimSR, type="l",
       xlab=expression(paste("Spawning biomass in year ",
           italic(t), "-1, ", italic(B)[t-1], " (t)", sep="")),
       ylab=expression(paste("Recruitment in year ",
           italic(t), ",", italic(R)[t], " (1000s)"), sep="") )
  #points(obj$B$SB, obj$B$R)
  #text(obj$B$SB, obj$B$R, labels=substring(as.character(years), 3), cex=0.5)
  text(obj$B[-length(years), "SB"], obj$B[-1, "R"], labels = substring(as.character(years), 3), cex = 0.5)
  dev.off()
  
  #windows()
  #plt.lengthResids( stdRes.CL( obj$CLs ),
  #  main=paste("Survey",mainTitle,"Series",i) )
  #if ( save )
  #  savePlot( "surveyLengthResids", type="png" )
  # plt.idx( obj$CPUE,  main="Commercial Fishery",save="fishResids" )
  # plt.idx( obj$Survey,main="Survey",save="surveyResids" )
  closeAllWin()
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^plt.mpdGraphs


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
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^plt.numR


#plt.quantBio---------------------------2011-08-31
# From popScapeRuns2.r  - AME now replaced yLim to force 0.
# This prints out tables (if run from command line), so be good to
#  use as template for decisions tables once we have MSY.
#----------------------------------------------AME
plt.quantBio <- function( obj, projObj=NULL, policy=NULL,
                  p=c(0.025,0.25,0.5,0.75,0.975),
                  xyType="lines",
                  lineType=c(3,2,1,2,3),
                  refLines=NULL, xLim=NULL, yLim=NULL,
                  userPrompt=FALSE, save=TRUE,
                  yaxis.lab="Spawning biomass" )
{
  plt.qB <- function( obj, xyType="lines", new=TRUE, xLim, yLim, line.col="black", ... )   #AME line.col="black", col is for filling rect
  {
    if ( new )
      plot( xLim,yLim, type="n", xlab="",ylab="" )

    yrs <- as.numeric(dimnames(obj)[[2]])

    # Connect the quantiles with lines.
    if ( xyType=="lines" )
    {
      for ( i in 1:nrow(obj) )
      {
        # Plot reconstructed biomass.
        lines( yrs,obj[i,], lty=lineType[i],... )
      }
    }

    # ARK vertical line-dot plot.
    # Assumes that five quantiles requested, with median as one.
    if ( xyType=="lineDot" )
    {
      points( yrs,obj[2,], pch=3, cex=0.5,... )
      points( yrs,obj[3,], pch=1,... )
      points( yrs,obj[4,], pch=3, cex=0.5,... )
      segments( yrs,obj[1,], yrs,obj[5,], lty=1,... )
    }

    # Quantile boxplots - assumes five quantiles.
    if ( xyType=="quantBox" )
    {
      delta <- 0.25
      # Draw the outer whiskers.
      segments( yrs,obj[1,], yrs,obj[5,], lty=1, col=line.col) #, ... )
             #AME col=1 removed, col=line.col ... so can have red for projs
      # Overlay the box.
      for ( i in 1:length(yrs) )
        rect( yrs[i]-delta,obj[2,i], yrs[i]+delta, obj[4,i],
             border=line.col, col="white")#AME border,col=NA (empty)
      # Add the median.
      segments( yrs-delta,obj[3,],yrs+delta,obj[3,],lty=1,col="black")     #line.col )   #AME black
    }
  }

  # Plot quantiles of biomass using the posterior densities.
  # If proj!=NULL then add the projections for all policies.
  # If policy!=NULL than plot only the specified policy.

  # Plotting ranges for reconstruction (1) and projection (2).
  # ***
  nCol <- min( length(policy),3 )
  if ( !is.null(policy) )
    nRow <- length( policy )/nCol

  nRow <- 3
  nCol <- 2

  par( oma=c(2,2,1,1), mar=c(2,2,1.5,1), mfrow=c(nRow,nCol), userPrompt )     # AME mfcol -> mfrow to fill left to right

  yrs1 <- NULL
  yrs2 <- NULL
  result1 <- NULL
  result2 <- NULL

  # Calculate the quantiles of the reconstructed biomass.
  result1 <- apply( obj,2,quantile,probs=p )
  yrs1 <- as.numeric(dimnames(result1)[[2]])

  if ( is.null(yLim) )
    # yLim <- range( result1 )
    yLim = c(0, max( result1 ))

  # Reconstructed biomass.
  if ( is.null(projObj) )
  {
    plt.qB( result1,xLim=range(yrs1),yLim=yLim, xyType=xyType )
    mtext( side=1, line=2, cex=1.0, "Year" )
    mtext( side=2, line=2, cex=1.0, "Biomass" )
  }

  # Reconstructed biomass and projections.
  if ( !is.null(projObj) )
  {
    # Get the policies to be plotted.
    if ( is.null(policy) )
      policyList <- names( projObj )
    else
      policyList <- policy

    # Loop over the policies.
    result2 <- as.list(policyList)
    names( result2 ) <- policyList

    iPage <- 1
    nPolicies <- length(policyList)
    for ( j in 1:nPolicies )
    {
      # Calculate quantiles of the projected biomass for policy.
      pol <- policyList[j]
      result2[[j]] <- apply( projObj[[pol]],2,quantile,probs=p )

      # cat( "\n\nQuantiles of projection for policy = ",
      #   policyList[j],"\n" )
      # print( result2[[j]] )

      yrs2 <- as.numeric( dimnames(result2[[j]])[[2]] )
      if ( is.null(xLim) )
        xLim <- range( yrs1,yrs2 )
      # yLim <- range( result1,result2[[j]] )  # AME to get same axes
      # yLim <- c(0,yLim[2])

      # Plot the quantiles of the biomass.
      if ( xyType=="quantBox" )
        plt.qB( result1, xyType=xyType, new=TRUE, xLim,yLim, col=NA, line.col="black", med.col="red" )      
      else
        plt.qB( result1, xyType=xyType, new=TRUE, xLim,yLim,col="red" )

      if ( !is.null(refLines) )
        abline( v=refLines,lty=4 )

      # Plot the quantiles of the projected biomass.
      if ( xyType=="quantBox" )
        plt.qB( result2[[j]], xyType=xyType, new=FALSE, xLim,yLim, line.col="red")  # AME: col fills in box, I want to change line cols 
      else
        plt.qB( result2[[j]], xyType=xyType, new=FALSE, xLim,yLim, col="red" )

      #for ( i in 1:nrow(result2[[j]]) )
      #  lines( yrs2,result2[[j]][i,],lty=lineType[i],lwd=2,col=2 )

      abline( v=yrs2[1]-0.5, lty=2 )

      mfg <- par( "mfg" )
      if ( mfg[1]==1 & mfg[2]==1 )
      {
        mtext( side=1, line=0.8, cex=1.0, outer=TRUE, "Year" )
                                       #AME line=0 changed
        mtext( side=2, line=0.8, cex=1.0, outer=TRUE, yaxis.lab,
           ) # "   and Spawning
      }
      mtext( side=3, line=0.25, cex=0.8,
        paste( "Catch strategy:",policyList[j]) )

      # Last panel on page or last policy.
      if ( mfg[1]==mfg[3] & mfg[2]==mfg[4] | j==nPolicies )
      {
        if ( save )
          savePlot( paste( "policyProj",iPage,sep=""),type="png" )
        iPage <- iPage + 1

        if ( j < nPolicies )
        {
          windows()
          par( oma=c(2,2,1,1), mar=c(2,2,1.5,1),
               mfcol=c(nRow,nCol), userPrompt )
        }
      }
    }
  }
  par( mfrow=c(1,1) )

  val <- list( recon=result1, proj=result2 )
  val
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^plt.quantBio


#plt.quantBioBB0------------------------2011-08-31
# This is for B/B0 for each run, with just one projection. Doing one
#  plot here instead of multiple, so taking some from plotBVBnorm
#  which did one for each run. Don't think used for biomasses.
#  Now using for single recruitment projection plot
#----------------------------------------------AME
plt.quantBioBB0 <- function( obj, projObj=NULL, policy=NULL,
                  p=c(0.025,0.25,0.5,0.75,0.975),
                  xyType="lines",
                  lineType=c(3,2,1,2,3),
                  refLines=NULL, xLim=NULL, yLim=NULL,
                  userPrompt=FALSE, save=TRUE, main="", cex.main="",
                  tcl.val=-0.2, xaxis.by = 1, yaxis.by=10000,
                  xaxis.lab = "Year", yaxis.lab= "Spawning biomass" )
{
  plt.qB <- function( obj, xyType="lines", new=TRUE, xLim, yLim, line.col="black", med.col="black", ... )   #AME line.col="black", col is for filling rect med.col is for median
  {
    if ( new )
      plot( xLim,yLim, type="n", xlab=xaxis.lab,ylab=yaxis.lab )

    yrs <- as.numeric(dimnames(obj)[[2]])

    # Connect the quantiles with lines.
    if ( xyType=="lines" )
    {
      for ( i in 1:nrow(obj) )
      {
        # Plot reconstructed biomass.
        lines( yrs,obj[i,], lty=lineType[i],... )
      }
    }

    # ARK vertical line-dot plot.
    # Assumes that five quantiles requested, with median as one.
    if ( xyType=="lineDot" )
    {
      points( yrs,obj[2,], pch=3, cex=0.5,... )
      points( yrs,obj[3,], pch=1,... )
      points( yrs,obj[4,], pch=3, cex=0.5,... )
      segments( yrs,obj[1,], yrs,obj[5,], lty=1,... )
    }

    # Quantile boxplots - assumes five quantiles.
    if ( xyType=="quantBox" )
    {
      delta <- 0.25
      # Draw the outer whiskers.
      segments( yrs,obj[1,], yrs,obj[5,], lty=1, col=line.col) #, ... )
             #AME col=1 removed, col=line.col ... so can have red for projs
      # Overlay the box.
      for ( i in 1:length(yrs) )
        rect( yrs[i]-delta,obj[2,i], yrs[i]+delta, obj[4,i],
             border=line.col, col="white")#AME border,col=NA (empty)
      # Add the median.
      # segments( yrs-delta,obj[3,],yrs+delta,obj[3,],lty=1,col=line.col )   #AME black
      segments( yrs-delta,obj[3,],yrs+delta,obj[3,],lty=1,col=med.col )   #AME black
    }
  }

  # Plot quantiles of biomass using the posterior densities.
  # If proj!=NULL then add the projections for all policies.
  # If policy!=NULL than plot only the specified policy.

  # Just one plot now. Deleted some.

  yrs1 <- NULL
  yrs2 <- NULL
  result1 <- NULL
  result2 <- NULL

  # Calculate the quantiles of the reconstructed biomass.
  result1 <- apply( obj,2,quantile,probs=p )
  yrs1 <- as.numeric(dimnames(result1)[[2]])

  if ( is.null(yLim) )
    yLim <- range( result1 )

  # Reconstructed biomass.
  if ( is.null(projObj) )
  {
    plt.qB( result1,xLim=range(yrs1),yLim=yLim, xyType=xyType )
    # mtext( side=1, line=2, cex=1.0, "Year" )
    # mtext( side=2, line=2, cex=1.0, "Biomass" )
  }

  # Reconstructed biomass and projections.
  if ( !is.null(projObj) )
  {
    # Get the policies to be plotted.
    if ( is.null(policy) )
      policyList <- names( projObj )
    else
      policyList <- policy

    # Loop over the policies.
    result2 <- as.list(policyList)
    names( result2 ) <- policyList

    iPage <- 1
    nPolicies <- length(policyList)
    for ( j in 1:nPolicies )
    {
      # Calculate quantiles of the projected biomass for policy.
      pol <- policyList[j]
      result2[[j]] <- apply( projObj[[pol]],2,quantile,probs=p )

      # cat( "\n\nQuantiles of projection for policy = ",
      #   policyList[j],"\n" )
      # print( result2[[j]] )

      yrs2 <- as.numeric( dimnames(result2[[j]])[[2]] )
      if ( is.null(xLim) )
        xLim <- range( yrs1,yrs2 )
      # yLim <- range( result1,result2[[j]] )  # AME to get same axes
      # yLim <- c(0,yLim[2])

      # Plot the quantiles of the biomass.
      if ( xyType=="quantBox" )
        plt.qB( result1, xyType=xyType, new=TRUE, xLim,yLim, col="black", line.col="black", med.col="red" )      
      else
        plt.qB( result1, xyType=xyType, new=TRUE, xLim,yLim, col="red") # line.col="red")

      if ( !is.null(refLines) )
        abline( v=refLines,lty=4 )

      # Plot the quantiles of the projected biomass.
      if ( xyType=="quantBox" )
        plt.qB( result2[[j]], xyType=xyType, new=FALSE, xLim,yLim, line.col="red", med.col="black")  # AME: col fills in box, I want to change line cols 
      else
        plt.qB( result2[[j]], xyType=xyType, new=FALSE, xLim,yLim, col="red" )

      #for ( i in 1:nrow(result2[[j]]) )
      #  lines( yrs2,result2[[j]][i,],lty=lineType[i],lwd=2,col=2 )

      abline( v=yrs2[1]-0.5, lty=2 )
      axis(1, at=seq(xLim[1], xLim[2], by=xaxis.by),
           tcl=tcl.val, labels=FALSE)
      axis(2, at = seq(0, yLim[2], by=yaxis.by),
           tcl=tcl.val, labels=FALSE)
      
      # mfg <- par( "mfg" )
      # if ( mfg[1]==1 & mfg[2]==1 )
      # {
      #  mtext( side=1, line=0.8, cex=1.0, outer=TRUE, "Year" )
      #                                 #AME line=0 changed
      #  mtext( side=2, line=0.8, cex=1.0, outer=TRUE,
      #        "Spawning biomass" ) # "   and Spawning
      #}
      #mtext( side=3, line=0.25, cex=0.8,
      #  paste( "Policy:",policyList[j]) )

      # Last panel on page or last policy.
      #if ( mfg[1]==mfg[3] & mfg[2]==mfg[4] | j==nPolicies )
      #{
      #  if ( save )
      #    savePlot( paste( "policyProj",iPage,sep=""),type="png" )
      #  iPage <- iPage + 1
      #
      #  if ( j < nPolicies )
      #  {
      #    windows()
      #    par( oma=c(2,2,1,1), mar=c(2,2,1.5,1),
      #         mfcol=c(nRow,nCol), userPrompt )
      #  }
      # }
    }
  }
  # par( mfrow=c(1,1) )

  val <- list( recon=result1, proj=result2 )
  val
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^plt.quantBioBB0


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
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^plt.ssbVbCatch


#plt.stdResids--------------------------2011-08-31
# Plot standardised residuals (AME adding xlim's for POP).
#----------------------------------------------AME
plt.stdResids <- function( obj, pct=c(5,25,50,75,95),
                           main=NULL, yLim=NULL, xLim=xLim )
{
  # Input is a data.frame with columns "Year", "stdRes", "Fit".
  # cex.val = 1.2       # size of text, if want to shrink down
  par( oma=c(1,1,2,1), mar=c(3,3,1,1), mfrow=c(3,1) )

  if ( is.null(yLim) )
    yLim <- range( obj$stdRes, na.rm=TRUE )
	ptcex=1; ptpch=19

  # Plot the standardised residuals against time.
  plot( obj$Year, obj$stdRes, type="n", xlab="", ylab="", ylim=yLim,
       xlim = xLim)
  abline( h=0, lty=2 )
  points( obj$Year, obj$stdRes, cex=ptcex, pch=ptpch) #, bg="orange" )
  mtext( side=1, line=2, cex=0.8, "Year" )

  # Plot the standardised residuals against predicted values.
  plot( log(obj$Fit), obj$stdRes, type="n", xlab="", ylab="", ylim=yLim )
  abline( h=0, lty=2 )
  points( log(obj$Fit), obj$stdRes, cex=ptcex, pch=ptpch) #, bg="orange" )
  mtext( side=1, line=2, cex=0.8, "Predicted" )

  # Plot the q-normal plot of the standardised residuals.
  wiggle=qqnorm( obj$stdRes,xlab="",ylab="",main="" , pch=ptpch, cex=ptcex)
  abline( a=0,b=1 )
  abline( h=quantile(obj$stdRes,p=pct/100,na.rm=TRUE),lty=2 )
  points(wiggle, pch=ptpch, cex=ptcex) #, bg="orange")
  mtext( side=1, line=2, cex=0.8, "Theoretical quantiles" )

  mtext( side=2, line=-0.5, cex=0.8, outer=TRUE, "Standardised residuals" )
  if ( !is.null(main) )
    mtext( side=3, line=0, cex=1.0, outer=TRUE, main )
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^plt.stdResids


#importMCMC.ddiff-----------------------2011-08-31
#  Import Functions for PJS Delay Difference Model
# Purpose is to make an scapeMCMC object identical in format to
# the result of importMCMC from PJS delay difference model output.
# The difference is that B is biomass defined by ddiff model.
#----------------------------------------------AME
importMCMC.ddiff <- function()
{
  # Get the likelihood MCMCs. PJS includes other derived parameters here.
  L <- read.table( "mcmclike.csv", header=TRUE, sep="," )

  # Get the parameter MCMCs.
  params <-read.table( "mcmcparams.csv", header=TRUE, sep="," )
  rdevID <- grep( "rdev",names(params) )
  P <- params[ ,setdiff(c(1:ncol(params)),rdevID ) ]

  # Get the biomass MCMCs and strip "biom" off to leave only "yyyy".
  B <- read.table( "mcmcbiom.csv", header=TRUE, sep="," )
  names( B ) <- substring( names( B ),5 )

  # Get the recruitments.  They are already in the mcmcparams file
  # with headers "rdev_yyyy".
  R <- params[ ,rdevID ]

  # Now rename with year only as per scapeMCMC projection object.
  names( R ) <- gsub( "rdev_","",names(R) )

  result <- list( L=L, P=P, B=B, R=R )
  result
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^importMCMC.ddiff


#importProj.ddiff-----------------------2011-08-31
# Purpose is to make an scapeMCMC object identical in format to
# the result of importProj from PJS delay difference model output.
# The difference is that B is biomass defined by ddiff model.
#----------------------------------------------AME
importProj.ddiff <- function( yrVal="2006" )
{
  # Get the biomass projection - PJS does 1 year ahead projection.
  # The column "X" appears as the last column because trailing "," exist in the mcmcprojbiom.csv file.
  # Note also that "cat=" in csv file becomes "cat." in read.table.
  projbiom <- read.table( "mcmcprojbiom.csv",header=TRUE,sep="," )
  projNames <- gsub( "cat.","",names( projbiom ) )
  names( projbiom ) <- projNames

  # Remove the last column of NA values created by trailing ",".
  projbiom <- projbiom[ ,names(projbiom)!="X" ]
  projNames <- names( projbiom )

  # Convert the columns of projbiom to matrices as list elements.
  B <- as.list( 1:length(projNames) )
  names( B ) <- projNames
  for ( j in 1:length(projNames) )
  {
    B[[j]] <- matrix( projbiom[,j],ncol=1 )
    dimnames( B[[j]] ) <- list( NULL,yrVal )
  }
  # Get the harvest rate projection - PJS does 1 year ahead projection.
  # These appear to be F's rather than annual harvest rates?
  # The column "X" appears as the last column because trailing ","
  # exist in the mcmcprojbiom.csv file.
  # Note also that "cat=" in csv file becomes "cat." in read.table.
  projhr <- read.table( "mcmcprojhr.csv",header=TRUE,sep="," )
  projNames <- gsub( "cat.","",names( projhr ) )
  names( projhr ) <- projNames

  # Remove the last column of NA values created by trailing ",".
  projhr <- projhr[ ,names(projhr)!="X" ]
  projNames <- names( projhr )

  # Convert instantaneous F's to annual rates.
  projU <- (1.0-exp(-projhr))

  # Calculate yield (catch) at annual rate.
  catch <- projbiom * projU

  # Convert the columns of projbiom to matrices as list elements.
  Y <- as.list( 1:length(projNames) )
  names( Y ) <- projNames
  for ( j in 1:length(projNames) )
  {
    Y[[j]] <- matrix( catch[,j],ncol=1 )
    dimnames( Y[[j]] ) <- list( NULL,yrVal )
  }

  result <- list( B=B, Y=Y )
  result
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^importProj.ddiff


#msyCalc--------------------------------2012-10-16
# To load in MSY.out and calculated the MSY.
# Call this function with msy = msyCalc().
#----------------------------------------------AME
msyCalc = 
function (dir = getwd(), error.rep = 1) 
{
# rewriting msyCalc here to report the convergence numbers:
	control = readLines(paste(dir, "/Yields.ctl", sep = ""))
	maxProj = as.numeric(control[3])
	tolerance = as.numeric(control[5])
	test = read.table(paste(dir, "/MSY.out", sep = ""), header = TRUE, sep = "\t")
	num.draws = dim(test)[1]
	nProjIndex = grep("nProj", names(test))
	nProjMat = test[, nProjIndex]   # matrix of number of projections
	nProjMatTF = rowSums(nProjMat > maxProj - 1)   # sum by row
	if (error.rep == 1 & (sum(nProjMatTF) > 0)) {
		stop(paste("Simulations reach maximum year for", sum(nProjMatTF), 
			"of the", num.draws * dim(nProjMat)[2], "simulations, so need to run for longer or reduce tolerance to reach equilibrium"))
	}
	yieldIndex = grep("Yield", names(test))
	yieldMat = test[, yieldIndex]
	uIndex = grep("U", names(test))
	uMat = test[, uIndex]
	VBIndex = grep("VB_", names(test))
	VBMat = test[, VBIndex]
	SBIndex = grep("SB_", names(test))
	SBMat = test[, SBIndex]
	imsy = apply(yieldMat, 1, which.max)
	if (error.rep == 1 & max(imsy) == dim(yieldMat)[2]) {
		stop("Need to use a higher max U to reach the MSY, for at least one draw")
	}
	msy = vector()
	umsy = vector()
	VBmsy = vector()
	Bmsy = vector()
	nProj = vector()
	for (i in 1:num.draws) {
		ind = imsy[i]
		msy[i] = yieldMat[i, ind]
		umsy[i] = uMat[i, ind]
		VBmsy[i] = VBMat[i, ind]
		Bmsy[i] = SBMat[i, ind]
		nProj[i] = nProjMat[i, ind]
	}
	return(list(yield = msy, u = umsy, VB = VBmsy, B = Bmsy, 
		nProj = nProj, uMin=uMat[,1], uMax=uMat[,dim(uMat)[2]], imsy = imsy, maxUind = rep(dim(yieldMat)[2], num.draws),  maxProj = rep(maxProj, num.draws), tolerance = rep(tolerance, num.draws), nProjMatTF=nProjMatTF))
	# imsy is index of tested u's for which you get umsy, so
	#  so if it's 1 or maxUind for an MCMC then that one 
	#  reached the bounds of u. Report that below in Sweave.
	# uMin, uMax are vectors of the min/max tested u's
	#  same for all MCMC samples, but need a vector to use
	#  sapply below (and same for the following variables).
	# maxProj is maximum number of projection years tried
	#  (from Yields.ctl file), so need to report if that's
	#  reached for any of the nProj. Again, need a vector.
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^msyCalc


#refPoints------------------------------2011-08-31
# Call from Sweave as  refPoints() or, in full:
# refPoints(currentMCMC, currentProj, currentMSY, refLevels=c(0.4, 0.8, 1))
#----------------------------------------------AME
refPoints <- function( mcmcObj=currentMCMC, projObj=currentProj,
                     msyObj=currentMSY, refLevels=c(0.4, 0.8, 1))
                     # refLevels are %age of msyObj
{
  refPlist = as.list(c("LRP", "URP", "Bmsy"))  # Can't have 0.4Bmsy
  names(refPlist) = c("LRP", "URP", "Bmsy")   # as numeric at start. '0.4Bmsy'
  for(i in 1:length(refLevels))
    {
    refPlist[[i]] = refLevels[i] * msyObj$B
    }
  return(refPlist)
}

refPointsB0 <- function( mcmcObj=currentMCMC, projObj=currentProj,
                     B0Obj=B0.MCMC, refLevels=B0refLevels,
                     refNames=B0refNames)
  {
  refPlist = as.list(refNames)
  names(refPlist) = c(refNames)
  for(i in 1:length(refLevels))
    {
    refPlist[[i]] = refLevels[i] * B0Obj
    }
  return(refPlist)
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^refPoints


#plotSnail------------------------------2012-08-17
# Plot snail-trail plots for MCMC analysis.
#   AME: replacing "2010" with as.character(currYear - 1)
#----------------------------------------------AME
plotSnail = function (BoverBmsy, UoverUmsy, p=c(0.1,0.9), xLim=NULL, yLim=NULL, Lwd=2)
{
    BoverBmsy.med = apply(BoverBmsy, 2, median)
    UoverUmsy.med = apply(UoverUmsy, 2, median)
    BoverBmsy.med = BoverBmsy.med[-length(BoverBmsy.med)]
    colPal = colorRampPalette(c("grey95", "grey30"))
    n = length(BoverBmsy.med)
    if (is.null(xLim)) {
        xLim = c(0, max(c(BoverBmsy.med, quantile(currentMCMC$BoverBmsy[, 
            as.character(currYear - 1)], p[2]), 1)))
    }
    if (is.null(yLim)) {
        yLim = c(0, max(c(UoverUmsy.med, quantile(currentMCMC$UoverUmsy[, 
            as.character(currYear - 1)], p[2]), 1)))
    }
    plot(BoverBmsy.med, UoverUmsy.med, type = "l", xlim = xLim, 
        ylim = yLim, xlab = expression(paste(B[t]/B[msy])), ylab = expression(paste(u[t]/u[msy])), 
        col = "grey", lwd = Lwd)
    points(BoverBmsy.med, UoverUmsy.med, type = "p", pch = 19, 
        col = colPal(n))
    points(BoverBmsy.med[1], UoverUmsy.med[1], pch = 19, col = "blue")
    points(BoverBmsy.med[as.character(currYear - 1)], UoverUmsy.med[as.character(currYear - 1)], pch = 19, 
        col = "red")
    segments(quantile(BoverBmsy[, as.character(currYear - 1)], p[1]), UoverUmsy.med[as.character(currYear - 1)], 
        quantile(BoverBmsy[, as.character(currYear - 1)], p[2]), UoverUmsy.med[as.character(currYear - 1)], 
        col = "red")
    segments(BoverBmsy.med[as.character(currYear - 1)], quantile(UoverUmsy[, as.character(currYear - 1)], 
        p[1]), BoverBmsy.med[as.character(currYear - 1)], quantile(UoverUmsy[, as.character(currYear - 1)], 
        p[2]), col = "red")
    abline(h = 1, col = "grey", lwd = Lwd)
    abline(v = 0.4, col = "grey", lwd = Lwd)
    abline(v = 0.8, col = "grey", lwd = Lwd)
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^plotSnail


#plotDensPOPparsPrior-------------------2011-08-31
# Adding the prior automatically.
#----------------------------------------------AME
plotDensPOPparsPrior <-
    function (mcmc, probs = c(0.025, 0.975), points = FALSE, axes = TRUE, 
    same.limits = FALSE, between = list(x = axes, y = axes), 
    div = 1, log = FALSE, base = 10, main = NULL, xlab = NULL, 
    ylab = NULL, cex.main = 1.2, cex.lab = 1, cex.strip = 0.8, 
    cex.axis = 0.7, las = 0, tck = 0.5, tick.number = 5, lty.density = 1, 
    lwd.density = 3, col.density = "black", lty.median = 2, lwd.median = 1, 
    col.median = "darkgrey", lty.outer = 3, lwd.outer = 1, col.outer = "darkgrey", 
    pch = "|", cex.points = 1, col.points = "black", plot = TRUE,
    MPD.height = 0.04, mpd=mcmc[1,], ...)    # MPD.height, how far up to put MPD
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
        panel.abline(v = median(x), lty = lty.median,
            lwd = lwd.median, 
            col = col.median)
        panel.xyplot(mpd[panel.number()], current.panel.limits()$ylim[2]*MPD.height,
            pch=19, col="red") # AME, MPD. 0.04 of way up
                                        #  assumes ..ylim[1]=0
        panel.xyplot(mpd[panel.number()], current.panel.limits()$ylim[2]*MPD.height,
            pch=1, col="black") #AME
        # panel.curve(priorDistList[[panel.number()]], min=-1, current.panel.limits()$xlim[1], current.panel.limits()$xlim[2], col="blue")
        # panel.curve(priorDistList[[1]], col="blue")
        #   panel.curve(priorDistList[[panel.number()]](x), from = max(priorBoundsList[[panel.number()]][1],
        #      current.panel.limits()$xlim[1]), to = min(priorBoundsList[[panel.number()]][2], 
        #      current.panel.limits()$xlim[2]), col = "blue") # need the bounds, from max of lower bound and panel xlim[1], to min of upper bound and panel xlim[2]
        panel.curve(priorDistList[[panel.number()]]
            (x, priorInput[panel.number(), ] ),
            from = max(priorInput[panel.number(), 2] , 
              current.panel.limits()$xlim[1]),
            to = min(priorInput[panel.number(), 3],
              current.panel.limits()$xlim[2]),
            col = "blue")
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
    require(grid, quietly = TRUE, warn.conflicts = FALSE)
    require(lattice, quietly = TRUE, warn.conflicts = FALSE)
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
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^plotDensPOPparsPrior


#findTarget-----------------------------2012-08-14
#  To derive decision tables for moving windows and find
#  the times to achieve recovery with given confidence.
#   Vmat   = matrix of projected B-values (MCMC projections x Year)
#   yrP    = user-specified projection years
#   yrG    = number of years for moving target window (e.g, 90y = 3 YMR generations). Might not work for all possibilities.
#   ratio  = recovery target ratio
#   target = recovery target values (e.g., B0, Bmsy).
#            = B0.MCMC for ratios of B0
#            = Bmsy.MCMC for ratios of Bmsy
#            = Bt.MCMC for moving window
#   conf   = confidence level required
#   plotit = logical to plot the probability of Bt/target
#   retVal = character name of object to return
#    retVal="N", look for the global object "Ttab" (number of years
#     to acheive target)
#    retVal="p.hi" gives global object "Ptab", a list of decision
#     tables where row is the catch option and column is the year
#     Values are probabilities of acheiving target.
#-----------------------------------------------RH
findTarget = function(Vmat, yrU=as.numeric(dimnames(Vmat)[[2]]), yrG=90, ratio=0.5, target=B0.MCMC,
    conf=0.95, plotit=FALSE, retVal="N") {
	
	# oldpar = par(no.readonly=TRUE);  on.exit(par(oldpar))
	yrA    = as.numeric(dimnames(Vmat)[[2]])   # years available
	yrP    = sort(intersect(yrA,yrU))          # years for proj
	yr0    = yrP[1]; yrN = rev(yrP)[1]         # 

	vmat = Vmat[,is.element(dimnames(Vmat)[[2]],as.character(yrP))]             # include only yrP years
	if (is.data.frame(target) || is.matrix(target)) {
		yrM   = yrP - yrG                                                        # moving target years
		yrM1  = intersect(as.numeric(dimnames(target)[[2]]),yrM)                 # available target years from MCMC
		if (length(yrM1)==0) {                                                   # projection not long enough for any overlap with 3 generations
			if (retVal=="N") return(NA)
			else {p.hi = rep(NA,length(yrP)); names(p.hi)=yrP }; return(p.hi) }
		yrMr  = range(yrM1)                                                      # range of years to use from MCMC
		targM = target[,as.character(yrM1)]                                      # target data from MCMC
		yrM2  = setdiff(yrM,yrM1)                                                # missing target years (can occur before and after the MCMC years)
#browser(); return()
		if (length(yrM2)>0) {
			nrow = dim(target)[1]
			if (any(yrM2<yrMr[1])) {
				yrMo  = yrM2[yrM2<yrMr[1]]                                         # years of data older than MCMCs
				ncol  = length(yrMo)
				targ0 = matrix(rep(target[,as.character(yrM1[1])],ncol),
					nrow=nrow, ncol=ncol, dimnames=list(1:nrow,yrMo))               # repeat B0 (first column)
				targM = cbind(as.data.frame(targ0),targM)                          # moving target
			}
			if (any(yrM2>yrMr[2])) {
				yrMn  = yrM2[yrM2>yrMr[2]]                                         # years of data newer than MCMCs
				ncol  = length(yrMn)
				targN = vmat[,as.character(yrMn)]                                  # start using projections
				targM = cbind(targM,targN)                                         # moving target
			}
		}
		rats = vmat/targM                                                        # matrix of ratios Bt/ moving target
	}
	else    # if it's a vector, so no moving window
		rats = apply(vmat,2,function(x,targ){x/targ},targ=target)                # matrix of ratios Bt/ target (B0 or Bmsy)
	p.hi = apply(rats,2,function(x,r){xhi=x>=r; sum(xhi)/length(xhi)},r=ratio)  # vector of probabilities Bt/B0 > target ratio for each year.
	# p.hi can become each row of a decision table (AME checked
	#  the numbers for 0.4 Bmsy match my existing
	#  independent calculations). Need to save this for moving window.
#browser(); return()
	z.hi = p.hi >= conf                                                         # logical: is p.hi >= confidence limit specified

	if (all(z.hi))       yrT = yr0                      # all p.hi exceed the confidence level
	else if (!any(z.hi)) yrT = yrN                      # no  p.hi exceed the confidence level
	else {
		pdif = diff(p.hi)                                # one-year change in trend
		z1 = diff(p.hi)>0                                # logical: trend increasing?
		z2 = c(pdif[-1],FALSE)>0                         # logical: trend one period later increasing?
		z3 = z.hi[-1]                                    # does the probability of equalling or exceeding the target ratio exceed the confidence level?
		z  = z1 & z2 & z3                                # logical: potential years when target reached
		if (!any(z)) yrT = yrN                           # target not reached within the projection period
		else         yrT = as.numeric(names(z)[z][1])    # first year when target reached
	}
	N = yrT - yr0                                       # number of years to reach target
	if (plotit) {
		par(mar=c(4,5,0.5,0.5))
		#ylim = c(0, max(p.hi,ratio))
		ylim = c(min(p.hi,0.5),1)
		plot(yr0:yrN,p.hi,type="n",ylim=ylim,ylab="",mgp=c(2.0,0.75,0))
		lines(yr0:yrN,p.hi,col="grey")
		points(yr0:yrN,p.hi,pch=20,col="orange",cex=1.2)
		mtext(text=expression(p~~frac(B[t],B[Target]) ), side=2, line=1.5, cex=1.5)
		abline(h=conf,v=yrT,col="dodgerblue")
	}
	eval(parse(text=paste("return(",retVal,")",sep=""))) 
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^findTarget


#importProjRec--------------------------2011-08-31
# Imports the projected recruitments
#  (actually what's saved is the N(0,1) random numbers,
#   which for a particular MCMC sample are the same for all the catch strategies),
#   which 'importProj' does not do. Need this for YMR. 4th July 2011.
#----------------------------------------------AME
importProjRec = function (dir, info = "", coda = FALSE, quiet = TRUE) 
{
	get.Policies <- function() {
		if (!quiet) cat("Policies  ")
		Policies <- read.table(paste(dir, "strategy.out", sep = "/"), skip = 1)
		if (!quiet) cat("file...")
		Policies <- unique(as.vector(as.matrix(Policies)))
		if (!quiet) cat("unique...OK\n")
		return(Policies)
	}
	get.Years <- function() {
		if (!quiet) cat("Years...")
		Years <- read.table(paste(dir, "strategy.out", sep = "/"), nrows = 1)
		if (!quiet) cat("file...")
		Years <- unlist(strsplit(as.matrix(Years), "_"))
		if (!quiet) cat("labels...")
		Years <- unique(matrix(Years, nrow = 3)[2, ])
		if (!quiet) cat("unique...OK\n")
		return(Years)
	}
	get.B <- function(Policies, Years) {
		if (!quiet) cat("Biomass   ")
		B <- read.table(paste(dir, "projspbm.out", sep = "/"), header = TRUE)[, -c(1, 2)]
		if (!quiet) cat("file...")
		Blist <- list()
		for (p in 1:length(Policies)) {
			from <- (p - 1) * length(Years) + 1
			to <- p * length(Years)
			Blist[[p]] <- B[, from:to]
			names(Blist[[p]]) <- Years
		}
		names(Blist) <- Policies
		B <- Blist
		if (!quiet) cat("list...OK\n")
		return(B)
	}
	get.Y <- function(Policies, Years) {
		if (!quiet) cat("Landings...")
		Y <- read.table(paste(dir, "procatch.out", sep = "/"), header = TRUE)
		if (!quiet) cat("file...")
		Ylist <- list()
		for (p in 1:length(Policies)) {
			from <- (p - 1) * length(Years) + 1
			to <- p * length(Years)
			Ylist[[p]] <- Y[, from:to]
			names(Ylist[[p]]) <- Years
		}
		names(Ylist) <- Policies
		Y <- Ylist
		if (!quiet) cat("list...OK\n")
		return(Y)
	}
	# AME adding to load in projected recruitment residuals.
	#  Bit different to the others as the file is saved
	#  differently. NOTE that this returns the 'tempdev' values
	#  from Awatea, which are just N(0,1) values, without
	#  multiplying by the sigma_R. Call this epsilon_t
	#  and maybe multiply here by sigmaR, which is
	#	obj$extra$residuals$p_log_RecDev[6]
	get.eps <- function(Policies, Years) {
		if (!quiet) cat("Recruitment...")
		eps <- read.table(paste(dir, "RecRes.out", sep = "/"), header = FALSE, skip=1)
		if (!quiet) cat("file...")
		nRow = dim(eps)[1] / length(Policies)
		epslist <- list()
		for (p in 1:length(Policies)) {
			rows <- (0:(nRow-1)) * length(Policies) + p 
			epslist[[p]] <- sigmaR * eps[rows, ]
			names(epslist[[p]]) <- Years
		}
		names(epslist) <- Policies
		eps <- epslist
		if (!quiet) cat("list...OK\n")
		return(eps)
	}
	files <- paste(dir, c("strategy.out", "projspbm.out", "procatch.out", "recres.out"), sep = "/")
	sapply(files, function(f) if (!file.exists(f)) 
		stop("File ", f, " does not exist. Please check the 'dir' argument.", call. = FALSE))
	if (!quiet) cat("\nParsing files in directory ", dir, ":\n\n", sep = "")
	Policies <- get.Policies()
	Years <- get.Years()
	B <- get.B(Policies, Years)
	Y <- get.Y(Policies, Years)
	eps <- get.eps(Policies, Years)
	if (!quiet) cat("\n")
	output <- list(B = B, Y = Y, eps = eps)
	if (coda) {
		require(coda, quietly = TRUE, warn.conflicts = FALSE)
		output <- lapply(output, function(x) lapply(x, mcmc))
	}
	attr(output, "call") <- match.call()
	attr(output, "info") <- info
	return(output)
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^importProjRec


#srFun----------------------------------2011-08-31
# Stock-recruitment function. From ProjRecCalcs.r
#  To input a vector of spawners in year t-1 and calculate recruits in year t.
#  Output for recruits is vector, each element corresponds to spawners the
#  the year before, so will usually want to shift the output by 1 so that 
#  recruits in year t are based on spawners in year t-1.
#  Can also have each input as a vector (used when calculating a single 
#  year but multiple MCMCs, as in first year of projections is based on 
#  penultimate year of MCMC calcualtions.
#----------------------------------------------AME
srFun = function(spawners, h=h.mpd, R0=R0.mpd, B0=B0.mpd) {
# to input a vector of spawners in year t-1 and calculate recruits in year t 
	4 * h * R0 * spawners / ( ( 1 - h) * B0 + (5 * h - 1) * spawners)
}

