**PBSawatea** provides an R interface for running ADMB Awatea software, which is a variant of the Coleraine fish population software. The added functionality here includes automation of:

1. reweighting abundance data (by weighting survey and commercial index CVs) and composition data (by weighting effective sample sizes *N* of proportion-at-ages) before choosing an MPD (mode of the posterior density) run,
2. launching an MCMC (Monte Carlo Markoff Chain) simulation,
3. calculating Bmsy (biomass at maximum sustainable yield) and Umsy (exploitation rate at MSY), and
4. customising Sweave files for individual runs and reweightings from various master Sweave files.

All the automation offers enormous time-saving when trying numerous model runs.

**PBSawatea** requires the R packages **PBSmodelling**, **scape**, and **scapeMCMC** (all posted on CRAN). It also borrows some functionality from the scape packages by adopting the code from a few functions and creating variants herein. We try to acknowledge the original source wherever possible.
