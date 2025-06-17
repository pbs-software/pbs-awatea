## PBSawatea: Tools for running Awatea and visualizing the results ##
&copy; Fisheries and Oceans Canada (2011-2025)

<b>PBSawatea</b> provides an R interface for running ADMB Awatea software, which is a variant of the Coleraine fish population software. The added functionality here includes automation of:

1. reweighting abundance data (by weighting survey and commercial index CVs) and composition data (by weighting effective sample sizes <i>N</i> of proportion-at-ages) before choosing an MPD (mode of the posterior distribution) run,
2. launching an MCMC (Monte Carlo Markoff Chain) simulation,
3. calculating <i>B</i><sub>MSY</sub> (biomass at maximum sustainable yield) and <i>u</i><sub>MSY</sub> (exploitation rate at MSY), and
4. customising Sweave files for individual runs and reweightings from various master Sweave files.

The automation offers substantial time-saving when trying numerous model runs.

<b>PBSawatea</b> requires the R packages <b>PBSmodelling</b>, <b>scape</b>, <b>plotMCMC</b>, <b>xtable</b>, <b>lattice</b>, <b>coda</b>, <b>gplots</b>, and <b>Hmisc</b> (all posted on <a href="https://cran.r-project.org/web/packages/">CRAN</a>).

<font color="red"><h3>Background</h3></font>

<b>PBSawatea</b> borrows some functionality from the <b>scape</b> package by adopting the code from a few functions and creating variants. We try to acknowledge the original source wherever possible.

<b>WARNING:</b> The reliance on so many packages is not ideal because each relies on other packages (e.g., Frank Harrell's <b>Hmisc</b> taps into a cascade of dependencies that rivals the machinations of Wizard Wickham). This means that a user likely has to install various additional R packages, which is a nuisance. Below you can see this maze of cross-dependencies in the <b>scape</b> package. Given enough time, we could extract ourselves from this morass by designing our own functions; however, extrication is slow.

<b>PBSawatea</b> depends on and imports <b>PBSmodelling</b>. Additionally <b>PBSawatea</b> imports the following packages:

<b>methods</b>, <b>lattice</b>, <b>scape</b>, <b>plotMCMC</b>, and <b>xtable</b>.

<b>PBSawatea</b> imports specific functions from the following:

import(methods, PBSmodelling, lattice, scape, plotMCMC, xtable)

importFrom("graphics", "abline", "arrows", "axis", "barplot", "box", "boxplot", "contour", "frame", "legend", "lines", "mtext", "pairs", "par", "plot", "points", "polygon", "rect", "segments", "text")

importFrom("grDevices", "boxplot.stats", "colorRampPalette", "dev.cur", "dev.list", "dev.off", "extendrange", "pdf", "png", "postscript", "savePlot", "win.metafile", "windows")

importFrom("Hmisc", "Cbind", "panel.xYplot")

importFrom("stats", "acf", "cor", "fitted", "lowess", "median", "na.pass", "qqnorm", "quantile", "residuals", "rnorm", "runif", "sd", "smooth.spline", "var")

importFrom("utils", "Sweave", "data", "flush.console", "help", "installed.packages", "menu", "read.table", "tail", "write.table")

importFrom("akima", "interp")

importFrom("PBStools", "changeLangOpts", "eop", ".flush.cat", "findRC", "inWord", ".su", "scaleVec")

In the end, the target audience is very small, comprising only those people who actually use Awatea. The platform was last used in 2024 for an update of the model used in the 2019 Bocaccio stock assessment, primarliy to follow up on the huge 2016 recruitment event. 

<b>PBSawatea</b> is not available on <a href="https://cran.r-project.org/">CRAN</a> (Comprehensive R Archive Network); however, the source code is available on <a href="https://github.com/pbs-software/pbs-awatea">GitHub</a>. Most versions undergo a CRAN-like check, via `R CMD check` on a <b>Windows</b> 64-bit system. The package has fallen out of use because the Offshore Rockfish Program switched to using Stock Synthesis 3 in 2021 (see <a href="https://github.com/pbs-software/pbs-synth">PBSsynth</a>). Regardless, the best way to use this code is to source it directly when performing stock assessments, so think of the repo as a depot.

<font color="red"><h3>Disclaimer</h3></font>

"Fisheries and Oceans Canada (DFO) GitHub project code is provided on an 'as is' basis and the user assumes responsibility for its use. DFO relinquishes control of the information and assumes no responsibility to protect the integrity, confidentiality, or availability of the information. Any claims against DFO stemming from the use of its GitHub project will be governed by all applicable Canadian Federal laws. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favouring by DFO. The Fisheries and Oceans Canada seal and logo, or the seal and logo of a DFO bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DFO or the Canadian Government.”

As with any freely available product, there is no warranty or promise that **PBSawatea** will perform adequately for all circumstances. Additionally, coding errors are possible, and users should contact the package maintainer if bugs are detected.

Maintainer: <a href="mailto:rowan.haigh@dfo-mpo.gc.ca">Rowan Haigh</a>

<p align="right"><img src="DFOlogo_small.jpg" alt="DFO logo" style="height:30px;"></p> 

