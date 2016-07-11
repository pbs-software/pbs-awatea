###Introduction###

**PBSawatea** provides an R interface for running ADMB Awatea software, which is a variant of the Coleraine fish population software. The added functionality here includes automation of:

1. reweighting abundance data (by weighting survey and commercial index CVs) and composition data (by weighting effective sample sizes *N* of proportion-at-ages) before choosing an MPD (mode of the posterior density) run,
2. launching an MCMC (Monte Carlo Markoff Chain) simulation,
3. calculating Bmsy (biomass at maximum sustainable yield) and Umsy (exploitation rate at MSY), and
4. customising Sweave files for individual runs and reweightings from various master Sweave files.

All the automation offers enormous time-saving when trying numerous model runs.

**PBSawatea** requires the R packages **PBStools**, **scape**, **plotMCMC**, **xtable**, **lattice**, **coda**, **gplots**, and **Hmisc** (all posted on CRAN except **PBStools**, available from GitHub: <https://github.com/pbs-software/pbs-tools>).

This package borrows some functionality from the **scape** package by adopting the code from a few functions and creating variants. We try to acknowledge the original source wherever possible.

The reliance on so many packages is not ideal because each relies on other packages. This means that a user likely has to install various additional R packages, which is a nuisance. Also, there appears to be redundant cross-dependencies (see below); however, **PBSawatea** only imports specific functions from the following:

importFrom("coda", "mcmc")

importFrom("gplots", "plotCI")

importFrom("Hmisc", "Cbind", "panel.xYplot")

In the end, the target audience is very small, comprising only those people who actually use Awatea.

###PBSawatea dependencies###

**PBStools** depends on **PBSmapping**, **PBSmodelling**, **PBSdata**, **RODBC**<br>
&emsp;**PBSmodelling** imports *methods*, *tcltk*, *XML*<br>
&emsp;**RODBC** imports *stats*<br>

**scape** imports **coda**, **Hmisc**, **lattice**<br>
&emsp; **coda** imports **lattice**<br>
&emsp;&emsp; **lattice** imports *grid*, *grDevices*, *graphics*, *stats*, *utils*<br>
&emsp; **Hmisc** depends on **lattice**, **survival**, **Formula**, **ggplot2**<br>
&emsp;&emsp; **lattice** imports *grid*, *grDevices*, *graphics*, *stats*, *utils*<br>
&emsp;&emsp; **survival** imports *graphics*, **Matrix**, *methods*, *splines*, *stats*, *utils*<br>
&emsp;&emsp;&emsp; **Matrix** imports *methods*, *graphics*, *grid*, *stats*, *utils*, **lattice**<br>
&emsp; **Hmisc** imports *methods*, **latticeExtra**, **cluster**, **rpart**, **nnet**, **acepack**, **foreign**, **gtable**, *grid*, **gridExtra**, **data.table**<br>
&emsp;&emsp; **latticeExtra** depends on **lattice**, **RColorBrewer**<br>
&emsp;&emsp; **cluster** imports *graphics*, *grDevices*, *stats*, *utils*<br>
&emsp;&emsp; **rpart** depends on *graphics*, *stats*, *grDevices*<br>
&emsp;&emsp; **nnet** depends on *stats*, *utils*<br>
&emsp;&emsp; **foreign** imports *methods*, *utils*, *stats*<br>
&emsp;&emsp; **gtable** imports *grid*<br>
&emsp;&emsp; **gridExtra** imports **gtable**, *grid*, *grDevices*, *graphics*, *utils*<br>
&emsp;&emsp; **data.table** imports *methods*, **chron**<br>
&emsp;&emsp;&emsp; **chron** imports *graphics*, *stats*<br>
&emsp; **lattice** imports *grid*, *grDevices*, *graphics*, *stats*, *utils*<br>

**plotMCMC** imports **coda**, **gplots**, **lattice**<br>
&emsp; **coda** imports **lattice**<br>
&emsp; **gplots** imports **gtools**, **gdata**, *stats*, **caTools**, **KernSmooth**<br>
&emsp;&emsp; **gdata** imports **gtools**, *stats*, *methods*, *utils*<br>
&emsp;&emsp; **caTools** imports **bitops**<br>
&emsp;&emsp; **KernSmooth** depends on *stats*<br>
&emsp; **lattice** imports *grid*, *grDevices*, *graphics*, *stats*, *utils*<br>

**xtable** imports *stats*, *utils*<br>

**lattice** imports *grid*, *grDevices*, *graphics*, *stats*, *utils*<br>

**coda** imports **lattice**<br>
&emsp; **lattice** imports *grid*, *grDevices*, *graphics*, *stats*, *utils*<br>

**gplots** imports **gtools**, **gdata**, *stats*, **caTools**, **KernSmooth**<br>
&emsp; **gdata** imports **gtools**, *stats*, *methods*, *utils*<br>
&emsp; **caTools** imports **bitops**<br>
&emsp; **KernSmooth** depends on *stats*<br>

**Hmisc** depends on **lattice**, **survival**, **Formula**, **ggplot2**<br>
&emsp; **lattice** imports *grid*, *grDevices*, *graphics*, *stats*, *utils*<br>
&emsp; **survival** imports *graphics*, **Matrix**, *methods*, *splines*, *stats*, *utils*<br>
&emsp;&emsp; **Matrix** imports *methods*, *graphics*, *grid*, *stats*, *utils*, **lattice**<br>
**Hmisc** imports *methods*, **latticeExtra**, **cluster**, **rpart**, **nnet**, **acepack**, **foreign**, **gtable**, *grid*, **gridExtra**, **data.table**<br>
&emsp; **latticeExtra** depends on **lattice**, **RColorBrewer**<br>
&emsp; **cluster** imports *graphics*, *grDevices*, *stats*, *utils*<br>
&emsp; **rpart** depends on *graphics*, *stats*, *grDevices*<br>
&emsp; **nnet** depends on *stats*, *utils*<br>
&emsp; **foreign** imports *methods*, *utils*, *stats*<br>
&emsp; **gtable** imports *grid*<br>
&emsp; **gridExtra** imports **gtable**, *grid*, *grDevices*, *graphics*, *utils*<br>
&emsp; **data.table** imports *methods*, **chron**<br>
&emsp;&emsp; **chron** imports *graphics*, *stats*<br>
