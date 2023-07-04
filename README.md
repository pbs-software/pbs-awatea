## PBSawatea: Tools for running Awatea and visualizing the results ##
&copy; Fisheries and Oceans Canada (2011-2022)

<b>PBSawatea</b> provides an R interface for running ADMB Awatea software, which is a variant of the Coleraine fish population software. The added functionality here includes automation of:

1. reweighting abundance data (by weighting survey and commercial index CVs) and composition data (by weighting effective sample sizes <i>N</i> of proportion-at-ages) before choosing an MPD (mode of the posterior distribution) run,
2. launching an MCMC (Monte Carlo Markoff Chain) simulation,
3. calculating <i>B</i><sub>MSY</sub> (biomass at maximum sustainable yield) and <i>u</i><sub>MSY</sub> (exploitation rate at MSY), and
4. customising Sweave files for individual runs and reweightings from various master Sweave files.

The automation offers substantial time-saving when trying numerous model runs.

<b>PBSawatea</b> requires the R packages <b>PBSmodelling</b>, <b>scape</b>, <b>plotMCMC</b>, <b>xtable</b>, <b>lattice</b>, <b>coda</b>, <b>gplots</b>, and <b>Hmisc</b> (all posted on <a href="https://cran.r-project.org/web/packages/">CRAN</a>).

<b>PBSawatea</b> borrows some functionality from the <b>scape</b> package by adopting the code from a few functions and creating variants. We try to acknowledge the original source wherever possible.

<b>WARNING:</b> The reliance on so many packages is not ideal because each relies on other packages (e.g., Frank Harrell's <b>Hmisc</b> taps into a cascade of dependencies that rivals the machinations of Wizard Wickham). This means that a user likely has to install various additional R packages, which is a nuisance. Below you can see this maze of cross-dependencies in the <b>scape</b> package. Given enough time, we could extract ourselves from this morass by designing our own functions; however, extrication is slow.

<b>PBSawatea</b> depends on and imports <b>PBSmodelling</b>. Additionally <b>PBSawatea</b> imports the following packages:

<b>methods</b>, <b>lattice</b>, <b>scape</b>, <b>plotMCMC</b>, and <b>xtable</b>.

<b>PBSawatea</b> imports specific functions from the following:

<!--- importFrom("coda", "mcmc") --->

<!--- importFrom("gplots", "plotCI") --->

importFrom("graphics", "abline", "arrows", "axis", "barplot", "box", "boxplot", "legend", "lines", "mtext", "pairs", "par", "plot", "points", "polygon", "rect", "segments", "text")

importFrom("grDevices", "boxplot.stats", "colorRampPalette", "dev.cur", "dev.list", "dev.off", "extendrange", "pdf", "png", "postscript", "savePlot", "win.metafile", "windows")

importFrom("Hmisc", "Cbind", "panel.xYplot")

importFrom("stats", "acf", "cor", "lowess", "median", "na.pass", "qqnorm", "quantile", "rnorm", "runif", "sd", "var")

importFrom("utils", "Sweave", "data", "flush.console", "help", "installed.packages", "menu", "read.table", "tail", "write.table")

In the end, the target audience is very small, comprising only those people who actually use Awatea.

<b>PBSawatea</b> is not available on <a href="https://cran.r-project.org/">CRAN</a> (Comprehensive R Archive Network); however, the source code is available on <a href="https://github.com/pbs-software/pbs-awatea">GitHub</a>. Additionally, both source tarball and Windows binary &#8212; built using R-devel and checked via `R CMD check --as-cran` routine on a <b>Windows 7</b> 64-bit system &#8212; are available on <a href="https://drive.google.com/drive/folders/0B2Bkic2Qu5LGOGx1WkRySVYxNFU?usp=sharing">Google Drive</a>. The best way to use this code is to source it directly when performing stock assessments, so think of the repo as a depot.

As with any freely available product, there is no warranty or promise that **PBSawatea** will perform adequately for all circumstances. Additionally, coding errors are possible, and users should contact the package maintainer if bugs are detected.

Maintainer: <a href="mailto:rowan.haigh@dfo-mpo.gc.ca">Rowan Haigh</a>

<p align="right"><img src="DFOlogo_small.jpg" alt="DFO logo" style="height:30px;"></p> 


### PBSawatea dependencies ###
 Packages in <b>bold</b> depend on other packages, packages in <i>italics</i> do not. (updated 2018-10-17)

<hr>

<ul style="margin-left:-1em; padding-bottom: 0;"><!-----start PBSmodelling package----->
<li><b>PBSmodelling</b> imports <b>methods</b>, <b>tcltk</b>, <b>XML</b>
  <ul><!-----start PBSmodelling import----->
  <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
    <ul><!-----start methods import----->
    <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
      <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
    </li></ul><!-----end methods----->
  <li><b>tcltk</b> imports <i>utils</i></li>
  <li><b>XML</b> depends on <b>methods</b>, <i>utils</i>
    <ul><!-----start XML depend----->
    <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
      <ul><!-----start methods import----->
      <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
        <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
      </li></ul><!-----end methods----->
    </li></ul><!-----end XML depend----->
  </li></ul><!-----end PBSmodelling import----->
</li></ul><!-----end PBSmodelling package----->

<hr>

<ul style="margin-left:-1em; padding-bottom: 0;"><!-----start scape package----->
<li><b>scape</b> imports <b>stats</b>, <i>utils</i>, <b>coda</b>, <b>Hmisc</b>, <b>lattice</b>
  <ul><!-----start scape import----->
  <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
    <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
    </li><!-----end stats----->
  <li><b>coda</b> imports <b>lattice</b>
    <ul><!-----start coda import----->
    <li><b>lattice</b> imports <b>grid</b>, <i>grDevices</i>, <i>graphics</i>, <b>stats</b>, <i>utils</i>
      <ul><!-----start lattice import----->
      <li><b>grid</b> imports <i>grDevices</i>, <i>utils</i></li>
      <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
        <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
      </li></ul><!-----end lattice----->
    </li></ul><!-----end coda----->
  <li><b>Hmisc</b> depends on <b>lattice</b>, <b>survival</b>, <b>Formula</b>, <b>ggplot2</b>, <b>htmlTable</b>, <b>viridis</b>, <b>htmltools</b>, <i>base64enc</i>
    <ul><!-----start Hmisc import----->
    <li><b>lattice</b> imports <b>grid</b>, <i>grDevices</i>, <i>graphics</i>, <b>stats</b>, <i>utils</i>
      <ul><!-----start lattice import----->
      <li><b>grid</b> imports <i>grDevices</i>, <i>utils</i></li>
      <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
        <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul></li>
      </li></ul><!-----end lattice----->
    <li><b>survival</b> imports <b>graphics</b>, <b>Matrix</b>, <b>methods</b>, <b>splines</b>, <b>stats</b>, <i>utils</i>
      <ul><!-----start survival import----->
      <li><b>graphics</b> imports <i>grDevices</i></li>
      <li><b>Matrix</b> imports <b>methods</b>, <b>graphics</b>, <b>grid</b>, <b>stats</b>, <i>utils</i>, <b>lattice</b>
        <ul><!-----start Matrix import----->
        <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
          <ul><!-----start methods import----->
          <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
            <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
          </li></ul><!-----end methods----->
        <li><b>graphics</b> imports <i>grDevices</i></li>
        <li><b>grid</b> imports <i>grDevices</i>, <i>utils</i></li>
        <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
          <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul></li>
        <li><b>lattice</b> imports <b>grid</b>, <i>grDevices</i>, <i>graphics</i>, <b>stats</b>, <i>utils</i>
          <ul><!-----start lattice import----->
          <li><b>grid</b> imports <i>grDevices</i>, <i>utils</i></li>
          <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
            <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul></li>
          </li></ul><!-----end lattice----->
        </li></ul><!-----end Matrix import----->
        <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
          <ul><!-----start methods import----->
          <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
            <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
          </li></ul><!-----end methods----->
        <li><b>splines</b> imports <b>graphics</b>, <b>stats</b>
          <ul><!-----start splines import----->
          <li><b>graphics</b> imports <i>grDevices</i></li>
          <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
            <ul><li><b>graphics</b> imports <i>grDevices</i></l></ul></li>
          </li></ul><!-----end splines----->
        <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
          <ul><b>graphics</b> imports <i>grDevices</i></ul></li>
      </li></ul><!-----end survival----->
    <li><b>Formula</b> depends on <b>stats</b>
      <ul><!-----start Formula import----->
      <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
        <ul><li><b>graphics</b> imports <i>grDevices</i></l></ul></li>
      </li></ul><!-----end Formula----->
    <li><b>ggplot2</b> imports <i>digest</i>, <b>grid</b>, <b>gtable</b>, <i>lazyeval</i>, <b>MASS</b>, <b>mgcv</b>, <b>plyr</b>, <b>reshape2</b>, <i>rlang</i>, <b>scales</b>, <b>stats</b>, <b>tibble</b>, <i>viridisLite</i>, <b>withr</b>
      <ul><!-----start ggplot2 import----->
      <li><b>grid</b> imports <i>grDevices</i>, <i>utils</i></li>
      <li><b>gtable</b> imports <b>grid</b>
        <ul><li><b>grid</b> imports <i>grDevices</i>, <i>utils</i></li></ul></li>
      <li><b>MASS</b> depends on <i>grDevices</i>, <b>graphics</b>, <b>stats</b>, <i>utils</i>
        <ul><!-----start MASS depend----->
        <li><b>graphics</b> imports <i>grDevices</i></li>
        <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
          <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul></li>
        </li></ul><!-----end MASS depend----->
      <li><b>MASS</b> imports <b>methods</b>
        <ul><!-----start MASS import----->
        <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
          <ul><!-----start methods import----->
          <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
            <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
          </li></ul><!-----end methods----->
        </li></ul><!-----end MASS import----->
      <li><b>mgcv</b> depends on <b>nlme</b>
        <ul><!-----start mgcv depend----->
        <li><b>nlme</b> imports <b>graphics</b>, <b>stats</b>, <i>utils</i>, <b>lattice</b>
        <li><b>graphics</b> imports <i>grDevices</i></li>
        <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
          <ul><li><b>graphics</b> imports <i>grDevices</i></l></ul></li>
        <li><b>lattice</b> imports <b>grid</b>, <i>grDevices</i>, <i>graphics</i>, <b>stats</b>, <i>utils</i>
          <ul><!-----start lattice import----->
          <li><b>grid</b> imports <i>grDevices</i>, <i>utils</i></li>
          <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
            <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul></li>
          </li></ul><!-----end lattice----->
        </li></ul><!-----end mgcv depend----->
      <li><b>mgcv</b> imports <b>methods</b>, <b>stats</b>, <b>graphics</b>, <b>Matrix</b>
        <ul><!-----start mgcv import----->
        <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
          <ul><!-----start methods import----->
          <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
            <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
          </li></ul><!-----end methods----->
        <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
          <ul><li><b>graphics</b> imports <i>grDevices</i></l></ul></li>
        <li><b>graphics</b> imports <i>grDevices</i></li>
        <li><b>Matrix</b> imports <b>methods</b>, <b>graphics</b>, <b>grid</b>, <b>stats</b>, <i>utils</i>, <b>lattice</b>
          <ul><!-----start Matrix import----->
          <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
            <ul><!-----start methods import----->
            <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
              <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
            </li></ul><!-----end methods----->
          <li><b>graphics</b> imports <i>grDevices</i></li>
          <li><b>grid</b> imports <i>grDevices</i>, <i>utils</i></li>
          <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
            <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul></li>
          <li><b>lattice</b> imports <b>grid</b>, <i>grDevices</i>, <i>graphics</i>, <b>stats</b>, <i>utils</i>
            <ul><!-----start lattice import----->
            <li><b>grid</b> imports <i>grDevices</i>, <i>utils</i></li>
            <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
              <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul></li>
            </li></ul><!-----end lattice----->
          </li></ul><!-----end Matrix import----->
        </li></ul><!-----end mgcv dimport----->
      <li><b>plyr</b> imports <b>Rcpp</b>
        <ul><!-----start plyr import----->
        <li><b>Rcpp</b> imports <b>methods</b>, <i>utils</i>
          <ul><!-----start Rcpp import----->
            <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
              <ul><!-----start methods import----->
              <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
                <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
              </li></ul><!-----end methods----->
          </li></ul><!-----end Rcpp import----->
        </li></ul><!-----end plyr import----->
      <li><b>reshape2</b> imports <b>plyr</b>, <b>Rcpp</b>, <b>stringr</b>
        <ul><!-----start reshape2 import----->
        <li><b>plyr</b> imports <b>Rcpp</b>
          <ul><!-----start plyr import----->
          <li><b>Rcpp</b> imports <b>methods</b>, <i>utils</i>
            <ul><!-----start Rcpp import----->
              <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
                <ul><!-----start methods import----->
                <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
                  <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
                </li></ul><!-----end methods----->
            </li></ul><!-----end Rcpp import----->
          </li></ul><!-----end plyr import----->
        <li><b>Rcpp</b> imports <b>methods</b>, <i>utils</i>
          <ul><!-----start Rcpp import----->
          <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
            <ul><!-----start methods import----->
            <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
              <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
            </li></ul><!-----end methods----->
          </li></ul><!-----end Rcpp import----->
        <li><b>stringr</b> imports <b>glue</b>, <i>magrittr</i>, <b>stringi</b>
          <ul><!-----start stringr import----->
          <li><b>glue</b> imports <b>methods</b>
            <ul><!-----start glue import----->
            <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
              <ul><!-----start methods import----->
              <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
                <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
              </li></ul><!-----end methods----->
            </li></ul><!-----end glue import----->
          <li><b>stringi</b> imports <i>tools</i>, <i>utils</i>, <b>stats</b>
            <ul><!-----start stringi import----->
            <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
            <ul><li><b>graphics</b> imports <i>grDevices</i></l></ul></li>
            </li></ul><!-----end stringi import----->
          </li></ul><!-----end stringr import----->
        </li></ul><!-----end reshape2 import----->
      <li><b>scales</b> imports <i>labeling</i>, <b>munsell</b>, <i>R6</i>, <i>RColorBrewer</i>, <b>Rcpp</b>, <i>viridisLite</i>
        <ul><!-----start scales import----->
        <li><b>munsell</b> imports <b>colorspace</b>, <b>methods</b>
          <ul><!-----start munsell import----->
          <li><b>colorspace</b> imports <b>graphics</b>, <i>grDevices</i>
            <ul><li><b>graphics</b> imports <i>grDevices</i></l></li></ul>
          <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
            <ul><!-----start methods import----->
            <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
              <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
            </li></ul><!-----end methods----->
          </li></ul><!-----end munsell import----->
        <li><b>Rcpp</b> imports <b>methods</b>, <i>utils</i>
          <ul><!-----start Rcpp import----->
          <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
            <ul><!-----start methods import----->
            <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
              <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
            </li></ul><!-----end methods----->
          </li></ul><!-----end Rcpp import----->
        </li></ul><!-----end scales import----->
      <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
        <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
      <li><b>tibble</b> imports <b>cli</b>, <b>crayon</b>, <b>methods</b>, <b>pillar</b>, <i>rlang</i>, <i>utils</i>
        <ul><!-----start tibble import----->
        <li><b>cli</b> imports <b>assertthat</b>, <b>crayon</b>, <b>methods</b>, <i>utils</i>
          <ul><!-----start cli import----->
          <li><b>assertthat</b> imports <i>tools</i></li>
          <li><b>crayon</b> imports <i>grDevices</i>, <b>methods</b>, <i>utils</i>
            <ul><!-----start crayon import----->
            <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
              <ul><!-----start methods import----->
              <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
                <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
              </li></ul><!-----end methods----->
            </li></ul><!-----end crayon import----->
            <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
              <ul><!-----start methods import----->
              <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
                <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
              </li></ul><!-----end methods----->
          </li></ul><!-----end cli import----->
        <li><b>crayon</b> imports <i>grDevices</i>, <b>methods</b>, <i>utils</i>
          <ul><!-----start crayon import----->
          <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
            <ul><!-----start methods import----->
            <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
              <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
            </li></ul><!-----end methods----->
          </li></ul><!-----end crayon import----->
        <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
          <ul><!-----start methods import----->
          <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
            <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
          </li></ul><!-----end methods----->
        <li><b>pillar</b> imports <b>cli</b>, <b>crayon</b>, <i>fansi</i>, <b>methods</b>, <i>rlang</i>, <i>utf8</i>
          <ul><!-----start pillar import----->
          <li><b>cli</b> imports <b>assertthat</b>, <b>crayon</b>, <b>methods</b>, <i>utils</i>
            <ul><!-----start cli import----->
            <li><b>assertthat</b> imports <i>tools</i></li>
            <li><b>crayon</b> imports <i>grDevices</i>, <b>methods</b>, <i>utils</i>
              <ul><!-----start crayon import----->
              <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
                <ul><!-----start methods import----->
                <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
                  <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
                </li></ul><!-----end methods----->
              </li></ul><!-----end crayon import----->
              <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
                <ul><!-----start methods import----->
                <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
                  <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
                </li></ul><!-----end methods----->
            </li></ul><!-----end cli import----->
          <li><b>crayon</b> imports <i>grDevices</i>, <b>methods</b>, <i>utils</i>
            <ul><!-----start crayon import----->
            <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
              <ul><!-----start methods import----->
              <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
                <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
              </li></ul><!-----end methods----->
            </li></ul><!-----end crayon import----->
          <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
            <ul><!-----start methods import----->
            <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
              <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
            </li></ul><!-----end methods----->
          </li></ul><!-----end pillat import----->
        </li></ul><!-----end tibble import----->
      <li><b>withr</b> imports <b>stats</b>, <b>graphics</b>, <i>grDevices</i>
        <ul><!-----start withr import----->
          <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
            <ul><li><b>graphics</b> imports <i>grDevices</i></l></ul></li>
          <li><b>graphics</b> imports <i>grDevices</i></li>
        </li></ul><!-----end withr import----->
      </li></ul><!-----end gglot2 import----->
    <li><b>htmlTable</b> imports <b>stringr</b>, <b>knitr</b>, <i>magrittr</i>, <b>methods</b>, <b>checkmate</b>, <b>htmlwidgets</b>, <b>htmltools</b>, <i>rstudioapi</i>
      <ul><!-----start htmlTable import----->
      <li><b>stringr</b> imports <b>glue</b>, <i>magrittr</i>, <b>stringi</b>
        <ul><!-----start stringr import----->
        <li><b>glue</b> imports <b>methods</b>
          <ul><!-----start glue import----->
          <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
            <ul><!-----start methods import----->
            <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
              <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
            </li></ul><!-----end methods----->
          </li></ul><!-----end glue import----->
        <li><b>stringi</b> imports <i>tools</i>, <i>utils</i>, <b>stats</b>
          <ul><!-----start stringi import----->
          <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
          <ul><li><b>graphics</b> imports <i>grDevices</i></l></ul></li>
          </li></ul><!-----end stringi import----->
        </li></ul><!-----end stringr import----->
      <li><b>knitr</b> imports <b>evaluate</b>, <i>highr</i>, <b>markdown</b>, <b>stringr</b>, <i>yaml</i>, <b>methods</b>, <i>tools</i>
        <ul><!-----start knitr import----->
        <li><b>evaluate</b> imports <b>methods</b>, <b>stringr</b>
          <ul><!-----start evaluate import----->
          <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
            <ul><!-----start methods import----->
            <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
              <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
            </li></ul><!-----end methods----->
          <li><b>stringr</b> imports <b>glue</b>, <i>magrittr</i>, <b>stringi</b>
            <ul><!-----start stringr import----->
            <li><b>glue</b> imports <b>methods</b>
              <ul><!-----start glue import----->
              <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
                <ul><!-----start methods import----->
                <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
                  <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
                </li></ul><!-----end methods----->
              </li></ul><!-----end glue import----->
            <li><b>stringi</b> imports <i>tools</i>, <i>utils</i>, <b>stats</b>
              <ul><!-----start stringi import----->
              <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
              <ul><li><b>graphics</b> imports <i>grDevices</i></l></ul></li>
              </li></ul><!-----end stringi import----->
            </li></ul><!-----end stringr import----->
          </li></ul><!-----end evaluate import----->
        <li><b>markdown</b> imports <i>utils</i>, <b>mime</b>
          <ul><!-----start markdown import----->
          <li><b>mime</b> imports <i>tools</i></li>
          </li></ul><!-----end markdown import----->
        <li><b>stringr</b> imports <b>glue</b>, <i>magrittr</i>, <b>stringi</b>
          <ul><!-----start stringr import----->
          <li><b>glue</b> imports <b>methods</b>
            <ul><!-----start glue import----->
            <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
              <ul><!-----start methods import----->
              <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
                <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
              </li></ul><!-----end methods----->
            </li></ul><!-----end glue import----->
          <li><b>stringi</b> imports <i>tools</i>, <i>utils</i>, <b>stats</b>
            <ul><!-----start stringi import----->
            <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
            <ul><li><b>graphics</b> imports <i>grDevices</i></l></ul></li>
            </li></ul><!-----end stringi import----->
          </li></ul><!-----end stringr import----->
        <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
          <ul><!-----start methods import----->
          <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
            <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
          </li></ul><!-----end methods----->
        </li></ul><!-----end knitr import----->
      <li><b>checkmate</b> imports <b>backports</b>, <i>utils</i>
        <ul><!-----start checkmate import----->
        <li><b>backports</b> imports <i>utils</i></li>
        </li></ul><!-----end checkmate import----->
      <li><b>htmlwidgets</b> imports <i>grDevices</i>, <b>htmltools</b>, <b>jsonlite</b>, <i>yaml</i>
        <ul><!-----start htmlwidgets import----->
        <li><b>htmltools</b> imports <i>utils</i>, <i>digest</i>, <b>Rcpp</b>
          <ul><!-----start htmltools import----->
          <li><b>Rcpp</b> imports <b>methods</b>, <i>utils</i>
            <ul><!-----start Rcpp import----->
            <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
              <ul><!-----start methods import----->
              <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
                <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
              </li></ul><!-----end methods----->
            </li></ul><!-----end Rcpp import----->
          </li></ul><!-----end htmltools import----->
        <li><b>jsonlite</b> depends on <b>methods</b>
          <ul><!-----start jsonlite import----->
          <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
            <ul><!-----start methods import----->
            <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
              <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
            </li></ul><!-----end methods----->
          </li></ul><!-----end jsonlite import----->
        </li></ul><!-----end htmlwidgets import----->
      <li><b>htmltools</b> imports <i>utils</i>, <i>digest</i>, <b>Rcpp</b>
        <ul><!-----start htmltools import----->
        <li><b>Rcpp</b> imports <b>methods</b>, <i>utils</i>
          <ul><!-----start Rcpp import----->
          <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
            <ul><!-----start methods import----->
            <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
              <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
            </li></ul><!-----end methods----->
          </li></ul><!-----end Rcpp import----->
        </li></ul><!-----end htmltools import----->
      </li></ul><!-----end htmlTable import----->
    <li><b>viridis</b> depends on <i>viridisLite</i></li>
    <li><b>viridis</b> imports <b>stats</b>, <b>ggplot2</b>, <b>gridExtra</b>
      <ul><!-----start viridis import----->
      <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
        <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
      <li><b>ggplot2</b> imports <i>digest</i>, <b>grid</b>, <b>gtable</b>, <i>lazyeval</i>, <b>MASS</b>, <b>mgcv</b>, <b>plyr</b>, <b>reshape2</b>, <i>rlang</i>, <b>scales</b>, <b>stats</b>, <b>tibble</b>, <i>viridisLite</i>, <b>withr</b>   (...<b>ggplot2</b> expansion withheld for humane reasons, see expansion above)
      <li><b>gridExtra</b> imports <b>gtable</b>, <b>grid</b>, <i>grDevices</i>, <b>graphics</b>, <i>utils</i>
        <ul><!-----start gridExtra import----->
        <li><b>gtable</b> imports <b>grid</b>
          <ul><!-----start gtable import----->
          <li><b>grid</b> imports <i>grDevices</i>, <i>utils</i></li>
          </li></ul><!-----end gtable import----->
          <li><b>grid</b> imports <i>grDevices</i>, <i>utils</i></li>
          <li><b>graphics</b> imports <i>grDevices</i></li>
        </li></ul><!-----end gridExtra import----->
      </li></ul><!-----end viridis import----->
    <li><b>htmltools</b> imports <i>utils</i>, <i>digest</i>, <b>Rcpp</b>
      <ul><!-----start htmltools import----->
      <li><b>Rcpp</b> imports <b>methods</b>, <i>utils</i>
        <ul><!-----start Rcpp import----->
        <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
          <ul><!-----start methods import----->
          <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
            <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
          </li></ul><!-----end methods----->
        </li></ul><!-----end Rcpp import----->
      </li></ul><!-----end htmltools import----->
    </li></ul><!-----end Hmisc import----->
  <li><b>lattice</b> imports <b>grid</b>, <i>grDevices</i>, <i>graphics</i>, <b>stats</b>, <i>utils</i>
    <ul><!-----start lattice import----->
    <li><b>grid</b> imports <i>grDevices</i>, <i>utils</i></li>
    <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
      <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul></li>
    </li></ul><!-----end lattice----->
  </li></ul><!-----end scape import----->
</li></ul><!-----end scape package----->

<hr>

<ul style="margin-left:-1em; padding-bottom: 0;"><!-----start plotMCMC package----->
<li><b>plotMCMC</b> imports <b>coda</b>, <b>gplots</b>, <b>lattice</b>
  <ul><!-----start plotMCMC import----->
  <li><b>coda</b> imports <b>lattice</b>
    <ul><!-----start coda import----->
    <li><b>lattice</b> imports <b>grid</b>, <i>grDevices</i>, <i>graphics</i>, <b>stats</b>, <i>utils</i>
      <ul><!-----start lattice import----->
      <li><b>grid</b> imports <i>grDevices</i>, <i>utils</i></li>
      <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
        <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul></li>
      </li></ul><!-----end lattice----->
    </li></ul><!-----end coda import----->
  <li><b>gplots</b> imports <b>gtools</b>, <b>gdata</b>, <b>stats</b>, <b>caTools</b>, <b>KernSmooth</b>
    <ul><!-----start gplots import----->
    <li><b>gtools</b> depends on <b>methods</b>, <b>stats</b>, <i>utils</i>
      <ul><!-----start gtools import----->
      <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
        <ul><!-----start methods import----->
        <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
          <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
        </li></ul><!-----end methods----->
        <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
          <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul></li>
      </li></ul><!-----end gtools import----->
    <li><b>gdata</b> imports <b>gtools</b>, <b>stats</b>, <b>methods</b>, <i>utils</i>
      <ul><!-----start gdata import----->
      <li><b>gtools</b> depends on <b>methods</b>, <b>stats</b>, <i>utils</i>
        <ul><!-----start gtools import----->
        <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
          <ul><!-----start methods import----->
          <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
            <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
          </li></ul><!-----end methods----->
          <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
            <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul></li>
        </li></ul><!-----end gtools import----->
      <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
        <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul></li>
      <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
        <ul><!-----start methods import----->
        <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
          <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
        </li></ul><!-----end methods----->
      </li></ul><!-----end gdata import----->
    <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
      <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
    <li><b>catTools</b> imports <i>bitops</i></li>
    <li><b>KernSmooth</b> depends on <b>stats</b>
      <ul><!-----start KernSmooth import----->
      <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
        <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
      </li></ul><!-----end KernSmooth import----->
    </li></ul><!-----end gplots import----->
  <li><b>lattice</b> imports <b>grid</b>, <i>grDevices</i>, <i>graphics</i>, <b>stats</b>, <i>utils</i>
    <ul><!-----start lattice import----->
    <li><b>grid</b> imports <i>grDevices</i>, <i>utils</i></li>
    <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
      <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul></li>
    </li></ul><!-----end lattice----->
  </li></ul><!-----end plotMCMC import----->
</li></ul><!-----end plotMCMC package----->

<hr>

<ul style="margin-left:-1em; padding-bottom: 0;"><!-----start lattice package----->
<li><b>lattice</b> imports  <b>grid</b>, <i>grDevices</i>, <i>graphics</i>, <b>stats</b>, <i>utils</i>
  <ul><!-----start lattice import----->
  <li><b>grid</b> imports <i>grDevices</i>, <i>utils</i></li>
  <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
    <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul></li>
  </li></ul><!-----end lattice import----->
</li></ul><!-----end lattice package----->

<hr>

<ul style="margin-left:-1em; padding-bottom: 0;"><!-----start xtable package----->
<li><b>xtable</b> imports  <b>stats</b>, <i>utils</i>
  <ul><!-----start xtable import----->
  <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
    <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul></li>
  </li></ul><!-----end xtable import----->
</li></ul><!-----end xtable package----->

<hr>

<ul style="margin-left:-1em; padding-bottom: 0;"><!-----start methods package----->
<li><b>methods</b> imports <i>utils</i>, <b>stats</b>
  <ul><!-----start methods import----->
  <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
    <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
   </li></ul><!-----end methods----->
</li></ul><!-----end xtable package----->

<hr>

### PBStools ###

Although <b>PBSawatea</b> does not depend on <b>PBStools</b>, this latter package adds functionality that sometimes makes life easier. <b>PBStools</b> is available from the <a href="https://github.com/pbs-software">PBS Software Repository</a> on GitHub.

<ul style="margin-left:-1em; padding-bottom: 0;"><!-----start PBStools package----->
<li><b>PBStools</b> depends on <i>PBSmapping</i>, <b>PBSmodelling</b>, <i>PBSdata</i>, <b>RODBC</b>
  <ul><!-----start PBStools depend----->
  <li><b>PBSmodelling</b> imports <b>methods</b>, <b>tcltk</b>, <b>XML</b>
    <ul><!-----start PBSmodelling import----->
    <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
      <ul><!-----start methods import----->
      <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
        <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
      </li></ul><!-----end methods----->
    <li><b>tcltk</b> imports <i>utils</i></li>
    <li><b>XML</b> depends on <b>methods</b>, <i>utils</i>
      <ul><!-----start XML depend----->
      <li><b>methods</b> imports <i>utils</i>, <b>stats</b>
        <ul><!-----start methods import----->
        <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
          <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul>
        </li></ul><!-----end methods----->
      </li></ul><!-----end XML depend----->
    </li></ul><!-----end PBSmodelling import----->
  <li><b>RODBC</b> imports <b>stats</b>
    <ul><!-----start RODBC import----->
    <li><b>stats</b> imports <i>utils</i>, <i>grDevices</i>, <b>graphics</b>
      <ul><li><b>graphics</b> imports <i>grDevices</i></li></ul></li>
    </li></ul><!-----end RODBC import----->
  </li></ul><!-----end PBStools depend----->
</li></ul><!-----end PBStools package----->

