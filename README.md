# metascreen: a modular pipeline for drug sensitivity screens

An R package with a modular library for the design and analysis of high-throughput drug combination screens


code name: metascreen   dev: alpha 0.1.0   author: roberthanes


<img src="https://github.com/Enserink-lab/metascreen/blob/main/doc/metascreen.png?raw=true" width="100%" align="left"></img><br />
<br>


The R-package consists of individual functions that can be used to (a) set-up a screen and generate a dispensing file, (b) to read raw measurements and consolidate different data sets, (c) perform a quality control, (d) analyze, (e) report and visualize the results of a drug combination screen.

#### INSTALLATION

```r
install.packages('devtools')
library(devtools)

# Install R-package from github
devtools::install_github('Enserink-lab/metascreen', build = TRUE, build_opts = c("--no-resave-data", "--no-build-vignettes"))
# Load 'metascreen' package
library(metascreen)

```

A set of files can be found in data/ and used as a reference for building data sets.

