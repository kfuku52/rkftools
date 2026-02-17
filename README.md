# Overview
[![R-CMD-check](https://github.com/kfuku52/rkftools/actions/workflows/r-cmd-check.yaml/badge.svg)](https://github.com/kfuku52/rkftools/actions/workflows/r-cmd-check.yaml)
[![Version](https://img.shields.io/badge/version-0.1.0-informational)](https://github.com/kfuku52/rkftools)
[![R](https://img.shields.io/badge/R-%3E%3D%203.4.1-276DC3?logo=r)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/license-BSD_3--Clause-blue.svg)](LICENSE)
[![Last commit](https://img.shields.io/github/last-commit/kfuku52/rkftools)](https://github.com/kfuku52/rkftools/commits/master)

This R package contains various tools to handle data in evolutionary biology.

# Dependency
Required:
* [R](https://www.r-project.org/)
* [ape](http://ape-package.ird.fr/)
* [phytools](https://github.com/liamrevell/phytools)

Optional (used by specific functions/modes):
* [PhylogeneticEM](https://github.com/pbastide/PhylogeneticEM)
* [Rphylopars](https://github.com/ericgoolsby/Rphylopars)

# Programs that depend on rkftools
* [RADTE](https://github.com/kfuku52/RADTE)

# Programs that rkftools can process its outputs
* [ape](http://ape-package.ird.fr/)
* [l1ou](https://github.com/khabbazian/l1ou) 
* [PhylogeneticEM](https://github.com/pbastide/PhylogeneticEM)
* [NOTUNG](http://www.cs.cmu.edu/~durand/Notung/)

# Installation
```
install.packages("devtools")
devtools::install_github(repo="kfuku52/rkftools", ref="master")
```

# Parallel tuning
`MAD_parallel()` automatically uses available CPU cores when `ncpu` is omitted.

To cap cores globally:
```r
options(rkftools.max_cores = 8)
```
