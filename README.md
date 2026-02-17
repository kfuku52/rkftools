# Overview
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
