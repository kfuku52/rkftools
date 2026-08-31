# Overview
[![R-CMD-check](https://github.com/kfuku52/rkftools/actions/workflows/r-cmd-check.yaml/badge.svg)](https://github.com/kfuku52/rkftools/actions/workflows/r-cmd-check.yaml)
[![Version](https://img.shields.io/badge/version-0.1.12-informational)](https://github.com/kfuku52/rkftools)
[![R](https://img.shields.io/badge/R-%3E%3D%204.1.0-276DC3?logo=r)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/license-MIT-yellow.svg)](LICENSE.md)
[![Last commit](https://img.shields.io/github/last-commit/kfuku52/rkftools)](https://github.com/kfuku52/rkftools/commits/master)
[![GitHub issues](https://img.shields.io/github/issues/kfuku52/rkftools)](https://github.com/kfuku52/rkftools/issues)
[![GitHub pull requests](https://img.shields.io/github/issues-pr/kfuku52/rkftools)](https://github.com/kfuku52/rkftools/pulls)
[![GitHub stars](https://img.shields.io/github/stars/kfuku52/rkftools)](https://github.com/kfuku52/rkftools)
[![GitHub forks](https://img.shields.io/github/forks/kfuku52/rkftools)](https://github.com/kfuku52/rkftools/network/members)

This R package contains various tools to handle data in evolutionary biology.

# Dependency
Required:
* [R (>= 4.1.0)](https://www.r-project.org/)
* [ape](https://github.com/emmanuelparadis/ape)
* [phytools](https://github.com/liamrevell/phytools)

Optional (used by specific functions/modes):
* [PhylogeneticEM](https://github.com/pbastide/PhylogeneticEM)
* [Rphylopars](https://github.com/ericgoolsby/Rphylopars)

# Supported upstream outputs
rkftools includes helpers for objects and outputs produced by:
* [ape](https://github.com/emmanuelparadis/ape) `phylo` objects (`table2phylo()`, `phylo2table()`, and tree utilities)
* [l1ou](https://github.com/khabbazian/l1ou) outputs (`get_tree_table()`, `get_regime_table()`, `get_leaf_table()`, and `get_bootstrap_table()`)
* [PhylogeneticEM](https://github.com/pbastide/PhylogeneticEM) outputs (`get_tree_table()`, `get_regime_table()`, and `get_leaf_table()`)
* [NOTUNG](https://www.cs.cmu.edu/~durand/Notung/) parsable output (`read_notung_parsable()`)

# Installation
```r
install.packages("remotes")
remotes::install_github("kfuku52/rkftools", ref = "master")
```

# Documentation

- [Usage and examples](docs/usage.md): branch tables, species-overlap scoring,
  parallel settings, and transformation contracts.
- [Development](docs/development.md): isolated setup, all tests, optional backend
  fixtures, compatibility checks, and release validation.
- [Performance](docs/performance.md): reproducible before/after measurements.
- [Changelog](NEWS.md)
