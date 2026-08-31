DEV_LIB ?= .local/R-library
export RKFTOOLS_DEV_LIBRARY := $(abspath $(DEV_LIB))
export R_LIBS := $(abspath $(DEV_LIB))$(if $(R_LIBS),:$(R_LIBS))
BENCHMARK_ARGS ?=

.PHONY: setup setup-minimal test coverage check check-full check-as-cran document figures benchmark release-check

setup:
	Rscript tools/setup.R

setup-minimal:
	Rscript tools/setup.R --minimal

test:
	Rscript tools/test.R

coverage:
	Rscript tools/coverage.R

check:
	R CMD build --no-build-vignettes .
	_R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual rkftools_$$(Rscript -e 'cat(read.dcf("DESCRIPTION")[1,"Version"])').tar.gz

check-full:
	R CMD build --no-build-vignettes .
	_R_CHECK_FORCE_SUGGESTS_=true R CMD check --no-manual rkftools_$$(Rscript -e 'cat(read.dcf("DESCRIPTION")[1,"Version"])').tar.gz

check-as-cran:
	R CMD build --no-build-vignettes .
	_R_CHECK_FORCE_SUGGESTS_=true R CMD check --as-cran --no-manual rkftools_$$(Rscript -e 'cat(read.dcf("DESCRIPTION")[1,"Version"])').tar.gz

document:
	Rscript -e 'roxygen2::roxygenise(".")'

figures:
	Rscript tools/generate-doc-figures.R

benchmark:
	Rscript tools/benchmark.R $(BENCHMARK_ARGS)

release-check:
	Rscript tools/release-check.R
