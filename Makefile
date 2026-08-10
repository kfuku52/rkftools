.PHONY: test coverage check check-as-cran document benchmark release-check

test:
	Rscript -e 'testthat::test_local(".")'

coverage:
	Rscript tools/coverage.R

check:
	R CMD build --no-build-vignettes .
	_R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual rkftools_$$(Rscript -e 'cat(read.dcf("DESCRIPTION")[1,"Version"])').tar.gz

check-as-cran:
	R CMD build --no-build-vignettes .
	_R_CHECK_FORCE_SUGGESTS_=false R CMD check --as-cran --no-manual rkftools_$$(Rscript -e 'cat(read.dcf("DESCRIPTION")[1,"Version"])').tar.gz

document:
	Rscript -e 'roxygen2::roxygenise(".")'

benchmark:
	Rscript tools/benchmark.R

release-check:
	Rscript tools/release-check.R
