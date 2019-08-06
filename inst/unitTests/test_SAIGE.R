#############################################################
#
# DESCRIPTION: test the SAIGE calculation
#

library(RUnit)
library(SAIGEgds)


# create an R object file for package checking
.create_test_file <- function()
{
	# open a GDS file
	fn <- system.file("extdata/grm1k_10k_snp.gds", package="SAIGEgds")
	gdsfile <- seqOpen(fn)
	on.exit(seqClose(gdsfile))

	# load phenotype
	phenofn <- system.file("extdata/pheno.txt.gz", package="SAIGEgds")
	pheno <- read.table(phenofn, header=TRUE, as.is=TRUE)

	# fit the null model
	glmm <- seqFitNullGLMM_SPA(y ~ x1 + x2, pheno, gdsfile)
	# save model
	save(glmm, file="saige_model.rda", compress="xz")

	# p-value calculation
	seqAssocGLMM_SPA(gdsfile, glmm, mac=4, res.savefn="saige_pval.rda")
}


test.saige_fit_null_model <- function()
{
	# 'Rounding' was the default in versions prior to R_3.6.0
	# it is used for reproduction of the results created by R (v3.5.2)
	suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rounding"))

	# load the previous model
	mod <- get(load(system.file("unitTests/saige_model.rda", package="SAIGEgds")))

	# open a GDS file
	fn <- system.file("extdata/grm1k_10k_snp.gds", package="SAIGEgds")
	gdsfile <- seqOpen(fn)
	on.exit(seqClose(gdsfile))

	# load phenotype
	phenofn <- system.file("extdata/pheno.txt.gz", package="SAIGEgds")
	pheno <- read.table(phenofn, header=TRUE, as.is=TRUE)

	# fit the null model
	glmm <- seqFitNullGLMM_SPA(y ~ x1 + x2, pheno, gdsfile)

	# check
	checkEquals(mod, glmm, "check the SAIGE parameters", tolerance=1e-4)
}


test.saige_pval <- function()
{
	# load the previous model
	mod <- get(load(system.file("unitTests/saige_model.rda", package="SAIGEgds")))
	pval <- get(load(system.file("unitTests/saige_pval.rda", package="SAIGEgds")))

	# open a GDS file
	fn <- system.file("extdata/grm1k_10k_snp.gds", package="SAIGEgds")
	gdsfile <- seqOpen(fn)
	on.exit(seqClose(gdsfile))

	# p-value calculation
	assoc <- seqAssocGLMM_SPA(gdsfile, mod, mac=4)

	# check
	checkEquals(pval, assoc, "check the SAIGE p-value output", tolerance=1e-7)
}
