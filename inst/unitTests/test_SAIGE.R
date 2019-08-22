#############################################################
#
# DESCRIPTION: test the SAIGE calculation
#

library(RUnit)
library(SAIGEgds)


# create an R object file for package checking
.create_test_file <- function()
{
	# 'Rounding' was the default in versions prior to R_3.6.0
	# it is used for reproduction of the results created by R (v3.5.2)
	tryCatch(suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rounding")),
		error=function(e) FALSE)

	# open a GDS file
	fn <- system.file("extdata/grm1k_10k_snp.gds", package="SAIGEgds")
	gdsfile <- seqOpen(fn)
	on.exit(seqClose(gdsfile))

	# load phenotype
	phenofn <- system.file("extdata/pheno.txt.gz", package="SAIGEgds")
	pheno <- read.table(phenofn, header=TRUE, as.is=TRUE)

	# fit the null model for binary outcomes
	glmm <- seqFitNullGLMM_SPA(y ~ x1 + x2, pheno, gdsfile)
	# save model
	save(glmm, file="saige_model.rda", compress="xz")
	# p-value calculation
	seqAssocGLMM_SPA(gdsfile, glmm, mac=4, res.savefn="saige_pval.rda")

	# fit the null model for quantitative outcomes
	glmm <- seqFitNullGLMM_SPA(yy ~ x1 + x2, pheno, gdsfile, trait.type="quantitative")
	# save model
	save(glmm, file="saige_model_quant.rda", compress="xz")
	# p-value calculation
	seqAssocGLMM_SPA(gdsfile, glmm, mac=4, res.savefn="saige_pval_quant.rda")
}


test.saige_fit_null_model <- function()
{
	# 'Rounding' was the default in versions prior to R_3.6.0
	# it is used for reproduction of the results created by R (v3.5.2)
	tryCatch(suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rounding")),
		error=function(e) FALSE)

	# load the previous models
	mod1 <- get(load(system.file("unitTests/saige_model.rda", package="SAIGEgds")))
	mod2 <- get(load(system.file("unitTests/saige_model_quant.rda", package="SAIGEgds")))

	# open a GDS file
	fn <- system.file("extdata/grm1k_10k_snp.gds", package="SAIGEgds")
	gdsfile <- seqOpen(fn)
	on.exit(seqClose(gdsfile))

	# load phenotype
	phenofn <- system.file("extdata/pheno.txt.gz", package="SAIGEgds")
	pheno <- read.table(phenofn, header=TRUE, as.is=TRUE)

	# fit the null model, check binary outcomes
	glmm <- seqFitNullGLMM_SPA(y ~ x1 + x2, pheno, gdsfile)
	checkEquals(mod1, glmm, "check the SAIGE parameters (binary outcomes)",
		tolerance=1e-4)

	# fit the null model, check quantitative outcomes
	glmm <- seqFitNullGLMM_SPA(yy ~ x1 + x2, pheno, gdsfile, trait.type="quantitative")
	checkEquals(mod2, glmm, "check the SAIGE parameters (quantitative outcomes)",
		tolerance=1e-4)
}


test.saige_pval <- function()
{
	# load the previous models and results
	mod1 <- get(load(system.file("unitTests/saige_model.rda", package="SAIGEgds")))
	mod2 <- get(load(system.file("unitTests/saige_model_quant.rda", package="SAIGEgds")))
	pval1 <- get(load(system.file("unitTests/saige_pval.rda", package="SAIGEgds")))
	pval2 <- get(load(system.file("unitTests/saige_pval_quant.rda", package="SAIGEgds")))

	# open a GDS file
	fn <- system.file("extdata/grm1k_10k_snp.gds", package="SAIGEgds")
	gdsfile <- seqOpen(fn)
	on.exit(seqClose(gdsfile))

	# p-value calculation, check binary outcomes
	assoc <- seqAssocGLMM_SPA(gdsfile, mod1, mac=4)
	checkEquals(pval1, assoc, "check the SAIGE p-value output (binary outcomes)",
		tolerance=1e-7)

	# p-value calculation, check quantitative outcomes
	assoc <- seqAssocGLMM_SPA(gdsfile, mod2, mac=4)
	checkEquals(pval2, assoc, "check the SAIGE p-value output (quantitative outcomes)",
		tolerance=1e-7)
}
