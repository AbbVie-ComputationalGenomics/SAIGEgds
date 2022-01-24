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
	fn <- system.file("extdata", "grm1k_10k_snp.gds", package="SAIGEgds")
	gdsfile <- seqOpen(fn)
	on.exit(seqClose(gdsfile))

	# load phenotype
	phenofn <- system.file("extdata", "pheno.txt.gz", package="SAIGEgds")
	pheno <- read.table(phenofn, header=TRUE, as.is=TRUE)

	# fit the null model for binary outcomes
	glmm <- seqFitNullGLMM_SPA(y ~ x1 + x2, pheno, gdsfile)
	# save model
	saveRDS(glmm, file="saige_model.rds", compress="xz")
	# p-value calculation
	seqAssocGLMM_SPA(gdsfile, glmm, mac=4, res.savefn="saige_pval.rds")

	# fit the null model for quantitative outcomes
	glmm <- seqFitNullGLMM_SPA(yy ~ x1 + x2, pheno, gdsfile,
		trait.type="quantitative")
	# save model
	saveRDS(glmm, file="saige_model_quant.rds", compress="xz")
	# p-value calculation
	seqAssocGLMM_SPA(gdsfile, glmm, mac=4, res.savefn="saige_pval_quant.rds")
}


test.saige_fit_null_model <- function()
{
	# 'Rounding' was the default in versions prior to R_3.6.0
	# it is used for reproduction of the results created by R (v3.5.2)
	tryCatch(suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rounding")),
		error=function(e) FALSE)

	# load the previous models
	mod1 <- readRDS(system.file("unitTests", "saige_model.rds",
		package="SAIGEgds"))
	mod2 <- readRDS(system.file("unitTests", "saige_model_quant.rds",
		package="SAIGEgds"))

	# open a GDS file
	fn <- system.file("extdata", "grm1k_10k_snp.gds", package="SAIGEgds")
	gdsfile <- seqOpen(fn)
	on.exit(seqClose(gdsfile))

	# load phenotype
	phenofn <- system.file("extdata", "pheno.txt.gz", package="SAIGEgds")
	pheno <- read.table(phenofn, header=TRUE, as.is=TRUE)

	# fit the null model, check binary outcomes
	glmm <- seqFitNullGLMM_SPA(y ~ x1 + x2, pheno, gdsfile)
	checkEquals(mod1, glmm, "check the SAIGE parameters (binary outcomes)",
		tolerance=1e-4)

	# fit the null model, check quantitative outcomes
	glmm <- seqFitNullGLMM_SPA(yy ~ x1 + x2, pheno, gdsfile,
		trait.type="quantitative")
	checkEquals(mod2, glmm,
		"check the SAIGE parameters (quantitative outcomes)", tolerance=1e-4)
}


test.saige_pval <- function()
{
	# load the previous models and results
	mod1 <- readRDS(system.file("unitTests", "saige_model.rds",
		package="SAIGEgds"))
	mod2 <- readRDS(system.file("unitTests", "saige_model_quant.rds",
		package="SAIGEgds"))
	pval1 <- readRDS(system.file("unitTests", "saige_pval.rds",
		package="SAIGEgds"))
	pval2 <- readRDS(system.file("unitTests", "saige_pval_quant.rds",
		package="SAIGEgds"))

	# open a GDS file
	fn <- system.file("extdata", "grm1k_10k_snp.gds", package="SAIGEgds")
	gdsfile <- seqOpen(fn)
	on.exit(seqClose(gdsfile))

	# p-value calculation, check binary outcomes
	assoc <- seqAssocGLMM_SPA(gdsfile, mod1, mac=4)
	checkEquals(pval1, assoc,
		"check the SAIGE p-value output (binary outcomes)", tolerance=1e-7)

	# p-value calculation, check quantitative outcomes
	assoc <- seqAssocGLMM_SPA(gdsfile, mod2, mac=4)
	checkEquals(pval2, assoc,
		"check the SAIGE p-value output (quantitative outcomes)",
		tolerance=1e-7)
}


test.saige_acta_o <- function()
{
	# load the prefit model
	mod <- readRDS(system.file("unitTests", "saige_model.rds",
		package="SAIGEgds"))

	# open a GDS file
	fn <- system.file("extdata", "grm1k_10k_snp.gds", package="SAIGEgds")
	gdsfile <- seqOpen(fn)
	on.exit(seqClose(gdsfile))

	# get a list of variant units for aggregate tests
	ut <- seqUnitSlidingWindows(gdsfile, win.size=200, win.shift=100)

	# run burden, ACAT-V & ACAT-O
	o <- seqAssocGLMM_spaACAT_O(gdsfile, mod, ut)
	v <- seqAssocGLMM_spaACAT_V(gdsfile, mod, ut)
	b <- seqAssocGLMM_spaBurden(gdsfile, mod, ut)

	# check p-value
	checkEquals(o$pval.b1_1, b$pval.b1_1, "ACAT-O vs Burden, beta(1,1)")
	checkEquals(o$pval.b1_25, b$pval.b1_25, "ACAT-O vs Burden, beta(1,25)")
	checkEquals(o$pval.v1_1, v$pval.v1_1, "ACAT-O vs ACAT-V, beta(1,1)")
	checkEquals(o$pval.v1_25, v$pval.v1_25, "ACAT-O vs ACAT-V, beta(1,25)")
}


test.pACAT <- function()
{
    # R implementation
    ps <- 10^-seq(1, 15, 0.1)
    A1 <- matrix(0, nrow=length(ps), ncol=length(ps))
    for (i in seq_along(ps))
    {
        for (j in seq_along(ps))
        {
            T <- mean(c(tanpi(0.5 - ps[i]), tanpi(0.5 - ps[j])))
            A1[i, j] <- 0.5 - atan(T)/pi
        }
    }

    # C implementation
    A2 <- matrix(0, nrow=length(ps), ncol=length(ps))
    for (i in seq_along(ps))
        for (j in seq_along(ps))
            A2[i, j] <- pACAT(c(ps[i], ps[j]))

    # check
	checkEquals(A1, A2, "R / C implementation of ACAT")
}
