#######################################################################
#
# Package name: SAIGEgds
#
# Description:
#     Scalable and accurate implementation of generalized mixed models
# using GDS files
#
# Copyright (C) 2021    Xiuwen Zheng / AbbVie-ComputationalGenomics
# License: GPL-3
#


.minor_allele_geno <- function(geno)
{
    if (anyNA(geno))
    {
        af <- mean(geno, na.rm=TRUE)
        if (is.na(af)) af <- 0L
        geno[is.na(geno)] <- af
    }
    af <- mean(geno)
    if (af > 1L) geno <- 2L - geno
    geno
}

.show_coeff <- function(fit0, indent="", initial=TRUE)
{
    if (initial)
        .cat(indent, "Initial fixed-effect coefficients:")
    else
        .cat(indent, "Fixed-effect coefficients:")
    v <- as.data.frame(t(fit0$coefficients))
    rownames(v) <- paste0("   ", indent)
    print(v)
    invisible()
}


#######################################################################
# Fit the null model
#

seqGLMM_GxG_spa <- function(formula, data, gds_grm, gds_assoc, snp_pair,
    trait.type=c("binary", "quantitative"), sample.col="sample.id", maf=0.005,
    missing.rate=0.01, max.num.snp=1000000L, variant.id=NULL, inv.norm=TRUE,
    X.transform=TRUE, tol=0.02, maxiter=20L, nrun=30L, tolPCG=1e-5,
    maxiterPCG=500L, tau.init=c(0,0), use_approx_tau=FALSE, glm_threshold=FALSE,
    traceCVcutoff=0.0025, ratioCVcutoff=0.001, geno.sparse=TRUE, num.thread=1L,
    model.savefn="", seed=200L, fork.loading=FALSE, verbose=TRUE,
    verbose.detail=TRUE)
{
    stopifnot(inherits(formula, "formula"))
    stopifnot(is.data.frame(data))
    stopifnot(is.character(gds_grm) | inherits(gds_grm, "SeqVarGDSClass"))
    stopifnot(is.null(gds_assoc) | is.character(gds_assoc) |
        inherits(gds_assoc, "SeqVarGDSClass") | is.matrix(gds_assoc))
    if (is.matrix(gds_assoc))
        stopifnot(is.matrix(gds_assoc), is.numeric(gds_assoc))
    stopifnot(is.data.frame(snp_pair), ncol(snp_pair)>=2L, nrow(snp_pair)>0L)
    trait.type <- match.arg(trait.type)
    stopifnot(is.character(sample.col), length(sample.col)==1L,
        !is.na(sample.col))
    stopifnot(is.numeric(maf), length(maf)==1L)
    stopifnot(is.numeric(missing.rate), length(missing.rate)==1L)
    stopifnot(is.numeric(max.num.snp), length(max.num.snp)==1L)
    stopifnot(is.logical(inv.norm), length(inv.norm)==1L)
    stopifnot(is.logical(X.transform), length(X.transform)==1L)
    stopifnot(is.numeric(tol), length(tol)==1L)
    stopifnot(is.numeric(maxiter), length(maxiter)==1L)
    stopifnot(is.numeric(nrun), length(nrun)==1L)
    stopifnot(is.numeric(tolPCG), length(tolPCG)==1L)
    stopifnot(is.numeric(maxiterPCG), length(maxiterPCG)==1L)
    stopifnot(is.numeric(tau.init), length(tau.init)==2L)
    stopifnot(is.logical(use_approx_tau), length(use_approx_tau)==1L,
        !is.na(use_approx_tau))
    stopifnot(is.logical(glm_threshold) | is.numeric(glm_threshold),
        length(glm_threshold)==1L)
    stopifnot(is.numeric(traceCVcutoff), length(traceCVcutoff)==1L)
    stopifnot(is.numeric(ratioCVcutoff), length(ratioCVcutoff)==1L)
    stopifnot(is.logical(geno.sparse), length(geno.sparse)==1L)
    stopifnot(is.numeric(num.thread), length(num.thread)==1L)
    stopifnot(is.character(model.savefn), length(model.savefn)==1L)
    stopifnot(is.numeric(seed), length(seed)==1L, is.finite(seed))
    stopifnot(is.logical(fork.loading), length(fork.loading)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    stopifnot(is.logical(verbose.detail), length(verbose.detail)==1L)
    if (!verbose) verbose.detail <- FALSE

    if (verbose)
    {
        .cat(.crayon_inverse("SAIGE association analysis on the GxG interaction:"))
        .cat(.crayon_underline(date()))
    }

    # check snp_pair
    if (anyNA(snp_pair))
        stop("'snp_pair' should not have missing values.")
    if (any(snp_pair[,1L] == snp_pair[,2L]))
        stop("'snp_pair' should not have the same variant in a pair.")

    if (is.character(gds_grm))
    {
        if (verbose)
            .cat("Open ", sQuote(gds_grm))
        gds_grm <- seqOpen(gds_grm, allow.duplicate=TRUE)
        on.exit(seqClose(gds_grm))
    } else {
        # save the filter on GDS file
        seqSetFilter(gds_grm, action="push", verbose=FALSE)
        on.exit(seqSetFilter(gds_grm, action="pop", verbose=FALSE))
    }

    if (is.character(gds_assoc))
    {
        if (verbose)
            .cat("Open ", sQuote(gds_assoc))
        gds_assoc <- seqOpen(gds_assoc, allow.duplicate=TRUE)
        on.exit(seqClose(gds_assoc), add=TRUE)
    } else if (is.null(gds_assoc))
    {
        gds_assoc <- gds_grm
    } else if (inherits(gds_assoc, "SeqVarGDSClass")) {
        # save the filter on GDS file
        seqSetFilter(gds_assoc, action="push", verbose=FALSE)
        on.exit(seqSetFilter(gds_assoc, action="pop", verbose=FALSE), add=TRUE)
    } else {
        # should be a matrix
        if (is.null(rownames(gds_assoc)))
            stop("rownames(gds_assoc) should be sample IDs, if gds_assoc is a matrix.")
        if (is.null(colnames(gds_assoc)))
            colnames(gds_assoc) <- seq_len(ncol(gds_assoc))
    }

    # show warnings immediately
    saveopt <- options(warn=1L)
    on.exit(options(warn=saveopt$warn), add=TRUE)
    set.seed(seed)

    # variables in the formula
    vars <- all.vars(formula)
    phenovar <- all.vars(formula)[1L]
    y <- data[[phenovar]]
    if (is.null(y))
        stop("There is no '", phenovar, "' in the input data frame.")
    if (!is.factor(y) && !is.numeric(y) && !is.logical(y))
        stop("The response variable should be numeric or a factor.")

    # check sample id
    if (sample.col %in% vars)
        stop(sprintf("'%s' should not be in the formula.", sample.col))
    if (!(sample.col %in% colnames(data)))
        stop(sprintf("'%s' should be one of the columns in 'data'.", sample.col))
    if (is.factor(data[[sample.col]]))
        stop(sprintf("'%s' should not be a factor variable.", sample.col))
    if (any(is.na(data[[sample.col]])))
        stop(sprintf("'%s' should not have any missing value.", sample.col))
    if (anyDuplicated(data[[sample.col]]))
        stop(sprintf("'%s' in data should be unique.", sample.col))

    # remove missing values
    data <- data[, c(sample.col, vars)]
    data <- na.omit(data)
    data <- droplevels(data)
    seqResetFilter(gds_grm, sample=TRUE, verbose=FALSE)
    sid <- seqGetData(gds_grm, "sample.id")
    i <- match(sid, data[[sample.col]])
    i <- i[!is.na(i)]
    data <- data[i, ]
    if (nrow(data) <= 0L)
        stop("No common sample.id between 'data' and the GDS file.")
    sid <- data[[sample.col]]
    seqSetFilter(gds_grm, sample.id=sid, verbose=FALSE)

    if (is.null(variant.id))
    {
        # filters of maf, mac, missing.rate
        if (verbose)
            cat("Filtering variants:\n")
        seqSetFilterCond(gds_grm, maf=maf, mac=NaN, missing.rate=missing.rate,
            parallel=num.thread, .progress=TRUE, verbose=FALSE)
    } else {
        seqSetFilter(gds_grm, variant.id=variant.id, verbose=FALSE)
    }

    # the max mumber of SNPs
    v <- seqGetFilter(gds_grm)$variant.sel
    n <- sum(v, na.rm=TRUE)
    if (max.num.snp>0L && n>max.num.snp)
    {
        set.seed(seed)
        seqSetFilter(gds_grm, variant.sel=sample(which(v), max.num.snp),
            verbose=FALSE)
    }

    # get the number of samples / variants
    dm <- seqSummary(gds_grm, "genotype", verbose=FALSE)$seldim
    n_samp <- dm[2L]
    n_var  <- dm[3L]
    if (verbose)
    {
        cat("Fit the null model:", format(formula), "+ var(GRM)\n")
        .cat("    # of samples: ", .pretty(n_samp))
        cat("    # of variants:", .pretty(n_var))
        if (n > max.num.snp)
            cat(" (randomly selected from ", .pretty(n), ")", sep="")
        cat("\n")
    }

    # set the number of internal threads
    if (is.na(num.thread) || num.thread < 1L)
        num.thread <- 1L
    # setThreadOptions(num.thread)  # no need to call setThreadOptions()
    if (verbose)
        .cat("    using ", num.thread, " thread", ifelse(num.thread>1L, "s", ""))

    # load SNP genotypes
    if (verbose)
        cat("Start loading SNP genotypes:\n")
    nfork <- 1L
    if (SeqArray:::.IsForking(num.thread) && isTRUE(fork.loading))
        nfork <- num.thread
    if (isTRUE(geno.sparse))
    {
        # sparse genotypes
        # integer genotypes or numeric dosages
        if (!is.null(index.gdsn(gds_grm, "genotype/data", silent=TRUE)))
        {
            nm <- "$dosage_alt"
        } else if (!is.null(index.gdsn(gds_grm, "annotation/format/DS/data", silent=TRUE)))
        {
            nm <- "annotation/format/DS"
            if (verbose) cat("    using 'annotation/format/DS'\n")
        } else
            stop("'genotype' and 'annotation/format/DS' are not available.")
        # internal buffer
        buffer <- integer(n_samp + 4L)
        .cfunction2("saige_init_sparse")(n_samp, buffer)
        fc <- .cfunction("saige_get_sparse")
        # run
        packed.geno <- seqParallel(nfork, gds_grm, FUN=function(f) {
            seqApply(f, nm, fc, as.is="list", .useraw=TRUE, .list_dup=FALSE,
                .progress=nfork==1L && verbose)
        }, .balancing=TRUE, .bl_size=5000L, .bl_progress=verbose)
        rm(buffer)
    } else {
        # 2-bit packed genotypes
        packed.geno <- SeqArray:::.seqGet2bGeno(gds_grm, verbose)
    }
    if (verbose)
    {
        .cat("    using ", .pretty_size(as.double(object.size(packed.geno))),
            " (", ifelse(isTRUE(geno.sparse), "sparse", "dense"), " matrix)")
    }

    # initialize internal variables and buffers
    buf_std_geno <- double(4L*n_var)
    buf_sigma <- double(n_samp)
    buf_crossprod <- matrix(0.0, nrow=n_samp, ncol=num.thread)
    if (isTRUE(geno.sparse))
    {
        .Call(saige_store_sp_geno, packed.geno, n_samp, buf_std_geno, buf_sigma,
            buf_crossprod)
    } else {
        .Call(saige_store_2b_geno, packed.geno, n_samp, buf_std_geno, buf_sigma,
            buf_crossprod)
    }

    # parameters for fitting the model
    param <- list(
        num.thread = num.thread,
        seed = seed,
        tol = tol, tolPCG = tolPCG,
        maxiter = maxiter, maxiterPCG = maxiterPCG, no_iteration = FALSE,
        nrun = nrun,
        num.marker = 1L,
        traceCVcutoff = traceCVcutoff,
        ratioCVcutoff = ratioCVcutoff,
        verbose = verbose.detail,
        indent = "    "
    )
    tau.init[is.na(tau.init)] <- 0
    tau.init[tau.init < 0] <- 0

    # check gds_assoc: SNP variant IDs
    ii <- unique(c(snp_pair[,1L], snp_pair[,2L]))
    if (inherits(gds_assoc, "SeqVarGDSClass"))
    {
        i <- seqSetFilter(gds_assoc, variant.id=ii, ret.idx=TRUE, warn=FALSE,
            verbose=FALSE)
    } else {
        i <- match(ii, colnames(gds_assoc))
    }
    if (anyNA(i))
        stop("No variant ID(s): ", paste(ii[is.na(i)], collapse=", "))

    # check gds_assoc: sample IDs
    if (inherits(gds_assoc, "SeqVarGDSClass"))
    {
        i_geno <- seqSetFilter(gds_assoc, sample.id=sid, ret.idx=TRUE,
            verbose=FALSE)$sample_idx
    } else {
        i_geno <- match(sid, rownames(gds_assoc))
    }
    if (anyNA(i_geno))
    {
        if (all(is.na(i_geno)))
            stop("No common samples in the association GDS file.")
        if (verbose)
        {
            cat(sprintf("Missing sample rate in the association GDS file: %.2f%%\n",
                sum(is.na(i_geno))/length(i_geno)))
        }
    }

    # show the distribution of outcomes
    if (verbose)
    {
        y <- data[[phenovar]]
        .show_outcome(trait.type, y, phenovar)
    }


    ####  Initial guess of tau  ####

    ori_X <- X <- model.matrix(formula, data)
    y <- data[[phenovar]]
    if (use_approx_tau)
    {
        if (verbose)
        {
            cat("Fitting the model without the SNP markers to find the initial tau:\n")
            cat("    Formula: ")
            print(formula)
        }
        if (NCOL(X) > 1L && isTRUE(X.transform))
        {
            if (verbose)
                cat("    Transform on the design matrix with QR decomposition:\n")
            # check multi-collinearity
            m <- lm(y ~ X - 1)
            i_na <- which(is.na(m$coefficients))
            if (length(i_na) > 0L)
            {
                X <- X[, -i_na]
                if (verbose)
                {
                    .cat("        exclude ", length(i_na), " covariates (",
                        paste(colnames(X)[i_na], collapse=", "),
                        ") to avoid multi collinearity.")
                }
            }
            X_name <- colnames(X)
            Xqr <- qr(X)  # QR decomposition
            X_new <- qr.Q(Xqr) * sqrt(nrow(X))
            X_qrr <- qr.R(Xqr)
            data <- data.frame(cbind(y, X_new))
            nm <- paste0("x_", seq_len(ncol(X_new))-1L)
            colnames(data) <- c("y", nm)
            formula <- as.formula(paste("y ~", paste(nm, collapse=" + "), "-1"))
            if (verbose)
                .cat("        new formula: ", format(formula))
        }

        # fit the null model without SNP markers
        fit0 <- glm(formula, data=data, family=binomial)
        if (verbose)
            .show_coeff(fit0, "    ")
        obj.noK <- SPAtest:::ScoreTest_wSaddleApprox_NULL_Model(formula, data)
        # initial tau
        tau <- fixtau <- c(0,0)
        if (fit0$family$family %in% c("binomial", "poisson"))
            tau[1L] <- fixtau[1L] <- 1
        if (sum(tau.init[fixtau==0]) == 0)
            tau[fixtau==0] <- 0.5
        else
            tau[fixtau==0] <- tau.init[fixtau==0]
        fixtau <- tau
        # iterate
        X <- model.matrix(fit0)
        glmm <- .Call(saige_fit_AI_PCG_binary, fit0, X, tau, param)
        tau.init <- glmm$tau
        if (verbose) cat("    ")
    } else {
        # no initial guess
        tau <- c(0,0)
        if (trait.type == "binary")
            tau[1L] <- 1
        if (sum(tau.init[tau==0]) == 0)
            tau[tau==0] <- 0.5
        else
            tau[tau==0] <- tau.init[tau==0]
        tau.init <- tau
    }

    if (verbose && use_approx_tau)
        .cat("Use tau for the interaction: (", tau.init[1L], ", ", tau.init[2L], ")")

    # check whether use glm p-value threshold
    if (is.na(glm_threshold)) glm_threshold <- FALSE
    if (isTRUE(glm_threshold)) glm_threshold <- 0.01
    if (verbose && !isFALSE(glm_threshold))
        .cat("GLM p-value threshold: ", glm_threshold)


    ####  Enumerate each SNP pair  ####

    param$no_iteration <- use_approx_tau
    rv_dat <- NULL  # output data
    if (verbose)
        .cat("Testing the interaction, # of SNP pairs: ", nrow(snp_pair))

    for (ii in seq_len(nrow(snp_pair)))
    {
        i1 <- snp_pair[ii, 1L]
        i2 <- snp_pair[ii, 2L]
        if (verbose)
        {
            s <- paste0("==> ", ii, ": SNP ", i1, " x SNP ", i2,
                "\t[", date(), "] <==")
            .cat(.crayon_inverse(s))
        }

        # first SNP
        if (inherits(gds_assoc, "SeqVarGDSClass"))
        {
            seqSetFilter(gds_assoc, variant.id=i1, verbose=FALSE)
            g1 <- .minor_allele_geno(seqGetData(gds_assoc, "$dosage_alt")[i_geno])
            s1 <- seqGetData(gds_assoc, "$chrom_pos_allele")
        } else {
            i <- match(i1, colnames(gds_assoc))
            g1 <- .minor_allele_geno(gds_assoc[i_geno, i])
            s1 <- i1
            i1 <- i
        }
        maf1 <- mean(g1)*0.5
        if (verbose)
            cat(sprintf("    SNP1 (%s), MAF: %.5g\n", s1, maf1))

        # second SNP
        if (inherits(gds_assoc, "SeqVarGDSClass"))
        {
            seqSetFilter(gds_assoc, variant.id=i2, verbose=FALSE)
            g2 <- .minor_allele_geno(seqGetData(gds_assoc, "$dosage_alt")[i_geno])
            s2 <- seqGetData(gds_assoc, "$chrom_pos_allele")
        } else {
            i <- match(i2, colnames(gds_assoc))
            g2 <- .minor_allele_geno(gds_assoc[i_geno, i])
            s2 <- i2
            i2 <- i
        }
        maf2 <- mean(g2)*0.5
        if (verbose)
            cat(sprintf("    SNP2 (%s), MAF: %.5g\n", s2, maf2))

        # calculate
        X <- cbind(ori_X, g1, g2)
        # check multi-collinearity
        if (verbose.detail)
            cat("    Transform on the design matrix with QR decomposition\n")
        m <- lm(y ~ X - 1)
        i_na <- which(is.na(m$coefficients))
        if (length(i_na) > 0L)
        {
            X <- X[, -i_na]
            if (verbose.detail)
            {
                .cat("    excluding ", length(i_na), " covariates (",
                    paste(colnames(X)[i_na], collapse=", "),
                    ") to avoid multi collinearity.")
            }
        }
        X_name <- colnames(X)
        Xqr <- qr(X)  # QR decomposition
        X_new <- qr.Q(Xqr) * sqrt(nrow(X))
        X_qrr <- qr.R(Xqr)
        new_dat <- data.frame(cbind(y, X_new))
        nm <- paste0("x_", seq_len(ncol(X_new))-1L)
        colnames(new_dat) <- c("y", nm)
        fm <- as.formula(paste("y ~", paste(nm, collapse=" + "), "-1"))
        if (verbose.detail)
            .cat("    New formula: ", format(fm))

        # fit the model
        if (trait.type == "binary")
        {
            fit0 <- glm(fm, data=new_dat, family=binomial)
            if (verbose.detail)
                .show_coeff(fit0, "    ", !param$no_iteration)
            obj.noK <- SPAtest:::ScoreTest_wSaddleApprox_NULL_Model(fm, new_dat)
            if (verbose.detail)
                .cat("    Calculate the interaction term:")
            X <- model.matrix(fit0)
            pv <- NA_real_
            run_glmm <- TRUE
            if (!isFALSE(glm_threshold))
            {
                # calculate the interaction term using glm
                param$no_iteration <- TRUE
                tau <- c(1, 0)
                glmm <- .Call(saige_fit_AI_PCG_binary, fit0, X, tau, param)
                d <- .Call(saige_GxG_snp_bin, fit0, glmm, g1*g2, obj.noK,
                    param, verbose)
                param$no_iteration <- use_approx_tau
                pv <- d$pval
                pv2 <- d$p.norm
                d$pval <- d$p.norm <- NaN
                d$p.glm <- pv
                d$p.glm.norm <- pv2
                run_glmm <- is.finite(pv) && (pv <= glm_threshold)
                if (verbose.detail)
                {
                    .cat("    glm p-value: ", pv, " ",
                        ifelse(run_glmm, "<= threshold", "> threshold (skip glmm)"))
                }
            }
            if (run_glmm)
            {
                # iterate
                tau <- tau.init
                glmm <- .Call(saige_fit_AI_PCG_binary, fit0, X, tau, param)
                d <- .Call(saige_GxG_snp_bin, fit0, glmm, g1*g2, obj.noK, param,
                    verbose)
                if (!is.na(pv))
                {
                    d$p.glm <- pv
                    d$p.glm.norm <- pv2
                }
            }

        } else if (trait.type == "quantitative")    
        {
            stop("Not implement yet.")
            # inverse normal transformation
            if (isTRUE(inv.norm))
            {
                if (isTRUE(X.transform)) phenovar <- "y"
                fit0 <- glm(formula, data=data)
                resid.sd <- sd(fit0$residuals)
                new.y <- .rank_norm(fit0$residuals) * resid.sd
                data[[phenovar]] <- new.y
                if (verbose.detail)
                {
                    .cat("Inverse normal transformation on residuals with standard deviation: ",
                        resid.sd)
                }
            }

            # fit the null model
            fit0 <- glm(formula, data=data)
            if (verbose.detail) .show_coeff(fit0, "    ")

            # ScoreTest_wSaddleApprox_NULL_Model_q
            obj.noK <- list()
            X1 <- model.matrix(fit0)
            X1 <- SPAtest:::ScoreTest_wSaddleApprox_Get_X1(X1)
            obj.noK$y <- fit0$y
            obj.noK$mu <- fit0$fitted.values
            obj.noK$res <- fit0$y - obj.noK$mu
            obj.noK$V <- rep(1, length(fit0$y))
            obj.noK$X1 <- X1
            obj.noK$XV <- t(X1)
            obj.noK$XXVX_inv <- X1 %*% solve(t(X1) %*% X1)
            class(obj.noK) <- "SA_NULL"

            # initial tau
            tau <- tau.init
            if (sum(tau) == 0) tau <- c(0.5, 0.5)
            yy <- fit0$y
            offset <- fit0$offset
            if (is.null(offset)) offset <- rep(0, length(yy))
            eta <- fit0$linear.predictors
            mu <- fit0$fitted.values
            mu.eta <- fit0$family$mu.eta(eta)
            Y <- eta - offset + (yy - mu)/mu.eta
            tau <- var(Y) * tau / sum(tau)

            # iterate
            glmm <- .Call(saige_fit_AI_PCG_quant, fit0, X1, tau, param)

            # calculate the variance ratio
            if (verbose)
                .cat("Calculate the interaction term:")
            set.seed(seed)
            var.ratio <- .Call(saige_calc_var_ratio_quant, fit0, glmm, obj.noK,
                param, sample.int(n_var, n_var))
            var.ratio <- var.ratio[order(var.ratio$id), ]
            var.ratio$id <- seqGetData(gds_grm, "variant.id")[var.ratio$id]
            rownames(var.ratio) <- NULL

            # ans$obj.glm.null <- fit0
            glmm$obj.noK <- obj.noK
            glmm$var.ratio <- var.ratio

        } else {
            stop("Invalid 'trait.type'.")
        }

        # combine results
        d <- cbind(id1=i1, snp1=s1, maf1=maf1, id2=i2, snp2=s2, maf2=maf2, d)
        rv_dat <- rbind(rv_dat, d)
        rv_ans <- rv_dat
        if (ncol(snp_pair) > 2L)
        {
            s <- c(colnames(rv_ans), colnames(snp_pair)[3:ncol(snp_pair)])
            rv_ans <- cbind(rv_ans, snp_pair[1:nrow(rv_ans), 3:ncol(snp_pair)])
            colnames(rv_ans) <- s
        }
        if (use_approx_tau)
            attr(rv_ans, "tau_G") <- tau.init[2L]
        # save to a file
        if (!is.na(model.savefn) && model.savefn!="")
        {
            if (verbose)
                .cat("    Save the results to ", sQuote(model.savefn))
            if (grepl("\\.(rda|RData)$", model.savefn, ignore.case=TRUE))
            {
                .x <- rv_ans
                save(.x, file=model.savefn)
            } else if (grepl("\\.rds$", model.savefn, ignore.case=TRUE))
            {
                saveRDS(rv_ans, file=model.savefn)
            } else if (grepl("\\.txt$", model.savefn, ignore.case=TRUE))
            {
                write.table(rv_ans, model.savefn, sep="\t", quote=FALSE, row.names=FALSE)
            } else if (grepl("\\.csv$", model.savefn, ignore.case=TRUE))
            {
                write.csv(rv_ans, model.savefn, row.names=FALSE)
            } else
                stop("Unknown format of the output file, and it should be RData or RDS.")
        }
    }

    if (verbose)
        .cat(.crayon_inverse(paste0("Done (", date(), ").")))

    if (!is.na(model.savefn) && model.savefn!="")
        return(invisible(rv_ans))
    else
        return(rv_ans)
}
