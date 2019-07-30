#######################################################################
#
# Package name: SAIGEgds
#
# Description:
#     Scalable and accurate implementation of generalized mixed models
# using GDS files
#
# Copyright (C) 2019    Xiuwen Zheng / AbbVie-ComputationalGenomics
# License: GPL-3
#


# Package-wide variable
.packageEnv <- new.env()


#######################################################################
# Internal functions
#

.cfunction <- function(name)
{
    fn <- function(x) NULL
    f <- quote(.Call(SEQ_ExternalName1, x))
    f[[1L]] <- .Call
    f[[2L]] <- getNativeSymbolInfo(name, "SAIGEgds")$address
    body(fn) <- f
    fn
}

.cfunction2 <- function(name)
{
    fn <- function(x, y) NULL
    f <- quote(.Call(SEQ_ExternalName2, x, y))
    f[[1L]] <- .Call
    f[[2L]] <- getNativeSymbolInfo(name, "SAIGEgds")$address
    body(fn) <- f
    fn
}

.pretty <- function(x) prettyNum(x, big.mark=",", scientific=FALSE)

SIMD <- function() .Call(saige_simd_version)

.crayon_inverse <- function(s)
{
    if (requireNamespace("crayon", quietly=TRUE))
        s <- crayon::inverse(s)
    s
}

.crayon_underline <- function(s)
{
    if (requireNamespace("crayon", quietly=TRUE))
        s <- crayon::underline(s)
    s
}


# Internal model checking
.check_saige_model <- function(obj)
{
    stopifnot(is.list(obj))
    for (nm in c("sample.id", "trait.type", "var.ratio"))
    {
        if (!(nm %in% names(obj)))
            stop("'", nm, "' should be stored in the SAIGE model.")
    }
    if (!(obj$trait.type %in% c("binary", "quantitative")))
        stop("'trait.type' should be binary or quantitative.")
    invisible()
}


# Write to GDS file
.write_gds <- function(out.gds, out.nm, in.gds, in.nm, cm)
{
    n <- add.gdsn(out.gds, out.nm, storage=index.gdsn(in.gds, in.nm), compress=cm)
    seqApply(in.gds, in.nm, `c`, as.is=n)
    readmode.gdsn(n)
    invisible()
}



#######################################################################
# Load the association p-values in a GDS file
#

seqSAIGE_LoadPval <- function(fn, varnm=NULL, index=NULL, verbose=TRUE)
{
    # check
    stopifnot(is.character(fn), length(fn)>0L, all(!is.na(fn)))
    stopifnot(is.null(varnm) || is.character(varnm))
    stopifnot(is.null(index) || is.numeric(index) || is.logical(index))
    stopifnot(is.logical(verbose), length(verbose)==1L)

    if (length(fn) == 1L)
    {
        if (verbose)
            cat("Loading '", fn, "' ...\n", sep="")
        if (grepl("\\.gds$", fn, ignore.case=TRUE))
        {
            f <- openfn.gds(fn)
            on.exit(closefn.gds(f))
            if (identical(get.attr.gdsn(f$root)$FileFormat, "SAIGE_OUTPUT"))
            {
                if (is.null(varnm))
                    varnm <- ls.gdsn(f$root)
                varnm <- setdiff(varnm, "sample.id")
                rv <- list()
                for (nm in varnm)
                    rv[[nm]] <- readex.gdsn(index.gdsn(f, nm), index)
                rv <- as.data.frame(rv, stringsAsFactors=FALSE)
            } else {
                stop("FileFormat should be 'SAIGE_OUTPUT'.")
            }
        } else if (grepl("\\.(rda|RData)$", fn, ignore.case=TRUE))
        {
            rv <- get(load(fn))
            if (!is.null(varnm)) rv <- rv[, varnm]
            if (!is.null(index)) rv <- rv[index, ]
        } else {
            stop(sprintf("Unknown format (%s), should be RData or gds.",
                basename(fn)))
        }
    } else {
        if (!is.null(index))
            stop("'index' should be NULL for multiple input files.")
        rv <- sapply(fn, function(nm) seqSAIGE_LoadPval(nm, varnm), simplify=FALSE)
        if (verbose) cat("Merging ...")
        rv <- Reduce(rbind, rv)
        if (verbose) cat(" [done]\n")
    }
    rv
}



#######################################################################
# Fit the null model
#

seqFitNullGLMM_SPA <- function(formula, data, gdsfile,
    trait.type=c("binary", "quantitative"), sample.col="sample.id", maf=0.005,
    missing.rate=0.01, max.num.snp=1000000L, variant.id=NULL, inv.norm=TRUE,
    X.transform=TRUE, tol=0.02, maxiter=20L, nrun=30L, tolPCG=1e-5, maxiterPCG=500L,
    num.marker=30L, tau.init=c(0,0), traceCVcutoff=0.0025, ratioCVcutoff=0.001,
    geno.sparse=TRUE, num.thread=1L, model.savefn="", seed=200L, no.fork.loading=FALSE,
    verbose=TRUE)
{
    stopifnot(inherits(formula, "formula"))
    stopifnot(is.data.frame(data))
    stopifnot(is.character(gdsfile) | inherits(gdsfile, "SeqVarGDSClass"))
    trait.type <- match.arg(trait.type)
    stopifnot(is.character(sample.col), length(sample.col)==1L, !is.na(sample.col))
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
    stopifnot(is.numeric(num.marker), length(num.marker)==1L)
    stopifnot(is.numeric(tau.init), length(tau.init)==2L)
    stopifnot(is.numeric(traceCVcutoff), length(traceCVcutoff)==1L)
    stopifnot(is.numeric(ratioCVcutoff), length(ratioCVcutoff)==1L)
    stopifnot(is.logical(geno.sparse), length(geno.sparse)==1L)
    stopifnot(is.numeric(num.thread), length(num.thread)==1L)
    stopifnot(is.character(model.savefn), length(model.savefn)==1L)
    stopifnot(is.numeric(seed), length(seed)==1L, is.finite(seed))
    stopifnot(is.logical(no.fork.loading), length(no.fork.loading)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    if (verbose)
        cat(.crayon_inverse("SAIGE association analysis:\n"))

    if (is.character(gdsfile))
    {
        if (verbose)
            cat("Open the genotype file '", gdsfile, "'\n", sep="")
        gdsfile <- seqOpen(gdsfile)
        on.exit(seqClose(gdsfile))
    } else {
        # save the filter on GDS file
        seqSetFilter(gdsfile, action="push", verbose=FALSE)
        on.exit(seqSetFilter(gdsfile, action="pop", verbose=FALSE))
    }

    # variables in the formula
    vars <- all.vars(formula)
    phenovar <- all.vars(formula)[1L]
    y <- data[[phenovar]]
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
    seqResetFilter(gdsfile, sample=TRUE, verbose=FALSE)
    sid <- seqGetData(gdsfile, "sample.id")
    i <- match(sid, data[[sample.col]])
    i <- i[!is.na(i)]
    data <- data[i, ]
    if (nrow(data) <= 0L)
        stop("No common sample.id between 'data' and the GDS file.")
    seqSetFilter(gdsfile, sample.id=data[[sample.col]], verbose=FALSE)

    if (is.null(variant.id))
    {
        # filters of maf, mac, missing.rate
        if (verbose)
            cat("Filtering variants:\n")
        seqSetFilterCond(gdsfile, maf=maf, mac=NaN, missing.rate=missing.rate,
            parallel=num.thread, .progress=TRUE, verbose=FALSE)
    } else {
        seqSetFilter(gdsfile, variant.id=variant.id, verbose=FALSE)
    }

    # the max mumber of SNPs
    v <- seqGetFilter(gdsfile)$variant.sel
    n <- sum(v, na.rm=TRUE)
    if (max.num.snp>0L && n>max.num.snp)
    {
        set.seed(seed)
        seqSetFilter(gdsfile, variant.sel=sample(which(v), max.num.snp),
            verbose=FALSE)
    }

    # get the number of samples / variants
    dm <- seqSummary(gdsfile, "genotype", verbose=FALSE)$seldim
    n_samp <- dm[2L]
    n_var  <- dm[3L]
    if (verbose)
    {
        cat("Fit the null model:", format(formula), "+ var(GRM)\n")
        cat("    # of samples: ", .pretty(n_samp), "\n", sep="")
        cat("    # of variants:", .pretty(n_var))
        if (n > max.num.snp)
            cat(" (randomly selected from ", .pretty(n), ")", sep="")
        cat("\n")
    }

    # set the number of internal threads
    if (is.na(num.thread) || num.thread < 1L)
        num.thread <- 1L
    setThreadOptions(num.thread)
    if (verbose)
    {
        cat("    using ", num.thread, " thread",
            ifelse(num.thread>1L, "s", ""), "\n", sep="")
    }

    if (isTRUE(X.transform))
    {
        if (verbose)
            cat("Transform on the design matrix with QR decomposition:\n")
        X <- model.matrix(formula, data)
        frm <- model.frame(formula, data)
        y <- model.response(frm, type="any")
        # check multi-collinearity
        m <- lm(y ~ X - 1)
        i_na <- which(is.na(m$coefficients))
        if (length(i_na) > 0L)
        {
            X <- X[, -i_na]
            if (verbose)
            {
                cat("    exclude ", length(i_na), " covariates (",
                    paste(colnames(X)[i_na], collapse=", "),
                    ") to avoid multi collinearity.\n", sep="")
            }
        }
        X_name <- colnames(X)
        Xqr = qr(X)  # QR decomposition
        X_new <- qr.Q(Xqr) * sqrt(nrow(X))
        X_qrr <- qr.R(Xqr)
        data <- data.frame(cbind(y, X_new))
        nm <- paste0("x", seq_len(ncol(X_new))-1L)
        colnames(data) <- c("y", nm)
        formula <- as.formula(paste("y ~", paste(nm, collapse=" + "), "-1"))
        if (verbose)
            cat("    new formula: ", format(formula), "\n", sep="")
    }

    # load SNP genotypes
    if (verbose)
        cat("Start loading SNP genotypes:\n")
    nfork <- 1L
    if (SeqArray:::.IsForking(num.thread) && !no.fork.loading)
        nfork <- num.thread
    if (isTRUE(geno.sparse))
    {
        # sparse genotypes
        buffer <- integer(n_samp + 4L)
        fc <- .cfunction2("saige_get_sparse")
        packed.geno <- seqParallel(nfork, gdsfile, FUN=function(f) {
            seqApply(f, "$dosage_alt", fc, as.is="list", y=buffer, .useraw=TRUE,
                .list_dup=FALSE, .progress=nfork==1L && verbose)
        }, .balancing=TRUE, .bl_size=5000L, .bl_progress=verbose)
        rm(buffer)
    } else {
        # 2-bit packed genotypes
        packed.geno <- SeqArray:::.seqGet2bGeno(gdsfile, verbose)
    }
    if (verbose)
    {
        cat("    using ")
        cat(SeqArray:::.pretty_size(as.double(object.size(packed.geno))))
        cat(ifelse(isTRUE(geno.sparse), " (sparse matrix)\n", " (dense matrix)\n"))
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
        maxiter = maxiter, maxiterPCG = maxiterPCG,
        nrun = nrun,
        num.marker = num.marker,
        traceCVcutoff = traceCVcutoff,
        ratioCVcutoff = ratioCVcutoff,
        verbose = verbose
    )

    tau.init[is.na(tau.init)] <- 0
    tau.init[tau.init < 0] <- 0

    # fit the model
    if (trait.type == "binary")
    {
        # binary outcome
        if (verbose)
        {
            cat("Binary outcome: ", phenovar, "\n", sep="")
            if (isTRUE(X.transform))
                y <- data$y
            else
                y <- data[[phenovar]]
            v <- table(y)
            v <- data.frame(v, as.numeric(prop.table(v)))
            v[, 1L] <- paste0("      ", v[, 1L])
            colnames(v) <- c(phenovar, "Number", "Proportion")
            print(v, row.names=FALSE)
        }

        # fit the null model
        fit0 <- glm(formula, data=data, family=binomial)
        if (verbose)
        {
            cat("Initial fixed-effect coefficients:\n")
            v <- as.data.frame(t(fit0$coefficients))
            rownames(v) <- "   "
            print(v)
        }
        obj.noK <- SPAtest:::ScoreTest_wSaddleApprox_NULL_Model(formula, data)

        # initial tau
        tau <- fixtau <- c(0,0)
        if (fit0$family$family %in% c("binomial", "poisson"))
            tau[1] <- fixtau[1] <- 1
        if (sum(tau.init[fixtau==0]) == 0)
            tau[fixtau==0] <- 0.5
        else
            tau[fixtau==0] <- tau.init[fixtau==0]

        # iterate
        X <- model.matrix(fit0)
        glmm <- .Call(saige_fit_AI_PCG_binary, fit0, X, tau, param)

        # calculate the variance ratio
        if (verbose)
            cat(.crayon_underline("Calculate the average ratio of variances:\n"))
        set.seed(seed)
        var.ratio <- .Call(saige_calc_var_ratio_binary, fit0, glmm, obj.noK,
            param, sample.int(n_var, n_var))
        var.ratio <- var.ratio[order(var.ratio$id), ]
        var.ratio$id <- seqGetData(gdsfile, "variant.id")[var.ratio$id]
        rownames(var.ratio) <- NULL

        # ans$obj.glm.null <- fit0
        glmm$obj.noK <- obj.noK
        glmm$var.ratio <- var.ratio

    } else if (trait.type == "quantitative")    
    {
        # quantitative outcome
        if (verbose)
        {
            cat("Quantitative outcome: ", phenovar, "\n", sep="")
            if (isTRUE(X.transform))
                y <- data$y
            else
                y <- data[[phenovar]]
            v <- data.frame(mean=mean(y), sd=sd(y), min=min(y), max=max(y))
            rownames(v) <- "   "
            print(v)
        }

        stop("Quantitative implementation is not ready!")

        # fit the null model
        fit0 <- glm(formula, data=data)
        if (verbose)
        {
            cat("Initial fixed-effect coefficients:\n")
            v <- as.data.frame(t(fit0$coefficients))
            rownames(v) <- "   "
            print(v)
        }
        # ScoreTest_wSaddleApprox_NULL_Model_q
        obj.noK <- list()
        X1 <- model.matrix(fit0)
        X1 <- SPAtest:::ScoreTest_wSaddleApprox_Get_X1(X1)
        obj.noK$y <- fit0$y
        obj.noK$mu <- fit0$fitted.values
        obj.noK$res <- fit0$y - obj.noK$mu
        obj.noK$V <- 1
        obj.noK$X1 <- X1
        obj.noK$XV <- t(X1)
        obj.noK$XVX_inv <- solve(t(X1) %*% X1)
        obj.noK$XXVX_inv <- X1 %*% obj.noK$XVX_inv
        class(obj.noK) <- "SA_NULL"

        # initial tau
        tau <- tau.init
        if (sum(tau) == 0)
        {
            y <- fit0$y
            offset <- fit0$offset
            if (is.null(offset)) offset <- rep(0, length(y))
            eta <- fit0$linear.predictors
            mu <- fit0$fitted.values
            mu.eta <- fit0$family$mu.eta(eta)
            Y <- eta - offset + (y - mu)/mu.eta
            tau[] <- var(Y)/2
        }

        # iterate
        glmm <- .Call(saige_fit_AI_PCG_quant, fit0, X1, tau, param)

        # calculate the variance ratio
        if (verbose)
            cat(.crayon_underline("Calculate the average ratio of variances:\n"))
        set.seed(seed)
        var.ratio <- .Call(saige_calc_var_ratio_quant, fit0, glmm, obj.noK,
            param, sample.int(n_var, n_var))
        var.ratio <- var.ratio[order(var.ratio$id), ]
        var.ratio$id <- seqGetData(gdsfile, "variant.id")[var.ratio$id]
        rownames(var.ratio) <- NULL

        # ans$obj.glm.null <- fit0
        glmm$obj.noK <- obj.noK
        glmm$var.ratio <- var.ratio

    } else {
        stop("Invalid 'trait.type'.")    
    } 

    if (verbose)
    {
        cat("    ratio avg. is ", mean(var.ratio$ratio),
            ", sd: ", sd(var.ratio$ratio), "\n", sep="")
    }

    # tweak the result
    if (!isTRUE(X.transform))
    {
        names(glmm$coefficients) <- colnames(obj.noK$X1)
    } else {
        coef <- solve(X_qrr, glmm$coefficients * sqrt(nrow(data)))
        names(coef) <- X_name
        glmm$coefficients <- coef
    }
    names(glmm$tau) <- c("Sigma_E", "Sigma_G")
    glmm$trait.type <- trait.type
    glmm$sample.id <- seqGetData(gdsfile, "sample.id")
    glmm$variant.id <- seqGetData(gdsfile, "variant.id")

    if (!is.na(model.savefn) && model.savefn!="")
    {
        cat("Save the model to '", model.savefn, "'\n", sep="")
        .glmm <- glmm
        save(.glmm, file=model.savefn)
    }
    if (verbose)
        cat(.crayon_inverse("Done."), "\n", sep="")

    if (!is.na(model.savefn) && model.savefn!="")
        return(invisible(glmm))
    else
        return(glmm)
}




#######################################################################
# SAIGE single variant analysis
#

seqAssocGLMM_SPA <- function(gdsfile, modobj, maf=NaN, mac=10,
    dsnode="", spa.pval=0.05, var.ratio=NaN,
    res.savefn="", res.compress=c("LZMA", "LZMA_RA", "ZIP", "ZIP_RA", "none"),
    parallel=FALSE, verbose=TRUE)
{
    stopifnot(inherits(gdsfile, "SeqVarGDSClass") | is.character(gdsfile))
    stopifnot(is.numeric(maf), length(maf)==1L)
    stopifnot(is.numeric(mac), length(mac)==1L)
    stopifnot(is.character(dsnode), length(dsnode)==1L, !is.na(dsnode))
    stopifnot(is.character(dsnode), length(dsnode)==1L, !is.na(dsnode))
    stopifnot(is.numeric(spa.pval), length(spa.pval)==1L)
    stopifnot(is.numeric(var.ratio), length(var.ratio)==1L)
    stopifnot(is.character(res.savefn), length(res.savefn)==1L)
    res.compress <- match.arg(res.compress)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    # check model
    if (is.character(modobj))
    {
        stopifnot(length(modobj)==1L)
        modobj <- get(load(modobj))
    }
    .check_saige_model(modobj)

    # GDS file
    if (is.character(gdsfile))
    {
        if (verbose)
            cat("Open '", gdsfile, "' ...\n", sep="")
        gdsfile <- seqOpen(gdsfile)
        on.exit(seqClose(gdsfile))
    } else {
        # save the filter on GDS file
        seqSetFilter(gdsfile, action="push", verbose=FALSE)
        on.exit(seqSetFilter(gdsfile, action="pop", verbose=FALSE))
    }

    # determine the GDS node for dosages
    if (dsnode == "")
    {
        n <- index.gdsn(gdsfile, "genotype/data", silent=TRUE)
        if (!is.null(n))
        {
            dsnode <- "$dosage_alt"
        } else {
            n <- index.gdsn(gdsfile, "annotation/format/DS", silent=TRUE)
            if (!is.null(n))
            {
                dsnode <- "annotation/format/DS"
            } else {
                stop("Dosages should be stored in genotype or annotation/format/DS.")
            }
        }
    }

    if (verbose)
        cat("SAIGE association analysis:\n")

    # check sample ID
    seqSetFilter(gdsfile, sample.id=modobj$sample.id, verbose=FALSE)
    sid <- seqGetData(gdsfile, "sample.id")
    if (length(sid) != length(modobj$sample.id))
        stop("Some of sample IDs are not available in the GDS file.")
    ii <- match(sid, modobj$sample.id)
    if (any(is.na(ii)))
        stop("Sample IDs do not match.")

    dm <- seqSummary(gdsfile, "genotype", verbose=FALSE)$seldim
    nVariant <- dm[3L]
    if (verbose)
    {
        cat("    # of samples: ", .pretty(dm[2L]), "\n", sep="")
        cat("    # of variants: ", .pretty(dm[3L]), "\n", sep="")
        cat("    p-value threshold for SPA adjustment: ", spa.pval, "\n", sep="")
    }

    if (!is.finite(var.ratio))
        var.ratio <- mean(modobj$var.ratio$ratio, na.rm=TRUE)
    if (verbose)
        cat("    variance ratio for approximation: ", var.ratio, "\n", sep="")

    if (dm[2L] <= 0) stop("No sample in the genotypic data set!")
    if (dm[3L] <= 0) stop("No variant in the genotypic data set!")

    # initialize the internal model parameters
    y <- unname(modobj$obj.noK$y)
    mu <- unname(modobj$fitted.values)
    X1 <- modobj$obj.noK$X1[ii, ]
    n <- length(ii)
    mobj <- list(
        maf = maf, mac = mac, spa.pval = spa.pval,
        tau = modobj$tau,
        y = y[ii], mu = mu[ii],
        y_mu = y[ii] - mu[ii],  # y - mu
        mu2 = (mu * (1 - mu))[ii],
        t_XXVX_inv = t(modobj$obj.noK$XXVX_inv[ii, ]),  # K x n_samp (K << n_samp, more efficient)
        XV = modobj$obj.noK$XV[, ii],  # K x n_samp
        t_XVX_inv_XV = t(modobj$obj.noK$XXVX_inv[ii, ] * modobj$obj.noK$V[ii]),  # K x n_samp
        t_X = t(X1),  # K x n_samp
        var.ratio = var.ratio,
        # buffer
        buf_dosage = double(n),
        buf_coeff = double(nrow(modobj$obj.noK$XV)),
        buf_adj_g = double(n),
        buf_index = integer(n),
        buf_B = double(n),
        buf_g_tilde = double(n),
        buf_spa = double(n+n),
        buf_tmp = double(ncol(X1))
    )
    mobj$XVX <- t(X1) %*% (X1 * mobj$mu2)  # a matrix: K x K
    mobj$S_a <- colSums(X1 * mobj$y_mu)    # a vector of size K

    if (!is.finite(mobj$var.ratio))
        stop("Invalid variance ratio in the SAIGE model.")

    # is forking or not?
    is_fork <- SeqArray:::.IsForking(parallel)
    njobs <- SeqArray:::.NumParallel(parallel)
    if (verbose && njobs>1L)
        cat("    # of processes: ", njobs, "\n", sep="")

    # initialize internally
    if (njobs<=1L || is_fork)
    {
        .Call(saige_score_test_init, mobj)
        initfun <- finalfun <- NULL
    } else {
        if (verbose)
            cat("Distribute the model parameters to the", njobs, "processes\n")
        # pass the model parameters to each process
        initfun <- function(proc_id, mobj)
        {
            eval(parse(text="library(Rcpp)"))
            eval(parse(text="library(SAIGEgds)"))
            .packageEnv$modobj <- mobj
            .Call(saige_score_test_init, mobj)
        }
        # clear when exit
        finalfun <- function(proc_id, param)
        {
            .packageEnv$modobj <- NULL
            remove(modobj, envir=.packageEnv)
        }
    }

    # scan all (selected) variants
    if (modobj$trait.type == "binary")
    {
        rv <- seqParallel(parallel, gdsfile, split="by.variant",
            .initialize=initfun, .finalize=finalfun, .initparam=mobj,
            .balancing=TRUE, .bl_size=50000L, .bl_progress=verbose,
            FUN = function(f, dsnode, pverbose)
            {
                seqApply(f, dsnode, .cfunction("saige_score_test_bin"), as.is="list",
                    parallel=FALSE, .progress=pverbose, .list_dup=FALSE, .useraw=NA)
            }, dsnode=dsnode, pverbose=verbose & (njobs==1L))
    } else if (modobj$trait.type == "quantitative")    
    {
        stop("Quantitative implementation is not ready.")

        rv <- seqParallel(parallel, gdsfile, split="by.variant",
            FUN = function(f, dsnode, verbose)
            {
                seqApply(f, dsnode, .cfunction("saige_score_test_quant"), as.is="list",
                    parallel=FALSE, .progress=verbose, .list_dup=FALSE, .useraw=NA)
            }, dsnode=dsnode, verbose=verbose)
    } else {
        stop("Invalid 'modobj$trait.type'.")
    }

    # if any maf/mac filter
    if (length(rv) != nVariant)
        stop("Internal error: seqParallel() returns a vector of wrong length.")
    x <- sapply(rv, is.null)
    if (any(x))
    {
        x <- !x
        seqSetFilter(gdsfile, variant.sel=x, action="intersect", verbose=FALSE)
        rv <- rv[x]
    }
    if (verbose)
    {
        cat("# of variants after filtering MAF/MAC: ",
            .pretty(length(rv)), "\n", sep="")
    }

    # output to a GDS file?
    isfn <- !is.na(res.savefn) && res.savefn!=""
    if (isfn && grepl("\\.gds$", res.savefn, ignore.case=TRUE))
    {
        if (verbose)
            cat("Save to '", res.savefn, "' ...\n", sep="")
        cm <- res.compress[1L]
        # create a GDS file
        f <- createfn.gds(res.savefn)
        on.exit(closefn.gds(f), add=TRUE)
        put.attr.gdsn(f$root, "FileFormat", "SAIGE_OUTPUT")
        put.attr.gdsn(f$root, "Version",
            paste0("SAIGEgds_", packageVersion("SAIGEgds")))
        # add sample IDs
        add.gdsn(f, "sample.id", seqGetData(gdsfile, "sample.id"),
            compress=cm, closezip=TRUE)
        # write data
        .write_gds(f, "id", gdsfile, "variant.id", cm)
        .write_gds(f, "chr", gdsfile, "chromosome", cm)
        .write_gds(f, "pos", gdsfile, "position", cm)
        # rs.id
        if (!is.null(index.gdsn(gdsfile, "annotation/id", silent=TRUE)))
        {
            .write_gds(f, "rs.id", gdsfile, "annotation/id", cm)
        }
        # ref and alt alleles
        add.gdsn(f, "ref", seqGetData(gdsfile, "$ref"), compress=cm, closezip=TRUE)
        add.gdsn(f, "alt", seqGetData(gdsfile, "$alt"), compress=cm, closezip=TRUE)
        # other data
        add.gdsn(f, "AF.alt", sapply(rv, `[`, i=1L), compress=cm, closezip=TRUE)
        add.gdsn(f, "AC.alt", sapply(rv, `[`, i=2L), compress=cm, closezip=TRUE)
        add.gdsn(f, "num", as.integer(sapply(rv, `[`, i=3L)), compress=cm, closezip=TRUE)
        add.gdsn(f, "beta", sapply(rv, `[`, i=4L), compress=cm, closezip=TRUE)
        add.gdsn(f, "SE", sapply(rv, `[`, i=5L), compress=cm, closezip=TRUE)
        add.gdsn(f, "pval", sapply(rv, `[`, i=6L), compress=cm, closezip=TRUE)
        if (modobj$trait.type == "binary")
        {
            add.gdsn(f, "pval.noadj", sapply(rv, `[`, i=7L), compress=cm, closezip=TRUE)
            add.gdsn(f, "converged", sapply(rv, `[`, i=8L)==1, compress=cm, closezip=TRUE)
        }
        if (verbose) cat("Done.\n")
        # output
        invisible()
    } else {
        # output
        ans <- data.frame(
            id  = seqGetData(gdsfile, "variant.id"),
            chr = seqGetData(gdsfile, "chromosome"),
            pos = seqGetData(gdsfile, "position"),
            stringsAsFactors = FALSE
        )
        # add RS IDs if possible
        if (!is.null(index.gdsn(gdsfile, "annotation/id", silent=TRUE)))
            ans$rs.id <- seqGetData(gdsfile, "annotation/id")
        ans$ref <- seqGetData(gdsfile, "$ref")
        ans$alt <- seqGetData(gdsfile, "$alt")
        ans$AF.alt <- sapply(rv, `[`, i=1L)
        ans$AC.alt <- sapply(rv, `[`, i=2L)
        ans$num  <- as.integer(sapply(rv, `[`, i=3L))
        ans$beta <- sapply(rv, `[`, i=4L)
        ans$SE   <- sapply(rv, `[`, i=5L)
        ans$pval <- sapply(rv, `[`, i=6L)
        if (modobj$trait.type == "binary")
        {
            ans$pval.noadj <- sapply(rv, `[`, i=7L)
            ans$converged <- as.logical(sapply(rv, `[`, i=8L))
        }

        # save file?
        if (isfn)
        {
            if (grepl("\\.(rda|RData)$", res.savefn, ignore.case=TRUE))
            {
                if (verbose)
                    cat("Save to '", res.savefn, "' ...\n", sep="")
                cm <- switch(res.compress, LZMA="xz", LZMA_RA="xz",
                    ZIP="gzip", ZIP_RA="gzip", none="gzip", TRUE)
                .res <- ans
                save(.res, file=res.savefn, compress=cm)
                if (verbose) cat("Done.\n")
                invisible()
            } else {
                stop("Unknown format of the output file, and it should be RData or gds.")
            }
        } else {
            if (verbose) cat("Done.\n")
            ans
        }
    }
}
