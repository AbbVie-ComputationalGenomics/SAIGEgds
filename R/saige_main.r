#######################################################################
#
# Package name: SAIGEgds
#
# Description:
#     Scalable and accurate implementation of generalized mixed models
# using GDS files
#
# Copyright (C) 2019-2020    Xiuwen Zheng / AbbVie-ComputationalGenomics
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

.pretty_size <- function(x)
{
    if (x >= 1024^4)
        sprintf("%.1fT", x/1024^4)
    else if (x >= 1024^3)
        sprintf("%.1fG", x/1024^3)
    else if (x >= 1024^2)
        sprintf("%.1fM", x/1024^2)
    else if (x >= 1024)
        sprintf("%.1fK", x/1024)
    else
        sprintf("%g bytes", x)
}

SIMD <- function() .Call(saige_simd_version)

.rank_norm <- function(x) qnorm((rank(x) - 0.5)/length(x))

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


# Write to GDS file
.write_gds <- function(out.gds, out.nm, in.gds, in.nm, cm)
{
    i <- index.gdsn(in.gds, in.nm)
    n <- add.gdsn(out.gds, out.nm, storage=i, compress=cm)
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
    geno.sparse=TRUE, num.thread=1L, model.savefn="", seed=200L, fork.loading=FALSE,
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
    stopifnot(is.logical(fork.loading), length(fork.loading)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    if (verbose)
    {
        cat(.crayon_inverse("SAIGE association analysis:\n"))
        cat(.crayon_underline(date()), "\n", sep="")
    }

    rand_seed <- eval(parse(text="set.seed"))
    rand_seed(seed)

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

    # show warnings immediately
    saveopt <- options(warn=1L)
    on.exit(options(warn=saveopt$warn), add=TRUE)

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
        rand_seed(seed)
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
        Xqr <- qr(X)  # QR decomposition
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
    if (SeqArray:::.IsForking(num.thread) && isTRUE(fork.loading))
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
        cat(.pretty_size(as.double(object.size(packed.geno))))
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
            n <- length(v) 
            v <- data.frame(v, as.numeric(prop.table(v)))
            v[, 1L] <- paste0("      ", v[, 1L])
            colnames(v) <- c(phenovar, "Number", "Proportion")
            print(v, row.names=FALSE)
            if (n != 2L)
                stop("The outcome variable has more than 2 categories!")
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
        {
            cat(.crayon_inverse("Calculate the average ratio of variances:\n"))
            cat(.crayon_underline(date()), "\n", sep="")
        }
        rand_seed(seed)
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

        # Inverse normal transformation
        if (isTRUE(inv.norm))
        {
            if (isTRUE(X.transform)) phenovar <- "y"
            fit0 <- glm(formula, data=data)
            resid.sd <- sd(fit0$residuals)
            new.y <- .rank_norm(fit0$residuals) * resid.sd
            data[[phenovar]] <- new.y
            if (verbose)
            {
                cat("Inverse normal transformation on residuals with standard deviation: ",
                    resid.sd, "\n", sep="")
            }
        }

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
        obj.noK$V <- rep(1, length(fit0$y))
        obj.noK$X1 <- X1
        obj.noK$XV <- t(X1)
        obj.noK$XXVX_inv <- X1 %*% solve(t(X1) %*% X1)
        class(obj.noK) <- "SA_NULL"

        # initial tau
        tau <- tau.init
        if (sum(tau) == 0) tau <- c(0.5, 0.5)

        y <- fit0$y
        offset <- fit0$offset
        if (is.null(offset)) offset <- rep(0, length(y))
        eta <- fit0$linear.predictors
        mu <- fit0$fitted.values
        mu.eta <- fit0$family$mu.eta(eta)
        Y <- eta - offset + (y - mu)/mu.eta
        tau <- var(Y) * tau / sum(tau)

        # iterate
        glmm <- .Call(saige_fit_AI_PCG_quant, fit0, X1, tau, param)

        # calculate the variance ratio
        if (verbose)
        {
            cat(.crayon_inverse("Calculate the average ratio of variances:\n"))
            cat(.crayon_underline(date()), "\n", sep="")
        }
        rand_seed(seed)
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
    class(glmm) <- "ClassSAIGE_NullModel"

    if (!is.na(model.savefn) && model.savefn!="")
    {
        cat("Save the model to '", model.savefn, "'\n", sep="")
        .glmm <- glmm
        save(.glmm, file=model.savefn)
    }
    if (verbose)
    {
        cat(.crayon_underline(date()), "\n", sep="")
        cat(.crayon_inverse("Done."), "\n", sep="")
    }

    if (!is.na(model.savefn) && model.savefn!="")
        return(invisible(glmm))
    else
        return(glmm)
}


print.ClassSAIGE_NullModel <- function(x, ...) str(x)
