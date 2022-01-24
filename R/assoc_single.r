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


#######################################################################
# Internal model initialization

.init_nullmod <- function(modobj, ii, maf, mac, missing, spa.pval, var.ratio,
    summac=NaN, wbeta=double(), num_wbuf=0L, acatv_mac=10)
{
    # check
    if (!is.finite(var.ratio))
        stop("Invalid variance ratio in the SAIGE model.")
    # initialize the internal model parameters
    y <- unname(modobj$obj.noK$y)
    mu <- unname(modobj$fitted.values)
    X1 <- modobj$obj.noK$X1[ii,, drop=FALSE]
    n <- length(ii)
    mobj <- list(
        maf = maf, mac = mac, missing = missing, spa.pval = spa.pval,
        tau = modobj$tau,
        y = y[ii], mu = mu[ii],
        y_mu = (y - mu)[ii],  # y - mu
        mu2 = (mu * (1 - mu))[ii],
        t_XXVX_inv = t(modobj$obj.noK$XXVX_inv[ii,, drop=FALSE]),  # K x n_samp (K << n_samp, more efficient)
        XV = modobj$obj.noK$XV[, ii, drop=FALSE],  # K x n_samp
        t_XVX_inv_XV = t(modobj$obj.noK$XXVX_inv[ii,, drop=FALSE] * modobj$obj.noK$V[ii]),  # K x n_samp
        t_X = t(X1),  # K x n_samp
        var.ratio = var.ratio,
        # buffer
        buf_dosage = double(n),
        buf_coeff = double(NROW(modobj$obj.noK$XV)),
        buf_adj_g = double(n),
        buf_index = integer(n),
        buf_B = double(n),
        buf_g_tilde = double(n),
        buf_X1 = double(NCOL(X1)),
        buf_spa = double(n+n)
    )
    # additional design matrix
    if (modobj$trait.type == "binary")
    {
        mobj$XVX <- t(X1) %*% (X1 * mobj$mu2)  # a matrix: K x K
        mobj$S_a <- colSums(X1 * mobj$y_mu)    # a vector of size K
    } else if (modobj$trait.type == "quantitative")
    {
        mobj$XVX <- t(X1) %*% X1               # a matrix: K x K
        mobj$S_a <- colSums(X1 * mobj$y_mu)    # a vector of size K
    } else
        stop("Invalid 'modobj$trait.type': ", modobj$trait.type, ".")
    # additional for aggregate tests
    mobj$summac <- summac
    mobj$acatv_mac <- acatv_mac
    mobj$buf_wbeta <- as.double(wbeta)
    mobj$buf_unitsz <- double(num_wbuf*5L)
    # output
    mobj
}

.dsnode <- function(gdsfile, nm)
{
    if (nm == "")
    {
        n <- index.gdsn(gdsfile, "genotype/data", silent=TRUE)
        if (!is.null(n))
        {
            nm <- "$dosage_alt"
        } else {
            nm <- getOption("seqarray.node_ds", "annotation/format/DS")
            n <- index.gdsn(gdsfile, nm, silent=TRUE)
            if (is.null(n))
                stop("Dosages should be stored in genotype or annotation/format/DS.")
        }
    }
    nm
}


#######################################################################
# SAIGE single variant analysis
#

seqAssocGLMM_SPA <- function(gdsfile, modobj, maf=NaN, mac=10, missing=0.1,
    dsnode="", spa.pval=0.05, var.ratio=NaN, res.savefn="", res.compress="LZMA",
    parallel=FALSE, verbose=TRUE)
{
    stopifnot(inherits(gdsfile, "SeqVarGDSClass") | is.character(gdsfile))
    stopifnot(is.numeric(maf), length(maf)==1L)
    stopifnot(is.numeric(mac), length(mac)==1L)
    stopifnot(is.numeric(missing), length(missing)==1L)
    stopifnot(is.character(dsnode), length(dsnode)==1L, !is.na(dsnode))
    stopifnot(is.character(dsnode), length(dsnode)==1L, !is.na(dsnode))
    stopifnot(is.numeric(spa.pval), length(spa.pval)==1L)
    stopifnot(is.numeric(var.ratio), length(var.ratio)==1L)
    stopifnot(is.character(res.savefn), length(res.savefn)==1L)
    if (!is.element(res.compress, c("LZMA", "LZMA_RA", "ZIP", "ZIP_RA", "none")))
        stop("`res.compress` should be one of LZMA, LZMA_RA, ZIP, ZIP_RA and none.")
    stopifnot(is.logical(verbose), length(verbose)==1L)

    if (verbose)
        cat(.crayon_inverse("SAIGE association analysis:\n"))

    # check model
    modobj <- .check_modobj(modobj, verbose)

    # GDS file
    if (is.character(gdsfile))
    {
        if (verbose)
            .cat("    open ", sQuote(gdsfile))
        gdsfile <- seqOpen(gdsfile, allow.duplicate=TRUE)
        on.exit(seqClose(gdsfile))
    } else {
        # save the filter on GDS file
        seqSetFilter(gdsfile, action="push", verbose=FALSE)
        on.exit(seqSetFilter(gdsfile, action="pop", verbose=FALSE))
    }

    # show warnings immediately
    saveopt <- options(warn=1L)
    on.exit(options(warn=saveopt$warn), add=TRUE)

    # determine the GDS node for dosages
    dsnode <- .dsnode(gdsfile, dsnode)

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
        .cat("    # of samples: ", .pretty(dm[2L]))
        .cat("    # of variants: ", .pretty(dm[3L]))
        .cat("    MAF threshold: ", maf)
        .cat("    MAC threshold: ", mac)
        .cat("    missing threshold for variants: ", missing)
        .cat("    p-value threshold for SPA adjustment: ", spa.pval)
    }

    if (!is.finite(var.ratio))
        var.ratio <- mean(modobj$var.ratio$ratio, na.rm=TRUE)
    if (verbose)
        .cat("    variance ratio for approximation: ", var.ratio)

    if (dm[2L] <= 0) stop("No sample in the genotypic data set!")
    if (dm[3L] <= 0) stop("No variant in the genotypic data set!")

    # initialize the internal model parameters
    mobj <- .init_nullmod(modobj, ii, maf, mac, missing, spa.pval, var.ratio)

    # is forking or not?
    is_fork <- SeqArray:::.IsForking(parallel)
    njobs <- SeqArray:::.NumParallel(parallel)
    if (verbose)
        .cat("    # of processes: ", njobs)

    # initialize internally
    if (njobs<=1L || is_fork)
    {
        # forking, no need to distribute model parameters
        .Call(saige_score_test_init, mobj)
        initfun <- finalfun <- NULL
    } else {
        # pass the model parameters to each process
        if (verbose)
            cat("Distribute the model parameters to the", njobs, "processes\n")
        # initialize child processes internally
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
                seqApply(f, dsnode, .cfunction("saige_score_test_bin"),
                    as.is="list", parallel=FALSE, .progress=pverbose,
                    .list_dup=FALSE, .useraw=NA)
            }, dsnode=dsnode, pverbose=verbose & (njobs==1L))
    } else if (modobj$trait.type == "quantitative")
    {
        rv <- seqParallel(parallel, gdsfile, split="by.variant",
            .initialize=initfun, .finalize=finalfun, .initparam=mobj,
            .balancing=TRUE, .bl_size=50000L, .bl_progress=verbose,
            FUN = function(f, dsnode, pverbose)
            {
                seqApply(f, dsnode, .cfunction("saige_score_test_quant"),
                    as.is="list", parallel=FALSE, .progress=pverbose,
                    .list_dup=FALSE, .useraw=NA)
            }, dsnode=dsnode, pverbose=verbose & (njobs==1L))
    } else
        stop("Invalid 'modobj$trait.type'.")

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
        .cat(
        "# of variants after filtering by MAF, MAC and missing thresholds: ",
            .pretty(length(rv)))
    }

    # output to a GDS file?
    isfn <- !is.na(res.savefn) && res.savefn!=""
    if (isfn && grepl("\\.gds$", res.savefn, ignore.case=TRUE))
    {
        if (verbose)
            .cat("Save to ", sQuote(res.savefn), " ...")
        cm <- res.compress[1L]
        # create a GDS file
        outf <- createfn.gds(res.savefn)
        on.exit(closefn.gds(outf), add=TRUE)
        put.attr.gdsn(outf$root, "FileFormat", "SAIGE_OUTPUT")
        put.attr.gdsn(outf$root, "Version",
            paste0("SAIGEgds_", packageVersion("SAIGEgds")))
        # add function
        Add <- function(varnm, val)
            add.gdsn(outf, varnm, val, compress=cm, closezip=TRUE)
        # add sample IDs
        Add("sample.id", seqGetData(gdsfile, "sample.id"))
        # write data
        .write_gds(outf, "id", gdsfile, "variant.id", cm)
        .write_gds(outf, "chr", gdsfile, "chromosome", cm)
        .write_gds(outf, "pos", gdsfile, "position", cm)
        # rs.id
        if (!is.null(index.gdsn(gdsfile, "annotation/id", silent=TRUE)))
            .write_gds(outf, "rs.id", gdsfile, "annotation/id", cm)
        # ref and alt alleles
        Add("ref", seqGetData(gdsfile, "$ref"))
        Add("alt", seqGetData(gdsfile, "$alt"))
        # other data
        Add("AF.alt", sapply(rv, `[`, i=1L))
        Add("mac", sapply(rv, `[`, i=2L))
        Add("num", as.integer(sapply(rv, `[`, i=3L)))
        Add("beta", sapply(rv, `[`, i=4L))
        Add("SE", sapply(rv, `[`, i=5L))
        Add("pval", sapply(rv, `[`, i=6L))
        if (modobj$trait.type == "binary")
        {
            Add("p.norm", sapply(rv, `[`, i=7L))
            Add("converged", sapply(rv, `[`, i=8L)==1)
        }
        if (verbose) .cat(.crayon_inverse("Done."))
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
        ans$mac <- sapply(rv, `[`, i=2L)
        ans$num  <- as.integer(sapply(rv, `[`, i=3L))
        ans$beta <- sapply(rv, `[`, i=4L)
        ans$SE   <- sapply(rv, `[`, i=5L)
        ans$pval <- sapply(rv, `[`, i=6L)
        if (modobj$trait.type == "binary")
        {
            ans$p.norm <- sapply(rv, `[`, i=7L)
            ans$converged <- as.logical(sapply(rv, `[`, i=8L))
        }

        # save file?
        if (isfn)
        {
            cm <- switch(res.compress, LZMA="xz", LZMA_RA="xz", ZIP="gzip",
                ZIP_RA="gzip", TRUE)
            if (verbose)
                .cat("Save to ", sQuote(res.savefn), " ...")
            if (grepl("\\.(rda|RData)$", res.savefn, ignore.case=TRUE))
            {
                .res <- ans
                save(.res, file=res.savefn, compress=cm)
            } else if (grepl("\\.rds$", res.savefn, ignore.case=TRUE))
            {
                saveRDS(ans, file=res.savefn, compress=cm)
            } else {
                stop("Unknown format of the output file, and it should be RData, RDS or gds.")
            }
            if (verbose) .cat(.crayon_inverse("Done."))
            invisible()
        } else {
            if (verbose) .cat(.crayon_inverse("Done."))
            ans
        }
    }
}



#######################################################################
# SAIGE single variant analysis using user-defined genotypes
#

.UserGLMM_SPA <- function(nVariant, modobj, fun, spa.pval=0.05, var.ratio=NaN,
    res.savefn="", res.compress="LZMA", parallel=FALSE, verbose=TRUE, ...)
{
    stopifnot(is.numeric(nVariant), length(nVariant)==1L, nVariant>0L)
    stopifnot(is.function(fun))
    stopifnot(is.numeric(spa.pval), length(spa.pval)==1L)
    stopifnot(is.numeric(var.ratio), length(var.ratio)==1L)
    stopifnot(is.character(res.savefn), length(res.savefn)==1L)
    if (!is.element(res.compress, c("LZMA", "LZMA_RA", "ZIP", "ZIP_RA", "none")))
        stop("`res.compress` should be one of LZMA, LZMA_RA, ZIP, ZIP_RA and none.")
    stopifnot(is.logical(verbose), length(verbose)==1L)

    # check model
    if (is.character(modobj))
    {
        stopifnot(length(modobj)==1L)
        modobj <- get(load(modobj))
    }
    stopifnot(inherits(modobj, "ClassSAIGE_NullModel"))

    # show warnings immediately
    saveopt <- options(warn=1L)
    on.exit(options(warn=saveopt$warn), add=TRUE)

    nSamp <- length(modobj$sample.id)
    if (nSamp <= 0) stop("No sample!")
    if (nVariant <= 0) stop("No variant!")

    if (verbose)
    {
        .cat(.crayon_inverse("SAIGE association analysis:"))
        .cat("    # of samples: ", .pretty(nSamp))
        .cat("    # of variants: ", .pretty(nVariant))
        .cat("    p-value threshold for SPA adjustment: ", spa.pval)
    }

    if (!is.finite(var.ratio))
        var.ratio <- mean(modobj$var.ratio$ratio, na.rm=TRUE)
    if (verbose)
        .cat("    variance ratio for approximation: ", var.ratio)

    # initialize the internal model parameters
    mobj <- .init_nullmod(modobj, ii, 0, 0, 1, spa.pval, var.ratio)

    # is forking or not?
    is_fork <- SeqArray:::.IsForking(parallel)
    njobs <- SeqArray:::.NumParallel(parallel)
    if (verbose) .cat("    # of processes: ", njobs)

    # initialize internally
    if (njobs<=1L || is_fork)
    {
        # forking, no need to distribute model parameters
        .Call(saige_score_test_init, mobj)
        initfun <- finalfun <- NULL
    } else {
        # pass the model parameters to each process
        if (verbose)
            cat("Distribute the model parameters to the", njobs, "processes\n")
        # initialize child processes internally
        initfun <- function(proc_id, mobj)
        {
            eval(.load_pkg)
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

    # split
    pssplit <- SeqArray:::.file_split(nVariant, njobs, multiple=1L)

    # scan all (selected) variants
    if (modobj$trait.type == "binary")
    {
        rv <- seqParallel(parallel, NULL, split="none",
            .initialize=initfun, .finalize=finalfun, .initparam=mobj,
            FUN = function(fun, pssplit, verbose)
            {
                i <- SeqArray:::process_index
                st <- pssplit[[1L]][i] - 1L
                n <- pssplit[[2L]][i]
                ans <- vector("list", n)
                f <- .cfunction("saige_score_test_bin")
                if (n > 0L)
                {
                    pg <- if (verbose && i==1L) SeqArray:::.seqProgress(n, 1L) else NULL
                    for (j in seq_len(n))
                    {
                        v <- f(fun(j+st))
                        if (!is.null(v)) ans[[j]] <- v
                        SeqArray:::.seqProgForward(pg, 1L)
                    }
                    remove(pg)
                }
                ans
            }, fun=fun, pssplit=pssplit, verbose=verbose)
    } else if (modobj$trait.type == "quantitative")
    {
        rv <- seqParallel(parallel, NULL, split="none",
            .initialize=initfun, .finalize=finalfun, .initparam=mobj,
            FUN = function(fun, pssplit, verbose)
            {
                i <- SeqArray:::process_index
                st <- pssplit[[1L]][i] - 1L
                n <- pssplit[[2L]][i]
                ans <- vector("list", n)
                f <- .cfunction("saige_score_test_quant")
                if (n > 0L)
                {
                    pg <- if (verbose && i==1L) SeqArray:::.seqProgress(n, 1L) else NULL
                    for (j in seq_len(n))
                    {
                        v <- f(fun(j+st))
                        if (!is.null(v)) ans[[j]] <- v
                        SeqArray:::.seqProgForward(pg, 1L)
                    }
                    remove(pg)
                }
                ans
            }, fun=fun, pssplit=pssplit, verbose=verbose)
    } else {
        stop("Invalid 'modobj$trait.type'.")
    }

    # if any maf/mac filter
    if (length(rv) != nVariant)
        stop("Internal error: seqParallel() returns a vector of wrong length.")
    x <- sapply(rv, is.null)
    ii <- which(!x)
    if (any(x)) rv <- rv[ii]
    if (verbose)
    {
        cat("# of variants after filtering by MAF, MAC and missing thresholds: ",
            .pretty(length(rv)), "\n", sep="")
    }

    # output to a GDS file?
    isfn <- !is.na(res.savefn) && res.savefn!=""
    if (isfn && grepl("\\.gds$", res.savefn, ignore.case=TRUE))
    {
        if (verbose)
            .cat("Save to ", sQuote(res.savefn), " ...")
        cm <- res.compress[1L]
        # create a GDS file
        outf <- createfn.gds(res.savefn)
        on.exit(closefn.gds(outf), add=TRUE)
        put.attr.gdsn(outf$root, "FileFormat", "SAIGE_OUTPUT")
        put.attr.gdsn(outf$root, "Version",
            paste0("SAIGEgds_", packageVersion("SAIGEgds")))
        # add function
        Add <- function(varnm, val)
            add.gdsn(outf, varnm, val, compress=cm, closezip=TRUE)
        # write data
        Add("sample.id", modobj$sample.id)
        Add("id", ii)
        Add("AF.alt", sapply(rv, `[`, i=1L))
        Add("mac", sapply(rv, `[`, i=2L))
        Add("num", as.integer(sapply(rv, `[`, i=3L)))
        Add("beta", sapply(rv, `[`, i=4L))
        Add("SE", sapply(rv, `[`, i=5L))
        Add("pval", sapply(rv, `[`, i=6L), compress=cm, closezip=TRUE)
        if (modobj$trait.type == "binary")
        {
            Add("p.norm", sapply(rv, `[`, i=7L))
            Add("converged", sapply(rv, `[`, i=8L)==1L)
        }
        if (verbose) cat(.crayon_inverse("Done.\n"))
        # output
        invisible()
    } else {
        # output
        ans <- data.frame(
            id = ii,
            AF.alt = sapply(rv, `[`, i=1L),
            mac = sapply(rv, `[`, i=2L),
            num = as.integer(sapply(rv, `[`, i=3L)),
            beta = sapply(rv, `[`, i=4L),
            SE   = sapply(rv, `[`, i=5L),
            pval = sapply(rv, `[`, i=6L))
        if (modobj$trait.type == "binary")
        {
            ans$p.norm <- sapply(rv, `[`, i=7L)
            ans$converged <- as.logical(sapply(rv, `[`, i=8L))
        }

        # save file?
        if (isfn)
        {
            if (grepl("\\.(rda|RData)$", res.savefn, ignore.case=TRUE))
            {
                if (verbose)
                    .cat("Save to ", sQuote(res.savefn), " ...")
                cm <- switch(res.compress, LZMA="xz", LZMA_RA="xz",
                    ZIP="gzip", ZIP_RA="gzip", none="gzip", TRUE)
                .res <- ans
                save(.res, file=res.savefn, compress=cm)
                if (verbose) cat(.crayon_inverse("Done.\n"))
                invisible()
            } else {
                stop("Unknown format of the output file, and it should be RData or gds.")
            }
        } else {
            if (verbose) cat(.crayon_inverse("Done.\n"))
            ans
        }
    }
}
