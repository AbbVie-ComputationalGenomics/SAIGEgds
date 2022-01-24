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
# Beta parameters in aggregate tests
#

AggrParamBeta <- structure(c(1,1,1,25), dim=c(2L,2L),
    dimnames=list(c("shape1", "shape2"), c("b_1_1", "b_1_25")))

.check_beta <- function(wbeta)
{
    stopifnot(is.numeric(wbeta), is.finite(wbeta))
    err <- "'wbeta' should be a length-two vector or a matrix with two rows."
    if (is.vector(wbeta))
    {
        if (length(wbeta) != 2L) stop(err)
    } else if (is.matrix(wbeta))
    {
        if (NCOL(wbeta) <= 0L) stop(err)
    }
    invisible()
}

.show_maf <- function(gdsfile, parallel)
{
    cat("Calculating minor allele frequencies (MAF):\n")
    maf <- seqAlleleFreq(gdsfile, minor=TRUE, verbose=TRUE, parallel=parallel)
    cat(sprintf(
        "    MAF: avg (%.5f), min (%.5f), max (%.5f), sd (%.5f)\n",
        mean(maf, na.rm=TRUE), min(maf, na.rm=TRUE),
        max(maf, na.rm=TRUE), sd(maf, na.rm=TRUE)))
    invisible()
}


#######################################################################
# SAIGE burden tests
#

seqAssocGLMM_spaBurden <- function(gdsfile, modobj, units, wbeta=AggrParamBeta,
    summac=3, dsnode="", spa.pval=0.05, var.ratio=NaN, res.savefn="",
    res.compress="LZMA", parallel=FALSE, verbose=TRUE, verbose.maf=TRUE)
{
    stopifnot(inherits(gdsfile, "SeqVarGDSClass") | is.character(gdsfile))
    stopifnot(inherits(units, "SeqUnitListClass"))
    .check_beta(wbeta)
    stopifnot(is.numeric(summac), length(summac)==1L, is.finite(summac))
    stopifnot(is.character(dsnode), length(dsnode)==1L, !is.na(dsnode))
    stopifnot(is.numeric(spa.pval), length(spa.pval)==1L)
    stopifnot(is.numeric(var.ratio), length(var.ratio)==1L)
    stopifnot(is.character(res.savefn), length(res.savefn)==1L)
    stopifnot(is.character(res.compress), length(res.compress)==1L)
    if (!is.element(res.compress, c("LZMA", "LZMA_RA", "ZIP", "ZIP_RA", "none")))
        stop("`res.compress` should be one of LZMA, LZMA_RA, ZIP, ZIP_RA and none.")
    stopifnot(is.logical(verbose), length(verbose)==1L)
    stopifnot(is.logical(verbose.maf), length(verbose.maf)==1L)

    if (verbose)
        .cat(.crayon_inverse("SAIGE burden analysis:"))

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

    # determine the GDS node for dosages
    dsnode <- .dsnode(gdsfile, dsnode)

    # check sample ID
    suppressWarnings(seqSetFilter(gdsfile, sample.id=modobj$sample.id,
        verbose=FALSE))
    sid <- seqGetData(gdsfile, "sample.id")
    if (length(sid) != length(modobj$sample.id))
        stop("Some of sample IDs are not available in the GDS file.")
    ii <- match(sid, modobj$sample.id)
    if (any(is.na(ii)))
        stop("Sample IDs do not match.")

    # show summary
    suppressWarnings(seqSetFilter(gdsfile, variant.sel=unlist(units$index),
        verbose=FALSE))
    mm <- seqSummary(gdsfile, "genotype", verbose=FALSE)
    dm <- mm$seldim
    v <- lengths(units$index)
    sz_wmax <- max(v)
    if (verbose)
    {
        .cat("    # of samples: ", .pretty(dm[2L]))
        .cat("    # of variants in total: ", .pretty(dm[3L]))
        .cat("    # of units: ", .pretty(length(units$index)))
        .cat("    avg. # of variants per unit: ", mean(v))
        .cat("    min # of variants in a unit: ", min(v))
        .cat("    max # of variants in a unit: ", sz_wmax)
        .cat("    sd  # of variants in a unit: ", sd(v))
        .cat("    p-value threshold for SPA adjustment: ", spa.pval)
    }
    remove(v)

    if (!is.finite(var.ratio))
        var.ratio <- mean(modobj$var.ratio$ratio, na.rm=TRUE)
    if (verbose)
        .cat("    variance ratio for approximation: ", var.ratio)

    if (!is.matrix(wbeta))
        wbeta <- matrix(wbeta, nrow=2L)
    wb_colnm <- sprintf("b%g_%g", wbeta[1L,], wbeta[2L,])
    if (verbose)
    {
        v <- apply(wbeta, 2L, function(x)
            paste0("beta(", x[1L], ",", x[2L], ")"))
        .cat("    variant weights: ", paste(v, collapse=", "))
    }

    if (dm[2L] <= 0) stop("No sample in the genotypic data set!")
    if (dm[3L] <= 0) stop("No variant in the genotypic data set!")

    # is forking or not?
    is_fork <- SeqArray:::.IsForking(parallel)
    njobs <- SeqArray:::.NumParallel(parallel)
    if (verbose)
        .cat("    # of processes: ", njobs)

    # get allele frequencies
    if (verbose && isTRUE(verbose.maf)) .show_maf(gdsfile, parallel)

    # initialize the internal model parameters
    mobj <- .init_nullmod(modobj, ii, 0, 0, 1, spa.pval, var.ratio,
        summac, wbeta, sz_wmax)

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
        # initialize
        seqParallel(parallel, NULL, split="none", .combine="none",
            FUN = function(mobj) {
                eval(parse(text="library(Rcpp)"))
                eval(parse(text="library(SAIGEgds)"))
                .packageEnv$mobj <- mobj
                .Call(saige_score_test_init, mobj)
            }, mobj=mobj)
        # finalize
        on.exit({
            seqParallel(parallel, NULL, split="none", .combine="none",
                FUN = function() { .packageEnv$mobj <- NULL })
        }, add=TRUE)
    }

    # scan all variant units
    if (verbose)
        cat("Calculating p-values:\n")
    if (modobj$trait.type == "binary")
    {
        rv <- seqUnitApply(gdsfile, units, dsnode,
            FUN=function(x) .Call(saige_burden_test_bin, x), as.is="list",
            parallel=parallel, .useraw=NA, .progress=verbose)
    } else if (modobj$trait.type == "quantitative")
    {
        rv <- seqUnitApply(gdsfile, units, dsnode,
            FUN=function(x) .Call(saige_burden_test_quant, x), as.is="list",
            parallel=parallel, .useraw=NA, .progress=verbose)
    } else
        stop("Invalid 'modobj$trait.type'.")

    if (length(rv) != length(units$index))
        stop("seqUnitApply() returns a vector of wrong length.")

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
        put.attr.gdsn(outf$root, "FileFormat", "SAIGE_OUTPUT_SET")
        put.attr.gdsn(outf$root, "Version",
            paste0("SAIGEgds_", packageVersion("SAIGEgds")))
        # add function
        Add <- function(varnm, val)
            add.gdsn(outf, varnm, val, compress=cm, closezip=TRUE)
        # add sample IDs
        Add("sample.id", seqGetData(gdsfile, "sample.id"))
        # write summary variant data
        for (nm in names(units$desp))
            Add(nm, units$desp[[nm]])
        Add("numvar", lengths(units$index))
        Add("maf.avg", sapply(rv, `[`, i=1L))
        Add("maf.sd",  sapply(rv, `[`, i=2L))
        Add("maf.min", sapply(rv, `[`, i=3L))
        Add("maf.max", sapply(rv, `[`, i=4L))
        Add("mac.avg", sapply(rv, `[`, i=5L))
        Add("mac.sd",  sapply(rv, `[`, i=6L))
        Add("mac.min", sapply(rv, `[`, i=7L))
        Add("mac.max", sapply(rv, `[`, i=8L))
        # write p-values
        k <- 8L
        for (i in seq_len(ncol(wbeta)))
        {
            s <- ""
            if (length(wb_colnm) > 1L)
                s <- paste0(".", wb_colnm[i])
            Add(paste0("summac", s), sapply(rv, `[`, i=k+1L))
            Add(paste0("beta", s),   sapply(rv, `[`, i=k+2L))
            Add(paste0("SE", s),     sapply(rv, `[`, i=k+3L))
            Add(paste0("pval", s),   sapply(rv, `[`, i=k+4L))
            k <- k + 4L
            if (modobj$trait.type == "binary")
            {
                Add(paste0("p.norm", s), sapply(rv, `[`, i=k+1L))
                Add(paste0("cvg", s), as.logical(sapply(rv, `[`, i=k+2L)))
                k <- k + 2L
            }
        }
        if (verbose) cat(.crayon_inverse("Done.\n"))
        # output
        invisible()

    } else {
        # output
        ans <- units$desp
        ans$numvar  <- lengths(units$index)
        ans$maf.avg <- sapply(rv, `[`, i=1L)
        ans$maf.sd  <- sapply(rv, `[`, i=2L)
        ans$maf.min <- sapply(rv, `[`, i=3L)
        ans$maf.max <- sapply(rv, `[`, i=4L)
        ans$mac.avg <- sapply(rv, `[`, i=5L)
        ans$mac.sd  <- sapply(rv, `[`, i=6L)
        ans$mac.min <- sapply(rv, `[`, i=7L)
        ans$mac.max <- sapply(rv, `[`, i=8L)
        k <- 8L
        for (i in seq_len(ncol(wbeta)))
        {
            s <- ""
            if (length(wb_colnm) > 1L)
                s <- paste0(".", wb_colnm[i])
            ans[[paste0("summac", s)]] <- sapply(rv, `[`, i=k+1L)
            ans[[paste0("beta", s)]]   <- sapply(rv, `[`, i=k+2L)
            ans[[paste0("SE", s)]]     <- sapply(rv, `[`, i=k+3L)
            ans[[paste0("pval", s)]]   <- sapply(rv, `[`, i=k+4L)
            k <- k + 4L
            if (modobj$trait.type == "binary")
            {
                ans[[paste0("p.norm", s)]] <- sapply(rv, `[`, i=k+1L)
                ans[[paste0("cvg", s)]] <- as.logical(sapply(rv, `[`, i=k+2L))
                k <- k + 2L
            }
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
            if (verbose) cat(.crayon_inverse("Done.\n"))
            ans
        }
    }
}



#######################################################################
# SAIGE ACAT-V tests
#

seqAssocGLMM_spaACAT_V <- function(gdsfile, modobj, units, wbeta=AggrParamBeta,
    burden.mac=10, burden.summac=3, dsnode="", spa.pval=0.05, var.ratio=NaN,
    res.savefn="", res.compress="LZMA", parallel=FALSE,
    verbose=TRUE, verbose.maf=TRUE)
{
    stopifnot(inherits(gdsfile, "SeqVarGDSClass") | is.character(gdsfile))
    stopifnot(inherits(units, "SeqUnitListClass"))
    .check_beta(wbeta)
    stopifnot(is.numeric(burden.mac), length(burden.mac)==1L,
        is.finite(burden.mac))
    stopifnot(is.numeric(burden.summac), length(burden.summac)==1L,
        is.finite(burden.summac))
    stopifnot(is.character(dsnode), length(dsnode)==1L, !is.na(dsnode))
    stopifnot(is.numeric(spa.pval), length(spa.pval)==1L)
    stopifnot(is.numeric(var.ratio), length(var.ratio)==1L)
    stopifnot(is.character(res.savefn), length(res.savefn)==1L)
    stopifnot(is.character(res.compress), length(res.compress)==1L)
    if (!is.element(res.compress, c("LZMA", "LZMA_RA", "ZIP", "ZIP_RA", "none")))
        stop("`res.compress` should be one of LZMA, LZMA_RA, ZIP, ZIP_RA and none.")
    stopifnot(is.logical(verbose), length(verbose)==1L)
    stopifnot(is.logical(verbose.maf), length(verbose.maf)==1L)

    if (verbose)
        .cat(.crayon_inverse("SAIGE ACAT-V analysis:"))

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

    # determine the GDS node for dosages
    dsnode <- .dsnode(gdsfile, dsnode)

    # check sample ID
    suppressWarnings(seqSetFilter(gdsfile, sample.id=modobj$sample.id,
        verbose=FALSE))
    sid <- seqGetData(gdsfile, "sample.id")
    if (length(sid) != length(modobj$sample.id))
        stop("Some of sample IDs are not available in the GDS file.")
    ii <- match(sid, modobj$sample.id)
    if (any(is.na(ii)))
        stop("Sample IDs do not match.")

    # show summary
    suppressWarnings(seqSetFilter(gdsfile, variant.sel=unlist(units$index),
        verbose=FALSE))
    mm <- seqSummary(gdsfile, "genotype", verbose=FALSE)
    dm <- mm$seldim
    v <- lengths(units$index)
    sz_wmax <- max(v)
    if (verbose)
    {
        .cat("    # of samples: ", .pretty(dm[2L]))
        .cat("    # of variants in total: ", .pretty(dm[3L]))
        .cat("    # of units: ", .pretty(length(units$index)))
        .cat("    avg. # of variants per unit: ", mean(v))
        .cat("    min # of variants in a unit: ", min(v))
        .cat("    max # of variants in a unit: ", sz_wmax)
        .cat("    sd  # of variants in a unit: ", sd(v))
        .cat("    p-value threshold for SPA adjustment: ", spa.pval)
    }
    remove(v)

    if (!is.finite(var.ratio))
        var.ratio <- mean(modobj$var.ratio$ratio, na.rm=TRUE)
    if (verbose)
        .cat("    variance ratio for approximation: ", var.ratio)

    if (!is.matrix(wbeta))
        wbeta <- matrix(wbeta, nrow=2L)
    wb_colnm <- sprintf("v%g_%g", wbeta[1L,], wbeta[2L,])
    if (verbose)
    {
        v <- apply(wbeta, 2L, function(x)
            paste0("beta(", x[1L], ",", x[2L], ")"))
        .cat("    variant weights: ", paste(v, collapse=", "))
        .cat("    mac threshold for burden test: ", burden.mac)
    }

    if (dm[2L] <= 0) stop("No sample in the genotypic data set!")
    if (dm[3L] <= 0) stop("No variant in the genotypic data set!")

    # is forking or not?
    is_fork <- SeqArray:::.IsForking(parallel)
    njobs <- SeqArray:::.NumParallel(parallel)
    if (verbose)
        .cat("    # of processes: ", njobs)

    # get allele frequencies
    if (verbose && isTRUE(verbose.maf)) .show_maf(gdsfile, parallel)

    # initialize the internal model parameters
    mobj <- .init_nullmod(modobj, ii, 0, 0, 1, spa.pval, var.ratio,
        burden.summac, wbeta, sz_wmax, burden.mac)

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
        # initialize
        seqParallel(parallel, NULL, split="none", .combine="none",
            FUN = function(mobj) {
                eval(parse(text="library(Rcpp)"))
                eval(parse(text="library(SAIGEgds)"))
                .packageEnv$mobj <- mobj
                .Call(saige_score_test_init, mobj)
            }, mobj=mobj)
        # finalize
        on.exit({
            seqParallel(parallel, NULL, split="none", .combine="none",
                FUN = function() { .packageEnv$mobj <- NULL })
        }, add=TRUE)
    }

    # scan all variant units
    if (verbose)
        cat("Calculating p-values:\n")
    if (modobj$trait.type == "binary")
    {
        rv <- seqUnitApply(gdsfile, units, dsnode,
            FUN=function(x) .Call(saige_acatv_test_bin, x), as.is="list",
            parallel=parallel, .useraw=NA, .progress=verbose)
    } else if (modobj$trait.type == "quantitative")
    {
        rv <- seqUnitApply(gdsfile, units, dsnode,
            FUN=function(x) .Call(saige_acatv_test_quant, x), as.is="list",
            parallel=parallel, .useraw=NA, .progress=verbose)
    } else
        stop("Invalid 'modobj$trait.type'.")

    if (length(rv) != length(units$index))
        stop("seqUnitApply() returns a vector of wrong length.")

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
        put.attr.gdsn(outf$root, "FileFormat", "SAIGE_OUTPUT_SET")
        put.attr.gdsn(outf$root, "Version",
            paste0("SAIGEgds_", packageVersion("SAIGEgds")))
        # add function
        Add <- function(varnm, val)
            add.gdsn(outf, varnm, val, compress=cm, closezip=TRUE)
        # add sample IDs
        Add("sample.id", seqGetData(gdsfile, "sample.id"))
        # write summary variant data
        for (nm in names(units$desp))
            Add(nm, units$desp[[nm]])
        Add("numvar", lengths(units$index))
        Add("maf.avg", sapply(rv, `[`, i=1L))
        Add("maf.sd",  sapply(rv, `[`, i=2L))
        Add("maf.min", sapply(rv, `[`, i=3L))
        Add("maf.max", sapply(rv, `[`, i=4L))
        Add("mac.avg", sapply(rv, `[`, i=5L))
        Add("mac.sd",  sapply(rv, `[`, i=6L))
        Add("mac.min", sapply(rv, `[`, i=7L))
        Add("mac.max", sapply(rv, `[`, i=8L))
        Add("n.single", as.integer(sapply(rv, `[`, i=9L)))
        Add("n.burden", as.integer(sapply(rv, `[`, i=10L)))
        # write p-values
        st <- 10L
        for (i in seq_len(ncol(wbeta)))
        {
            s <- ""
            if (length(wb_colnm) > 1L)
                s <- paste0(".", wb_colnm[i])
            Add(paste0("pval", s),  sapply(rv, `[`, i=st+1L))
            Add(paste0("p.med", s), sapply(rv, `[`, i=st+2L))
            Add(paste0("p.min", s), sapply(rv, `[`, i=st+3L))
            Add(paste0("p.max", s), sapply(rv, `[`, i=st+4L))
            st <- st + 4L
        }
        if (verbose) cat(.crayon_inverse("Done.\n"))
        # output
        invisible()

    } else {
        # output
        ans <- units$desp
        ans$numvar  <- lengths(units$index)
        ans$maf.avg <- sapply(rv, `[`, i=1L)
        ans$maf.sd  <- sapply(rv, `[`, i=2L)
        ans$maf.min <- sapply(rv, `[`, i=3L)
        ans$maf.max <- sapply(rv, `[`, i=4L)
        ans$mac.avg <- sapply(rv, `[`, i=5L)
        ans$mac.sd  <- sapply(rv, `[`, i=6L)
        ans$mac.min <- sapply(rv, `[`, i=7L)
        ans$mac.max <- sapply(rv, `[`, i=8L)
        ans$n.single <- as.integer(sapply(rv, `[`, i=9L))
        ans$n.burden <- as.integer(sapply(rv, `[`, i=10L))
        st <- 10L
        for (i in seq_len(ncol(wbeta)))
        {
            s <- ""
            if (length(wb_colnm) > 1L)
                s <- paste0(".", wb_colnm[i])
            ans[[paste0("pval", s)]]  <- sapply(rv, `[`, i=st+1L)
            ans[[paste0("p.med", s)]] <- sapply(rv, `[`, i=st+2L)
            ans[[paste0("p.min", s)]] <- sapply(rv, `[`, i=st+3L)
            ans[[paste0("p.max", s)]] <- sapply(rv, `[`, i=st+4L)
            st <- st + 4L
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
            if (verbose) cat(.crayon_inverse("Done.\n"))
            ans
        }
    }
}



#######################################################################
# SAIGE ACAT-O tests
#

seqAssocGLMM_spaACAT_O <- function(gdsfile, modobj, units, wbeta=AggrParamBeta,
    burden.mac=10, burden.summac=3, dsnode="", spa.pval=0.05, var.ratio=NaN,
    res.savefn="", res.compress="LZMA", parallel=FALSE,
    verbose=TRUE, verbose.maf=TRUE)
{
    stopifnot(inherits(gdsfile, "SeqVarGDSClass") | is.character(gdsfile))
    stopifnot(inherits(units, "SeqUnitListClass"))
    .check_beta(wbeta)
    stopifnot(is.numeric(burden.mac), length(burden.mac)==1L,
        is.finite(burden.mac))
    stopifnot(is.numeric(burden.summac), length(burden.summac)==1L,
        is.finite(burden.summac))
    stopifnot(is.character(dsnode), length(dsnode)==1L, !is.na(dsnode))
    stopifnot(is.numeric(spa.pval), length(spa.pval)==1L)
    stopifnot(is.numeric(var.ratio), length(var.ratio)==1L)
    stopifnot(is.character(res.savefn), length(res.savefn)==1L)
    stopifnot(is.character(res.compress), length(res.compress)==1L)
    if (!is.element(res.compress, c("LZMA", "LZMA_RA", "ZIP", "ZIP_RA", "none")))
        stop("`res.compress` should be one of LZMA, LZMA_RA, ZIP, ZIP_RA and none.")
    stopifnot(is.logical(verbose), length(verbose)==1L)
    stopifnot(is.logical(verbose.maf), length(verbose.maf)==1L)

    if (verbose)
        .cat(.crayon_inverse("SAIGE ACAT-O analysis (combining burden & ACAT-V):"))

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

    # determine the GDS node for dosages
    dsnode <- .dsnode(gdsfile, dsnode)

    # check sample ID
    suppressWarnings(seqSetFilter(gdsfile, sample.id=modobj$sample.id,
        verbose=FALSE))
    sid <- seqGetData(gdsfile, "sample.id")
    if (length(sid) != length(modobj$sample.id))
        stop("Some of sample IDs are not available in the GDS file.")
    ii <- match(sid, modobj$sample.id)
    if (any(is.na(ii)))
        stop("Sample IDs do not match.")

    # show summary
    suppressWarnings(seqSetFilter(gdsfile, variant.sel=unlist(units$index),
        verbose=FALSE))
    mm <- seqSummary(gdsfile, "genotype", verbose=FALSE)
    dm <- mm$seldim
    v <- lengths(units$index)
    sz_wmax <- max(v)
    if (verbose)
    {
        .cat("    # of samples: ", .pretty(dm[2L]))
        .cat("    # of variants in total: ", .pretty(dm[3L]))
        .cat("    # of units: ", .pretty(length(units$index)))
        .cat("    avg. # of variants per unit: ", mean(v))
        .cat("    min # of variants in a unit: ", min(v))
        .cat("    max # of variants in a unit: ", sz_wmax)
        .cat("    sd  # of variants in a unit: ", sd(v))
        .cat("    p-value threshold for SPA adjustment: ", spa.pval)
    }
    remove(v)

    if (!is.finite(var.ratio))
        var.ratio <- mean(modobj$var.ratio$ratio, na.rm=TRUE)
    if (verbose)
        .cat("    variance ratio for approximation: ", var.ratio)

    if (!is.matrix(wbeta))
        wbeta <- matrix(wbeta, nrow=2L)
    wb_colnm <- sprintf("%g_%g", wbeta[1L,], wbeta[2L,])
    if (verbose)
    {
        v <- apply(wbeta, 2L, function(x)
            paste0("beta(", x[1L], ",", x[2L], ")"))
        .cat("    variant weights: ", paste(v, collapse=", "))
        .cat("    mac threshold for burden test: ", burden.mac)
    }

    if (dm[2L] <= 0) stop("No sample in the genotypic data set!")
    if (dm[3L] <= 0) stop("No variant in the genotypic data set!")

    # is forking or not?
    is_fork <- SeqArray:::.IsForking(parallel)
    njobs <- SeqArray:::.NumParallel(parallel)
    if (verbose)
        .cat("    # of processes: ", njobs)

    # get allele frequencies
    if (verbose && isTRUE(verbose.maf)) .show_maf(gdsfile, parallel)

    # initialize the internal model parameters
    mobj <- .init_nullmod(modobj, ii, 0, 0, 1, spa.pval, var.ratio,
        burden.summac, wbeta, sz_wmax, burden.mac)

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
        # initialize
        seqParallel(parallel, NULL, split="none", .combine="none",
            FUN = function(mobj) {
                eval(parse(text="library(Rcpp)"))
                eval(parse(text="library(SAIGEgds)"))
                .packageEnv$mobj <- mobj
                .Call(saige_score_test_init, mobj)
            }, mobj=mobj)
        # finalize
        on.exit({
            seqParallel(parallel, NULL, split="none", .combine="none",
                FUN = function() { .packageEnv$mobj <- NULL })
        }, add=TRUE)
    }

    # scan all variant units
    if (verbose)
        cat("Calculating p-values:\n")
    if (modobj$trait.type == "binary")
    {
        rv <- seqUnitApply(gdsfile, units, dsnode,
            FUN=function(x) .Call(saige_acato_test_bin, x), as.is="list",
            parallel=parallel, .useraw=NA, .progress=verbose)
    } else if (modobj$trait.type == "quantitative")
    {
        rv <- seqUnitApply(gdsfile, units, dsnode,
            FUN=function(x) .Call(saige_acato_test_quant, x), as.is="list",
            parallel=parallel, .useraw=NA, .progress=verbose)
    } else
        stop("Invalid 'modobj$trait.type'.")

    if (length(rv) != length(units$index))
        stop("seqUnitApply() returns a vector of wrong length.")

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
        put.attr.gdsn(outf$root, "FileFormat", "SAIGE_OUTPUT_SET")
        put.attr.gdsn(outf$root, "Version",
            paste0("SAIGEgds_", packageVersion("SAIGEgds")))
        # add function
        Add <- function(varnm, val)
            add.gdsn(outf, varnm, val, compress=cm, closezip=TRUE)
        # add sample IDs
        Add("sample.id", seqGetData(gdsfile, "sample.id"))
        # write summary variant data
        for (nm in names(units$desp))
            Add(nm, units$desp[[nm]])
        Add("numvar", lengths(units$index))
        Add("maf.avg", sapply(rv, `[`, i=1L))
        Add("maf.sd",  sapply(rv, `[`, i=2L))
        Add("maf.min", sapply(rv, `[`, i=3L))
        Add("maf.max", sapply(rv, `[`, i=4L))
        Add("mac.avg", sapply(rv, `[`, i=5L))
        Add("mac.sd",  sapply(rv, `[`, i=6L))
        Add("mac.min", sapply(rv, `[`, i=7L))
        Add("mac.max", sapply(rv, `[`, i=8L))
        Add("pval",    sapply(rv, `[`, i=9L))
        # write p-values
        for (i in seq_len(ncol(wbeta)))
        {
            s <- wb_colnm[i]
            Add(paste0("pval.b", s), sapply(rv, `[`, i=8L+i*2L))
            Add(paste0("pval.v", s), sapply(rv, `[`, i=9L+i*2L))
        }
        if (verbose) cat(.crayon_inverse("Done.\n"))
        # output
        invisible()

    } else {
        # output
        ans <- units$desp
        ans$numvar  <- lengths(units$index)
        ans$maf.avg <- sapply(rv, `[`, i=1L)
        ans$maf.sd  <- sapply(rv, `[`, i=2L)
        ans$maf.min <- sapply(rv, `[`, i=3L)
        ans$maf.max <- sapply(rv, `[`, i=4L)
        ans$mac.avg <- sapply(rv, `[`, i=5L)
        ans$mac.sd  <- sapply(rv, `[`, i=6L)
        ans$mac.min <- sapply(rv, `[`, i=7L)
        ans$mac.max <- sapply(rv, `[`, i=8L)
        ans$pval    <- sapply(rv, `[`, i=9L)
        for (i in seq_len(ncol(wbeta)))
        {
            s <- wb_colnm[i]
            ans[[paste0("pval.b", s)]] <- sapply(rv, `[`, i=8L+i*2L)
            ans[[paste0("pval.v", s)]] <- sapply(rv, `[`, i=9L+i*2L)
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
            if (verbose) cat(.crayon_inverse("Done.\n"))
            ans
        }
    }
}
