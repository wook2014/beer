#' Derive point estimates for c, pi, phi, and Z for a particular sample
#'
#' Posterior means are used as point estimates for \eqn{c}, \eqn{\pi},
#' \eqn{\phi}, and \eqn{Z}. As super-enriched peptides are tossed out before
#' MCMC sampling, super-enriched peptides return \code{NA} for the \eqn{\phi}
#' and \eqn{Z} point estimates. Indices corresponding to a particular peptide in
#' the MCMC sampler are mapped back to the original peptide names.
#'
#' @param object a \code{\link[PhIPData]{PhIPData}} object
#' @param file path to rds file
#' @param se.matrix logical matrix indicating which peptides were identified as
#' super-enriched peptides
#' @param burn.in number of iterations to be burned
#' @param post.thin thinning parameter
#'
#' @return list of point estimates for c, pi, phi and Z
summarizeRunOne <- function(object, file, se.matrix,
    burn.in = 0, post.thin = 1, run.matrix = NULL) {
    sample <- sub("\\.rds$", "", basename(file))
    if (is.null(run.matrix)) run.matrix <- !se.matrix

    rds <- readRDS(file)

    ## translate peptide indices
    pep_ind <- rep(NA, nrow(object))
    run_ind <- run.matrix[, sample] & !se.matrix[, sample]
    pep_ind[which(run_ind)] <- seq(sum(run_ind))
    names(pep_ind) <- rownames(object)

    point_phi_value <- rep(NA_real_, nrow(object))
    point_phi_enriched <- rep(NA_real_, nrow(object))
    point_Z_value <- rep(NA_real_, nrow(object))
    point_phi_value[!se.matrix[, sample] & !run_ind] <- 1
    point_Z_value[!se.matrix[, sample] & !run_ind] <- 0

    if (!is.list(rds) || !isTRUE(rds$no_mcmc)) {
        mcmc_matrix <- as.matrix(rds)
        iter_ind <- seq(burn.in + 1, nrow(mcmc_matrix), by = post.thin)
        mcmc_matrix <- mcmc_matrix[iter_ind, , drop = FALSE]

        # for convenience extract parameter specific samples
        samples_c <- mcmc_matrix[, grepl("c", colnames(mcmc_matrix)),
            drop = FALSE
        ]
        samples_pi <- mcmc_matrix[, grepl("pi", colnames(mcmc_matrix)),
            drop = FALSE
        ]
        samples_phi <- mcmc_matrix[, grepl("phi\\[", colnames(mcmc_matrix)),
            drop = FALSE
        ]
        samples_Z <- mcmc_matrix[, grepl("Z\\[", colnames(mcmc_matrix)),
            drop = FALSE
        ]

        point_phi_value[run_ind] <- unname(colMeans(samples_phi))
        point_phi_enriched[run_ind] <- unname(
            colSums(samples_phi * samples_Z) / colSums(samples_Z)
        )
        point_Z_value[run_ind] <- unname(colMeans(samples_Z))
        point_c_value <- mean(samples_c)
        point_pi_value <- mean(samples_pi)
    } else {
        point_c_value <- NA_real_
        point_pi_value <- 0
    }

    # summarize info
    point_c <- data.frame(
        parameter = "c",
        sample = sample,
        est_value = point_c_value
    )
    point_pi <- data.frame(
        parameter = "pi",
        sample = sample,
        est_value = point_pi_value
    )
    point_phi <- data.frame(
        parameter = "phi",
        sample = sample,
        peptide = rownames(object),
        est_value = point_phi_value,
        est_enriched = point_phi_enriched
    )
    point_Z <- data.frame(
        parameter = "Z",
        sample = sample,
        peptide = rownames(object),
        est_value = point_Z_value
    )

    list(
        point_c = point_c,
        point_pi = point_pi,
        point_phi = point_phi,
        point_Z = point_Z
    )
}

#' Summarize MCMC chain and return point estimates for BEER parameters
#'
#' Posterior means are used as point estimates for \eqn{c}, \eqn{\pi},
#' \eqn{\phi}, and \eqn{Z}. As super-enriched peptides are tossed out before
#' MCMC sampling, super-enriched peptides return \code{NA} for the \eqn{\phi}
#' and \eqn{Z} point estimates. Indices corresponding to a particular peptide in
#' the MCMC sampler are mapped back to the original peptide names.
#'
#' @param object a \code{\link[PhIPData]{PhIPData}} object
#' @param jags.files list of files containing MCMC sampling results
#' @param se.matrix logical matrix indicating which peptides were identified as
#' super-enriched peptides
#' @param burn.in number of iterations to be burned
#' @param post.thin thinning parameter
#' @param run.matrix logical matrix indicating which peptides were run in JAGS
#' @param assay.names named vector of specifying where to store point estimates
#' @param BPPARAM \code{[BiocParallel::BiocParallelParam]} passed to
#' BiocParallel functions.
#'
#' @return PhIPData object with point estimates stored in the assays specified
#' by `assay.names`.
#'
#' @importFrom progressr handlers progressor
#' @importFrom BiocParallel bplapply
#' @import PhIPData SummarizedExperiment
summarizeRun <- function(object, jags.files, se.matrix,
    burn.in = 0, post.thin = 1,
    assay.names = c(
        phi = NULL, phi_Z = "logfc", Z = "prob",
        c = "sampleInfo", pi = "sampleInfo"
    ),
    run.matrix = NULL,
    BPPARAM = BiocParallel::bpparam()) {
    if (is.null(run.matrix)) run.matrix <- !se.matrix

    ## Check that all files are present
    if (!all(file.exists(jags.files))) {
        stop(paste0(
            "Cannot find the following files: ",
            paste0(jags.files[!file.exists(jags.files)],
                collapse = ", "
            )
        ))
    }

    samples <- sub("\\.rds$", "", basename(jags.files))
    names(jags.files) <- samples

    ## Pre-allocate containers
    point_c <- if (assay.names["c"] %in% colnames(sampleInfo(object))) {
        sampleInfo(object)[[assay.names["c"]]]
    } else {
        rep(NA, ncol(object))
    }
    point_pi <- if (assay.names["pi"] %in% colnames(sampleInfo(object))) {
        sampleInfo(object)[[assay.names["pi"]]]
    } else {
        rep(NA, ncol(object))
    }
    point_phi <- if (!is.null(assay.names["phi"]) &
        assay.names["phi"] %in% assayNames(object)) {
        assay(object, assay.names["phi"])
    } else {
        matrix(NA, nrow = nrow(object), ncol = ncol(object))
    }
    point_phi <- if (!is.null(assay.names["phi"]) &
        assay.names["phi"] %in% assayNames(object)) {
        assay(object, assay.names["phi"])
    } else {
        matrix(NA, nrow = nrow(object), ncol = ncol(object))
    }
    point_phi_Z <- if (!is.null(assay.names["phi_Z"]) &
        assay.names["phi_Z"] %in% assayNames(object)) {
        assay(object, assay.names["phi_Z"])
    } else {
        matrix(NA, nrow = nrow(object), ncol = ncol(object))
    }
    point_Z <- if (!is.null(assay.names["Z"]) &
        assay.names["Z"] %in% assayNames(object)) {
        assay(object, assay.names["Z"])
    } else {
        matrix(NA, nrow = nrow(object), ncol = ncol(object))
    }

    names(point_c) <- names(point_pi) <- colnames(point_phi) <-
        colnames(point_phi_Z) <- colnames(point_Z) <-
        colnames(object)
    rownames(point_phi) <- rownames(point_phi_Z) <- rownames(point_Z) <-
        rownames(object)

    ## Summarize each file, use bplapply here to enable reading/import
    ## to be parallelized
    progressr::handlers("txtprogressbar")
    p <- progressr::progressor(along = jags.files)
    .summarizeRunOne <- summarizeRunOne
    files_out <- bplapply(jags.files, function(file) {
        file_counter <- paste0(
            which(file == jags.files), " of ",
            length(jags.files)
        )
        p(file_counter, class = "sticky", amount = 1)
        .summarizeRunOne(object, file, se.matrix, burn.in, post.thin,
            run.matrix = run.matrix)
    }, BPPARAM = BPPARAM)

    for (out in files_out) {
        sample <- out$point_c$sample

        point_c[sample] <- out$point_c$est_value
        point_pi[sample] <- out$point_pi$est_value
        point_phi[, sample] <- out$point_phi$est_value
        point_phi_Z[, sample] <- out$point_phi$est_enriched
        point_Z[, sample] <- out$point_Z$est_value
    }

    ## Assign c and pi to sampleInfo
    if (!is.na(assay.names["c"])) object$c <- point_c
    if (!is.na(assay.names["pi"])) object$pi <- point_pi

    ## Assign phi, phi_Z, and Z to assays
    assay <- c("phi", "phi_Z", "Z")[!is.na(assay.names[c("phi", "phi_Z", "Z")])]
    assays(object)[assay.names[assay]] <- list(
        phi = point_phi, phi_Z = point_phi_Z,
        Z = point_Z
    )[assay]
    object
}
