#' output a log message starting with the current system time.
#' @noRd
#' @param ... Strings which should be part of the message.
#' @param verbose Logical. Whether the message should be output at all. Allows for control of verbosity of any functions that call on this function
.timestamped_msg <- function(..., verbose = TRUE) {
    if (verbose) {
        cat("[ ", format(Sys.time()), " ] ", ..., "\n", sep = "")
    }
}

#' error if packages not installed
#' @noRd
#' @param ... package names to check for installation
#' @param fxn string naming the user-requested action that requires the packages
.check_packages <- function(..., fxn = "the requested functionality") {

    missing <- c()
    for (i in ...length) {
        this_package <- ...elt(i)
        if (!requireNamespace(this_package, quietly = TRUE)) {
            missing <- c(missing, this_package)
        }
    }

    if (length(missing) > 0) {
        stop(
            "Missing package(s) required for", fxn, ":",
            paste(missing, collapse = ", ")
        )
    }
}

#' Minor wrapper on top of saveRDS that automatically sets whether to compress the object based on cell number.
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param sce_file_name String. A filepath at which to save the \code{sce}.
#' @param compress Logical or NA. Whether to compress the saved object per like-named input of \code{\link[base]{readRDS}}. When left as \code{NA}, compression is performed if the \code{sce} contains more than 1,000,000 cells.
#' @param verbose Logical. Whether to output timestamped log messages during running.
#' @param ... Additional arguments passed to \code{\link[base]{readRDS}}
#' @return None
#' @author Daniel Bunis
#' @export
save_sce <- function(sce, sce_file_name, compress = NA, verbose, ...) {
    if (is.na(compress)) {
        compress <- ncol(sce) > 1000000
    }
    .timestamped_msg(
        ifelse(compress, "Saving SCE, compressed", "Saving SCE"), " to: ",
        sce_file_name, verbose = verbose)
    saveRDS(sce, file = sce_file_name, compress = compress, ...)
}

#' A function to automatically import cyclone outputs into a SingleCellExperiment structure OR load one that has already been generated and saved.
#' @param sce_file_name String. File path of either the sce to load in, or where to save to when \code{save = TRUE}.
#' @param checkpoint1,checkpoint8 Strings. File paths pointing to cyclone outputs 'checkpoint_1.RData' and 'checkpoint_8.RData' to use.
#' @param load_checkpoints Logical. Whether or not to load in checkpoint1 and checkpoint8 data.
#' @param save Logical. Whether to save the SCE object if newly creating it here from checkpoint data
#' @param make_clusters_factors Logical. When creating the SCE object, whether to ensure 'cluster' metadata are factors with levels in numeric order.
#' Doing so ensures clusters appear in order from 1,2,3,4,5, etc. when plotting.
#' @param verbose Logical. Whether to output timestamped log messages during running.
#' @param verbose_checkpoint_load = Logical. Whether to set 'verbose = TRUE' in \code{\link[base]{load}} calls for loading 'checkpoint1' and 'checkpoint8' data.
#' @return a \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @details The function starts by loading in primary outputs of the cyclone pipeline if \code{load_checkpoints = TRUE}.
#' When doing so, it will read in \code{checkpoint1} before \code{checkpoint8},
#' a necessary ordering because the cell_metadata element of \code{checkpoint1} is updated to include clustering information by \code{checkpoint8}.
#' Next, if a file exists at the specified \code{sce_file_name} location, that file will be loaded in and is assumed to be an SCE previously created from the cyclone data.
#' Otherwise, a new SingleCellExperiment (SCE) object is created where: \itemize{
#' \item The 'trans_exp' output from checkpoint1 is used to fill a 'transformed' assay.
#' \item The 'cell_metadata' output from checkpoint8 is used to fill both  colData and a reducedDim named 'umap'.
#' }
#' When \code{make_clusters_factors} is \code{TRUE}, all colData columns whose names start with 'column' are turned into factors with levels ordered from min value to max value.
#' When \code{save} is \code{TRUE} and an SCE was newly created, \code{\link{save_sce}} is then used to save the newly created SCE to the \code{sce_file_name} file path.
#' Finally, the SCE is returned.
#' @author Daniel Bunis
#' @export
make_or_load_full_sce <- function(
    sce_file_name,
    checkpoint1,
    checkpoint8,
    load_checkpoints = TRUE,
    save = FALSE,
    make_clusters_factors = TRUE,
    verbose = TRUE,
    verbose_checkpoint_load = TRUE) {

    if (load_checkpoints) {
        .timestamped_msg("Loading Checkpoint 1", verbose = verbose)
        load(checkpoint1,
             envir = .GlobalEnv,
             verbose = verbose_checkpoint_load)

        .timestamped_msg("Loading Checkpoint 8", verbose = verbose)
        load(checkpoint8,
             envir = .GlobalEnv,
             verbose = verbose_checkpoint_load)
    }

    if ( file.exists(sce_file_name) ) {
        .timestamped_msg("Reading in previously made Rds file, ", sce_file_name, verbose = verbose)
        sce <- readRDS(sce_file_name)
    } else {
        if (!load_checkpoints) {
            stop("No file at 'sce_file_name', but checkpoint data was not read in.")
        }
        .timestamped_msg("Making SCE.", verbose = verbose)
        .check_packages(
            "SingleCellExperiment", "S4Vectors", # S4Vectors is dep of SCE
            fxn = "creating a SingleCellExperiment object")
        sce <- SingleCellExperiment::SingleCellExperiment(
            assays = list(
                transformed=t(trans_exp)
            ),
            colData = S4Vectors::DataFrame(
                cell_metadata[, !grepl("UMAP",colnames(cell_metadata))]
            ),
            reducedDims = list(
                umap=cell_metadata[, grepl("UMAP",colnames(cell_metadata))])
            )
        created_sce <- TRUE
    }

    if (make_clusters_factors) {
        .check_packages(
            "SummarizedExperiment", # Another dep of SCE
            fxn = "making cluster metadata into factors")
        for (res in grep("^cluster", SummarizedExperiment::colData(sce), value = TRUE)) {
            this_clusts <- as.numeric(as.character(sce[[res, drop = TRUE]]))
            sce[[res]] <- factor(
                sce[[res, drop = TRUE]],
                levels = min(this_clusts):max(this_clusts)
            )
        }
    }

    if (created_sce && save) {
        save_sce(sce, sce_file_name)
    }

    .timestamped_msg("Done.", verbose = verbose)
    sce
}

#' A function to create a down-sampled SingleCellExperiment object from a full one OR to load one in that has already been generated and saved.
#' @section NOTE:
#' performs a simple downsample that does NOT attempt to pull equally from each
#' sample or cluster.  The purpose here is assumed to simply be rapid testing of
#' visualizations which can take minutes longer to produce from millions of
#' cells than from thousands of cells.
#' @param down_sce_file_name String. File path of either the sce to load in, or where to save to when \code{save = TRUE}.
#' @param full_sce the \code{\link[SingleCellExperiment]{SingleCellExperiment}} object to downsample
#' @param n_keep Positive integer. The number of cells to retain.
#' @param save Logical. Whether to save the SCE object if newly creating it.
#' @param verbose Logical. Whether to output timestamped log messages during running.
#' @return a \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @details If a file exists at the specified \code{down_sce_file_name} location, that file will be loaded in and returned.
#' Otherwise, \code{full_sce} will be subset to a randomly selected \code{n_keep} number of cells.
#' Then, if \code{save} is \code{TRUE}, \code{\link{save_sce}} is then used to save the downsampled SCE to the \code{down_sce_file_name} file path.
#' Finally, the downsampled SCE is returned.
#' @author Daniel Bunis
#' @export
make_or_load_downsample_sce <- function(
    down_sce_file_name,
    full_sce,
    n_keep = 100000,
    save = TRUE,
    verbose = TRUE
    ) {

    if (file.exists(down_sce_file_name)) {
        .timestamped_msg("Reading in previously made Rds file, ", down_sce_file_name, verbose = verbose)
        sce_down <- readRDS(down_sce_file_name)
    } else {
        .timestamped_msg("Creating downsampled SCE", verbose = verbose)
        kept_for_downsample <- sample(ncol(full_sce), min(ncol(full_sce), 100000))
        sce_down <- full_sce[,kept_for_downsample]
        if (save) {
            save_sce(sce_down, down_sce_file_name)
        }
        .timestamped_msg("Done.", verbose = verbose)
    }
    sce_down
}
