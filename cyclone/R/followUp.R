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
    for (i in seq_len(...length())) {
        this_package <- ...elt(i)
        if (!requireNamespace(this_package, quietly = TRUE)) {
            missing <- c(missing, this_package)
        }
    }

    if (length(missing) > 0) {
        stop(
            "Missing package(s) required for ", fxn, ": ",
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
    .timestamped_msg("Saving Complete.", verbose = verbose)
}

#' A function to automatically import cyclone outputs into a SingleCellExperiment structure OR load one that has already been generated and saved.
#' @param sce_file_name String. File path of the SingleCellExperiment (SCE) to load in, if previously created.
#' @param checkpoint1,checkpoint8 Strings. File paths pointing to cyclone outputs 'checkpoint_1.RData' and 'checkpoint_8.RData' to use when \code{load_checkpoints = TRUE}.
#' @param load_checkpoints Logical. Whether or not to load in checkpoint1 and checkpoint8 data.
#' @param make_clusters_factors Logical. When creating the SCE object, whether to ensure numeric 'cluster_#x#' metadata are factors with levels in numeric order.
#' Doing so ensures clusters appear in order from 1,2,3,4,5, etc. when plotting.
#' @param verbose Logical. Whether to output timestamped log messages during running.
#' @param verbose_checkpoint_load Logical. Whether to set 'verbose = TRUE' in \code{\link[base]{load}} calls for loading 'checkpoint1' and 'checkpoint8' data.
#' @return a \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @details The function starts by loading in primary outputs of the cyclone pipeline if \code{load_checkpoints = TRUE}.
#' When doing so, it will read in \code{checkpoint1} before \code{checkpoint8},
#' a necessary ordering because the cell_metadata element of \code{checkpoint1} is updated to include clustering information by \code{checkpoint8}.
#' Next, if a file exists at the specified \code{sce_file_name} location, that file will be loaded in and is assumed to be an SCE previously created from the cyclone data.
#' Otherwise, a new SingleCellExperiment (SCE) object is created where: \itemize{
#' \item The 'trans_exp' output from checkpoint1 is used to fill a 'transformed' assay.
#' \item The 'cell_metadata' output from checkpoint8 is used to fill both  colData and a reducedDim named 'umap'.
#' \item If \code{make_clusters_factors} is \code{TRUE}, all colData columns whose names match with 'cluster_#x#' are turned into factors with levels ordered from min value to max value.
#' }
#' Finally, the SCE is returned.
#' @author Daniel Bunis
#' @export
make_or_load_full_sce <- function(
    sce_file_name = NULL,
    checkpoint1,
    checkpoint8,
    load_checkpoints = TRUE,
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

    if ( !identical(sce_file_name, NULL) && file.exists(sce_file_name) ) {
        .timestamped_msg("Reading in previously made Rds file, ", sce_file_name, verbose = verbose)
        sce <- readRDS(sce_file_name)
    } else {
        if (!load_checkpoints) {
            stop("No file at 'sce_file_name', but checkpoint data was not read in to make one.")
        }
        .timestamped_msg("No file at 'sce_file_name', Making SCE.", verbose = verbose)
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

        if (make_clusters_factors) {
            .check_packages(
                "SummarizedExperiment", # Another dep of SCE
                fxn = "turning numeric cluster metadata into factors")
            .timestamped_msg("Turning numeric cluster metadata into factors", verbose = verbose)
            for (res in grep("^cluster_(\\d)+x(\\d)+$", colnames(SummarizedExperiment::colData(sce)), value = TRUE)) {
                this_clusts <- sce[[res, drop = TRUE]]
                sce[[res]] <- factor(
                    sce[[res, drop = TRUE]],
                    levels = min(this_clusts, na.rm = TRUE):max(this_clusts, na.rm = TRUE)
                )
            }
        }
    }

    .timestamped_msg("Done.", verbose = verbose)
    sce
}

#' A function to create a down-sampled SingleCellExperiment object from a full one.
#' @section NOTE:
#' performs a simple downsample that does NOT attempt to pull equally from each
#' sample or cluster.  The purpose here is assumed to simply be rapid testing of
#' visualizations which can take minutes longer to produce from millions of
#' cells than from thousands of cells.
#' @param full_sce the \code{\link[SingleCellExperiment]{SingleCellExperiment}} object to downsample
#' @param n_keep Positive integer. The number of cells to retain.
#' @param verbose Logical. Whether to output timestamped log messages during running.
#' @return a \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @details \code{full_sce} is subset to a randomly selected \code{n_keep} number of cells, and then returned.
#' @author Daniel Bunis
#' @export
downsample_sce <- function(
    full_sce,
    n_keep = 100000,
    verbose = TRUE
    ) {

    .timestamped_msg("Creating downsampled SCE", verbose = verbose)
    kept_for_downsample <- sample(ncol(full_sce), min(ncol(full_sce), n_keep))
    sce_down <- full_sce[,kept_for_downsample]
    .timestamped_msg("Done.", verbose = verbose)
    sce_down
}

#' Calculate per-sample frequencies of clusters or cell annotations, and compare them across group.
#' @param object A \code{\link[SingleCellExperiment]{SingleCellExperiment}} (or Seurat) object
#' @param cell.by String name of a per-cell metadata (a column of \code{colData(object)}) containing the cluster or cell annotation identities to quantify and assess.
#' @param group.by String name of a per-cell metadata (a column of \code{colData(object)}) containing sample-group identities.
#' @param group.1,group.2 Strings naming the 2 groups within the \code{group.by} metadata which you aim to compare.
#' @param sample.by String name of a per-cell metadata (a column of \code{colData(object)}) containing which sample each cell belongs to.
#' Recommendations for cyclone data, (standardiazed because they are required elements of the file_metadata input!): \itemize{
#' \item 'file_name': holds which original fcs file each cell came from.
#' \item 'donor_id': holds which patient/mouse each cell came from.
#' \item somthing else: sometimes, your data might both break up samples' data acquisition accross multiple .fcs files (so 'file_name' would then be too specific) & contain multiple timepoints or conditions per sample (so 'donor_id' is not specific enough).
#' In such a case, the burden lies on the user to create a viable metadata (\code{<object>$<metadata-name> <- <properly-uniqued-values>}, and then use \code{sample.by = <metadata-name>})
#' \item 'donor_id' & data subsetting with `cells.use`: As an alternative to creating a new metadata to use for \code{sample.by}, subsetting to only cells from a specific timepoint or condition might serve a dual purpose of achieving your specific analysis goal && allowing 'donor_id' to properly scope to individual samples.
#' See \code{cells.use} input description for further details.
#' }
#' @param cell.targs (Optional) Single string or a string vector naming which cell groups of the \code{cell.by} metadata which should be targetted.
#' When not provided, the function will loop through all cell groups in the \code{cell.by} metadata.
#' @param cells.use Logical vector, the same length as the number of cells in the object, which sets which cells to include (TRUE) versus ignore (FALSE).
#' @param pseudocount Number, an ideally small value relative to the lowest expected cell frequencies of the data, which is added to both \code{group.1} and \code{group.2} median frequencies to prevent division by zero in fold_change calculation.
#' @param p.adjust.method String, passed along to the \code{method} input of \code{\link[stats]{p.adjust}}, any valid option for that input will work. "fdr" by default.
#' @param data.out Logical. When set to \code{TRUE}, changes the output from the stats data.frame alone to a named list containing both the stats ("stats") and the underlying per-sample frequency calculations ("data").
#' @return a data.frame, or if \code{data.out} was set to \code{TRUE}, a named list containing 2 data.frames, 'stats' and the underlying 'data'.
#' @details The function starts by utilizing \code{\link[dittoSeq]{dittoFreqPlot}} for
#' \code{cell.by}-cell frequency calculation within \code{sample.by}-samples,
#' percent normalization,
#' marking which \code{sample.by}-samples belong to which \code{group.by}-groups,
#' and trimming to only: 1. requested \code{cell.targs}, 2. \code{group.by}-groups \code{group.1} and \code{group.2}, and 3. cells matching the \code{cells.use} requirements if any were given.
#' It then removes some unnecessary columns from the data.frame returned by \code{\link[dittoSeq]{dittoFreqPlot}}. (Set \code{data.out = TRUE} to see what this cleaned return looks like!)
#'
#' Afterwards, it loops through all \code{cell.targs}, building a row of the eventual stats return for each.
#' Of note, a small \code{pseudocount} is introduced in median fold change calculation to prevent division by zero errors. This \code{pseudocount} has no effect on p-values and only a nominal effect on fold_changes for most cell types, but can be made smaller to decrease the effect on differences among very rare (fraction less than 0.0001) populations.
#' Lastly, \code{p.adjust.method} correction, FDR by default, is applied to the 'p' column and added as a 'padj' column before data is returned.
#' @section The stats data.frame return:
#' Each row holds statistics for an individual comparison.
#' The columns represent:
#' \itemize{
#' \item cell_group: this row's cluster or cell-annotation
#' \item comparison: this groups of \code{group.by} compared in this row, formatted \code{<group.1>_vs_<group.2>}.
#' (For compatibility with running the function multiple times, each targwtting distinct groups, and then concatenating all outputs together!)
#' \item median_g1: the median frequency for the given cell_group within samples from \code{group.1}
#' \item median_g2: the median frequency for the given cell_group within samples from \code{group.2}
#' \item median_fold_change: \code{(median_g1 + pseudocount) / (median_g2 + pseudocount)}. A small \code{pseudocount} is used here to prevent division by zero errors.
#' \item median_log2_fold_chang: \code{log2( median_fold_change )}
#' \item positive_fc_means_up_in: Value = \code{group.1}, just a minor note to help remember the directionality of these fold changes!
#' \item p: The p-value associated with comparison of cell_group percent frequencies of group.1 samples versus group.2 samples using a Mann Whitney U Test / wilcoxon rank sum test (\code{\link[stats]{wilcox.test}}).
#' \item padj: p-values corrected by the chosen \code{p.adjust.method}, FDR by default, built from running \code{p.adjust(stats$p, method = p.adjust.method)} per all hypotheses tested in this call to the \code{freq_stats} function.
#' }
#' @author Daniel Bunis
#' @export
#' @importFrom stats wilcox.test
#' @importFrom stats p.adjust
#' @importFrom stats median
freq_stats <- function(
        object,
        sample.by,
        cell.by,
        group.by, group.1, group.2,
        cell.targs = NULL,
        cells.use = TRUE,
        pseudocount = 1e-6,
        p.adjust.method = "fdr",
        data.out = FALSE
) {

    .check_packages(
        "dittoSeq", # S4Vectors is dep of SCE
        fxn = "this frequency calculation function")

    if (is.null(cell.targs)) {
        cell.targs <- dittoSeq::metaLevels(cell.by, object)
    }

    # Collect stats with dittoSeq
    data <- dittoSeq::dittoFreqPlot(
        object,
        var = cell.by,
        vars.use = cell.targs,
        sample.by = sample.by,
        group.by = group.by,
        cells.use = cells.use & object[[group.by, drop=TRUE]] %in% c(group.1, group.2),
        data.out = TRUE
    )$data

    # Clean
    data$var.data <- NULL # Column only needed for making the plot
    data$grouping <- NULL # Column only needed for making the plot, it's included in a column with the metadata's own name
    data$label.count.total.per.facet <- NULL # Not needed once used for percent calculation
    names(data)[1] <- "cell_group"

    # Here, we loop through all the cell_groups being targeted, 1- calculating stats and 2- building a data.frame during each iteration.
    #  The lapply call performs the iteration, and gathers the data.frames output by each iteration into a list.
    #  That list of data.frames created in our lapply is then 'rbind'ed into a single data.frame.
    stats <- do.call(
        rbind,
        lapply(
            unique(data$cell_group),
            function(clust) {
                data_use <- data[data$cell_group==clust,]
                g1s <- as.vector(data_use[[group.by]]==group.1)
                g2s <- as.vector(data_use[[group.by]]==group.2)
                new <- data.frame(
                    cell_group = clust,
                    comparison = paste0(group.1, "_vs_", group.2),
                    median_g1 = median(data_use$percent[g1s]),
                    median_g2 = median(data_use$percent[g2s]),
                    stringsAsFactors = FALSE
                )
                new$median_fold_change <- (new$median_g1 + pseudocount) / (new$median_g2 + pseudocount)
                new$median_log2_fold_change <- log2(new$median_fold_change)
                new$positive_fc_means_up_in <- group.1
                new$p <- wilcox.test(x=data_use$percent[g1s],
                                     y=data_use$percent[g2s])$p.value
                new
            })
    )

    # Apply FDR correction
    stats$padj <- p.adjust(stats$p, method = p.adjust.method)

    # Output
    if (data.out) {
        list(stats = stats, data = data)
    } else {
        stats
    }
}
