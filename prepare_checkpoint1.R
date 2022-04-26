prepare_checkpoint1 <- function(raw_exp=NULL, trans_exp=NULL, file_metadata=NULL, marker_metadata=NULL, cell_metadata=NULL, out_dir=NULL, arcsinh_cofactor=5) {
  if(is.null(trans_exp)){
    if(is.null(raw_exp)) {
      stop("At least one of raw expression matrix and arcsinh-transformed expression matrix is required.")
    } else {
      trans_exp <- asinh( raw_exp / arcsinh_cofactor )
    }
  }
  if(is.null(file_metadata)){
    stop("File matadata is required.")
  }
  if(is.null(marker_metadata)){
    stop("Marker matadata is required.")
  }
  if(is.null(cell_metadata)){
    stop("Cell matadata is required.")
  }
  if(!all(rownames(cell_metadata) == rownames(trans_exp))) {
    stop("Rownames of cell_metadata, raw_exp and trans_exp must match.")
  }
  if(sum(marker_metadata$marker_name %in% colnames(trans_exp)) > 0) {
    overlap_markers <- marker_metadata$marker_name[ marker_metadata$marker_name %in% colnames(trans_exp) ]
    cat("Following", length(overlap_markers), "markers overlap between expression data and marker_metadata\n")
    cat(paste0(overlap_markers, collapse = ", "))
  } else {
    stop("No marker_name match column names of the expression matrix")
  }
  CHECKPOINT <- 1
  save(CHECKPOINT, raw_exp, trans_exp, cell_metadata, file_metadata, marker_metadata, file=file.path(out_dir,"checkpoint_1.RData"))
}

