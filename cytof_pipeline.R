# Objective: process CyTOF data, make UMAP, run clustering, an export various statistics for downstream analysis

 

### NOTES:
# Assumes that:
#      the FCS data doesn't need compensation
#      the FCS files have live and singlets, i.e. normalization beads, multiplets and dead cells are already gated out.
#      the out_dir, data_dir, scaffold_dir contain absolute paths.
#      all FCS files have the same set of channels and the corresponding markers.
#      
# 


################ Define & parse command-line options ###################
suppressMessages({
library(optparse)
})

# Define command-line options
option_list <- list(
  make_option(c("-c", "--config_file"), type="character", default=NULL,
              help="YAML configuration file", metavar="character")
)

# Parse command-line options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
config_file <- opt$config_file

if( is.null(config_file) ) {
  cat("No input is found. Use -h/--help option to read the usage.\n")
  quit(save = "no", status=1)
}
############## END: Define & parse command-line options ################



########################### Loading packages ###########################
suppressMessages({
  ## CyTOF data analysis
  library(flowCore) # FCS I/O
  library(cluster) # Clustering
  library(FlowSOM) # Clustering
  library(uwot) # UMAP
  
  ## Modify tables and text
  library(tidyverse)

  ## Plotting data
  library(viridis)
  library(RColorBrewer)
  library(gridExtra) # plot grid of graphs
  library(pals) #colors
  library(ggrepel)
  library(pheatmap)
  library(ggpubr)
  #library(ggExtra) # When ggExtra is loaded here, the "FORK" mode of doParallel (not PSOCK) fail. So will load ggExtra later while making plots.
  
  library(optparse)
  library(configr)
  library(igraph)
  
  library(doParallel)
  library(foreach)
  library(clusterSim)
})
########################### END: Loading packages ######################



########## Printing a list of important variable names for Joel ########
important_df_names <- c("raw_exp","trans_exp","cell_metadata","file_metadata","marker_metadata")
important_list_names <- c("param_list")
####### END: Printing a list of important variable names for Joel ######







########################### Functions ##################################
scaffold_save <- function(obj, f_name) {
  con <- file(f_name, "wb")
  serialize(obj, con, ascii = F)
  close(con)
}

scaffold_load <- function(f_name) {
  con <- file(f_name, "rb")
  retval <- unserialize(con)
  close(con)
  return(retval)
}

prep_param_list <- function() {
  return( list(
    seed = seed,
    #### OS related
    isWindows = isWindows,
    #### Data processing related
    data_dir = data_dir,
    out_dir = out_dir,
    rerun_from_step = rerun_from_step,
    arcsinh_cofactor = arcsinh_cofactor,
    subsample = subsample,
    subsample_n = subsample_n,
    #### UMAP related
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    spread = spread,
    learning_rate = learning_rate,
    init = init,
    #### Clustering related
    clustering_method = clustering_method,
    clara_params = clara_params,
    flowsom_params = flowsom_params,
    nthreads = nthreads,
    #### SCAFFoLD related
    scaffold_dir = scaffold_dir,
    exclude_controls = exclude_controls
  ) )
}

print_param_list <- function(param_list) {
  cat("Parameters:\n{\n")
  for(i in names(param_list) ) {
    cat(paste("\t", i, "=", param_list[[i]], "\n" ))
  }
  cat("}\n")
}

get_checkpoint_filename <- function(out_dir, CHECKPOINT) {
  return( paste0(out_dir, "/checkpoint_", CHECKPOINT, ".RData") )
}

print_message <- function(message) {
  cat("[", format(Sys.time()), "]", message, "\n")
}

validate_path_and_exit <- function(path, dir=FALSE) {
  path_exists <- FALSE
  if(dir)
    path_exists <- dir.exists(path)
  else
    path_exists <- file.exists(path)
  if( ! path_exists ) {
    print_message(paste(path, "does not exist. Exiting..."))
    quit(save = "no", status=1)
  }
}



#' grid_optimization_plots
#' @description Generate plots for grid optimization results
#' @param flowsom_out a list containing results of FlowSOM clustering for each grid size
grid_optimization_plots <- function(flowsom_out) {
  my_df = data.frame()
  for(i in 1:length(flowsom_out)) {
    obj = flowsom_out[[i]]
    my_df = rbind(my_df, data.frame(xdim = obj$xdim, ydim = obj$ydim, DBI = obj$DBI, ctime = obj$ctime, nclust=obj$xdim * obj$ydim))
  }
  dbi_plot = ggplot(my_df, aes(x=nclust, y=DBI, label=paste0(xdim,"x",ydim))) + 
    geom_line() + geom_point() + 
    geom_text_repel(nudge_x=5, color="grey") + 
    theme_classic() + xlab("Cluster counts") + ylab("DBI")
  time_plot = ggplot(my_df, aes(x=nclust, y=ctime, label=paste0(xdim,"x",ydim))) + 
    geom_line() + geom_point() + 
    geom_text_repel(nudge_x=5, color="grey") + 
    theme_classic() + xlab("Cluster counts") + ylab("Time (min)")
  pdf(file.path(out_dir, "clustering_param_optimization.pdf"), width=15, height=6)
  print(dbi_plot)
  print(time_plot)
  dev.off()
}




#' run_flowsom
#' @description Perform clustering 
#' @param trans_exp_submarkers data frame of transformed expression of markers (columns) to be used for clustering.
#' @param clust_params list containing clustering parameters. See the function call below to understand how this list is made (using prep_param_list function).
#' @param xdim a single xdim value or a vector of xdim values
#' @param xdim a single ydim value or a vector of ydim values, there has to be same number of xdim and ydim values.
#' @param optimize_grid a boolean indicating whether the current function call is for optimizing the grid sizes. If TRUE, the DBI values are calculated for each grid size and a plot of cluster number vs DBI is generated. The clustering labels are not saved. If FALSE, the cluster labels are saved for each grid size; the DBI is not calculated.
run_flowsom <- function(trans_exp_submarkers, clust_params, xdim, ydim, optimize_grid) {

  if(length(xdim) != length(ydim)) {
    print_message("Fatal error: The number of xdim and ydim values are not the same. Make sure the input is correct and rerun the pipeline. Exiting.")
    quit(save = "no", status = 1)
  }
  
  print_message("Preparing the flowFrame for input to FlowSOM: ")
  my_flowframe = flowCore::flowFrame(as.matrix(trans_exp_submarkers))
  
  if(nthreads == -1)
    nthreads = parallel::detectCores() - 1
  
  if(nthreads > 1) {
    my.cluster <- parallel::makeCluster(
      nthreads, 
      type = "FORK"
    )
    # Register the parallel processes. 
    doParallel::registerDoParallel(cl = my.cluster)
  }
  
  flowsom_out <- foreach(i = 1:length(xdim)) %dopar% {
    cstart = Sys.time()
    print_message(paste0("FlowSOM clustering begins for grid: ", xdim[i], "x", ydim[i] ))
    fSOM <- FlowSOM(my_flowframe,
                    # Input options:
                    compensate = FALSE,
                    transform = FALSE,
                    scale = FALSE,
                    # SOM options:
                    xdim = xdim[i], ydim = ydim[i],
                    # Metaclustering options:
                    nClus = clust_params$k,
                    silent = FALSE)
    cend = Sys.time()
    print_message(paste0("Clustering completed for grid: ", xdim[i], "x", ydim[i] ))
    ctime = as.numeric(difftime(cend, cstart, units="mins"))
    
    fSOM$data = NULL
    
    if(clust_params$meta_cluster) {
      clusters = as.vector(GetMetaclusters(fSOM))
    }
    else {
      clusters = as.vector(GetClusters(fSOM))
    }
    
    if(optimize_grid) {
      DBI_obj = index.DB(trans_exp_submarkers, as.numeric(clusters) )
      DBI = DBI_obj$DB
      list("xdim" = xdim[i], "ydim" = ydim[i], "ctime" = ctime, "DBI" = DBI)
    } else {
      list("clusters" = clusters, "clust_obj" = fSOM, "xdim" = xdim[i], "ydim" = ydim[i], "ctime" = ctime)
    }
  }
  
  # Add the grid ids to each element of flowsom_out
  names(flowsom_out) = paste0("cluster_", xdim, "x", ydim)
  
  if(optimize_grid) {
    grid_optimization_plots(flowsom_out = flowsom_out)
    return()
  } else {
    return(flowsom_out)
  }
}









old_run_flowsom <- function(trans_exp_submarkers, clust_params) {
  optimize_grid = clust_params$optimize_grid
  if(optimize_grid & file.exists(clust_params$grid_sizes_file)) {
    grid_sizes = read.csv(grid_sizes, header=T)
    if( is.null(grid_sizes$xdim) | is.null(grid_sizes$ydim) ) {
      print_message("Warning: cannot find xdim and/or ydim columns in the grid sizes file. Using the user-provided xdim and ydim values instead of optimizing the grid size")
      optimize_grid = FALSE
    } else{ 
      xdim = as.numeric(grid_sizes$xdim)
      ydim = as.numeric(grid_sizes$ydim)
      optimize_grid = TRUE
    }
  } else {
    print_message("Using the user-provided xdim and ydim values.")
    xdim = as.numeric(clust_params$xdim)
    ydim = as.numeric(clust_params$ydim)
    optimize_grid = FALSE
  }
  
  # Reset the optimize_grid in clust_params
  clust_params$optimize_grid = optimize_grid
  
  # Add the list of grids being tested in clust_params.
  clust_params$xdim = paste0(xdim, collapse = ",")
  clust_params$ydim = paste0(ydim, collapse = ",")
  
  print_message("Preparing the flowFrame for input to FlowSOM: ")
  my_flowframe = flowCore::flowFrame(as.matrix(trans_exp_submarkers))
  
  if(nthreads == -1)
    nthreads = parallel::detectCores() - 1
  
  if(nthreads > 1) {
    my.cluster <- parallel::makeCluster(
      nthreads, 
      type = "FORK"
    )
    # Register the parallel processes. 
    doParallel::registerDoParallel(cl = my.cluster)
  }
  
  flowsom_out <- foreach(i = 1:length(xdim)) %dopar% {
    cstart = Sys.time()
    print_message(paste0("FlowSOM clustering begins for grid: ", xdim[i], "x", ydim[i] ))
    fSOM <- FlowSOM(my_flowframe,
                    # Input options:
                    compensate = FALSE,
                    transform = FALSE,
                    scale = FALSE,
                    # SOM options:
                    xdim = xdim[i], ydim = ydim[i],
                    # Metaclustering options:
                    nClus = clust_params$k,
                    silent = FALSE)
    cend = Sys.time()
    print_message(paste0("Clustering completed for grid: ", xdim[i], "x", ydim[i] ))
    ctime = as.numeric(difftime(cend, cstart, units="mins"))
    
    fSOM$data = NULL
    
    if(clust_params$meta_cluster) {
      clusters = as.vector(GetMetaclusters(fSOM))
    }
    else {
      clusters = as.vector(GetClusters(fSOM))
    }
    
    DBI_obj = index.DB(trans_exp_submarkers, as.numeric(cell_metadata$cluster) )
    DBI = DBI_obj$DB
    
    if(optimize_grid) {
      list("xdim" = xdim, "ydim" = ydim, "ctime" = ctime, "DBI" = DBI)
    } else {
      list("clusters" = clusters, "clust_obj" = fSOM, "xdim" = xdim, "ydim" = ydim, "ctime" = ctime, "DBI" = DBI)
    }
  }  
  
  if(optimize_grid) {
    grid_optimization_plots(flowsom_out = flowsom_out)
    return()
  } else {
    flowsom_out$clusters = flowsom_out[[1]][["clusters"]]
    flowsom_out$clust_obj = flowsom_out[[1]][["clust_obj"]]
    return(flowsom_out)
  }
}



run_clara <- function(trans_exp_submarkers, clust_params) {
  print_message("Clara clustering begins")
  clara_res <- clara(trans_exp_submarkers, k = clust_params$k, metric = clust_params$metric, samples = clust_params$samples, rngR = TRUE)
  clusters = clara_res$clustering
  
  # Deleting the redundant data from the clara_res object.
  clara_res$data <- NULL  # Deleting the input data from the clara_res object
  clara_res$clustering <- NULL  # Deleting the clustering info from the clara_res object

  return(list("clusters" = clusters, "clust_obj" = clara_res))
}


# This function returns nrow and ncol the plots need to be arranged for the no_of_plots
get_plot_grid_layout <- function(no_of_plots) {
  sq = sqrt(no_of_plots)
  fl = floor(sq)
  if(sq > fl+0.5) {
    nrow = fl + 1
    ncol = fl + 1
  } else if(sq == fl) {
    nrow = fl
    ncol = fl
  } else {
    nrow = fl + 1
    ncol = fl
  }
  if(nrow * ncol < no_of_plots)
    print("The grid size doesn't fit all plots.")
  return(list("nrow"=nrow, "ncol"=ncol))
}


print_step_startup_msg <- function() {
  cat(paste0("\n\n\n##########################################################\nStarting step#: ", CHECKPOINT+1, "\n"))
}


########################### END: Functions #############################




# Load config parameters
config <- read.config(file = config_file)
for( x in names(config) ) {
  assign(x, config[[x]])
}

isWindows = Sys.info()["sysname"] == "Windows"
nthreads = ifelse(isWindows, 1, nthreads)

# Create out_dir if it doesn't exist
if( ! dir.exists( out_dir ) )
  dir.create( out_dir )



########################### Global variables ###########################
# Default step id.
STEP <- 1
# Default checkpoint id.
CHECKPOINT <- 0

# A list of parameter values to be stored at each checkpoints.
param_list <- list()

########################### END: Global variables ######################




############ Load checkpoint data if exists in the out_dir ##############
# Prepare checkpoint data file name.
checkpoint_files <- list.files(out_dir, pattern="checkpoint_\\d+.RData", full.names=T)

# Assign the new checkpoint code based on the number of checkpoint.RData files. If there are two checkpoint RData files, that must mean that the analysis should resume from CHECKPOINT=2
# If there are no checkpoint.RData files, the value of 0 will be assigned to CHECKPOINT
CHECKPOINT <- length( checkpoint_files )

# Update the CHECKPOINT to the user provided value if the user has provided one.
# If the user has provided a checkpoint value for which there are no checkpoint.RData files, we will use the smallest of the user defined checkpoint and the number of checkpoint.RData files.
if( rerun_from_step != -1 ) CHECKPOINT <- min( rerun_from_step - 1, CHECKPOINT )

# If the checkpoint data file exists in the out_dir, load the data.
if( CHECKPOINT > 0 ) {
  # Load all checkpoint files
  for( my_checkpoint in 1:CHECKPOINT ) {
    checkpoint_file = get_checkpoint_filename(out_dir, my_checkpoint)

    if(! file.exists(checkpoint_file) ) {
      # If the checkpoint file corresponding to the current my_checkpoint doesn't exist, that means there is a missing checkpoint file. Ask users to delete the checkpoint files that come after the current checkpoint file and rerun the secript.
      print_message(paste("Fatal error: The", checkpoint_file, "is missing. Delete all checkpoint_*.RData files following the missing file and rerun the script. Exiting."))
      quit(save = "no", status = 1)
    }
    print_message(paste("Loading", checkpoint_file))
    load(checkpoint_file)
  }
}



print_message(paste0( "Starting the analysis from checkpoint=", CHECKPOINT ))

# In case of resuming the analysis, the param_list will be non-empty. In that case, print message to the console warning the user if there are parameters for which new values are provided by the user, compared to what was stored in the last checkpoint.RData file.
old_new_val_msg <- ""
for(x in names(param_list)) {
  if( exists( x ) ) {
    old_val <- param_list[[x]]
    new_val <- get(x)
    if( ! identical(old_val, new_val) ) {
      msg <- paste0("\t", x, ": old value = \"", old_val, "\"; new value = \"", new_val, "\"\n" )
      old_new_val_msg <- paste0( old_new_val_msg, msg )
    }
  }
}
if( old_new_val_msg != "") {
  print_message( paste0("New values were supplied for the following parameters:") )
  print_message( old_new_val_msg )
}

######### END: Load checkpoint data if exist in the out_dir ############



################ Check if the input files exist if the CHECKPOINT is 0 ##################
if(CHECKPOINT == 0) {
  # Exit if the input files don't exist
  validate_path_and_exit(data_dir, dir = TRUE)
  
  # Marker metadata
  validate_path_and_exit(marker_metadata_csvfile, dir = FALSE)
  marker_metadata <- read.csv(marker_metadata_csvfile) 
  
  # File metadata
  validate_path_and_exit(file_metadata_csvfile, dir = FALSE)
  file_metadata <- read.csv(file_metadata_csvfile)
  if(exclude_controls)
    file_metadata <- file_metadata %>% filter( ! control_sample )
}
################ END: Check if the input files exist ##################



################ Set up SCAFFoLD analysis directories ##################
# Define and create scaffold_dir
scaffold_dir <- file.path( out_dir, "scaffold" )
if( dir.exists( scaffold_dir ) & CHECKPOINT < 5 ) {
  print_message( paste0("SCAFFoLD directory exists. Deleting ", scaffold_dir) )
  unlink( scaffold_dir , recursive=T)
  dir.create( scaffold_dir, showWarnings = TRUE )
} else {
  dir.create( scaffold_dir, showWarnings = TRUE )
}

# Don't make SCAFFoLD map, unless the gated populations are provided.
make_scaffold_map <- FALSE
if( dir.exists(gated_fcs_dir) ) {
  dir.create(file.path(scaffold_dir, "gated"), showWarnings = TRUE )
  gated_fcss <- list.files(gated_fcs_dir, full.names = TRUE, pattern = ".fcs$")
  copy_gated_fcss <- file.copy(gated_fcss, file.path(scaffold_dir, "gated"), overwrite = TRUE)
  if( sum(copy_gated_fcss) > 0 ) {
    make_scaffold_map <- TRUE
  } else {
    print_message("Could not copy any FCS files from the gated population directory. SCAFFoLD map generation will be skipped.")
  }
}
############# END: Set up SCAFFoLD analysis directories ################



###### Begin the analysis: read, preprocess and transform the data.
if(CHECKPOINT == 0) {

  print_step_startup_msg()
  print_message("Starting to read FCS file and transform data using following parameters: ")
  print_param_list( prep_param_list() )
  
  # Setting a seed for the randomization used in downsampling.
  set.seed(seed)

  # Read in FCS files, ASINH transform, concatenate data from all files and prepare cell metadata df
  
    fcs_files <- list.files(path = data_dir, pattern = "*.fcs$", recursive = TRUE) # List all FCS files recursively in the data_dir
    fcs_files <- fcs_files[ basename(fcs_files) %in% file_metadata$file_name ]  # Keep only those files that are included in file_metadata

    raw_exp <- data.frame() # Initialize a data frame for storing raw expression data.
    trans_exp <- data.frame() # Initialize a data frame for storing transformed expression data.
    cell_metadata <- data.frame() # Initialize a data frame for storing cell metadata

    # A variable indicating the FCS first file.
    is_first_file = TRUE
    
    for (f in fcs_files) { # Loop through the FCS files in data_dir 
      
      print_message(paste0("Reading in: ", data_dir, "/",f)) # Print file name
      tmp_fcs_data <- read.FCS(paste0(data_dir, "/",f)) # Read FCS file and store in a temporary variable

      # Finding the "desc" values from the first FCS object for the user-specified channel_name in marker_metadata. The "desc" column is stored in marker_metada for using it as column names for the data exported for SCAFFoLD analysis. 
      if(is_first_file) {
        fcs_param_data = tmp_fcs_data@parameters@data %>% mutate( desc = ifelse(is.na(desc), name, desc) )
        marker_metadata = marker_metadata %>% 
          mutate( desc = fcs_param_data[ match( marker_metadata$channel_name, fcs_param_data$name), ]$desc )
        is_first_file = FALSE
      }
      tmp_raw_exp <- tmp_fcs_data@exprs %>% data.frame() %>% dplyr::select(marker_metadata$channel_name) # Convert @exprs to a data frame
      # Remove cells containing Inf values
      tmp_raw_exp <- tmp_raw_exp[!is.infinite(rowSums(tmp_raw_exp)),]

      
      # Subsampling if subsample == TRUE
      if(subsample) {
        tmp_raw_exp <- sample_n(tmp_raw_exp, subsample_n, replace=TRUE)
      }
      
        
      # arcsinh transformation
      tmp_raw_exp_minimal = tmp_raw_exp %>%
        dplyr::select( marker_metadata$channel_name ) %>% # Select markers included in the marker_metadata
        setNames( marker_metadata$marker_name ) # Rename the marker names using the marker_name column
      
      tmp_trans_exp <- asinh( tmp_raw_exp_minimal / arcsinh_cofactor )

      
      # Count total number of cells
      no_of_cells <- nrow(tmp_trans_exp)
      
      # Generate row names by merging the filename with unique cell ids and assign them to the row.names of the raw and transformed expression df
      f_wo_ext <- str_replace(basename(f), ".fcs", "")     # file name without extension
      row_names <- paste0(f_wo_ext, "_", 1:no_of_cells)
      row.names(tmp_raw_exp) <- row_names
      row.names(tmp_trans_exp) <- row_names

      # Append the data from the current FCS file to the raw_exp and trans_exp dataframes.  
      raw_exp <- rbind( raw_exp, tmp_raw_exp)
      trans_exp <- rbind( trans_exp, tmp_trans_exp)
      
      
      # Extract metadata for the file in f and append that to cell_metadata dataframe.
      file_info <- file_metadata %>% filter( file_name == basename(f) )

      cell_metadata <- rbind( cell_metadata, 
                              data.frame( 
                                          file_info[ rep(1, no_of_cells), ], # Associating the file_metadata info for the file to all cells of that file
                                          row.names = row_names
                                          )
                              )

      # Cleaning up large objects that are not needed anymore
      rm(tmp_raw_exp)
      rm(tmp_raw_exp_minimal)
      rm(tmp_trans_exp)
      rm(tmp_fcs_data)
    }
  print_message("Checkpoint #1 reached. Saving the newly generated data.")
  CHECKPOINT <- 1
  param_list <- prep_param_list()
  checkpoint_file = get_checkpoint_filename(out_dir, CHECKPOINT)
  save(CHECKPOINT, raw_exp, trans_exp, cell_metadata, marker_metadata, file_metadata, param_list, file=checkpoint_file )
  print_message(paste("Data is saved in", checkpoint_file) )
}
##############################




###### UMAP calculation
if(CHECKPOINT == 1) {

  print_step_startup_msg()
  print_message("Starting UMAP calculation using following parameters: ")
  print_param_list( prep_param_list() )

  # Setting a seed for the randomization used in UMAP calculation.
  set.seed(seed)

  # Get the used_for_UMAP status for the markers included in trans_exp df.
  use_markers_for_UMAP <- marker_metadata[ match( colnames(trans_exp), marker_metadata$marker_name), ]$used_for_UMAP

  # Subset the columns of trans_exp to retain only those markers that should be used_for_UMAP.
  trans_exp_submarkers <- trans_exp[, use_markers_for_UMAP ]

  print_message("Using the following markers for UMAP: ")
  print(paste( colnames(trans_exp_submarkers), collapse = "," ))
  
  # UMAP calculation
  print_message("Dimensionality Reduction Begins")
  UMAP <- umap(trans_exp_submarkers, n_neighbors = n_neighbors, learning_rate = learning_rate, min_dist = min_dist, spread = spread, init = init)
  print_message("Dimensionality Reduction Completed")

  cell_metadata$UMAP1 <- UMAP[,1]
  cell_metadata$UMAP2 <- UMAP[,2]

  print_message("Checkpoint #2 reached. Saving the newly generated data.")
  CHECKPOINT <- 2
  param_list <- prep_param_list()
  checkpoint_file = get_checkpoint_filename(out_dir, CHECKPOINT)
  save(CHECKPOINT, cell_metadata, param_list, file=checkpoint_file )
  print_message(paste("Data is saved in", checkpoint_file) )
}
##############################



###### (Clustering parameter optimization)
if(CHECKPOINT == 2) {
  
  print_step_startup_msg()
  print_message("Starting the optimization of clustering parameters using following parameters: ")
  print_param_list( prep_param_list() )
  
  # Setting a seed for the randomization used in clustering.
  set.seed(seed)
  
  if( clustering_method == "flowsom") {
    # Check if the grid size file exists, otherwise exit with an error message.
    validate_path_and_exit(flowsom_params$grid_sizes_file)
  
    # Read the grid sizes from the input CSV file and save the sizes in xdim and ydim variables.
    grid_sizes <- read.csv(flowsom_params$grid_sizes_file, header=T)
    if( is.null(grid_sizes$xdim) | is.null(grid_sizes$ydim) ) {
      print_message(paste0("Fatal error: cannot find xdim and/or ydim columns in the grid sizes file. Make sure that ", flowsom_params$grid_sizes_file, " has two columns, xdim and ydim and rerun the pipeline. Exiting." ))
      quit(save = "no", status=1)
    
    } else{ 
      xdim <- as.numeric(grid_sizes$xdim)
      ydim <- as.numeric(grid_sizes$ydim)
    }
  }
  
  # Get the used_for_clustering status for the markers included in trans_exp df.
  use_markers_for_clustering <- marker_metadata[ match( colnames(trans_exp), marker_metadata$marker_name), ]$used_for_clustering
  
  # Subset the columns of trans_exp to retain only those markers that should be used_for_clustering.
  trans_exp_submarkers <- trans_exp[, use_markers_for_clustering ]
  
  print_message("Using the following markers for optimization: ")
  print(paste( colnames(trans_exp_submarkers), collapse = "," ))
  
  # Optimization
  if( clustering_method == "flowsom") {
    print_message("Optimization Begins:")
    optimization_res = run_flowsom(trans_exp_submarkers, flowsom_params, xdim, ydim, optimize_grid = TRUE)
    print_message("Optimization completed.")
  }
  else if( clustering_method == "clara") {
    print_message("There is no optimization for clara clustering.")
  }
  
  print_message("Checkpoint #3 reached. Saving the newly generated data.")
  CHECKPOINT <- 3
  param_list <- prep_param_list()
  checkpoint_file = get_checkpoint_filename(out_dir, CHECKPOINT)
  save(CHECKPOINT, optimization_res, file=checkpoint_file )
  print_message(paste("Data is saved in", checkpoint_file) )
  print_message("NOTE::: The output of optimization has been generated in clustering_param_optimization.pdf in the output directory. Determine the optimal parameters based on the plots in the PDF, input them in the config.yaml file and rerun the script.")
  quit(save = "no", status=0)
}
##############################


###### Clustering
if(CHECKPOINT == 3) {
  
  print_step_startup_msg()
  print_message("Starting the clustering using following parameters: ")
  print_param_list( prep_param_list() )

  # Setting a seed for the randomization used in clustering.
  set.seed(seed)
  
  # Get the used_for_clustering status for the markers included in trans_exp df.
  use_markers_for_clustering <- marker_metadata[ match( colnames(trans_exp), marker_metadata$marker_name), ]$used_for_clustering
  
  # Subset the columns of trans_exp to retain only those markers that should be used_for_clustering.
  trans_exp_submarkers <- trans_exp[, use_markers_for_clustering ]
  
  print_message("Using the following markers for clustering: ")
  print(paste( colnames(trans_exp_submarkers), collapse = "," ))
  
  # Clustering
  print_message("Clustering Begins:")
  if( clustering_method == "flowsom") {
    if( grepl(",", flowsom_params$xdim) ) {
      xdim <- as.numeric(unlist(strsplit(flowsom_params$xdim, ",")))
      ydim <- as.numeric(unlist(strsplit(flowsom_params$ydim, ",")))
    } else {
      xdim <- flowsom_params$xdim
      ydim <- flowsom_params$ydim
    }
    clust_res = run_flowsom(trans_exp_submarkers, flowsom_params, xdim, ydim, optimize_grid = FALSE)
    # Storing clustering results in variables.
    cell_metadata$cluster <- as.character(clust_res[[1]]$clusters)
    cell_metadata = cbind(cell_metadata, 
                          sapply(clust_res, function(x) x$clusters) %>% as.data.frame() 
                          )
    clust_obj = clust_res[[1]]$clust_obj
    cluster_levels = unique( cell_metadata$cluster )
    cluster_metadata <- data.frame(cluster = cluster_levels, row.names = cluster_levels )
  }
  else if( clustering_method == "clara") {
    clust_res = run_clara(trans_exp_submarkers, clara_params)
    # Storing clustering results in variables.
    cell_metadata$cluster <- as.character(clust_res$clusters)
    clust_obj = clust_res$clust_obj
    cluster_levels = unique( cell_metadata$cluster )
    cluster_metadata <- data.frame(cluster = cluster_levels, row.names = cluster_levels )
  }
  print_message("Clustering Completed")


  print_message("Checkpoint #4 reached. Saving the newly generated data.")
  CHECKPOINT <- 4
  param_list <- prep_param_list()
  checkpoint_file = get_checkpoint_filename(out_dir, CHECKPOINT)
  save(CHECKPOINT, cell_metadata, cluster_metadata, param_list, clust_obj, file=checkpoint_file )
  print_message(paste("Data is saved in", checkpoint_file) )
}
##############################




if(CHECKPOINT == 4) {
  print_step_startup_msg()
  ###### Calculate frequency matrices
  file_by_cluster_freq <- cell_metadata %>% 
                          dplyr::select(file_name, cluster) %>% 
                          table() %>% data.frame()
  file_by_cluster_freq_norm  <- cell_metadata %>% 
                                dplyr::select(file_name, cluster) %>% 
                                table() %>% prop.table(margin = 1) %>% 
                                data.frame() %>% mutate(Freq = Freq * 100)
  ##############################
  
  
  
  ###### Calculate median expression matrices
  file_by_cluster_median_exp <- trans_exp %>% 
                                aggregate( by = list(file_name = cell_metadata$file_name, 
                                                     cluster = cell_metadata$cluster), 
                                           median 
                                           )
  cluster_median_exp <- trans_exp %>% 
                                aggregate( by = list(cluster = cell_metadata$cluster), median ) %>% column_to_rownames("cluster")

  file_median_exp <- trans_exp %>% 
    aggregate( by = list(file = cell_metadata$file_name), median ) %>% column_to_rownames("file") %>%
    arrange( match( row.names(.), file_metadata$file_name  ) )   # Arrange the rows (files) by the order of files in file_metadata
  
  ##############################
  CHECKPOINT <- 5
  print_message("Checkpoint #5 reached. Saving the newly generated data.")
  param_list <- prep_param_list()
  checkpoint_file = get_checkpoint_filename(out_dir, CHECKPOINT)
  save(CHECKPOINT, file_by_cluster_freq, file_by_cluster_freq_norm, file_by_cluster_median_exp, cluster_median_exp, file_median_exp, param_list, file=checkpoint_file )
  print_message(paste("Data is saved in", checkpoint_file) )
}





###### Prep files for Scaffold
if(CHECKPOINT == 5) {
  print_step_startup_msg()
  if( make_scaffold_map ) {  
    # Changing the column names of cluster_median_exp to the "desc" column from the parameter slot of the FCS objects, these column names must match the "desc" column in the FCS files from the "gated" populations.
    cluster_median_exp_scaff = cluster_median_exp %>% 
      setNames( marker_metadata[ match( colnames(cluster_median_exp), marker_metadata$marker_name ) , ]$desc )
    
    file_by_cluster_exp_freq <- cluster_median_exp_scaff %>% rownames_to_column("cluster") %>% 
                                full_join(file_by_cluster_freq, by = c("cluster" = "cluster") ) %>%
                                dplyr::rename(popsize = Freq, cellType = cluster) %>%   # Need to explicitely specify dplyr:: because I think there is another rename function in scaffold code(?) or some where else.
                                relocate(file_name,  .before = popsize)
    
    for( f in unique( as.vector(file_by_cluster_freq$file_name)) ) {
      exp_freq_f <- file_by_cluster_exp_freq %>% filter(file_name == f) %>% dplyr::rename( sample = file_name)
  
      cells_f <- cell_metadata$file_name == f
      cell_metadata_f <- cell_metadata %>% filter( cells_f )
      #trans_exp_f <- trans_exp %>% filter( cells_f ) %>% mutate( cellType = cell_metadata_f$cluster )
      
      trans_exp_f <- trans_exp %>% filter( cells_f ) %>% 
        setNames( marker_metadata[ match( colnames(trans_exp), marker_metadata$marker_name ) , ]$desc ) %>% 
        mutate( cellType = cell_metadata_f$cluster )
      
      write.table( exp_freq_f, paste0(scaffold_dir, "/", f, ".clustered.txt"), sep = "\t", quote = FALSE, row.names = F )
      scaffold_save( trans_exp_f, paste0(scaffold_dir, "/", f, ".clustered.all_events.RData") )
      file.create(paste0(scaffold_dir, "/", f))
      
    }
  }
  else {
    print_message("Gated populations are not provided, therefore, skipping the step of generating files for SCAFFoLD.")
  }
  CHECKPOINT <- 6
  print_message("Checkpoint #6 reached. SCAFFoLD input files are generated")
  param_list <- prep_param_list()
  checkpoint_file = get_checkpoint_filename(out_dir, CHECKPOINT)
  save(CHECKPOINT, param_list, file=checkpoint_file )
  print_message(paste("Data is saved in", checkpoint_file) )
}








###### Make SCAFFoLD map
if(CHECKPOINT == 6) {
  print_step_startup_msg()
  if( make_scaffold_map ) {
    
    #source all R files of scaffold
    #scaffold_Rcode_dir <- "/krummellab/data1/rpatel5/software/scaffold_tmp/statisticalScaffold-master/R/"
    #scaffold_Rfiles <- list.files(scaffold_Rcode_dir)
    #tmp_msgs = sapply( paste0(scaffold_Rcode_dir,scaffold_Rfiles), source )
    #library(scaffold)
    
    # Assumes that the "gated" directory exists in the scaffold_dir and that the first FCS file is used as a reference set.
    # Yet to add more validations to confirm the required inputs exist.
    ref_dataset_file <-
      paste0(file_by_cluster_freq$file_name[1], ".clustered.txt")
    scaffold:::run_analysis_gated(scaffold_dir, ref_dataset_file, marker_metadata[marker_metadata$used_for_scaffold,]$desc, 5)
    
    # Get the closest landmark populations and store them in cluster_metadata, and calculate various statistics and add them to .scaffold file.
    scaf_data <-
      scaffold_load(file.path(scaffold_dir, paste0(ref_dataset_file, ".scaffold")))
    
    
    # Get the closest landmark population and store in cluster_metadata
    one_graph = scaf_data$graphs[[1]]
    
    groups = vertex_attr(one_graph, "groups")
    high_scoring_edges = vertex_attr(one_graph, "highest_scoring_edge")
    high_scoring_edges_ends = ends(one_graph, high_scoring_edges[!is.na(groups)]) %>% 
      data.frame() %>% 
      setNames(c("Landmarks","Index")) %>% 
      mutate(Index = as.numeric(Index), cluster = groups[Index])
    cluster_metadata <-
      cluster_metadata %>% 
      mutate(landmarks = high_scoring_edges_ends[match(cluster, high_scoring_edges_ends$cluster),]$Landmarks)
      
    
    cov <-
      file_by_cluster_freq_norm %>% group_by(cluster) %>% summarise(val = var(Freq) / mean(Freq))
    var <-
      file_by_cluster_freq_norm %>% group_by(cluster) %>% summarise(val = var(Freq))
    mean <-
      file_by_cluster_freq_norm %>% group_by(cluster) %>% summarise(val = mean(Freq))
    median <-
      file_by_cluster_freq_norm %>% group_by(cluster) %>% summarise(val = median(Freq))
    mad <-
      file_by_cluster_freq_norm %>% group_by(cluster) %>% summarise(val = mad(Freq))
    norm_mad <-
      file_by_cluster_freq_norm %>% group_by(cluster) %>% summarise(val = mad(Freq)/median(Freq))
    
    stat_woOutliers <-
      function(d, stat = "cov") {
        # Removes Tukey's (boxplot) outliers (values beyond Q1 +/- 1.5 IQR) upto maximum of 10% of data points and returns CoV.
        outliers <- boxplot.stats(d)$out
        tenP_count <- round(length(d) * 0.1)
        if(length(outliers) > tenP_count)
          outliers <- tail(sort(outliers), tenP_count)
        d = d[!d %in% outliers]
        if (stat == "cov")
          return(var(d) / mean(d))
        if (stat == "var")
          return(var(d))
        if (stat == "mean")
          return(mean(d))
      }
    
    trimmed_cov <-
      file_by_cluster_freq_norm %>% group_by(cluster) %>% summarise(val = stat_woOutliers(Freq, "cov"))

    
    all_stats <- list(
      "CoV" = cov,
      "variance" = var,
      "mean_size" = mean,
      "median_size" = median,
      "MAD" = mad,
      "norm_MAD" = norm_mad,
      "trimmed_CoV" = trimmed_cov
    )
    
    
    # Saving the cluster-level statistics in cluster_metadata
    cluster_metadata <-
      cluster_metadata %>% 
      mutate(CoV = cov[match(cluster, cov$cluster),]$val) %>% 
      mutate(variance = var[match(cluster, var$cluster),]$val) %>% 
      mutate(mean_size = mean[match(cluster, mean$cluster),]$val) %>% 
      mutate(median_size = median[match(cluster, median$cluster),]$val) %>% 
      mutate(trimmed_CoV = trimmed_cov[match(cluster, trimmed_cov$cluster),]$val) %>%
      mutate(mad = mad[match(cluster, mad$cluster),]$val) %>%
      mutate(norm_mad = norm_mad[match(cluster, norm_mad$cluster),]$val)
      
    
    
        
    groups <- get.vertex.attribute( scaf_data$graphs[[1]], "groups" )
    index <- which( ! is.na( groups ) )
    
    # Preparing various statistics to be stored as vertex attributes.
    # Prepare values: add value 0 for landmark population and replace any NA/NaN values to 0
    all_stats_values <- list()
    for( n in names(all_stats) ) {
      all_stats_values[[n]] = c(
        rep(0, sum(is.na(groups))), 
        (all_stats[[n]] %>% slice(match(groups[index], cluster)) %>% 
           pull(val) %>% as.numeric())
      ) %>% ifelse(is.na(.), 0, .)
    }

    
    # Storing same statistic values as vertex attributes to all graphs.
    for( n in names(scaf_data$graphs) ) {
      for( s in names(all_stats_values) ) {
        values = all_stats_values[[s]]
        scaf_data$graphs[[n]] <-
          set.vertex.attribute(
            graph = scaf_data$graphs[[n]],
            name = s,
            value = values
          )
        scaf_data$dataset.statistics$max.marker.vals[[s]] <- max(values)
      }
    }

    scaffold_save(scaf_data, file.path(scaffold_dir, paste0(ref_dataset_file, ".scaffold")))
    
    setwd(out_dir)
  }
  else {
    print_message("Gated populations are not provided, therefore, skipping the SCAFFoLD map generation step.")
  }
  CHECKPOINT <- 7
  print_message("Checkpoint #7 reached.")
  param_list <- prep_param_list()
  checkpoint_file = get_checkpoint_filename(out_dir, CHECKPOINT)
  save(CHECKPOINT, cluster_metadata, param_list, file=checkpoint_file )
  print_message(paste("Data is saved in", checkpoint_file) )
}



# The following chunk must be executed at last to store the final processed data objects.
if(CHECKPOINT == 7) {
  print_step_startup_msg()
  CHECKPOINT <- 8
  print_message("Checkpoint #8 reached.")
  param_list <- prep_param_list()
  #checkpoint_file = file.path(out_dir,"processed_data.RData")
  checkpoint_file = get_checkpoint_filename(out_dir, CHECKPOINT)
  save(file_metadata, marker_metadata, cell_metadata, cluster_metadata, 
       file_by_cluster_freq, file_by_cluster_freq_norm, file_by_cluster_median_exp, 
       cluster_median_exp, file_median_exp, 
       param_list, file=checkpoint_file )
  print_message(paste("Final processed data objects are stored in:", checkpoint_file) )
}







library(ggExtra)
######### Prelim. Plotting commands

subsmpl = sample(1:nrow(trans_exp),10000)
myf = function(col) {
  exp = as.data.frame(trans_exp[subsmpl,col])[,1]
  as_tibble(cell_metadata[subsmpl,c("UMAP1","UMAP2")]) %>%
    ggplot(aes(x= UMAP1, y = UMAP2, col=exp)) +
    geom_point(size=1, alpha=0.3) +
    scale_color_viridis_c(option = "A", direction=-1) + 
    labs(x = "UMAP-1", 
         y = "UMAP-2", 
         title = col) + 
    coord_fixed() +
    theme_classic() +
    theme(plot.title = element_text(size = 12))
}
allPlots = lapply(as.vector(marker_metadata$marker_name[ marker_metadata$used_for_UMAP ]), myf)
png(file.path(out_dir, "feature_plots.png"), width=20, height=15, units = "in", res=600, pointsize = 4)
ggarrange(plotlist = allPlots, nrow = 6, ncol = 7)
dev.off()




## UMAPs split by clusters
set.seed(1234)
cell_metadata_sub = cell_metadata %>% sample_n(min(nrow(.),100000)) 
cell_metadata_sub$col_to_use = cell_metadata_sub$cluster
x_range = range(cell_metadata_sub$UMAP1) %>% scales::expand_range(add = 0.2)
y_range = range(cell_metadata_sub$UMAP2) %>% scales::expand_range(add = 0.2)
cm_tmp = cell_metadata_sub %>% group_by(col_to_use) %>% sample_n(33, replace=T)
plot_list = list()
groups = sort(unique(cell_metadata_sub$col_to_use))
for(i in groups ) {
  p = cell_metadata_sub %>% filter(col_to_use == i) %>% ggplot(aes(UMAP1, UMAP2)) + geom_rect(data=cm_tmp, aes(xmin=UMAP1-0.05,xmax=UMAP1+0.05,ymin=UMAP2-0.05,ymax=UMAP2+0.05), color="gray95", fill="gray95") + geom_point(size=0.1, alpha=0.3) + theme_classic() + xlim(x_range) + ylim(y_range) + ggtitle(i)
  plot_list[[i]] = ggMarginal(p, type = "histogram", fill = "red", binwidth = 0.1)
}
grid_size = get_plot_grid_layout(length(groups))
png(file.path(out_dir, "split_umap_by_cluster.png"), width=30, height=30, units = "in", res=600, pointsize = 4)
ggarrange(plotlist=plot_list, ncol=grid_size$ncol, nrow=grid_size$nrow)
dev.off()




batch_levels = file_metadata$pool_id %>% as.character() %>% unique()
batch_colors = kelly( length(batch_levels) + 1)[-1] %>% setNames(batch_levels)

pdf(file.path(out_dir, "plots.pdf"))
  # Downsampled UMAP plot with all samples
  set.seed(seed)
  down_size_for_plotting <- min( nrow( cell_metadata ), 100000)
  cell_metadata %>% sample_n(down_size_for_plotting) %>%
    ggplot( aes(x=UMAP1, y=UMAP2, color=as.factor(pool_id) ) ) + 
    geom_point( alpha=0.25, size=0.1 ) + 
    scale_color_manual(values = batch_colors) + 
    theme_classic() + 
    guides(color = guide_legend(override.aes = list(size = 3, alpha=1)))

  plotlist = list()  
  for( b in sort(batch_levels) ) {
    cell_metadata_sub = cell_metadata %>% filter( pool_id == b )
    down_size_for_plotting <- min( nrow( cell_metadata_sub ), 10000)
    plotlist[[ b ]] <-
      cell_metadata_sub %>% sample_n(down_size_for_plotting) %>%
        ggplot( aes(x=UMAP1, y=UMAP2 ) ) + 
        geom_point( alpha=0.25, size=0.1, color = batch_colors[b] ) + 
        theme_classic() + 
        ggtitle(b)
        guides(color = FALSE)
  }
  plot_grid_size = get_plot_grid_layout(length(batch_levels))
  ggarrange(plotlist = plotlist, nrow=plot_grid_size$nrow, ncol=plot_grid_size$ncol)
  
  
  if( any(file_metadata$control_sample) ) {
    # Downsampled UMAP plot with only control samples
    down_size_for_plotting <- min( nrow( cell_metadata[cell_metadata$control_sample,] ), 100000)
    p <-
      cell_metadata[cell_metadata$control_sample,] %>% sample_n(down_size_for_plotting) %>% 
      ggplot( aes(x=UMAP1, y=UMAP2, color=as.factor(pool_id) ) ) + 
      geom_point( alpha=0.25, size=0.1 ) + 
      scale_color_manual(values = batch_colors) + 
      theme_classic() + 
      guides(color = guide_legend(override.aes = list(size = 3, alpha=1))) +
      ggtitle("UMAP (Controls only)")
    print(p)
    
    # Downsampled UMAP plot with samples other than the control samples
    down_size_for_plotting <- min( nrow( cell_metadata[ ! cell_metadata$control_sample,] ), 100000)
    p <- 
      cell_metadata[ ! cell_metadata$control_sample,] %>% sample_n(down_size_for_plotting) %>% 
      ggplot( aes(x=UMAP1, y=UMAP2, color=as.factor(pool_id) ) ) + 
      geom_point( alpha=0.25, size=0.1 ) + 
      scale_color_manual(values = batch_colors) + 
      theme_classic() + 
      guides(color = guide_legend(override.aes = list(size = 3, alpha=1))) +
      ggtitle("UMAP (Non-controls only)")
    print(p)
  } 
  
  

  my_palette_reds = colorRampPalette(brewer.pal(n = 9, name = "Reds"))
  breaksList = seq(0, max(log1p(file_by_cluster_freq_norm$Freq)), by = 0.1)
  
  pheatmap( file_by_cluster_freq_norm %>% 
              pivot_wider(names_from="file_name",values_from=Freq) %>% 
              column_to_rownames("cluster") %>%
              log1p(),
            annotation_col =
              data.frame(
                batch = as.character(file_metadata$pool_id),
                control = as.character(file_metadata$control_sample),
                row.names = file_metadata$file_name
              ),
            annotation_colors = list(
              batch = batch_colors,
              control = c("TRUE" = "black", "FALSE" = "lightgrey")
            ),
            clustering_method = "ward.D2", show_colnames=F,
            breaks = breaksList, 
            color = my_palette_reds(length(breaksList)), 
            main = "File x Cluster: cell freq. (norm-log1p)",
            border_color = NA
            )


  pheatmap( file_by_cluster_freq %>% 
               mutate(pool_id = file_metadata[match( (.)$file_name, file_metadata$file_name),]$pool_id ) %>% 
               mutate(file_name = NULL) %>% 
               group_by(cluster, pool_id) %>% summarise(Freq = sum(Freq)) %>%
               pivot_wider(names_from=pool_id,values_from=Freq) %>% 
               column_to_rownames("cluster") %>% 
               as.matrix() %>% prop.table(2) %>% 
               data.frame(check.names = FALSE) %>% 
               mutate(. * 100) %>% log1p, 
            breaks = breaksList, 
            color = my_palette_reds(length(breaksList)), 
            main="Batch x Cluster: cell freq. (norm-log1p)",
            border_color = NA
            )

  if( any( file_metadata$control_sample ) ) {
    control_file_names <- file_metadata %>% filter(control_sample) %>% pull(file_name)
    
    
    pheatmap( file_by_cluster_freq %>% filter( file_name %in% control_file_names ) %>% 
                mutate(pool_id = file_metadata[match( (.)$file_name, file_metadata$file_name),]$pool_id ) %>% 
                mutate(file_name = NULL) %>% 
                group_by(cluster, pool_id) %>% summarise(Freq = sum(Freq)) %>%
                pivot_wider(names_from=pool_id,values_from=Freq) %>% 
                column_to_rownames("cluster") %>% 
                as.matrix() %>% prop.table(2) %>% 
                data.frame(check.names = FALSE) %>% 
                mutate(. * 100) %>% log1p, 
              breaks = breaksList, 
              color = my_palette_reds(length(breaksList)), 
              main="Batch x Cluster: cell freq. (norm-log1p) (Controls only)",
              border_color = NA
    )
    
    pheatmap( file_by_cluster_freq %>% filter( ! file_name %in% control_file_names ) %>% 
                mutate(pool_id = file_metadata[match( (.)$file_name, file_metadata$file_name),]$pool_id ) %>% 
                mutate(file_name = NULL) %>% 
                group_by(cluster, pool_id) %>% summarise(Freq = sum(Freq)) %>%
                pivot_wider(names_from=pool_id,values_from=Freq) %>% 
                column_to_rownames("cluster") %>% 
                as.matrix() %>% prop.table(2) %>% 
                data.frame(check.names = FALSE) %>% 
                mutate(. * 100) %>% log1p, 
              breaks = breaksList, 
              color = my_palette_reds(length(breaksList)), 
              main="Batch x Cluster: cell freq. (norm-log1p) (Non-controls only)",
              border_color = NA
    )
  }
  
  
  
  my_palette_greens = colorRampPalette(brewer.pal(n = 9, name = "Greens"))
  breaksList = seq(0, max(log1p(file_median_exp)), by = 0.1)
  
  pheatmap(
            file_median_exp %>% dplyr::select(which(marker_metadata$used_for_clustering)) %>% t() %>% log1p(),
            annotation_col =
              data.frame(
                batch = as.character(file_metadata$pool_id),
                control = as.character(file_metadata$control_sample),
                row.names = file_metadata$file_name
              ),
            annotation_colors = list(
              batch = batch_colors,
              control = c("TRUE" = "black", "FALSE" = "lightgrey")
            ),
            clustering_method = "ward.D2", show_colnames=F,
            breaks = breaksList, 
            color = my_palette_greens(length(breaksList)), 
            main = "File x Cluster: expression (arcsinh-log1p)",
            border_color = NA
  )
  
  

  
  dev.off()
  
  
  
  
  


  # DON'T RUN
  if(FALSE) {

    
    file_median_exp_z = apply(file_median_exp, 2, scale) %>% t() %>% data.frame() %>% 
      setNames(row.names(file_median_exp)) %>% na.omit()
    
    pal_set3 = colorRampPalette(brewer.pal(12,"Set3"))
    pal_set1 = colorRampPalette(brewer.pal(9,"Set1"))
    pool_levels = file_metadata$pool_id %>% unique()
    donor_levels = file_metadata$donor_id %>% unique()
    ha = HeatmapAnnotation(
      pool = file_metadata$pool_id, 
      donor = file_metadata$donor_id, 
      col = list(pool = pal_set1( length(pool_levels) ) %>% setNames( pool_levels ),
                 donor = pal_set3( length(donor_levels) ) %>% setNames( donor_levels )
      ),
      na_col = "black"
    )
    
    
    Heatmap(file_median_exp_z, name = "", top_annotation = ha)
    
    

gex_sub = inner_join(rownames_to_column(cell_metadata), rownames_to_column(trans_exp), by=c("rowname" = "rowname")) %>% 
  sample_n( 10000 )

myf = function(col, gex_sub) {
  print(col)
  # Mute the extremely high values.
  # Get the 0.1%th value and replace all values higher than that by that value and store the values to a new column called "mod_val"
  limit_val = gex_sub %>% pull(get(col)) %>% quantile(0.999)
  gex_sub = gex_sub %>% mutate( mod_val = ifelse( get(col) > limit_val, limit_val, get(col) )  )
  
  gex_sub %>% 
    ggplot( aes(x=UMAP1, y=UMAP2, color= mod_val ) ) + 
    geom_point( alpha=0.25, size=0.1 ) + 
    scale_color_viridis_c(option="inferno", direction = -1) +
    labs(x = "UMAP-1", 
         y = "UMAP-2", 
         title = col) + 
    coord_fixed() +
    theme_classic() +
    theme(plot.title = element_text(size = 12))
}
allPlots = lapply(marker_metadata$marker_name, myf, gex_sub = gex_sub)
png(file.path(out_dir, "feature_plots.png"), width=20, height=15, units = "in", res=600, pointsize = 4)
ggarrange(plotlist = allPlots, nrow = 6, ncol = 7)
dev.off()


V(t$graphs[[1]])$tt = 1:62
t$dataset.statistics$max.marker.vals[["tt"]] = 62



}











