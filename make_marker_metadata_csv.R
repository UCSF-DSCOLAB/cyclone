########################### Loading packages ###########################
suppressMessages({
  library(flowCore) # FCS I/O
  library(tidyverse)
  library(optparse)
  library(stringr)
})
########################### END: Loading packages ######################



################ Define & parse command-line options ###################
# Define command-line options
option_list <- list(
  make_option(c("-f", "--fcs_file"), type="character", default=NULL,
              help="Path of one of the FCS files that will be used by cytof.R.", metavar="character"),
  make_option(c("-o", "--output_csv_file"), type="character", default=NULL,
              help="Output file name. The marker metadata will be stored in this file in CSV format.", metavar="character")
)

# Parse command-line options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
fcs_file <- opt$fcs_file
output_csv_file <- opt$output_csv_file
############## END: Define & parse command-line options ################

if( is.null(fcs_file) & is.null(output_csv_file) ) {
  cat("No input is found. Use -h/--help option to read the usage.\n")
  quit(save = "no", status=1)
}

if( is.null(fcs_file) | is.null(output_csv_file) ) {
  cat("Please provide both, the FCS file and the output file name.\n")
  quit(save = "no", status=1)
}

if(file.exists(fcs_file)) {
  fcs_data <-
    read.FCS(fcs_file) # Read FCS file and store in a temporary variable
  fcs_param_data <- fcs_data@parameters@data  # Extract parameter data
  marker_metadata <- fcs_param_data %>%
    mutate(desc = ifelse(is.na(desc), name, desc)) %>%   # Convert NA values to character "NA"
    mutate(marker_name_short = gsub("^[^_]+_", "", desc) ) %>%
    select(name, desc, marker_name_short) %>% setNames(c("channel_name", "marker_name", "marker_name_short"))  # Select channel name and marker name and rename the column names
  # Make the marker_name_short column unique, but adding _1, _2, etc to the duplicated values (e.g. DEAD_1, DEAD_2, etc)
  duplicated_levels <- marker_metadata$marker_name_short[ duplicated(marker_metadata$marker_name_short) ] %>% unique()
  for(l in duplicated_levels) {
    dup_ind <- marker_metadata$marker_name_short == l
    marker_metadata$marker_name_short[ dup_ind ] = paste(l,1:sum(dup_ind), sep="_")
  }
  
  marker_metadata$used_for_UMAP <- FALSE
  marker_metadata$used_for_clustering <- FALSE
  marker_metadata$used_for_scaffold <- FALSE
  write.csv( marker_metadata, file = output_csv_file )  # Write out the output file
} else {
  cat("Cannot read ", fcs_file, ".\n", sep = "")
}
