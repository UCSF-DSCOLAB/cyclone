# Running the CyTOF Pipeline

<p class="author-name">Daniel Bunis<span class="affil-mark">1†</span>, Rebecca Jaszczak<span class="affil-mark">1‡</span> and Ravi Patel<span class="affil-mark">1§</span></p>
<p class="author-affiliation"><span class="affil-mark">1</span>UCSF DSCoLab</p>
<p class="author-email"><span class="affil-mark">†</span><a href="mailto:daniel.bunis@ucsf.edu">daniel.bunis@ucsf.edu</a><br><span class="affil-mark">‡</span><a href="mailto:rebecca.jaszczak@ucsf.edu">rebecca.jaszczak@ucsf.edu</a><br><span class="affil-mark">§</span><a href="mailto:ravi.patel2@ucsf.edu">ravi.patel2@ucsf.edu</a></p>

## Contents

-   [Introduction](#introduction)
    -   [Starting point](#starting-point)
    -   [What this pipeline does do](#what-this-pipeline-does-do)
    -   [What this pipeline does NOT
        do](#what-this-pipeline-does-not-do)
-   [Setup](#setup)
    -   [Installing FlosSOM and other
        dependencies](#installing-flossom-and-other-dependencies)
    -   [Choosing an output directory](#choosing-an-output-directory)
    -   [File Metadata](#file-metadata)
    -   [Marker Metadata](#marker-metadata)
    -   [Config.yml File](#config.yml-file)
    -   [Notes on Parallelization](#notes-on-parallelization)
-   [Running the Pipeline](#running-the-pipeline)
-   [Checkpoints](#checkpoints)
    -   [Output Checkpoints 1-8](#output-checkpoints-1-8)
        -   [Plots output:](#plots-output)

## Introduction

This pipeline was developed collaboratively by a small team of wet lab
“bench” scientists and data scientists of [the
CoLabs](https://colabs.ucsf.edu/) and associated labs at UCSF. We aimed
to develop a pipeline for CyTOF analysis which scaled well for datasets
with large numbers of events and provided consistency of clusters when
down-sampling. We also want the pipeline to be easily implemented by
“bench” scientists to analyze their own CyTOF data, bypassing the need
for programming expertise. Cyclone is a stand-alone tool for the
analysis of full single-cell high dimensional CyTOF data. It can
integerate with SCAFFoLD workflow, and perform custom statistics and
visualizations.

### Starting point

To get started with this pipeline, you should have:

-   A folder containing fcs files from cytof data
-   file metadata (described below)
-   marker metadata (described below)
-   csv of grid sizes (`CyTOF_pipeline/grid_sizes.csv`)

### What this pipeline does do

-   Chunks work of the steps below across 8 ‘checkpoints’
-   Calculates UMAP embeddings
-   Performs clustering with FlowSOM (first, over an optimization space
    to allow users to pick their ideal parameters)
-   Then, with user-chosen parameters, calculates and outputs various
    metrics (detailed in the Checkpoints sections below)

### What this pipeline does NOT do

The pipeline assumes you have already performed any of the below
pre-processing steps which may be necessary for your data:

-   subset populations
-   debarcoding
-   bead normalization
-   batch correction
-   compensation
-   special transformations (beyond the arcsinh transformation which the
    pipeline does)

Optionally, if your data has been batch corrected by a method that does
not output corrected FCS files, and only outputs a corrected matrix, it
is still possible to use the pipeline. Feel free to contact us if that
is your use-case.

## Setup

To run the pipeline, you first need to:

-   Clone the repo
-   Install all required packages
-   Choose or create an intended output directory
-   Create a set of metadata files:
    -   `file_metadata.csv`
    -   `marker_metadata.csv`
-   Finalize your `config.yml` file
-   Submit your pipeline call as an `Rscript`.

### Installing FlosSOM and other dependencies

Simply installing this <package-name> counterpart will take care of most
required installations. You can install from github with:

    if (!requireNamespace("remotes")) {
        install.packages("remotes")
    }
    remotes::install_github("ravipatel4/CyTOF_pipeline", subdir = "cyclone")

Some packages are not absolutely required, so we leave to the user
whether to install them. If you want them, you will need to install them
yourself:

-   Scaffold

<!-- -->

    remotes::install_github("nolanlab/scaffold")

### Choosing an output directory

This is the directory where all your checkpoints will be saved to, as
well as any plots that result from running the pipeline.

### File Metadata

Required columns:

-   `file_name`: the name (not including full path) of the FCS files
-   `donor_id`: specifies the sample origin
-   `pool_id`: use any ID mechanism to indicate samples from the same
    CyTOF pool
-   `control_sample`: ‘control’ means the batch control used in pools
    (NOT experimental control)

In the event that your data doesn’t have pools, or different donors, you
can simply populate these columns with the same value, e.g. `1`.

How it’s used: To direct the pipeline to the location of the FCS files,
and associate some metadata with them (donor, pool).

How to create: You can make this in R/python/Excel.

Recommended format: `CSV (Comma delimited)`

### Marker Metadata

Required columns:

-   `channel_name`: from the FCS files
-   `marker_name`: from the FCS files
-   `used_for_UMAP`: Boolean, whether the marker should be used in UMAP
    calculations
-   `used_for_clustering`: Boolean, whether the marker should be used in
    clustering calculations
-   `used_for_scaffold`: Boolean, whether the marker should be used in
    scaffold analysis

In most cases, `used_for_UMAP`, `used_for_clustering`,
`used_for_scaffold` can be the same, but there may be cases where you
might use a marker for clustering, but not for UMAP calculation, for
example.

How to create: Use `CyTOF_pipeline/make_marker_metadata_csv.R` script.

To invoke this script, it is run from the terminal/command line prompt, NOT
in Rstudio. 

`$ Rscript make_marker_metadata_csv.R -f fcs_file.fcs -o marker_metadata.csv`

Above is an example invocation. You must update your path to wherever your
copy of `make_marker_metadata_csv.R` is located, and replace `fcs_file.fcs`
with the full path and name of one of your own FCS files to run. 
`marker_metadata.csv` is the suggested name output.

### Config.yml File

The file directs the pipeline to find all inputs, and is the main point
of control for how the pipeline will be run.

Contents

-   path of directory where FCS files are; will be searched recursively
-   path for where you want the output files saved
-   where `file_metadata.csv` is located
-   where `marker_metadata.csv` is located
-   path for `gated_fcs_dir`, the template for scaffold analysis (can be
    left blank if not running scaffold)
-   `rerun_from_step`: if you need to restart the pipeline at an
    intermediate step
-   data processing related: can vary per user needs or leave as default
-   UMAP related: can vary per user needs or leave as default
-   clustering: `flowsom` recommended, as it is faster, although `clara`
    is an alternative
-   nthreads: the option to parallelize the optimization step (see Notes
    on Parallelization; not available on Windows)

Template: `CyTOF_pipeline/config.yml`

How to create: Fill in the `config.yml` file, saving it with whatever
name you would like. We will use `my_config.yml`.

### Notes on Parallelization

Parallelization is available on Mac and Linux operating systems. It is
not available for Windows users.

If you are submitting on a shared compute cluster, here is a suggested
workflow:

-   Structure your initial job call to have one more thread than you
    will assign in `my_config.yml`.
    -   E.g., if you intend to parallel process over 4 cores, request 5
        cores from your scheduler.
    -   In the same example, ensure that `nthreads: 4` value is assigned
        in `my_config.yml`
-   Larger datasets (over 20 million cells) may need to allocate scratch
    directories, e.g. `--gres=scratch:500G`, pending cluster
    requirements.
-   Request enough RAM for the number of cells contained in your
    dataset.
    -   As structured currently, R will split the total amount of RAM
        passed in your job to each of the processes requested in
        `my_config.yml`
    -   We have used ~100Gb RAM without parallaization with ~20million
        cells (config 'nthreads: 1'), ~420Gb RAM for 4x parallelization
        with ~40million cells (config 'nthreads: 4', scheduler ntasks=5),
        and ~520Gb RAM for 3x parallelization with ~50 million cells
        (config 'nthreads: 3', scheduler ntasks=4). This may vary based on
        your computer/cluster specs.

If processing in parallel, please note that there will be no log
messages `FlowSOM clustering begins for grid...` after you have received
the `Starting the FlowSOM clustering...` message.

## Running the Pipeline

Once everything is set up,run the pipeline, passing as a command line
variable your modified `my_config.yml`.

    Rscript cytof_pipeline.R -c my_config.yml

## Checkpoints

The pipeline is broken down into multiple chunks, where each chunk saves
its output as a checkpoint.Rdata file.

### Output Checkpoints 1-8

-   Read, preprocess and transform the data. Output:
    `checkpoint_1.RData`
-   Calculate UMAP. Output: `checkpoint_2.RData`
-   Optimize clustering parameters. Output: `heckpoint_3.RData` &
    `clustering_param_optimization.pdf`
-   Perform clustering. Output: `checkpoint_4.RData`
-   Generate various matrices. Output: `checkpoint_5.RData`
-   Prepare input files for SCAFFoLD analysis. Output:
    `checkpoint_6.RData`
-   Generate SCAFFoLD map. Output: `checkpoint_7.RData`
-   Produce final set of visualizations. Output: `checkpoint_8.RData`
    and plots below

#### Plots output:

-   `clustering_param_optimization.pdf` - DBI x cluster counts graph
-   `feature_plots.png` - a set list of markers helpful for cell type
    identification
-   `split_umap_by_cluster.png` - histogram of each clusters density in
    UMAP space dimensions
-   `Rplots.pdf` per cluster, histogram of median expression of all
    markers
-   `plots.pdf` - umap colored by pool ID to survey potential batch
    effects
