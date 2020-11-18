if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("CEMiTool")
library(CEMiTool)
library(ggplot2)

# Import GMT file
gmt_in <- read_gmt("FILE_LOCATION/Human_geneset.gmt")

# run cemitool
#Apply vst = FALSE ------ Microarray
#Apply vst = TRUE ------ RNASeq


# Filter FALSE
cem <- cemitool(expr0, sample_annot, gmt_in, int_df, 
                  filter=FALSE, plot=TRUE, verbose=TRUE, apply_vst=TRUE)

#Filter TRUE
cem <- cemitool(expr0, sample_annot, gmt_in, int_df, 
                filter=TRUE, plot=TRUE, verbose=TRUE, apply_vst=TRUE, min_ngen=15, n_genes=4000)

# create report as html document
generate_report(cem, directory="./Report")

# write analysis results into files
write_files(cem, directory="./Tables")
force=TRUE

# save all plots
save_plots(cem, "all", directory="./Plots")
