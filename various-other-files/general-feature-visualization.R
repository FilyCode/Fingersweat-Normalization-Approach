source("X:/Studenten_Schueler/Philipp_Trollmann/R Project/MS-data-analysis_functions.R")
source("X:/Studenten_Schueler/Philipp_Trollmann/R Project/feature-visualization-script.R")


# set all the experiment and file names
wd = "X:/Studenten_Schueler/Philipp_Trollmann/Experiments/"
#wd = "X:/Studenten_Schueler/Ulrike Jiras/Experiments/"
setwd(wd)

#exp = "Weinstudie_Chrom-Kim"
#exp = "FiS-45min-Abstandsmessung"
exp = "MTX1-8"
datadir = paste0(wd, exp, "/data/")
resultsdir = paste0(wd, exp, "/results/")
figuredir  = paste0(wd, "Figures/", exp)

# a name that is unique for each sample type, to distinguish between different sample types (ex. "Chrom", "Kim)
distinguishing_sample_names <- c("Chrom", "Kim")
#distinguishing_sample_names <- c("231107", "231213", "240103", "240208", "240214", "240314", "240418")

# choose 2 names of your distinguishing_sample_names to compare in scatter plots
#sample_names_for_comparison <- c("240314", "240418")
sample_names_for_comparison <- c("Chrom", "Kim")



### Calculations and result plots (dont need to change anything here)

# get all significant and abundant features
significant_abundant_features_object <- get_significant_abundant_features(datadir, resultsdir = resultsdir, figuredir = figuredir, exp = exp)

# plot the mean of every sample in various plots to pdf
get_all_mean_plots(significant_abundant_features_object$significant_abundant_features,
               distinguish_sample_names = distinguishing_sample_names, resultsdir = resultsdir, figuredir = figuredir)

# plots each sample separatly in various plots to pdf
plot_feature_ranking_each_sample_to_pdf(significant_abundant_features_object$significant_abundant_features, resultsdir, figuredir)

# plots scatter plots for comparison of two sample types and gives back featrues that differ strongly
varying_features <- generate_scatter_plots_and_get_filterd_featrues(significant_abundant_features_object$significant_abundant_features, 
                                                                    sample_names = sample_names_for_comparison, resultsdir)
# write features that vary strongly to a csv file  
write.csv(varying_features$mean_comparison, file = paste0(resultsdir,
                                               paste0(sample_names_for_comparison[1],"vs", sample_names_for_comparison[2],
                                                                 "_varying-mean-Features-", exp, ".csv")))
write.csv(varying_features$rank_comparison, file = paste0(resultsdir, 
                                               paste0(sample_names_for_comparison[1],"vs", sample_names_for_comparison[2],
                                                                  "_varying-rank-Features-", exp, ".csv")))

