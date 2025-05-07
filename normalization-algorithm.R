source("//prot-data/MS Storage 3/tmp_PT/R Project/MS-data-analysis_functions.R")
source("//prot-data/MS Storage 3/tmp_PT/R Project/feature-visualization-script.R")

# set all the experiment and file names
#wd = "X:/Studenten_Schueler/Philipp_Trollmann/Experiments/"
wd = "//prot-data/MS Storage 3/tmp_PT/Experiments/"
setwd(wd)

#exp = "Julia_Coffein"
exp = "Weinstudie-neu"

datadir = paste0(wd, exp, "/data/")
resultsdir = paste0(wd, exp, "/results/")

# When finished files exist use this
  #targeted_experiment_data <- read.csv(paste0(resultsdir, 'all_targeted_experiment_data_filtered.csv'), row.names = 1)
  #norm_feature_list <- read.csv(paste0(resultsdir, 'all_norm_feature_list_filtered.csv'), row.names = 1)
  #pqn_vsn_list <- read.csv(paste0(resultsdir, 'all_pqn_vsn_list_filtered.csv'), row.names = 1)
  #full_prep_data <- read.csv(paste0(resultsdir, 'all_full_prep_data_filtered.csv'), row.names = 1)

# Otherwise create needed dataframes here
tlists = paste0(wd, "Transition-Lists/")
tfile <- "WinestudyTransitionList_curves_adjusted.csv"
clustered_features <- 'WinestudyTransitionList_clustered-molecules.csv'
# normalization_list  = paste0(wd, "normalization/biological_normalization_molecule_list_for-combination-test-for-wine-study_bigger3.csv")
# normalization_list  = paste0(wd, "normalization/biological_normalization_molecule_list_for-combination-test-for-wine-study_random-test.csv")
# normalization_list  = paste0(wd, "normalization/biological_normalization_molecule_list_for-combination-test-for-wine-study_filtered-good-sweat-features.csv")
# normalization_list  = paste0(wd, "normalization/filtered_good_sweat-features_for-norm.csv")
normalization_list  = paste0(wd, "normalization/filtered_selected-features_with-filtered-good_sweat-features_for-norm.csv")



info_file <- "/Wine_names.csv"
add_info_file <- paste0(wd, exp, info_file)
add_info <- read.csv(add_info_file)
if (length(add_info) < 2) {
  add_info <- read.csv2(add_info_file)
}

exp_data <- process_merged_files(datadir)
all_data <- exp_data$data
df_ids <- exp_data$df_ids

transitionFile <- paste0(tlists, tfile)
targeted_experiment_data <- extract_feature_list(all_data, df_ids, transitionFile)
# Standardize the rownames and colnames to ensure compatibility
normalize_names <- function(names) {
  # Remove 6-digit numbers followed by an underscore (if present)
  gsub("^[0-9]{6}_", "", names)
}
colnames(targeted_experiment_data) <- normalize_names(colnames(targeted_experiment_data))

transitionFile_cluster <- paste0(tlists, clustered_features)
clustered_data <- extract_feature_list(all_data, df_ids, transitionFile_cluster)
# Standardize the rownames and colnames to ensure compatibility
colnames(clustered_data) <- normalize_names(colnames(clustered_data))


info_file_dir <- paste0(wd, exp, info_file)
significant_abundant_features_object <- get_significant_abundant_features(datadir, exp_data, info_file_dir, resultsdir, figuredir, exp, tlists)
all_needed_features <- significant_abundant_features_object$significant_abundant_features

norm_feature_names <- read.csv(normalization_list)
if (ncol(norm_feature_names) < 2) {
  norm_feature_names <- read.csv2(normalization_list)
}
norm_feature_list <- extract_feature_list(all_data, df_ids, normalization_list)
#norm_feature_list <- extract_feature_list(all_data, df_ids, normalization_list, ppm_range = 1, rt_range = 0.01)
norm_feature_list <- norm_feature_list %>% distinct(id, .keep_all = TRUE) # remove duplicates
colnames(norm_feature_list) <- normalize_names(colnames(norm_feature_list))


# Remove outlier from dataframes
targeted_experiment_data <- targeted_experiment_data[, c(colnames(targeted_experiment_data[c(1:4)]), intersect(colnames(targeted_experiment_data[-c(1:4, ncol(targeted_experiment_data))]), 
                                                                        colnames(all_needed_features)), colnames(targeted_experiment_data[ncol(targeted_experiment_data)]))]
clustered_data <- clustered_data[, c(colnames(clustered_data[c(1:4)]), intersect(colnames(clustered_data[-c(1:4, ncol(clustered_data))]), 
                                                                                                               colnames(all_needed_features)), colnames(clustered_data[ncol(clustered_data)]))]
norm_feature_list <- norm_feature_list[, c(colnames(norm_feature_list[c(1:4)]), intersect(colnames(norm_feature_list[-c(1:4, ncol(norm_feature_list))]), 
                                                                                                        colnames(all_needed_features)), colnames(norm_feature_list[ncol(norm_feature_list)]))]
full_prep_data <- prepare_data_for_plot(targeted_experiment_data, add_info_file)
full_prep_data <- full_prep_data[!is.na(full_prep_data$Donor),]

plot_molecule_curves(resultsdir, full_prep_data, Group = 'Donor')


full_prep_data <- full_prep_data %>% # Adjust time and process the data
  # Remove donors with negControl as TRUE
  filter(!negControl) %>%
  # Adjust time based on Molecule.Name
  mutate(
    time = case_when(
      Molecule.Name %in% c("Tryptanthrin", 
                           "3,7-Dihydro-1-butyl-7-(5,6-dihydroxyhexyl)-3-methyl-1H-purine-2,6-dione") ~ time - 120,
      TRUE ~ time - 180
    )
  ) %>%
  # Remove rows with negative times
  filter(time >= 0) %>%
  # Set Area to 0 for time = 0
  mutate(Area = ifelse(time == 0, 300, Area)) 


# Choose 15 random features from the ones with known formula
  # set.seed(456)  # Set seed for reproducibility
  # random_rows <- all_needed_features[!is.na(all_needed_features$Formula), ][sample(nrow(all_needed_features[!is.na(all_needed_features$Formula), ]), min(15, nrow(all_needed_features[!is.na(all_needed_features$Formula), ]))), ]

# create pure random normalization list
  # Define the log-scale range
    min_exp <- 3   # log10(1e3)
    max_exp <- 8   # log10(1e8)
    set.seed(456)
  # Get the columns to update
    cols_to_update <- 5:(ncol(norm_feature_list)-1)
  # Fill with log-uniform values
    norm_feature_list[ , cols_to_update] <- replicate(
                                              length(cols_to_update),
                                              10^runif(n = nrow(norm_feature_list), min = min_exp, max = max_exp))


# Call optimization algorithm for normalization combination
  results <- calculate_normalization_combinations_check_fits_for_biological_kinetics_with_optimization_algorithm(
    all_needed_features, all_data, targeted_experiment_data, norm_feature_list, full_prep_data, add_info, clustered_data, resultsdir, Group = 'Donor', df_ids, transitionFile, transitionFile_cluster, normalization_list)

# Call bruteforce algorithm for normalization combination
  results <- calculate_normalization_combinations_check_fits_for_biological_kinetics_with_optimization_algorithm(
    all_needed_features, all_data, targeted_experiment_data, norm_feature_list, full_prep_data, add_info, clustered_data, paste0(resultsdir,"2compartment-bateman-normalization-bruteforce-combinatoric_selected-with-sweat-features/"), Group = 'Donor', df_ids, transitionFile, transitionFile_cluster, normalization_list, bayesoptim = FALSE)
  

# Analyse results from file (from single run)
  #results <- as.data.table(read.csv(paste0(resultsdir, 'bayesOpt_normalising_scoreSummaryTable.csv'), row.names = 1))
  resultsdir <- paste0(wd, exp, "/results/",'2compartment-bateman-normalization-bruteforce-combinatoric_sweat-linear-features/')
  results1 <- as.data.table(fread(file=paste0(resultsdir, 'bayesOpt_normalising_scoreSummaryTable_1to6-from-selected_all_combinations.csv'))) %>% filter(!is.na(significant_features_norm)) %>% select(-V1)
  results1$group <- 'selected features'
  results2 <- as.data.table(fread(file=paste0(resultsdir, 'bayesOpt_normalising_scoreSummaryTable_1to6-from-random-features_all_combinations.csv'))) %>% filter(!is.na(significant_features_norm))
  results2$group <- 'random features'
  results3 <-  as.data.table(fread(file=paste0(resultsdir, 'bayesOpt_normalising_scoreSummaryTable_1to6-from-linear-sweat_all_combinations.csv'))) %>% filter(!is.na(significant_features_norm))
  results3$group <- 'linear sweat features'
  results4 <-  as.data.table(fread(file=paste0(resultsdir, 'bayesOpt_normalising_scoreSummaryTable_1to11-from-contaminations_all_combinations.csv'))) %>% filter(!is.na(significant_features_norm))
  results4$group <- 'linear contamination features'
  results5 <-  as.data.table(fread(file=paste0(resultsdir, 'bayesOpt_normalising_scoreSummaryTable_1to6-from-random-values_all_combinations.csv'))) %>% filter(!is.na(significant_features_norm))
  results5$group <- 'pure random values'
  results <- rbind(results1, results2, results3, results4, results5)
  fwrite(results, paste0(resultsdir, 'scoreSummaryTable_all_combinations-1to6combinatoric-from_selected_randomFeat_linSweat_linCont_pureRand.csv'))
  all_results <- results
  
  resultsdir <- paste0(wd, exp, "/results/",'2compartment-bateman-normalization-bruteforce-combinatoric_selected-with-sweat-features/1to8/')
  results1 <- as.data.table(fread(file=paste0(resultsdir, 'bayesOpt_normalising_scoreSummaryTable_all_combinations_1to7-combinatorial.csv'))) %>% filter(!is.na(significant_features_norm))
  results2 <- as.data.table(fread(file=paste0(resultsdir, 'bayesOpt_normalising_scoreSummaryTable_all_combinations_8-combinatorial.csv'))) %>% filter(!is.na(significant_features_norm))
  results <- rbind(results1, results2)
  fwrite(results, paste0(resultsdir, 'scoreSummaryTable_all_combinations-1to8combinatoric-from_best-selected-and-linSweat.csv'))
  all_results <- results
  
  resultsdir <- paste0(wd, exp, "/results/",'2compartment-bateman-normalization-bruteforce-combinatoric_selected-with-sweat-features/')
  results <- as.data.table(fread(file=paste0(resultsdir, 'bayesOpt_normalising_scoreSummaryTable_all_combinations_1to7-combinatorial.csv'))) %>% filter(!is.na(significant_features_norm))
  all_results <- results
  
  resultsdir <- paste0(wd, exp, "/results/",'2compartment-bateman-normalization-bruteforce-combinatoric-random-reference-test/')
  results <- as.data.table(read.csv(paste0(resultsdir, 'bayesOpt_normalising_scoreSummaryTable_all_combinations.csv'))) %>% filter(!is.na(significant_features_norm))
  all_results <- results
  
  resultsdir <- paste0(wd, exp, "/results/",'2compartment-bateman-normalization-bruteforce-combinatoric_sweat-linear-features/')
  results <- as.data.table(read.csv(paste0(resultsdir, 'bayesOpt_normalising_scoreSummaryTable_all_combinations.csv'))) %>% filter(!is.na(significant_features_norm))
  all_results <- results
  
  resultsdir <- paste0(wd, exp, "/results/",'2compartment-bateman-normalization-bruteforce-combinatoric_selected-with-sweat-features/')
  results <- as.data.table(read.csv(paste0(resultsdir, 'bayesOpt_normalising_scoreSummaryTable_all_combinations_1to7-combinatorial.csv'))) %>% filter(!is.na(significant_features_norm))
  all_results <- results
  
  
  print_norm_molecule_network_to_pdf(results, cutoff = NA, remove_duplicates = FALSE)
  print_norm_molecule_network_to_pdf(results, cutoff = NA, remove_duplicates = TRUE)
  print_norm_molecule_network_to_pdf(results, cutoff = 2, remove_duplicates = FALSE)
  print_norm_molecule_network_to_pdf(results, cutoff = 2, remove_duplicates = TRUE)
  print_norm_molecule_network_to_pdf(results, cutoff = 1.5, remove_duplicates = FALSE)
  print_norm_molecule_network_to_pdf(results, cutoff = 1.5, remove_duplicates = TRUE)
  print_norm_molecule_network_to_pdf(results, cutoff = 1.3, remove_duplicates = FALSE)
  print_norm_molecule_network_to_pdf(results, cutoff = 1.3, remove_duplicates = TRUE)
  print_norm_molecule_network_to_pdf(results, cutoff = 1.2, remove_duplicates = FALSE)
  print_norm_molecule_network_to_pdf(results, cutoff = 1.2, remove_duplicates = TRUE)
  print_norm_molecule_network_to_pdf(results, cutoff = 1.1, remove_duplicates = FALSE)
  print_norm_molecule_network_to_pdf(results, cutoff = 1.1, remove_duplicates = TRUE)


# Run the optimization several times as different local minima can be found (let run 10 times for max. of 14 epochs or 24h)
  resultsdir <- paste0(resultsdir,'2compartment-bateman-normalization-even-further-advanced-14-epochs-10runs/')
  for (run in c(8:10)) {
    results <- calculate_normalization_combinations_check_fits_for_biological_kinetics_with_optimization_algorithm(
      all_needed_features, all_data, targeted_experiment_data, norm_feature_list, full_prep_data, add_info, clustered_data, resultsdir, Group = 'Donor', df_ids, transitionFile, transitionFile_cluster, normalization_list, run = run)
    
    results$run <- run
    
    print_norm_molecule_network_to_pdf(results, cutoff = NA, remove_duplicates = FALSE, run = run)
    print_norm_molecule_network_to_pdf(results, cutoff = NA, remove_duplicates = TRUE, run = run)
    print_norm_molecule_network_to_pdf(results, cutoff = 2, remove_duplicates = FALSE, run = run)
    print_norm_molecule_network_to_pdf(results, cutoff = 2, remove_duplicates = TRUE, run = run)
    print_norm_molecule_network_to_pdf(results, cutoff = 1.5, remove_duplicates = FALSE, run = run)
    print_norm_molecule_network_to_pdf(results, cutoff = 1.5, remove_duplicates = TRUE, run = run)
    print_norm_molecule_network_to_pdf(results, cutoff = 1.3, remove_duplicates = FALSE, run = run)
    print_norm_molecule_network_to_pdf(results, cutoff = 1.3, remove_duplicates = TRUE, run = run)
    print_norm_molecule_network_to_pdf(results, cutoff = 1.2, remove_duplicates = FALSE, run = run)
    print_norm_molecule_network_to_pdf(results, cutoff = 1.2, remove_duplicates = TRUE, run = run)
    print_norm_molecule_network_to_pdf(results, cutoff = 1.1, remove_duplicates = FALSE, run = run)
    print_norm_molecule_network_to_pdf(results, cutoff = 1.1, remove_duplicates = TRUE, run = run)
  }

# Get a total overview of the data of all runs
  for (run in c(1:10)) {
    results <- read.csv(paste0(resultsdir, 'run', run, '/bayesOpt_normalising_scoreSummaryTable_run', run, '.csv'), row.names = 1)
    results$run <- run
    
    if(run==1) {
      all_results <- read.csv(paste0(resultsdir, 'run0/bayesOpt_normalising_scoreSummaryTable.csv'), row.names = 1)
      all_results$run <- 0
    } else {
      all_results <- rbind(all_results, results)
    }
    
    if (run == 10) {
      write.csv(all_results, paste0(resultsdir, 'all_results_combined.csv'))
      all_results <- as.data.table(all_results)
      
      print_norm_molecule_network_to_pdf(all_results, cutoff = NA, remove_duplicates = FALSE, run = '-all')
      print_norm_molecule_network_to_pdf(all_results, cutoff = NA, remove_duplicates = TRUE, run = '-all')
      print_norm_molecule_network_to_pdf(all_results, cutoff = 2, remove_duplicates = FALSE, run = '-all')
      print_norm_molecule_network_to_pdf(all_results, cutoff = 2, remove_duplicates = TRUE, run = '-all')
      print_norm_molecule_network_to_pdf(all_results, cutoff = 1.5, remove_duplicates = FALSE, run = '-all')
      print_norm_molecule_network_to_pdf(all_results, cutoff = 1.5, remove_duplicates = TRUE, run = '-all')
      print_norm_molecule_network_to_pdf(all_results, cutoff = 1.3, remove_duplicates = FALSE, run = '-all')
      print_norm_molecule_network_to_pdf(all_results, cutoff = 1.3, remove_duplicates = TRUE, run = '-all')
      print_norm_molecule_network_to_pdf(all_results, cutoff = 1.2, remove_duplicates = FALSE, run = '-all')
      print_norm_molecule_network_to_pdf(all_results, cutoff = 1.2, remove_duplicates = TRUE, run = '-all')
      print_norm_molecule_network_to_pdf(all_results, cutoff = 1.1, remove_duplicates = FALSE, run = '-all')
      print_norm_molecule_network_to_pdf(all_results, cutoff = 1.1, remove_duplicates = TRUE, run = '-all')
      
    }
  }
  

# Get an overview of the effect of the different methods and molecules
  all_results_distinct <- all_results %>% distinct(pqn_vsn_used, method_name, comb, .keep_all = TRUE) %>%
                                          mutate(
                                            pqn_vsn_used = ifelse(pqn_vsn_used == "", "None", pqn_vsn_used),
                                            pqn_vsn_used = replace_na(pqn_vsn_used, "None")
                                          ) %>%
                                          mutate(
                                            pqn_vsn_used = factor(pqn_vsn_used),
                                            method_name = factor(method_name)
                                          )
  all_results_distinct$comb <- iconv(all_results_distinct$comb, from = "latin1", to = "UTF-8", sub = "")
  all_results_distinct_split <- all_results_distinct %>% mutate(Count = as.factor(str_count(comb, pattern = ", ")+1)) %>%
    separate_rows(comb, sep = ", ") %>% mutate(Molecule = as.factor(comb))
  
  all_results_distinct_split$full_Molecule_name <- all_results_distinct_split$Molecule
  all_results_distinct_split$Molecule <- abbreviate(all_results_distinct_split$Molecule, minlength = 10)
  
  fwrite(all_results_distinct_split, paste0(resultsdir, "all_combinations-1to8-combinatoric-from_best-selected-and-linSweat_split.csv"))
  # all_results_distinct_split_all_groups <- fread(paste0(resultsdir, 'all_combinations-combinatoric-from_selected_randomFeat_linSweat_linCont_pureRand_split.csv'))
  
  all_results_distinct_split_all_groups <- all_results_distinct_split
  # all_results_distinct_split_all_groups <- rbind(all_results_distinct_split, all_results_distinct_split_all_groups)
  # 
  # fwrite(all_results_distinct_split_all_groups, paste0(resultsdir, 'scoreSummaryTable_1to6-from-selected-randomFeature-pure-random-1to11-linear-sweat_all_combinations.csv'))
  
  
  all_results_distinct_filter <- all_results_distinct %>% 
    mutate(Count = str_count(comb, pattern = ", ")+1)  %>%
    filter(Count > 2, !is.na(Score))
  
  molecule_list <- read.csv(paste0(resultsdir, "filtered_selected-features_with-filtered-good_sweat-features_for-norm.csv")) %>% select(Molecule.Name, origin)
  
  # all_results_distinct_filter$Score <- -all_results_distinct_filter$norm_curve_error * (abs(floor(min(all_results_distinct_filter$corr_factor_diff))) + all_results_distinct_filter$corr_factor_diff)
  # resultsdir <- paste0(resultsdir, 'curve-score/')
  
  print_norm_molecule_network_to_pdf(all_results_distinct_filter, cutoff = NA, remove_duplicates = TRUE, abbrevate = TRUE, abbrevate_len = 15, molecule_list = molecule_list)
  print_norm_molecule_network_to_pdf(all_results_distinct_filter, cutoff = 2, remove_duplicates = TRUE, abbrevate = TRUE, abbrevate_len = 15, molecule_list = molecule_list)
  print_norm_molecule_network_to_pdf(all_results_distinct_filter, cutoff = 1.5, remove_duplicates = TRUE, abbrevate = TRUE, abbrevate_len = 15, molecule_list = molecule_list)
  print_norm_molecule_network_to_pdf(all_results_distinct_filter, cutoff = 1.3, remove_duplicates = TRUE, abbrevate = TRUE, abbrevate_len = 15, molecule_list = molecule_list)
  print_norm_molecule_network_to_pdf(all_results_distinct_filter, cutoff = 1.2, remove_duplicates = TRUE, abbrevate = TRUE, abbrevate_len = 15, molecule_list = molecule_list)
  print_norm_molecule_network_to_pdf(all_results_distinct_filter, cutoff = 1.1, remove_duplicates = TRUE, abbrevate = TRUE, abbrevate_len = 15, molecule_list = molecule_list)
  print_norm_molecule_network_to_pdf(all_results_distinct_filter, cutoff = 1.05, remove_duplicates = TRUE, abbrevate = TRUE, abbrevate_len = 15, molecule_list = molecule_list)
  
  
  
  plot <- ggplot(all_results_distinct, aes(x = pqn_vsn_used, y = log10(-Score))) +
    geom_boxplot() +
    labs(x = "Variance Stabilizing Method", y = "log10(Score)", title = "Boxplot of tested Variance Stabilizing Methods") +
    theme_minimal()
  ggsave(paste0(resultsdir, "pqn-vsn-method-boxplot.png"), plot = plot, width = 8, height = 6, dpi = 1000)
  
  plot <- ggplot(all_results_distinct, aes(x = method_name, y = log10(-Score))) +
    geom_boxplot() +
    labs(x = "Weighting Method", y = "log10(Score)", title = "Boxplot of tested Weighting Methods") +
    theme_minimal()
  ggsave(paste0(resultsdir, "weighting-method-boxplot.png"), plot = plot, width = 8, height = 6, dpi = 1000)
  
  plot <- ggplot(all_results_distinct_split, aes(x = Molecule, y = log10(-Score))) +
    geom_boxplot(outlier.color = NA, width = 0.7) +
    geom_jitter(aes(color = Count), width = 0.1, height = 0, alpha = 0.4, size = 0.6) +  # Color only in jitter
    scale_color_viridis_d() +  # Use discrete color scale
    labs(x = "Molecules", y = "log10(Score)", title = "Boxplot of Tested Molecules for Normalization\ncolored by Molecules in Combination") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste0(resultsdir, "normalisation-molecules-boxplot-colored-by-count.png"), plot = plot, width = 8, height = 6, dpi = 1000)
  
  plot <- ggplot(all_results_distinct_split, aes(x = Molecule, y = log10(-Score))) +
    geom_boxplot(outlier.color = NA, width = 0.7) +
    geom_jitter(aes(color = pqn_vsn_used), width = 0.1, height = 0, alpha = 0.4, size = 0.6) +  # Color only in jitter
    scale_color_viridis_d() +  # Use discrete color scale
    labs(x = "Molecules", y = "log10(Score)", title = "Boxplot of Tested Molecules for Normalization\ncolored by used Variance Stabilizing Method") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste0(resultsdir, "normalisation-molecules-boxplot-colored-by-VarStabMeth.png"), plot = plot, width = 8, height = 6, dpi = 1000)
  
  plot <- ggplot(all_results_distinct_split, aes(x = Molecule, y = log10(-Score))) +
    geom_boxplot(outlier.color = NA, width = 0.8) +
    geom_jitter(aes(color = method_name), width = 0.1, height = 0, alpha = 0.4, size = 0.6) +  # Color only in jitter
    scale_color_viridis_d() +  # Use discrete color scale
    labs(x = "Molecules", y = "log10(Score)", title = "Boxplot of Tested Molecules for Normalization\ncolored by used Weighting Method") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste0(resultsdir, "normalisation-molecules-boxplot-colored-by-WeightMeth.png"), plot = plot, width = 8, height = 6, dpi = 1000)
  
  
# plot good selected and good linear sweat feature combinations  
  
  all_results_annotated <- merge(all_results_distinct_split, molecule_list,
                                 by.x = "comb", by.y = "Molecule.Name", all.x = TRUE)
  
  plot <- ggplot(all_results_annotated, aes(x = Molecule, y = log10(-Score), fill = origin)) +
    geom_boxplot(outlier.color = NA, width = 0.9) +
    scale_fill_manual(values = c("selected feature" = "blue", "linear sweat feature" = "green4")) +  # Custom fill colors
    geom_jitter(aes(color = Count), width = 0.15, height = 0, alpha = 0.4, size = 0.6) +  # Color only in jitter
    scale_color_viridis_d() +  # Use discrete color scale
    labs(x = "Molecules", y = "log10(Score)", title = "Boxplot of Tested Molecules for Normalization\ncolored by Molecules in Combination") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
  ggsave(paste0(resultsdir, "normalisation-molecules-boxplot-colored-by-count-all-groups.jpg"), plot = plot, width = 14, height = 6, dpi = 1000)
  
  all_results_annotated$only_curve_error <- all_results_annotated$norm_curve_error * (abs(floor(min(all_results_annotated$corr_factor_diff))) + all_results_annotated$corr_factor_diff)
  
  
  plot <- ggplot(all_results_annotated, aes(x = Molecule, y = log10(only_curve_error), fill = origin)) +
    geom_boxplot(outlier.color = NA, width = 0.9) +
    scale_fill_manual(values = c("selected feature" = "blue", "linear sweat feature" = "green4")) +  # Custom fill colors
    geom_jitter(aes(color = Count), width = 0.15, height = 0, alpha = 0.4, size = 0.6) +  # Color only in jitter
    scale_color_viridis_d() +  # Use discrete color scale
    labs(x = "Molecules", y = "log10(Curve Score)", title = "Boxplot of Tested Molecules for Normalization\ncolored by Molecules in Combination") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
  ggsave(paste0(resultsdir, "normalisation-molecules-boxplot-colored-by-count-curve-score-all-groups.jpg"), plot = plot, width = 14, height = 6, dpi = 1000)
  
  
  
# plot selected and random molecules
  resultsdir <- paste0(wd, exp, "/results/",'2compartment-bateman-normalization-bruteforce-combinatoric/')
  all_results_distinct_split_all_groups <- as.data.frame(fread(paste0(resultsdir, 'bayesOpt_normalising_scoreSummaryTable_1to6-from-selected-and-random-molecules_all_combinations.csv')))
  all_results_distinct_split_all_groups <- all_results_distinct_split_all_groups %>%
    arrange(group, Molecule) %>%
    mutate(Molecule = factor(Molecule, levels = unique(Molecule)),
           Count = factor(Count, levels = unique(Count)))
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = pqn_vsn_used, y = log10(-Score), fill = group)) +
    geom_boxplot() +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red")) +  # Custom fill colors
    labs(x = "Variance Stabilizing Method", y = "log10(Score)", title = "Boxplot of tested Variance Stabilizing Methods") +
    theme_minimal()
  ggsave(paste0(resultsdir, "pqn-vsn-method-boxplot-all-groups.png"), plot = plot, width = 8, height = 6, dpi = 1000)
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = method_name, y = log10(-Score), fill = group)) +
    geom_boxplot() +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red")) +  # Custom fill colors
    labs(x = "Weighting Method", y = "log10(Score)", title = "Boxplot of tested Weighting Methods") +
    theme_minimal()
  ggsave(paste0(resultsdir, "weighting-method-boxplot-all-groups.png"), plot = plot, width = 8, height = 6, dpi = 1000)
  
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = Molecule, y = log10(-Score), fill = group)) +
    geom_boxplot(outlier.color = NA, width = 0.9) +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red")) +  # Custom fill colors
    geom_jitter(aes(color = Count), width = 0.2, height = 0, alpha = 0.4) +  # Color only in jitter
    scale_color_viridis_d() +  # Use discrete color scale
    labs(x = "Molecules", y = "log10(Score)", title = "Boxplot of Tested Molecules for Normalization\ncolored by Molecules in Combination") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
  ggsave(paste0(resultsdir, "normalisation-molecules-boxplot-colored-by-count-all-groups.jpg"), plot = plot, width = 12, height = 6, dpi = 1000)
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = Molecule, y = log10(-Score), fill = group)) +
    geom_boxplot(outlier.color = NA, width = 0.9) +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red")) +  # Custom fill colors
    geom_jitter(aes(color = pqn_vsn_used), width = 0.2, height = 0, alpha = 0.4) +  # Color only in jitter
    scale_color_viridis_d() +  # Use discrete color scale
    labs(x = "Molecules", y = "log10(Score)", title = "Boxplot of Tested Molecules for Normalization\ncolored by used Variance Stabilizing Method") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
  ggsave(paste0(resultsdir, "normalisation-molecules-boxplot-colored-by-VarStabMeth-all-groups.jpg"), plot = plot, width = 12, height = 6, dpi = 1000)
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = Molecule, y = log10(-Score), fill = group)) +
    geom_boxplot(outlier.color = NA, width = 0.9) +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red")) +  # Custom fill colors
    geom_jitter(aes(color = method_name), width = 0.2, height = 0, alpha = 0.4) +  # Color only in jitter
    scale_color_viridis_d() +  # Use discrete color scale
    labs(x = "Molecules", y = "log10(Score)", title = "Boxplot of Tested Molecules for Normalization\ncolored by used Weighting Method") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
  ggsave(paste0(resultsdir, "normalisation-molecules-boxplot-colored-by-WeightMeth-all-groups.jpg"), plot = plot, width = 12, height = 6, dpi = 1000)
  

  
  
  # plot selected + random + linear sweat molecules
  resultsdir <- paste0(wd, exp, "/results/",'2compartment-bateman-normalization-bruteforce-combinatoric_sweat-linear-features/')
  all_results_distinct_split_all_groups <- as.data.frame(fread(paste0(resultsdir, 'scoreSummaryTable_1to6-from-selected-random-1to11-linear-sweat_all_combinations.csv')))
  all_results_distinct_split_all_groups <- all_results_distinct_split_all_groups %>%
    arrange(group, Molecule) %>%
    mutate(Molecule = factor(Molecule, levels = unique(Molecule)),
           Count = factor(Count, levels = unique(Count)))
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = pqn_vsn_used, y = log10(-Score), fill = group)) +
    geom_boxplot() +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4")) +  # Custom fill colors
    labs(x = "Variance Stabilizing Method", y = "log10(Score)", title = "Boxplot of tested Variance Stabilizing Methods") +
    theme_minimal()
  ggsave(paste0(resultsdir, "pqn-vsn-method-boxplot-all-groups.png"), plot = plot, width = 8, height = 6, dpi = 1000)
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = method_name, y = log10(-Score), fill = group)) +
    geom_boxplot() +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4")) +  # Custom fill colors
    labs(x = "Weighting Method", y = "log10(Score)", title = "Boxplot of tested Weighting Methods") +
    theme_minimal()
  ggsave(paste0(resultsdir, "weighting-method-boxplot-all-groups.png"), plot = plot, width = 8, height = 6, dpi = 1000)
  
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = Molecule, y = log10(-Score), fill = group)) +
    geom_boxplot(outlier.color = NA, width = 0.9) +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4")) +  # Custom fill colors
    geom_jitter(aes(color = method_name), width = 0.1, height = 0, alpha = 0.4, size = 0.6) +  # Color only in jitter
    scale_color_viridis_d() +  # Use discrete color scale
    labs(x = "Molecules", y = "log10(Score)", title = "Boxplot of Tested Molecules for Normalization\ncolored by Molecules in Combination") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
  ggsave(paste0(resultsdir, "normalisation-molecules-boxplot-colored-by-count-all-groups.jpg"), plot = plot, width = 14, height = 6, dpi = 1000)
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = Molecule, y = log10(-Score), fill = group)) +
    geom_boxplot(outlier.color = NA, width = 0.9) +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4")) +  # Custom fill colors
    geom_jitter(aes(color = pqn_vsn_used), width = 0.1, height = 0, alpha = 0.4, size = 0.6) +  # Color only in jitter
    scale_color_viridis_d() +  # Use discrete color scale
    labs(x = "Molecules", y = "log10(Score)", title = "Boxplot of Tested Molecules for Normalization\ncolored by used Variance Stabilizing Method") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
  ggsave(paste0(resultsdir, "normalisation-molecules-boxplot-colored-by-VarStabMeth-all-groups.jpg"), plot = plot, width = 14, height = 6, dpi = 1000)
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = Molecule, y = log10(-Score), fill = group)) +
    geom_boxplot(outlier.color = NA, width = 0.9) +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4")) +  # Custom fill colors
    geom_jitter(aes(color = method_name), width = 0.1, height = 0, alpha = 0.4, size = 0.6) +  # Color only in jitter
    scale_color_viridis_d() +  # Use discrete color scale
    labs(x = "Molecules", y = "log10(Score)", title = "Boxplot of Tested Molecules for Normalization\ncolored by used Weighting Method") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
  ggsave(paste0(resultsdir, "normalisation-molecules-boxplot-colored-by-WeightMeth-all-groups.jpg"), plot = plot, width = 14, height = 6, dpi = 1000)
  
  
  
  # plot selected + random + linear sweat molecules + linear contaminations
  resultsdir <- paste0(wd, exp, "/results/",'2compartment-bateman-normalization-bruteforce-combinatoric_sweat-linear-features/')
  all_results_distinct_split_all_groups <- as.data.frame(fread(paste0(resultsdir, 'all_combinations-combinatoric-from_selected_randomFeat_linSweat_linCont_pureRand_split.csv')))
  all_results_distinct_split_all_groups <- all_results_distinct_split_all_groups %>% filter(group != "pure random values") %>%
    arrange(group, Molecule) %>%
    mutate(Molecule = factor(Molecule, levels = unique(Molecule)),
           Count = factor(Count, levels = unique(Count)))
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = pqn_vsn_used, y = log10(-Score), fill = group)) +
    geom_boxplot() +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4", "linear contamination features" = "orange2")) +  # Custom fill colors
    labs(x = "Variance Stabilizing Method", y = "log10(Score)", title = "Boxplot of tested Variance Stabilizing Methods") +
    theme_minimal()
  ggsave(paste0(resultsdir, "pqn-vsn-method-boxplot-all-groups.png"), plot = plot, width = 8, height = 6, dpi = 1000)
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = method_name, y = log10(-Score), fill = group)) +
    geom_boxplot() +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4", "linear contamination features" = "orange2")) +  # Custom fill colors
    labs(x = "Weighting Method", y = "log10(Score)", title = "Boxplot of tested Weighting Methods") +
    theme_minimal()
  ggsave(paste0(resultsdir, "weighting-method-boxplot-all-groups.png"), plot = plot, width = 8, height = 6, dpi = 1000)
  
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = Molecule, y = log10(-Score), fill = group)) +
    geom_boxplot(outlier.color = NA, width = 0.9) +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4", "linear contamination features" = "orange2")) +  # Custom fill colors
    geom_jitter(aes(color = Count), width = 0.1, height = 0, alpha = 0.4, size = 0.6) +  # Color only in jitter
    scale_color_viridis_d() +  # Use discrete color scale
    labs(x = "Molecules", y = "log10(Score)", title = "Boxplot of Tested Molecules for Normalization\ncolored by Molecules in Combination") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
  ggsave(paste0(resultsdir, "normalisation-molecules-boxplot-colored-by-count-all-groups.jpg"), plot = plot, width = 14, height = 6, dpi = 1000)
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = Molecule, y = log10(-Score), fill = group)) +
    geom_boxplot(outlier.color = NA, width = 0.9) +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4", "linear contamination features" = "orange2")) +  # Custom fill colors
    geom_jitter(aes(color = pqn_vsn_used), width = 0.1, height = 0, alpha = 0.4, size = 0.6) +  # Color only in jitter
    scale_color_viridis_d() +  # Use discrete color scale
    labs(x = "Molecules", y = "log10(Score)", title = "Boxplot of Tested Molecules for Normalization\ncolored by used Variance Stabilizing Method") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
  ggsave(paste0(resultsdir, "normalisation-molecules-boxplot-colored-by-VarStabMeth-all-groups.jpg"), plot = plot, width = 14, height = 6, dpi = 1000)
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = Molecule, y = log10(-Score), fill = group)) +
    geom_boxplot(outlier.color = NA, width = 0.9) +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4", "linear contamination features" = "orange2")) +  # Custom fill colors
    geom_jitter(aes(color = method_name), width = 0.1, height = 0, alpha = 0.4, size = 0.6) +  # Color only in jitter
    scale_color_viridis_d() +  # Use discrete color scale
    labs(x = "Molecules", y = "log10(Score)", title = "Boxplot of Tested Molecules for Normalization\ncolored by used Weighting Method") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
  ggsave(paste0(resultsdir, "normalisation-molecules-boxplot-colored-by-WeightMeth-all-groups.jpg"), plot = plot, width = 14, height = 6, dpi = 1000)
  
  
  
  
  
  # plot selected + random + linear sweat molecules
  resultsdir <- paste0(wd, exp, "/results/",'2compartment-bateman-normalization-bruteforce-combinatoric_sweat-linear-features/')
  all_results_distinct_split_all_groups <- as.data.frame(fread(paste0(resultsdir, 'all_combinations-combinatoric-from_selected_randomFeat_linSweat_linCont_pureRand_split.csv')))
  all_results_distinct_split_all_groups <- all_results_distinct_split_all_groups %>%
    arrange(group, Molecule) %>%
    mutate(Molecule = factor(Molecule, levels = unique(Molecule)),
           Count = factor(Count, levels = unique(Count)))
  resultsdir <- paste0(resultsdir, "including_random_values/")
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = pqn_vsn_used, y = log10(-Score), fill = group)) +
    geom_boxplot() +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4", "linear contamination features" = "orange2", "pure random values" = "violet")) +  # Custom fill colors
    labs(x = "Variance Stabilizing Method", y = "log10(Score)", title = "Boxplot of tested Variance Stabilizing Methods") +
    theme_minimal()
  ggsave(paste0(resultsdir, "pqn-vsn-method-boxplot-all-groups.png"), plot = plot, width = 8, height = 6, dpi = 1000)
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = method_name, y = log10(-Score), fill = group)) +
    geom_boxplot() +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4", "linear contamination features" = "orange2", "pure random values" = "violet")) +  # Custom fill colors
    labs(x = "Weighting Method", y = "log10(Score)", title = "Boxplot of tested Weighting Methods") +
    theme_minimal()
  ggsave(paste0(resultsdir, "weighting-method-boxplot-all-groups.png"), plot = plot, width = 8, height = 6, dpi = 1000)
  
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = Molecule, y = log10(-Score), fill = group)) +
    geom_boxplot(outlier.color = NA, width = 0.9) +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4", "linear contamination features" = "orange2", "pure random values" = "violet")) +  # Custom fill colors
    geom_jitter(aes(color = Count), width = 0.1, height = 0, alpha = 0.4, size = 0.6) +  # Color only in jitter
    scale_color_viridis_d() +  # Use discrete color scale
    labs(x = "Molecules", y = "log10(Score)", title = "Boxplot of Tested Molecules for Normalization\ncolored by Molecules in Combination") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
  ggsave(paste0(resultsdir, "normalisation-molecules-boxplot-colored-by-count-all-groups.jpg"), plot = plot, width = 14, height = 6, dpi = 1000)
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = Molecule, y = log10(-Score), fill = group)) +
    geom_boxplot(outlier.color = NA, width = 0.9) +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4", "linear contamination features" = "orange2", "pure random values" = "violet")) +  # Custom fill colors
    geom_jitter(aes(color = pqn_vsn_used), width = 0.1, height = 0, alpha = 0.4, size = 0.6) +  # Color only in jitter
    scale_color_viridis_d() +  # Use discrete color scale
    labs(x = "Molecules", y = "log10(Score)", title = "Boxplot of Tested Molecules for Normalization\ncolored by used Variance Stabilizing Method") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
  ggsave(paste0(resultsdir, "normalisation-molecules-boxplot-colored-by-VarStabMeth-all-groups.jpg"), plot = plot, width = 14, height = 6, dpi = 1000)
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = Molecule, y = log10(-Score), fill = group)) +
    geom_boxplot(outlier.color = NA, width = 0.9) +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4", "linear contamination features" = "orange2", "pure random values" = "violet")) +  # Custom fill colors
    geom_jitter(aes(color = method_name), width = 0.1, height = 0, alpha = 0.4, size = 0.6) +  # Color only in jitter
    scale_color_viridis_d() +  # Use discrete color scale
    labs(x = "Molecules", y = "log10(Score)", title = "Boxplot of Tested Molecules for Normalization\ncolored by used Weighting Method") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
  ggsave(paste0(resultsdir, "normalisation-molecules-boxplot-colored-by-WeightMeth-all-groups.jpg"), plot = plot, width = 14, height = 6, dpi = 1000)
  
  
  
  
  
  
  # plot only curve scores for selected + random + linear sweat molecules + linear contamination
  resultsdir <- paste0(wd, exp, "/results/",'2compartment-bateman-normalization-bruteforce-combinatoric_sweat-linear-features/')
  all_results_distinct_split_all_groups <- as.data.frame(fread(paste0(resultsdir, 'all_combinations-combinatoric-from_selected_randomFeat_linSweat_linCont_pureRand_split.csv')))
  
  # get only curve error, get norm_curve fit and add corr_factor into it, negative corr_factor means the sd_norm is smaller then sd_raw
  all_results_distinct_split_all_groups$only_curve_error <- all_results_distinct_split_all_groups$norm_curve_error * (abs(floor(min(all_results_distinct_split_all_groups$corr_factor_diff))) + all_results_distinct_split_all_groups$corr_factor_diff)
  
  all_results_distinct_split_all_groups <- all_results_distinct_split_all_groups %>% filter(group != "pure random values") %>%
    arrange(group, Molecule) %>%
    mutate(Molecule = factor(Molecule, levels = unique(Molecule)),
           Count = factor(Count, levels = unique(Count)))
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = pqn_vsn_used, y = log10(only_curve_error), fill = group)) +
    geom_boxplot() +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4", "linear contamination features" = "orange2")) +  # Custom fill colors
    labs(x = "Variance Stabilizing Method", y = "log10(Curve Score)", title = "Boxplot of tested Variance Stabilizing Methods") +
    theme_minimal()
  ggsave(paste0(resultsdir, "pqn-vsn-method-boxplot-curve-score-all-groups.png"), plot = plot, width = 8, height = 6, dpi = 1000)
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = method_name, y = log10(only_curve_error), fill = group)) +
    geom_boxplot() +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4", "linear contamination features" = "orange2")) +  # Custom fill colors
    labs(x = "Weighting Method", y = "log10(Curve Score)", title = "Boxplot of tested Weighting Methods") +
    theme_minimal()
  ggsave(paste0(resultsdir, "weighting-method-boxplot-curve-score-all-groups.png"), plot = plot, width = 8, height = 6, dpi = 1000)
  
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = Molecule, y = log10(only_curve_error), fill = group)) +
    geom_boxplot(outlier.color = NA, width = 0.9) +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4", "linear contamination features" = "orange2")) +  # Custom fill colors
    geom_jitter(aes(color = method_name), width = 0.1, height = 0, alpha = 0.4, size = 0.6) +  # Color only in jitter
    scale_color_viridis_d() +  # Use discrete color scale
    labs(x = "Molecules", y = "log10(Curve Score)", title = "Boxplot of Tested Molecules for Normalization\ncolored by Molecules in Combination") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
  ggsave(paste0(resultsdir, "normalisation-molecules-boxplot-colored-by-count-curve-score-all-groups.jpg"), plot = plot, width = 14, height = 6, dpi = 1000)
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = Molecule, y = log10(only_curve_error), fill = group)) +
    geom_boxplot(outlier.color = NA, width = 0.9) +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4", "linear contamination features" = "orange2")) +  # Custom fill colors
    geom_jitter(aes(color = pqn_vsn_used), width = 0.1, height = 0, alpha = 0.4, size = 0.6) +  # Color only in jitter
    scale_color_viridis_d() +  # Use discrete color scale
    labs(x = "Molecules", y = "log10(Curve Score)", title = "Boxplot of Tested Molecules for Normalization\ncolored by used Variance Stabilizing Method") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
  ggsave(paste0(resultsdir, "normalisation-molecules-boxplot-colored-by-VarStabMeth-curve-score-all-groups.jpg"), plot = plot, width = 14, height = 6, dpi = 1000)
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = Molecule, y = log10(only_curve_error), fill = group)) +
    geom_boxplot(outlier.color = NA, width = 0.9) +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4", "linear contamination features" = "orange2")) +  # Custom fill colors
    geom_jitter(aes(color = method_name), width = 0.1, height = 0, alpha = 0.4, size = 0.6) +  # Color only in jitter
    scale_color_viridis_d() +  # Use discrete color scale
    labs(x = "Molecules", y = "log10(Curve Score)", title = "Boxplot of Tested Molecules for Normalization\ncolored by used Weighting Method") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
  ggsave(paste0(resultsdir, "normalisation-molecules-boxplot-colored-by-WeightMeth-curve-score-all-groups.jpg"), plot = plot, width = 14, height = 6, dpi = 1000)
  
  
  
  # plot only curve scores for selected + random + linear sweat + linear contamination + random values
  resultsdir <- paste0(wd, exp, "/results/",'2compartment-bateman-normalization-bruteforce-combinatoric_sweat-linear-features/')
  all_results_distinct_split_all_groups <- as.data.frame(fread(paste0(resultsdir, 'all_combinations-combinatoric-from_selected_randomFeat_linSweat_linCont_pureRand_split.csv')))
  resultsdir <- paste0(resultsdir, "including_random_values/")
  
  
  # get only curve error, get norm_curve fit and add corr_factor into it, negative corr_factor means the sd_norm is smaller then sd_raw
  all_results_distinct_split_all_groups$only_curve_error <- all_results_distinct_split_all_groups$norm_curve_error * (abs(floor(min(all_results_distinct_split_all_groups$corr_factor_diff))) + all_results_distinct_split_all_groups$corr_factor_diff)
  
  all_results_distinct_split_all_groups <- all_results_distinct_split_all_groups %>%
    arrange(group, Molecule) %>%
    mutate(Molecule = factor(Molecule, levels = unique(Molecule)),
           Count = factor(Count, levels = unique(Count)))
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = pqn_vsn_used, y = log10(only_curve_error), fill = group)) +
    geom_boxplot() +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4", "linear contamination features" = "orange2", "pure random values" = "violet")) +  # Custom fill colors
    labs(x = "Variance Stabilizing Method", y = "log10(Curve Score)", title = "Boxplot of tested Variance Stabilizing Methods") +
    theme_minimal()
  ggsave(paste0(resultsdir, "pqn-vsn-method-boxplot-curve-score-all-groups.png"), plot = plot, width = 8, height = 6, dpi = 1000)
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = method_name, y = log10(only_curve_error), fill = group)) +
    geom_boxplot() +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4", "linear contamination features" = "orange2", "pure random values" = "violet")) +  # Custom fill colors
    labs(x = "Weighting Method", y = "log10(Curve Score)", title = "Boxplot of tested Weighting Methods") +
    theme_minimal()
  ggsave(paste0(resultsdir, "weighting-method-boxplot-curve-score-all-groups.png"), plot = plot, width = 8, height = 6, dpi = 1000)
  
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = Molecule, y = log10(only_curve_error), fill = group)) +
    geom_boxplot(outlier.color = NA, width = 0.9) +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4", "linear contamination features" = "orange2", "pure random values" = "violet")) +  # Custom fill colors
    geom_jitter(aes(color = method_name), width = 0.1, height = 0, alpha = 0.4, size = 0.6) +  # Color only in jitter
    scale_color_viridis_d() +  # Use discrete color scale
    labs(x = "Molecules", y = "log10(Curve Score)", title = "Boxplot of Tested Molecules for Normalization\ncolored by Molecules in Combination") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
  ggsave(paste0(resultsdir, "normalisation-molecules-boxplot-colored-by-count-curve-score-all-groups.jpg"), plot = plot, width = 14, height = 6, dpi = 1000)
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = Molecule, y = log10(only_curve_error), fill = group)) +
    geom_boxplot(outlier.color = NA, width = 0.9) +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4", "linear contamination features" = "orange2", "pure random values" = "violet")) +  # Custom fill colors
    geom_jitter(aes(color = pqn_vsn_used), width = 0.1, height = 0, alpha = 0.4, size = 0.6) +  # Color only in jitter
    scale_color_viridis_d() +  # Use discrete color scale
    labs(x = "Molecules", y = "log10(Curve Score)", title = "Boxplot of Tested Molecules for Normalization\ncolored by used Variance Stabilizing Method") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
  ggsave(paste0(resultsdir, "normalisation-molecules-boxplot-colored-by-VarStabMeth-curve-score-all-groups.jpg"), plot = plot, width = 14, height = 6, dpi = 1000)
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = Molecule, y = log10(only_curve_error), fill = group)) +
    geom_boxplot(outlier.color = NA, width = 0.9) +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4", "linear contamination features" = "orange2", "pure random values" = "violet")) +  # Custom fill colors
    geom_jitter(aes(color = method_name), width = 0.1, height = 0, alpha = 0.4, size = 0.6) +  # Color only in jitter
    scale_color_viridis_d() +  # Use discrete color scale
    labs(x = "Molecules", y = "log10(Curve Score)", title = "Boxplot of Tested Molecules for Normalization\ncolored by used Weighting Method") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
  ggsave(paste0(resultsdir, "normalisation-molecules-boxplot-colored-by-WeightMeth-curve-score-all-groups.jpg"), plot = plot, width = 14, height = 6, dpi = 1000)
  
  
  
  
  # filter the methods of all groups
  # plot selected + random + linear sweat molecules
  resultsdir <- paste0(wd, exp, "/results/",'2compartment-bateman-normalization-bruteforce-combinatoric_sweat-linear-features/')
  all_results_distinct_split_all_groups <- as.data.frame(fread(paste0(resultsdir, 'all_combinations-combinatoric-from_selected_randomFeat_linSweat_linCont_pureRand_split.csv')))
  all_results_distinct_split_all_groups <- all_results_distinct_split_all_groups %>%
                                            filter(!method_name == 'Log2', !method_name == 'Sqrt') %>% filter(pqn_vsn_used == 'post-VSN') %>% filter(Count > 2, Count < 7)
  
  
  all_results_distinct_split_all_groups <- all_results_distinct_split_all_groups %>%
    arrange(group, Molecule) %>%
    mutate(Molecule = factor(Molecule, levels = unique(Molecule)),
           Count = factor(Count, levels = unique(Count)))
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = pqn_vsn_used, y = log10(-Score), fill = group)) +
    geom_boxplot() +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4", "linear contamination features" = "orange2", "pure random values" = "violet")) +  # Custom fill colors
    labs(x = "Variance Stabilizing Method", y = "log10(Score)", title = "Boxplot of tested Variance Stabilizing Methods") +
    theme_minimal()
  ggsave(paste0(resultsdir, "pqn-vsn-method-boxplot-all-groups-filtered-data2.png"), plot = plot, width = 8, height = 6, dpi = 2000)
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = method_name, y = log10(-Score), fill = group)) +
    geom_boxplot() +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4", "linear contamination features" = "orange2", "pure random values" = "violet")) +  # Custom fill colors
    labs(x = "Weighting Method", y = "log10(Score)", title = "Boxplot of tested Weighting Methods") +
    theme_minimal()
  ggsave(paste0(resultsdir, "weighting-method-boxplot-all-groups-filtered-data2.png"), plot = plot, width = 8, height = 6, dpi = 2000)
  
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = Molecule, y = log10(-Score), fill = group)) +
    geom_boxplot(outlier.color = NA, width = 0.9) +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4", "linear contamination features" = "orange2", "pure random values" = "violet")) +  # Custom fill colors
    geom_jitter(aes(color = method_name), width = 0.1, height = 0, alpha = 0.4, size = 0.6) +  # Color only in jitter
    scale_color_viridis_d() +  # Use discrete color scale
    labs(x = "Molecules", y = "log10(Score)", title = "Boxplot of Tested Molecules for Normalization\ncolored by Molecules in Combination") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
  ggsave(paste0(resultsdir, "normalisation-molecules-boxplot-colored-by-count-all-groups-filtered-data2.jpg"), plot = plot, width = 14, height = 6, dpi = 2000)
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = Molecule, y = log10(-Score), fill = group)) +
    geom_boxplot(outlier.color = NA, width = 0.9) +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4", "linear contamination features" = "orange2", "pure random values" = "violet")) +  # Custom fill colors
    geom_jitter(aes(color = pqn_vsn_used), width = 0.1, height = 0, alpha = 0.4, size = 0.6) +  # Color only in jitter
    scale_color_viridis_d() +  # Use discrete color scale
    labs(x = "Molecules", y = "log10(Score)", title = "Boxplot of Tested Molecules for Normalization\ncolored by used Variance Stabilizing Method") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
  ggsave(paste0(resultsdir, "normalisation-molecules-boxplot-colored-by-VarStabMeth-all-groups-filtered-data2.jpg"), plot = plot, width = 14, height = 6, dpi = 2000)
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = Molecule, y = log10(-Score), fill = group)) +
    geom_boxplot(outlier.color = NA, width = 0.9) +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4", "linear contamination features" = "orange2", "pure random values" = "violet")) +  # Custom fill colors
    geom_jitter(aes(color = method_name), width = 0.1, height = 0, alpha = 0.4, size = 0.6) +  # Color only in jitter
    scale_color_viridis_d() +  # Use discrete color scale
    labs(x = "Molecules", y = "log10(Score)", title = "Boxplot of Tested Molecules for Normalization\ncolored by used Weighting Method") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
  ggsave(paste0(resultsdir, "normalisation-molecules-boxplot-colored-by-WeightMeth-all-groups-filtered-data2.jpg"), plot = plot, width = 14, height = 6, dpi = 2000)
  
  
  
  # get only curve error, get norm_curve fit and add corr_factor into it, negative corr_factor means the sd_norm is smaller then sd_raw
  # filter the methods of all groups
  # plot selected + random + linear sweat molecules
  resultsdir <- paste0(wd, exp, "/results/",'2compartment-bateman-normalization-bruteforce-combinatoric_sweat-linear-features/')
  all_results_distinct_split_all_groups <- as.data.frame(fread(paste0(resultsdir, 'all_combinations-combinatoric-from_selected_randomFeat_linSweat_linCont_pureRand_split.csv')))
  all_results_distinct_split_all_groups <- all_results_distinct_split_all_groups %>%
    filter(!method_name == 'Log2', !method_name == 'Sqrt') %>% filter(pqn_vsn_used == 'post-VSN') %>% filter(Count > 2, Count < 7)
  
  all_results_distinct_split_all_groups$only_curve_error <- all_results_distinct_split_all_groups$norm_curve_error * (abs(floor(min(all_results_distinct_split_all_groups$corr_factor_diff))) + all_results_distinct_split_all_groups$corr_factor_diff)
  
  all_results_distinct_split_all_groups <- all_results_distinct_split_all_groups %>%
    arrange(group, Molecule) %>%
    mutate(Molecule = factor(Molecule, levels = unique(Molecule)),
           Count = factor(Count, levels = unique(Count)))
  
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = pqn_vsn_used, y = log10(-Score), fill = group)) +
    geom_boxplot() +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4", "linear contamination features" = "orange2", "pure random values" = "violet")) +  # Custom fill colors
    labs(x = "Variance Stabilizing Method", y = "log10(Score)", title = "Boxplot of tested Variance Stabilizing Methods") +
    theme_minimal()
  ggsave(paste0(resultsdir, "pqn-vsn-method-boxplot-all-groups-curve-score-filtered-data2.png"), plot = plot, width = 8, height = 6, dpi = 2000)
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = method_name, y = log10(-Score), fill = group)) +
    geom_boxplot() +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4", "linear contamination features" = "orange2", "pure random values" = "violet")) +  # Custom fill colors
    labs(x = "Weighting Method", y = "log10(Score)", title = "Boxplot of tested Weighting Methods") +
    theme_minimal()
  ggsave(paste0(resultsdir, "weighting-method-boxplot-all-groups-curve-score-filtered-data2.png"), plot = plot, width = 8, height = 6, dpi = 2000)
  
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = Molecule, y = log10(-Score), fill = group)) +
    geom_boxplot(outlier.color = NA, width = 0.9) +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4", "linear contamination features" = "orange2", "pure random values" = "violet")) +  # Custom fill colors
    geom_jitter(aes(color = method_name), width = 0.1, height = 0, alpha = 0.4, size = 0.6) +  # Color only in jitter
    scale_color_viridis_d() +  # Use discrete color scale
    labs(x = "Molecules", y = "log10(Score)", title = "Boxplot of Tested Molecules for Normalization\ncolored by Molecules in Combination") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
  ggsave(paste0(resultsdir, "normalisation-molecules-boxplot-colored-by-count-all-groups-curve-score-filtered-data2.jpg"), plot = plot, width = 14, height = 6, dpi = 2000)
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = Molecule, y = log10(-Score), fill = group)) +
    geom_boxplot(outlier.color = NA, width = 0.9) +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4", "linear contamination features" = "orange2", "pure random values" = "violet")) +  # Custom fill colors
    geom_jitter(aes(color = pqn_vsn_used), width = 0.1, height = 0, alpha = 0.4, size = 0.6) +  # Color only in jitter
    scale_color_viridis_d() +  # Use discrete color scale
    labs(x = "Molecules", y = "log10(Score)", title = "Boxplot of Tested Molecules for Normalization\ncolored by used Variance Stabilizing Method") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
  ggsave(paste0(resultsdir, "normalisation-molecules-boxplot-colored-by-VarStabMeth-all-groups-curve-score-filtered-data2.jpg"), plot = plot, width = 14, height = 6, dpi = 2000)
  
  plot <- ggplot(all_results_distinct_split_all_groups, aes(x = Molecule, y = log10(-Score), fill = group)) +
    geom_boxplot(outlier.color = NA, width = 0.9) +
    scale_fill_manual(values = c("selected features" = "blue", "random features" = "red", "linear sweat features" = "green4", "linear contamination features" = "orange2", "pure random values" = "violet")) +  # Custom fill colors
    geom_jitter(aes(color = method_name), width = 0.1, height = 0, alpha = 0.4, size = 0.6) +  # Color only in jitter
    scale_color_viridis_d() +  # Use discrete color scale
    labs(x = "Molecules", y = "log10(Score)", title = "Boxplot of Tested Molecules for Normalization\ncolored by used Weighting Method") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
  ggsave(paste0(resultsdir, "normalisation-molecules-boxplot-colored-by-WeightMeth-all-groups-curve-score-filtered-data2.jpg"), plot = plot, width = 14, height = 6, dpi = 2000)
  
  
  
  
  
  
  # filter the methods of all groups
  # plot selected + random + linear sweat molecules
  resultsdir <- paste0(wd, exp, "/results/",'2compartment-bateman-normalization-bruteforce-combinatoric_sweat-linear-features/')
  all_results_distinct_all_groups <- as.data.frame(fread(file = paste0(resultsdir, 'scoreSummaryTable_all_combinations-1to6combinatoric-from_selected_randomFeat_linSweat_linCont_pureRand.csv')))
  all_results_distinct_all_groups <- all_results_distinct_all_groups %>% 
    mutate(Count = str_count(comb, pattern = ", ")+1)  %>%
    filter(!method_name == 'Log2', !method_name == 'Sqrt') %>% filter(pqn_vsn_used == 'post-VSN') %>% filter(Count > 2, Count < 7, !is.na(Score))
  all_results_distinct_all_groups <- as.data.table(all_results_distinct_all_groups)
  
  print_norm_molecule_network_to_pdf(all_results_distinct_all_groups, cutoff = NA, remove_duplicates = TRUE, abbrevate = TRUE)
  print_norm_molecule_network_to_pdf(all_results_distinct_all_groups, cutoff = 2, remove_duplicates = TRUE, abbrevate = TRUE)
  print_norm_molecule_network_to_pdf(all_results_distinct_all_groups, cutoff = 1.5, remove_duplicates = TRUE, abbrevate = TRUE)
  print_norm_molecule_network_to_pdf(all_results_distinct_all_groups, cutoff = 1.3, remove_duplicates = TRUE, abbrevate = TRUE)
  print_norm_molecule_network_to_pdf(all_results_distinct_all_groups, cutoff = 1.2, remove_duplicates = TRUE, abbrevate = TRUE)
  print_norm_molecule_network_to_pdf(all_results_distinct_all_groups, cutoff = 1.1, remove_duplicates = TRUE, abbrevate = TRUE)
  print_norm_molecule_network_to_pdf(all_results_distinct_all_groups, cutoff = 1.05, remove_duplicates = TRUE, abbrevate = TRUE)
  
  
  resultsdir <- paste0(wd, exp, "/results/",'2compartment-bateman-normalization-bruteforce-combinatoric_sweat-linear-features/')
  all_results_distinct_all_groups <- as.data.frame(fread(file = paste0(resultsdir, 'scoreSummaryTable_all_combinations-1to6combinatoric-from_selected_randomFeat_linSweat_linCont_pureRand.csv')))
  all_results_distinct_all_groups <- all_results_distinct_all_groups %>% 
    mutate(Count = str_count(comb, pattern = ", ")+1)  %>%
    filter(!method_name == 'Log2', !method_name == 'Sqrt') %>% filter(pqn_vsn_used == 'post-VSN') %>% filter(Count > 2, Count < 7, !is.na(Score))
  
  resultsdir <- paste0(resultsdir, "filtered-selected/")
  all_results_distinct_all_groups <- as.data.table(all_results_distinct_all_groups %>% filter(group == "selected features"))
  
  print_norm_molecule_network_to_pdf(all_results_distinct_all_groups, cutoff = NA, remove_duplicates = TRUE)
  print_norm_molecule_network_to_pdf(all_results_distinct_all_groups, cutoff = 2, remove_duplicates = TRUE)
  print_norm_molecule_network_to_pdf(all_results_distinct_all_groups, cutoff = 1.5, remove_duplicates = TRUE)
  print_norm_molecule_network_to_pdf(all_results_distinct_all_groups, cutoff = 1.3, remove_duplicates = TRUE)
  print_norm_molecule_network_to_pdf(all_results_distinct_all_groups, cutoff = 1.2, remove_duplicates = TRUE)
  print_norm_molecule_network_to_pdf(all_results_distinct_all_groups, cutoff = 1.1, remove_duplicates = TRUE)
  print_norm_molecule_network_to_pdf(all_results_distinct_all_groups, cutoff = 1.05, remove_duplicates = TRUE)
  
  
  resultsdir <- paste0(wd, exp, "/results/",'2compartment-bateman-normalization-bruteforce-combinatoric_sweat-linear-features/')
  all_results_distinct_all_groups <- as.data.frame(fread(file = paste0(resultsdir, 'scoreSummaryTable_all_combinations-1to6combinatoric-from_selected_randomFeat_linSweat_linCont_pureRand.csv')))
  all_results_distinct_all_groups <- all_results_distinct_all_groups %>% 
    mutate(Count = str_count(comb, pattern = ", ")+1)  %>%
    filter(!method_name == 'Log2', !method_name == 'Sqrt') %>% filter(pqn_vsn_used == 'post-VSN') %>% filter(Count > 2, Count < 7, !is.na(Score))
  
  resultsdir <- paste0(resultsdir, "filtered-sweat/")
  all_results_distinct_all_groups <- as.data.table(all_results_distinct_all_groups %>% filter(group == "linear sweat features"))
  
  print_norm_molecule_network_to_pdf(all_results_distinct_all_groups, cutoff = NA, remove_duplicates = TRUE, abbrevate = TRUE, abbrevate_len = 15)
  print_norm_molecule_network_to_pdf(all_results_distinct_all_groups, cutoff = 2, remove_duplicates = TRUE, abbrevate = TRUE, abbrevate_len = 15)
  print_norm_molecule_network_to_pdf(all_results_distinct_all_groups, cutoff = 1.5, remove_duplicates = TRUE, abbrevate = TRUE, abbrevate_len = 15)
  print_norm_molecule_network_to_pdf(all_results_distinct_all_groups, cutoff = 1.3, remove_duplicates = TRUE, abbrevate = TRUE, abbrevate_len = 15)
  print_norm_molecule_network_to_pdf(all_results_distinct_all_groups, cutoff = 1.2, remove_duplicates = TRUE, abbrevate = TRUE, abbrevate_len = 15)
  print_norm_molecule_network_to_pdf(all_results_distinct_all_groups, cutoff = 1.1, remove_duplicates = TRUE, abbrevate = TRUE, abbrevate_len = 15)
  print_norm_molecule_network_to_pdf(all_results_distinct_all_groups, cutoff = 1.05, remove_duplicates = TRUE, abbrevate = TRUE, abbrevate_len = 15)
  
  