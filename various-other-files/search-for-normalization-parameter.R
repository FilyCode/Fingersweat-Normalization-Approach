#source("X:/Studenten_Schueler/Philipp_Trollmann/R Project/MS-data-analysis_functions.R")
#source("X:/Studenten_Schueler/Philipp_Trollmann/R Project/feature-visualization-script.R")
source("W:/tmp_PT/R Project/MS-data-analysis_functions.R")
source("W:/tmp_PT/R Project/feature-visualization-script.R")

# set all the experiment and file names
  #wd = "X:/Studenten_Schueler/Philipp_Trollmann/Experiments/"
  wd = "D:/tmp_PT/Experiments/"
  #wd = "E:/Experiments/"
  setwd(wd)
  
  exp = "FiS-45min-Abstandsmessung"
  #exp = "MTX1-8"
  #exp = "Coffein-Samplings_SimonBAProject"
  #exp = "Julia_Coffein"
  #exp = "CoffeeMatcha-Sampling"
  #exp = "Weinstudie_Chrom-Kim"
  #exp = "Weinstudie-all-samples"
  #exp = "Weinstudie-Chrom"
  
  datadir = paste0(wd, exp, "/data/")
  resultsdir = paste0(wd, exp, "/results/")
  
  tlists = paste0(wd, "Transition-Lists/")
  tfile <- "CaffeineTransitionList_CoffeeOnly.csv"
  #tfile <- "CaffeineTransitionList_biological-norm-feature-check.csv"
  #tfile <- "MTX1-7_targeted-features.csv"
  #tfile <- "CaffeineTransitionList_CoffeeOnly_shift-Simon-Caffein.csv"
  #tfile <- "CaffeineTransitionList_CoffeeOnly_shift-Julia-Caffein.csv"
  #tfile <- "CaffeineTransitionList_biological-norm-feature-check_shift-Simon-Data.csv"
  #tfile <- "WinestudyTransitionList.csv"
  #tfile <- "Screening_TransitionList.csv"
  
  info_file <- "/FiS_45min-times_long-names.csv"
  #info_file <- "/MTX1-7_Samplelist.csv"
  info_file <- "/FiS_Caffein-times.csv"
  #info_file <- "/FiS_Caffein-Matcha-times.csv"
  #info_file <- "/FiS_Caffein-Matcha-times_coffee-matcha-split.csv"
  #info_file <- "/Weinstudie-samplelist.csv"
  
  #normalization_list  = paste0(wd, "normalization/featurelist-normalization_52-correlating-features.csv")
  #normalization_list  = paste0(wd, "normalization/biological_normalization_molecule_list.csv")
  #normalization_list  = paste0(wd, "normalization/biological_normalization_molecule_list_filtered.csv")
  #normalization_list  = paste0(wd, "normalization/biological_normalization_molecule_list_filtered_größer_0.6rt.csv")
  #normalization_list  = paste0(wd, "normalization/biological_normalization_molecule_list_filtered_größer_0.6rt_without-niacinamid.csv")
  #normalization_list  = paste0(wd, "normalization/biological_normalization_molecule_list_filtered_shift-Simon-Caffein.csv")
  #normalization_list  = paste0(wd, "normalization/biological_normalization_molecule_list_filter-for-combination-test_shiftSimon.csv")
  #normalization_list  = paste0(wd, "normalization/biological_normalization_molecule_list_filter-for-combination-test_shiftJulia.csv")
  #normalization_list  = paste0(wd, "normalization/biological_normalization_molecule_list_filter-for-combination-test.csv")
  normalization_list  = paste0(wd, "normalization/biological_normalization_molecule_list_filter-for-combination-test_reduced.csv")
  normalization_list  = paste0(wd, "normalization/biological_normalization_molecule_list_filter-for-combination-test_shiftSimon_reduced.csv")
  normalization_list  = paste0(wd, "normalization/biological_normalization_molecule_list_filter-for-combination-test_shiftJulia_reduced.csv")
  #normalization_list  = paste0(wd, "normalization/Tyrosine-only.csv")
  
###############################
## Start of Data Processing  ##
###############################

# Get the whole processed data set from file (already prepared csv file with mzMine)
  exp_data <- process_merged_files(datadir)
  all_data <- exp_data$data
  df_ids <- exp_data$df_ids


# Define the molecules to look at from transition list and extract them from the whole data set
  info_file_dir <- paste0(wd, exp, info_file)
  transitionFile <- paste0(tlists, tfile)
  targeted_experiment_data <- extract_feature_list(all_data, df_ids, transitionFile)
  significant_abundant_features_object <- get_significant_abundant_features(datadir, exp_data, info_file_dir, resultsdir, figuredir, exp, tlists)
  all_needed_features <- significant_abundant_features_object$significant_abundant_features
  se <- significant_abundant_features_object$se
  
  
  pqn_data <- as.data.frame(assay(se, i = "log2_pqn"))
  vsn_data <- as.data.frame(assay(se, i = "vsn"))
  pqn_data <- 2**(pqn_data-20) # make from log2 to "normal" again
  vsn_data <- 2**vsn_data # make from log2 to "normal" again
  targeted_experiment_data_pqn <- extract_feature_list(pqn_data, df_ids, transitionFile)
  targeted_experiment_data_vsn <- extract_feature_list(vsn_data, df_ids, transitionFile)
  
  # remove duplicates
  targeted_experiment_data <- targeted_experiment_data[c(1,3,5),]
  targeted_experiment_data_pqn <- targeted_experiment_data_pqn[c(1,3,5),]
  targeted_experiment_data_vsn <- targeted_experiment_data_vsn[c(1,3,5),]
  
  
  # if you want to do it with the SIRIUS annotation
    norm_feature_names <- read.csv(normalization_list)
    if (ncol(norm_feature_names) < 2) {
      norm_feature_names <- read.csv2(normalization_list)
    }    
    norm_feature_list <- all_needed_features[all_needed_features$Annotation %in% norm_feature_names$Molecule.Name,]
    normalization_molecules <- norm_feature_names$Molecule.Name
  
  # if you want to look for the targeted molecules
    norm_feature_names <- read.csv(normalization_list)
    if (ncol(norm_feature_names) < 2) {
      norm_feature_names <- read.csv2(normalization_list)
    }
    norm_feature_list <- extract_feature_list(all_data, df_ids, normalization_list)
    normalization_molecules <- norm_feature_names$Molecule.Name
    normalization_molecules_id <- norm_feature_list$id
  
    norm_feature_list_pqn <- extract_feature_list(pqn_data, df_ids, normalization_list)
    norm_feature_list_vsn <- extract_feature_list(vsn_data, df_ids, normalization_list)
    
    # remove duplicates
    norm_feature_list <- norm_feature_list %>% distinct(id, .keep_all = TRUE)
    norm_feature_list_pqn <- norm_feature_list_pqn %>% distinct(id, .keep_all = TRUE)
    norm_feature_list_vsn <- norm_feature_list_vsn %>% distinct(id, .keep_all = TRUE)
    
    
  # make PCA for normalization molecules
    # Function to replace NAs with the median of the column
    replace_na_with_median <- function(x) {
      if (is.numeric(x)) {
        x[is.na(x)] <- median(x, na.rm = TRUE)
      }
      return(x)
    }
    
    norm_mol_df <- all_needed_features %>%
      mutate(Group = ifelse(id %in% normalization_molecules_id, "Normalization Feature", "General Feature"))
    norm_mol_df <-  norm_mol_df[-c(1:9)]
    norm_mol_df <- as.data.frame(lapply(norm_mol_df, replace_na_with_median))
    #rownames(norm_mol_df_for_pca) <- norm_mol_df_for_pca$Molecule.Name
    norm_mol_df_for_pca <- norm_mol_df
    norm_mol_df_for_pca <- norm_mol_df[-length(norm_mol_df_for_pca)]
    norm_mol_df_for_pca <- scale(norm_mol_df_for_pca)
    pca_result <- prcomp(norm_mol_df_for_pca, center = TRUE, scale. = TRUE)
    pca_data <- data.frame(pca_result$x, Group = norm_mol_df$Group)
    summary(pca_result)
    # Visualize the PCA results
      fviz_eig(pca_result)
      fviz_pca_ind(pca_result, geom.ind = "point", col.ind = norm_mol_df$Group, 
                   palette = c("#00AFBB", "#FC4E07"), 
                   addEllipses = TRUE, 
                   legend.title = "Groups",
                   repel = TRUE)
      #fviz_pca_var(pca_result, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)
      autoplot(pca_result, data = norm_mol_df,
               color = "Group",
               label = TRUE, label.size = 3,
               frame = TRUE, frame.type = 'norm',
               #loadings = TRUE, loadings.colour = 'blue',
               #loadings.label = TRUE, loadings.label.size = 3
               )
    
  unique_molecules <- unique(targeted_experiment_data$Molecule.Name)
  add_info_file <- paste0(wd, exp, info_file)
  add_info <- read.csv2(add_info_file)
  
# Plot all molecules in transition list and check various normalization strategies
  # for 52 correlating features
  plot_to_pdf_normalization_methods_check(targeted_experiment_data, unique_molecules, normalization_molecules, all_needed_features, add_info_file, resultsdir)
  plot_to_pdf_normalization_methods_check_per_person(targeted_experiment_data, unique_molecules, normalization_molecules, all_needed_features, add_info_file, resultsdir)
  # for targeted biological features
  plot_to_pdf_normalization_methods_check(targeted_experiment_data, unique_molecules, normalization_molecules, all_needed_features, add_info_file, resultsdir, norm_feature_list)
  plot_to_pdf_normalization_methods_check_per_person(targeted_experiment_data, unique_molecules, normalization_molecules, all_needed_features, add_info_file, resultsdir, norm_feature_list)


# Plot all molecules in transition list against time and normalize against various molecules to check which of them have potential to be 
  targeted_exp_data <- targeted_experiment_data
  # If Urocanic Acid is used in list than need to sum up cis and trans forms
    # Sum up cis and trans
    summed_row <- targeted_exp_data %>%
      filter(Molecule.Name %in% c('trans-Urocanic acid', 'cis-Urocanic acid')) %>%
      summarise(across(5:(ncol(targeted_exp_data)-1), sum)) %>%
      mutate(Molecule.Name = 'Urocanic acid sum')
    
    # Extract the metadata from trans-Urocanic acid
    trans_metadata <- targeted_exp_data %>%
      filter(Molecule.Name == 'trans-Urocanic acid') %>%
      select(id, rt, mz, charge) 
    
    # Combine the metadata and the summed values
    summed_row <- bind_cols(trans_metadata, summed_row)
  
    # Remove the original rows and add the new row
    targeted_exp_data <- targeted_exp_data %>%
      filter(!Molecule.Name %in% c('trans-Urocanic acid', 'cis-Urocanic acid')) %>%
      bind_rows(summed_row)
  
    row.names(targeted_exp_data) <- targeted_exp_data$id
    
  # Transform data to plot them
  experiment_data_plots <- targeted_exp_data %>%
    pivot_longer(cols = -c(id, rt, mz, charge, Molecule.Name), names_to = "Sample", values_to = "Area") %>%
    mutate(Sample = factor(Sample, levels = unique(Sample))) 
  
  # Add additional information (like sample times, different condition like different paper etc)
  add_info_file <- paste0(wd, exp, info_file)
  add_info <- read.csv(add_info_file)
  if (length(add_info) < 2) {
    add_info <- read.csv2(add_info_file)
  }
  
  experiment_data_plots <- merge(experiment_data_plots, add_info, by.x = "Sample", by.y = "Sample", all.x = TRUE)
  
  
  # change the data values to the according format
  experiment_data_plots$date <- as.Date(experiment_data_plots$date, format = "%d.%m.%Y")
  
  # Convert time column to minutes between sampling
  experiment_data_plots$time <- sapply(experiment_data_plots$time, convert_to_minutes)
  
  unique_molecules <- unique(experiment_data_plots$Molecule.Name)

  #normalization_molecules <- c('Tyrosine', 'Nicotinamide', 'Hypoxanthine', 'Oleamide', 'Glutaric Acid', 
  #                             'Creatine', 'Creatinine', 'Uric Acid', 'Serin')

  normalization_molecules <- c("DL-Tyrosine", "DL-Tryptophan", "Hippuric acid", "Urocanic acid sum", "Niacinamide",
                               "Succinic acid", "Glutaric Acid", "Creatinine", "Carnitine", "Uric Acid", "Uracil", "Hypoxanthine",
                               "Serine", "Glutamic acid", "Iso-Leucine", "DL-Phenylalanine", "Suberic acid", "Carnosine", "Xanthine",
                               "Adenosine", "Citrulline", "Glutamine", "Asparagine", "Aspartate", "Histidine", "Arginine", "Alanine",
                               "Glycine", "Betaine1", "Betaine2", "Allantoin")
  
  
  final_raw_results <- readRDS("D:/tmp_PT/Experiments/Coffein-Samplings_SimonBAProject/results/normalization combination tests/with_raw_data/all_results_of_all_normalization_combinations_with_curve_data.rds")
  #ranked_results <- readRDS("D:/tmp_PT/Experiments/Coffein-Samplings_SimonBAProject/results/normalization combination tests/top_20_normalization_results_with_curve_data.rds")
  final_pqn_results <- readRDS("D:/tmp_PT/Experiments/Coffein-Samplings_SimonBAProject/results/normalization combination tests/with_pqn_normalized_data/all_results_of_all_normalization_combinations_with_curve_data.rds")
  final_vsn_results <- readRDS("D:/tmp_PT/Experiments/Coffein-Samplings_SimonBAProject/results/normalization combination tests/with_vsn_normalized_data/all_results_of_all_normalization_combinations_with_curve_data.rds")
  final_results <- readRDS("D:/tmp_PT/Experiments/FiS-45min-Abstandsmessung/results/normalization combination tests (followup after simon data)/second try - vsn data - after further filtering of Simon combinations/all_results_of_all_normalization_combinations_with_curve_data.rds")
  
  
  top1000_raw <- rank_top1000_combinations(final_raw_results)
  top1000_pqn <- rank_top1000_combinations(final_pqn_results)
  top1000_vsn <- rank_top1000_combinations(final_vsn_results)
  
  plot_histogram_of_residuals_of_different_normalization_strategies(top1000_raw, top1000_pqn, top1000_vsn, resultsdir, molecule)
  
  
  
  # Test different normalization features with objective methods
  plot_to_pdf_normalization_check(resultsdir, experiment_data_plots, unique_molecules, normalization_molecules)
  plot_to_pdf_normalization_check_per_person_caffein(resultsdir, experiment_data_plots, unique_molecules, normalization_molecules)
  
  plot_to_pdf_normalization_check_objective_tests(resultsdir, norm_feature_list, all_needed_features, add_info, Group = "paper")
  plot_to_pdf_normalization_check_objective_tests(resultsdir, norm_feature_list, all_needed_features, add_info, Group = "Donor")
  plot_to_pdf_normalization_check_objective_tests(resultsdir, norm_feature_list, all_needed_features, add_info, Group = "Intervention")
  
  plot_to_pdf_normalization_check_objective_tests(resultsdir, norm_feature_list, targeted_experiment_data, add_info, Group = "Donor")
  plot_to_pdf_normalization_check_objective_tests(resultsdir, norm_feature_list, targeted_experiment_data, add_info, Group = "Intervention")
  
  
  # Test different normalization methods with objective methods
  # Normalization methods
  norm_methods <- c("Tyrosine", "Total_Sum", "non-weighted", "Sqrt", "Rank2", "Median", "Med.Std.dev", "Mean")
  
  plot_to_pdf_normalization_check_objective_tests(resultsdir, norm_feature_list, all_needed_features, add_info, Group = "paper", norm_methods = norm_methods)
  plot_to_pdf_normalization_check_objective_tests(resultsdir, norm_feature_list, all_needed_features, add_info, Group = "Donor", norm_methods = norm_methods)
  plot_to_pdf_normalization_check_objective_tests(resultsdir, norm_feature_list, all_needed_features, add_info, Group = "Intervention", norm_methods = norm_methods)
  
  plot_to_pdf_normalization_check_objective_tests(resultsdir, norm_feature_list, targeted_experiment_data, add_info, Group = "Donor", norm_methods = norm_methods)
  plot_to_pdf_normalization_check_objective_tests(resultsdir, norm_feature_list, targeted_experiment_data, add_info, Group = "Intervention", norm_methods = norm_methods)
  
  
  # Test different normalization features by showing the goodness of fit to kinetik of biological molecule like caffeine
  plot_to_pdf_various_normalization_fits_for_biological_kinetics(targeted_experiment_data, norm_feature_list, resultsdir, add_info_file, molecule = "Caffein")
  
  
  
  
  
  
# calculate and evaluate CV to check for possible Housekeeping Molecules
  # Calculate the coefficient of variation (CV) for each molecule within each paper
  consistently_stable_molecules <- evaluate_stability_of_molecule_paper_list(experiment_data_plots)#, stable_threshold = 0.25, min_sample_count = 10)
  # Outputs the names and average CV of molecules considered consistently stable
  print(consistently_stable_molecules)
  # Extract only the names of stable molecules
  stable_molecule_names_across_papers <- consistently_stable_molecules$Molecule.Name
  
  
  
# check for correlations in the untargeted dataset
  #identify_correlating_features(all_data, df_ids, datadir, resultsdir, threshold = 0.90,at_least_n_samples = 100, 
  #                              p_value_threshold = 0.05, filename = "corr_features_t0.9_p0.05.csv")
  
  identify_correlating_features_fixed_thresholds(datadir, resultsdir)
  


  
  
# Check change in CV for normalization (PQN, VSN and Tyrosine)
  normalization_list  = paste0(wd, "normalization/Tyrosine-only.csv")
  significant_abundant_features_object <- get_significant_abundant_features(datadir, exp_data, info_file_dir, resultsdir, figuredir, exp, tlists)
  raw_data <- significant_abundant_features_object$significant_abundant_features[,-c(1:9)]
  se <- significant_abundant_features_object$se
  df_f <- as.data.frame(rowData(se))
  pqn_data <- as.data.frame(assay(se, i = "log2_pqn"))
  vsn_data <- as.data.frame(assay(se, i = "vsn"))
  pqn_data <- 2**(pqn_data-20) # make from log2 to "normal" again
  vsn_data <- 2**vsn_data # make from log2 to "normal" again
  PQN.Tyr.norm.data <- normalize_data(pqn_data, normalization_list, df_f)[,-c(1:9)]
  VSN.Tyr.norm.data <- normalize_data(vsn_data, normalization_list, df_f)[,-c(1:9)]
  
  significant_abundant_features_object <- get_significant_abundant_features(datadir, exp_data, info_file_dir, resultsdir, figuredir, exp, tlists, normalization_list)
  all_needed_features <- significant_abundant_features_object$significant_abundant_features
  se <- significant_abundant_features_object$se
  Tyr.PQN.norm.data <- as.data.frame(assay(se, i = "log2_pqn"))
  Tyr.VSN.norm.data <- as.data.frame(assay(se, i = "vsn"))
  Tyr.PQN.norm.data <- 2**(Tyr.PQN.norm.data-20) # make from log2 to "normal" again
  Tyr.VSN.norm.data <- 2**Tyr.VSN.norm.data # make from log2 to "normal" again
  Tyr.norm.data <- significant_abundant_features_object$significant_abundant_features[,-c(1:9)]

  
  raw.data.cv <- get_CV_overview_of_dataset(raw_data) %>% mutate(Typ = 'raw')
  pqn.data.cv <- get_CV_overview_of_dataset(pqn_data) %>% mutate(Typ = 'PQN')
  vsn.data.cv <- get_CV_overview_of_dataset(vsn_data) %>% mutate(Typ = 'VSN')
  pqn.tyr.data.cv <- get_CV_overview_of_dataset(PQN.Tyr.norm.data) %>% mutate(Typ = 'PQN+Tyr')
  vsn.tyr.data.cv <- get_CV_overview_of_dataset(VSN.Tyr.norm.data) %>% mutate(Typ = 'VSN+Tyr')
  tyr.data.cv <- get_CV_overview_of_dataset(Tyr.norm.data) %>% mutate(Typ = 'Tyr')
  tyr.pqn.data.cv <- get_CV_overview_of_dataset(Tyr.PQN.norm.data) %>% mutate(Typ = 'Tyr+PQN')
  tyr.vsn.data.cv <- get_CV_overview_of_dataset(Tyr.VSN.norm.data) %>% mutate(Typ = 'Tyr+VSN')
  full.norm.plot.data <- rbind(raw.data.cv, pqn.data.cv, vsn.data.cv, pqn.tyr.data.cv, vsn.tyr.data.cv,
                              tyr.data.cv, tyr.pqn.data.cv, tyr.vsn.data.cv)
  
  lowest_10_percent_features <- rownames(raw.data.cv[raw.data.cv$cv <= quantile(raw.data.cv$cv, 0.1, na.rm = TRUE), ])
  
  ggplot(full.norm.plot.data, aes(x = factor(Typ), y = cv, fill = factor(Typ))) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +  # Suppress outliers in boxplot for clarity
    geom_jitter(aes(color = Feature %in% lowest_10_percent_features), 
                width = 0.4, size = 1, alpha = 0.3) +  # Conditional color for points
    labs(title = "CV Distribution by Normalization Combination", 
         x = 'Normalization Typ', 
         y = "Coefficient of Variation (CV)", 
         fill = 'Typ', 
         color = "Highlighted Features") +  # Add legend for conditional coloring
    theme_minimal() +
    scale_color_manual(values = c("gray1", "red"), 
                       labels = c("Other", "Lowest 10%"))  # Custom legend labels
  ggsave(paste0(resultsdir, 'CV-comparison-of-norm-methods.png'))
  
  
  
  
  
# get PQN and VSN factors
  data_matrix <- t(as.matrix(all_needed_features[-c(1:9)]))
  reference_spectrum <- apply(data_matrix, 2, median, na.rm = TRUE) # Calculate the reference spectrum (median of all samples)
  pqn_factors <- apply(data_matrix, 1, function(sample) { # Compute normalization factors for each sample
    median(sample / reference_spectrum, na.rm = TRUE)
  })
  pqn_normalization_factors <- t(as.data.frame(pqn_factors))
  
  
  # Ensure non-negative values in the data (VSN requires non-negative data)
  data_matrix[data_matrix < 0] <- 0
  # Apply VSN normalization
  vsn_fit <- vsn2(data_matrix, minDataPointsPerStratum = 5)
  vsn_normalized_data <- predict(vsn_fit, data_matrix)  # Normalize the data
  vsn_normalization_factors <- rowMeans(vsn_normalized_data)   # compute normalization factors (e.g., mean or scaling for each sample)

  pqn_vsn_none <- rbind(pqn_normalization_factors, vsn_normalization_factors)

  
# put several experiments together
  add_info_file <- paste0(wd, exp, info_file)
  add_info <- read.csv(add_info_file)
  if (length(add_info) < 2) {
    add_info <- read.csv2(add_info_file)
  }
  
  norm_feature_list <- norm_feature_list %>% distinct(id, .keep_all = TRUE)
  #targeted_experiment_data <- targeted_experiment_data[-c(2,4,6,8),]
  full_prep_data <- prepare_data_for_plot(targeted_experiment_data, add_info_file)
  
  # first experiment
  all_full_prep_data <- full_prep_data
  all_norm_feature_list <- norm_feature_list
  all_targeted_experiment_data <- targeted_experiment_data
  all_pqn_vsn_list <- pqn_vsn_none
  
  # additional experiments
  all_full_prep_data <- rbind(all_full_prep_data[c("Sample","id","rt","mz","charge","Molecule.Name","Area",'Donor',"Timepoint","time")], full_prep_data[c("Sample","id","rt","mz","charge","Molecule.Name","Area",'Donor',"Timepoint","time")])
  all_norm_feature_list <- cbind(all_norm_feature_list[, -which(names(all_norm_feature_list) == "Molecule.Name")], norm_feature_list[-6, -c(1:4)])
  all_targeted_experiment_data <- cbind(all_targeted_experiment_data[, -which(names(all_targeted_experiment_data) == "Molecule.Name")], targeted_experiment_data[, -c(1:4)])
  all_pqn_vsn_list <- cbind(all_pqn_vsn_list, pqn_vsn_none)
  
  # save dataframes
  write.csv(all_targeted_experiment_data, paste0(resultsdir, 'all_targeted_experiment_data_filtered.csv'))
  write.csv(all_norm_feature_list, paste0(resultsdir, 'all_norm_feature_list_filtered.csv'))
  write.csv(all_pqn_vsn_list, paste0(resultsdir, 'all_pqn_vsn_list_filtered.csv'), row.names = TRUE)
  write.csv(all_full_prep_data, paste0(resultsdir, 'all_full_prep_data_filtered.csv'))
  
  # remove columns/samples that were removed in the filtering of the data
  all_targeted_experiment_data <- all_targeted_experiment_data %>%
    select(c(1:4, intersect(colnames(all_pqn_vsn_list), colnames(all_targeted_experiment_data)[-c(1:4, ncol(all_targeted_experiment_data))])), ncol(all_targeted_experiment_data))
  all_norm_feature_list <- all_norm_feature_list %>% 
    select(c(1:4, intersect(colnames(all_pqn_vsn_list), colnames(all_norm_feature_list)[-c(1:4, ncol(all_norm_feature_list))])), ncol(all_norm_feature_list))
  all_full_prep_data <- all_full_prep_data %>% 
    filter(Sample %in% colnames(all_pqn_vsn_list))
  all_pqn_vsn_list <- as.data.frame(all_pqn_vsn_list) %>%
    select(intersect(colnames(all_targeted_experiment_data)[-c(1:4, ncol(all_targeted_experiment_data))], colnames(all_pqn_vsn_list)))
  