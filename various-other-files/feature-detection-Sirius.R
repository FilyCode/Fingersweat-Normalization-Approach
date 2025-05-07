source("Y:/Studenten_Schueler/Philipp_Trollmann/R Project/MS-data-analysis_functions.R")

# set all the experiment and file names
  wd = "Y:/Studenten_Schueler/Philipp_Trollmann/Experiments/"
  setwd(wd)
  
  exp = "FiS-45min-Abstandsmessung"
  datadir = paste0(wd, exp, "/data/")
  resultsdir = paste0(wd, exp, "/results/")
  tlists = paste0(wd, "Transition-Lists/")




# read in SIRIUS table, reduce df to important columns
# parse "name" column for "null" entry and replace with "InChIkey2D" entry, preceded by "InChIKey:"
  SIRIUS_raw <- read.delim(paste0(datadir, "/SIRIUS/compound_identifications.tsv"), header = TRUE)
  SIRIUS_raw <- SIRIUS_raw %>%
    mutate(pubchemids = str_split(pubchemids, ";") %>% sapply(`[`, 1))
  KeepColumns <- c("id", "featureId", "molecularFormula", "name", "InChIkey2D", "ConfidenceScore", "ionMass", "retentionTimeInSeconds", "smiles", "pubchemids")
  SIRIUS <- SIRIUS_raw[,KeepColumns]
  for(m in 1:dim(SIRIUS)[1]){
    if(SIRIUS$name[m] == "null"){
      SIRIUS$name[m] <- paste0("InChIKey: ", SIRIUS$InChIkey2D[m])
    }
  }

# reordering df for export
  neworder <- c("featureId", "name", "retentionTimeInSeconds", "ionMass", "molecularFormula", "ConfidenceScore", "smiles", "pubchemids")
  newnames <- c("id", "Annotation", "RT", "MZ", "Formula", "Confidence", "SMILES", "CID")
  SIRIUS <- SIRIUS[,neworder]
  colnames(SIRIUS) <- newnames
  SIRIUS$RT <- round(SIRIUS$RT / 60,3)
  SIRIUS$Confidence <- as.numeric(SIRIUS$Confidence)

# webchem package for better annotation
  pc.properties <- c("Title", "XLogP", "IUPACName")
  pc.info <- pc_prop(SIRIUS$CID, properties = pc.properties) %>%
    mutate(Title = ifelse(is.na(Title), IUPACName, Title)) %>%
    select(-IUPACName) %>%
    dplyr::rename(Name = Title)
  
  SIRIUS <- SIRIUS %>%
    mutate(CID = as.numeric(CID)) %>%
    left_join(pc.info, by = "CID", multiple = "first") %>%
    mutate(Annotation = Name) %>%
    select(-Name)

# Expressions Sets require Feature data row names to be of equal identity and length so I'm adding the missing feature from the original assay data
# Get the whole processed data set from file (already prepared csv file with mzMine)
  exp_data <- process_merged_files(datadir)
  all_data <- exp_data$data
  df_ids <- exp_data$df_ids
  
  df_data <- all_data

# Feature data from Sirius merged with data ids
  df_features <- merge(df_ids, SIRIUS, by = "id", all.x =TRUE)


# Creating meta data straight from column names
  df_phenodata <- data.frame(ID = colnames(df_data)) %>%
                            separate_wider_delim(cols = ID, delim = "_", names = c("Experiment", "Identifier", "Timepoint"),
                                                 cols_remove = FALSE, too_few = "align_start") %>%
                            select(-Experiment) %>%
                            mutate(Intervention = case_when(
                              Timepoint == 't0' ~ "pre_Intervention",
                              Timepoint != 't0' ~ "post_Intervention"
                            )) %>%
                            mutate(Intervention = factor(Intervention, levels = c("pre_Intervention", "post_Intervention")),
                                   Identifier = factor(Identifier, levels = unique(Identifier))) %>%
                            column_to_rownames("ID")



# SummarizedExperiment: Assembly
  se <- SummarizedExperiment(assays = list("raw" = df_data),
                             colData = df_phenodata,
                             rowData = df_features)


# Feature Reduction, Imputation, Normalization and Clustering
  se <- removeFeatures(se, i = "raw", group = se$Intervention, cut = 0.80)

  assay(se, i = "norm_log2") <- assay(se, i = "raw") %>%
    mutate_all(~log2(.)+20)
  se <- imputeIntensity(se, i = "norm_log2", name = "norm_MinProb", method = "MinProb")
  se <- normalizeIntensity(se, i = "norm_MinProb", name = "norm_MinProb_pqn", method = "pqn")
  se <- clusterFeatures(se,i = "norm_MinProb_pqn", rtime_var = "rt", rt_cut = 0.003)

# get significant features
  feature_groups <- as.data.frame(rowData(se)[,c("id", "rt", "mz", "rtime_group", "feature_group", "Annotation", "Confidence", "SMILES")])
  significant_features <- feature_groups[!is.na(feature_groups$Confidence) & feature_groups$Confidence > 0.25, ]
  
  columns_to_append <- df_data[significant_features$id, ]
  all_significant_feature_data <- cbind(significant_features, columns_to_append)
  data_length <- ncol(all_significant_feature_data)
  all_significant_feature_data[, 9:data_length][is.na(all_significant_feature_data[, 9:data_length])] <- 1 #make Na 1
  num_cols_abundant <- rowSums(all_significant_feature_data[, 9:data_length] > 1e5, na.rm = TRUE)
  significant_abundant_features <- all_significant_feature_data[num_cols_abundant >= 2, ]
  
  mean_abundance_each_col <- colMeans(significant_abundant_features[, 9:data_length], na.rm = TRUE)





# rank the features by their abundance, plot each sample and safe to file
  pdf(file = paste0(resultsdir, "feature_ranking_each-sample.pdf"), width = 150, height = 80) 
  
  # plot raw abundance value
    plot_list <- list()
    data_length <- ncol(significant_abundant_features)
    
    
    # Loop over columns 9 to data_length (all samples from study)
    for (col_index in seq(9, data_length)) {
      
      df_rank <- data.frame(name = significant_abundant_features$Annotation, sample_area = significant_abundant_features[, col_index])
      df_rank <- na.omit(df_rank)
      
      # Rank the rows based on the values in the current column
      df_rank$rank <- rank(-df_rank$sample_area, na.last = NA,  ties.method = 'first')  # Use -df_rank[, col_index] for descending order
      
      # Identify the top 20 rows
      top20 <- head(df_rank[order(-df_rank$sample_area), ], 20)
      
      # Create the plot for the current column
      p <- ggplot(df_rank, aes(x = rank, y = sample_area)) +
        geom_point() +
        geom_text_repel(data = top20, aes(label = name, color = "Top20"), max.overlaps = Inf, min.segment.length = 0.01) + 
        labs(x = "Rank", y = "Abundance",
             title = paste("Feature Ranking of", colnames(significant_abundant_features[col_index])))
      
      # Add the plot to the list
      plot_list[[col_index - 8]] <- p
    }
    
    # Arrange the plots in a grid
    do.call(grid.arrange, plot_list)
  
  
  # plot log2 abundance value
    plot_list <- list()
    
    # make raw abundance to log abundance
    significant_abundant_features_log <- significant_abundant_features[, 1:data_length]
    significant_abundant_features_log[, 9:data_length] <- apply(significant_abundant_features_log[, 9:data_length], 2, log2)
    
    
    # Loop over columns 9 to data_length (all samples from study)
    for (col_index in seq(9, data_length)) {
      
      df_rank <- data.frame(name = significant_abundant_features_log$Annotation, sample_area = significant_abundant_features_log[, col_index])
      df_rank <- na.omit(df_rank)
      
      # Rank the rows based on the values in the current column
      df_rank$rank <- rank(-df_rank$sample_area, na.last = NA,  ties.method = 'first')  # Use -df_rank[, col_index] for descending order
      
      # Identify the top 20 rows
      top20 <- head(df_rank[order(-df_rank$sample_area), ], 20)
      
      # Create the plot for the current column
      p <- ggplot(df_rank, aes(x = rank, y = sample_area)) +
        geom_point() +
        geom_text_repel(data = top20, aes(label = name, color = "Top20"), max.overlaps = Inf, min.segment.length = 0.01) + 
        labs(x = "Rank", y = "Abundance",
             title = paste("Feature Ranking of", colnames(significant_abundant_features_log[col_index])))
      
      # Add the plot to the list
      plot_list[[col_index - 8]] <- p
    }
    
    # Arrange the plots in a grid
    do.call(grid.arrange, plot_list)
  
  dev.off() # Close the PDF device


  

  


  
  
  
# normalize the areas of each column with the sum of the according column
  # make dataframes for raw and normalized feature data with all data from both papers
    raw_sign_abund_feature <- significant_abundant_features
    norm_sign_abund_feature <- raw_sign_abund_feature
    data_length <- ncol(raw_sign_abund_feature)
    column_sums <- colSums(norm_sign_abund_feature[, 9:data_length], na.rm = TRUE)
  
  # make dataframes for raw and normalized feates data from only Chrom OR Kim
    # Chrom
      raw_sign_abund_feature_chrom <- raw_sign_abund_feature
      info_columns <- raw_sign_abund_feature[, 1:8]
      raw_sign_abund_feature_chrom <- cbind(info_column, raw_sign_abund_feature_chrom[, grep("Chrom", colnames(raw_sign_abund_feature_chrom))])
      norm_sign_abund_feature_chrom <- raw_sign_abund_feature_chrom
      data_length_chrom <- ncol(raw_sign_abund_feature_chrom)
      column_sums_chrom <- colSums(norm_sign_abund_feature_chrom[, 9:data_length_chrom], na.rm = TRUE)
    
    # Kim
      raw_sign_abund_feature_kim <- raw_sign_abund_feature
      info_columns <- raw_sign_abund_feature[, 1:8]
      raw_sign_abund_feature_kim <- cbind(info_column, raw_sign_abund_feature_kim[, grep("Kim", colnames(raw_sign_abund_feature_kim))])
      norm_sign_abund_feature_kim <- raw_sign_abund_feature_kim
      data_length_kim <- ncol(raw_sign_abund_feature_kim)
      column_sums_kim <- colSums(norm_sign_abund_feature_kim[, 9:data_length_kim], na.rm = TRUE)
    
  
  # normalization step
    # for data with bot papers
      for (col_index in seq(9,data_length)) {
        col_name <- colnames(norm_sign_abund_feature)[col_index]
        # get percentage of abundance per feature
        norm_sign_abund_feature[col_index] <- (norm_sign_abund_feature[[col_name]] / column_sums[col_index-8]) * 100 
        colnames(norm_sign_abund_feature)[col_index] <- paste0(col_name, "_norm")
      }
    
    # for data from Chrom paper
      for (col_index in seq(9,data_length_chrom)) {
        col_name <- colnames(norm_sign_abund_feature_chrom)[col_index]
        # get percentage of abundance per feature
        norm_sign_abund_feature_chrom[col_index] <- (norm_sign_abund_feature_chrom[[col_name]] / column_sums_chrom[col_index-8]) * 100 
        colnames(norm_sign_abund_feature_chrom)[col_index] <- paste0(col_name, "_norm")
      }
    
    # for data from Kim paper
      for (col_index in seq(9,data_length_kim)) {
        col_name <- colnames(norm_sign_abund_feature_kim)[col_index]
        # get percentage of abundance per feature
        norm_sign_abund_feature_kim[col_index] <- (norm_sign_abund_feature_kim[[col_name]] / column_sums_kim[col_index-8]) * 100 
        colnames(norm_sign_abund_feature_kim)[col_index] <- paste0(col_name, "_norm")
      }

# rank the features by their abundance, plot the mean and std. dist. over all samples (raw and normalized) and safe to file
  # get mean, std dist, rank for whole data with both papers
    # mean, std dist, rank for raw data
      raw_sign_abund_feature$mean <- rowMeans(raw_sign_abund_feature[, 9:data_length], na.rm = TRUE)
      raw_sign_abund_feature$std_dev <- apply(raw_sign_abund_feature[, 9:data_length], 1, sd, na.rm = TRUE)
      raw_sign_abund_feature$rel_std_dev <- (raw_sign_abund_feature$std_dev / raw_sign_abund_feature$mean) * 100
      raw_sign_abund_feature$mean_rank <- rank(-raw_sign_abund_feature$mean)
      mean_top_20 <- head(raw_sign_abund_feature[order(-raw_sign_abund_feature$mean), ], 20)
    # apply log2 and get mean, std dist, rank for raw data
      raw_sign_abund_feature_log <- raw_sign_abund_feature
      raw_sign_abund_feature_log[, 9:data_length] <- apply(raw_sign_abund_feature_log[, 9:data_length], 2, log2)
      colnames(raw_sign_abund_feature_log)[9:data_length] <- paste0(colnames(raw_sign_abund_feature_log)[9:data_length], "_log")
      raw_sign_abund_feature_log$mean_log <- rowMeans(raw_sign_abund_feature_log[, 9:data_length], na.rm = TRUE)
      raw_sign_abund_feature_log$std_dev_log <- apply(raw_sign_abund_feature_log[, 9:data_length], 1, sd, na.rm = TRUE)
      raw_sign_abund_feature_log$rel_std_dev_log <- (raw_sign_abund_feature_log$std_dev_log / raw_sign_abund_feature_log$mean_log) * 100
      raw_sign_abund_feature_log$mean_log_rank <- rank(-raw_sign_abund_feature_log$mean_log)
      mean_log_top_20 <- head(raw_sign_abund_feature_log[order(-raw_sign_abund_feature_log$mean_log), ], 20)
      
    # mean, std dist, rank for normalized data
      norm_sign_abund_feature$mean_norm <- rowMeans(norm_sign_abund_feature[, 9:data_length], na.rm = TRUE)
      norm_sign_abund_feature$std_dev_norm <- apply(norm_sign_abund_feature[, 9:data_length], 1, sd, na.rm = TRUE)
      norm_sign_abund_feature$rel_std_dev_norm <- (norm_sign_abund_feature$std_dev / norm_sign_abund_feature$mean) * 100
      norm_sign_abund_feature$mean_norm_rank <- rank(-norm_sign_abund_feature$mean_norm)
      norm_mean_top_20 <- head(norm_sign_abund_feature[order(-norm_sign_abund_feature$mean_norm), ], 20)
    # apply log2 and get mean, std dist, rank for normalized data
      norm_sign_abund_feature_log <- norm_sign_abund_feature
      norm_sign_abund_feature_log[, 9:data_length] <- apply(norm_sign_abund_feature_log[, 9:data_length], 2, log2)
      colnames(norm_sign_abund_feature_log)[9:data_length] <- paste0(colnames(norm_sign_abund_feature_log)[9:data_length], "_log")
      norm_sign_abund_feature_log$mean_norm_log <- rowMeans(norm_sign_abund_feature_log[, 9:data_length], na.rm = TRUE)
      norm_sign_abund_feature_log$std_dev_norm_log <- apply(norm_sign_abund_feature_log[, 9:data_length], 1, sd, na.rm = TRUE)
      norm_sign_abund_feature_log$rel_std_dev_norm_log <- (norm_sign_abund_feature_log$std_dev_norm_log / norm_sign_abund_feature_log$mean_norm_log) * 100
      norm_sign_abund_feature_log$mean_norm_log_rank <- rank(-norm_sign_abund_feature_log$mean_norm_log)
      norm_mean_log_top_20 <- head(norm_sign_abund_feature_log[order(-norm_sign_abund_feature_log$mean_norm_log), ], 20)
    
    
  # get mean, std dist, rank for data from Chrom paper
    # mean, std dist, rank for raw data
      raw_sign_abund_feature_chrom$mean <- rowMeans(raw_sign_abund_feature_chrom[, 9:data_length_chrom], na.rm = TRUE)
      raw_sign_abund_feature_chrom$std_dev <- apply(raw_sign_abund_feature_chrom[, 9:data_length_chrom], 1, sd, na.rm = TRUE)
      raw_sign_abund_feature_chrom$rel_std_dev <- (raw_sign_abund_feature_chrom$std_dev / raw_sign_abund_feature_chrom$mean) * 100
      raw_sign_abund_feature_chrom$mean_rank <- rank(-raw_sign_abund_feature_chrom$mean)
      mean_top_20_chrom <- head(raw_sign_abund_feature_chrom[order(-raw_sign_abund_feature_chrom$mean), ], 20)
    # apply log2 and get mean, std dist, rank for raw data
      raw_sign_abund_feature_chrom_log <- raw_sign_abund_feature_chrom
      raw_sign_abund_feature_chrom_log[, 9:data_length_chrom] <- apply(raw_sign_abund_feature_chrom_log[, 9:data_length_chrom], 2, log2)
      colnames(raw_sign_abund_feature_chrom_log)[9:data_length_chrom] <- paste0(colnames(raw_sign_abund_feature_chrom_log)[9:data_length_chrom], "_log")
      raw_sign_abund_feature_chrom_log$mean_log <- rowMeans(raw_sign_abund_feature_chrom_log[, 9:data_length_chrom], na.rm = TRUE)
      raw_sign_abund_feature_chrom_log$std_dev_log <- apply(raw_sign_abund_feature_chrom_log[, 9:data_length_chrom], 1, sd, na.rm = TRUE)
      raw_sign_abund_feature_chrom_log$rel_std_dev_log <- (raw_sign_abund_feature_chrom_log$std_dev_log / raw_sign_abund_feature_chrom_log$mean_log) * 100
      raw_sign_abund_feature_chrom_log$mean_log_rank <- rank(-raw_sign_abund_feature_chrom_log$mean_log)
      mean_log_top_20_chrom <- head(raw_sign_abund_feature_chrom_log[order(-raw_sign_abund_feature_chrom_log$mean_log), ], 20)
      mean_log_top_10_chrom <- head(raw_sign_abund_feature_chrom_log[order(-raw_sign_abund_feature_chrom_log$mean_log), ], 10)
      
    # mean, std dist, rank for normalized data
      norm_sign_abund_feature_chrom$mean_norm <- rowMeans(norm_sign_abund_feature_chrom[, 9:data_length_chrom], na.rm = TRUE)
      norm_sign_abund_feature_chrom$std_dev_norm <- apply(norm_sign_abund_feature_chrom[, 9:data_length_chrom], 1, sd, na.rm = TRUE)
      norm_sign_abund_feature_chrom$rel_std_dev_norm <- (norm_sign_abund_feature_chrom$std_dev / norm_sign_abund_feature_chrom$mean) * 100
      norm_sign_abund_feature_chrom$mean_norm_rank <- rank(-norm_sign_abund_feature_chrom$mean_norm)
      norm_mean_top_20_chrom <- head(norm_sign_abund_feature_chrom[order(-norm_sign_abund_feature_chrom$mean_norm), ], 20)
    # apply log2 and get mean, std dist, rank for normalized data
      norm_sign_abund_feature_chrom_log <- norm_sign_abund_feature_chrom
      norm_sign_abund_feature_chrom_log[, 9:data_length_chrom] <- apply(norm_sign_abund_feature_chrom_log[, 9:data_length_chrom], 2, log2)
      colnames(norm_sign_abund_feature_chrom_log)[9:data_length_chrom] <- paste0(colnames(norm_sign_abund_feature_chrom_log)[9:data_length_chrom], "_log")
      norm_sign_abund_feature_chrom_log$mean_norm_log <- rowMeans(norm_sign_abund_feature_chrom_log[, 9:data_length_chrom], na.rm = TRUE)
      norm_sign_abund_feature_chrom_log$std_dev_norm_log <- apply(norm_sign_abund_feature_chrom_log[, 9:data_length_chrom], 1, sd, na.rm = TRUE)
      norm_sign_abund_feature_chrom_log$rel_std_dev_norm_log <- (norm_sign_abund_feature_chrom_log$std_dev_norm_log / norm_sign_abund_feature_chrom_log$mean_norm_log) * 100
      norm_sign_abund_feature_chrom_log$mean_norm_log_rank <- rank(-norm_sign_abund_feature_chrom_log$mean_norm_log)
      norm_mean_log_top_20_chrom <- head(norm_sign_abund_feature_chrom_log[order(-norm_sign_abund_feature_chrom_log$mean_norm_log), ], 20)
      norm_mean_log_top_10_chrom <- head(norm_sign_abund_feature_chrom_log[order(-norm_sign_abund_feature_chrom_log$mean_norm_log), ], 10)
  
      
  # get mean, std dist, rank for data from Kim paper
    # mean, std dist, rank for raw data
      raw_sign_abund_feature_kim$mean <- rowMeans(raw_sign_abund_feature_kim[, 9:data_length_kim], na.rm = TRUE)
      raw_sign_abund_feature_kim$std_dev <- apply(raw_sign_abund_feature_kim[, 9:data_length_kim], 1, sd, na.rm = TRUE)
      raw_sign_abund_feature_kim$rel_std_dev <- (raw_sign_abund_feature_kim$std_dev / raw_sign_abund_feature_kim$mean) * 100
      raw_sign_abund_feature_kim$mean_rank <- rank(-raw_sign_abund_feature_kim$mean)
      mean_top_20_kim <- head(raw_sign_abund_feature_kim[order(-raw_sign_abund_feature_kim$mean), ], 20)
    # apply log2 and get mean, std dist, rank for raw data
      raw_sign_abund_feature_kim_log <- raw_sign_abund_feature_kim
      raw_sign_abund_feature_kim_log[, 9:data_length_kim] <- apply(raw_sign_abund_feature_kim_log[, 9:data_length_kim], 2, log2)
      colnames(raw_sign_abund_feature_kim_log)[9:data_length_kim] <- paste0(colnames(raw_sign_abund_feature_kim_log)[9:data_length_kim], "_log")
      raw_sign_abund_feature_kim_log$mean_log <- rowMeans(raw_sign_abund_feature_kim_log[, 9:data_length_kim], na.rm = TRUE)
      raw_sign_abund_feature_kim_log$std_dev_log <- apply(raw_sign_abund_feature_kim_log[, 9:data_length_kim], 1, sd, na.rm = TRUE)
      raw_sign_abund_feature_kim_log$rel_std_dev_log <- (raw_sign_abund_feature_kim_log$std_dev_log / raw_sign_abund_feature_kim_log$mean_log) * 100
      raw_sign_abund_feature_kim_log$mean_log_rank <- rank(-raw_sign_abund_feature_kim_log$mean_log)
      mean_log_top_20_kim <- head(raw_sign_abund_feature_kim_log[order(-raw_sign_abund_feature_kim_log$mean_log), ], 20)
      mean_log_top_10_kim <- head(raw_sign_abund_feature_kim_log[order(-raw_sign_abund_feature_kim_log$mean_log), ], 10)
      
    # mean, std dist, rank for normalized data
      norm_sign_abund_feature_kim$mean_norm <- rowMeans(norm_sign_abund_feature_kim[, 9:data_length_kim], na.rm = TRUE)
      norm_sign_abund_feature_kim$std_dev_norm <- apply(norm_sign_abund_feature_kim[, 9:data_length_kim], 1, sd, na.rm = TRUE)
      norm_sign_abund_feature_kim$rel_std_dev_norm <- (norm_sign_abund_feature_kim$std_dev / norm_sign_abund_feature_kim$mean) * 100
      norm_sign_abund_feature_kim$mean_norm_rank <- rank(-norm_sign_abund_feature_kim$mean_norm)
      norm_mean_top_20_kim <- head(norm_sign_abund_feature_kim[order(-norm_sign_abund_feature_kim$mean_norm), ], 20)
    # apply log2 and get mean, std dist, rank for normalized data
      norm_sign_abund_feature_kim_log <- norm_sign_abund_feature_kim
      norm_sign_abund_feature_kim_log[, 9:data_length_kim] <- apply(norm_sign_abund_feature_kim_log[, 9:data_length_kim], 2, log2)
      colnames(norm_sign_abund_feature_kim_log)[9:data_length_kim] <- paste0(colnames(norm_sign_abund_feature_kim_log)[9:data_length_kim], "_log")
      norm_sign_abund_feature_kim_log$mean_norm_log <- rowMeans(norm_sign_abund_feature_kim_log[, 9:data_length_kim], na.rm = TRUE)
      norm_sign_abund_feature_kim_log$std_dev_norm_log <- apply(norm_sign_abund_feature_kim_log[, 9:data_length_kim], 1, sd, na.rm = TRUE)
      norm_sign_abund_feature_kim_log$rel_std_dev_norm_log <- (norm_sign_abund_feature_kim_log$std_dev_norm_log / norm_sign_abund_feature_kim_log$mean_norm_log) * 100
      norm_sign_abund_feature_kim_log$mean_norm_log_rank <- rank(-norm_sign_abund_feature_kim_log$mean_norm_log)
      norm_mean_log_top_20_kim <- head(norm_sign_abund_feature_kim_log[order(-norm_sign_abund_feature_kim_log$mean_norm_log), ], 20)  
      norm_mean_log_top_10_kim <- head(norm_sign_abund_feature_kim_log[order(-norm_sign_abund_feature_kim_log$mean_norm_log), ], 10)  
  
  
  # dateframe to make paper comparison
    # Combine columns for Chrom data
    chrom_combined <- cbind(
      raw_sign_abund_feature_chrom_log[, c("mean_log", "std_dev_log", "mean_log_rank")],
      norm_sign_abund_feature_chrom_log[, c("mean_norm_log", "std_dev_norm_log", "mean_norm_log_rank")]
    )
    # Add postfix "_chrom" to column names from Chrom data except for "id"
    names(chrom_combined) <- gsub("$", "_chrom", names(chrom_combined))
    # Add id and Annotation
    chrom_combined <- cbind(raw_sign_abund_feature_chrom_log[, c("id", "Annotation")], chrom_combined)
    
    # Combine columns for Kim data
    kim_combined <- cbind(
      raw_sign_abund_feature_kim_log[, c("mean_log", "std_dev_log", "mean_log_rank")],
      norm_sign_abund_feature_kim_log[, c("mean_norm_log", "std_dev_norm_log", "mean_norm_log_rank")]
    )
    # Add postfix "_kim" to column names from Kim data
    names(kim_combined) <- gsub("$", "_kim", names(kim_combined))
    # Add id
    kim_combined <- cbind(raw_sign_abund_feature_kim_log$id, kim_combined)
    names(kim_combined)[1] <- "id" # need to set column name again as it gets lost with cbind if only 1 column is combined
    
    # Merge Chrom and Kim data based on the common column ID
    paper_comparison_data <- merge(chrom_combined, kim_combined, by = "id", all = TRUE)
      
    
    
  # Features to label additionally
    names_to_label <- c("Paraxanthine", "Theobromine", "Theophylline")
  
  # get the number of features
    num_features <- nrow(raw_sign_abund_feature)
  
  # Calculate the sum of mean_norm for the top 20 and the rest
    sum_mean_norm_top20 <- sum(norm_mean_top_20$mean_norm)
    sum_mean_norm_rest <- sum(norm_sign_abund_feature$mean_norm) - sum_mean_norm_top20
    sum_mean_log_norm_top20 <- sum(norm_mean_log_top_20$mean_norm_log)
    sum_mean_log_norm_rest <- sum(norm_sign_abund_feature_log$mean_norm_log) - sum_mean_log_norm_top20
    sum_mean_norm_top20_chrom <- sum(norm_mean_top_20_chrom$mean_norm)
    sum_mean_norm_rest_chrom <- sum(norm_sign_abund_feature_chrom$mean_norm) - sum_mean_norm_top20_chrom
    sum_mean_log_norm_top20_chrom <- sum(norm_mean_log_top_20_chrom$mean_norm_log)
    sum_mean_log_norm_rest_chrom <- sum(norm_sign_abund_feature_chrom_log$mean_norm_log) - sum_mean_log_norm_top20_chrom
    sum_mean_norm_top20_kim <- sum(norm_mean_top_20_kim$mean_norm)
    sum_mean_norm_rest_kim <- sum(norm_sign_abund_feature_kim$mean_norm) - sum_mean_norm_top20_kim
    sum_mean_log_norm_top20_kim <- sum(norm_mean_log_top_20_kim$mean_norm_log)
    sum_mean_log_norm_rest_kim <- sum(norm_sign_abund_feature_kim_log$mean_norm_log) - sum_mean_log_norm_top20_kim
  
  
  pdf(file = paste0(resultsdir, "feature_ranking_sample-mean.pdf"), width = 15, height = 15) 
  
  ### BOTH PAPERS
    # plot the absolute values of both papers
      plot_list <- list()

      # Plot mean and standard deviation for each row
      # plot raw data
      plot_list[[1]] <- ggplot(raw_sign_abund_feature, aes(x = mean_rank, y = mean)) +
                          geom_ribbon(aes(ymin = mean - std_dev, ymax = mean + std_dev, fill = "Standard Deviation"), alpha = 0.4) +  # Add ribbon for standard deviation
                          geom_point(aes(color = "Mean"), size = 1) +  # Add points for mean with bigger size
                          geom_text_repel(data = mean_top_20, aes(label = Annotation), color = "black", max.overlaps = Inf, min.segment.length = 0.01) +  # Add text labels for top 20 rows
                          geom_text_repel(data = raw_sign_abund_feature[raw_sign_abund_feature$Annotation %in% names_to_label, ],
                                          aes(label = Annotation), color = "darkgrey", max.overlaps = Inf, vjust = 5, min.segment.length = 0.01) +  # Add text labels for specific points with line always
                          labs(x = "Mean Rank", y = "raw Abundance", color = NULL, fill = NULL,  # Remove color and fill legends from default mapping
                          title = "Feature Ranking of raw Abundance") +  # Add plot title              
                          scale_color_manual(values = c("Mean" = "darkblue"), guide = guide_legend(title = "Mean")) +  # Manual color legend for mean
                          scale_fill_manual(values = "skyblue", guide = guide_legend(title = "Standard Deviation")) +  # Manual fill legend for standard deviation
                          theme(legend.position = "top") +  # Position legend at the top of the plot
                          annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, 
                                   label = paste("Number of Features:", num_features))  # Add number of features as text to top right corner
      
      # plot normalized data
      plot_list[[2]] <- ggplot(norm_sign_abund_feature, aes(x = mean_norm_rank, y = mean_norm)) +
                      geom_ribbon(aes(ymin = mean_norm - std_dev_norm, ymax = mean_norm + std_dev_norm, fill = "Standard Deviation"), alpha = 0.4) +  # Add ribbon for standard deviation
                      geom_point(aes(color = "Mean"), size = 1) +  # Add points for mean with bigger size
                      geom_text_repel(data = norm_mean_top_20, aes(label = Annotation), color = "black", max.overlaps = Inf, min.segment.length = 0.01) +  # Add text labels for top 20 rows
                      geom_text_repel(data = norm_sign_abund_feature[norm_sign_abund_feature$Annotation %in% names_to_label, ],
                                      aes(label = Annotation), color = "darkgrey", max.overlaps = Inf, vjust = 5, min.segment.length = 0.01) +  # Add text labels for specific points with line always
                      labs(x = "Mean Rank", y = "norm Abundance", color = NULL, fill = NULL,  # Remove color and fill legends from default mapping
                      title = "Feature Ranking of normalized Abundance") +  # Add plot title              
                      scale_color_manual(values = c("Mean" = "darkblue"), guide = guide_legend(title = "Mean")) +  # Manual color legend for mean
                      scale_fill_manual(values = "skyblue", guide = guide_legend(title = "Standard Deviation")) +  # Manual fill legend for standard deviation
                      theme(legend.position = "top") +  # Position legend at the top of the plot
                      annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, 
                               label = paste("Number of Features:", num_features, "\n",
                                             "Sum of mean_norm (Top 20):", sum_mean_norm_top20, "\n",
                                             "Sum of mean_norm (Rest):", sum_mean_norm_rest))  # Add text with number of features and sum of mean_norm to top right corner
        
      # Arrange the plots in a grid
      do.call(grid.arrange, plot_list)
      
    
    # plot the log2 values of both papers
      plot_list <- list()
      
      # Plot mean and standard deviation for each row
      # plot raw data
      plot_list[[1]] <- ggplot(raw_sign_abund_feature_log, aes(x = mean_log_rank, y = mean_log)) +
        geom_ribbon(aes(ymin = mean_log - std_dev_log, ymax = mean_log + std_dev_log, fill = "Standard Deviation"), alpha = 0.4) +  # Add ribbon for standard deviation
        geom_point(aes(color = "Mean"), size = 1) +  # Add points for mean with bigger size
        geom_text_repel(data = mean_log_top_20, aes(label = Annotation), color = "black", max.overlaps = Inf, min.segment.length = 0.01) +  # Add text labels for top 20 rows
        geom_text_repel(data = raw_sign_abund_feature_log[raw_sign_abund_feature_log$Annotation %in% names_to_label, ],
                        aes(label = Annotation), color = "darkgrey", max.overlaps = Inf, vjust = 5, min.segment.length = 0.01) +  # Add text labels for specific points with line always
        labs(x = "Mean Rank", y = "raw log2(Abundance)", color = NULL, fill = NULL,  # Remove color and fill legends from default mapping
             title = "Feature Ranking of raw log2 Abundance") +  # Add plot title              
        scale_color_manual(values = c("Mean" = "darkblue"), guide = guide_legend(title = "Mean")) +  # Manual color legend for mean
        scale_fill_manual(values = "skyblue", guide = guide_legend(title = "Standard Deviation")) +  # Manual fill legend for standard deviation
        theme(legend.position = "top") +  # Position legend at the top of the plot
        annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, 
                 label = paste("Number of Features:", num_features))  # Add number of features as text to top right corner
      
      # plot normalized data
      plot_list[[2]] <- ggplot(norm_sign_abund_feature_log, aes(x = mean_norm_log_rank, y = mean_norm_log)) +
        geom_ribbon(aes(ymin = mean_norm_log - std_dev_norm_log, ymax = mean_norm_log + std_dev_norm_log, fill = "Standard Deviation"), alpha = 0.4) +  # Add ribbon for standard deviation
        geom_point(aes(color = "Mean"), size = 1) +  # Add points for mean with bigger size
        geom_text_repel(data = norm_mean_log_top_20, aes(label = Annotation), color = "black", max.overlaps = Inf, min.segment.length = 0.01) +  # Add text labels for top 20 rows
        geom_text_repel(data = norm_sign_abund_feature_log[norm_sign_abund_feature_log$Annotation %in% names_to_label, ],
                        aes(label = Annotation), color = "darkgrey", max.overlaps = Inf, vjust = 5, min.segment.length = 0.01) +  # Add text labels for specific points with line always
        labs(x = "Mean Rank", y = "norm log2(Abundance)", color = NULL, fill = NULL,  # Remove color and fill legends from default mapping
             title = "Feature Ranking of normalized log2 Abundance") +  # Add plot title              
        scale_color_manual(values = c("Mean" = "darkblue"), guide = guide_legend(title = "Mean")) +  # Manual color legend for mean
        scale_fill_manual(values = "skyblue", guide = guide_legend(title = "Standard Deviation")) +  # Manual fill legend for standard deviation
        theme(legend.position = "top") +  # Position legend at the top of the plot
        annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, 
                 label = paste("Number of Features:", num_features, "\n",
                               "Sum of mean_norm (Top 20):", sum_mean_log_norm_top20, "\n",
                               "Sum of mean_norm (Rest):", sum_mean_log_norm_rest))  # Add text with number of features and sum of mean_norm to top right corner
      
      # Arrange the plots in a grid
      do.call(grid.arrange, plot_list)
    
  
  ### CHROM
    # plot the absolute values of Chrom paper
      plot_list <- list()
      
      # Plot mean and standard deviation for each row
      # plot raw data
      plot_list[[1]] <- ggplot(raw_sign_abund_feature_chrom, aes(x = mean_rank, y = mean)) +
        geom_ribbon(aes(ymin = mean - std_dev, ymax = mean + std_dev, fill = "Standard Deviation"), alpha = 0.4) +  # Add ribbon for standard deviation
        geom_point(aes(color = "Mean"), size = 1) +  # Add points for mean with bigger size
        geom_text_repel(data = mean_top_20_chrom, aes(label = Annotation), color = "black", max.overlaps = Inf, min.segment.length = 0.01) +  # Add text labels for top 20 rows
        geom_text_repel(data = raw_sign_abund_feature_chrom[raw_sign_abund_feature_chrom$Annotation %in% names_to_label, ],
                        aes(label = Annotation), color = "darkgrey", max.overlaps = Inf, vjust = 5, min.segment.length = 0.01) +  # Add text labels for specific points with line always
        labs(x = "Mean Rank", y = "raw Abundance", color = NULL, fill = NULL,  # Remove color and fill legends from default mapping
             title = "Feature Ranking of raw Abundance of Chrom Paper") +  # Add plot title              
        scale_color_manual(values = c("Mean" = "darkblue"), guide = guide_legend(title = "Mean")) +  # Manual color legend for mean
        scale_fill_manual(values = "skyblue", guide = guide_legend(title = "Standard Deviation")) +  # Manual fill legend for standard deviation
        theme(legend.position = "top") +  # Position legend at the top of the plot
        annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, 
                 label = paste("Number of Features:", num_features))  # Add number of features as text to top right corner
      
      # plot normalized data
      plot_list[[2]] <- ggplot(norm_sign_abund_feature_chrom, aes(x = mean_norm_rank, y = mean_norm)) +
        geom_ribbon(aes(ymin = mean_norm - std_dev_norm, ymax = mean_norm + std_dev_norm, fill = "Standard Deviation"), alpha = 0.4) +  # Add ribbon for standard deviation
        geom_point(aes(color = "Mean"), size = 1) +  # Add points for mean with bigger size
        geom_text_repel(data = norm_mean_top_20_chrom, aes(label = Annotation), color = "black", max.overlaps = Inf, min.segment.length = 0.01) +  # Add text labels for top 20 rows
        geom_text_repel(data = norm_sign_abund_feature_chrom[norm_sign_abund_feature_chrom$Annotation %in% names_to_label, ],
                        aes(label = Annotation), color = "darkgrey", max.overlaps = Inf, vjust = 5, min.segment.length = 0.01) +  # Add text labels for specific points with line always
        labs(x = "Mean Rank", y = "norm Abundance", color = NULL, fill = NULL,  # Remove color and fill legends from default mapping
             title = "Feature Ranking of normalized Abundance of Chrom Paper") +  # Add plot title              
        scale_color_manual(values = c("Mean" = "darkblue"), guide = guide_legend(title = "Mean")) +  # Manual color legend for mean
        scale_fill_manual(values = "skyblue", guide = guide_legend(title = "Standard Deviation")) +  # Manual fill legend for standard deviation
        theme(legend.position = "top") +  # Position legend at the top of the plot
        annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, 
                 label = paste("Number of Features:", num_features, "\n",
                               "Sum of mean_norm (Top 20):", sum_mean_norm_top20_chrom , "\n",
                               "Sum of mean_norm (Rest):", sum_mean_norm_rest_chrom))  # Add text with number of features and sum of mean_norm to top right corner
      
      # Arrange the plots in a grid
      do.call(grid.arrange, plot_list)
      
      
      
    # plot the log2 values of Chrom Paper
      plot_list <- list()
      
      # Plot mean and standard deviation for each row
      # plot raw data
      plot_list[[1]] <- ggplot(raw_sign_abund_feature_chrom_log, aes(x = mean_log_rank, y = mean_log)) +
        geom_ribbon(aes(ymin = mean_log - std_dev_log, ymax = mean_log + std_dev_log, fill = "Standard Deviation"), alpha = 0.4) +  # Add ribbon for standard deviation
        geom_point(aes(color = "Mean"), size = 1) +  # Add points for mean with bigger size
        geom_text_repel(data = mean_log_top_20_chrom, aes(label = Annotation), color = "black", max.overlaps = Inf, min.segment.length = 0.01) +  # Add text labels for top 20 rows
        geom_text_repel(data = raw_sign_abund_feature_chrom_log[raw_sign_abund_feature_chrom_log$Annotation %in% names_to_label, ],
                        aes(label = Annotation), color = "darkgrey", max.overlaps = Inf, vjust = 5, min.segment.length = 0.01) +  # Add text labels for specific points with line always
        labs(x = "Mean Rank", y = "raw log2(Abundance)", color = NULL, fill = NULL,  # Remove color and fill legends from default mapping
             title = "Feature Ranking of raw log2 Abundance of Chrom Paper") +  # Add plot title              
        scale_color_manual(values = c("Mean" = "darkblue"), guide = guide_legend(title = "Mean")) +  # Manual color legend for mean
        scale_fill_manual(values = "skyblue", guide = guide_legend(title = "Standard Deviation")) +  # Manual fill legend for standard deviation
        theme(legend.position = "top") +  # Position legend at the top of the plot
        annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, 
                 label = paste("Number of Features:", num_features))  # Add number of features as text to top right corner
      
      # plot normalized data
      plot_list[[2]] <- ggplot(norm_sign_abund_feature_chrom_log, aes(x = mean_norm_log_rank, y = mean_norm_log)) +
        geom_ribbon(aes(ymin = mean_norm_log - std_dev_norm_log, ymax = mean_norm_log + std_dev_norm_log, fill = "Standard Deviation"), alpha = 0.4) +  # Add ribbon for standard deviation
        geom_point(aes(color = "Mean"), size = 1) +  # Add points for mean with bigger size
        geom_text_repel(data = norm_mean_log_top_20_chrom, aes(label = Annotation), color = "black", max.overlaps = Inf, min.segment.length = 0.01) +  # Add text labels for top 20 rows
        geom_text_repel(data = norm_sign_abund_feature_chrom_log[norm_sign_abund_feature_chrom_log$Annotation %in% names_to_label, ],
                        aes(label = Annotation), color = "darkgrey", max.overlaps = Inf, vjust = 5, min.segment.length = 0.01) +  # Add text labels for specific points with line always
        labs(x = "Mean Rank", y = "norm log2(Abundance)", color = NULL, fill = NULL,  # Remove color and fill legends from default mapping
             title = "Feature Ranking of normalized log2 Abundance of Chrom Paper") +  # Add plot title              
        scale_color_manual(values = c("Mean" = "darkblue"), guide = guide_legend(title = "Mean")) +  # Manual color legend for mean
        scale_fill_manual(values = "skyblue", guide = guide_legend(title = "Standard Deviation")) +  # Manual fill legend for standard deviation
        theme(legend.position = "top") +  # Position legend at the top of the plot
        annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, 
                 label = paste("Number of Features:", num_features, "\n",
                               "Sum of mean_norm (Top 20):", sum_mean_log_norm_top20_chrom, "\n",
                               "Sum of mean_norm (Rest):", sum_mean_log_norm_rest_chrom ))  # Add text with number of features and sum of mean_norm to top right corner
      
      # Arrange the plots in a grid
      do.call(grid.arrange, plot_list)
    
      
    
  ### KIM
    # plot the absolute values of Kim paper
      plot_list <- list()
      
      # Plot mean and standard deviation for each row
      # plot raw data
      plot_list[[1]] <- ggplot(raw_sign_abund_feature_kim, aes(x = mean_rank, y = mean)) +
        geom_ribbon(aes(ymin = mean - std_dev, ymax = mean + std_dev, fill = "Standard Deviation"), alpha = 0.4) +  # Add ribbon for standard deviation
        geom_point(aes(color = "Mean"), size = 1) +  # Add points for mean with bigger size
        geom_text_repel(data = mean_top_20_kim, aes(label = Annotation), color = "black", max.overlaps = Inf, min.segment.length = 0.01) +  # Add text labels for top 20 rows
        geom_text_repel(data = raw_sign_abund_feature_kim[raw_sign_abund_feature_kim$Annotation %in% names_to_label, ],
                        aes(label = Annotation), color = "darkgrey", max.overlaps = Inf, vjust = 5, min.segment.length = 0.01) +  # Add text labels for specific points with line always
        labs(x = "Mean Rank", y = "raw Abundance", color = NULL, fill = NULL,  # Remove color and fill legends from default mapping
             title = "Feature Ranking of raw Abundance of Kim Paper") +  # Add plot title              
        scale_color_manual(values = c("Mean" = "darkblue"), guide = guide_legend(title = "Mean")) +  # Manual color legend for mean
        scale_fill_manual(values = "skyblue", guide = guide_legend(title = "Standard Deviation")) +  # Manual fill legend for standard deviation
        theme(legend.position = "top") +  # Position legend at the top of the plot
        annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, 
                 label = paste("Number of Features:", num_features))  # Add number of features as text to top right corner
      
      # plot normalized data
      plot_list[[2]] <- ggplot(norm_sign_abund_feature_kim, aes(x = mean_norm_rank, y = mean_norm)) +
        geom_ribbon(aes(ymin = mean_norm - std_dev_norm, ymax = mean_norm + std_dev_norm, fill = "Standard Deviation"), alpha = 0.4) +  # Add ribbon for standard deviation
        geom_point(aes(color = "Mean"), size = 1) +  # Add points for mean with bigger size
        geom_text_repel(data = norm_mean_top_20_kim, aes(label = Annotation), color = "black", max.overlaps = Inf, min.segment.length = 0.01) +  # Add text labels for top 20 rows
        geom_text_repel(data = norm_sign_abund_feature_kim[norm_sign_abund_feature_kim$Annotation %in% names_to_label, ],
                        aes(label = Annotation), color = "darkgrey", max.overlaps = Inf, vjust = 5, min.segment.length = 0.01) +  # Add text labels for specific points with line always
        labs(x = "Mean Rank", y = "norm Abundance", color = NULL, fill = NULL,  # Remove color and fill legends from default mapping
             title = "Feature Ranking of normalized Abundance of Kim Paper") +  # Add plot title              
        scale_color_manual(values = c("Mean" = "darkblue"), guide = guide_legend(title = "Mean")) +  # Manual color legend for mean
        scale_fill_manual(values = "skyblue", guide = guide_legend(title = "Standard Deviation")) +  # Manual fill legend for standard deviation
        theme(legend.position = "top") +  # Position legend at the top of the plot
        annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, 
                 label = paste("Number of Features:", num_features, "\n",
                               "Sum of mean_norm (Top 20):", sum_mean_norm_top20_kim, "\n",
                               "Sum of mean_norm (Rest):", sum_mean_norm_rest_kim))  # Add text with number of features and sum of mean_norm to top right corner
      
      # Arrange the plots in a grid
      do.call(grid.arrange, plot_list)
      
      
    # plot the log2 values of Kim paper
      plot_list <- list()
      
      # Plot mean and standard deviation for each row
      # plot raw data
      plot_list[[1]] <- ggplot(raw_sign_abund_feature_kim_log, aes(x = mean_log_rank, y = mean_log)) +
        geom_ribbon(aes(ymin = mean_log - std_dev_log, ymax = mean_log + std_dev_log, fill = "Standard Deviation"), alpha = 0.4) +  # Add ribbon for standard deviation
        geom_point(aes(color = "Mean"), size = 1) +  # Add points for mean with bigger size
        geom_text_repel(data = mean_log_top_20_kim, aes(label = Annotation), color = "black", max.overlaps = Inf, min.segment.length = 0.01) +  # Add text labels for top 20 rows
        geom_text_repel(data = raw_sign_abund_feature_kim_log[raw_sign_abund_feature_kim_log$Annotation %in% names_to_label, ],
                        aes(label = Annotation), color = "darkgrey", max.overlaps = Inf, vjust = 5, min.segment.length = 0.01) +  # Add text labels for specific points with line always
        labs(x = "Mean Rank", y = "raw log2(Abundance)", color = NULL, fill = NULL,  # Remove color and fill legends from default mapping
             title = "Feature Ranking of raw log2 Abundance of Kim Paper") +  # Add plot title              
        scale_color_manual(values = c("Mean" = "darkblue"), guide = guide_legend(title = "Mean")) +  # Manual color legend for mean
        scale_fill_manual(values = "skyblue", guide = guide_legend(title = "Standard Deviation")) +  # Manual fill legend for standard deviation
        theme(legend.position = "top") +  # Position legend at the top of the plot
        annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, 
                 label = paste("Number of Features:", num_features))  # Add number of features as text to top right corner
      
      # plot normalized data
      plot_list[[2]] <- ggplot(norm_sign_abund_feature_kim_log, aes(x = mean_norm_log_rank, y = mean_norm_log)) +
        geom_ribbon(aes(ymin = mean_norm_log - std_dev_norm_log, ymax = mean_norm_log + std_dev_norm_log, fill = "Standard Deviation"), alpha = 0.4) +  # Add ribbon for standard deviation
        geom_point(aes(color = "Mean"), size = 1) +  # Add points for mean with bigger size
        geom_text_repel(data = norm_mean_log_top_20_kim, aes(label = Annotation), color = "black", max.overlaps = Inf, min.segment.length = 0.01) +  # Add text labels for top 20 rows
        geom_text_repel(data = norm_sign_abund_feature_kim_log[norm_sign_abund_feature_kim_log$Annotation %in% names_to_label, ],
                        aes(label = Annotation), color = "darkgrey", max.overlaps = Inf, vjust = 5, min.segment.length = 0.01) +  # Add text labels for specific points with line always
        labs(x = "Mean Rank", y = "norm log2(Abundance)", color = NULL, fill = NULL,  # Remove color and fill legends from default mapping
             title = "Feature Ranking of normalized log2 Abundance of Kim Paper") +  # Add plot title              
        scale_color_manual(values = c("Mean" = "darkblue"), guide = guide_legend(title = "Mean")) +  # Manual color legend for mean
        scale_fill_manual(values = "skyblue", guide = guide_legend(title = "Standard Deviation")) +  # Manual fill legend for standard deviation
        theme(legend.position = "top") +  # Position legend at the top of the plot
        annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, 
                 label = paste("Number of Features:", num_features, "\n",
                               "Sum of mean_norm (Top 20):", sum_mean_log_norm_top20_kim, "\n",
                               "Sum of mean_norm (Rest):", sum_mean_log_norm_rest_kim))  # Add text with number of features and sum of mean_norm to top right corner
      
      # Arrange the plots in a grid
      do.call(grid.arrange, plot_list)
  
      
  
  # COMPARE Chrom and Kim Abundance in same plot
      plot_list <- list()
      
      # Plot mean and standard deviation for each row
      # plot raw data
      plot_list[[1]] <- ggplot(paper_comparison_data) +
        geom_ribbon(aes(x = mean_log_rank_chrom, ymin = mean_log_chrom - std_dev_log_chrom, ymax = mean_log_chrom + std_dev_log_chrom, fill = "Standard Deviation Chrom"), alpha = 0.4) +  # Add ribbon for standard deviation
        geom_point(aes(x = mean_log_rank_chrom, y = mean_log_chrom, color = "Mean Chrom"), size = 1) +  # Add points for mean with bigger size
        geom_ribbon(aes(x = mean_log_rank_kim, ymin = mean_log_kim - std_dev_log_kim, ymax = mean_log_kim + std_dev_log_kim, fill = "Standard Deviation Kim"), alpha = 0.4) +  # Add ribbon for standard deviation
        geom_point(aes(x = mean_log_rank_kim, y = mean_log_kim, color = "Mean Kim"), size = 1) +  # Add points for mean with bigger size
        labs(x = "Mean Rank", y = "raw log2(Abundance)", color = NULL, fill = NULL,  # Remove color and fill legends from default mapping
             title = "Feature Ranking of raw log2 Abundance (Chrom vs Kim)") +  # Add plot title              
        scale_color_manual(values = c("Mean Chrom" = "darkblue", "Mean Kim" = "darkred"), guide = guide_legend(title = "Mean")) +  # Manual color legend for mean
        scale_fill_manual(values = c("Standard Deviation Chrom" = "lightblue", "Standard Deviation Kim" = "pink2"), guide = guide_legend(title = "Standard Deviation")) +  # Manual fill legend for standard deviation
        theme(legend.position = "top") +  # Position legend at the top of the plot
        annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, 
                 label = paste("Number of Features:", num_features))  # Add number of features as text to top right corner
      
      # plot normalized data
      plot_list[[2]] <- ggplot(paper_comparison_data) +
        geom_ribbon(aes(x = mean_norm_log_rank_chrom, ymin = mean_norm_log_chrom - std_dev_norm_log_chrom, ymax = mean_norm_log_chrom + std_dev_norm_log_chrom, fill = "Standard Deviation Chrom"), alpha = 0.4) +  # Add ribbon for standard deviation
        geom_point(aes(x = mean_norm_log_rank_chrom, y = mean_norm_log_chrom, color = "Mean Chrom"), size = 1) +  # Add points for mean with bigger size
        geom_ribbon(aes(x = mean_norm_log_rank_kim, ymin = mean_norm_log_kim - std_dev_norm_log_kim, ymax = mean_norm_log_kim + std_dev_norm_log_kim, fill = "Standard Deviation Kim"), alpha = 0.4) +  # Add ribbon for standard deviation
        geom_point(aes(x = mean_norm_log_rank_kim, y = mean_norm_log_kim, color = "Mean Kim"), size = 1) +  # Add points for mean with bigger size
        labs(x = "Mean Rank", y = "raw log2(Abundance)", color = NULL, fill = NULL,  # Remove color and fill legends from default mapping
             title = "Feature Ranking of raw log2 Abundance (Chrom vs Kim)") +  # Add plot title              
        scale_color_manual(values = c("Mean Chrom" = "darkblue", "Mean Kim" = "darkred"), guide = guide_legend(title = "Mean")) +  # Manual color legend for mean
        scale_fill_manual(values = c("Standard Deviation Chrom" = "lightblue", "Standard Deviation Kim" = "pink2"), guide = guide_legend(title = "Standard Deviation")) +  # Manual fill legend for standard deviation
        theme(legend.position = "top") +  # Position legend at the top of the plot
        annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, 
                 label = paste("Number of Features:", num_features))  # Add number of features as text to top right corner
      
      # Arrange the plots in a grid
      do.call(grid.arrange, plot_list)
  
      
  # COMPARE Chrom and Kim with Scatter plot
    names_for_mean_comparison <- c("Paraxanthine", "Theobromine", "Theophylline", "Caffeine", "DL-Tyrosine", "Uric Acid")
    condition_mean <- with(paper_comparison_data, mean_log_chrom > 2 * mean_log_kim | 
                             mean_log_chrom * 2 < mean_log_kim)
    annot <- paper_comparison_data$Annotation[condition_mean]
    filtered_annot <- annot[!is.na(annot) & annot != ""]
    names_for_mean_comparison <- c(names_for_mean_comparison, filtered_annot)
    
    names_for_std_dev_comparison <- c("Paraxanthine", "Theobromine", "Theophylline", "Caffeine", "DL-Tyrosine", "Uric Acid")
    condition_std_dev <- with(paper_comparison_data, std_dev_log_chrom > 4 * std_dev_log_kim | 
                                std_dev_log_chrom * 4 < std_dev_log_kim)
    annot <- paper_comparison_data$Annotation[condition_std_dev]
    filtered_annot <- annot[!is.na(annot) & annot != ""]
    names_for_std_dev_comparison <- c(names_for_std_dev_comparison, filtered_annot)
      
    names_for_rank_comparison <- c("Paraxanthine", "Theobromine", "Theophylline", "Caffeine", "DL-Tyrosine", "Uric Acid")
    condition_rank <- with(paper_comparison_data, mean_log_rank_chrom > 4 * mean_log_rank_kim | 
                             mean_log_rank_chrom * 4 < mean_log_rank_kim)
    annot <- paper_comparison_data$Annotation[condition_rank]
    filtered_annot <- annot[!is.na(annot) & annot != ""]
    names_for_rank_comparison <- c(names_for_rank_comparison, filtered_annot)
    
      
    names_for_rank_comparison <- c("Paraxanthine", "Theobromine", "Theophylline", "Caffeine", "DL-Tyrosine", "Uric Acid")
    condition_rank <- with(paper_comparison_data, mean_log_rank_chrom > 3 * mean_log_rank_kim | 
                                              mean_log_rank_chrom * 3 < mean_log_rank_kim)
    annot <- paper_comparison_data$Annotation[condition]
    filtered_annot <- annot[!is.na(annot) & annot != ""]
    names_for_rank_comparison <- c(names_for_comparison, filtered_annot)
    
    
    # plot the log2 mean values of chrom and kim as comparison
      plot_list <- list()
      
      # plot raw log2 mean data
      plot_list[[1]] <- ggplot(paper_comparison_data, aes(x = mean_log_chrom, y = mean_log_kim)) +
        geom_point() +  # Add points for the scatterplot
        #geom_text_repel(data = paper_comparison_data[paper_comparison_data$Annotation %in% names_for_mean_comparison, ],
        #                aes(label = Annotation), color = "black", min.segment.length = 0.01) +
        labs(x = "log2(Abundance) Chrom", y = "log2(Abundance) Kim",
             title = "Comparison of Feature Abundance (Chrom vs Kim)") +
        coord_fixed() +  # Make the plot square
        xlim(c(min(paper_comparison_data$mean_log_chrom, paper_comparison_data$mean_log_kim),
               max(paper_comparison_data$mean_log_chrom, paper_comparison_data$mean_log_kim))) +  # Set same x-axis limits
        ylim(c(min(paper_comparison_data$mean_log_chrom, paper_comparison_data$mean_log_kim),
               max(paper_comparison_data$mean_log_chrom, paper_comparison_data$mean_log_kim)))  # Set same y-axis limits
      
      # plot log2 normalized mean data
      plot_list[[2]] <- ggplot(paper_comparison_data, aes(x = mean_norm_log_chrom, y = mean_norm_log_kim)) +
        geom_point() +  # Add points for the scatterplot
        #geom_text_repel(data = paper_comparison_data[paper_comparison_data$Annotation %in% names_for_mean_comparison, ],
        #                aes(label = Annotation), color = "black", min.segment.length = 0.01) +
        labs(x = "log2(norm. Abundance) Chrom", y = "log2(norm. Abundance) Kim",
             title = "Comparison of normalized Feature Abundance (Chrom vs Kim)") +
        coord_fixed() +  # Make the plot square
        xlim(c(min(paper_comparison_data$mean_norm_log_chrom, paper_comparison_data$mean_norm_log_kim),
               max(paper_comparison_data$mean_norm_log_chrom, paper_comparison_data$mean_norm_log_kim))) +  # Set same x-axis limits
        ylim(c(min(paper_comparison_data$mean_norm_log_chrom, paper_comparison_data$mean_norm_log_kim),
               max(paper_comparison_data$mean_norm_log_chrom, paper_comparison_data$mean_norm_log_kim)))  # Set same y-axis limits
      
      # Arrange the plots in a grid
      do.call(grid.arrange, plot_list)    
      
      
    
    # plot the log2 std dev values of chrom and kim as comparison
      plot_list <- list()
      
      # plot raw log2 std_dev data
      plot_list[[1]] <- ggplot(paper_comparison_data, aes(x = std_dev_log_chrom, y = std_dev_log_kim)) +
        geom_point() +  # Add points for the scatterplot
        #geom_text_repel(data = paper_comparison_data[paper_comparison_data$Annotation %in% names_for_std_dev_comparison, ],
        #                aes(label = Annotation), color = "black", min.segment.length = 0.01) +
        labs(x = "log2(Abundance) Chrom", y = "log2(Abundance) Kim",  
             title = "Comparison of Std. dev. of Feature Abundance (Chrom vs Kim)") +
        coord_fixed() +  # Make the plot square
        xlim(c(min(paper_comparison_data$std_dev_log_chrom, paper_comparison_data$std_dev_log_kim),
               max(paper_comparison_data$std_dev_log_chrom, paper_comparison_data$std_dev_log_kim))) +  # Set same x-axis limits
        ylim(c(min(paper_comparison_data$std_dev_log_chrom, paper_comparison_data$std_dev_log_kim),
               max(paper_comparison_data$std_dev_log_chrom, paper_comparison_data$std_dev_log_kim)))  # Set same y-axis limits
      
      
      # plot log2 normalized std_dev data
      plot_list[[2]] <- ggplot(paper_comparison_data, aes(x = std_dev_norm_log_chrom, y = std_dev_norm_log_kim)) +
        geom_point() +  # Add points for the scatterplot
        #geom_text_repel(data = paper_comparison_data[paper_comparison_data$Annotation %in% names_for_std_dev_comparison, ],
        #                aes(label = Annotation), color = "black", min.segment.length = 0.01) +
        labs(x = "log2(norm. Abundance) Chrom", y = "log2(norm. Abundance) Kim",  
             title = "Comparison of Std. dev. of normalized Feature Abundance (Chrom vs Kim)") +
        coord_fixed() +  # Make the plot square
        xlim(c(min(paper_comparison_data$std_dev_norm_log_chrom, paper_comparison_data$std_dev_norm_log_kim),
               max(paper_comparison_data$std_dev_norm_log_chrom, paper_comparison_data$std_dev_norm_log_kim))) +  # Set same x-axis limits
        ylim(c(min(paper_comparison_data$std_dev_norm_log_chrom, paper_comparison_data$std_dev_norm_log_kim),
               max(paper_comparison_data$std_dev_norm_log_chrom, paper_comparison_data$std_dev_norm_log_kim)))  # Set same y-axis limits
      
      
      # Arrange the plots in a grid
      do.call(grid.arrange, plot_list)  
      
      
    
    # plot the log2 rank values of chrom and kim as comparison
      plot_list <- list()
      
      # plot raw log2 rank data
      plot_list[[1]] <- ggplot(paper_comparison_data, aes(x = mean_log_rank_chrom, y = mean_log_rank_kim)) +
        geom_point() +  # Add points for the scatterplot
        #geom_text_repel(data = paper_comparison_data[paper_comparison_data$Annotation %in% names_for_rank_comparison, ],
        #                aes(label = Annotation), color = "black", min.segment.length = 0.01) +
        labs(x = "rank of log2(mean Abundance) Chrom", y = "rank of log2(mean Abundance) Kim",  
             title = "Comparison of Rank of Feature Abundance (Chrom vs Kim)") +
        coord_fixed() +  # Make the plot square
        xlim(c(min(paper_comparison_data$mean_log_rank_chrom, paper_comparison_data$mean_log_rank_kim),
               max(paper_comparison_data$mean_log_rank_chrom, paper_comparison_data$mean_log_rank_kim))) +  # Set same x-axis limits
        ylim(c(min(paper_comparison_data$mean_log_rank_chrom, paper_comparison_data$mean_log_rank_kim),
               max(paper_comparison_data$mean_log_rank_chrom, paper_comparison_data$mean_log_rank_kim)))  # Set same y-axis limits
      
      
      # plot log2 normalized rank data
      plot_list[[2]] <- ggplot(paper_comparison_data, aes(x = mean_norm_log_rank_chrom, y = mean_norm_log_rank_kim)) +
        geom_point() +  # Add points for the scatterplot
        #geom_text_repel(data = paper_comparison_data[paper_comparison_data$Annotation %in% names_for_rank_comparison, ],
        #                aes(label = Annotation), color = "black", min.segment.length = 0.01) +
        labs(x = "rank of log2(norm. mean Abundance) Chrom", y = "rank of log2(norm. mean Abundance) Kim",  
             title = "Comparison of Rank of normalized Feature Abundance (Chrom vs Kim)") +
        coord_fixed() +  # Make the plot square
        xlim(c(min(paper_comparison_data$mean_norm_log_rank_chrom, paper_comparison_data$mean_norm_log_rank_kim),
               max(paper_comparison_data$mean_norm_log_rank_chrom, paper_comparison_data$mean_norm_log_rank_kim))) +  # Set same x-axis limits
        ylim(c(min(paper_comparison_data$mean_norm_log_rank_chrom, paper_comparison_data$mean_norm_log_rank_kim),
               max(paper_comparison_data$mean_norm_log_rank_chrom, paper_comparison_data$mean_norm_log_rank_kim)))  # Set same y-axis limits
      
      
      # Arrange the plots in a grid
      do.call(grid.arrange, plot_list)  
      
        
  dev.off() # Close the PDF device
  
  










# try some normalization plots by normalising witht he mean abundance, but not really needed

experiment_data_plots <- raw_sign_abund_feature %>%
  pivot_longer(cols = -c(id, rt, mz, rtime_group, feature_group, Annotation, Confidence, SMILES), names_to = "Sample", values_to = "Area") %>%
  mutate(Sample = factor(Sample, levels = unique(Sample))) 

# Add additional information (like sample times, different condition like different paper etc)
add_info_file <- paste0(wd, exp,"/FiS_45min-times.csv")
add_info <- read.csv2(add_info_file)
add_info$mean_abundance <- mean_abundance_each_col

experiment_data_plots <- merge(experiment_data_plots, add_info, by.x = "Sample", by.y = "Sample", all.x = TRUE)


# change the data values to the according format
experiment_data_plots$date <- as.Date(experiment_data_plots$date, format = "%d.%m.%Y")

# Convert time column to minutes between sampling
experiment_data_plots$time <- sapply(experiment_data_plots$time, convert_to_minutes)



unique_molecules <- unique(raw_sign_abund_feature$Annotation)
normalization_molecules <- c('DL-Tyrosine', 'Uric Acid', 'mean abundance')

# Set up PDF device
pdf(file = paste0(resultsdir, "45min-Papercomparison-feature_data3.pdf"), width = 15, height = 8) 

# Loop over each unique molecule
for (molecule in unique_molecules) {
  plot_list <- list()
  
  # subsets for the current molecule
  filtered_data <- subset(experiment_data_plots, Annotation == molecule)
  
  # raw data plot
  raw_plot <- ggplot(data = filtered_data, 
                     aes(x = time, y = Area, color = interaction(paper, date, hand))) +
    geom_line() +
    geom_point(size = 3) +
    labs(x = "Time (minutes)", y = "Area", color = "Paper", 
         title = paste(molecule, "over time (Raw Data)")) +
    theme_minimal()
  
  plot_list[[length(plot_list) + 1]] <- raw_plot
  
  # loop over each normalization variant
  for (norm_mol in normalization_molecules) {
    if (norm_mol %in% unique_molecules) {
      # subsets for normalization variant
      norm_data <- subset(experiment_data_plots, Annotation == norm_mol)
      
      # normalization
      filtered_data_norm <- merge(filtered_data, norm_data, by = "Sample", suffixes = c("", "_norm_mol"))
      filtered_data_norm$Normalized_Area <- filtered_data_norm$Area / filtered_data_norm$Area_norm_mol
      
      # normalized data plot
      normalized_plot <- ggplot(data = filtered_data_norm, 
                                aes(x = time, y = Normalized_Area, color = interaction(paper, date, hand))) +
        geom_line() +
        geom_point(size = 3) +
        labs(x = "Time (minutes)", y = "Normalized Area", color = "Paper", 
             title = paste(norm_mol, "normalized", molecule, "over time")) +
        theme_minimal()
      
      plot_list[[length(plot_list) + 1]] <- normalized_plot
    } else {
      # normalization
      filtered_data_norm <- filtered_data
      filtered_data_norm$Normalized_Area <- filtered_data_norm$Area / filtered_data_norm$mean_abundance
      
      # normalized data plot
      normalized_plot <- ggplot(data = filtered_data_norm, 
                                aes(x = time, y = Normalized_Area, color = interaction(paper, date, hand))) +
        geom_line() +
        geom_point(size = 3) +
        labs(x = "Time (minutes)", y = "Normalized Area", color = "Paper", 
             title = paste(norm_mol, "normalized", molecule, "over time")) +
        theme_minimal()
      
      plot_list[[length(plot_list) + 1]] <- normalized_plot
    }
  }
  
  # combine all plots for the current molecule into one plot
  do.call(grid.arrange, plot_list)
}
# Close PDF device
dev.off()




