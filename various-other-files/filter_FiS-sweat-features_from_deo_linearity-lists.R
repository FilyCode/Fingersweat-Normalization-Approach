source("//prot-data/MS Storage 3/tmp_PT/R Project/MS-data-analysis_functions.R")
source("//prot-data/MS Storage 3/tmp_PT/R Project/feature-visualization-script.R")

library(fuzzyjoin)

# set all the experiment and file names
wd = "//prot-data/MS Storage 3/tmp_PT/Experiments/"
setwd(wd)
exp = "FiS-good-sweat-features/"

resultsdir = paste0(wd, exp)

sweat_features <- read.csv(paste0(resultsdir, "GOOD_SWEAT_FEATURES.csv"))
contamination_features <- read.csv(paste0(resultsdir, "FiS_limma-filter_aludeo-effect.csv"))
good_features <- read.csv(paste0(resultsdir, "250411_good_features_procblank_dilution.csv"))
good_features_all <- read.csv(paste0(resultsdir, "compressed_dilution_data.csv"))


# Select required columns from sweat_features
sweat_selected <- sweat_features %>%
  select(rt, mz, charge, Annotation, Formula, Confidence, n, L1_L2, L2_L3, R1_R2, R2_R3)
sweat_selected$mean_sweat_change <- rowMeans(sweat_selected[, c("L1_L2", "R2_R3")], na.rm = TRUE)
sweat_selected$mean_sd_sweat_stability <- apply(sweat_selected, 1, function(row) {
  sd1 <- sd(c(row["L2_L3"], row["R1_R2"]), na.rm = TRUE)
  sd2 <- sd(c(row["L1_L2"], row["R2_R3"]), na.rm = TRUE)
  mean(c(sd1, sd2), na.rm = TRUE)
})

# Add linearity data to good feature subset 
good_features <- good_features %>% mutate(rt_round = round(rt, 3), mz_round = round(mz, 3))
good_features_all <- good_features_all %>% mutate(rt_round = round(rt, 3), mz_round = round(mz, 3))

good_features_merged <- good_features %>%
  left_join(good_features_all %>% select(rt_round, mz_round, charge, r_squared, avg_residuals, sd_residuals), 
            by = c("rt_round" = "rt_round", "mz_round" = "mz_round", "charge" = "charge"))

# Perform fuzzy join on mz and rt with conditions
merged_df <- fuzzy_inner_join(
  sweat_selected, good_features_merged,
    by = c("mz" = "mz", "rt" = "rt", "charge" = "charge"),
    match_fun = list(
      function(x, y) abs(x - y) <= 0.005,  # MZ within ±0.005
      function(x, y) abs(x - y) <= 0.1,  # RT will be handled later
      `==`  # Charge must be identical
    )
  ) %>%
  group_by(rt.x, mz.x, charge.x) %>%  # Group by original features
  slice_min(abs(rt.x - rt.y), n = 1) %>%  # Pick the closest RT
  ungroup() %>%
  select(-rt.y, -mz.y, -charge.y) %>%  # Remove duplicate RT and MZ columns
  dplyr::rename(rt = rt.x, mz = mz.x, charge = charge.x)  # Rename back

# Choose best annotation based on highest confidence
merged_df$Confidence <- pmax(merged_df$Confidence.x, merged_df$Confidence.y, na.rm = TRUE)
merged_df$Annotation <- ifelse(
  merged_df$Confidence.x >= merged_df$Confidence.y,
  merged_df$Annotation.x,
  merged_df$Annotation.y
)
merged_df$Formula <- ifelse(
  merged_df$Confidence.x >= merged_df$Confidence.y,
  merged_df$Formula.x,
  merged_df$Formula.y
)

# For sweat
merged_df_filtered <- merged_df %>% filter(r_squared > 0.99, Confidence > 0.2, mean_sweat_change < -4, mean_sd_sweat_stability < 1.5)
write.csv(merged_df_filtered, paste0(resultsdir, "filtered_good_sweat-features.csv"), row.names = FALSE)


filtered_good_sweat_norm_mols <- read.csv("//prot-data/MS Storage 3/tmp_PT/Experiments/normalization/filtered_good_sweat-features_for-norm.csv")

plot <- ggplot(good_features_merged, aes(x = rt, y = r_squared, label = Annotation)) +
  geom_point(alpha = 0.6) +
  labs(title = "R^2 vs. RT of all Molecules",
       x = "RT (min)", y = "R^2") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18), 
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14) 
  )
ggsave(paste0(resultsdir, "R2-vs-RT-scatter-plot-all.png"), plot = plot, width = 8, height = 6, dpi = 300)


plot <- ggplot(sweat_selected, aes(x = rt, y = mean_sweat_change, label = Annotation)) +
  geom_point(alpha = 0.6) +
  labs(title = "mean Sweat Change vs. RT of Sweat Molecules",
       x = "RT (min)", y = "mean Sweat Change") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18), 
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14) 
  )
ggsave(paste0(resultsdir, "meanSweatChange-vs-RT-scatter-plot-all.png"), plot = plot, width = 8, height = 6, dpi = 300)


plot <- ggplot(merged_df, aes(x = rt, y = r_squared, label = Annotation)) +
  geom_point(alpha = 0.6) +
  geom_point(data = filtered_good_sweat_norm_mols, aes(x = rt, y = r_squared), color = "red") +
  labs(title = "R^2 vs. RT of merged Molecules",
       x = "RT (min)", y = "R^2") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18), 
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14) 
  )
ggsave(paste0(resultsdir, "R2-vs-RT-scatter-plot-merged.png"), plot = plot, width = 8, height = 6, dpi = 300)


plot <- ggplot(merged_df, aes(x = rt, y = mean_sweat_change, label = Annotation)) +
  geom_point(alpha = 0.6) +
  geom_point(data = filtered_good_sweat_norm_mols, aes(x = rt, y = mean_sweat_change), color = "red") +
  labs(title = "mean Sweat Change  vs. RT of merged Molecules",
       x = "RT (min)", y = "mean Sweat Change ") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18), 
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14) 
  )
ggsave(paste0(resultsdir, "meanSweatChange-vs-RT-scatter-plot-merged.png"), plot = plot, width = 8, height = 6, dpi = 300)



plot <- ggplot(merged_df_filtered, aes(x = rt, y = r_squared, label = Annotation)) +
  geom_point(alpha = 0.6) +
  geom_point(data = filtered_good_sweat_norm_mols, aes(x = rt, y = r_squared), color = "red") +
  labs(title = "R^2 vs. RT of merged Molecules",
       x = "RT (min)", y = "R^2") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18), 
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14) 
  )
ggsave(paste0(resultsdir, "R2-vs-RT-scatter-plot-merged-filtered.png"), plot = plot, width = 8, height = 6, dpi = 300)


plot <- ggplot(merged_df_filtered, aes(x = rt, y = mean_sweat_change, label = Annotation)) +
  geom_point(alpha = 0.6) +
  geom_point(data = filtered_good_sweat_norm_mols, aes(x = rt, y = mean_sweat_change), color = "red") +
  labs(title = "mean Sweat Change  vs. RT of merged Molecules",
       x = "RT (min)", y = "mean Sweat Change ") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18), 
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14) 
  )
ggsave(paste0(resultsdir, "meanSweatChange-vs-RT-scatter-plot-merged-filtered.png"), plot = plot, width = 8, height = 6, dpi = 300)




##### FOR CONTAMINATION
exp = "FiS-good-sweat-features/deoFeatures/"

resultsdir = paste0(wd, exp)

contamination_features <- read.csv(paste0(resultsdir, "FiS_limma-filter_aludeo-effect.csv"))

contamination_selected <- contamination_features %>%
  select(rt, mz, charge, Annotation, Formula, Confidence, logFC, adj.P.Val)

# Add linearity data to good feature subset 
good_features <- good_features %>% mutate(rt_round = round(rt, 3), mz_round = round(mz, 3))
good_features_all <- good_features_all %>% mutate(rt_round = round(rt, 3), mz_round = round(mz, 3))

good_features_all_filtered <- good_features_all %>% filter(r_squared > 0.99, avg_residuals > -0.01, avg_residuals < 0.002, sd_residuals < 0.035)

good_features_merged <- good_features %>%
  left_join(good_features_all %>% select(rt_round, mz_round, charge, r_squared, avg_residuals, sd_residuals), 
            by = c("rt_round" = "rt_round", "mz_round" = "mz_round", "charge" = "charge"))

# Perform fuzzy join on mz and rt with conditions
merged_df <- fuzzy_inner_join(
  contamination_selected, good_features_all_filtered,
  by = c("mz" = "mz", "rt" = "rt", "charge" = "charge"),
  match_fun = list(
    function(x, y) abs(x - y) <= 0.005,  # MZ within ±0.005
    function(x, y) abs(x - y) <= 0.1,  # RT will be handled later
    `==`  # Charge must be identical
  )
) %>%
  group_by(rt.x, mz.x, charge.x) %>%  # Group by original features
  slice_min(abs(rt.x - rt.y), n = 1) %>%  # Pick the closest RT
  ungroup() %>%
  select(-rt.y, -mz.y, -charge.y) %>%  # Remove duplicate RT and MZ columns
  dplyr::rename(rt = rt.x, mz = mz.x, charge = charge.x)  # Rename back

# Choose best annotation based on highest confidence
merged_df$Confidence <- pmax(merged_df$Confidence.x, merged_df$Confidence.y, na.rm = TRUE)
merged_df$Annotation <- ifelse(
  merged_df$Confidence.x >= merged_df$Confidence.y,
  merged_df$Annotation.x,
  merged_df$Annotation.y
)


# For contamination
merged_df_filtered <- merged_df %>% filter(r_squared > 0.99, Confidence > 0.2, logFC > 8, adj.P.Val < 0.01)
merged_df_filtered$Abbrevation <- abbreviate(merged_df_filtered$Annotation, minlength = 15)

write.csv(merged_df_filtered, paste0(resultsdir, "deofiltered_linear-contamination-features.csv"), row.names = FALSE)


filtered_contamination_norm_mols <- read.csv("//prot-data/MS Storage 3/tmp_PT/Experiments/normalization/biological_normalization_molecule_list_for-combination-test-for-wine-study_filtered-linear-contamination-features.csv")

plot <- ggplot(good_features_all_filtered, aes(x = rt, y = r_squared, label = Annotation)) +
  geom_point(alpha = 0.6) +
  labs(title = "R^2 vs. RT of all Molecules",
       x = "RT (min)", y = "R^2") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18), 
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14) 
  )
ggsave(paste0(resultsdir, "R2-vs-RT-scatter-plot-all.png"), plot = plot, width = 8, height = 6, dpi = 300)


plot <- ggplot(contamination_selected, aes(x = rt, y = logFC, label = Annotation)) +
  geom_point(alpha = 0.6) +
  labs(title = "logFC vs. RT of Sweat Molecules",
       x = "RT (min)", y = "logFC") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18), 
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14) 
  )
ggsave(paste0(resultsdir, "logFC-vs-RT-scatter-plot-all.png"), plot = plot, width = 8, height = 6, dpi = 300)


plot <- ggplot(merged_df, aes(x = rt, y = r_squared)) +
  geom_point(alpha = 0.6) +
  geom_point(data = filtered_contamination_norm_mols, aes(x = rt, y = r_squared), color = "red") +
  labs(title = "R^2 vs. RT of merged Molecules",
       x = "RT (min)", y = "R^2") +
  scale_x_continuous(limits = c(0, 5)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18), 
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14) 
  )
ggsave(paste0(resultsdir, "R2-vs-RT-scatter-plot-merged.png"), plot = plot, width = 8, height = 6, dpi = 300)


plot <- ggplot(merged_df, aes(x = rt, y = logFC)) +
  geom_point(alpha = 0.6) +
  geom_point(data = filtered_contamination_norm_mols, aes(x = rt, y = logFC), color = "red") +
  labs(title = "logFC vs. RT of merged Molecules",
       x = "RT (min)", y = "logFC ") +
  scale_x_continuous(limits = c(0, 5)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18), 
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14) 
  )
ggsave(paste0(resultsdir, "logFC-vs-RT-scatter-plot-merged.png"), plot = plot, width = 8, height = 6, dpi = 300)

