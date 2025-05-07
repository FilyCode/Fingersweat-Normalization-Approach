library(ggplot2)
library(dplyr)
library(tidyr)
library(fmsb)

wd = "X:/Studenten_Schueler/Philipp_Trollmann/Experiments/"
setwd(wd)

exp = "Dataset-comparison"

datafile = paste0(wd, exp, "/dataset-comparison-biological-features.csv")
resultsdir = paste0(wd, exp)
data <- read.csv(datafile)


# Function to create boxplots and save to pdf
create_plots_for_molecule_comparison <- function(df, sample_patterns, donor_pattern = NULL) {
  # Get unique molecules
  unique_molecules <- unique(df$Molecule.Name)
  
  # Check if "Caffein-D9" is in the unique molecules
  if (!"Caffein-D9" %in% unique_molecules) {
    stop("Caffein-D9 is not found in the Molecule.Name column.")
  }
  
  # Create a new column to classify the samples based on patterns
  df$Sample.Group <- sapply(df$Sample, function(sample) {
    match <- sapply(sample_patterns, function(pattern) grepl(pattern, sample))
    if (any(match)) {
      return(sample_patterns[which(match)])
    } else {
      return("Other")
    }
  })
  
  # Create a column to store the normalized area
  df$Normalized.Area <- NA
  
  # Loop through each sample to normalize
  for (sample in unique(df$Sample)) {
    # Get the Caffein-D9 value for the current sample
    caffein_d9_value <- df$Area[df$Sample == sample & df$Molecule.Name == "Caffein-D9"]
    
    if (length(caffein_d9_value) == 0) {
      stop(paste("No Caffein-D9 value found for sample", sample))
    }
    
    # Normalize the Area values by the Caffein-D9 value
    df$Normalized.Area[df$Sample == sample] <- df$Area[df$Sample == sample] / caffein_d9_value
  }
  
  df$Normalized.Area[!is.finite(df$Normalized.Area)] <- NA
  
  # Create a PDF to save the plots
  pdf(file = paste0(resultsdir, "/biological-features-comparision_boxplots_Caffein-D9-normalized.pdf"), width = 15, height = 8)
  
  # Loop through each molecule (excluding "Caffein-D9")
  for (molecule in unique_molecules[unique_molecules != "Caffein-D9"]) {
    # Subset the dataframe for the current molecule
    molecule_df <- subset(df, Molecule.Name == molecule)
    
    # Count the number of samples in each group
    sample_counts <- table(molecule_df$Sample.Group)
    sample_labels <- paste0(names(sample_counts), "\n(n = ", sample_counts, ")")
    names(sample_labels) <- names(sample_counts)
    
    # Create the boxplot
    p <- ggplot(molecule_df, aes(x = Sample.Group, y = Normalized.Area)) +
      geom_boxplot() +
      scale_x_discrete(labels = sample_labels) +
      labs(title = paste("Comparison of", molecule, "between different studies"),
           x = "Sample Group",
           y = "Normalized Area") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45))
    
    # Print the plot to the PDF
    print(p)
  }
  
  p <- ggplot(df[df$Molecule.Name != "Caffein-D9", ], aes(x = Molecule.Name, y = Normalized.Area)) +
    geom_boxplot() +
    labs(title = paste("Comparison of biological molecules (data includes various studies)"),
         x = "Molecules",
         y = "Normalized Area") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90))
  
  # Print the plot to the PDF
  print(p)
  
  
  # Create a new plot with each molecule scaled from 0 to 1
  df$Scaled.Normalized.Area <- NA
  for (molecule in unique_molecules) {
    molecule_values <- df$Normalized.Area[df$Molecule.Name == molecule]
    max_value <- max(molecule_values, na.rm = TRUE)
    df$Scaled.Normalized.Area[df$Molecule.Name == molecule] <- df$Normalized.Area[df$Molecule.Name == molecule] / max_value
  }
  
  p <- ggplot(df[df$Molecule.Name != "Caffein-D9", ], aes(x = Molecule.Name, y = Scaled.Normalized.Area)) +
    geom_boxplot() +
    labs(title = paste("Comparison of biological molecules (scaled from 0 to 1)"),
         x = "Molecules",
         y = "Scaled normalized Area") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90))
  
  # Print the plot to the PDF
  print(p)
  
  
  # Plot the non Caffein-D9 normalized data for comparison
  p <- ggplot(df, aes(x = Molecule.Name, y = Area)) +
    geom_boxplot() +
    labs(title = paste("Comparison of non-normalized biological molecules (data includes various studies)"),
         x = "Molecules",
         y = "Area") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90))
  
  # Print the plot to the PDF
  print(p)
  
  
  # Calculate Relative Standard Deviation (RSD) for each molecule in each sample group
  rsd_df <- df %>%
    group_by(Sample.Group, Molecule.Name) %>%
    summarize(RSD = sd(Area) / mean(Area) * 100, .groups = 'drop') %>%
    pivot_wider(names_from = Molecule.Name, values_from = RSD)
  
  # Prepare the data for radar plot
  #rsd_df <- rsd_df[complete.cases(rsd_df), ]  # Remove rows with NA
  rsd_df <- as.data.frame(rsd_df)
  rownames(rsd_df) <- rsd_df$Sample.Group
  rsd_df <- rsd_df[, -1]  # Remove Sample.Group column
  
  # Add max and min values for the radar plot
  rsd_df <- rbind("max" = rep(max(rsd_df, na.rm = TRUE) + 10, ncol(rsd_df)),  # Adding +10 to make scaling more appropiate
                  "min" = rep(0, ncol(rsd_df)),
                  rsd_df)
  
  # Calculate the maximum value and step size for the axis labels
  max_val <- max(rsd_df, na.rm = TRUE)
  num_intervals <- 4
  step_size <- floor(max_val / num_intervals)
  
  # Create the radar plot
  par(mar = c(2, 2, 2, 7), xpd = TRUE)
  radarchart(rsd_df, axistype = 1,
             pcol = rainbow(nrow(rsd_df)-2),
             plwd = 1, plty = 1,
             cglcol = "grey", cglty = 1, cglwd = 0.8, vlcex = 0.8, 
             axislabcol = "grey", caxislabels = seq(0, max_val, by = step_size),
             title = "Relative Standard Deviation (RSD) of Molecules Across Experiments")
  
  legend(
    x = "topright", legend = rownames(rsd_df)[3:nrow(rsd_df)], horiz = FALSE,
    bty = "n", pch = 20, col = rainbow(nrow(rsd_df)-2),
    text.col = "black", cex = 1.2, pt.cex = 1.5
  )
  
 
  # Loop through each experiment in donor_pattern
  for (experiment in names(donor_pattern)) {
    # Initialize storage for RSD values for donors within the current experiment
    donor_rsd_list <- list()
    
    for (donor in donor_pattern[[experiment]]) {
      # Filter data for each donor within the current experiment
      donor_df <- df %>%
        filter(grepl(donor, Sample) & Sample.Group == experiment)
      
      # Calculate RSD for the donor within the current experiment
      donor_rsd <- donor_df %>%
        group_by(Molecule.Name) %>%
        summarize(RSD = sd(Area) / mean(Area) * 100, .groups = 'drop')
      
      # Store RSD values with donor name as row identifier
      donor_rsd_list[[donor]] <- donor_rsd$RSD
      names(donor_rsd_list[[donor]]) <- donor_rsd$Molecule.Name
    }
    
    # Combine RSD values of all donors into a dataframe
    donor_rsd_df <- as.data.frame( t( do.call(cbind, donor_rsd_list))) # make to dataframe and transpose

    # Add max and min values for the radar plot
    max_val <- max(donor_rsd_df, na.rm = TRUE) + 10
    min_val <- 0
    donor_rsd_df <- rbind("max" = rep(max_val, ncol(donor_rsd_df)), 
                          "min" = rep(min_val, ncol(donor_rsd_df)), 
                          donor_rsd_df)
    
    # Calculate the step size for the axis labels
    step_size <- floor(max_val / num_intervals)
    
   # Create the radar plot for the current experiment
    par(mar = c(2, 2, 2, 7), xpd = TRUE)
    radarchart(donor_rsd_df, axistype = 1,
               pcol = rainbow(nrow(donor_rsd_df)-2),
               plwd = 1, plty = 1,
               cglcol = "grey", cglty = 1, cglwd = 0.8, vlcex = 0.8, 
               axislabcol = "grey", caxislabels = seq(0, max_val, by = step_size),
               title = paste("Relative Standard Deviation (RSD) of Molecules for", experiment))
    
    legend(
      x = "topright", legend = rownames(donor_rsd_df)[3:nrow(donor_rsd_df)], horiz = FALSE,
      bty = "n", pch = 20, col = rainbow(nrow(donor_rsd_df)-2),
      text.col = "black", cex = 1.2, pt.cex = 1.5
    )
    
    
  }
  
  # Close the PDF
  dev.off()
}


sample_patterns <- c("MTX", "PhilKoff", "Weinstudie", "KaMa", "Granat", "LongCovid", "Kepler", "AttZol", "ProcBlank")
donor_pattern <- list("PhilKoff" = c("Chrom", "Kim"), 
                      "Weinstudie" = c("CG", "SMM", "MW", "YF"),
                      "KaMa" = c("P1", "P2", "P3", "P4", "P5", "P6"),
                      "Granat" = c("AH", "DW", "GH", "GrH", "KK", "MW"))

# Call the function
create_plots_for_molecule_comparison(data, sample_patterns, donor_pattern)

