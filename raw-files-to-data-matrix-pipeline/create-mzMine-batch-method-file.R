source(mzMine_batch-method-steps.R)


generate_xml_file <- function(output_file_path, beginning, ms_peaks_substraction, ms_peak_join_pick, merge, used_MS_methods, methods_parameters, sample_name_pattern) {
  # Get mzMine Batch method steps as strings to variable for file generation
  all_strings <- get_mzMine_batch_method_steps()
  beginning <- all_strings$beginning
  ms_peaks_substraction <- all_strings$ms_peaks_substraction
  ms_peak_join_pick <- all_strings$ms_peak_join_pick
  merge <- all_strings$merge
  
  
  # Define XML content
  xml_content <- c(beginning)

# Generate XML content for each MS method and sample name pattern
for (ms_method in used_MS_methods) {
  for (sample_pattern in sample_name_pattern) {
    # Create name pattern to insert
    pattern <- paste0('*', sample_pattern, '*', ms_method, '*')
    
    # Replace placeholder with name pattern
    ms_peaks_substraction <- gsub('insertNamehere', pattern, ms_peaks_substraction)
    
    # Replace minMass and maxMass with method parameters
    params <- methods_parameters[[ms_method]]
    ms_peaks_substraction <- gsub('insertMinMass', params['minMass'], ms_peaks_substraction)
    ms_peaks_substraction <- gsub('insertMaxMass', params['maxMass'], ms_peaks_substraction)
    
    # Append MS method content to XML content
    xml_content <- c(xml_content, ms_peaks_substraction)
  }
  # Append MS method content to XML content
  xml_content <- c(xml_content, ms_peak_join_pick)
}

# Append the rest of the XML content
xml_content <- c(xml_content, merge)

# Write to file
writeLines(xml_content, output_file_path)
}

# Example usage:
output_file_path <- "output.xml"
used_MS_methods <- c('hpos', 'npos', 'neg')
methods_parameters <- list(hpos = c(minMass = 180, maxMass = 580), npos =  c(minMass = 80, maxMass = 600), neg = c(minMass = 180, maxMass = 580))
sample_name_pattern <- c('Chrom', 'Kim')

