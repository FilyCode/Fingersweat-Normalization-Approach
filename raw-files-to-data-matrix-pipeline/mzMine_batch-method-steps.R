# A file that contains all the steps for the mzMine batch method


get_mzMine_batch_method_steps <- function() {
  beginning <- c()
  
  ms_peaks_substraction <- c()
  
  ms_peak_join_pick <- c()
  
  merge <- c()
  
  return(c(beginning = beginning, ms_peaks_substraction = ms_peaks_substraction, ms_peak_join_pick = ms_peak_join_pick, merge = merge))
}