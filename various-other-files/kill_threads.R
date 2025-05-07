
# Function to clean up unnecessary R processes
cleanup_R_processes <- function(process_kill) {
  current_pid <- Sys.getpid()
  
  # List all processes
  all_processes <- ps::ps()
  
  # Filter R processes
  all_such_processes <- all_processes[grepl(process_kill, all_processes$name), ]
  #all_such_processes <- all_processes[is.na(all_processes$name),]
  
  # Filter out the current R process
  processes_to_kill <- all_such_processes[all_such_processes$pid != current_pid, ]
  
  # Kill each unnecessary process
  sapply(processes_to_kill$pid, function(pid) {
    cat("Killing process:", pid, "\n")
    ps::ps_kill(ps::ps_handle(pid))
  })
}

# Call the cleanup function
cleanup_R_processes("R")

