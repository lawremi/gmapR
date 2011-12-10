##' This function parallelizes gsnap across a cluster.
##'
##' @title Parallelized gsnap
##' @param num_machines number of cluster machines to run gsnap across
##' @param procs_per_machine number of gsnap jobs to run on each cluster node
##' @param gsnap location of gsnap executable
##' @param gsnap_params a string of the parameters to pass to gsnap. Note "--split-output" is not allowed.
##' @param input_file location of the input fastq to align
##' @param input_file_2 location of the second input fastq to align (for paired ends)
##' @param paired_end logical indicating if gsnap should be run in paired-end mode
##' @param output_dir directory for gsnap to write its output
##' @param test_mode used for debugging. turns off parallelization.
##' @param intern indicates whether to have STDOUT the return
##' value. identical to 'intern' argument to the system() function
##' @param multifile_out tells gsnap to write to 3 output files for
##' single-end data, or 7 output files for paired-end data
##' @param record_sys_call_dir If a directory is supplied, will record the system call made to gsnap here
##' @return a list of the return values of the system function for
##' each parallelized piece
##' @author Cory Barr
##' @export
parallelized_gsnap <- function(num_machines,
                               procs_per_machine,
                               gsnap,
                               gsnap_params,
                               input_file,
                               input_file_2=NULL,
                               paired_end=F,
                               output_dir,
                               test_mode=FALSE,
                               intern=FALSE,
                               multifile_out=TRUE,
                               record_sys_call_dir=NA) {
  
  if (length(grep("split-output ", gsnap_params) > 0))
    stop("This function cannot handle gsnap parameters with the --split-output switch enabled")

  if(procs_per_machine == 1) {
    apply_func <- lapply
  } else {
    apply_func <- mclapply
  }
  
  .gsnap_per_machine <- function(i,
                                 procs_per_machine,
                                 gsnap,
                                 gsnap_params,
                                 paired_end,
                                 input_file,
                                 input_file_2,
                                 output_dir,
                                 total_gsnap_parts,
                                 intern,
                                 multifile_out,
                                 record_sys_call_dir) {

    ##library calls needed so slaves have access
    library("multicore")
    
    gsnap_parts <-
      (i - 1) * procs_per_machine + (seq_len(procs_per_machine) - 1)
    gsnap_part_params <- paste("--part=",
                               gsnap_parts,
                               "/",
                               total_gsnap_parts,
                               sep="")

    if (paired_end)
      gsnap_params <- paste(gsnap_params, "-a paired")
    
    ##make the prefix of gsnaps 'split-output' outs be unique
    if (multifile_out) {
      gsnap_multi_out_params <- paste("gsnap_out.", gsnap_parts, sep="")
      gsnap_multi_out_params <- paste(output_dir, gsnap_multi_out_params, sep="/")
      sys_commands <- paste(gsnap,
                            gsnap_params,
                            gsnap_part_params,
                            "--split-output", gsnap_multi_out_params,
                            input_file)
    } else {
      stop("Not tested when not using gsnap's --split-output flag")
    }
    
    if (paired_end) {
      sys_commands <- paste(sys_commands, input_file_2)
    }
    
    sys_commands <- paste("nice", sys_commands)
    sys_commands <- paste(sys_commands, "2>&1")

    if(!is.na(record_sys_call_dir)) {
      if(!(file.exists(record_sys_call_dir)))
        dir.create(record_sys_call_dir)

      sys_call_file <- file.path(record_sys_call_dir,
                                 paste("aligner.sys_call",
                                       i - 1,
                                       sep="."))
      writeLines(sys_commands,
                 con=sys_call_file)
    }
    
    ##TODO: this awkward bit is about to 
    if (intern) { ##parsing gsnap output. also sending multimaps out through STDOUT
      stop("intern=TRUE is not working yet")
    } else { #do not alter gsnap output  
      results <- apply_func(sys_commands, system, intern=TRUE)
    }

    ##check results to see if each process finished
    successful <- sapply(seq_len(length(results)),
                         function(index) {
                           returned_lines <- results[[index]]
                           res_len <- length(returned_lines)
                           grepl("^Processed \\d+ queries in [\\d\\.]+ seconds",
                                 returned_lines[res_len],
                                 perl=TRUE)
                         })
    if(!all(successful))
      stop(paste("aligner did not return successfully. System call is "),
           head(sys_commands[successful == FALSE], n=1))
        
    return(as.integer(successful) - 1)
  }
  
  if (test_mode) {
    results <- lapply(seq_len(num_machines),
                      .gsnap_per_machine,
                      procs_per_machine=procs_per_machine,
                      gsnap=gsnap,
                      gsnap_params=gsnap_params,
                      paired_end=paired_end,
                      input_file=input_file,
                      input_file_2=input_file_2,
                      output_dir=output_dir,
                      total_gsnap_parts=num_machines * procs_per_machine,
                      intern=intern,
                      multifile_out=multifile_out,
                      record_sys_call_dir=record_sys_call_dir)
  } else {
    if (num_machines > 1) {
      library("snow")
      cl <- makeCluster(num_machines, "MPI")
      results <- parLapply(cl,
                           seq_len(num_machines),
                           .gsnap_per_machine,
                           procs_per_machine=procs_per_machine,
                           gsnap=gsnap,
                           gsnap_params=gsnap_params,
                           paired_end=paired_end,
                           input_file=input_file,
                           input_file_2=input_file_2,
                           output_dir=output_dir,
                           total_gsnap_parts=num_machines * procs_per_machine,
                           intern=intern,
                           multifile_out=multifile_out,
                           record_sys_call_dir=record_sys_call_dir)
      stopCluster(cl)
    } else {      
      results <- apply_func(seq_len(num_machines),
                            .gsnap_per_machine,
                            procs_per_machine=procs_per_machine,
                            gsnap=gsnap,
                            gsnap_params=gsnap_params,
                            paired_end=paired_end,
                            input_file=input_file,
                            input_file_2=input_file_2,
                            output_dir=output_dir,
                            total_gsnap_parts=num_machines * procs_per_machine,
                            intern=intern,
                            multifile_out=multifile_out,
                            record_sys_call_dir=record_sys_call_dir)
    }
  }  
  
  return(results)
}
