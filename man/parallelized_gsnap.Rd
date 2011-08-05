\name{parallelized_gsnap}
\alias{parallelized_gsnap}
\title{Parallelized gsnap}
\usage{parallelized_gsnap(num_machines, procs_per_machine, gsnap,
    gsnap_params, input_file, input_file_2, paired_end=F, output_dir,
    test_mode=FALSE, intern=FALSE, multifile_out=TRUE,
    record_sys_call_dir=NA)}
\description{This function parallelizes gsnap across a cluster.}
\value{a list of the return values of the system function for
each parallelized piece}
\author{Cory Barr}
\arguments{\item{num_machines}{number of cluster machines to run gsnap across}
\item{procs_per_machine}{number of gsnap jobs to run on each cluster node}
\item{gsnap}{location of gsnap executable}
\item{gsnap_params}{a string of the parameters to pass to gsnap. Note "--split-output" is not allowed.}
\item{input_file}{location of the input fastq to align}
\item{input_file_2}{location of the second input fastq to align (for paired ends)}
\item{paired_end}{logical indicating if gsnap should be run in paired-end mode}
\item{output_dir}{directory for gsnap to write its output}
\item{test_mode}{used for debugging. turns off parallelization.}
\item{intern}{indicates whether to have STDOUT the return
value. identical to 'intern' argument to the system() function}
\item{multifile_out}{tells gsnap to write to 3 output files for
single-end data, or 7 output files for paired-end data}}
