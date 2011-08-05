\name{getGsnapFileCounts}
\alias{getGsnapFileCounts}
\title{Return number of reads in each of gsnap's output files}
\usage{getGsnapFileCounts(aligner_dir, parallelize=TRUE)}
\description{Return number of reads in each of gsnap's output files}
\details{With the --split-output argument, gsnap can output multiple
files. If these files are BAM files, this function counts the
number of reads in each file.}
\value{named list of read counts per gsnap output file}
\author{Cory Barr}
\arguments{\item{aligner_dir}{Directory containing output bam files from gsnap}
\item{parallelize}{Boolean indicating if counting should be done in parallel}}
