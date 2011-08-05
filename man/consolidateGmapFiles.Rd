\name{consolidateGmapFiles}
\alias{consolidateGmapFiles}
\title{Consolidate gmap's output files}
\usage{consolidateGmapFiles(sam_file_dir, paired_end=FALSE,
    remove_merged=FALSE)}
\description{Consolidates all pieces from parallelized gsnap output into
appropriately merged files}
\details{If gmap was run in single-end mode, this will consolidate all
parallelized into the 3 output files. If gmap was run in
paired-end mode, consolidates to the 7 output files.}
\value{list of the names of the files created from consolidation}
\author{Cory Barr}
\arguments{\item{sam_file_dir}{directory where gmap's sam file output is stored}
\item{paired_end}{indicated whether gmap was run in paired_end mode}
\item{remove_merged}{remove the individual pieces after they are merged}}
