\name{consolidateGsnapFiles}
\alias{consolidateGsnapFiles}
\title{Consolidate gsnap's parallelized output files}
\usage{consolidateGsnapFiles(sam_file_dir, remove_merged=FALSE,
    parallelized=TRUE)}
\description{Consolidates all pieces from parallelized gsnap output into
appropriately merged files}
\details{If gsnap was run in single-end mode, this will consolidate all
parallelized into the 3 output files. If gsnap was run in
paired-end mode, consolidates to the 7 output files.}
\value{list of the names of the files created from consolidation}
\author{Cory Barr}
\arguments{\item{sam_file_dir}{directory where gsnap's sam file output is stored}
\item{remove_merged}{remove the individual pieces after they are merged}}
