\name{consolidateSAMFiles}
\alias{consolidateSAMFiles}
\title{Consolidate parallelized gsnap SAM files}
\usage{consolidateSAMFiles(sam_files, outfile, remove_merged)}
\description{Consolidate SAM files created by gsnap}
\details{This function takes uses samtools to sort and merge SAM files into a single BAM file}
\value{character vector indicating the name of the resulting BAM file}
\author{Cory Barr}
\arguments{\item{sam_files}{character vector of the SAM files}
\item{outfile}{the name of the output BAM file. Will be place in the
same directory as the first SAM file}
\item{remove_merged}{logical indication whether to delete the
contents of sam_files after merging}}
