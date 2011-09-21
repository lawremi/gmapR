\name{buildGmapIndex}
\alias{buildGmapIndex}
\title{Create Gmap Reference Genome Index}
\usage{buildGmapIndex(genome, gmap_data_dir)}
\description{Gsnap index of hg19 (with haplotypes and unassembled regions excluded)}
\details{a reference genome into gmap/gsnap by downloading a
genome and constructing an IIT. Currently only hg19 is supported.}
\value{0}
\author{Cory Barr}
\arguments{\item{genome}{The name of the genome, often an assembly (e.g., hg19)}
\item{gmap_data_dir}{The directory to store the Gmap IIT}}
