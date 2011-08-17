\name{buildGmapIITFromFasta}
\alias{buildGmapIITFromFasta}
\title{Build GMAP/GSNAP Index from a FASTA file}
\usage{buildGmapIITFromFasta(genome, fasta, gmap_data_dir)}
\description{Given a FASTA file, will make the FASTA a reference sequence for GMAP/GSNAP}
\value{0 if successful}
\author{Cory Barr}
\arguments{\item{genome}{the name of the genome you want to pass to gsnap to align against your FASTA}
\item{fasta}{A fasta of the reference sequence you wish to align against}
\item{gmap_data_dir}{The location to save the GMAP/GSNAP index into}}
