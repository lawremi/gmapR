\name{gsnapBuildRefseqSplicesites}
\alias{gsnapBuildRefseqSplicesites}
\title{Build RefSeq Splice Site IIT for gsnap}
\usage{gsnapBuildRefseqSplicesites(genome)}
\description{add refSeq splice sites to gsnap}
\details{downloads refSeq table from UCSC, parses out splice site donors
and acceptors, creates IIT for gsnap, then copies IIT to genome's
shared gsnap dir}
\value{0}
\author{Cory Barr}
\arguments{\item{genome}{Genome to construct split sites for (only hg19 supported at the moment)}}
