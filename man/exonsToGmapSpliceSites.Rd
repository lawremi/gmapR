\name{exonsToGmapSpliceSites}
\alias{exonsToGmapSpliceSites}
\title{Convert a data frame of exon info into input for incorporating known splice sites into Gmap/Gsnap}
\usage{exonsToGmapSpliceSites(exons)}
\description{Convert a data frame of exon information into input to build a splice site index for Gmap/Gsnap}
\details{the input data frame needs to have columns named "chrom,"
"strand," "exonStarts", "exonEnds", and "name."  A common way to
obtain this is from UCSC via rtracklayer:

library(rtracklayer)
session <- browserSession()
genome(session) <- genome
query <- ucscTableQuery(session, "refGene")
refGeneTable <- getTable(query)

To create a splice-site index for Gmap/Gsnap, take the following
output (here called spliceSites) and pass it to buildGmapSpliceSites:

spliceSites <- exonsToSpliceSites(refGeneTable)}
\value{a character vector}
\author{Cory Barr}
\arguments{\item{exons}{A data frame containing the fields described in the details section.}}
