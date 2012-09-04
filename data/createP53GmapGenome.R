##This script was made to create a small GMAP genome to use for
##development and testing. Later, p53 was used.

##How to get the p53 sequence
##p53 location obtained via UCSC link to NM_000546: http://genome.ucsc.edu/cgi-bin/hgTracks?hgHubConnect.destUrl=..%2Fcgi-bin%2FhgTracks&clade=mammal&org=Human&db=hg19&position=tp53&hgt.positionInput=tp53&hgt.suggestTrack=knownGene&Submit=submit&hgsid=296149047
##chr17:7,571,720-7,590,868
##added 1 megabase to each side:

p53Genome <- local({library("Biostrings")
                    library("BSgenome.Hsapiens.UCSC.hg19")
                    p53Seq <- getSeq(x=Hsapiens, names="chr17",
                                     start=7571720 - 1000000,
                                     end=7590868 + 1000000,
                                     strand="-")
                    p53Seq <- as(p53Seq, "DNAStringSet")
                    names(p53Seq) <- "p53Genome"
                    p53Genome <- GmapGenome(genome=p53Seq,
                                            name="p53Genome", create=TRUE)
                    return(p53Genome)
                  })
