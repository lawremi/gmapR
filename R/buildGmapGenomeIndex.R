##' Gsnap index of hg19 (with haplotypes and unassembled regions excluded)
##' 
##'Incorporate a reference genome into gmap/gsnap by downloading a
##' genome and constructing an IIT. Currently only hg19 is supported.
##' @title Create Gmap Reference Genome Index
##' @param genome 
##' @return 0
##' @author Cory Barr
##' @export
buildGmapIndex <- function(genome) {
  
  retriveGenomeFasta <- function(genome) {
    if(genome != 'hg19')
      stop("only hg19 supported at the moment")

    tmp_dir <- tempdir()
    if(file.exists(tmp_dir))
      unlink(tmp_dir, recursive=TRUE)
    dir.create(tmp_dir)

    on.exit(setwd(file_path_as_absolute(getwd())))
    setwd(tmp_dir)
    sys_command <- paste('wget http://hgdownload.cse.ucsc.edu/goldenPath/',
                         genome,
                         '/bigZips/chromFa.tar.gz',
                         sep='')
    system(sys_command)
    system('tar xzf chromFa.tar.gz')

    genome_fasta_file <- paste(genome, 'fa', sep='.')
    ##excluding haplotypes, cat genome seqs together into one fasta:
    sys_call <- paste('ls *fa | grep -v "_gl" | grep -v "_hap" | xargs -i cat {} > ',
                      genome_fasta_file,
                      sep='')
    system(sys_call)
    
    file_path_as_absolute(genome_fasta_file)     
  }

  genome_fa <- retrieveGenomeFasta(genome)

  buildGmapIITFromFasta <- function(genome, fasta) {

    ##currently, hg19 for the Genie pipelines is names hg19_ucsc to
    ##distinguish origin. This will be phased out
    genome <- 'hg19_ucsc'
    
    gmap_setup <- file.path(globals()['gmap_bin'],
                              "gmap_setup")  
    sys_command <- paste(gmap_setup,
                         "-d",
                         genome,
                         fasta)
    system(sys_command)
      
    sys_command <- paste("make -f",
                         paste('Makefile.', genome, sep=''),
                         "coords")
    system(sys_command)
      
    sys_command <- paste("make -f",
                         paste('Makefile.', genome, sep=''),
                         "gmapdb")
    system(sys_command)
      
    sys_command <- paste("make -f",
                         paste('Makefile.', genome, sep=''),
                         "install")
    system(sys_command)
  }
  
  buildGmapIITFromFasta(genome, genome_fa)
  return(0)
}
