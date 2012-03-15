##' title tallyInDel uses bam_tally to get the indels from a bam file
##' param bamfiles bam file to be tallied 
##' param genome gsnap genome that was used to build the bam file
##' param genome_dir directory the genome lives in
##' param genomeGR A GRanges object with the chromosome information 
##' param variant_strand how many strands should the variant be seen on
##' param min_cov what is the minimum depth to look at
##' param map_qual minimum map quality for reads
##' param mc.cores how many cores to run the mclapply over
##' return A InDel variant GRanges
##' author Jeremiah Degenhardt
##'

tallyInDel <- function(bamfiles,
                       genome = "hg19_ucsc",
                       genome_dir = "/gnet/is2/data/bioinfo/gmap/data/genomes",
                       genomeGR = getChrGRangesFromGmap(genome=genome, genome_dir=genome_dir),
                       variant_strand=1,
                       min_cov =1,
                       map_qual=0,
                       mc.cores =1)
{
  options(scipen=500)
  chr_ga <-  readBamGappedAlignments(bam_file)
  chr_grl <- grglist(chr_ga, drop.D.ranges=TRUE)
  chr_ids <- as.list(as.character(seqnames(genomeGR)))
  list_of_gr <- mclapply(chr_ids, mc.cores = mc.cores, function(chr_name){
    tmp <- tempfile(file = chr_name)
    on.exit(unlink(tmp))
    tally <- pipe(paste("bam_tally --block-format=0 --cycles --indels -q", map_qual,
                        " --variants=", variant_strand, " --min-depth=", min_cov,
                        " --db=", genome, " --dir=", genome_dir, " ", bamfiles,
                        paste(" '", chr_name, ":' 2> /dev/null; echo $? >", tmp, sep = ""), sep =""))
    tab <- read.table(tally, colClasses = c("character", "integer",
                               "character"), sep = "\t",
                      col.names = c("chrom", "position", "cycles"))
    if(readLines(tmp) == '139'){
      cat(paste("bam_tally returned Segmentation fault on chromosome:", chr_name, "\n", sep =""))
      stop(paste("bam_tally returned Segmentation fault on chromosome:", chr_name))
    }
    if(dim(tab)[1] <1){
      gr <- GRanges()
    } else { 
      
      ins <- tab[grepl(tab[,1], pattern = "^\\^"),]
      del <- tab[grepl(tab[,1], pattern = "^_"),]
    
    
#######################
### Process insertions
#######################
      if(dim(ins)[1]<1){
        ins_gr <- GRanges()
      }else{
        counts <- strsplit(ins[[3]], " ", fixed=TRUE)
        counts_part <- PartitioningByWidth(elementLengths(counts))
        chr <- rep(substr(ins[[1]], 2, nchar(ins[[1]])), width(counts_part)/2)
        pos <- rep(ins[[2]], width(counts_part)/2)
        ##cov_base <- as.integer(cov[[chrom]][pos])
        location <- paste(chr, pos, sep = ":")
        location_dup <- duplicated(location)
        strand <- rep("*", length(location))
        counts_flat <- unlist(counts, use.names=FALSE)
        groups <- as.data.frame(cbind(location[togroup(counts)], counts_flat))
        base <- groups[seq(1,dim(groups)[1], by = 2),]
        info <- groups[seq(2,dim(groups)[1], by =2),]
        cycles <- strsplit(sub(".*\\((.*?)\\).*", "\\1", info[,2]), ",",
                           fixed=TRUE)
        ncycles <- elementLengths(cycles)
        cycles_flat <- unlist(cycles, use.names=FALSE)
        cycles_counts <- unlist(strsplit(cycles_flat, "@", fixed=TRUE),
                                use.names=FALSE)
        cycles_mat <- matrix(suppressWarnings(as.integer(cycles_counts)), nrow=2)
        ## insane code for counting up reads by strand
        neg_cycle <- cycles_mat[2,] < 0
        coord_ind <- rep(seq_along(cycles), ncycles)
        strand_key <- paste(coord_ind, neg_cycle)
        strand_key_uniq <- !duplicated(strand_key)
        strand_key <- factor(strand_key, strand_key[strand_key_uniq])
        neg_cycle_uniq <- neg_cycle[strand_key_uniq]
        coord_ind_uniq <- coord_ind[strand_key_uniq]
        strand_count <- rowsum(cycles_mat[1,], strand_key)
        count.pos <- count.neg <- rep(0L, length(cycles))
        ## have to use which() here to drop the NA values
        count.pos[coord_ind_uniq[which(!neg_cycle_uniq)]] <-
          strand_count[which(!neg_cycle_uniq)]
        count.neg[coord_ind_uniq[which(neg_cycle_uniq)]] <-
          strand_count[which(neg_cycle_uniq)]
        base_counts <- as.integer(sub("(\\d+).*", "\\1", info[,2]))
        over <- GRanges(seqnames = chr, IRanges(pos, pos+1))
        total_cov<-countOverlaps(query=over, subject=chr_grl, type = "within")
        count.ref <- total_cov - base_counts
        
        ## Getting the ref base for the position to the left of the event
        tmp_pos <- tempfile(file = paste(chr_name, "positions", sep = ""))
        on.exit(unlink(tmp_pos))
        writeLines(location, con = tmp_pos)
        command <- pipe(paste("cat ", tmp_pos, " | get-genome -d ", genome, " -D ", genome_dir, sep = "")) 
        ref_base <- read.table(command, sep = "\t", colClasses='character')
        ref_base <- ref_base[seq(2, length(ref_base[,1]), by = 2),]
        
        ## making the granges
        
        ins_gr <- GRanges(chr, IRanges(pos, width=1L), strand, location,
                          ref = DNAStringSet(ref_base),
                          read = DNAStringSet(paste(ref_base, as.character(base[,2]), sep = "")),
                          ncycles,
                          count = base_counts, count.ref, count.total=total_cov,
                          count.pos, count.neg, type = "insertion")
      }
      
#######################
### Process deletions
#######################
      if(dim(del)[1] <1){
        del_gr <- GRanges()
      }else{
        counts <- strsplit(del[[3]], " ", fixed=TRUE)
        counts_part <- PartitioningByWidth(elementLengths(counts))
        chr <- rep(substr(del[[1]], 2,nchar(del[[1]])), width(counts_part)/2)
        pos <- rep(del[[2]], width(counts_part)/2)
        pos <- pos-1
        ##cov_base <- as.integer(cov[[chrom]][pos])
        location <- paste(chr, pos, sep = ":")
        location_dup <- duplicated(location)
        strand <- rep("*", length(location))
        counts_flat <- unlist(counts, use.names=FALSE)
        groups <- as.data.frame(cbind(location[togroup(counts)], counts_flat))
        base <- groups[seq(1,dim(groups)[1], by = 2),]
        info <- groups[seq(2,dim(groups)[1], by =2),]
        cycles <- strsplit(sub(".*\\((.*?)\\).*", "\\1", info[,2]), ",",
                           fixed=TRUE)
        ncycles <- elementLengths(cycles)
        cycles_flat <- unlist(cycles, use.names=FALSE)
        cycles_counts <- unlist(strsplit(cycles_flat, "@", fixed=TRUE),
                                use.names=FALSE)
        cycles_mat <- matrix(suppressWarnings(as.integer(cycles_counts)), nrow=2)
        ## insane code for counting up reads by strand
        neg_cycle <- cycles_mat[2,] < 0
        coord_ind <- rep(seq_along(cycles), ncycles)
        strand_key <- paste(coord_ind, neg_cycle)
        strand_key_uniq <- !duplicated(strand_key)
        strand_key <- factor(strand_key, strand_key[strand_key_uniq])
        neg_cycle_uniq <- neg_cycle[strand_key_uniq]
        coord_ind_uniq <- coord_ind[strand_key_uniq]
        strand_count <- rowsum(cycles_mat[1,], strand_key)
        count.pos <- count.neg <- rep(0L, length(cycles))
        ## have to use which() here to drop the NA values
        count.pos[coord_ind_uniq[which(!neg_cycle_uniq)]] <-
          strand_count[which(!neg_cycle_uniq)]
        count.neg[coord_ind_uniq[which(neg_cycle_uniq)]] <-
          strand_count[which(neg_cycle_uniq)]
        base_counts <- as.integer(sub("(\\d+).*", "\\1", info[,2]))
        over <- GRanges(seqnames = chr, IRanges(pos, width=nchar(as.character(base[,2]))))
        count.ref<-countOverlaps(query=over, subject=chr_grl)
        total_cov <- count.ref + base_counts
        
        tmp_pos <- tempfile(file = paste(chr_name, "positions", sep = ""))
        on.exit(unlink(tmp_pos))
        writeLines(location, con = tmp_pos)
        command <- pipe(paste("cat ", tmp_pos, " | get-genome -d ", genome, " -D ", genome_dir, sep = ""))
        ref_base <- read.table(command, sep = "\t", colClasses='character')
        ref_base <- ref_base[seq(2, length(ref_base[,1]), by = 2),] 
        
        
        
        del_gr <- GRanges(chr, IRanges(pos, width=1L), strand, location,
                          ref = DNAStringSet(paste(ref_base, as.character(base[,2]), sep = "")),
                          read = DNAStringSet(ref_base),
                          ncycles,
                          count = base_counts, count.ref, count.total=total_cov,
                          count.pos, count.neg, type = "deletion")
      }
######################################
### combine ins and dels and return gr
######################################
      
      gr <- c(ins_gr, del_gr)
      gr <- sort(gr)
    }
    seqlevels(gr) <- unlist(chr_ids)
    
    return(gr)
    
  })
  GR_full <- do.call(c, list_of_gr)
  seqlevels(GR_full) <- seqlevels(genomeGR)
  seqlengths(GR_full) <- seqlengths(genomeGR)
  return(GR_full)
}
                         
