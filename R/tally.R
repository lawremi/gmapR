##' Tallies the bases in a BAM alignment by position, strand and
##' cycle. Essentially, this just calls bam_tally -A and stores the
##' result in an IIT.
##'
##' @title Tally bases in a BAM
##' @param bamfiles The input bamfiles
##' @param iitfiles The output iit files, without the .iit suffix
##' @param genome The gsnap genome identifier
##' @param genome_dir The gsnap genome directory
##' @return The paths to the iit files
##' @author Michael Lawrence
##' @importMethodsFrom IRanges as.vector

tally2iit <- function(bamfiles, iitfiles = sub("\\.bam$", "", bamfiles),
                      genome = "hg19_ucsc",
                      genome_dir = "/gnet/is2/data/bioinfo/gmap/data/genomes",
                      positions = c("all", "variant", "both-strands"))
{
  positions <- match(positions, eval(formals(sys.function())$positions)) - 1L
  bam_tally <- paste("bam_tally -C -X", positions, "-T -d", genome,
                     "-D", genome_dir, bamfiles,
                     "| iit_store -o", iitfiles)
  sapply(bam_tally, system)
  structure(paste(iitfiles, "iit", sep = "."), names = iitfiles)
}

gsnapChromosomes <- function(genome, genome_dir) {
  read.table(pipe(paste("get-genome -L -d", genome, "-D", genome_dir)))[,1]
}

##' Reads tally information stored in an IIT generated via sam_tally -A
##'
##' @title Read tally
##' @param x The name of a tally IIT file
##' @param regions The regions of interest, as a
##' \code{\link[IRanges]{RangedData}}
##' @param breaks Integer vector indicating breaks for binning of
##' cycles. Specifying this argument yields columns in the result for
##' the counts in each bin.
##' @param cycleDetails Whether to output a record in the result for
##' every cycle position with a count. This can yield a very large and
##' complex result. If this is \code{TRUE}, the \code{breaks} argument
##' is ignored. Specifying \code{breaks} instead may be a better idea.
##' @param variantsOnly Filter out all reference rows (usually reduces
##' data size dramatically)
##' @param mask A \code{RangedData} or \code{GRanges} with regions to
##' drop from the result. If there is a column named \code{alt}, a
##' tally is only dropped if its read base matches the value in
##' \code{alt}. This is useful for ignoring known variants.
##' @return A \code{GRanges} with counts for each position within a
##' region of interest. Common columns include \code{ref} base,

##' ##' \code{read} base, and \code{count}. If \code{cycleDetails} is
##' \code{TRUE}, the counts are broken down by read \code{cycle} and
##' \code{strand}. If \code{FALSE}, the number of unique cycles
##' (according to cycle number and strand), \code{ncycles}, is
##' returned for each variant. There is also a \code{count.ref} and
##' \code{count.total} for the reference and total counts,
##' respectively. Analogously, \code{ncycles.ref} gives the number of
##' unique cycles matching the reference. If \code{breaks} is given,
##' there are additional columns with the counts for each bin.
##' @author Michael Lawrence
readTally <- function(x, regions, breaks, cycleDetails = FALSE,
                      variantsOnly = TRUE, mask = NULL)
{
  ## annoyingly, the iit_get waits for me to read the result of a query
  ## before processing the next one, so we redirect to a file
  tmpfile <- tempfile(x, "/local")
  message("iit_get from ", x, " to ", tmpfile)
  iit_pipe <- pipe(paste("iit_get -T", x, ">", tmpfile))
  regions <- as(regions, "GRanges")
  regions <- reduce(regions)
  if (!is.null(mask))
    regions <- setdiff(regions, as(mask, "GRanges"))
  query <- paste(seqnames(regions), ":", start(regions), "-", end(regions),
                 sep="")
  writeLines(query, iit_pipe)
  close(iit_pipe)
  message("reading ", tmpfile)
  firstLine <- readLines(tmpfile, 4)
  message(paste("read.table on pipe from ", query))
  if (length(firstLine)>3) {
    tab <- read.table(tmpfile, colClasses = c("character", "integer", "integer",
                                 "character"), sep = "\t")
  } else tab <- data.frame(character(), integer(), integer(), character())
  unlink(tmpfile)
  message("generating GRanges for ", x)
  counts <- strsplit(tab[[4]], " ", fixed=TRUE)
  
  counts_part <- PartitioningByWidth(elementLengths(counts))

  chr <- rep(tab[[1]], width(counts_part))
  pos <- rep(tab[[2]], width(counts_part))
  count.total <- rep(tab[[3]], width(counts_part))
  ## we only get duplicated locations if there is overlap across strands
  ## thus, we force the duplicated locations onto different strands

  location <- paste(tab[[1]], tab[[2]], sep = ":")
  location_dup <- duplicated(location)
  strand <- rep("*", length(location))
  
  strand[match(location[location_dup], location)] <- "+"
  strand[location_dup] <- "-"
  strand <- rep(strand, width(counts_part))
  location <- rep(location, width(counts_part))
  counts_flat <- unlist(counts, use.names=FALSE)
  bases <- sub("(\\D).*", "\\1", counts_flat)
  ref_rows <- start(counts_part)
  ref <- rep(bases[ref_rows], width(counts_part))
  zero_count <- !grepl("(", counts_flat, fixed=TRUE)
### FIXME: obscene hack to simplify code
  counts_flat[zero_count] <- paste(counts_flat[zero_count], "(0@NA)", sep="")
  cycles <- strsplit(sub(".*\\((.*?)\\).*", "\\1", counts_flat), ",", 
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
##  
  if (cycleDetails) {
    strand <- factor(ifelse(cycles_mat[2,] < 0, "-", "+"), 
                     levels=levels(strand()))
    gr <- GRanges(Rle(chr, ncycles), IRanges(rep(pos, ncycles), width=1L),
                  rep(strand, ncycles),
                  location = rep(location, ncycles),
                  ref = DNAStringSet(rep(ref, ncycles)),
                  read = DNAStringSet(rep(bases, ncycles)),
                  cycle = abs(cycles_mat[2,]), strand = strand, 
                  count = cycles_mat[1,])
  } else {
    if (!missing(breaks)) {
### NOTE: overuse of rep() here might lead to over-long vectors            
      cycle_bins <- cut(rep(abs(cycles_mat[2,]), cycles_mat[1,]), breaks)      
      cycle_i <- factor(rep(rep(seq(length(cycles)), ncycles), cycles_mat[1,]),
                        seq(length(cycles)))
      cycle_tab <- table(cycle_i, cycle_bins)
      colnames(cycle_tab) <- paste(head(breaks, -1), tail(breaks, -1),
                                   sep = ".")
      ref_cycle_tab <- cycle_tab[rep(ref_rows, width(counts_part)),]
      colnames(ref_cycle_tab) <- paste(colnames(cycle_tab), "ref", sep = ".")
      cycle_tab <- cbind(cycle_tab, ref_cycle_tab)
    } else cycle_tab <- matrix(nrow = length(pos), ncol = 0)
    base_counts <- as.integer(sub("\\D(\\d+).*", "\\1", counts_flat))
    ncycles[zero_count] <- 0L # otherwise 1 due to our obscene hack above
    ncycles.ref <- rep(ncycles[ref_rows], width(counts_part))
    count.ref <- rep(base_counts[ref_rows], width(counts_part))
    count.pos[zero_count] <- 0L
    count.pos.ref <- rep(count.pos[ref_rows], width(counts_part))
    count.neg[zero_count] <- 0L
    count.neg.ref <- rep(count.neg[ref_rows], width(counts_part))
    gr <- GRanges(chr, IRanges(pos, width=1L), strand, location, 
                  ref = DNAStringSet(ref), read = DNAStringSet(bases),
                  ncycles, ncycles.ref,
                  count = base_counts, count.ref, count.total,
                  count.pos, count.pos.ref, count.neg, count.neg.ref,
                  cycleCount = unclass(cycle_tab))
  }
  if (variantsOnly) # could do this earlier, but let's not optimize yet
    gr <- gr[as.character(values(gr)$ref) != as.character(values(gr)$read)]
  region_ol <- findOverlaps(gr, regions)
  region_strand <- as.vector(strand(regions))[subjectHits(region_ol)]
  strand(gr) <- region_strand
  rc <- region_strand == "-"
  values(gr)$ref[rc] <- reverseComplement(values(gr)$ref[rc])
  values(gr)$read[rc] <- reverseComplement(values(gr)$read[rc])
  values(gr)$location <- paste(values(gr)$location, strand(gr), sep = ":")
  message("finished ", x)
  gr
}

##' tally2GR
##' uses Tom's bam_tally to produce variant calls from a BAM file


##' @export
tally2GR<- function(bamfiles,
                    genome = "hg19_ucsc",
                    genome_dir = "/gnet/is2/data/bioinfo/gmap/data/genomes",
                    chr_gr = getChrGRangesFromGmap(genome, genome_dir),
                    regions,
                    variant_strand=1,
                    min_cov =1,
                    breaks=NULL,
                    bqual_thresh=0,
                    map_qual=0,
                    mc.cores =1)
{
  chr_ids <-paste(seqnames(chr_gr), ":", start(chr_gr), "-", end(chr_gr), sep = "")
  chr_ids <- as.list(as.character(chr_ids))
  has_regions <- !missing(regions)
  list_of_gr <- mclapply(chr_ids, mc.cores = mc.cores, function(chr_name){
    tmp <- tempfile(file = chr_name)
    on.exit(unlink(tmp))
    tally <- pipe(paste("bam_tally --block-format=0 --cycles --quality-scores -q", map_qual,
                        " --variants=", variant_strand, " --min-depth=", min_cov,
                        " --totals --db=", genome, " --dir=", genome_dir, " ", bamfiles,
                        paste(" '", chr_name, "' 2> /dev/null; echo $? >", tmp, sep = ""), sep =""))
    tab <- read.table(tally, colClasses = c("character", "integer", "integer",
                               "character"), sep = "\t",
                      col.names = c("chrom", "position", "count", "cycles"))

    if(readLines(tmp) == '139'){
      cat(paste("bam_tally returned Segmentation fault on chromosome:", chr_name, "\n", sep =""))
      stop(paste("bam_tally returned Segmentation fault on chromosome:", chr_name))
    }
    if(dim(tab)[1] <1){
      gr <- GRanges()
      message(paste("finished chr ", chr_name))
      seqlevels(gr) <- seqlevels(chr_gr)
      return(gr)
    } else {      
      gr <- GRanges()
      counts <- strsplit(tab[[4]], " ", fixed=TRUE)
      lens <- unlist(lapply(counts, function(x) length(x)))
      invar <- counts[lens==1]
      inv_tab <- tab[lens==1,]
      counts <- counts[lens>1]
      var_tab <- tab[lens>1,]
      if(length(counts)>=1){
        gr <-gmapR:::.processReads(counts,
                                   var_tab,
                                   type = "var",
                                   breaks=breaks,
                                   bqual_thresh=bqual_thresh)
      }
      if(length(invar)>=1){
        gr2<- gmapR:::.processReads(invar,
                                    inv_tab,
                                    type = "invar",
                                    breaks=breaks,
                                    bqual_thresh=bqual_thresh)
        gr<- c(gr, gr2)
        gr<- sort(gr)
      }
      
      if (has_regions) {
          region_ol <- findOverlaps(gr, regions)
          region_strand <- as.vector(strand(regions))[subjectHits(region_ol)]
          strand(gr) <- region_strand
          rc <- region_strand == "-"
          values(gr)$ref[rc] <- reverseComplement(values(gr)$ref[rc])
          values(gr)$read[rc] <- reverseComplement(values(gr)$read[rc])
        }
      values(gr)$location <- paste(values(gr)$location, strand(gr), sep = ":")
    }
    message(paste("finished chr ", chr_name))
    seqlevels(gr) <- seqlevels(chr_gr)
    return(gr)
  })
  GR_full <- do.call(c, list_of_gr)
  return(GR_full)
}


      
.processReads <-function(counts, tab, type="var", breaks, bqual_thresh){ 
  counts_part <- PartitioningByWidth(elementLengths(counts))
  chr <- rep(tab[[1]], width(counts_part))
  pos <- rep(tab[[2]], width(counts_part))
  count.total <- rep(tab[[3]], width(counts_part))
  ## we only get duplicated locations if there is overlap across strands
  ## thus, we force the duplicated locations onto different strands
  
  location <- paste(tab[[1]], tab[[2]], sep = ":")
  location_dup <- duplicated(location)
  strand <- rep("*", length(location))
  
  strand[match(location[location_dup], location)] <- "+"
  strand[location_dup] <- "-"
  strand <- rep(strand, width(counts_part))
  location <- rep(location, width(counts_part))
  counts_flat <- unlist(counts, use.names=FALSE)
  bases <- sub("(\\D).*", "\\1", counts_flat)
  ref_rows <- start(counts_part)
  ref <- rep(bases[ref_rows], width(counts_part))
  zero_count <- !grepl("(", counts_flat, fixed=TRUE)
### FIXME: obscene hack to simplify code
  counts_flat[zero_count] <- paste(counts_flat[zero_count], "(0@NA)(0Q-1)", sep="")
  cycles <- strsplit(sub(".*\\((.*?)\\).*", "\\1", counts_flat), ",", 
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
  ##message("working on breaks") 
  if (!is.null(breaks)) {
    ## NOTE: overuse of rep() here might lead to over-long vectors            
    cycle_bins <- cut(rep(abs(cycles_mat[2,]), cycles_mat[1,]), breaks)      
    cycle_i <- factor(rep(rep(seq(length(cycles)), ncycles), cycles_mat[1,]),
                      seq(length(cycles)))
    cycle_tab <- table(cycle_i, cycle_bins)
    colnames(cycle_tab) <- paste(head(breaks, -1), tail(breaks, -1),
                                 sep = ".")
    ref_cycle_tab <- cycle_tab[rep(ref_rows, width(counts_part)),]
    colnames(ref_cycle_tab) <- paste(colnames(cycle_tab), "ref", sep = ".")
    cycle_tab <- cbind(cycle_tab, ref_cycle_tab)
  } else  cycle_tab <- matrix(nrow = length(pos), ncol = 0)
  ##message("completed breaks")
  base_counts <- as.integer(sub("\\D(\\d+).*", "\\1", counts_flat))
  ncycles[zero_count] <- 0L # otherwise 1 due to our obscene hack above
  ncycles.ref <- rep(ncycles[ref_rows], width(counts_part))
  count.ref <- rep(base_counts[ref_rows], width(counts_part))
  count.pos[zero_count] <- 0L
  count.pos.ref <- rep(count.pos[ref_rows], width(counts_part))
  count.neg[zero_count] <- 0L
  count.neg.ref <- rep(count.neg[ref_rows], width(counts_part))
  ##parcing the quality info from the tally out put
  quals <- strsplit(sub(".*\\)\\((.*?)\\).*", "\\1", counts_flat), ",", 
                    fixed=TRUE)
  nquals <- elementLengths(quals)
  
  count_above_thresh <- lapply(quals, function(x){
    mat <-matrix(unlist(strsplit(x, "Q", fixed=TRUE), 
                        use.names=FALSE), nrow=2)
    logical <- (as.numeric(mat[2,]))>bqual_thresh
    sum(as.numeric(mat[1,logical]), na.rm=T)
  })
  mean_qual <- lapply(quals, function(x){
    mat <-matrix(unlist(strsplit(x, "Q", fixed=TRUE),
                        use.names=FALSE), nrow=2)
    logical <- (as.numeric(mat[2,]))>bqual_thresh
    sum((as.numeric(mat[2,logical])*as.numeric(mat[1,logical])), na.rm=T) / sum(as.numeric(mat[1,logical]), na.rm=T)
  }) 
  high_qual <- unlist(count_above_thresh)      
  high_qual_ref <- rep(high_qual[ref_rows], width(counts_part))
  mean_qual <- unlist(mean_qual)
  mean_qual_ref <- rep(mean_qual[ref_rows], width(counts_part))
  
  gr <- GRanges(chr, IRanges(pos, width=1L), strand, location, 
                ref = DNAStringSet(ref), read = DNAStringSet(bases),
                ncycles, ncycles.ref,
                count = base_counts, count.ref, count.total,
                high.quality = high_qual,high.quality.ref=high_qual_ref,
                mean.quality = mean_qual,mean.quality.ref=mean_qual_ref,
                count.pos, count.pos.ref, count.neg, count.neg.ref,
                cycleCount = unclass(cycle_tab))
  
  ##adding ref counts back in as counts per break segment
  ##there seems to be a bug somewhere above that is causing the ref positions to be listed twice some times. Removing these for now
  if(type=="var"){
    gr <- gr[as.character(values(gr)$ref) != as.character(values(gr)$read)]
  }else{
    values(gr)$ncycles <- NA
    values(gr)$count <- NA
    values(gr)$high.quality <- NA
    values(gr)$mean.quality <- NA
    values(gr)$count.pos <- NA
    values(gr)$count.neg <- NA
  }
  return(gr)
}
