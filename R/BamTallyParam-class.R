### =========================================================================
### BamTallyParam class
### -------------------------------------------------------------------------
###
### Collects parameters for the bam_tally command
###


setClass("BamTallyParam",
         representation(genome = "GmapGenome",
                        which = "GenomicRanges",
                        desired_read_group = "characterORNULL",
                        minimum_mapq = "integer",
                        concordant_only = "logical",
                        unique_only = "logical",
                        primary_only = "logical",
                        ignore_duplicates = "logical",
                        min_depth = "integer",
                        variant_strand = "integer",
                        ignore_query_Ns = "logical",
                        indels = "logical",
                        include_soft_clips = "integer",
                        exon_iit = "characterORNULL",
                        xs = "logical",
                        read_pos = "logical",
                        min_base_quality = "integer",
                        noncovered = "logical",
                        nm = "logical"))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

normArgWhich <- function(x, genome) {
  if (!is(x, "GenomicRanges"))
    stop("'which' must be a GenomicRanges")
  si <- seqinfo(genome)
  seqinfo(x, new2old = match(seqlevels(si), seqlevels(x))) <-
    merge(si, seqinfo(x))
  x
}

setGeneric("normArgCdsIIT", function(exon_iit, genome, ...)
           standardGeneric("normArgCdsIIT"))

setMethod("normArgCdsIIT", "ANY", function(exon_iit, genome, BPPARAM) {
    if(!is.null(exon_iit)
       && (!is(exon_iit, "character") || length(exon_iit) > 1))
      stop("Invalid exon_iit value. See ?BamTallyParam for acceptable values")
    if(!is.null(exon_iit) && (nchar(exon_iit) && !file.exists(exon_iit)))
        stop("exon_iit is non-empty and does not point to an existing file")
    exon_iit

})

setMethod("normArgCdsIIT", "TxDb", function(exon_iit, genome, BPPARAM) {
    exon_iit = iit_store(exonsBy(exon_iit, "gene"),
      dest = tempfile(pattern = "cds", fileext = ".iit"), BPPARAM = BPPARAM)
    normArgCdsIIT(exon_iit)
})

setMethod("normArgCdsIIT", "GRangesList", function(exon_iit, genome, BPPARAM) {
    exon_iit = iit_store(exon_iit, dest = tempfile(pattern = "cds",
                                     fileext = ".iit"), BPPARAM = BPPARAM)
    normArgCdsIIT(exon_iit)
})

setMethod("normArgCdsIIT", "GmapGenome", function(exon_iit, genome, BPPARAM) {
    if(!identical(exon_iit, genome))
        stop("You cannot specify different GmapGenomes for exon_iit and genome")
    exon_iit = list.files(mapsDirectory(exon_iit), pattern = "genes.iit")
    if(!length(exon_iit))
        stop("No map matching the pattern 'genes.iit' found for this GmapGenome")
    if(length(exon_iit) > 1)
        stop("Multipe map matching the pattern 'genes.iit' found for this GmapGenome")
    normArgCdsIIT(exon_iit)
})

.makeCdsIIT <- function(exons,
                        filename = tempfile(pattern="exon_iit", fileext=".iit"))
{
    exonsFlat <- unlist(exons, use.names=FALSE)
    exonsPart <- PartitioningByWidth(exons)
    exonsHead <- exonsFlat[-end(exonsPart)]
    donors <- flank(exonsHead, 1L, start = FALSE)
    exonsTail <- exonsFlat[-start(exonsPart)]
    acceptors <- flank(exonsTail, 1L, start = TRUE)
    sites <- c(resize(donors, 2L, fix = "end"),
               resize(acceptors, 2L, fix = "start"))
    names(sites) <- values(sites)$exon_id
    info <- rep(c("donor", "acceptor"), each = length(donors))
    intronWidths <- abs(start(acceptors) - start(donors)) + 1L
    info <- paste(info, intronWidths)
    values(sites) <- DataFrame(info)
    iit_store(sites, filename)
    filename
}




BamTallyParam <- function(genome, which = GRanges(),
                          desired_read_group = NULL, minimum_mapq = 0L,
                          concordant_only = FALSE, unique_only = FALSE,
                          primary_only = FALSE, ignore_duplicates = FALSE,
                          min_depth = 0L, variant_strand = 0L,
                          ignore_query_Ns = FALSE,
                          indels = FALSE, include_soft_clips = 0L,
                          exon_iit = NULL, IIT_BPPARAM = NULL,
                          xs = FALSE, read_pos = FALSE,
                          min_base_quality = 0L, noncovered = FALSE,
                          nm = FALSE)
{
  if (!is.null(desired_read_group) && !isSingleString(desired_read_group))
    stop("'desired_read_group' must be NULL or a single, non-NA string")
  if (!isSingleNumber(minimum_mapq) || minimum_mapq < 0)
    stop("minimum_mapq must be a single, non-negative, non-NA number")
  if (!isTRUEorFALSE(concordant_only))
    stop("concordant_only must be TRUE or FALSE")
  if (!isTRUEorFALSE(unique_only))
    stop("unique_only must be TRUE or FALSE")
  if (!isTRUEorFALSE(primary_only))
    stop("primary_only must be TRUE or FALSE")
  if (!isTRUEorFALSE(ignore_duplicates))
    stop("ignore_duplicates must be TRUE or FALSE")
  if (!isSingleNumber(min_depth) || min_depth < 0)
    stop("min_depth must be a single, non-negative, non-NA number")
  if (!variant_strand %in% c(0, 1, 2))
    stop("variant_strand must be one of 0, 1, or 2")
  if (!isTRUEorFALSE(ignore_query_Ns))
    stop("ignore_query_Ns must be TRUE or FALSE")
  if (!isTRUEorFALSE(indels))
    stop("indels must be TRUE or FALSE")
  if (include_soft_clips < 0)
    stop("include_soft_clips must be non-negative")
  if (!isTRUEorFALSE(xs))
    stop("xs must be TRUE or FALSE")
  if (!isTRUEorFALSE(read_pos))
    stop("read_pos must be TRUE or FALSE")
  if (!isSingleNumber(min_base_quality) || min_base_quality < 0)
    stop("min_base_quality must be a single, non-negative, non-NA number")
  if (!isTRUEorFALSE(noncovered))
    stop("noncovered must be TRUE or FALSE")
  if (!isTRUEorFALSE(nm))
    stop("nm must be TRUE or FALSE")
  args <- names(formals(sys.function()))
  params <- mget(args, environment())
  params$genome <- as(genome, "GmapGenome")
  params$which <- normArgWhich(which, params$genome)
  params$exon_iit = normArgCdsIIT(params$exon_iit, BPPARAM = IIT_BPPARAM)
  integer_params <- c("minimum_mapq", "min_depth", "variant_strand",
                      "include_soft_clips")
  params[integer_params] <- lapply(params[integer_params], as.integer)
  params = params[names(params) != "IIT_BPPARAM"]
  do.call(new, c("BamTallyParam", params))  
}



### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

setAs("BamTallyParam", "list", function(from) {
  sapply(slotNames(from), slot, object = from, simplify = FALSE)
})

setMethod("as.list", "BamTallyParam", function(x) as(x, "list"))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

showSlot <- function(name, value, ...) {
  S4Vectors:::labeledLine(name, showAsCell(value), ...)
}

showSlots <- function(object, exclude = character(), ...) {
  snames <- setdiff(slotNames(object), exclude)
  slots <- sapply(snames, slot, object = object, simplify = FALSE)
  mapply(showSlot, names(slots), slots, MoreArgs = list(...))
}

setMethod("show", "BamTallyParam", function(object) {
  cat("A", class(object), "object\n", sep = " ")
  cat(showSlots(object, count = FALSE), sep = "")
})
