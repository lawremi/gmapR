## High-level interface to gsnap 

## We probably want a parameters object. The question is whether we
## really want to maintain accessors, etc for so many parameters.
## Choices:
## - Just have people use a list if they want to group parameters
## - Have an S4 class with every parameter; with accessors
## - Have an S4 class with slots for common parameters, and an extra list

## If we chose the latter, how many arguments would we need?
## --part, so people can use multiple nodes
## --batch, so people can take advantage of their memory
## --max_mismatches, default is reasonable, but people will want to change this
## --use_snps, for SNP tolerance
## --mode, for methylation and RNA editing support
## --nthreads, so people can use their cores
## --novelsplicing, for RNA-seq
## --splicing, known splice sites, for RNA-seq
## --npaths, to limit multi-mappers
## --quiet-if-excessive, to limit multi-mappers
## --nofails, to ignore non-mappers
## --format, to get BAM
## --split_output, to organize the output
## --read_group_id, read_group_name for identifying the reads

## So that is 14 parameters. We definitely need a class to hold
## those. Accessors for 'format' and 'mode' will collide with
## fundamental methods in R. What to do about that? Probably by
## prefixing gsnap on to all of them.

## parallelize_gsnap would also need to take this object,
## overriding the 'part' parameter.

.valid_GsnapParam_part <- function(x) {
  ## validate the part string (i/n)
}
.valid_GsnapParam_batch <- function(x) {
  ## one of 0, 1, 2, 3, 4
}
.valid_GsnapParam_max_mismatches <- function(x) {
  ## non-negative number or NULL
}
.valid_GsnapParam_use_snps <- function(x) {
  ## a file that exists (not sure)
}
.valid_GsnapParam_mode <- function(x) {
  ## one of standard, cmet, atoi
}
.valid_GsnapParam_nthreads <- function(x) {
  ## positive number
}
.valid_GsnapParam_novelsplicing <- function(x) {
  ## TRUE or FALSE
}
.valid_GsnapParam_splicing <- function(x) {
  ## a file that exists (not sure) or NULL
}
.valid_GsnapParam_npaths <- function(x) {
  ## positive number
}
.valid_GsnapParam_quiet_if_excessive <- function(x) {
  ## TRUE or FALSE
}
.valid_GsnapParam_nofails <- function(x) {
  ## TRUE or FALSE
}
.valid_GsnapParam_format <- function(x) {
  ## sam or NULL (this is a low-level check!)
}
.valid_GsnapParam_split_output <- function(x) {
  ## TRUE or FALSE
}
.valid_GsnapParam_read_group_id <- function(x) {
  ## single string (any character constraints?) or NULL
}
.valid_GsnapParam_read_group_name  <- function(x) {
  ## single string or NULL
}

.valid_gmap_parameter <- function(name, params) {
  validator <- get(paste(".valid_GsnapParam_", name, sep = ""))
  if (!is.null(validator))
    validator(params)
  else NULL
}

.valid_GsnapParam <- function(x) {
  x <- as.list(x) # converts to low-level parameter list
  do.call(c, lapply(names(x), .valid_gmap_parameter, x))
}

setClassUnion("IntegerORNULL", c("integer", "NULL"))

setClass("GsnapParam",
         representation(part = "characterORNULL", # used by parallelized_gsnap
                        batch = "character", # weird "0", "1", ... 
                        max_mismatches = "integerORNULL",
                        use_snps = "characterORNULL",
                        mode = "character",
                        nthreads = "integer",
                        novelsplicing = "logical",
                        splicing = "characterORNULL",
                        npaths = "integer",
                        quiet_if_excessive = "logical",
                        nofails = "logical", 
                        format = "character", # bam, sam, gsnap
                        split_output = "logical",
                        read_group_id = "characterORNULL",
                        read_group_name = "characterORNULL"),
         prototype = list(batch = "2", mode = "standard", npaths = 100L,
           quiet_if_excessive = FALSE, nofails = FALSE, format = "bam",
           split_output = FALSE),
         validity = .valid_GsnapParam)

setGeneric("gsnap", function(input_a, input_b = NULL, params, ...)
           standardGeneric("gsnap"))

setMethod("gsnap", c("character", "character", "GsnapParams"),
          function(input_a, input_b, params, ...) {
### TODO: vectorize over input_a and input_b, with recycling
            params <- as.list(initialize(params, ...))
            do.call(.gsnap, c(input_a = input_a, input_b = input_b, params))
### We will end up with a character vector of BAM files
### (or CharacterList when splitting output)
          })

## Low-level interface to gnsap, with all the params
.gsnap <- function(db, part = NULL, input_buffer_size = 1000L,
                   barcode_length = 0L, pc_linefeeds = FALSE,
                   orientation = c("FR", "RF", "FF"), gunzip = FALSE,
                   batch = c("2", "0", "1", "3", "4"), max_mismatches = NULL,
                   query_unk_mismatch = FALSE, genome_unk_mismatch = TRUE,
                   terminal_penalty = 3L, indel_penalty = 2L,
                   indel_endlength = 4L, max_middle_insertions = 9L,
                   max_middle_deletions = 30L, max_end_insertions = 3L,
                   max_end_deletions = 6L, suboptimal_levels = 0L,
                   masking = c("2", "1", "3", "4"), adapter_string = NULL,
                   trim_mismatch_score = -3L, snpsdir = NULL, use_snps = NULL,
                   cmetdir = NULL, atoidir = NULL,
                   mode = c("standard", "cmet", "atoi"),
                   tallydir = NULL, use_tally = NULL, nthreads = NULL,
                   novelsplicing = FALSE, splicing = NULL,
                   novel_doublesplices = FALSE, localsplicedist = 200000L,
                   local_splice_penalty = 0L, distant_splice_penalty = 3L,
                   distant_splice_endlength = 16L,
                   shortend_splice_endlength = 2L,
                   distant_splice_identity = 0.95, pairmax_dna = 1000L,
                   pairmax_rna = 200000L, pairexpect = 200L,
                   quality_protocol = c("sanger", "illumina"),
                   quality_zero_score = 33L, quality_print_shift = 0L,
                   mapq_unique_score = NULL, npaths = 100L,
                   quiet_if_excessive = FALSE, ordered = FALSE,
                   show_rediff = FALSE, clip_overlap = FALSE,
                   print_snps = FALSE, failsonly = FALSE,
                   nofails = FALSE, fails_as_input = FALSE,
                   format = NULL, split_output = NULL, no_sam_headers = FALSE,
                   sam_headers_batch = NULL, read_group_id = NULL,
                   read_group_name = NULL, input_a = NULL, input_b = NULL)
{
  formals <- formals(sys.function())
  problems <-
    do.call(c, lapply(names(formals), .valid_gmap_parameter, formals))
  if (!is.null(problems))
    stop("validation failed:\n  ", paste(problems, collapse = "\n  "))
  
### TODO: if input_a is NULL, return a pipe()
  .system(commandLine("gsnap"))
### TODO: return the bam file path as character vector
### This means interpreting split_output correctly
}
