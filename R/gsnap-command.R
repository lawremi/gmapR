### =========================================================================
### gsnap command
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### High-level interface
###

setGeneric("gsnap", function(input_a, input_b = NULL, params, ...)
           standardGeneric("gsnap"))

setMethod("gsnap", c("character", "characterORNULL", "GsnapParam"),
          function(input_a, input_b, params,
                   output = file_path_sans_ext(input_a, TRUE),
                   consolidate = TRUE, ...)
          {
            output_dir <- dirname(output)
            if (!file.exists(output_dir))
              dir.create(output_dir, recursive = TRUE)

            params <- initialize(params, ...)
            params_list <- as.list(params)
            if (gsnap_split_output(params)) {
              params_list$split_output <- basename(output)
              output_path <- output_dir
            } else {
              output_path <- paste0(output, ".sam")
              params_list$.redirect <- paste(">", output_path)
            }

            if (is.null(params_list$gunzip)) {
              if (file_ext(input_a) == "gz" && file_ext(input_b) == "gz") {
                params_list$gunzip <- TRUE
              }
            }
            
            res <- do.call(.gsnap,
                         c(list(.input_a = input_a, .input_b = input_b,
                                format = "sam"),
                           params_list))
            gsnap_output <- GsnapOutput(path = output_path,
                                        version = gsnapVersion(),
                                        param = params)
            gsnap_output <- asBam(gsnap_output)
            if (consolidate)
              consolidate(gsnap_output)
            return(gsnap_output)
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level interface
###

.gsnap <- function(db = NULL, dir = NULL, part = NULL,
                   input_buffer_size = 1000L,
                   barcode_length = 0L,
                   orientation = c("FR", "RF", "FF"),
                   fastq_id_start = NULL, fastq_id_end = NULL,
                   filter_chastity = c("off", "either", "both"),
                   gunzip = FALSE,
                   batch = c("2", "0", "1", "3", "4"), max_mismatches = NULL,
                   query_unk_mismatch = 0L, genome_unk_mismatch = 1L,
                   terminal_threshold = 2L,
                   indel_penalty = 2L,
                   indel_endlength = 4L, max_middle_insertions = 9L,
                   max_middle_deletions = 30L, max_end_insertions = 3L,
                   max_end_deletions = 6L, suboptimal_levels = 0L,
                   adapter_string = NULL,
                   trim_mismatch_score = -3L, trim_indel_score = -4L,
                   snpsdir = NULL, use_snps = NULL,
                   cmetdir = NULL, atoidir = NULL,
                   mode = c("standard", "cmet-stranded", "cmet-nonstranded",
                     "atoi-stranded", "atoi-nonstranded"),
                   tallydir = NULL, use_tally = NULL,
                   use_runlength = NULL, nthreads = 1L,
                   gmap_mode = "pairsearch,terminal,improve",
                   trigger_score_for_gmap = 5L, max_gmap_pairsearch = 3L,
                   max_gmap_improvement = 3L, microexon_spliceprob = 0.90,
                   novelsplicing = 0L, splicingdir = NULL,
                   use_splicing = NULL, ambig_splice_noclip = FALSE,
                   localsplicedist = 200000L,
                   local_splice_penalty = 0L, distant_splice_penalty = 3L,
                   distant_splice_endlength = 16L,
                   shortend_splice_endlength = 2L,
                   distant_splice_identity = 0.95,
                   antistranded_penalty = 0L, merge_distant_samechr = FALSE,
                   pairmax_dna = 1000L,
                   pairmax_rna = 200000L, pairexpect = 200L, pairdev = 25L,
                   quality_protocol = c("sanger", "illumina"),
                   quality_zero_score = 33L, quality_print_shift = 0L,
                   mapq_unique_score = NULL, npaths = 100L,
                   quiet_if_excessive = FALSE, ordered = FALSE,
                   show_rediff = FALSE, clip_overlap = FALSE,
                   print_snps = FALSE, failsonly = FALSE,
                   nofails = FALSE, fails_as_input = FALSE,
                   format = NULL, split_output = NULL,
                   output_buffer_size = 1000L,
                   no_sam_headers = FALSE,
                   sam_headers_batch = NULL, read_group_id = NULL,
                   read_group_name = NULL, read_group_library = NULL,
                   read_group_platform = NULL, version = FALSE,
                   .input_a = NULL, .input_b = NULL,
                   .redirect = NULL ## e.g., > gsnap.sam
                   )
{
  formals <- formals(sys.function())
  problems <- do.call(c, lapply(names(formals), .valid_gmap_parameter, formals))
  if (!is.null(problems))
    stop("validation failed:\n  ", paste(problems, collapse = "\n  "))  

  orientation <- match.arg(orientation)
  batch <- match.arg(batch)
  mode <- match.arg(mode)
  quality_protocol <- match.arg(quality_protocol)
  filter_chastity <- match.arg(filter_chastity)
  
### TODO: if input_a is NULL, or split_output and .redirect are NULL:
###       return a pipe()
  .system_gsnap(commandLine("gsnap"))
}

..gsnap <- function(args, path = NULL) {
  .system_gsnap(.commandLine("gsnap", args, path))
}

.system_gsnap <- function(command) {
  stderr <- tempfile("gsnap-stderr", fileext = "log")
  command <- paste(command, "2>", stderr)
  .system(command)
  output <- readLines(stderr)
  if (!gsnapSucceeded(command, output))
    stop("Execution of gsnap failed, command-line: '", command,
         "'; last output line: '", tail(output, 1), "'")
  else output
}

gsnapSucceeded <- function(command, output) {
  if (!grepl("--version", command))
    any(grepl("^Processed \\d+ queries", output))
  else TRUE
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity Checking
###

.valid_GsnapParam_part <- function(x) {
  ## validate the part string (i/n)
}
.valid_GsnapParam_batch <- function(x) {
  ## one of 0, 1, 2, 3, 4
}
.valid_GsnapParam_max_mismatches <- function(x) {
  ## non-negative number or NULL
}
.valid_GsnapParam_suboptimal_levels <- function(x) {
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
  ## single string or NULL
}
.valid_GsnapParam_read_group_id <- function(x) {
  ## single string (any character constraints?) or NULL
}
.valid_GsnapParam_read_group_name  <- function(x) {
  ## single string or NULL
}

.valid_gmap_parameter <- function(name, params) {
  res <- TRUE ##if we're not validating the param, assume it is good

  validatorName <- paste(".valid_GsnapParam_", name, sep = "")
  if (exists(validatorName)) {
    validator <- get(validatorName)
    res <- validator(params)
  } else {
    res <- NULL
  }

  return(res)
}

.valid_GsnapParam <- function(object) {
  x <- as.list(object) # converts to low-level parameter list
  do.call(c, lapply(names(x), .valid_gmap_parameter, x))
}

setValidity("GsnapParam", .valid_GsnapParam)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities
###

gsnapVersion <- function() {
  output <- .gsnap(version = TRUE)
  version_text <- sub("GSNAP version (.*?) .*", "\\1", output[1])
  parseGsnapVersion(version_text)
}

parseGsnapVersion <- function(x) {
  as.POSIXlt(x, format = "%Y-%m-%d")
}
