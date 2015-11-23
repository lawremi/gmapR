### =========================================================================
### gmap command
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### High-level interface
###

setGeneric("gmap", function(input, params, ...) standardGeneric("gmap"))

setMethod("gmap", c("ANY", "GmapParam"),
          function(input, params, ...) {
              tmpfile <- tempfile()
              export(input, tmpfile, format="fasta")
              input <- tmpfile
              callGeneric()
          })

setMethod("gmap", c("FastaFile", "GmapParam"),
          function(input, params, ...) {
              input <- path(input)
              callGeneric()
          })

setMethod("gmap", c("character", "GmapParam"),
          function(input, params,
                   output = file.path(getwd(),
                       file_path_sans_ext(basename(input), TRUE)), ...)
              {
                  if (any(is.na(input)))
                      stop("'input' must not contain NA's")
                  if (length(input) > 1L) {
                      return(GmapOutputList(mapply(gmap, input,
                                                   MoreArgs =
                                                       list(params, output,
                                                            ...))))
                  }
                  
                  output_dir <- dirname(output)
                  if (!file.exists(output_dir))
                      dir.create(output_dir, recursive = TRUE)

                  params <- initialize(params, ...)
                  params_list <- as.list(params)
                  if (gsnap_split_output(params)) {
                      params_list$split_output <- output
                      output_path <- output_dir
                  } else {
                      output_path <- paste0(output, ".",
                                            formatToExt(params$format))
                      params_list$.redirect <- paste(">", output_path)
                  }

                  res <- do.call(.gmap,
                                 c(list(.input = input), params_list))
                  gmap_output<- GmapOutput(path = output_path,
                                           version = gmapVersion(),
                                           param = params)

                  if (params_list$format %in% c("samse", "sampe")) {
                      gmap_output <- asBam(gmap_output)
                  }
                  gmap_output
              })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level interface
###

.gmap <- function(db = NULL, dir = NULL,
                  kmer = NULL, basesize = NULL, sampling = NULL,
                  genomefull = FALSE, gseg = NULL,
                  selfalign = FALSE, pairalign = FALSE,
                  part = NULL, input_buffer_size = 1000L,
                  batch = c("0", "1", "2", "3", "4", "5"),
                  expand_offsets = FALSE, nosplicing = FALSE,
                  min_intronlength = 9L, intronlength = 1000000L,
                  localsplicedist = 2000000L, totallength = 2400000L,
                  chimera_margin = 40L, no_chimeras = FALSE,
                  nthreads = NULL, chrsubsetfile = NULL,
                  direction = c("auto", "sense_force", "antisense_force",
                      "sense_filter", "antisense_filter"),
                  trimendexons = 12L, canonical_mode = c("1", "0", "2"),
                  cross_species = FALSE, allow_close_indels = c("1", "0", "2"),
                  microexon_spliceprob = 0.90, cmetdir = NULL, atoidir = NULL,
                  mode = c("standard", "cmet-stranded", "cmet-nonstranded",
                      "atoi-stranded", "atoi-nonstranded"),
                  prunelevel = c("0", "1", "2", "3"),
                  format = NULL, npaths = 5L, quiet_if_excessive = FALSE,
                  suboptimal_score = NULL, ordered = FALSE, md5 = FALSE,
                  chimera_overlap = FALSE, failsonly = FALSE, nofails = FALSE,
                  fails_as_input = FALSE, snpsdir = NULL, use_snps = NULL,
                  split_output = NULL, append_output = FALSE,
                  output_buffer_size = 1000L, fulllength = FALSE,
                  cdsstart = NULL, truncate = FALSE, tolerant = FALSE,
                  no_sam_headers = FALSE, sam_use_0M = FALSE,
                  force_xs_dir = FALSE, md_lowercase_snp = FALSE,
                  read_group_id = NULL, read_group_name = NULL,
                  read_group_library = NULL, read_group_platform = NULL,
                  quality_protocol = c("sanger", "illumina"),
                  quality_print_shift = 0L, mapdir = NULL, map = NULL,
                  mapexons = FALSE, mapboth = FALSE, version = FALSE,
                  .input = NULL, .redirect = NULL)
{
    formals <- formals(sys.function())

    expand_offsets <- as.integer(expand_offsets)
    batch <- match.arg(batch)
    direction <- match.arg(direction)
    canonical_mode <- match.arg(canonical_mode)
    allow_close_indels <- match.arg(allow_close_indels)
    mode <- match.arg(mode)
    prunelevel <- match.arg(prunelevel)
    quality_protocol <- match.arg(quality_protocol)
    
    if (version) {
        .redirect <- ">/dev/null"
    }

### TODO: if input_a is NULL, or split_output and .redirect are NULL:
###       return a pipe()
    .system_gsnap(commandLine("gmap"))
}

..gmap <- function(args, path = NULL) {
    .system_gsnap(.commandLine("gmap", args, path))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities
###

gmapVersion <- function() {
    output <- .gmap(version = TRUE)
    version_text <- sub("GMAP version (.*?) .*", "\\1", output[1])
    parseGsnapVersion(version_text)
}

