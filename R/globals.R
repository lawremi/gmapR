###
### Reading global HTSeqGenieBase settings
###

.PARSE_COLS <- c("gsnap_parser_col.names", "gsnap_parser_colClasses")

.removeCommentsFromFile <- function(infile, outfile, comment.char="#")
{
    if (!isSingleString(infile))
        stop("'infile' must be a single string")
    if (!isSingleString(outfile))
        stop("'outfile' must be a single string")
    if (file.exists(outfile))
        stop("file '", outfile, "' already exists")
    infile <- file(infile, "r")
    #on.exit(close(infile))
    outfile <- file(outfile, "w")
    #on.exit(close(outfile)) # doesn't seem to work
    while (TRUE) {
        text <- readLines(infile, n=1L)
        if (length(text) == 0L)
            break
        if (substr(text, 1L, nchar(comment.char)) != comment.char)
            writeLines(text, outfile)
    }
    close(outfile)
    close(infile)
}

.read.globals <- function(file)
{
    tmp_file <- tempfile()
    .removeCommentsFromFile(file, tmp_file)
    on.exit(file.remove(tmp_file))
    m <- read.dcf(tmp_file)
    getLastNonNA <- function(x)
    {
        ii <- which(!is.na(x))
        unname(x[ii[length(ii)]])
    }
    ans <- lapply(seq_len(ncol(m)),
                  function(j) {
                      val <- getLastNonNA(m[ , j])
                      if (colnames(m)[j] %in% .PARSE_COLS)
                          val <- eval(parse(text=val))
                      val
                  })
    names(ans) <- colnames(m)
    ans
}

##' System-specific locations of files needed for HTSeq pipelines
##'
##' This function returns a named list of locations on disk to data
##' files and executables needed by pipelines
##' @title System-specific pipeline executables and file locations
##' @return a named list
##' @author Cory Barr
##' @export
globals <- function()
{
    default_file <- system.file("gmapR_globals-default.dcf",
                                package="gmapR")
    ans <- .read.globals(default_file)
    user_file <- "~/.gmapR_globals.dcf"
    if (file.exists(user_file)) {
        user_globals <- .read.globals(user_file)
        ans[names(user_globals)] <- user_globals
    }
    gmap_tools <- c("gmap", "gsnap", "gmap_build")
    gmap_tools_list <- as.list(paste(ans$gsnap_bin_dir, gmap_tools, sep="/"))
    names(gmap_tools_list) <- gmap_tools
    append(ans, gmap_tools_list)
}
