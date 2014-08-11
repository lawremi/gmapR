### =========================================================================
### iit_store command
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### High-level wrapper
###

setGeneric("iit_store", function(x, dest, BPPARAM=MulticoreParam(1), ...) standardGeneric("iit_store"))

gmapRange <- function(x) {
  pos <- strand(x) == "+"
  range <- paste0(seqnames(x), ":", ifelse(pos, start(x), end(x)))
  if (!all(width(x) == 1L)) {
    range <- paste0(range, "..", ifelse(pos, end(x), start(x)))
  }
  range
}



.make_iit_interval = function(grange, name) {
    pos = runValue(strand(grange))[1] == "+"
    strts = start(grange)
    ends = end(grange)
    o = order(strts, decreasing = !pos)
    
    strtpos = if(pos) min(strts) else max(ends)
    endpos = if(!pos) min(strts) else max(ends)
    hdr = paste0(">", name, " ", runValue(seqnames(grange))[1], ":", strtpos, "..", endpos )

    datline = paste(name, name, "NA")
    rngst = if(pos) strts[o] else ends[o]
    rngend = if(!pos) strts[o] else ends[o]
    rnglines = paste0(rngst , " ", rngend)
    c(hdr, datline, rnglines)
}
    
    

setMethod("iit_store", c("GenomicRangesList"),
          function(x, dest =  tempfile(pattern="iit", fileext=".iit"), BPPARAM= MulticoreParam(1)) {
              nms = gsub("( |:|\\.)", "_", names(x))
              lines = unlist(bpmapply(.make_iit_interval, x, nms, BPPARAM = BPPARAM), use.names = FALSE)
              p <- .iit_store(sort = "none", output = dest)
              #writeLines(lines, p)
              cat(paste(lines, collapse="\n"), file =p)
              close(p)
              dest
          })


setMethod("iit_store", c("GenomicRanges"),
          function(x, dest =  tempfile(pattern="iit", fileext=".iit"),
                   info = colnames(values(x))[1]) {
            lines <- paste0(">", names(x), " ", gmapRange(x), " ",
                            values(x)[[info]])
            p <- .iit_store(sort = "none", output = dest)
            writeLines(lines, p)
            close(p)
            dest
          })


setMethod("iit_store", c("character"),
          function(x, dest =  tempfile(pattern="iit", fileext=".iit"),
                   gff = file_ext(x) == "gff",
                   label = if (gff) "ID" else NULL)
          {
            .iit_store(gff = gff, label = label, sort = "none",
                       output = dest, inputfile = x)
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level wrapper
###

.iit_store <- function(fields = FALSE, gff = FALSE,
                       label = if (gff) "ID" else NULL,
                       sort = c("chrom", "none", "alpha", "numeric-alpha"),
                       output, .inputfile = NULL)
{
### TODO: assertions
  sort <- match.arg(sort)
  cl <- commandLine("iit_store")
  if (is.null(.inputfile))
    pipe(cl, open = "w")
  else .system(cl)
}
