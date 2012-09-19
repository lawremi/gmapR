test_gsnap_syscall_unique_only_TRUE <- function() {
  gmapGenome <- TP53Genome()
  fastqs <- LungCancerLines::LungCancerFastqFiles()
  gsnapParam <- GsnapParam(genome = gmapGenome,
                           unique_only = TRUE,
                           gunzip=TRUE)
  systemCall <- gmapR:::asSystemCall(gsnap(input_a=fastqs["H1993.first"],
                                           input_b=fastqs["H1993.last"],
                                           params=gsnapParam))
  ##"unique_only=TRUE" means npath=1, no_fails is true, and
  ##split_output is false
  pieces <- unlist(strsplit(systemCall, " "))
  checkTrue(any(pieces == "--npaths=1"))
  checkTrue(any(pieces == "--nofails"))
  checkTrue(!any(pieces == "split-output"))
}

test_gsnap_syscall_unique_only_FALSE <- function() {
  gmapGenome <- TP53Genome()
  fastqs <- LungCancerLines::LungCancerFastqFiles()
  gsnapParam <- GsnapParam(genome = gmapGenome,
                           unique_only = FALSE,
                           gunzip=TRUE)
  systemCall <- gmapR:::asSystemCall(gsnap(input_a=fastqs["H1993.first"],
                                           input_b=fastqs["H1993.last"],
                                           params=gsnapParam))
  ##"unique_only=FALSE" means npath=100L, no_fails is false, and
  ##split_output is true
  pieces <- unlist(strsplit(systemCall, " "))
  ##since npaths=100 is a default, "npaths" should not be in the
  ##command line
  checkTrue(!any(grepl("--npaths=", pieces)))
  checkTrue(!any(pieces == "--nofails"))
  checkTrue(sum(grepl("^--split-output=", pieces)) == 1)
}

test_gsnap_syscall_max_mismatches <- function() {
  gmapGenome <- TP53Genome()
  fastqs <- LungCancerLines::LungCancerFastqFiles()
  gsnapParam <- GsnapParam(genome = gmapGenome,
                           unique_only = FALSE,
                           gunzip=TRUE,
                           max_mismatches=8)
  systemCall <- gmapR:::asSystemCall(gsnap(input_a=fastqs["H1993.first"],
                                           input_b=fastqs["H1993.last"],
                                           params=gsnapParam))
  pieces <- unlist(strsplit(systemCall, " "))
  checkTrue(any(pieces == "--max-mismatches=8"))
}

test_gsnap_syscall_suboptimal_levels <- function() {
  gmapGenome <- TP53Genome()
  fastqs <- LungCancerLines::LungCancerFastqFiles()
  gsnapParam <- GsnapParam(genome = gmapGenome,
                           unique_only = FALSE,
                           gunzip=TRUE,
                           suboptimal_levels=3)
  systemCall <- gmapR:::asSystemCall(gsnap(input_a=fastqs["H1993.first"],
                                           input_b=fastqs["H1993.last"],
                                           params=gsnapParam))
  pieces <- unlist(strsplit(systemCall, " "))
  checkTrue(any(pieces == "--suboptimal-levels=3"))
}

test_gsnap_syscall_suboptimal_levels_default <- function() {
  gmapGenome <- TP53Genome()
  fastqs <- LungCancerLines::LungCancerFastqFiles()
  gsnapParam <- GsnapParam(genome = gmapGenome,
                           unique_only = FALSE,
                           gunzip=TRUE,
                           suboptimal_levels=0)
  systemCall <- gmapR:::asSystemCall(gsnap(input_a=fastqs["H1993.first"],
                                           input_b=fastqs["H1993.last"],
                                           params=gsnapParam))
  pieces <- unlist(strsplit(systemCall, " "))
  checkTrue(!any(grepl("^--suboptimal-levels=", pieces)))
}

test_gsnap_mode_incorrect <- function() {
  gmapGenome <- TP53Genome()
  fastqs <- LungCancerLines::LungCancerFastqFiles()
  gsnapParam <- GsnapParam(genome = gmapGenome,
                           gunzip=TRUE,
                           mode="cmet")
  checkException(gsnap(input_a=fastqs["H1993.first"],
                       input_b=fastqs["H1993.last"],
                       params=gsnapParam))
}

test_gsnap_syscall_mode <- function() {
  gmapGenome <- TP53Genome()
  fastqs <- LungCancerLines::LungCancerFastqFiles()
  gsnapParam <- GsnapParam(genome = gmapGenome,
                           gunzip=TRUE,
                           mode="cmet-stranded")
  systemCall <- gmapR:::asSystemCall(gsnap(input_a=fastqs["H1993.first"],
                                           input_b=fastqs["H1993.last"],
                                           params=gsnapParam))
  pieces <- unlist(strsplit(systemCall, " "))
  checkTrue(any(pieces == "--mode=cmet-stranded"))
}

test_gsnap_syscall_novelsplicing_FALSE <- function() {
  gmapGenome <- TP53Genome()
  fastqs <- LungCancerLines::LungCancerFastqFiles()
  gsnapParam <- GsnapParam(genome = gmapGenome,
                           novelsplicing = FALSE,
                           gunzip=TRUE)  
  systemCall <- gmapR:::asSystemCall(gsnap(input_a=fastqs["H1993.first"],
                                           input_b=fastqs["H1993.last"],
                                           params=gsnapParam))
  pieces <- unlist(strsplit(systemCall, " "))
  checkTrue(sum(grepl("^--novelsplicing=", pieces)) == 0)
}


test_gsnap_syscall_novelsplicing_TRUE <- function() {
  gmapGenome <- TP53Genome()
  fastqs <- LungCancerLines::LungCancerFastqFiles()
  gsnapParam <- GsnapParam(genome = gmapGenome,
                           unique_only = FALSE,
                           max_mismatches = NULL,
                           suboptimal_levels = 0, mode = "standard",
                           npaths = 10,
                           novelsplicing = TRUE,
                           splicing = NULL, 
                           nthreads = 1,
                           batch = 2L,
                           gunzip=TRUE)
  
  systemCall <- gmapR:::asSystemCall(gsnap(input_a=fastqs["H1993.first"],
                                           input_b=fastqs["H1993.last"],
                                           params=gsnapParam))
  pieces <- unlist(strsplit(systemCall, " "))
  checkTrue("--novelsplicing=1" %in% pieces)
}
