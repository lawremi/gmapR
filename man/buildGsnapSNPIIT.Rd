\name{buildGsnapSNPIIT}
\alias{buildGsnapSNPIIT}
\title{Build SNP-tolerant gsnap alignment}
\usage{buildGsnapSNPIIT(index_file, gsnap_bin, save_dir, gsnap_genome_dir,
    snp_file)}
\description{This function helps create the files needed for SNP-tolerant gsnap
alignment.}
\value{0}
\arguments{\item{index_file}{The fasta-like file containing the SNPs}
\item{gsnap_bin}{location of the gsnap executable. Defaults to
location at gsnap:::globals()$gsnap}
\item{save_dir}{location of the default gmap directory. Defaults
to location at gsnap:::globals()$gsnap_save_dir}
\item{gsnap_genome_dir}{name of the genome, which corresponds to
the subdirectory of save_dir. The IIT files go here. Typical
examples are "hg19" and "mm9"}
\item{snp_file}{This argument will specify the name of the
resulting final IIT of SNPs and consequently the name of the
database specified via gsnap's -v option. For example, if snp_file
is db131, then gsnap will align via gsnap -d genome -v db131
input.fastq}}
