\name{buildGsnapSNPIIT}
\alias{buildGsnapSNPIIT}
\title{Build SNP-tolerant gsnap alignment}
\usage{buildGsnapSNPIIT(snp_file, gsnap_bin, gsnap_data_dir, genome,
    snps_name)}
\description{This function helps create the files needed for SNP-tolerant gsnap
alignment.}
\value{0}
\arguments{\item{snp_file}{The fasta-like file containing the SNPs}
\item{gsnap_bin}{location of the gsnap executable. Defaults to
location at gsnap:::globals()$gsnap}
\item{gsnap_data_dir}{location of the gmap data directory where indexes are stored. Defaults
to location at gsnap:::globals()$gsnap_save_dir}
\item{genome}{name of the genome, which corresponds to
the subdirectory of save_dir. The IIT files go here. Typical
examples are "hg19" and "mm9"}
\item{snps_name}{This argument specifies the name of the
argument to call via gsnap's -v option for snp-tolerant alignment. For example, if snp_name
is db131, then gsnap will align via gsnap will incorporate your injected snps through the command 'gsnap -d genome -v db131'}}
