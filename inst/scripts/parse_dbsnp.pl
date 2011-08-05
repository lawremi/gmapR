#Hey Cory,
#
#I think we should be including class 'het' and the SNPs from class 'mixed' as well. The het class is easy, parsing SNPs out of the mixed class is going to be a bit more of a trick...
#
#Jeremiah
#
#On Fri, Nov 5, 2010 at 5:36 PM, Cory Barr <barr.cory@gene.com> wrote:
#
#I’m just keeping class=’single’ and locType=’exact’.  Too exclusive?

use warnings;
use strict;

sub revcomp {
    my @nucs = @_;
    @nucs = map(uc, @nucs);

    for (my $i=0; $i<@nucs; $i++) {
	if ($nucs[$i] eq "A") {
	    $nucs[$i] = "T";
	}
	elsif ($nucs[$i] eq "C") {
	    $nucs[$i] = "G";
	}
	elsif ($nucs[$i] eq "G") {
	    $nucs[$i] = "C";
	}
	elsif ($nucs[$i] eq "T") {
	    $nucs[$i] = "A";
	}
    }
    
    return(@nucs);
}

while (my $input=<>) {
    chomp($input);
    my @pieces = split("\t", $input);
    my $chr = $pieces[1];
    my $start = $pieces[2] + 1;
    my $end = $pieces[3];
    my $id = $pieces[4];
    my $strand = $pieces[6];
    my $ref = $pieces[8];
    my $seen = $pieces[9];
    my $type = $pieces[11];

    if ($type eq 'het') { $seen =~ s|\(HETEROZYGOUS\)/?||; }
    elsif ($type ne 'single' && $type ne 'mixed') {
	next;
    }

    if (not $ref =~ m/^[ACGT]$/) { next; } #cannot be a SNP unless matches this

    my @seen_nucs = split("/", $seen);
    if ($type eq 'mixed') {
	@seen_nucs = grep(/^[ACGT]$/, @seen_nucs);
    }

    if ($strand eq "-") {
	@seen_nucs = revcomp(@seen_nucs);
    }

    @seen_nucs = grep(!/^$ref/, @seen_nucs);
        
#all entries must be in this format: >rs62211261 21:14379270..14379270 CG
    foreach my $sn (@seen_nucs) {
	print 
	    ">" .
	    $id .
	    " " .
	    $chr .
	    ":" .
	    $start .
	    ".." .
	    $end .
	    " " .
	    $ref .
	    $sn .
	    "\n";
    }
}
