#!/usr/bin/perl

use warnings;
use strict;
use Array::Transpose;

my @dat;
&read_datN($ARGV[0], \@dat);
@dat=transpose(\@dat);

for( my $i = 1; $i < scalar(@dat); $i++ ) {
    &output_tdat(sprintf("%s_%03d.dat", $ARGV[1], $i-1), \@{$dat[0]}, \@{$dat[$i]});
}

#print STDERR sprintf("\nstatus: wrote all data items to %d files\n", scalar(@dat)-1);

exit 0;

sub read_datN($$) 
{
    my $fname = $_[0]; # filename for reading data
    my $a_ref = $_[1]; # reference for AoA to store acts

    open IN_F, "< $fname" or die "can't open data file \"$fname\"\n";
    #print STDERR "status: reading data from $fname\n"; 

    my $ntot = 0;
    while( <IN_F> ) {
        chomp;
        my @data_read = split /\s+/, $_;
        push @$a_ref, [ @data_read ];
        $ntot += scalar(@data_read);
    }
    close(IN_F);

    #print STDERR sprintf("status: read $ntot values in %d rows from $fname\n", scalar(@$a_ref)); 
}

sub output_tdat ($$$) 
{
    my $fname = $_[0]; # filename for output
    my $a_ref = $_[1]; # reference to array (t)
    my $b_ref = $_[2]; # reference to array (v)
    
    open OUTF, "> $fname" or
        die "error: can't open output file \"$fname\"\n";

    for( my $i = 0; $i < scalar(@$a_ref); $i++ ) {
        print OUTF "$a_ref->[$i] $b_ref->[$i]\n";
    }
    
    close OUTF;
    #print STDERR sprintf("status: wrote %d data items to file: %s\r", scalar(@$a_ref), $fname);
}

