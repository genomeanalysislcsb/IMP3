#!/usr/bin/env perl
## split_by_length.pl
#AHB, Feb 2019
#
use strict;
use warnings;
use Getopt::Long;

my ($midlen,$infile,$shortfile,$longfile);
GetOptions('m=i' => \$midlen,
           'i=s' => \$infile,
           's=s' => \$shortfile,
           'l=s' => \$longfile);

open (IN, $infile);
open (LONG, ">", $longfile); 
open (SHORT, ">", $shortfile);
local $/=">";
while(<IN>) {
	chomp;
        next unless /\w/;
        s/>$//gs;
        my @chunk = split /\n/;
        my $header = shift @chunk;
        my $seqlen = length join "", @chunk;
        if($seqlen >= $midlen){
            print LONG ">$_";
        }else{
	    print SHORT ">$_";
	}	
}
close (IN);
close (LONG);
close (SHORT);



