#!/usr/bin/env perl

use strict;

my $file=$ARGV[0];
my $name=$ARGV[1];

die "# ERROR new string is missing" if length($name) eq 0;
my $ct=0;

open(FILE,$file) or die "# ERROR: Cannot open $file\n";;
while(my $str=<FILE>){

	chomp($str);

	next if length($str)==0;

	if ($str =~ />(.+$)/){
		$ct++;
		print ">$name"."_"."$ct\n";
	}else{
		print $str,"\n";
	}
}
close(FILE);
