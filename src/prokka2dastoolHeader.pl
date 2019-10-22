#!/usr/bin/perl

use strict;

my $prokka=$ARGV[0];
my $contigs=$ARGV[1];

my %ids=();
my %contigs=();

open(CONTIGS,$contigs) or die $!;
while( my $str=<CONTIGS>){
	chomp($str);
	next if $str!~/\S/;
	my ($contig,$id,$num)=split(/\t/,$str);
	$ids{$id."_".$num}=$contig;	
	$contigs{$contig}{$num}=0;;
}
close(CONTIGS);

foreach my $contig (keys(%contigs)){
	my $ct=1;
	foreach my $num ( sort {$a <=> $b} keys %{$contigs{$contig}}){
		$contigs{$contig}{$num}=$ct;
		$ct++;
	}			
}

open(PROKKA,$prokka) or die $!;
while(my $str=<PROKKA>){
	chomp($str);
	next if $str!~/\S/;
	if($str=~/^>(.+)_(\d+)/){
		my $id=$1;
		my $num=$2;
		if (exists($ids{$id."_".$num})){
			my $contig=$ids{$id."_".$num};
			$str=">".$contig."_".$contigs{$contig}{$num};
		}else{
			print STDERR "# ERROR: $str: $id"."_"."$num does not exist\n";
			exit(0);
	        }
	}elsif($str=~/^>/){
		print STDERR "# ERROR: wrong header: $str\n";
	}
	print $str."\n";
}
close(PROKKA);
