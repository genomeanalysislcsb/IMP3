#!/usr/bin/perl

use strict;
# Enforces some good programming practices, such as requiring variables to be declared before use.

my $input=$ARGV[0];
my $contigs=$ARGV[1] if $ARGV[1];
# Retrieves command-line arguments. $input is the first argument (expected to be a file), and $contigs is the second, optional argument (also expected to be a file).

my %ids=();
my @ids=();
# Initializes a hash (%ids) to store ID keys with their associated values and an array (@ids) to keep a list of unique IDs.

if ($contigs){
    # If a second file was provided, process it as contig file.
    open(CONTIGS,$contigs) or die $!;
    # Opens the contigs file. If the file can't be opened, the script dies and prints the error message.
    while (my $line=<CONTIGS>){
        # Reads the contigs file line by line.
        chomp($line);
        # Removes the newline character from the end of the line.
        next if $line !~/\S/;
        # Skips the line if it contains only whitespace.
        if ($line=~/^>(.+)$/){
            # If the line starts with a ">", it is considered an ID line, and the ID is captured.
            my $id=$1;
            push @ids, $id;
            # The ID is then added to the @ids array.
        } 
    }
    close(CONTIGS);
    # Closes the contigs file.
}

open(FILE,$input) or die $!;
# Opens the first input file. If it can't be opened, the script dies and prints the error message.
while(my $line=<FILE>){
    # Reads the input file line by line.
	chomp($line);
    # Removes the newline character from the end of the line.
    next if $line !~/\S/;
    # Skips the line if it contains only whitespace.
	my @F=split(/\t/,$line);
    # Splits the line by tabs into an array @F.
	if(exists($ids{$F[0]})){
	    # If the first element of @F (an ID) already exists in %ids, its value is updated.
		$ids{$F[0]}+=$F[1]*$F[4];
	}else{
	    # If the ID doesn't exist in %ids, it is added with an initial value.
		$ids{$F[0]}=$F[1]*$F[4];
        push @ids,$F[0] if !$contigs;
        # If contigs file was not provided, the ID is also added to the @ids array.
	}
}
close(FILE);
# Closes the input file.

foreach my $id (@ids){
    # Iterates over each unique ID.
    if(exists($ids{$id})){
        # If the ID has an associated value in %ids, it is printed along with the value.
        printf("%s\t%.5f\n",$id,$ids{$id});
    }else{
        # If the ID doesn't have an associated value in %ids, it is printed with a value of 0.
        printf("%s\t%.5f\n",$id,0);    
    }
}
