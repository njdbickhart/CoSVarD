#!/usr/bin/perl
# simple file merger script designed to run prior to the start of other scripts

use strict;

my $out = pop(@ARGV);
my @input = @ARGV;
chomp @input;
chomp $out;

open(OUT, "> $out");
foreach my $i (@input){
	chomp $i;
	open(IN, "< $i") || print "Error opening $i!\n";
	while(<IN>){
		print OUT "$_";
	}
	close IN;
}
exit;