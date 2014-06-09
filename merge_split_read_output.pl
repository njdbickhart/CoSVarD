#!/usr/bin/perl
# This script takes the split read output and reformats it into bed format output
# split bed: chr, start, end, del/inv, balanced support, unbalanced support, anchor edit, split edit

use strict;
use Getopt::Std;

my $input = $ARGV[0];
my $output = $ARGV[1];
chomp $input;
chomp $output;

open (IN, "< $input") || die "Could not open $input\n";
open (OUT, "> $output") || die "Could not create $output!\n";
while (my $line = <IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	my @event = split(//, $segs[4]);
	my $len = $segs[2] - $segs[1];
	my $type;
	if($event[0] eq "D" && $len > 10000){	# too large to be a deletion; likely a MEI
		$type = "ins";
	}elsif ($event[0] eq "D"){
		$type = "del";
	}elsif($event[0] eq "I"){
		$type = "inv";
	}else{
		print "unknown type: $segs[4]\n";
		$type = "unk";
	}
	print OUT "$segs[0]\t$segs[1]\t$segs[2]\t$type\t$segs[5]\t$segs[6]\t$segs[7]\t$segs[8]\n";
}
close IN;
close OUT;
exit;

