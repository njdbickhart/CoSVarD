#!/usr/bin/perl
# This script will take the four DOC files and format them into a final bed format.
# doc bed: chr, start, end, gain/loss, copynumber

use strict;
use Getopt::Std;

my $usage = "$0 -i \<auto dups\> -d \<auto dels\> -x \<x dups\> -y \<x dels\> -c \<cn file\> -o \<output\>\n";
my %opts;
getopt('idxyco', \%opts);
unless(defined($opts{'i'}) && defined($opts{'d'}) && defined($opts{'x'}) && defined($opts{'y'})){
	print $usage;
	exit;
}

my %finaldata; #{chr}->[start, end, gain/loss, cn]
my %cnhash; # {chr}->[start, end, value]

# Store cn data for future lookup
open (IN, "< $opts{c}") || die "could not open cnfile: $opts{c}!\n";
while(my $line = <IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	push(@{$cnhash{$segs[0]}}, [$segs[1], $segs[2], $segs[3]]); 
}
close IN;

# Start working on input files
open_file_produce_output($opts{'i'}, "gain", \%cnhash, \%finaldata);
open_file_produce_output($opts{'d'}, "loss", \%cnhash, \%finaldata);
open_file_produce_output($opts{'x'}, "gain", \%cnhash, \%finaldata);
open_file_produce_output($opts{'y'}, "loss", \%cnhash, \%finaldata);

open (OUT, "> $opts{o}"); 

foreach my $chrs (sort {my ($achrs) = $a =~ m/chr(.+)/; my ($bchrs) = $b =~ m/chr(.+)/; if($achrs =~ /X/){$achrs = 500;} if($bchrs =~ /X/){$bchrs = 500;} $achrs <=> $bchrs} keys(%finaldata)){
	foreach my $rows (@{$finaldata{$chrs}}){
		print OUT "$chrs\t" . join("\t", @{$rows}) . "\n";
	}
}
close OUT;

exit;

sub open_file_produce_output {
	my($file, $type, $cndata, $finalhash) = @_;
	open(IN, "< $file") || print "could not open $file!\n";
	while (my $line = <IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		my @cnvals = ();
		foreach my $rows (@{$cndata->{$segs[0]}}){
			if($rows->[1] > $segs[1] && $rows->[0] < $segs[2]){
				push(@cnvals, $rows->[2]);
			}
		}
		if(scalar(@cnvals) == 0){next;} #Filter out anomalous entries that are due to gap filtration
		my $avgcn = average(\@cnvals);
		push(@{$finalhash->{$segs[0]}}, [$segs[1], $segs[2], $type, $avgcn]);
	}
	close IN;
}

sub average {
	my ($aref) = @_;
	my ($sum, $count);
	if(scalar(@{$aref}) == 0){ return 0;}
	foreach my $a (@{$aref}){
		$sum += $a;
		$count++;
	}
	if ($count == 0){ return 0;}
	return $sum / $count;
}