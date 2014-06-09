#!/usr/bin/perl
# This script identifies all of the mrsfast sam files in a directory and turns them into bam files to save space

use strict;
use Getopt::Std;
use FileHandle;

my $usage = "$0 -i <directory> -r <reference fai file> -b <bin directory> -c <create cleanup file>\n";
my %opts;
getopt('irbc', \%opts);

opendir(DIR, $opts{'i'}) || die "[Convert Sam to Bam] Could not open dir: $opts{'i'}!\n";
my $cf;
if($opts{'c'}){
	$cf = FileHandle->new("> $opts{i}/cleanup.sh");
}
while (my $file = readdir(DIR)){
	my $path = $opts{'i'} . "/$file";
	unless(-f $path && $file =~ m/.+mrsfast.sam/){
		next;
	}
	
	my @fsegs = split(/\./, $file);
	pop(@fsegs);
	my $outbam = $opts{'i'} . "/" . join(".", @fsegs) . ".bam";
	
	my $fai = $opts{'r'};
	my $bindir = $opts{'b'};
	system("$bindir/samtools view -bS -t $fai -o $outbam $path");
	
	if($opts{'c'}){
		print $cf "rm $path\n";
	}else{
		if(-s $outbam){
			system("rm $path");
		}
	}
}
closedir(DIR);