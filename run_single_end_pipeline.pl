#!/usr/bin/perl
# This is the subscript that will take input reads from the main pipeline and take it through the altered
# splitread pipeline I have created

use strict;
use Getopt::Std;
#use localtime;

my $usage = "$0 -i \<input interleaved fastq\> -o \<output directory\> -r \<reference\>
\t-e \<edit distance\> -n \<fq num\> -x \<outbase name\>\n";
my $mrsfast = "/ibfs7/asg2/bickhartd/bin/mrsfast";
my $samtools = "/ibfs7/asg2/bickhartd/bin/samtools";


my %opts;
getopt('iorenx', \%opts);

unless(defined($opts{'i'}) && defined($opts{'o'}) && defined($opts{'r'})){
	print "Improper number of arguments\n$usage";
	exit;
}

my $refindex = "$opts{r}.fai";
my @base = split(/\//, $opts{'i'});
my @fqname = split(/\./, $base[-1]);
#my $outbase = $fqname[0];
my $outbase = $opts{'x'};
# file names
my $sam_out = "$opts{o}/$outbase/$fqname[0].$fqname[1].$opts{n}.sam";
my $unmap = "$opts{o}/$outbase/$fqname[0].$fqname[1].$opts{n}.nohit";
my $bam_out = "$opts{o}/$outbase/$fqname[0].$fqname[1].$opts{n}.bam";
my $bam_sorted_out = "$opts{o}/$outbase/$fqname[0].$fqname[1].$opts{n}.single.sorted";

mkdir("$opts{o}") || print $! . "\n";
mkdir("$opts{o}/$outbase") || print $! . "\n";

my $initial_time = time;

my $tstr = return_time_str($initial_time);
print "Running MrsFAST alignment: $tstr\n";
system("$mrsfast --search $opts{r} --seq $opts{i} -e $opts{e} -o $sam_out -u $unmap");

my $tstr = return_time_str($initial_time);
print "Converting output to bam; saving: $tstr\n";
system("$samtools view -bt $refindex $sam_out > $bam_out");
system("$samtools sort $bam_out $bam_sorted_out");
system("rm $sam_out $unmap $bam_out");
print "removing single end fastq $opts{i}\n";
system("rm $opts{i}");
exit;

sub return_time_str {
	#require localtime;
	my ($initial_time) = @_;
	#my @time = localtime(time);
	#$time[5] += 1900;
	my $secs = time - $initial_time;
	#my $tstr = "$time[4]/$time[3]/$time[5] $time[2]:$time[1]:$time[0]. Secs from start: $secs";
	my $tstr = "Secs from start: $secs";
	return $tstr;
}