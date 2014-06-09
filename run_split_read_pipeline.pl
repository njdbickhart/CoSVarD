#!/usr/bin/perl
# This is the subscript that will take input reads from the main pipeline and take it through the altered
# splitread pipeline I have created

use strict;
use Getopt::Std;
#use localtime;

my $usage = "$0 -1 \<input sam1\> -2 \<input sam2\> -a \<fq1\> -z \<fq2\> -o \<output directory\> -r \<reference\>
\t-e \<edit distance\> -l \<low bound\> -m \<max bound\> -u \<upper\> -x \<read len\> -n \<fq num\>
\t-c \<outbase name\> -b \<bin dir\> -j \<java dir\>\n";



my %opts;
getopt('12orelmuxncbjaz', \%opts);
my $mrsfast = "$opts{b}/mrsfast";
my $samtools = "$opts{b}/samtools";
my $bindir = "$opts{b}";
my $javadir = "$opts{j}";

unless(defined($opts{'o'}) && defined($opts{'r'})){
	print "Improper number of arguments\n$usage";
	exit;
}

my $refindex = "$opts{r}.fai";
my @base = split(/\//, $opts{'i'});
my @fqname = split(/\./, $base[-1]);
#my $outbase = $fqname[0];
my $outbase = $opts{'c'};
# file names
my $sam_out = "$opts{o}/$outbase/$fqname[0].$fqname[1].$opts{n}.sam";
my $unmap = "$opts{o}/$outbase/$fqname[0].$fqname[1].$opts{n}.nohit";
my $bam_out = "$opts{o}/$outbase/$fqname[0].$fqname[1].$opts{n}.bam";
my $bam_sorted_out = "$opts{o}/$outbase/$fqname[0].$fqname[1].$opts{n}.sorted";
my $d_out = "$opts{o}/$outbase/divet";
my $sr_out = "$opts{o}/$outbase/split_read";
my $initial = "$opts{o}/$outbase/$fqname[0].$fqname[1].$opts{n}.sr.pe.data.sam";
my $working = "$opts{o}/$outbase/split_read/$fqname[0].$fqname[1].$opts{n}.sr";
my $divet = "$initial.divet.vh";

mkdir("$d_out") || print $! . "\n";
mkdir("$sr_out") || print $! . "\n";

my $initial_time = time;
# Creating split read fastqs



my $tstr = return_time_str($initial_time);
print "Running mrsfast on OEAs: $tstr\n";
system("$mrsfast --search $opts{r} --seq $working.split.fq -e $opts{e} -o $working.split.map -u $working.split.unmapped");
system("$bindir/split_match_sam -i $working.split.map -o $working.split.match -l $opts{x} -u $opts{x} -m $opts{m}");


my $tstr = return_time_str($initial_time);
print "Matching pairs; matching OEAs: $tstr\n";
print "$bindir/pair_match_v5 -i $initial.single.txt -fq $opts{i} -split $working.split.match -o $working.pair -l $opts{l} -u $opts{u} -m $opts{m}\n";
system("$bindir/pair_match_v5 -i $initial.single.txt -fq $opts{i} -split $working.split.match -o $working.pair -l $opts{l} -u $opts{u} -m $opts{m}");
system("cat $initial.single.txt $working.split.match.single.txt > $working.all.oea");
system("$bindir/SplitOEA_match -i $working.all.oea -o $working.all.oea.match -l $opts{l} -u $opts{u} -m $opts{m}");


my $tstr = return_time_str($initial_time);
print "Finalizing output to bams: $tstr\n";
system("$samtools view -bt $refindex $working.split.match > $working.split.match.bam");
system("$samtools sort $working.split.match.bam $working.split.match.sorted");
system("$samtools view -bt $refindex $working.split.match.single.txt > $working.split.match.single.bam");
system("$samtools sort $working.split.match.single.bam $working.split.match.single.sorted");
system("$samtools view -bt $refindex $initial.single.txt > $working.single.bam");
system("$samtools sort $working.single.bam $working.single.sorted");
system("rm $sam_out $unmap $bam_out");
system("rm $initial.inv.evert.txt $working.inv.evert.bam");
system("rm $initial.transchr.txt $working.transchr.bam");
system("rm $working.oea.name $working.oea.fq $working.split.map $working.split.unmapped");
system("rm $working.pair.single.txt $working.all.oea $working.all.oea.match.single.txt");
system("rm $working.split.match $working.split.match.bam $working.split.match.single.txt $working.split.match.single.bam $initial.single.txt $working.single.bam");
system("rm $working.all.oea.match.maxevert.txt $working.all.oea.match.trans.txt $working.pair.concordand.txt $working.pair.discordand.txt");
system("rm $working.inv.evert.sorted.bam $working.single.sorted.bam $working.split.fq $working.split.match.single.sorted.bam $working.split.match.sorted.bam");
system("rm $initial $working.pair.trans.txt");
if(-s "$working.all.oea.match" && -s "$working.pair"){
print "removing input interleaved fastq $opts{i}\n";
system("rm $opts{i}"); 
}else{
	print "$working.all.oea.match and $working.pair have size of zero!\n";
}
my $tstr = return_time_str($initial_time);
print "Finished with pipeline: $tstr\n";
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