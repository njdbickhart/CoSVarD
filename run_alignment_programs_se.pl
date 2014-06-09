#!/usr/bin/perl
# This script runs the BWA and MrsFAST alignment programs on the single end and mate pair split fastq files
# This will also use the Picard suite of programs in order to clean and filter the files

use strict;
use Getopt::Std;

my %opts;
my $usage = "$0 -i \<input fastq\>  -o \<output directory\> -r \<reference\> -j \<java path\> -v \<java memory\>
\t-e \<edit distance\> -n \<fq num\> -a \<unmasked reference\> -c \<outbase name\> -b \<bin dir\> -g \<read group number\>\n";

getopt('iojrenacbgv', \%opts);

if(scalar(keys(%opts)) < 11 || scalar(keys(%opts)) > 11){
	print "Incorrect number of options; All fields required!\n";
	print $usage;
	exit;
}

# Necessary variables
my $bindir = $opts{'b'};
my $outbase = $opts{'c'};
my $bnames = $opts{'c'};
my $fqrnum = $opts{'n'};
my $jopt = "-$opts{'v'}";
my $javapath = $opts{'j'};
my $splitedit = $opts{'e'} / 2;
my $rgnum = $opts{'g'};

# Output files
my $tsai = "$opts{o}/$outbase/$bnames.$fqrnum.se.sai";
my $reffai = "$opts{r}.fai";
my $mrsfastnohit = "$opts{o}/$outbase/$bnames.$fqrnum.mrsfast.nohit";
my $mrsfastsam = "$opts{o}/$outbase/$bnames.$fqrnum.se.mrsfast.sam";
my $mrsfastbam = "$opts{o}/$outbase/$bnames.$fqrnum.se.mrsfast.bam";
my $bwasam = "$opts{o}/$outbase/$bnames.$rgnum.$fqrnum.bwa.sam";
my $bwabam = "$opts{o}/$outbase/$bnames.$rgnum.$fqrnum.bwa.bam";
my $sortpre = "$opts{o}/$outbase/$bnames.$rgnum.$fqrnum.clean.sort";
my $sortbam = "$opts{o}/$outbase/$bnames.$rgnum.$fqrnum.clean.sort.bam";
my $cleaned = substr($bwabam, 0, length($bwabam) - 4) . ".clean";
my $reordered = "$cleaned.reorder.bam";
my $dup = "$opts{o}/$outbase/$bnames.$rgnum.$fqrnum.clean.sort.nodup.bam";
my $tmpdir = "$opts{o}/$outbase/tmp";

# Read group Information
my $sm = "$bnames.$opts{g}";
my $plat = "ILLUMINA";
my $rgstr = qq{\'\@RG\tID:$bnames\tPL:$plat\tSM:$sm\'};

# Ensuring that necessary directories are present
mkdir("$opts{o}") || print $! . "\n";
mkdir("$opts{o}/$outbase") || print $! . "\n";
mkdir("$tmpdir") || print "$tmpdir: $!\n";

# Running MrsFAST
#print "MrsFAST: $bnames $opts{1} $opts{2}\n";
system("$bindir/mrsfast --search $opts{r} --seq $opts{i} -e $opts{e} -o $mrsfastsam -u $mrsfastnohit");
system("rm $mrsfastnohit");
system("$bindir/samtools view -bS -t $reffai -o $mrsfastbam $mrsfastsam");
system("rm $mrsfastsam");

# Processing BWA bams 
#print "BWA: $bnames $opts{1} $opts{2}\n";
system("$bindir/bwa aln $opts{a} $opts{i} > $tsai");

system("$bindir/bwa samse -r $rgstr $opts{a} $tsai $opts{i} > $bwasam");

system("$bindir/samtools view -bS -o $bwabam $bwasam");

system("$javapath/java $jopt -jar $bindir/CleanSam.jar I=$bwabam O=$cleaned VALIDATION_STRINGENCY=LENIENT");
system("$javapath/java $jopt -jar $bindir/ReorderSam.jar I=$cleaned O=$reordered R=$opts{a} VALIDATION_STRINGENCY=LENIENT");

system("$bindir/samtools sort -m 2000000000 $reordered $sortpre");
system("rm $tsai $bwasam $bwabam $cleaned $reordered");

system("$javapath/java $jopt -jar $bindir/MarkDuplicates.jar I=$sortbam O=$dup M=metric TMP_DIR=$tmpdir REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT");
system("rm $sortbam");

# Check if the bwa and MrsFAST bams are present; if so, delete the initial fastq
if(-e $mrsfastbam && -e $dup){
	system("rm $opts{i}");
}