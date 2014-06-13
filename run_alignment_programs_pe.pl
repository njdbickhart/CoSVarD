#!/usr/bin/perl
# This script runs the BWA and MrsFAST alignment programs on the split fastq files
# This will also use the Picard suite of programs in order to clean and filter the files

use strict;
use Getopt::Std;

my %opts;
my $usage = "$0 -1 \<input fastq1\> -2 \<input fastq2\> -o \<output directory\> -r \<reference\> -j \<java path\> -v \<java memory\>
\t-e \<edit distance\> -n \<fq num\> -a \<unmasked reference\> -c \<outbase name\> -b \<bin dir\> -g \<read group number\>
\t-l \<low bound\> -m \<max bound\> -u \<upper\>\n";

getopt('12ojrenacbgvlmu', \%opts);

if(scalar(keys(%opts)) < 15 || scalar(keys(%opts)) > 15){
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

# Output files
my $tsai1 = "$opts{o}/$outbase/$bnames.$fqrnum.1.sai";
my $tsai2 = "$opts{o}/$outbase/$bnames.$fqrnum.2.sai";
my $reffai = "$opts{r}.fai";
my $mrsfastnohit = "$opts{o}/$outbase/$bnames.$fqrnum.mrsfast.nohit";
my $mrsfastsam1 = "$opts{o}/$outbase/$bnames.$fqrnum.1.mrsfast.sam";
my $mrsfastsam2 = "$opts{o}/$outbase/$bnames.$fqrnum.2.mrsfast.sam";
my $mrsfastbam1 = "$opts{o}/$outbase/$bnames.$fqrnum.1.mrsfast.bam";
my $mrsfastbam2 = "$opts{o}/$outbase/$bnames.$fqrnum.2.mrsfast.bam";
my $bwasam = "$opts{o}/$outbase/$bnames.$fqrnum.bwa.sam";
my $bwabam = "$opts{o}/$outbase/$bnames.$fqrnum.bwa.bam";
my $sortpre = "$opts{o}/$outbase/$bnames.$fqrnum.clean.sort";
my $sortbam = "$opts{o}/$outbase/$bnames.$fqrnum.clean.sort.bam";
my $d_out = "$opts{o}/$outbase/divet";
my $sr_out = "$opts{o}/$outbase/split_read";
my $working = "$sr_out/$bnames.$fqrnum.sr";
my $cleaned = substr($bwabam, 0, length($bwabam) - 4) . ".clean";
my $reordered = "$cleaned.reorder.bam";
my $dup = "$opts{o}/$outbase/$bnames.$fqrnum.clean.sort.nodup.bam";
my $tmpdir = "$opts{o}/$outbase/tmp";
my $splitfq1 = "$working.split.fq";
my $mrsfastsplitsam1 = "$working.mrsfast.sam";
my $divet = "$working.divet.vh";
my $single = "$working.single.txt";
my $dfinal = "$d_out/$bnames.$fqrnum.divet.vh.gz";

# Read group Information
my $sm = "$bnames.$opts{g}";
my $plat = "ILLUMINA";
my $rgstr = qq{\'\@RG\tID:$sm\tPL:$plat\tLB:$sm\tSM:$bnames\'};

# Ensuring that necessary directories are present
mkdir("$opts{o}") || print $! . "\n";
mkdir("$opts{o}/$outbase") || print $! . "\n";
mkdir("$d_out") || print $! . "\n";
mkdir("$sr_out") || print $! . "\n";
mkdir("$tmpdir") || print "$tmpdir: $!\n";

# Running first round of MrsFAST
#print "MrsFAST: $bnames $opts{1} $opts{2}\n";
system("$bindir/mrsfast --search $opts{r} --seq $opts{1} -e $opts{e} -o $mrsfastsam1 -u $mrsfastnohit");
#system("rm $mrsfastnohit");
system("$bindir/mrsfast --search $opts{r} --seq $opts{2} -e $opts{e} -o $mrsfastsam2 -u $mrsfastnohit");
system("rm $mrsfastnohit");

# convert MrsFAST sams to bams
system("$bindir/samtools view -bS -t $reffai -o $mrsfastbam1 $mrsfastsam1");
system("$bindir/samtools view -bS -t $reffai -o $mrsfastbam2 $mrsfastsam2");


# Running read spliting pipeline
my $upperb = int($opts{u});
system("$javapath/java $jopt -jar $bindir/pairMatchMrsfastSam.jar -i1 $mrsfastsam1 -i2 $mrsfastsam2 -f1 $opts{1} -f2 $opts{2} -o $working -m $opts{m} -l $opts{l} -u $upperb");
system("rm $mrsfastsam1 $mrsfastsam2");

# Running last round of MrsFAST
system("$bindir/mrsfast --search $opts{r} --seq $splitfq1 -e $splitedit -o $mrsfastsplitsam1 -u $mrsfastnohit");
system("rm $mrsfastnohit");
system("rm splitfq1");
system("gzip $divet");
system("gzip $single");

# Processing BWA bams 
#print "BWA: $bnames $opts{1} $opts{2}\n";
system("$bindir/bwa aln $opts{a} $opts{1} > $tsai1");
system("$bindir/bwa aln $opts{a} $opts{2} > $tsai2");

system("$bindir/bwa sampe -r $rgstr $opts{a} $tsai1 $tsai2 $opts{1} $opts{2} > $bwasam");

system("$bindir/samtools view -bS -o $bwabam $bwasam");

system("$javapath/java $jopt -jar $bindir/CleanSam.jar I=$bwabam O=$cleaned VALIDATION_STRINGENCY=LENIENT");
system("$javapath/java $jopt -jar $bindir/ReorderSam.jar I=$cleaned O=$reordered R=$opts{a} VALIDATION_STRINGENCY=LENIENT");

system("$bindir/samtools sort -m 2000000000 $reordered $sortpre");
system("rm $tsai1 $tsai2 $bwasam $bwabam $cleaned $reordered");

system("$javapath/java $jopt -jar $bindir/MarkDuplicates.jar I=$sortbam O=$dup M=metric TMP_DIR=$tmpdir REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT");
system("rm $sortbam");

# Check for necessary files; if present, delete initial fastq; otherwise the checkpoint program will know something's wrong
if(-e "$single.gz" && -e $mrsfastbam1 && -e $mrsfastbam2 && -e "$divet.gz" && -e $dup){
	system("rm $opts{1} $opts{2}");
	system("mv $divet.gz $dfinal");
}