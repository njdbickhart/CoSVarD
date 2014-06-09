#!/usr/bin/perl
# This script is designed to take account of bam files in a folder and then merge and sort them 
# Given the limits on file entries for samtools merge, I will have to delegate a merge limit of 1000

use strict;
use Forks::Super;
use IPC::Open2;
use Getopt::Std;
use POSIX ":sys_wait_h";

my %opts;
getopt('libopf', \%opts);
my $usage = "$0 -l \<number of files to merge\> -i \<input path to search\> -b \<bin dir\> -o \<output name\> -p \<max processors\> -f \<BOOLEAN: is threaded?\>";
unless(defined($opts{'l'}) && defined($opts{'i'}) && defined($opts{'b'}) && defined($opts{'o'})){
	print $usage;
	exit;
}
my $limit = $opts{'l'};
my $samtools = "$opts{b}/samtools";
my $output = $opts{'o'};

# Adding wildcard
$opts{'i'} .= "*nodup.bam";

my @bams = <$opts{'i'}>;
chomp(@bams);

my @merger_strings;
my @removal_list;
my $merge_it = 0;
my $limit_it = 0;

my $skip = 0;
for (my $x = 0; $x < scalar(@bams); $x++){
	my ($wtr, $rdr, $err, $pid);
	$skip = 0;
	system("$samtools view -H $bams[$x] > tmp.out 2> tmp.err");
	#$pid = open2($wtr, $rdr, "$samtools view -H $bams[$x] | head");
	my $counter = 0;
	open(ERR, "< tmp.err"); 
	while (my $line = <ERR>){
		
		if ($line =~ /invalid BAM binary header/ || $line =~ /EOF marker is absent/){
			$skip = 1;
			last;
		}
	}
	close ERR;
	#if($x % 100 == 0){
		#sleep(2);
	#}
	#close $wtr;
	#close $rdr;
	if($skip){
		next;
	}
	if($x == 0 || $merger_strings[$merge_it] eq ''){
		$merger_strings[$merge_it] = $bams[$x];
		$removal_list[$merge_it] = $bams[$x];
		$limit_it = 1;
	}
	elsif($x != 0 && $limit_it >= $limit){
		$merge_it++;
		$merger_strings[$merge_it] = $bams[$x];
		$removal_list[$merge_it] = $bams[$x];
		$limit_it = 1;
	}else{
		$merger_strings[$merge_it] .= " " . $bams[$x];
		$removal_list[$merge_it] .= " " . $bams[$x];
		$limit_it++;
	}
	# if(waitpid($pid, WNOHANG) > 0){
		# my $exitstatus = $? >> 8;
		# print "bam file exit status: $pid\r";
	# }
}
#print "\n";
system("rm tmp.out tmp.err");

my @final_bams;
for (my $i = 0; $i < scalar(@merger_strings); $i++){
	my $tout = "$output.tmp.$i";
	#print "$samtools merge $tout $merger_strings[$i]\n";
	if(!($opts{'f'})){
		fork {max_proc => $opts{'p'} , cmd => "bsub -J bammerge$i -K -R \'rusage[mem=10000]\' -oo $output.bammerge.$i.OUT $samtools merge -f $tout $merger_strings[$i]"};
	}else{
		fork {max_proc => $opts{'p'}, cmd => "$samtools merge -f $tout $merger_strings[$i]"};
	}
	push(@final_bams, $tout);
}

my $pid = waitall();
if($pid > 0){
	print "Waited for $pid threads of samtools merger\n";
}

if(scalar(@final_bams) > 1){
	my $f_output_str = join(" ", @final_bams);
	print "$samtools merge -f $output $f_output_str\n";
	system("$samtools merge -f $output $f_output_str");
	print "$samtools index $output\n";
	system("$samtools index $output");
	foreach my $fb (@final_bams){
		push(@removal_list, $fb);
	}
}else{
	print "Only one file; moving, indexing...\n";
	system("mv $final_bams[0] $output");
	print "$samtools index $output\n";
	system("$samtools index $output");
}
open (OUT, "> $output.rm");
foreach my $r (@removal_list){
	print OUT "rm $r\n";
}
close OUT;
print "Run $output.rm to remove previous bams\n";
print "Finished\n";
exit;