#!/usr/bin/perl
# This is the remake of the Alkan wssd pipeline C wrapper
# It is designed to be modular, and to omit the initial r^2 statistical analysis done by alkan
# The goals are to calculate statistics for some of the files, then pipe them into
# the shell script. 

use strict;
use Getopt::Long;

my $usage = " 
Usage: $0 [files]
	
  Files:
	--File1		= 5kb absolute sliding windows with hits
	--File2		= 1kb absolute sliding windows with hits
	--File3		= 1kb non-overlapping windows with hits
	--bed		= required with --hits option; output base name for bedfiles (same as Files 1-3)
	--hits		= required with --bed option; folder containing bedfiles
	--ref		= required with --hits option; reference genome fasta
	\n\n";
my $file1;

my $file2;
my $file3;
my $intersect = 0;
my $hits;
my $bedtools = "/mnt/gliu1_usb/dbickhart/BEDTools-Version-2.10.0/bin";
my $template1 = "/mnt/gliu1_usb/dbickhart/alkan_files/umd3/umd3_windows/umd3_template_file1.bed";
my $template2 = "/mnt/gliu1_usb/dbickhart/alkan_files/umd3/umd3_windows/umd3_template_file2.bed";
my $template3 = "/mnt/gliu1_usb/dbickhart/alkan_files/umd3/umd3_windows/umd3_template_file3.bed";
my $bedout;
my $reference; my $java; my $bin;

my $result = GetOptions(
	'File1=s' => \$file1,
	'File2=s' => \$file2,
	'File3=s' => \$file3,
	'bed=s' => \$bedout,
	'hits=s' => \$hits,
	'ref=s' => \$reference,
	'java=s' => \$java,
	'bin=s' => \$bin,
	"help|?|h" => sub{print $usage; exit;}
);
unless(defined($file1)&& defined($file2) && defined($file3)){
	print $usage;
	exit;
}
chomp $file1;
chomp $file2;
chomp $file3;
chomp $hits;
if (defined($bedout) && defined($hits) && defined($reference)){
	print "Reading files\n";
	system("$java/java -cp $bin/lib/sam-1.53.jar -jar $bin/mergeDocWindows.jar -I $hits -O $bedout -R $reference -n 8");
}

my $file1_c = substr($file1, 0, (length($file1) - 4));
my $file3_c = substr($file3, 0, (length($file3) - 4));

check_files_exist($file1, $file2, $file3);
calculate_pipe("$file1_c\_c.bed", $file1, $file2, "$file3_c\_c.bed", $file3);

exit;
sub check_files_exist{
	my @files = @_;
	foreach my $f (@files){
		open(IN, "< $f") || die "Could not open $f!\n";
		close IN;
	}
	return;
}
sub calculate_pipe (){
	my ($f1c, $f1, $f2, $f3c, $f3) = @_;
	#$file1_c should be used for the calculations
	open (IN, "< $f1") || die "Could not open $f1!\n";
	print "Loading numbers from $f1c into array...\n";
	my @hit_numbers;
	while (my $line = <IN>){
		my @segments = split(/\t/, $line);
		chomp $segments[4];
		push (@hit_numbers, $segments[4]);
	}
	close IN;
	print "Calculating average...\n";
	my $num_2 = average(\@hit_numbers);
	print "Calculating stdev...\n";
	my $st_dev = standard_deviation(\@hit_numbers);
	my $num_3 = $num_2 + ($st_dev * 3);
	my $num_4 = $num_3 / 5;
	my $num_7 = $num_2 - ($st_dev * 2);
	my $num_8 = $num_7 / 5;
	my $num_9 = $num_2 + ($st_dev * 2);
	my $num_10 = $num_2 / 5;
	print "Average was: $num_2\n";
	print "Stdev was; $st_dev\n";
	#$1 BAC doc file name
	#$2 unique bacs average DOC
	#$3 unique bacs avg+3stdev
	#$4 unique bacs (avg+3stdev)/5
	#$5 HG17 5K DOC file name
	#$6 HG17 1K DOC file name
	#$7 unique bacs avg-2stdev
	#$8 unique bacs (avg-2stdev)/5
	#$9 unique bacs avg+2stdev
	#$10 unique bacs avg/5
	#$11 BAC nonoverlapping 1k doc file name
	#$12 HG17 nonoverlapping 1k doc file name
	print "Running following command:\n";
	print "sh $bin/cattle_umd3_pipeline_bob_placeholder.sh $f1c $num_2 $num_3 $num_4 $f1 $f2 $num_7 $num_8 $num_9 $num_10 $f3c $f3\n\n";
	system("sh $bin/cattle_umd3_pipeline_bob_placeholder.sh $f1c $num_2 $num_3 $num_4 $f1 $f2 $num_7 $num_8 $num_9 $num_10 $f3c $f3");
}

sub standard_deviation {
	my $array_ref = shift(@_);
	#Prevent division by 0 error in case you get junk data
	my @numbers = @{$array_ref};
	return undef unless(scalar(@numbers));

	# Step 1, find the mean of the numbers
	my $total1 = 0;
	foreach my $num (@numbers) {
	$total1 += $num;
	}
	my $mean1 = $total1 / (scalar @numbers);

	# Step 2, find the mean of the squares of the differences
	# between each number and the mean
	my $total2 = 0;
	foreach my $num (@numbers) {
	$total2 += ($mean1-$num)**2;
	}
	my $mean2 = $total2 / (scalar @numbers);

	# Step 3, standard deviation is the square root of the
	# above mean
	my $std_dev = sqrt($mean2);
	return $std_dev;
}
sub average {
	my $array_r = shift(@_);
	my @numb = @{$array_r};
	my $total3 = 0;
	foreach my $num1 (@numb) {
	$total3 += $num1;
	}
	my $mean3 = $total3 / (scalar @numb);
	return $mean3;
}
