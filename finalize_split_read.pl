#!/usr/bin/perl
# this script is designed to finish the split_read pipeline and generate the final results
# It concatenates the pair and all.oea.match files and then uses the split read included scripts to 
# generate the final output.

use strict;
use Getopt::Std;

my %opts;
my $usage = "$0 -i \<split file folder\> -o \<final file output\> -b \<base name\> -B \<bin dir\> -j \<java dir\>\n";  

getopt('iobBj', \%opts);
unless(defined($opts{'i'}) && defined($opts{'o'})){
	print $usage;
	exit;
}

my $bindir = "$opts{B}";
my $javadir = "$opts{j}";
# Concatenate the files

my $paircat = "$opts{i}/$opts{b}.sr.pair";
my $matchcat = "$opts{i}/$opts{b}.sr.all.oea.match";
system("rm $paircat");
system("rm $matchcat");

my @pair_files = `ls $opts{i}/*.sr.pair`;
my @match_files = `ls $opts{i}/*.sr.all.oea.match`;

chomp(@pair_files);
chomp(@match_files);

my @tpair_files = grep{$_ ne $paircat} @pair_files;
my @tmatch_files = grep{$_ ne $matchcat} @match_files;

@pair_files = @tpair_files;
@match_files = @tmatch_files;

my $pairstr = join(" ", @pair_files);
my $matchstr = join(" ", @match_files);

my $pairformat = "$opts{i}/$opts{b}.sr.pair.format";
my $matchformat = "$opts{i}/$opts{b}.sr.all.oea.match.format";
my $eventsout = "$opts{i}.total.$opts{b}.sr.events";

print "Concatenating files\n";


system("$bindir/simple_file_merger.pl $pairstr $paircat");
system("$bindir/simple_file_merger.pl $matchstr $matchcat");

# Format them using the supplied scripts
print "Formating files $paircat and $matchcat\n";
system("$bindir/pairToSC.sh $paircat $pairformat $matchcat $matchformat");

# Run the "events" program; will have to optimize this later, I think
print "Running the events5 program\n";
#system("$bindir/events5 $pairformat $matchformat $eventsout 50");
system("$javadir/java -Xmx9G -jar $bindir/splitReadFinalEvents.jar $pairformat $matchformat $eventsout 50 4");

# run the final output script to process the end results
print "Final processing...\n";
system("$bindir/getFinalevents.sh $eventsout $opts{o}.raw");

open(IN, "< $opts{o}.raw");
open(OUT, "> $opts{o}");

# sensitivity and specificity filter proposed in split read supplemental
while(my $line = <IN>){
	chomp $line;
	my @segs = split(/\t/, $line); 
	if($segs[5] >= 2 && $segs[6] >= 2){
		print OUT join("\t", @segs) . "\n";
	}
}
close IN;
close OUT;

exit;  