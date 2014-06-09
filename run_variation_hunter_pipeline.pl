#!/usr/bin/perl
# This is the subscript that will take input divets from the main pipeline and run them through variation hunter

use strict;
use Getopt::Std;

my $usage = "$0 -i \<location of divet files\> -n \<base name for ls command\> -m \<min length of insert\>
\t-M \<max length of insert\> -p \<preprocess pruning prob.\> -t \<event threshold\>
\t-c \<chr length file\> -g \<gap file\> -I \<iteration name\> -b \<bin_dir\>\n";



my %opts;
getopt('inmMptcgIb', \%opts);
my $bindir = "$opts{b}";

unless(defined($opts{'i'}) && defined($opts{'n'})){
	print "Improper number of arguments\n$usage";
	exit;
}

my @dfiles = `ls $opts{i}/$opts{n}*.$opts{I}*.vh`;
my $full_divet = "$opts{i}/$opts{n}.$opts{I}.full.cat.divet";

open(CHR, "< $opts{c}") || die "Could not open chr len file $opts{c}!\n";
my %chrlens;
while(my $line = <CHR>){
	chomp $line;
	my @segs = split(/\t/, $line);
	$chrlens{$segs[0]} = $segs[1];
}
close CHR;


# combine all divet files into a massive divet
print "Found: " . scalar(@dfiles) . " files\n";
open(OUT, "> $full_divet") || die "Could not open $full_divet for output!\n";
foreach my $d (@dfiles){
	open(IN, "< $d") || print "Could not open $d!\n";
	while(my $line = <IN>){
		print OUT $line;
	}
	close IN;
}
close OUT;
my $readname_u = "$opts{i}/$opts{n}.$opts{I}.readname.u";
my $readname_sort = "$opts{i}/$opts{n}.$opts{I}.readName.sort";
my $divet_name_edit_prob = "$opts{i}/$opts{n}.$opts{I}.name.edit.prob";
my $divet_minedit_totalp = "$opts{i}/$opts{n}.$opts{I}.name.minedit.totalprob";
my $divet_prob_full = "$opts{i}/$opts{n}.$opts{I}.prob.full";
my $divet_prob_clean = "$opts{i}/$opts{n}.$opts{I}.prob.clean";
my $full_chr_file = "$opts{i}/$opts{n}.$opts{I}.full.chr.rr_ff";
my $full_sv_events = "$opts{i}/$opts{n}.$opts{I}.final.SV";

# create a list of all the unique read names, sort them and print them sequentially in a text file
system("cut -f1 $full_divet | sort -u > $readname_u");
system("$bindir/sortString $readname_u > $readname_sort");

# create a text file containing just the name, edit distance and prob ($s[0], $s[9], $s[11]) from divet
system("cut -f1,10,12 $full_divet > $divet_name_edit_prob");

print "Read name sorted, starting probability filtering\n";
# calProbMinEditRead
system("$bindir/calProbMinEditRead $readname_sort $divet_name_edit_prob > $divet_minedit_totalp");

# convertToInvLRProbMinEditDist
system("$bindir/convertToInvLRProbMinEditDist $full_divet $divet_minedit_totalp $opts{p} > $divet_prob_full");

# Take output of last program and print select segments ($s[0-6], $s[9])
system("cut -f1,2,3,4,5,6,7,10 $divet_prob_full > $divet_prob_clean");

# Cycle through chromosomes and run createSetsDelAsInsNoGapInvRLProb.alpha and combineRR_FF.Prob on them
my @removal;
my @RR_FF;
print "Cycling through chromosomes to calculate SVs on per chromosome basis\n";
foreach my $chrs (keys(%chrlens)){
	open(OUT, "> $opts{i}/$chrs.$opts{n}.$opts{I}.tmp");
	print OUT "$chrs\t$chrlens{$chrs}\n";
	close OUT;
	print ("$bindir/createSetsDelAsInsNoGapInvRLProb.alpha $opts{i}/$chrs.$opts{n}.$opts{I}.tmp $divet_prob_clean $opts{i}/$opts{n}.$opts{I}.$chrs $opts{M} $opts{m} $opts{g}\n");
	system("$bindir/createSetsDelAsInsNoGapInvRLProb.alpha $opts{i}/$chrs.$opts{n}.$opts{I}.tmp $divet_prob_clean $opts{i}/$opts{n}.$opts{I}.$chrs $opts{M} $opts{m} $opts{g}");
	print ("$bindir/combineRR_FF.Prob $readname_sort $opts{i}/$opts{n}.$opts{I}.$chrs $opts{M} > $opts{i}/$opts{n}.$opts{I}.$chrs.RR_FF\n");
	system("$bindir/combineRR_FF.Prob $readname_sort $opts{i}/$opts{n}.$opts{I}.$chrs $opts{M} > $opts{i}/$opts{n}.$opts{I}.$chrs.RR_FF");
	push(@removal, "$opts{i}/$chrs.$opts{n}.$opts{I}.tmp", "$opts{i}/$opts{n}.$opts{I}.$chrs");
	push(@RR_FF, "$opts{i}/$opts{n}.$opts{I}.$chrs.RR_FF");
}

# Combine all individual chr output into a large file
open(OUT, "> $full_chr_file");
foreach my $rf (@RR_FF){
	open(IN, "< $rf") || print "could not open chr rf file $rf!\n";
	while(my $line = <IN>){
		print OUT $line;
	}
	close IN;
}
close OUT;

# Removing unnecessary files
push(@removal, @RR_FF);
foreach my $r (@removal){
	system("rm $r");
}
@removal = ();

# setCoverProbWeighted
print "Making final calls...\n";
system("$bindir/setCoverProbWeighted $readname_sort $full_chr_file $opts{t} > $full_sv_events");

# Grep out the different SV types by number
# Maybe I will skip that step until later

# Remove remaining unnecessary files
print "Removing unnecessary files\n";
system("rm $readname_u");
system("rm $readname_sort");
system("rm $divet_name_edit_prob");
system("rm $divet_minedit_totalp");
system("rm $divet_prob_full");
system("rm $divet_prob_clean");
system("rm $full_chr_file");


exit;