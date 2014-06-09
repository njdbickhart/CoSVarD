#!/usr/bin/perl
# This script takes the formatted output of the previous programs and processes them to generate a final bed file
# I will also try to incorporate a means to pass the output from this program into my drawing script
# output1: chr outsidestart outsideend cnvtypestr score + insidestart insideend color
# output2: chr, insidestart, insideend, outsidestart, outsideend, cnvtypstr, copynumber, confidence, docconsis, splitprob, vhprob

# cnvtypstr: [gain/loss/ins/inv/del]_[dsv]_cnvr#
# Bed color: 10+ red, 5-9 purple, 3-4 blue, 1 grey, 0 green, n/a black

# Inputs:
# doc bed: chr, start, end, gain/loss, copynumber
# vh bed: chr, insidestart, insideend, outsidestart, outsideend, ins/del/inv/, support, sumprob
# split bed: chr, start, end, del/inv, balanced support, unbalanced support, anchor edit, split edit

use strict;
use Getopt::Std;

#my $bedtools = "/ibfs7/asg2/bickhartd/bin/bedtools merge";
my $bedtools = "mergeBed";
my $usage = "$0 -d \<doc events\> -p \<paired end events\> -s \<split read events\> -c \<copynumber file\> -o \<output\>\n";
my %opts;
getopt('dpsoc', \%opts);

unless(defined($opts{'d'}) && defined($opts{'p'}) && defined($opts{'s'}) && defined($opts{'c'})){
	print $usage;
	exit;
}
our %colors = (
	'red' => '255,0,0',
	'purple' => '148,0,211',
	'blue' => '0,0,205',
	'grey' => '112,138,144',
	'green' => '0,255,0',
	'black' => '0,0,0',
	'orange' => '242,176,102');

my %docdata;
my %vhdata;
my %splitdata;
my %cnhash;

my %vhsplit_merger; #{chr}->[r][insidestart, insideend, outsidestart, outsideend, split call, vh call, splitsupport, vhsupport, vhsplitcn]
my %final_data; #{chr}->[r][insidestart, insideend, outsidestart, outsideend, mergecall, doccn, confidence, docconsist, $splitprob, $vhprob]

# Read in files
print "Reading files\n";
store_file_in_hash(\%cnhash, $opts{'c'});
store_file_in_hash(\%docdata, $opts{'d'});
store_file_in_hash(\%vhdata, $opts{'p'});
store_file_in_hash(\%splitdata, $opts{'s'});

# merge vh and split read; 
print "merging split and vh data\n";
merge_vh_split(\%vhdata, \%splitdata, \%vhsplit_merger, \%cnhash);

# merge doc with the vh/sr data; mark any outliers
print "incorporating doc data\n";
incorporate_doc(\%docdata, \%vhsplit_merger, \%final_data);

# Send outliers to drawing program
# To be implemented

# create merger files using inside information (just parity checking for now)
print "Checking merger integrity and parity\n";
open (OUT, "> $opts{o}.tmp");
foreach my $chr (keys(%final_data)){
	foreach my $rows (@{$final_data{$chr}}){
		print OUT "$chr\t$rows->[0]\t$rows->[1]\t$rows->[4]\n";
	}
}
close OUT;
open (OUT, "> $opts{o}.problem");
open (IN, "$bedtools -i $opts{o}.tmp -nms | ") || print "$!\n";
while(my $line = <IN>){
	chomp $line;
	my @segs = split(/\t/, $line); 
	my @nms = split(/;/, $segs[3]);
	if (scalar(@nms) > 1){
		print "Found problem in this region: $line\n";
		print OUT "$line\n";
	}
}
close IN;
close OUT;

# system("rm $opts{o}.tmp");

# print final output formats
# output1: chr outsidestart outsideend cnvtypestr score + insidestart insideend color
# output2: chr, insidestart, insideend, outsidestart, outsideend, cnvtypstr, copynumber, confidence, docconsis, splitprob, vhprob
print "Printing final output\n";
open (OBED, "> $opts{o}.bed");
open (OTAB, "> $opts{o}.tab");
my $i = 0; # iterator for cnv number
foreach my $chr (sort {my ($achrs) = $a =~ m/chr(.+)/; my ($bchrs) = $b =~ m/chr(.+)/; if($achrs =~ /X/){$achrs = 500;} if($bchrs =~ /X/){$bchrs = 500;} $achrs <=> $bchrs} keys(%final_data)){
	foreach my $rows (sort {$a->[2] <=> $b->[2]} @{$final_data{$chr}}){
		$i++;
		my $cn; # Bed color: 10+ red, 5-9 purple, 3-4 blue, 1 grey, 0 green, n/a black
		if ($rows->[5] eq "null"){
			$cn = $colors{'black'}; 
		}elsif($rows->[5] >= 10){
			$cn = $colors{'red'};
		}elsif($rows->[5] >= 5 && $rows->[5] < 10){
			$cn = $colors{'purple'};
		}elsif($rows->[5] >= 3 && $rows->[5] < 5){
			$cn = $colors{'blue'};
		}elsif($rows->[5] < 2 && $rows->[5] >= 1){
			$cn = $colors{'grey'};
		}elsif($rows->[5] < 1){
			$cn = $colors{'green'};
		}else{
			$cn = $colors{'black'};
		}
		my $sper = int($rows->[6] * 1000);
		print OBED "$chr\t$rows->[2]\t$rows->[3]\t$rows->[4]\_$i\t$sper\t+\t$rows->[0]\t$rows->[1]\t$cn\n";
		print OTAB "$chr\t$rows->[0]\t$rows->[1]\t$rows->[2]\t$rows->[3]\t$rows->[4]\_$i\t$rows->[5]\t$rows->[6]\t$rows->[7]\t$rows->[8]\t$rows->[9]\n";
	}
}
close OBED;
close OTAB;

exit;
# Here is my idea: I have several probability ratios from split read and vh, lets combine them into a unified metric
# for assessing the reliability of each call. 
sub vh_calc_prob {
	# returns the average probability of all paired end 
	# This is the probability that the read would have mapped here based on differences in sequence
	my ($support, $sumprob) = @_;
	if ($support == 0){return 0;}
	return $sumprob / $support;
	
}
sub split_calc_prob {
	my ($bsupport, $usupport, $aedit, $sedit) = @_;
	# To incorporate the bsupport and usupport variables in the probability calculation, I will have to 
	# read from the read depth estimates; to be incorporated later

	# Also, I will have to incorporate a probability-based phred value in splitread to make better anchor edit and
	# split edit calculations

	# This is a very crude, and likely useless approximation in the meantime
	return (1 - (($aedit / 50) + ($sedit / 25)));
	
}
sub placeholder_confidence_calc {
	# While I wait to implement the Bayes method, here is a placeholder method
	my ($docconsistent, $splitprob, $vhprob) = @_;
	if ($splitprob ne "null" && $vhprob ne "null"){
		# If the doc is inconsistent, return 0; if it is consistent, return the average of the splitprob and vhprob (one-based percentages)
		return (($docconsistent + ($splitprob + $vhprob)) / 3);
	}elsif ($splitprob ne "null" && $vhprob eq "null"){
		return (($docconsistent + ($splitprob)) / 2);
	}elsif ($splitprob eq "null" && $vhprob ne "null"){
		return (($docconsistent + ($vhprob)) / 2);
	}else{
		return 0;
	}
	
}
sub bayes_doc_vs_splitvh{
	# This is an adaptation of Bayes theorem to attempt to calculate a confidence score for 
	# any DOC related incompatibilites with the split-read and vh data

	# It takes all of the DOC ins/del events and tries to calculate the likelihood of a DOC event
	# in the presence of split/vh calls that are similar or different
}

sub merge_call_string {
	my ($doc, $sr, $vh) = @_;
	my $mergecall;
	if($doc ne "null" && $sr eq "null" && $vh ne "null"){
		$mergecall = $doc . "." . $vh . "_dv";
	}elsif ($doc ne "null" && $sr ne "null" && $vh eq "null"){
		$mergecall = $doc . "." . $sr . "_ds";
	}elsif ($doc ne "null" && $sr eq "null" && $vh eq "null"){
		$mergecall = $doc . "_d";
	}elsif ($doc eq "null" && $sr ne "null" && $vh ne "null"){
		$mergecall = $sr . "." . $vh . "_sv";
	}elsif ($doc eq "null" && $sr eq "null" && $vh ne "null"){
		$mergecall = $vh . "_v";
	}elsif ($doc eq "null" && $sr ne "null" && $vh eq "null"){
		$mergecall = $sr . "_s";
	}else{
		$mergecall = $doc . "." . $sr . "." . $vh . "_dsv";
	}
	return $mergecall;
}

sub incorporate_doc {
	my ($docref, $vsmerge, $vsdfinal) = @_;
	# doc bed: chr, start, end, gain/loss, copynumber
	# vsmerge: {chr}->[row]->[insidestart, insideend, outsidestart, outsideend, split call, vh call, splitsupport, vhsupport, vhcn]
	# vsdfinal: {chr}->[insidestart, insideend, outsidestart, outsideend, mergecall, doccn, confidence, docconsist, $splitprob, $vhprob]
	# cnvtypstr: [gain/loss/ins/inv/del]_[dsv]_cnvr#
	# Store all of the chrs from both files
	my %chrs;
	foreach my $c (keys(%{$docref})){
		$chrs{$c} += 1;
	}
	foreach my $c (keys(%{$vsmerge})){
		$chrs{$c} += 1;
	}
	foreach my $c (keys(%chrs)){
		# if found only in one chromosome, push onto the final dataset 5
		if ($chrs{$c} == 1){
			if(exists($docref->{$c})){
				foreach my $row (@{$docref->{$c}}){
					push(@{$vsdfinal->{$c}}, [$row->[0], $row->[1], $row->[0], $row->[1], $row->[2] . "_d", $row->[3], 1, 1, "null", "null"]);
				}
			}elsif(exists($vsmerge->{$c})){
				foreach my $row(@{$vsmerge->{$c}}){
					my $callavg = placeholder_confidence_calc(1, $row->[6], $row->[7]);
					my $mergecall = merge_call_string("null", $row->[4], $row->[5]);
					push(@{$vsdfinal->{$c}}, [$row->[0], $row->[1], $row->[2], $row->[3], $mergecall, $row->[8], $callavg, "null", $row->[6], $row->[7]]);
				}
			}else{
				# Shouldn't be here
				print "Error in chromosome reference: $c!\n";
			}
		}else{
			foreach my $rows (@{$vsmerge->{$c}}){
				my $found = 0;
				foreach my $srows (@{$docref->{$c}}){
					# vsdfinal: {chr}->[insidestart, insideend, outsidestart, outsideend, mergecall, doccn, confidence, docconsist, $splitprob, $vhprob]
					if ($rows->[2] < $srows->[1] && $rows->[3] > $srows->[0]){
						# Determine if doc data overlaps with vhmerge data based on outside edges
						my $doconsistency;
						my $ovlpbp = overlap($srows->[0], $srows->[1], $rows->[2], $rows->[3]);
						if(($srows->[2] eq "gain" && (!($rows->[4] =~ /del/) && ($rows->[5] =~ /ins/ || !($rows->[5] =~ /del/)))) ||
						($srows->[2] eq "loss" && ($rows->[4] =~ /del/ && $rows->[5] =~ /del/)) ||
						($rows->[4] =~ /^inv$/ && $rows->[5] =~ /^inv$/)){
							$doconsistency = 1;
						}else{
							$doconsistency = 0;
						}
						my $mergecall = merge_call_string($srows->[2], $rows->[4], $rows->[5]);
						my $dlen = $srows->[1] - $srows->[0];
						my $vslen = $rows->[3] - $rows->[2];
						my $avgprob = placeholder_confidence_calc($doconsistency, $rows->[6], $rows->[7]);
						# svh data completely overlaps doc data
						if($rows->[2] <= $srows->[0] && $rows->[3] >= $srows->[1]){
							my($nostart, $noend) = expand_outside_edges($rows->[2], $srows->[0], $rows->[3], $srows->[1]);
							push(@{$vsdfinal->{$c}}, [$rows->[0], $rows->[1], $nostart, $noend, $mergecall, $srows->[3], $avgprob, $doconsistency, $rows->[6], $rows->[7]]);
							push(@{$srows}, 1); # Found this one, so we do not have to print it separately
							$found = 1;
							last;
						}
						# Determine if doc data overlaps by at least 50% reciprocally, instead
						elsif( ($ovlpbp / $dlen >= 0.50) && ($ovlpbp / $vslen >= 0.50) && !($rows->[4] =~ /ins/)) {
							my($nostart, $noend) = expand_outside_edges($rows->[2], $srows->[0], $rows->[3], $srows->[1]);
							push(@{$vsdfinal->{$c}}, [$rows->[0], $rows->[1], $nostart, $noend, $mergecall, $srows->[3], $avgprob, $doconsistency, $rows->[6], $rows->[7]]);
							push(@{$srows}, 1); # Found this one, so we do not have to print it separately
							$found = 1;
							last;
						}
						
					}
				}
				if(!($found)){
					my $avgprob = placeholder_confidence_calc(1 , $rows->[6], $rows->[7]);
					my $mergecall = merge_call_string("null", $rows->[4], $rows->[5]);
					push(@{$vsdfinal->{$c}}, [$rows->[0], $rows->[1], $rows->[2], $rows->[3], $mergecall, $rows->[8], $avgprob, "null", $rows->[6], $rows->[7]]);
				}
			}
			foreach my $rows(@{$docref->{$c}}){
				# Now, go through and place all of the "unfound" doc calls into the final array
				if ($rows->[7] != 1 || $rows->[7] eq ""){
					push(@{$vsdfinal->{$c}}, [$rows->[0], $rows->[1], $rows->[0], $rows->[1], $rows->[2] . "_d", $rows->[3], 1, 1, "null", "null"]);
				}
			}
		}
	}
}
sub merge_vh_split {
	my ($vhref, $sref, $vsfinal, $cnhash) = @_;
	# vh bed: chr, insidestart, insideend, outsidestart, outsideend, ins/del/inv/, support, sumprob
	# split bed: chr, start, end, del/inv, balanced support, unbalanced support, anchor edit, split edit
	# vsfinal: {chr}->[row]->[insidestart, insideend, outsidestart, outsideend, split call, vh call, splitsupport, vhsupport
	# Store all of the chrs from both files
	my %chrs;
	foreach my $c (keys(%{$vhref})){
		$chrs{$c} += 1;
	}
	foreach my $c (keys(%{$sref})){
		$chrs{$c} += 1;
	}

	# cycle through each chr, start with vh data (should be more limited) and confirm with split read data
	foreach my $c (keys(%chrs)){
		# If only found in one dataset, push onto $vsfinal
		if ($chrs{$c} == 1){
			if(exists($vhref->{$c})){
				foreach my $row (@{$vhref->{$c}}){
					my $avgprob = vh_calc_prob($row->[5], $row->[6]);
					push(@{$vsfinal->{$c}}, [$row->[0], $row->[1], $row->[2], $row->[3], "null", $row->[4], "null", $avgprob]);
				}
			}elsif(exists($sref->{$c})){
				foreach my $row(@{$sref->{$c}}){
					my $slen = $row->[1] - $row->[0];
					if($slen > 500000){$row->[7] = 1; next;}
					my $splitprob = split_calc_prob($row->[3], $row->[4], $row->[5], $row->[6]);
					push(@{$vsfinal->{$c}}, [$row->[0], $row->[1], $row->[0], $row->[1], $row->[2], "null", $splitprob, "null"]);
				}
			}else{
				# Shouldn't be here
				print "Error in chromosome reference: $c!\n";
			}
		}else{
			foreach my $rows (@{$vhref->{$c}}){
				my $found = 0;
				foreach my $srows (@{$sref->{$c}}){
					if ($srows->[0] < $rows->[1] && $srows->[1] > $rows->[0]){
						# Determine if split read data overlaps completely with vh data
						my $slen = $srows->[1] - $srows->[0];
						if($slen > 500000){$srows->[7] = 1; next;}
						my $vlen = $rows->[1] - $rows->[0];
						if($srows->[0] >= $rows->[0] && $srows->[1] <= $rows->[1]){
							my $avgprob = vh_calc_prob($rows->[5], $rows->[6]);
							my $splitprob = split_calc_prob($srows->[3], $srows->[4], $srows->[5], $srows->[6]);
							my $vhevent = $rows->[4];
							my $spevent = $srows->[2];
							my($nistart, $niend) = refine_inside_edges($rows->[0], $srows->[0], $rows->[1], $srows->[1]);
							my $vhcn = determine_avg_cn($cnhash->{$c}, $rows->[2], $rows->[3]);
							push(@{$vsfinal->{$c}}, [$nistart, $niend, $rows->[2], $rows->[3], $spevent, $vhevent, $splitprob, $avgprob, $vhcn]);
							push(@{$srows}, 1); # Found this one, so we do not have to print it separately
							$found = 1;
							last;
						}
						# Determine if split read data overlaps by at least 50%, instead
						elsif((($srows->[0] <= $rows->[1] - ($vlen * 0.50)) && ($srows->[1] >= $rows->[0] + ($vlen * 0.50))) ||
						(($rows->[0] <= $srows->[1] - ($slen * 0.25)) && ($rows->[1] >= $srows->[0] + ($slen * 0.25)))){
							my $avgprob = vh_calc_prob($rows->[5], $rows->[6]);
							my $splitprob = split_calc_prob($srows->[3], $srows->[4], $srows->[5], $srows->[6]);
							my $vhevent = $rows->[4];
							my $spevent = $srows->[2];
							my($nistart, $niend) = refine_inside_edges($rows->[0], $srows->[0], $rows->[1], $srows->[1]);
							my $vhcn = determine_avg_cn($cnhash->{$c}, $rows->[2], $rows->[3]);
							push(@{$vsfinal->{$c}}, [$nistart, $niend, $rows->[2], $rows->[3], $spevent, $vhevent, $splitprob, $avgprob, $vhcn]);
							push(@{$srows}, 1); # Found this one, so we do not have to print it separately
							$found = 1;
							last;
						}
						
					}
				}
				if(!($found)){
					my $avgprob = vh_calc_prob($rows->[5], $rows->[6]);
					my $vhcn = determine_avg_cn($cnhash->{$c}, $rows->[2], $rows->[3]);
					push(@{$vsfinal->{$c}}, [$rows->[0], $rows->[1], $rows->[2], $rows->[3], "null", $rows->[4], "null", $avgprob, $vhcn]);
				}
			}
			foreach my $rows(@{$sref->{$c}}){
				# Now, go through and place all of the "unfound" split reads into the final array
				if ($rows->[7] != 1 || $rows->[7] eq ""){
					my $splitprob = split_calc_prob($rows->[3], $rows->[4], $rows->[5], $rows->[6]);
					my $vhcn = determine_avg_cn($cnhash->{$c}, $rows->[0], $rows->[1]);
					push(@{$vsfinal->{$c}}, [$rows->[0], $rows->[1], $rows->[0], $rows->[1], $rows->[2], "null", $splitprob, "null", $vhcn]);
				}
			}
		}
	}
}
sub determine_avg_cn{
	my ($cnregions, $start, $end) = @_;
	my @cns;
	foreach my $reg (@{$cnregions}){
		if(overlap($reg->[0], $reg->[1], $start, $end) > 0){
			push(@cns, $reg->[2]);
		}
	}
	if(scalar(@cns) == 0){return "null";}
	return average(\@cns);
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
sub least{
	my ($a, $b) = @_;
	return ($a > $b)? $b : $a;
}
sub most {
	my ($a, $b) = @_;
	return ($a > $b)? $a : $b;
}
# Shamelessly copied from bedtools
# Returns a number with negative numbers indicating no overlap and positive numbers indicating overlap
sub overlap {
	my ($s1, $e1, $s2, $e2) = @_;
	return least($e1, $e2) - most($s1, $s2);
}
sub return_longer_string {
	my ($st1, $st2) = @_;
	if(length($st1) >= length($st2)){
		return $st1;
	}else{
		return $st2;
	}
}
sub expand_outside_edges{
	my($s1, $s2, $e1, $e2) = @_;
	my $start = ($s1 <= $s2) ? $s1 : $s2;
	my $end = ($e1 >= $e2) ? $e1 : $e2;
	return $start, $end;
	
}
sub refine_inside_edges{
	my ($s1, $s2, $e1, $e2) = @_;
	my $start = ($s1 >= $s2) ? $s1 : $s2;
	my $end = ($e1 <= $e2) ? $e1 : $e2;
	if ($start > $end){
		# Shouldn't be here; problem!
		print "Warning, inside start $start greater than end $end\nFlipping\n";
		my $t = $end;
		$end = $start;
		$start = $t;
	}
	return $start, $end;
}
sub store_file_in_hash{
	my ($href, $file) = @_;
	chomp $file;
	open(IN, "< $file") || die "Could not open $file for initial hash!\n";
	while(my $line = <IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		if (scalar(@segs) < 3){next;}
		my $chr = shift(@segs);
		push(@{$href->{$chr}}, [@segs]);
	}
	close IN;
}

sub save_for_simple_drawing {
	my ($doc, $vh, $split, $cnvr, $outbase) = @_;
	# Not yet implemented
}