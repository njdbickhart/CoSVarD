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
use lib '/ibfs7/asg2/bickhartd/bin';
use lib '/home/dbickhart/share/bob_missou_data/lib';
use merger_utils;
use base_cnv_struct;
use doc_cnv_struct;
use pem_cnv_struct;
use split_cnv_struct;
use final_cnv_struct;
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

# Hashes of the structs. Hash based on chromosome
my %docdata;
my %vhdata;
my %splitdata;
# Hash of the copynumber values from input file
my %cnhash;

my $utils = merger_utils->new();

#my %vhsplit_merger; #{chr}->[r][insidestart, insideend, outsidestart, outsideend, split call, vh call, splitsupport, vhsupport, vhsplitcn]
my %final_data; #{chr}->[r]->final_cnv_struct instance

# Read in files
print "Reading files\n";
$utils->store_file_in_hash(\%cnhash, $opts{'c'});
store_file_object(\%docdata, $opts{'d'}, "DOC");
store_file_object(\%vhdata, $opts{'p'}, "VH");
store_file_object(\%splitdata, $opts{'s'}, "SPLIT");

# merge vh and split read; 
print "merging split and vh data\n";
merge_vh_split(\%vhdata, \%splitdata, \%final_data, \%cnhash);

# merge doc with the vh/sr data; mark any outliers
print "incorporating doc data\n";
incorporate_doc(\%docdata, \%final_data);

# Send outliers to drawing program
# To be implemented

# create merger files using inside information (just parity checking for now)
print "Checking merger integrity and parity\n";
open (OUT, "> $opts{o}.tmp");
foreach my $chr (keys(%final_data)){
	foreach my $rows (@{$final_data{$chr}}){
		my $mergecall = merge_call_string($rows->docCall(), $rows->splitCall(), $rows->vhCall());
		print OUT "$chr\t" . $rows->insStart() . "\t" . $rows->insEnd() . "\t$mergecall\n";
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
foreach my $chr (sort {
	my ($achrs) = $a =~ m/chr(.+)/; 
	my ($bchrs) = $b =~ m/chr(.+)/; 
	if($achrs =~ /X/){$achrs = 500;} 
	if($bchrs =~ /X/){$bchrs = 500;} 
	$achrs <=> $bchrs} keys(%final_data)){
		
	foreach my $rows (sort {$a->outStart() <=> $b->outStart()} @{$final_data{$chr}}){
		$i++;
		my $cn; # Bed color: 10+ red, 5-9 purple, 3-4 blue, 1 grey, 0 green, n/a black
		if ($rows->cn() == -1){
			$cn = $colors{'black'}; 
		}elsif($rows->cn() >= 10){
			$cn = $colors{'red'};
		}elsif($rows->cn() >= 5 && $rows->cn() < 10){
			$cn = $colors{'purple'};
		}elsif($rows->cn() >= 3 && $rows->cn() < 5){
			$cn = $colors{'blue'};
		}elsif($rows->cn() < 2 && $rows->cn() >= 1){
			$cn = $colors{'grey'};
		}elsif($rows->cn() < 1){
			$cn = $colors{'green'};
		}else{
			$cn = $colors{'black'};
		}
		$rows->eventNum($i);
		my $conf = placeholder_confidence_calc($rows->docConsist(), $rows->splitSupport(), $rows->vhSupport());
		my $sper = int( $conf * 1000);
		my $mergecall = merge_call_string($rows->docCall(), $rows->splitCall(), $rows->vhCall());
		my $tcn = ($rows->cn() == -1)? 'null' : $rows->cn();
		my $tsplit = ($rows->splitSupport() == -1)? 'null' : $rows->splitSupport();
		my $tvh = ($rows->vhSupport() == -1)? 'null' : $rows->vhSupport();
		print OBED "$chr\t" . $rows->outStart() . "\t" . $rows->outEnd() . "\t$mergecall\_$i\t$sper\t+\t" . $rows->insStart() . "\t" . $rows->insEnd() . "\t$cn\n";
		print OTAB "$chr\t" . $rows->insStart() . "\t" . $rows->insEnd() . "\t" . $rows->outStart() . "\t" . $rows->outEnd() . "\t$mergecall\_$i\t" . $tcn 
		. "\t$conf\t" . $rows->docConsist() . "\t" . $tsplit . "\t" . $tvh ."\n";
	}
}
close OBED;
close OTAB;

exit;
# Here is my idea: I have several probability ratios from split read and vh, lets combine them into a unified metric
# for assessing the reliability of each call. 
sub store_file_object{
	my ($href, $file, $obtype) = @_;
	# Inputs:
	# doc bed: chr, start, end, gain/loss, copynumber
	# vh bed: chr, insidestart, insideend, outsidestart, outsideend, ins/del/inv/, support, sumprob
	# split bed: chr, start, end, del/inv, balanced support, unbalanced support, anchor edit, split edit
	require base_cnv_struct;
	require doc_cnv_struct;
	require pem_cnv_struct;
	require split_cnv_struct;
	
	open (IN, "< $file") || die "Could not open $file $obtype!\n";
	if ($obtype eq "DOC"){
		while (my $line = <IN>){
			chomp $line;
			my @segs = split(/\t/, $line);
			my $call = ($segs[3] eq "loss")? "del" : "gain";
			if ($segs[4] eq "null" || $segs[4] eq '' || !defined($segs[4])){
				$segs[4] = -1;
			}
			push(@{$href->{$segs[0]}}, doc_cnv_struct->new(
				'chr' => $segs[0], 
				'insStart' => $segs[1], 
				'insEnd' => $segs[2],
				'type' => 'd',
				'call' => $call,
				'cn' => $segs[4]));
		}
	}elsif($obtype eq "VH"){
		while (my $line = <IN>){
			chomp $line;
			my @segs = split(/\t/, $line);
			push(@{$href->{$segs[0]}}, pem_cnv_struct->new(
				'chr' => $segs[0], 
				'insStart' => $segs[1], 
				'insEnd' => $segs[2],
				'type' => 'v',
				'call' => $segs[5],
				'outStart' => $segs[3],
				'outEnd' => $segs[4],
				'support' => $segs[6],
				'sumProb' => $segs[7]));
		}
	}elsif($obtype eq "SPLIT"){
		while (my $line = <IN>){
			chomp $line;
			my @segs = split(/\t/, $line);
			push(@{$href->{$segs[0]}}, split_cnv_struct->new(
				'chr' => $segs[0], 
				'insStart' => $segs[1], 
				'insEnd' => $segs[2],
				'type' => 's',
				'call' => $segs[3],
				'balSupport' => $segs[4],
				'unbalSupport' => $segs[5],
				'anchorEdit' => $segs[6],
				'splitEdit' => $segs[7]));
		}
	}
	close IN;
}
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
	if ($splitprob != -1 && $vhprob != -1){
		# If the doc is inconsistent, return 0; if it is consistent, return the average of the splitprob and vhprob (one-based percentages)
		return (($docconsistent + ($splitprob + $vhprob)) / 3);
	}elsif ($splitprob != -1 && $vhprob != -1){
		return (($docconsistent + ($splitprob)) / 2);
	}elsif ($splitprob != -1 && $vhprob != -1){
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
	my ($docref, $vsdfinal) = @_;
	# doc bed: chr, start, end, gain/loss, copynumber
	# vsmerge: {chr}->[row]->[insidestart, insideend, outsidestart, outsideend, split call, vh call, splitsupport, vhsupport, vhcn]
	# vsdfinal: {chr}->[insidestart, insideend, outsidestart, outsideend, mergecall, doccn, confidence, docconsist, $splitprob, $vhprob]
	# cnvtypstr: [gain/loss/ins/inv/del]_[dsv]_cnvr#
	# Store all of the chrs from both files
	my %chrs;
	require merger_utils;
	require base_cnv_struct;
	require doc_cnv_struct;
	require pem_cnv_struct;
	require split_cnv_struct;
	require final_cnv_struct;

	my $utils = merger_utils->new();
	foreach my $c (keys(%{$docref})){
		$chrs{$c} += 1;
	}
	foreach my $c (keys(%{$vsdfinal})){
		$chrs{$c} += 1;
	}
	foreach my $c (keys(%chrs)){
		# if found only in one chromosome, push onto the final dataset 5
		if ($chrs{$c} == 1){
			if(exists($docref->{$c})){
				foreach my $row (@{$docref->{$c}}){
					push(@{$vsdfinal->{$c}}, final_cnv_struct->new(
						'chr' => $c,
						'insStart' =>$row->insStart(), 
						'insEnd' => $row->insEnd(), 
						'outStart' => $row->outStart(), 
						'outEnd' => $row->outEnd(), 
						'call' => $row->call(), 
						'type' => $row->type(),
						'cn' => $row->cn(),
						'docCall' => $row->call()));
				}
			}else{
				# Shouldn't be here
				print "Error in chromosome reference: $c!\n";
			}
		}else{
			foreach my $rows (@{$vsdfinal->{$c}}){
				my $found = 0;
				foreach my $srows (@{$docref->{$c}}){
					if($srows->used()){next;} # Found this one already; skipping
					if ($rows->outStart() < $srows->outEnd() && $rows->outEnd() > $srows->outStart()){
						# Determine if doc data overlaps with vhmerge data based on outside edges
						my $doconsistency;
						my $outovlpbp = $utils->overlap($srows->outStart(), $srows->outEnd(), $rows->outStart(), $rows->outEnd());
						my $inovlpbp = $utils->overlap($srows->insStart(), $srows->insEnd(), $rows->insStart(), $rows->insEnd());
						if(($srows->call() eq "gain" && (!($rows->splitCall() =~ /del/) && ($rows->splitCall() =~ /ins/ || !($rows->vhCall() =~ /del/)))) ||
						($srows->call() eq "loss" && ($rows->splitCall() =~ /del/ && $rows->vhCall() =~ /del/)) ||
						($rows->splitCall() =~ /inv/ && $rows->vhCall() =~ /inv/)){
							$doconsistency = 1;
						}else{
							$doconsistency = 0;
						}
						#my $mergecall = merge_call_string($srows->[2], $rows->[4], $rows->[5]);
						my $dolen = $srows->outEnd() - $srows->outStart();
						my $vsolen = $rows->outEnd() - $rows->outStart();
						my $dilen = $srows->insEnd() - $srows->insStart();
						my $vsilen = $rows->insEnd() - $rows->insStart();
						#my $avgprob = placeholder_confidence_calc($doconsistency, $rows->[6], $rows->[7]);
						# svh data overlaps doc data on inside edges at a greater than 50% rate reciprocally
						if(($inovlpbp / $dilen >= 0.50) && ($inovlpbp / $vsilen >= 0.50)){
							my($nistart, $niend) = $utils->refine_inside_edges($rows->insStart(), $srows->insStart(), $rows->insEnd(), $srows->insEnd());
							my($nostart, $noend) = $utils->expand_outside_edges($rows->outStart(), $srows->outStart(), $rows->outEnd(), $srows->outEnd());
							my $typestr = ($rows->type() eq 'sv')? 'dsv' : ($rows->type() eq 's')? 'ds' : 'dv';
								$rows->chr($c);
								$rows->insStart($nistart); 
								$rows->insEnd($niend); 
								$rows->outStart($nostart); 
								$rows->outEnd($noend); 
								$rows->call('mult'); 
								$rows->type($typestr);
								$rows->cn($srows->cn());
								$rows->docConsist($doconsistency);
								$rows->docCall($srows->call());
							$srows->used(1); # Found this one, so we do not have to print it separately
							$found = 1;
							last;
						}
						# Determine if doc data overlaps by at least 50% reciprocally
						elsif( ($outovlpbp / $dolen >= 0.50) && ($outovlpbp / $vsolen >= 0.50) && !($rows->splitCall()=~ /ins/)) {
							my($nostart, $noend) = $utils->expand_outside_edges($rows->outStart(), $srows->outStart(), $rows->outEnd(), $srows->outEnd());
							my $typestr = ($rows->type() eq 'sv')? 'dsv' : ($rows->type() eq 's')? 'ds' : 'dv';
								$rows->chr($c);
								$rows->insStart($rows->insStart()); 
								$rows->insEnd($rows->insEnd()); 
								$rows->outStart($nostart); 
								$rows->outEnd($noend); 
								$rows->call('mult'); 
								$rows->type($typestr);
								$rows->cn($srows->cn());
								$rows->docConsist($doconsistency);
								$rows->docCall($srows->call());
								
							$srows->used(1); # Found this one, so we do not have to print it separately
							$found = 1;
							last;
						}
						
					}
				}
				
			}
			foreach my $rows(@{$docref->{$c}}){
				# Now, go through and place all of the "unfound" doc calls into the final array
				if (!($rows->used())){
					push(@{$vsdfinal->{$c}}, final_cnv_struct->new(
						'chr' => $c,
						'insStart' =>$rows->insStart(), 
						'insEnd' => $rows->insEnd(), 
						'outStart' => $rows->outStart(), 
						'outEnd' => $rows->outEnd(), 
						'call' => $rows->call(), 
						'type' => $rows->type(),
						'cn' => $rows->cn(),
						'docCall' => $rows->call()));
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
	require merger_utils;
	require base_cnv_struct;
	require doc_cnv_struct;
	require pem_cnv_struct;
	require split_cnv_struct;
	require final_cnv_struct;
	
	foreach my $c (keys(%{$vhref})){
		$chrs{$c} += 1;
	}
	foreach my $c (keys(%{$sref})){
		$chrs{$c} += 1;
	}
	my $utils = merger_utils->new();
	my %vsfinal; # {chr}->[row]->final_cnv_struct instance
	# cycle through each chr, start with vh data (should be more limited) and confirm with split read data
	foreach my $c (keys(%chrs)){
		# If only found in one dataset, push onto $vsfinal
		if ($chrs{$c} == 1){
			if(exists($vhref->{$c})){
				foreach my $row (@{$vhref->{$c}}){
					my $avgprob = vh_calc_prob($row->support(), $row->sumProb());
					my $vhcn = $utils->determine_avg_cn($cnhash->{$c}, $row->outStart(), $row->outEnd());
					push(@{$vsfinal->{$c}}, final_cnv_struct->new(
						'chr' => $c,
						'insStart' =>$row->insStart(), 
						'insEnd' => $row->insEnd(), 
						'outStart' => $row->outStart(), 
						'outEnd' => $row->outEnd(), 
						'call' => $row->call(), 
						'type' => $row->type(),
						'cn' => $vhcn,
						'vhCall' => $row->call(),
						'vhSupport' => $avgprob));
				}
			}elsif(exists($sref->{$c})){
				foreach my $row(@{$sref->{$c}}){
					my $slen = $row->insEnd() - $row->insStart();
					if($slen > 500000){$row->used(1); next;}
					my $splitprob = split_calc_prob($row->balSupport(), $row->unbalSupport(), $row->anchorEdit(), $row->splitEdit());
					push(@{$vsfinal->{$c}}, final_cnv_struct->new(
						'chr' => $c,
						'insStart' =>$row->insStart(), 
						'insEnd' => $row->insEnd(), 
						'outStart' => $row->insStart(), 
						'outEnd' => $row->insEnd(), 
						'call' => $row->call(), 
						'type' => $row->type(),
						'splitCall' => $row->call(),
						'splitSupport' => $splitprob));
				}
			}else{
				# Shouldn't be here
				print "Error in chromosome reference: $c!\n";
			}
		}else{
			foreach my $rows (@{$vhref->{$c}}){
				my $found = 0;
				foreach my $srows (@{$sref->{$c}}){
					if ($srows->insStart() < $rows->insEnd() && $srows->insEnd() > $rows->insStart()){
						# Determine if split read data overlaps completely with vh data
						my $slen = $srows->insEnd() - $srows->insStart();
						if($slen > 500000){$srows->used(1); next;}
						my $vlen = $rows->insEnd() - $rows->insStart();
						if($srows->insStart() >= $rows->insStart() && $srows->insEnd() <= $rows->insEnd() && !($srows->used())){
							my $avgprob = vh_calc_prob($rows->support(), $rows->sumProb());
							my $splitprob = split_calc_prob($srows->balSupport(), $srows->unbalSupport(), $srows->anchorEdit(), $srows->splitEdit());
							my($nistart, $niend) = $utils->refine_inside_edges($rows->insStart(), $srows->insStart(), $rows->insEnd(), $srows->insEnd());
							my $vhcn = $utils->determine_avg_cn($cnhash->{$c}, $rows->outStart(), $rows->outEnd());
							push(@{$vsfinal->{$c}}, final_cnv_struct->new(
								'chr' => $c,
								'insStart' =>$nistart, 
								'insEnd' => $niend, 
								'outStart' => $rows->outStart(), 
								'outEnd' => $rows->outEnd(), 
								'call' => 'mult', 
								'type' => 'sv',
								'cn' => $vhcn,
								'vhSupport' => $avgprob,
								'splitSupport' => $splitprob,
								'splitCall' => $srows->call(),
								'vhCall' => $rows->call()));
							$srows->used(1); # Found this one, so we do not have to print it separately
							$found = 1;
							last;
						}
						# Determine if split read data overlaps by at least 50%, instead
						elsif((($srows->insStart() <= $rows->insEnd() - ($vlen * 0.50)) && ($srows->insEnd() >= $rows->insStart() + ($vlen * 0.50))) ||
						(($rows->insStart() <= $srows->insEnd() - ($slen * 0.25)) && ($rows->insEnd() >= $srows->insStart() + ($slen * 0.25))) && !($srows->used())){
							my $avgprob = vh_calc_prob($rows->support(), $rows->sumProb());
							my $splitprob = split_calc_prob($srows->balSupport(), $srows->unbalSupport(), $srows->anchorEdit(), $srows->splitEdit());
							#my $vhevent = $rows->[4];
							#my $spevent = $srows->[2];
							my($nistart, $niend) = $utils->refine_inside_edges($rows->insStart(), $srows->insStart(), $rows->insEnd(), $srows->insEnd());
							my $vhcn = $utils->determine_avg_cn($cnhash->{$c}, $rows->outStart(), $rows->outEnd());
							push(@{$vsfinal->{$c}}, final_cnv_struct->new(
								'chr' => $c,
								'insStart' =>$nistart, 
								'insEnd' => $niend, 
								'outStart' => $rows->outStart(), 
								'outEnd' => $rows->outEnd(), 
								'call' => 'mult', 
								'type' => 'sv',
								'cn' => $vhcn,
								'vhSupport' => $avgprob,
								'splitSupport' => $splitprob,
								'splitCall' => $srows->call(),
								'vhCall' => $rows->call()));
							$srows->used(1); # Found this one, so we do not have to print it separately
							$found = 1;
							last;
						}
						
					}
				}
				if(!($found)){
					my $avgprob = vh_calc_prob($rows->support(), $rows->sumProb());
					my $vhcn = $utils->determine_avg_cn($cnhash->{$c}, $rows->outStart(), $rows->outEnd());
					push(@{$vsfinal->{$c}}, final_cnv_struct->new(
						'chr' => $c,
						'insStart' =>$rows->insStart(), 
						'insEnd' => $rows->insEnd(), 
						'outStart' => $rows->outStart(), 
						'outEnd' => $rows->outEnd(), 
						'call' => $rows->call(), 
						'type' => 'v',
						'cn' => $vhcn,
						'vhSupport' => $avgprob,
						'vhCall' => $rows->call()));
				}
			}
			foreach my $rows(@{$sref->{$c}}){
				# Now, go through and place all of the "unfound" split reads into the final array
				if (!($rows->used())){
					my $splitprob = split_calc_prob($rows->balSupport(), $rows->unbalSupport(), $rows->anchorEdit(), $rows->splitEdit());
					push(@{$vsfinal->{$c}}, final_cnv_struct->new(
						'chr' => $c,
						'insStart' =>$rows->insStart(), 
						'insEnd' => $rows->insEnd(), 
						'outStart' => $rows->insStart(), 
						'outEnd' => $rows->insEnd(), 
						'call' => $rows->call(), 
						'type' => 's',
						'splitSupport' => $splitprob,
						'splitCall' => $rows->call()));
				}
			}
		}
	}
}