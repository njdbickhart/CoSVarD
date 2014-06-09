#!/usr/bin/perl
# This script takes one or more merge_tab files and converts them into the tigra-sv input format
# Tigra SV format:
# [1]Chr [2]outside_start [3]inside_start [4]inside_end [5]outside_end [6]event_type [7]outside_size [8]aligner [9]comma_separated_sample_names [10]computational_approach [11]Group
# [1], [2],[3], [4], [5], [6] and [9] are essential

use strict;

if(scalar(@ARGV) < 1){
	print "$0 \<File 1\> \<File 2\> ...\n";
}
chomp (@ARGV);
my %event_holder;
foreach my $f (@ARGV){
	open (IN, "< $f") || die "could not open input file $f!\n";
	my ($sample_name) = $f =~ m/final_.+merge_(.+)\.tab/;
	while (my $line = <IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		my $size = $segs[4] - $segs[3];
		if($segs[3] > $segs[1] || $segs[4] < $segs[2] || $segs[5] =~ /inv/ || $size > 1000000){
			next; # These events are likely outside the reach of tigra-sv so will skip them
		}
		my @tags = split(/_/, $segs[5]);
		my @base_tags = split(/\./, $tags[0]);
		my @final_tags;
		foreach my $b (@base_tags){
			my $tmp;
			if($b =~ /gain/){
				$tmp = "ITX";
			}elsif($b =~ /del/){
				$tmp = "DEL";
			}elsif($b =~ /ins/){
				$tmp = "INS";
			}
			if(length($tmp) > 0){
				push(@final_tags, $tmp);
			}
		}
		my $last;
		my $skip = 0;
		for (my $x = 0; $x < scalar(@final_tags); $x++){
			if($x == 0){
				$last = $final_tags[$x];
				next;
			}else{
				if($final_tags[$x] eq $last){
					$skip = 1;
					last;
				}
			}
		}
		unless($skip){
			push(@{$event_holder{$sample_name}}, [$segs[0], $segs[3], $segs[1], $segs[2], $segs[4], $final_tags[0], $size, "MRSFAST", "ILLUMINA", $sample_name, $tags[1], "BFGL"]);
		}
		
	}
	close IN;
}

# Now to merge entries that overlap by at least 50%
my %final_events; #{chr}->[rows]->[rest...]
foreach my $sample (sort {$a cmp $b} keys(%event_holder)){
	foreach my $row (@{$event_holder{$sample}}){
		if(exists($final_events{$row->[0]})){
			my $found = 0;
			for (my $x = 0; $x < scalar(@{$final_events{$row->[0]}}); $x++){
				my $srow = $final_events{$row->[0]}->[$x];
				my $ovlp_bp = overlap($row->[1], $row->[4], $srow->[0], $srow->[3]);
				my $rolen = $row->[4] - $row->[1];
				my $solen = $srow->[3] - $srow->[0];
				my $rilen = $row->[3] - $row->[2];
				my $silen = $srow->[2] - $row->[1];
				if($ovlp_bp >= $rolen * 0.5 && $ovlp_bp >= $solen * 0.5 && $row->[5] eq $srow->[4]){
					my $iovlp_bp = overlap($row->[2], $row->[3], $srow->[1], $srow->[2]);
					if($iovlp_bp >= 0){
						$srow->[0] = least($row->[1], $srow->[0]);
						$srow->[1] = most($row->[2], $srow->[1]);
						$srow->[2] = least($row->[3], $srow->[2]);
						$srow->[3] = most($row->[4], $srow->[3]);
						$srow->[8] .= "," . $sample;
						$found = 1;
						last;
					}else{
						# This is probably a big mistake, but I am just going to take a dumb merger of the outer coords of the inner region
						$srow->[0] = least($row->[1], $srow->[0]);
						$srow->[1] = least($row->[2], $srow->[1]);
						$srow->[2] = most($row->[3], $srow->[2]);
						$srow->[3] = most($row->[4], $srow->[3]);
						$srow->[8] .= "," . $sample;
						$found = 1;
						last;
					}
				}
			}
			if(!$found){
				push(@{$final_events{$row->[0]}}, [@{$row}[1 .. 11]]);
			}
		}else{
			push(@{$final_events{$row->[0]}}, [@{$row}[1 .. 11]]);
		}
	}
}

# Printing out the input file
open(OUT, "> tigra_sv_input_file.tab") || die "could not create output!\n";
foreach my $k (sort{$a cmp $b} keys(%final_events)){
	foreach my $row (@{$final_events{$k}}){
		print OUT "$k\t" . join("\t", @{$row}) . "\n";
	}
}
close OUT;

exit;
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