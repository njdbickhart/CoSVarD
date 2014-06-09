#!/usr/bin/perl
# This script is designed to take a series of final variationhunter output files and merge the coordinates
# It will produce the final output indicated in the main alignment pipeline
# vh bed: chr, insidestart, insideend, outsidestart, outsideend, ins/del/inv/, support, sumprob

use strict;

my $out = pop(@ARGV);
my @input = @ARGV;
chomp @input;
our %svtypes = (
	1 => "ins",
	2 => "del",
	3 => "inv",
	4 => "inv",
	5 => "inv");

if(scalar(@input) == 1){
	print "Found only one file, generating vh merger output\n";
	single_generate_final($input[0], $out);
}else{
	print "Found " . scalar(@input) . " files\n";
	merge_multiple(\@input, $out);
}

exit;
sub merge_multiple {
	my ($fileref, $out) = @_;
	my %storage; # {chr}->[row]->[is, ie, os, oe, sv, sup, prob]
	foreach my $f (@{$fileref}){
		open(IN, "< $f") || print "Could not open $f!\n";			
		while (my $line = <IN>){
			if(!($line =~ /Inside/)){next;} # clearing the header junk from the file
			chomp $line;
			my ($os, $is, $ie, $oe, $c, $sv, $sup, $avg, $prob) = $line =~ m/Inside_Start:(\d+) Inside_End:(\d+) OutSide_Start:(\d+) Oustide_End:(\d+) chro:(.+) SVtype:(\d+) sup:(\d+) Avg_Span:(\d+) sumProb:(.+)/;
			my $struct = $svtypes{$sv};
			if($os > $ie && $oe 
			if(!exists($storage{$c})){
				push(@{$storage{$c}}, [$is, $ie, $os, $oe, $struct, $sup, $prob]);
			}else{
				# Keep the smallest inside coordinates and only expand outside coords if at all
				my $found = 0;
				foreach my $row (@{$storage{$c}}){
					# outside edges overlap, inside edges overlap and sv matches
					if (($row->[2] < $oe && $row->[3] > $os) && ($row->[0] < $ie && $row->[1] > $is) && $struct eq $row->[4]){
						$found = 1;
						my ($istart, $iend) = refine_inside_edges($row->[0], $is, $row->[1], $ie);
						my ($estart, $eend) = expand_outside_edges($row->[2], $os, $row->[3], $oe);
						# Replace values
						$row->[0] = $istart; $row->[1] = $iend;
						$row->[2] = $estart; $row->[3] = $eend;
						$row->[5] += $sup; $row->[6] += $prob;
						last;
					} #outside edges overlap, inside edges overlap and sv does not match
					elsif(($row->[2] < $oe && $row->[3] > $os) && ($row->[0] < $ie && $row->[1] > $is) && $struct ne $row->[4]){
						$found = 1;
						my ($istart, $iend) = refine_inside_edges($row->[0], $is, $row->[1], $ie);
						my ($estart, $eend) = expand_outside_edges($row->[2], $os, $row->[3], $oe);
						# Replace values
						$row->[0] = $istart; $row->[1] = $iend;
						$row->[2] = $estart; $row->[3] = $eend;
						$row->[5] += $sup; $row->[6] += $prob;
						if (($row->[4] eq "ins" && $struct eq "inv") || ($row->[4] eq "inv" && $struct eq "ins")){
							$row->[4] = "insinv";
						}elsif(($row->[4] eq "del" && $struct eq "inv") || ($row->[4] eq "inv" && $struct eq "del")){
							$row->[4] = "delinv";
						}elsif(($row->[4] eq "ins" && $struct eq "del") || ($row->[4] eq "del" && $struct eq "ins")){
							# Hopefully this does not happen often!
							$row->[4] = "insdel";
						}elsif(($row->[4] eq "insinv" || $row->[4] eq "insdel" || $row->[4] eq "delinv") && !($row->[4] =~ m/$struct/)){
							# Again, another worst case scenario that we need to catch
							$row->[4] = "insdelinv";
						}elsif(($row->[4] eq "insinv" || $row->[4] eq "insdel" || $row->[4] eq "delinv") && ($row->[4] =~ m/$struct/)){
							# Everything is ok; no change
						}
						last;
					} # outside edges overlap, inside edges do not overlap
					elsif(($row->[2] < $oe && $row->[3] > $os) && !($row->[0] < $ie && $row->[1] > $is)){
						# Do nothing for now
					}
					
					
				}
				if(!($found)){
					push(@{$storage{$c}}, [$is, $ie, $os, $oe, $struct, $sup, $prob]);
				}
			}
		}
		close IN;
	}
	open (OUT, "> $out");
	foreach my $chr (sort {$a cmp $b} keys(%storage)){
		foreach my $rows (@{$storage{$chr}}){
			print OUT "$chr\t" . join("\t", @{$rows}) . "\n";
		}
	}
	close OUT;
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

sub single_generate_final {
	my ($file, $out) = @_;
	open(IN, "< $file") || die "could not open $file!\n";
	open(OUT, "> $out");
	my $here = <IN>;	# clearing the header junk from the file
	$here = <IN>;
	$here = <IN>;
	$here = <IN>;
	while (my $line = <IN>){
		chomp $line;
		my ($is, $ie, $os, $oe, $c, $sv, $sup, $avg, $prob) = $line =~ m/Inside_Start:(\d+) Inside_End:(\d+) OutSide_Start:(\d+) Oustide_End:(\d+) chro:(.+) SVtype:(\d+) sup:(\d+) Avg_Span:(\d+) sumProb:(.+)/;
		my $struct = $svtypes{$sv};
		print OUT "$c\t$ie\t$os\t$is\t$oe\t$struct\t$sup\t$prob\n";
	}
	close IN;
	close OUT;
}

sub most {
	my($a, $b) = @_;
	return ($a > $b)? $a : $b;
}
sub least {
	my ($a, $b) = @_;
	return ($a < $b) ? $a : $b;
}
