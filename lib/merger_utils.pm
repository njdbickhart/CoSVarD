#!/usr/bin/perl
# This is a collection of the subroutines used in the final assembly sorting and merger pipeline scripts

use strict;

package merger_utils;

sub new {
	my ($class,%arg) = @_;
	my $self={};
	bless($self,$class || ref($class));
	return $self;
}

sub exonerate_intersect {
	my ($self, $exfile, $confile, $chr, $start, $end, $svstr) = @_;
}

sub filterMatches {
	my ($self,$matches) = @_;
	my %coords=();
	for (@$matches) {
		next if $_ !~ /^\S+\s+\d+\s+\d+\s+\d+\s+/;
		my $line = $_;
		if ($line =~ /INV/ || $line eq 'DEL' || $line eq 'DELINS') {
			$coords{$1}{"$2-$3"}=$line if $line =~ /^(\S+)\s+(\d+)\s+(\d+)\s+\d+\s+/;
		}
	}
	my @pass;
	# this part filters out false INS from parseNeighbor
	if (keys %coords) {
		foreach my $match (@$matches) {
			my ($c,$s,$e) = ($1,$2,$3) if $match =~ /^(\S+)\s+(\d+)\s+(\d+)\s+\d+\s+/;
			next if !$s || !$e;
			if ($coords{$c}{"$s-$e"} && $coords{$c}{"$s-$e"} eq $match ) {
				push @pass, $match;
				next;
			}
			elsif ($match =~ /internal/) {
				push @pass, $match;
				next;
			}
			my $over=0;
			foreach my $span (keys %{$coords{$c}}) {
				my ($a,$b) = split "-", $span;
				foreach ($a,$b) {
					my $span1 = ($_-50)."-".($_+50);
					my $int = $self->spanCheck($span1,"$s-$s");
					my $int2 = $self->spanCheck($span1,"$e-$e");
					$over++ if $int>0 || $int2>0;
				}
			}
			if (!$over) {
				push @pass, $match;
			}
		}
		return \@pass;	
	}
	else { return $matches }
}

sub parseIndels {
	my ($self,$hits) = @_;
	my @res;
	for (@$hits) {
		next if /no hit/ || /^\s*$/;
		my @cols = $self->exonerateHitParse($_);
		($cols[5],$cols[6]) = ($cols[6],$cols[5]) if $cols[5]>$cols[6];
		$cols[11] = $self->reverseCIGAR if $cols[7] eq '-';
		my $gcoord=$cols[5];
		my %indels=();
		my ($merged,$change) = $self->mergeCIGAR($cols[11]);
		while ($merged =~ m/(\S+)\s+(\d+)/g) {
			my $type = 'internal';
			$type .= "m" if $change;
			if ($1 eq 'M') {
				$gcoord+=$2;
			}
			elsif ($1 eq 'I') {
				push @{$indels{$cols[0]}}, [$gcoord-1,$gcoord+1,$2,"INS$type"];
			}
			elsif ($1 eq 'D') {
				push @{$indels{$cols[0]}}, [$gcoord,$gcoord+$2,$2,"DEL$type"];
				$gcoord+=$2;
			}
		}
		for (keys %indels) {
			my $contig = $_;
			for (0..$#{$indels{$contig}}) {
				my $print = "$cols[4] ";
				$print .= join(" ", @{$indels{$contig}[$_]});
				$print .= " $contig";
				push @res, $print;
			}
		}
	}
	return \@res;

}

sub classifyBestHits {
	my ($self,$hits,$coords,$classes) = @_;
	my $topHitInRange='';
	my $exact=0;
	my @ranges;
	my @dir;
	my %seen;
	my $contiglen=0;
	my $topdir;
	for (@$hits) {
		my @col = split /\s+/, $_;
		next if $seen{"$col[1]-$col[2]"};
		my $set = $self->getSet($col[1],$col[2],$hits); # get sets with same qcoords (default q is always +; t is + or -)
		my @best = $self->getBestHits($set,$coords); # find best hits 
		
		if (!$topHitInRange && $best[0]) {
			$topHitInRange = "$col[1]-$col[2]";
			$exact = 1 if $col[10] == $col[2] && $col[1] == 0;
			$topdir = $col[7];
			if ($exact==1) {
				push @{$classes->{full}}, $best[0];
				$best[0] =~ s/\n/ FULL/;
			}
			else {
				push @{$classes->{part}}, $best[0];
				$best[0] =~ s/\n/ PART/;
				$contiglen=$col[10];
			}
		}
		push @{$seen{"$col[1]-$col[2]"}}, @best;
		push @ranges, "$col[1]-$col[2]" if $best[0]; # keep best scoring hits WITHIN EXPECTED RANGE in order
		push @dir, $col[7] if $best[0];
	}
	return($topHitInRange,$topdir,\@ranges,\@dir,\%seen,$exact,$contiglen,$classes);
}

sub parseContigs {
	my ($self,$file,$type) = @_;
	$type = lc $type;
	my $contigList = ();
	my $command = '';
	if ($type eq 'list') {
		$command = "cat $file";
	}
	elsif ($type eq 'fasta') {
		$command = qq(awk '/^>/' $file | cut -f 1 | sed 's/^>//' | sort -u);
	}
	if (!$command) {
		print STDERR "Can't parse format $type\n";
		return '0';
	}
	foreach my $item (`$command`) {
		chomp $item;
		$contigList->{$item} = 1;
	}
	return $contigList;
}

# get list of hits having the same query region

sub getSet {
	my ($self,$start,$end,$hits) = @_;
	my @set;
	foreach my $hit (@$hits) {
		my @col = split /\s+/, $hit;
		if ( $col[1]==$start && $col[2]==$end ) {
			push @set, $hit;
		}
	}
	return \@set;
}

# Input a set of hits and a set of coordinates to match to
# and get the best hit that overlaps the region

sub getBestHits {
	my ($self,$set,$coords) = @_;
	my $best = '';
	my $bestRegionHit = '';
	my $bestScore = 0;
	foreach my $hit (@$set) {
		my @col = split /\s+/, $hit;
		my @coords2 = @col[4..6];
		($coords2[1],$coords2[2]) = ($coords2[2],$coords2[1]) if $coords2[1]>$coords2[2];
		if (!$bestScore) { # top scoring hit
			$best = $hit;
			$bestScore = $col[8];
			next unless $coords->[0] eq $coords2[0]; # check Chr first
			my $coordCheck = $self->coordCheck($coords,\@coords2);
			if ($coordCheck) {
				$bestRegionHit = $hit;
				$best = '';
			}
		}
		# Look for best hit(s) to Expected region;
		# if top hit is not within expected region, return 
		# results as well
		elsif ( ($bestRegionHit && $coords->[0] eq $coords2[0] && $bestScore == $col[6]) ||
			(!$bestRegionHit && $coords->[0] eq $coords2[0]) )
		{
			my $coordCheck = $self->coordCheck($coords,\@coords2);
			if ($coordCheck) {
				$bestScore = $col[6] if !$bestRegionHit;
				$bestRegionHit .= "$hit";
			}
		}
	}
	return ($bestRegionHit,$best);
}

sub checkOverlap {
	my ($self,$obj,$range,$dir,$topdir,$rangedir) = @_;
	my $count=0;
	my @range = @$range;
	foreach my $i (0..$#range) {
		$count++;
		my $r = $range[$i];
		my ($min,$max) = split "-", $r;
		my $size = (max $obj - min $obj) < ($max-$min) ? (max $obj - min $obj) : ($max-$min);
		# overlap in opposite strand; allow more overlap
		if ($rangedir->[$i] ne $topdir) {
			if ($dir eq 'rev' && $max < (max $obj) && $max>(min $obj) && $min < (min $obj-20)) {
				my $diff = abs($max-(min $obj));
				my $ovlp = $max - min $obj;
				# ignore if overlap is more than 95% of the smallest hit
				return ($r,$count,$diff) if $ovlp/$size < 0.95; # exiting after the first match chooses the highest scoring hit
			}
			if ($dir eq 'for' && $min > (min $obj) && $min<(max $obj) && $max > 20+(max $obj)) {
				my $diff = abs((max $obj)-$min);
				my $ovlp = max $obj - $min;
				# ignore if overlap is more than 75%
				return ($r,$count,$diff) if $ovlp/$size < 0.95; # exiting after the first match chooses the highest scoring hit 
			}

		}
		else {	
			if ($dir eq 'rev' && $max <= 500+(min $obj) && $max>(min $obj) && $min < (min $obj)) {
				my $diff = abs($max-(min $obj));
				my $ovlp = $max - min $obj;
				# ignore if overlap is more than 75% of the smallest hit
				return ($r,$count,$diff) if $ovlp/$size < 0.75; # exiting after the first match chooses the highest scoring hit within 50bp
			}
			if ($dir eq 'for' && $min >= (max $obj)-500 && $min<(max $obj) &&  $max > (max $obj)) {
				my $diff = abs((max $obj)-$min);
				my $ovlp = max $obj - $min;
				# ignore if overlap is more than 75%
				return ($r,$count,$diff) if $ovlp/$size < 0.75; # exiting after the first match chooses the highest scoring hit within 50bp
			}
		}
	}
	return '';
}

sub checkBlunt {
	my ($self,$obj,$range,$dir) = @_;
	foreach my $r (@$range) {
		my ($min,$max) = split "-", $r;
		if ($dir eq 'rev' && $max == min $obj) {
			return $r; # exiting after the first match chooses the highest scoring hit
		}
		elsif ($dir eq 'for' && $min == max $obj) {
			return $r;
		}
	}
	return '0';
}

sub checkClosest {
	my ($self,$obj,$range,$dir) = @_;
	my $count=0;
	my @lowest=();
	my $diff='';
	foreach my $r (@$range) {
		$count++;
		my ($min,$max) = split "-", $r;
		if ($dir eq 'rev') {
			if ($max < min $obj) {
				if (!@lowest || (@lowest && $max > $lowest[1])) {
					@lowest = ($min,$max);
					$diff = abs((min $obj)-$max);
				}
			}
		}
		if ($dir eq 'for') {
			if ($min > max $obj) {
				if (!@lowest || (@lowest && $min < $lowest[0])) {
					@lowest = ($min,$max);
					$diff = abs($min+1-(max $obj));
				}
			}
		}
	}
	return ("$lowest[0]-$lowest[1]",$count,$diff) if @lowest;
	return '0' if !@lowest;

}

# check chrom and span
sub coordCheck {
	my ($self,$set1,$set2) = @_;
	# Check chromosome first
	if ($set1->[0] ne $set2->[0]) {
		return '0';
	}
	my $span1;
	# Check for intersect of coordinates
	#my $span1 = new Set::IntSpan "$set1->[1]-$set1->[2]";
	my $span2 = "$set2->[1]-$set2->[2]";
	if (intersect $span1 $span2) {
		return '1';
	}
	else {
		return '0';
	}
}

# change coords back to relative to chr start if using a subsequence of a chr

sub coordfix {
	my ($self,$offset,@hits) = @_;
	my @res = ();
	for (@hits) {
		my @cols = split /\s+/,$_;
		my $line = join " ", (@cols[0..4],($cols[5]+$offset),($cols[6]+$offset),@cols[7..$#cols]);
		push @res, "$line\n";
	}
	return @res;
}

# check span only
sub spanCheck {
	my ($self,$span1,$span2) = @_;
	my @span1 = split "-", $span1;
	my @span2 = split "-", $span2;
	my $int;
	#my $set = new Set::IntSpan "$span1[0]-$span1[1]";
	#my $int = intersect $set "$span2[0]-$span2[1]";
	my $intsize = $int=~/(\d+)-(\d+)/ ? $2-$1 : '0';
	if (!$intsize) {
		return '0';
	}
	else {
		return $intsize;
	}
}

sub exonerateHitParse {
	my ($self,$line) = @_;
	my @cols = split /\s+/, $line;
	my $rest = join " ", @cols[11..$#cols];
	return (@cols[0..10],$rest);
}

sub reverseCIGAR {
	my ($self,$line) = @_;
	my  @rev = '';
	while (m/(\S+\s+\d+)/g) {
		push @rev, $1;
	}
	@rev = reverse @rev if @rev;
	my $rev = join " ", @rev;
	return $rev;
}

sub mergeCIGAR {
	my ($self,$line) = @_;
	my $oldline = $line;
	foreach my $type ('D','I') {
		while ($line =~ m/($type\s+(\d+)\s+M\s+(\d+)\s+$type\s+(\d+))/g) {
			my $match = $1;
			my $pos = pos($line);
			pos($line) = $pos-(length$match)+1;
			if ($3<50) {
				#my $size = $2+$4;
				#my $msize = $3;
				my $size = $2+$3+$4; # don't count spurious matches in gappy regions
				#$line =~ s/$match/$type $size M $msize/;
				$line =~ s/$match/$type $size/;
				$line = $self->mergeMatches($line);
			}
		}
	}
	if ($oldline ne $line) {
		return ($line,1);
	}
	else {
		return $line;
	}
}

sub mergeMatches {
	my ($self,$line) = @_;
	while ($line =~ m/(M\s+(\d+)\s+M\s+(\d+))/g) {
		my $match = $1;
		my $pos = pos($line);
		pos($line) = $pos-(length$match)+1;
		if ($2<20||$3<20) {
			my $size = $2+$3;
			$line =~ s/$match/M $size/;
		}
	}
	return $line;
}

# adjust breakpoints for split alignments, if there
# are repeats flanking. Keep the 5' coord max, and shift the 3'coord
#
sub breakpointCheck {
	my ($self,$overlap,$c1,$c2,$first) = @_;
		if ($c1>$c2) {
			$c1+=$overlap;
			return ($c2,$c1);
		}
		elsif ($c1<$c2) {
			$c2+=$overlap;
			return ($c1,$c2);
		}
		elsif ($c1==$c2) {
			$c2+=$overlap if $first==1;
			$c1+=$overlap if $first==2;
			return ($c1,$c2);
		}
}
#Input: $line_from_exonerate
#Output: @array_of_cols
sub exonerateHitParse {
	my ($self,$line) = @_;
	if($line =~ /Query:/){next;}
	my @cols = split /\s+/, $line;
	my $rest = join " ", @cols[11..$#cols];
	return (@cols[0..10],$rest);
}
# Input: $cnregions_hash_ref, $start, $end
# Output: average of CN from region bounds
sub determine_avg_cn{
	my ($self, $cnregions, $start, $end) = @_;
	my @cns;
	foreach my $reg (@{$cnregions}){
		if($self->overlap($reg->[0], $reg->[1], $start, $end) > 0){
			push(@cns, $reg->[2]);
		}
	}
	if(scalar(@cns) == 0){return -1;}
	return $self->average(\@cns);
}
sub average {
	my ($self, $aref) = @_;
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
	my ($self, $a, $b) = @_;
	return ($a > $b)? $b : $a;
}
sub most {
	my ($self, $a, $b) = @_;
	return ($a > $b)? $a : $b;
}
# Shamelessly copied from bedtools
# Returns a number with negative numbers indicating no overlap and positive numbers indicating overlap
sub overlap {
	my ($self, $s1, $e1, $s2, $e2) = @_;
	return ($self->least($e1, $e2) - $self->most($s1, $s2));
}
sub return_longer_string {
	my ($self, $st1, $st2) = @_;
	if(length($st1) >= length($st2)){
		return $st1;
	}else{
		return $st2;
	}
}
sub expand_outside_edges{
	my($self, $s1, $s2, $e1, $e2) = @_;
	my $start = ($s1 <= $s2) ? $s1 : $s2;
	my $end = ($e1 >= $e2) ? $e1 : $e2;
	return $start, $end;
	
}
sub refine_inside_edges{
	my ($self, $s1, $s2, $e1, $e2) = @_;
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
# Input: $hashref_for_storage, $file_string
# Output: first argument as stored hash (chr as key)
sub store_file_in_hash{
	my ($self, $href, $file) = @_;
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
1;