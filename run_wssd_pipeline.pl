#!/usr/bin/perl
# This script takes input from mergeDocWindows.jar and processes it using Alkan's pipeline
# This repurposes Alkan's shell scripts into a more dynamic and less piecemeal wrapper, though not as streamlined as MrCaNNAVAR

use strict;
use Getopt::Std;

my $usage = "$0 -l \<file 1\> -s \<file 2\> -n \<file 3\>\n";
my %opts;
getopt('lsngb', \%opts);

# Essential files
my $bindir = "/ibfs7/asg2/bickhartd/bin";
my $gapfile = "/ibfs7/asg2/bickhartd/";
my $wssd_wgac = "";

my @f1_store; # Each index contains the reference to another array containing the data
my @f2_store;
my @f3_store;
my @f1c_store;
my @f3c_store;

# Calculate f1 and f3 average and stdev values

# Crop out the values that are above the stdev cutoff in f1 and f3 (remove wssd/wgac) and make control files

# Filter files to ensure they match line for line

# calculate control average and stdev values




exit;

sub filter_exact_bed {
	my ($a, $b, $o) = @_;
	open(AF, "< $a") || die "Could not open $a!\n";
	my %data;
	my $ls = 0;
	print "Putting $a into memory\n";
	while (my $line = <AF>){
		chomp $line;
		my @segs = split(/\t/, $line);
		$data{$segs[0]}->{$segs[1]} = $segs[2];
		$ls++;
	}
	close AF;

	open(OUT, "> $o");
	open(BF, "< $b") || die "Could not open $opts{b}!\n";
	my $rows = 0;
	print "Comparing $b to $a\n";
	while (my $l = <BF>){
		chomp $l;
		my @s = split(/\t/, $l);
		if ($data{$s[0]}->{$s[1]} == $s[2]){
			print OUT "$l\n";
			$rows++;
		} 
	}
	close BF;
	close OUT;
	my $err = abs($rows - $ls);
	print "The number of lines that did not match: $err\n";
}
sub twoPartOvp_mgsrt{
	my($input, $gapsfile) = @_;
	#my $left  = $opts{'L'};
	#my $right = $opts{'R'};
	my %fstlist = ();
	my %sndlist = ();
	#my $iCond   = $opts{'c'};
	#my $jCond   = $opts{'d'};
	my @output;

	    foreach my $row (@{$input})
	      {
		my @data = @{$row};
		push @{ $fstlist{$data[0]} }, [$data[1], $data[2]]
		  if(defined $fstlist{$data[0]} );
		
		$fstlist{$data[0]} = [ [$data[1], $data[2]] ]
		  if(!defined $fstlist{$data[0]} );
	      }

	    foreach my $row (@{$gapsfile})
	      {
		my @data = @{$row};
		push @{ $sndlist{$data[0]} }, [$data[1], $data[2]]
		  if(defined $sndlist{$data[0]} );
		
		$sndlist{$data[0]} = [ [$data[1], $data[2]] ]
		  if(!defined $sndlist{$data[0]} );
	      }


	my $union = join(' ', sort keys %fstlist);
	foreach my $key (sort keys %sndlist)
	  {
	    next if( $union =~ /\Q $key \E/ );
	    next if( $union =~ /^\Q$key\E/ );
	    next if( $union =~ /\Q$key\E$/ );
	    $union .= ' '.$key;
	  }
	my @allChrom = split(/ /, $union);

	foreach my $chrom (@allChrom)
	  {

	    @{$fstlist{$chrom}} = sort {$a->[0] <=> $b->[0] } @{$fstlist{$chrom}}
			if defined $fstlist{$chrom};
	    @{$sndlist{$chrom}} = sort {$a->[0] <=> $b->[0] } @{$sndlist{$chrom}}
			if defined $sndlist{$chrom};


	    # what is in i but not in j

		next if( !defined $fstlist{$chrom} );
		if( !defined $sndlist{$chrom} )
		  {
		    outputList( $fstlist{$chrom}, $chrom, *OUT);
		    next;
		  }
	      

	    # chrom is defined for both i and j
	    my $i = 0;
	    my $j = 0;
	    my $left = 1;
	    my $right = 0;
	    my ($leftStart, $rightStart) = (-1, -1);
	    while( $i < scalar( @{ $fstlist{$chrom} } ) && $j < scalar( @{ $sndlist{$chrom} } ))
	      {
		my $fstRoA = $fstlist{$chrom}->[$i];
		my $sndRoA = $sndlist{$chrom}->[$j];


		
		if( $fstRoA->[1] <= $sndRoA->[0] )
		  {
		    push(@output, [$chrom, $fstRoA->[0], $fstRoA->[1]]) if($left && $leftStart == -1 && $fstRoA->[1] > $fstRoA->[0] );
		    push(@output, [$chrom, $leftStart, $fstRoA->[1]]) if($left && $leftStart != -1 && $fstRoA->[1] > $leftStart);
		    $leftStart = -1;
		    $i++;
		    next;
		  }
		if( $sndRoA->[1] <= $fstRoA->[0] )
		  {
		    push(@output, [$chrom, $sndRoA->[0], $sndRoA->[1]]) if($right && $rightStart == -1 && $sndRoA->[1] > $sndRoA->[0]);
		    push(@output, [$chrom, $rightStart, $sndRoA->[1]]) if($right && $rightStart != -1 && $sndRoA->[1] > $rightStart);
		    $rightStart = -1;
		    $j++;
		    next;
		  }
		
		# pick the greatest start and least end
		my $start = ($fstRoA->[0] > $sndRoA->[0]) ? $fstRoA->[0] : $sndRoA->[0];
		my $end   = ($fstRoA->[1] > $sndRoA->[1]) ? $sndRoA->[1] : $fstRoA->[1];


		push(@output, [$chrom, $start, $end]) if( !$left && !$right && $end > $start );
		if( $left && $fstRoA->[0] < $start )
		  {
		    push(@output, [$chrom, $fstRoA->[0], $start]) if( $leftStart == -1 && $start > $fstRoA->[0] );
		    push(@output, [$chrom, $leftStart, $start]) if( $leftStart != -1 && $start > $leftStart );
		  }
		if( $right && $sndRoA->[0] < $start )
		  {
		    push(@output, [$chrom, $sndRoA->[0], $start]) if( $rightStart == -1 && $start > $sndRoA->[0] );
		    push(@output, [$chrom, $rightStart, $start]) if( $rightStart != -1 && $start > $rightStart  );
		  }
		

		# move the one whose end is smaller
		$i++ if( $fstRoA->[1] == $end );
		$j++ if( $sndRoA->[1] == $end );
		$leftStart  = ($fstRoA->[1] > $end) ? $end : -1;
		$rightStart = ($sndRoA->[1] > $end) ? $end : -1;
	      }

	    # append those that are left from -i table
	    while($left && $i < scalar( @{ $fstlist{$chrom} } ) )
	      {
		if( $leftStart != -1 )
		  {
		    push(@output, [$chrom, $leftStart, $fstlist{$chrom}->[$i]->[1]]) if( $fstlist{$chrom}->[$i]->[1] > $leftStart );
		    $leftStart = -1;
		  }
		else
		  {
		    push(@output, [$chrom, $fstlist{$chrom}->[$i]->[0], $fstlist{$chrom}->[$i]->[1]]);
		  }
		$i++;
	      }
	return @output;
	
}
sub coordsMerger_sort{
	my ($input) = @_;
	my @kcols = ("0");
	my @orderKey = ();   # to stored keys as their original order
	my %all_data = ();
	my @output;

	foreach my $row ( @{$input} )
	{
	    my @cols = @{$row};
	    my $key = '';
	  
		$key = $cols[0];

	    $cols[1] = $cols[1] -1 ;

	    my @addARRAY = ($cols[1], $cols[2]);

	    push @{$all_data{ $key } }, [@addARRAY]
	      if(defined $all_data{ $key } );

	    if(! defined $all_data{ $key } )
	      {
		$all_data{ $key } = [ [@addARRAY] ];
		push @orderKey, $key;
	      }
	  }


	foreach my $key (@orderKey)
	  {
	    my @sort_data = sort { $a->[0] <=> $b->[0] || $b->[1] <=> $a->[1] } @{$all_data{$key} };
	    my ($lastStart, $lastEnd, $lastV, $start, $end, $val);


	    foreach my $roA (@sort_data)
	      {
		($start, $end) = ($roA->[0], $roA->[1]);
		
		# 1st row
		if( !defined $lastStart )
		  {
		    ($lastStart, $lastEnd) = ($roA->[0], $roA->[1]);
		    next;
		  }


		# accual merge within the same name
		if( $start <= $lastEnd )
		  {
		    $lastEnd = $end if( $end > $lastEnd );
		    next;
		  }
		
		# $start > $lastEnd or name changed
		if(defined $lastStart)
		  {
			push(@output, [$key, $lastStart, $lastEnd]);
		   
		  }

		($lastStart, $lastEnd) = ($roA->[0], $roA->[1]);
	      }

	    if(defined $lastStart)
	      {
		push(@output, [$key, $lastStart, $lastEnd]);      
	      
	      }

	  }
	return @output;
}
sub wssd_picker {
	my ($f1, $f2, $ctgw, $highw, $newW, $oriW, $valCol) = @_;
	#my $ctgw  = $opts{'w'};
	#my $highw = $opts{'s'};
	#my $newW  = $opts{'i'};
	#my $oriW  = $opts{'n'};
	#my $valCol= $opts{'b'};
	
	my @output; 
	my %wins = ();
	foreach my $rows ( @{$f1})
	{

		my @data = @{$rows};
		next if( scalar(@data) < $valCol || $data[1] !~ /^\d+$/ );

		push @{ $wins{$data[0]} }, [$data[1], $data[2], $data[$valCol]]      if( defined  $wins{$data[0]} );
		$wins{$data[0]} =        [ [$data[1], $data[2], $data[$valCol]] ]    if( !defined $wins{$data[0]} );
		# the %wins hash is now a hash of arrays. It contains the start, end and value of interest
	}

	my %wssd_hash = ();
	foreach my $chrom (sort keys %wins)
	{

		my $winlist   = $wins{$chrom};	     # $winlist is a reference to the array of values
		my $wssdFlag  = 0;
		my $wssdChrom = '';
		my $wssdStart = 0;
		my $wssdLastStart = 0;                   # variable to store start of last win in WSSD region
		my $wssdEnd   = 0;                       # variable to store end of wssd region
		my $i         = 0;

		while( $i <= scalar( @$winlist ) - $ctgw )
		{

		# for each contiguous window wide
			my $highCt    = 0;    # how many wins in $ctgw with value > threshold
			my $start     = -1;   # to store the start pos of first win with value > threshold
			my $lastStart = -1;   # to store the start pos of last win with value > threshold
			my $end       = -1;   # to store the end pos of last win with value > threshold
	
			for(my $j = $i; $j < $i + $ctgw; $j++)
			{
				# ggcmp determines if the first argument and the third argument are less than the 
				# second argument. If they both are, then it returns 1. If not, it returns 0.
				# IF OPTION 'm' IS NOT CHOSEN:
				# gccmp will return a 1 if the first argument is greater than the second.
				# Otherwise, it will return a 0. 
				if( ggcmp($winlist->[$j]->[2], $opts{'c'}, $opts{'m'}) )
				{
					$highCt++;
					$start     = $winlist->[$j]->[0] if($start == -1);
					$lastStart = $winlist->[$j]->[0];
					$end       = $winlist->[$j]->[1];
				}
			}

			if( $highCt >= $highw )
			{
				if( !$wssdFlag )
				{
					($wssdChrom, $wssdStart) = ($chrom, $start);
					$wssdFlag = 1;
				}
				$wssdLastStart = $lastStart;
				$wssdEnd       = $end;
			}

			if( $wssdFlag && ( $i == scalar( @$winlist ) - $ctgw || $highCt < $highw ) )
			{
				push @{ $wssd_hash{$wssdChrom} }, [$wssdStart, $wssdEnd, $wssdLastStart] if( defined $wssd_hash{$wssdChrom}  );
				$wssd_hash{$wssdChrom} =      [ [$wssdStart, $wssdEnd, $wssdLastStart] ] if( !defined $wssd_hash{$wssdChrom} );
				$wssdFlag = 0;
			}

	
			$i++;
		}
	
		# take into consideration when total window # = highw and all of them are high enough
		if( scalar( @$winlist ) == $highw && $highw < $ctgw )
		{
			my $takeWSSD = 1;
			($wssdStart, $wssdEnd) = ( $winlist->[0]->[0],  $winlist->[scalar( @$winlist ) -1]->[1] );
			for(my $j = 0; $j < scalar( @$winlist ); $j++)
			{
				if( ! ggcmp($winlist->[$j]->[2], $opts{'c'}, $opts{'m'}) ) {  $takeWSSD = 0; last;	  }
			}
			if( $takeWSSD ) { $wssd_hash{$chrom} = [ [$wssdStart, $wssdEnd, $wssdStart] ]; }
		}
	}

	if(1)
	{
		%wins = ();
		foreach my $rows ( @{$f2})
		{
			chomp;
			my @data = @{$f2};
			next if( $data[1] !~ /^\d+$/ );

			push @{ $wins{$data[0]} }, [$data[1], $data[2], $data[$valCol]]   if(  defined $wins{$data[0]} );
			$wins{$data[0]} =        [ [$data[1], $data[2], $data[$valCol]] ] if( !defined $wins{$data[0]} );
		}
		

		foreach my $key (sort keys %wssd_hash)
		{
			my $refineWin = $wins{$key};
			my $wssdRoA   = $wssd_hash{$key};
			my $refineIdx = 0;
			my @wssdSort  = sort { $a->[0] <=> $b->[0] } @$wssdRoA;
			my @rWinSort  = sort { $a->[0] <=> $b->[0] } @$refineWin;
			my $wssd_prc  = 0;
			my @refWSSD   = ();

			foreach my $cdRoA (@wssdSort)
			{
				my $j = $refineIdx;
				my @ovalWSSD = ();             # WSSD no matter can be refined or not

				while( $j < scalar(@rWinSort) )
				{
		
					if( $rWinSort[$j]->[1] <= $cdRoA->[0] )
					{
						$j++;
						next;
					}
					last if( $rWinSort[$j]->[0] >= $cdRoA->[1] );
		
					if( !$wssd_prc && $rWinSort[$j]->[0] >= $cdRoA->[0] )
					{
						$refineIdx = $j;           # where to start in refined win for next WSSD record
						my $foundHigh = 0;

						while( $j < scalar(@rWinSort) && $j < $refineIdx + $oriW/$newW  )
						{
							if( $rWinSort[$j]->[2] > $opts{'t'} )
							{
								push @ovalWSSD, $rWinSort[$j]->[0];
								$foundHigh = 1;
								last;
							}
							$j++;
						}

						push @ovalWSSD, $cdRoA->[0] if( !$foundHigh );
						$wssd_prc = 1;
					}

					if( $rWinSort[$j]->[0] >= $cdRoA->[2] )
					{
						my $curIdx      = $j;
						my $lastHighEnd = -1;

						while( $j < scalar(@rWinSort) && $rWinSort[$j]->[1] <= $cdRoA->[1] )
						{
							if( ggcmp($rWinSort[$j]->[2], $opts{'t'}, $opts{'m'} ) )
							{
								$lastHighEnd = $rWinSort[$j]->[1];
							}
							$j++;
						}
						push @ovalWSSD, $lastHighEnd   if( $lastHighEnd > -1 );
						push @ovalWSSD, $cdRoA->[1]    if( $lastHighEnd == -1 );

						$wssd_prc = 0;
						last;    # jump out of refine window loop as cur. WSSD is finished
					}

					$j++;
				} # loop of refine window

				push @refWSSD, [@ovalWSSD] if ( scalar(@ovalWSSD) == 2);
			}# loop of raw WSSD
	
			@refWSSD = sort { $a->[0] <=> $b->[0] }  @refWSSD;
			foreach my $wssdRoA (@refWSSD)
			{
				push(@output, [$key, $wssdRoA->[0], $wssdRoA->[1]]);
				
			}
		}
    
	}
	return @output;

}
sub ggcmp{
    my ($val, $cut, $rev) = @_;

    return 1 if ($rev  && $val < $cut );
    return 1 if (!$rev && $val > $cut );

    return 0;
}
sub part_gc_depth {
	my ($input, $output) = @_;
	my @totdepth;
	my @count;
	my $gc;
	my $depth;
	my $i;
	my $j;
	for ($i=0;$i<1000;$i++) { $totdepth[$i]=0.0; $count[$i]=0; }
	
	open (OUT, "> $output") || die "Could not create $output!\n";
	
	foreach my $rows (@{$input}){
		my @seg = @{$rows};
		chomp $seg[4];
		$gc = $seg[3];
		$depth = $seg[4];
		$i = int($gc * 1000);
		$count[$i]++;
		$totdepth[$i]+= $depth;
	}

	for ($i=0;$i<1000;$i++){
		if ($count[$i]==0 || $totdepth[$i] == 0){
			$j=$i-1;
		while ($j>=0 && $totdepth[$j]==0){
			$j--;
		}
		if ($totdepth[$j]==0 || $count[$j]==0){
			$j=$i+1;
			while ($j<1000 && $totdepth[$j]==0){ 
				$j++;
			}
		}

		$totdepth[$i]=$totdepth[$j];
		$count[$i]=$count[$j];
		}
		if ($totdepth[$i] != 0.0){
			my $value = $totdepth[$i] / $count[$i];
			print OUT "$i\t$value\n";
		}
		else{
			print OUT "$i\tEMPTY\n";
		}  
	}
}
sub depth_loess_avg {
	my ($input, $avgname, $expect, $divide) = @_;
	$expect = $expect / $divide;
	open(AVG, "< $avgname") || die "Could not open $avgname!\n";
	
	my @return;
	my @avg;
	my ($i, $depth);
	while(<AVG>){
		my @seg = split(/\t/);
		$i = $seg[0];
		$depth = $seg[1];
		chomp $depth;
		$avg[$i]=$depth / $divide;
	}
	close AVG;  
	foreach my $rows (@{$input}){
		my @segments = @{$rows};
		my $chr = $segments[0];
		my $start = $segments[1];
		my $end = $segments[2];
		my $gc = $segments[3];
		my $depth = $segments[4];
		chomp $depth;
		$i = $gc * 1000;
		my $fx = $avg[$i]; 
		my $newcn = $depth - ($fx - $expect);
		if ($newcn < 0){ 
			$newcn = 0.0;
		}
		push(@return, ($chr, $start, $end, $gc, $newcn));
	}
	return @return;
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