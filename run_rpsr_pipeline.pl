#!/usr/bin/perl
# This is the pipeline that runs the read pair and split read caller
# It is just a glorified perl wrapper for my Java program, but it separates out the chromosomes into different threads

use Forks::Super;
use strict;
use Getopt::Std;
use Class::Struct;

my %opts;
my $usage = "$0 -s \<input variant map file\> -o \<outbase name\> -c \<basename of animal\> -l \<boolean: is threaded?\>
\t -g \<gap file\> -j \<java path\> -b \<bin dir\> -f \<chr lens file\> -m \<command dir\>
\t -p \<max proc\>";
getopt('soclgjbfmp', \%opts);

struct(events => {
	'insertion' => '$',
	'deletion' => '$',
	'tandem' => '$',
	'inversion' => '$',
});

unless(scalar(@{keys(%opts)}) == 10){
	print "All arguments are required!\n";
	print $usage;
	exit;
}

open(IN, "< $opts{f}") || die "$0: could not open chr lengths file - $opts{f}!\n";
my %chrs; # chr lengths container
while(my $line = <IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	$chrs{$segs[0]} = $segs[1];
}
close IN;

my %eventfiles; #{chr} = events
foreach my $c (sort {$a cmp $b} keys(%chrs)){
	my $command;
	if($opts{'l'}){
		$command = "$opts{j}/java -Xmx4g -jar setWeightCoverVHSRDiscovery.jar -s $opts{s} -o $opts{o}.$c -g $opts{g} -c $c";
	}else{
		$command = "bsub -J V_$opts{c}_$c\_RPSR -K ";
		$command .= "-o $opts{m}/out/V_$opts{c}_$c\_RPSR.out ";
		$command .= "-e $opts{m}/err/V_$opts{c}_$c\_RPSR.err ";
		$command .= qq{-R "rusage[mem=4000]" };
		$command .= "-q normal ";
		$command .= "\'$opts{j}/java -Xmx4g -jar setWeightCoverVHSRDiscovery.jar -s $opts{s} -o $opts{o}.$c -g $opts{g} -c $c\'";
	}
	
	open(OUT, "> $opts{m}/V_$opts{c}_$c\_RPSR.sh") || die "$0: could not create command outfile for chr: $c!\n";
	print OUT $command;
	close OUT;
	
	$eventfiles{$c} = events->new(
		'insertion' => "$opts{o}.$c.vhsr.insertions",
		'inversion' => "$opts{o}.$c.vhsr.inversions",
		'tandem' => "$opts{o}.$c.vhsr.tand",
		'deletion' => "$opts{o}.$c.vhsr.deletions"
	);
	
	fork { timeout => 1800, cmd => $command, max_proc => $opts{p} };
}

my $pidwait = waitall();
if($pidwait > 0){
	# Catch any remaining jobs
	print "$0: Had to wait for $pidwait child processes\n";
}

# Now to merge all of the entries into one file
open(INS, "> $opts{o}.cat.vhsr.insertions");
open(INV, "> $opts{o}.cat.vhsr.inversions");
open(TAN, "> $opts{o}.cat.vhsr.tand");
open(DEL, "> $opts{o}.cat.vhsr.deletions");
open(ALL, "> $opts{o}.cat.vhsr.all");
open(FILT, "> $opts{o}.cat.vhsr.allfilt");
foreach my $c ( map {$_->[0]} 
		sort {$a->[1] cmp $b->[1]} 
		map{[$_, returnChrNum($_)]} keys(%eventfiles)){
	my $e = $eventfiles{$c};
	my $ins = $e->insertion();
	my $inv = $e->inversion();
	my $tan = $e->tandem();
	my $del = $e->deletion();
	
	open(IN, "< $ins");
	while(my $line = <IN>){
		print INS $line;
		print ALL $line;
		if(donotfilter($line)){
			print FILT $line;
		}
	}
	close IN;
	
	open(IN, "< $inv");
	while(my $line = <IN>){
		print INV $line;
		print ALL $line;
		if(donotfilter($line)){
			print FILT $line;
		}
	}
	close IN;
	
	open(IN, "< $tan");
	while(my $line = <IN>){
		print TAN $line;
		print ALL $line;
		if(donotfilter($line)){
			print FILT $line;
		}
	}
	close IN;
	
	open(IN, "< $del");
	while(my $line = <IN>){
		print DEL $line;
		print ALL $line;
		if(donotfilter($line)){
			print FILT $line;
		}
	}
	close IN;
	
	system("rm $ins $inv $tan $del");
}
close DEL;
close TAN;
close ALL;
close INS;
close INV;
close FILT;
exit;
sub donotfilter{
	my($line) = @_;
	my @segs = split(/\t/, $line);
	chomp $segs[9];
	my $sum = $segs[6] + $segs[7] + $segs[8];
	if($sum == 1 || $segs[9] == 1){
		return 0;
	}else{
		return 1;
	}
}
sub returnChrNum{
	my ($chr) = @_;
	if($chr =~ m/chr(\d+)/){
		return $1;
	}elsif($chr =~ m/chr([XYZWxyzw])/){
		if($1 =~ m/[XxZz]/){
			return 500;
		}elsif($1 =~ m/[WwYy]/){
			return 501;
		}else{
			return 502;
		}
	}else{
		return 502;
	}	
}