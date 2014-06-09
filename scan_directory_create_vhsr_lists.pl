#!/usr/bin/perl
# This script is designed to fill the .list files for the pipeline 
#  BIBR07.280.sr.concordant.vh  BIBR07.340.sr.split.fq       BIBR07.401.sr.single.txt     BIBR07.98.sr.mrsfast.bam
#BIBR07.158.sr.single.txt     BIBR07.219.sr.mrsfast.bam    BIBR07.280.sr.divet.vh       BIBR07.341.sr.concordant.vh  BIBR07.401.sr.split.fq       BIBR07.98.sr.single.txt
#BIBR07.158.sr.split.fq       BIBR07.219.sr.single.txt     BIBR07.280.sr.mrsfast.bam    BIBR07.341.sr.divet.vh       BIBR07.402.sr.concordant.vh  BIBR07.98.sr.split.fq
#BIBR07.159.sr.concordant.vh  BIBR07.219.sr.split.fq       BIBR07.280.sr.single.txt 

use strict;
use Class::Struct;
use Getopt::Std;

struct(list => {
	'splitbam' => '$',
	'divet' => '$',
	'anchor' => '$',
	'insert' => '$',
	'stdev' => '$',
});

struct(conf => {
	'numlist' => '@',
	'insert' => '$',
	'stdev' => '$',
});

my $usage = "$0 -i <input config file> -d <animal name> -c <checkpoint file> -o <output list file name> -j <original dir>\n";
my %opts;
getopt('idcoj', \%opts);

my @orderedconf; 
open(CNF, "< $opts{i}") || die "Could not open config file!\n$usage";
my $isfiles = 0;
while(my $line = <CNF>){
	if($isfiles){
		chomp $line;
		my @segs = split(/\t/, $line);
		if($segs[5] eq $opts{'d'}){
			push(@orderedconf, conf->new('insert' => $segs[2], 'stdev' => $segs[3], 'numlist' => []));
		}
	}elsif($line =~ /\[files\]/){
		$isfiles = 1;
	}elsif($line =~ m;\[/files\];){
		$isfiles = 0;
	}
}
close CNF;

open(CHK, "< $opts{c}") || die "Could not open checkpoint file!\n$usage";
my $iterator = 0;
while(my $line = <CHK>){
	chomp $line;
	my @segs = split(/\t/, $line);
	if($segs[0] eq $opts{'d'} && $segs[1] eq "pfq1"){
		shift(@segs);
		shift(@segs);
		shift(@segs);
		foreach my $fq (@segs){
			my @fsegs = split(/[_\.]/, $fq);
			push(@{$orderedconf[$iterator]->numlist()}, $fsegs[-2]);
		}
		$iterator++;
		if($iterator >= scalar(@orderedconf)){
			print "Something wrong: $iterator " . scalar(@orderedconf) . "\n";
		}
	}
}
close CHK;

opendir(DIR, "$opts{d}/split_read") || die "Could not open directory: $opts{d}/split_read!\n$usage";
my %initialholder;
my %listholder; # {num} -> list
while(my $f = readdir(DIR)){
	chomp $f;
	my @fsegs = split(/\./, $f);
	push(@{$initialholder{$fsegs[1]}}, $f);
}
closedir(DIR);

foreach my $conflist (@orderedconf){
	my $nums = $conflist->numlist();
	foreach my $n (@{$nums}){
		my $files = $initialholder{$n};
		my $l = list->new();
		foreach my $f (@{$files}){
			if($f =~ /divet.vh/){
				$l->divet($f);
			}elsif($f =~ /mrsfast.bam/){
				$l->splitbam($f);
			}elsif($f =~ /single.txt/){
				$l->anchor($f);
			}
		}
		$l->insert($conflist->insert());
		$l->stdev($conflist->stdev());
		$listholder{$n} = $l;
	}
}

open(OUT, "> $opts{o}") || die "Could not open $opts{o}!\n$usage";
foreach my $nums (sort {$a <=> $b} keys(%listholder)){
	my $row = $listholder{$nums};
	my $splitbam = $row->splitbam();
	my $divet = $row->divet();
	my $anchor = $row->anchor();
	my $insert = $row->insert();
	my $stdev = $row->stdev();
	print OUT "$opts{j}/$opts{d}/split_read/$splitbam\t$opts{j}/$opts{d}/split_read/$divet\t$opts{j}/$opts{d}/split_read/$anchor\t$insert\t$stdev\n"; 
}
close OUT;
exit;