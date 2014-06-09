#!/usr/bin/perl
# This script takes input from the end of the lewis pipeline and generates spreadsheet files for use in publications
# The goal will be to intersect files from the various gene databases, snp coords and such
# Spreadsheet tab 1: CNVR designation with CNVR numbers based on numerical sorting of chromosomes with gene intersection data
# Spreadsheet tab 2: Individual CNV calls per animal with copynumber and gene intersection data (cross referenced with CNVRs)
# Spreadsheet tabs 3 - x: Gene copynumber estimations from database files input
#
# This uses the "tab" output from the merger program and uses the outside ends of the CNV call for analysis

use strict;
use Spreadsheet::WriteExcel;
use Getopt::Std;
use IO::Handle;
STDOUT->autoflush(1);

#my $bedtools = "/ibfs7/asg2/bickhartd/bin/bedtools";
my $bedtools = "";
my %opts;
my $usage = "$0 -m \<merged cnv calls file list\<tab files\>\> -c \<list of cn file\> -d \<file with gene db paths\> -o \<output file name\>\n";

getopt('mdco', \%opts);

unless(defined($opts{'m'}) && defined($opts{'d'}) && defined($opts{'o'})){
	print $usage;
	exit;
}

open(DB, "< $opts{d}") || die "could not open db file: $opts{d}\n";

# major genedb holder
my %file_hash; #{filebase}->{chr}->[row]->[start, end, gene]
while(my $l = <DB>){
	chomp $l;
	my @base = split(/\//, $l);
	my @bname = split(/\_/, $base[-1]);
	$file_hash{$bname[0]} = place_db_into_hash($l);
}
close DB;

# cn hash
my %cn_hash; #{animal base name}->{chr}->[row]->[start, end, cn#]
open(CN, "< $opts{c}") || die "could not open cn file: $opts{c}\n";
while(my $line = <CN>){
	my %thash;
	chomp $line;
	my @base = split(/\//, $line);
	my @bname = split(/[\_\.]/, $base[-1]);
	store_file_in_hash(\%thash, $line);
	$cn_hash{$bname[0]} = \%thash;
}
close CN;


my @dbfile_base = sort{$a cmp $b} keys(%file_hash);
my @header1 = ('cnvr#', 'calls', 'chr', 'start', 'end', 'cn#', '#animals', 'animal_list');

foreach my $db (@dbfile_base){
	push(@header1, "$db");
	push(@header1, "#$db");
}

my @header2 = ('chr', 'start', 'end','length','animal', 'cn#', 'cnvr#', 'calls');
foreach my $db(@dbfile_base){
	push(@header2, "$db");
	push(@header2, "#$db");
	push(@header2, "$db\_bp");
	push(@header2, "$db\%");
}

my @headerx = ('name', 'chr', 'start', 'end', 'length', '#cnvs','%cnv_ovlp', '#animals', 'animal_list');


# Worksheet line arrays
my @wks1data;
my @wks2data;
my %wksdbdata;

my @mergercnv_files;

# major individual cnv holder
my %indiv_cnvs; #{animal}->{chr}->[row]->[start, end, cn, callstr]
open(IN, "< $opts{m}") || die "Could not open merger files $opts{m}\n";
open(OUT, "> tempcat.bed") || die "Could not create temp bed merge output\n";
while(my $line = <IN>){
	chomp $line;
	my @anbase = split(/[\_\.]/, $line);
	push(@mergercnv_files, $line);
	open(AN, "< $line") || die "could not open individual merger file: $line\n";
	while(my $l = <AN>){
		chomp $l;
		my @segs = split(/\t/, $l);
		if($segs[4] - $segs[3] == 0){next;}
		push(@{$indiv_cnvs{$anbase[-3]}->{$segs[0]}}, [$segs[3], $segs[4], $segs[6], $segs[5]]);
		print OUT "$segs[0]\t$segs[3]\t$segs[4]\t$anbase[-3]\n";
	}
	close AN;
}
close IN;
close OUT;


# Now, working on sheet one data
print "Working on sheet one data\n";
# ('cnvr#', 'chr', 'start', 'end', 'cn#', '#animals', 'animal_list', 'dbstring','#db'
my %cnvrholder; # {chr}->[row]->[start, end, num]
#open(BED, "$bedtools merge -i tempcat.bed -nms |") || print "$!\n";
open(BED, "mergeBed -i tempcat.bed -nms |") || print "$!\n";
my $cnvrit = 1;
while (my $line = <BED>){
	chomp $line;
	my @segs = split(/\t/, $line);
	
	my @tans = split(/;/, $segs[3]);
	my @uans = retrieve_unique_ans(\@tans);
	my @callstrs;
	foreach my $ua (@uans){
		my $tstr = intersect_cnvcallstr($indiv_cnvs{$ua}->{$segs[0]}, $segs[1], $segs[2]);
		push(@callstrs, $tstr);
	}
	my @tmpwksline = ($cnvrit, join("/", @callstrs), $segs[0], $segs[1], $segs[2]);
	push(@{$cnvrholder{$segs[0]}}, [$segs[1], $segs[2], $cnvrit, join("/", @callstrs)]);
	my @cnavg;
	foreach my $t (@uans){
		my $tcn = determine_avg_cn($cn_hash{$t}->{$segs[0]}, $segs[1], $segs[2]); 
		push(@cnavg, $tcn);
	}
	push (@tmpwksline, average(\@cnavg));
	push(@tmpwksline, scalar(@tans));
	push(@tmpwksline, $segs[3]);
	foreach my $db (@dbfile_base){
		my @tempdb = intersect_gene_db($file_hash{$db}->{$segs[0]}, $segs[1], $segs[2]);
		push(@tmpwksline, ($tempdb[0], $tempdb[1]));
	}
	push(@wks1data, \@tmpwksline);
	$cnvrit++;
	print "$cnvrit..\r";
}
close BED;

# Now working on sheet two data
print "Working on sheet two data\n";
# ('chr', 'start', 'end','length','animal', 'cn#', 'cnvr#', 'db', '#db', 'db_bp', 'db%'
my %gene_cnv_hash; # {genename} = number individual hits
my %gene_avg_ovlp; # {genename}->[ovlp perc, ovlp perc...]
my %gene_animals; # {genename}->[an1, an2, ...]
my %gene_calls; # {genename}->[callstr1, callstr2, ...]
foreach my $canimals (sort{$a cmp $b} keys(%indiv_cnvs)){
	foreach my $chr (sort{my ($achrs) = $a =~ m/chr(.+)/; my ($bchrs) = $b =~ m/chr(.+)/; if($achrs =~ /X/){$achrs = 500;} if($bchrs =~ /X/){$bchrs = 500;} $achrs <=> $bchrs} keys(%{$indiv_cnvs{$canimals}})){
		foreach my $rref (@{$indiv_cnvs{$canimals}->{$chr}}){
			my @tmpwksline = ($chr, $rref->[0], $rref->[1], $rref->[1] - $rref->[0], $canimals, $rref->[2]);
			my ($cnvrnum, $callstr) = intersect_cnvr($cnvrholder{$chr}, $rref->[0], $rref->[1]);
			push(@tmpwksline, $cnvrnum, $callstr);
			foreach my $db (@dbfile_base){
				my @tempdb = intersect_gene_db($file_hash{$db}->{$chr}, $rref->[0], $rref->[1]);
				my @genes = split(/;/, $tempdb[0]);
				my @ovpercs = split(/;/, $tempdb[3]);
				for (my $x = 0; $x < scalar(@genes); $x++){
					$gene_cnv_hash{$genes[$x]} += 1;
					
					push(@{$gene_avg_ovlp{$genes[$x]}}, $ovpercs[$x]);
					push(@{$gene_animals{$genes[$x]}}, $canimals);
				}
				
				push(@tmpwksline, @tempdb);
			}
			push(@wks2data, \@tmpwksline);
		}
	}
}

# Now, doing the individual genedb entries
print "Now doing the individual genedb entries\n";
# 'name', 'chr', 'start', 'end', 'length', '#cnvs','%cnv_ovlp', '#animals', 'animal_list', animal cn1 ...);
push(@headerx, sort{$a cmp $b}(keys(%cn_hash)));
foreach my $db (@dbfile_base){
	foreach my $chr (sort{my ($achrs) = $a =~ m/chr(.+)/; my ($bchrs) = $b =~ m/chr(.+)/; if($achrs =~ /X/){$achrs = 500;} if($bchrs =~ /X/){$bchrs = 500;} $achrs <=> $bchrs} keys(%{$file_hash{$db}})){
		foreach my $row (@{$file_hash{$db}->{$chr}}){
			my @tmpwksline;
			my @ancnholder;
			foreach my $uanimals(sort{$a cmp $b}keys(%cn_hash)){
				my $tcn = determine_avg_cn($cn_hash{$uanimals}->{$chr}, $row->[0], $row->[1]);
				push(@ancnholder, $tcn);
			}
			if(exists($gene_cnv_hash{$row->[2]}) && exists($gene_avg_ovlp{$row->[2]})){
				my %unique;
				foreach my $v (@{$gene_animals{$row->[2]}}){ $unique{$v} = 1;}
				
				@tmpwksline = ($row->[2], $chr, $row->[0], $row->[1], $row->[1] - $row->[0], $gene_cnv_hash{$row->[2]}, average($gene_avg_ovlp{$row->[2]}), scalar(@{$gene_animals{$row->[2]}}), join(";", keys(%unique)), @ancnholder);
			}else{
				@tmpwksline = ($row->[2], $chr, $row->[0], $row->[1], $row->[1] - $row->[0], "", "", "", "", @ancnholder);
			}
			push(@{$wksdbdata{$db}}, \@tmpwksline);
		}
	}
}

print "Printing out to worksheet...\n";
my $workbook = Spreadsheet::WriteExcel->new($opts{'o'});
# Now for the CNVR sheet
my $cnvrsheet = $workbook->add_worksheet('CNVR listings');
my $headformat = $workbook->add_format();
$headformat->set_bold();
$headformat->set_bottom();
$cnvrsheet->freeze_panes(1,0); # freezes top header row
for (my $cols = 0; $cols < scalar(@header1); $cols++){
	$cnvrsheet->write(0, $cols, $header1[$cols], $headformat);
}
my $rown = 1;
foreach my $rows (@wks1data){
	for (my $cols = 0; $cols < scalar(@{$rows}); $cols++){
		$cnvrsheet->write($rown, $cols, $rows->[$cols]);
	}
	$rown++;
}

# Now for the indiv CNV sheet
my $indivsheet = $workbook->add_worksheet('Individual CNVs');
$indivsheet->freeze_panes(1,0);
for (my $cols = 0; $cols < scalar(@header2); $cols++){
	$indivsheet->write(0, $cols, $header2[$cols], $headformat);
}
my $rown = 1;
foreach my $rows (@wks2data){
	for (my $cols = 0; $cols < scalar(@{$rows}); $cols++){
		$indivsheet->write($rown, $cols, $rows->[$cols]);
	}
	$rown++;
}

# Finally, the individual geneDB sheets
foreach my $db (sort {$a cmp $b} keys(%wksdbdata)){
	my $nwks = $workbook->add_worksheet("$db\_data");
	$nwks->freeze_panes(1,0);
	for (my $cols = 0; $cols < scalar(@headerx); $cols++){
		$nwks->write(0, $cols, $headerx[$cols], $headformat);
	}
	my $rown = 1;
	foreach my $rows (@{$wksdbdata{$db}}){
		for (my $cols = 0; $cols < scalar(@{$rows}); $cols++){
			$nwks->write($rown, $cols, $rows->[$cols]);
		}
		$rown++;
	}
}
$workbook->close();
print "Finished. File is $opts{o}\n";
exit;
sub intersect_cnvcallstr{
	# CNVref is just "start, end, cn#, callstr at this point
	my ($cnvref, $start, $end) = @_;
	my @callstrs;
	foreach my $rows (@{$cnvref}){
		my $ovlp = overlap($rows->[0], $rows->[1], $start, $end);
		if($ovlp >= 1){
			push(@callstrs, $rows->[3]);
		}
	}
	return join(";", @callstrs);
}
sub retrieve_unique_ans{
	my ($aref) = @_;
	my %thash;
	foreach my $t (@{$aref}){
		$thash{$t} = 1;
	}
	return sort{$a cmp $b}(keys(%thash));
}
sub intersect_cnvr{
	# cnvrref is just "start, end, cnvr#, callstr" at this point
	my ($cnvrref, $start, $end) = @_;
	foreach my $rows (@{$cnvrref}){
		my $ovlp = overlap($rows->[0], $rows->[1], $start, $end);
		if($ovlp >= 1){
			return $rows->[2], $rows->[3];
		}
	}
	return 0, "";
}
sub intersect_gene_db{
	# dbregions is just "start, end and gene" at this point
	my ($dbregions, $start, $end) = @_;
	my @return; # ('genestring', #genes, bpoverlap, percent)
	my (@genes, @bpovlps, @percs);
	foreach my $rows (@{$dbregions}){
		my $ovlp = overlap($rows->[0], $rows->[1], $start, $end);
		if($ovlp >= 1){
			my $glen = $rows->[1] - $rows->[0];
			if($glen == 0){next;}
			my $perc = $ovlp / $glen;
			push(@genes, $rows->[2]);
			push(@bpovlps, $ovlp);
			push(@percs, $perc);
		}
	}
	$return[0] = join(";", @genes);
	$return[1] = scalar(@genes);
	$return[2] = join(";", @bpovlps);
	$return[3] = join(";", @percs);
	return @return;
}

sub place_db_into_hash{
	# dbfile must be in the bed style format: chr\tstart\tend\tname
	my ($dbfile) = @_;
	my %int_hash;
	open(IN, "< $dbfile") || die "Could not place DB listing into hash: $dbfile\n";
	while(my $line = <IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		push(@{$int_hash{$segs[0]}}, [$segs[1], $segs[2], $segs[3]]);
	}
	close IN;
	return \%int_hash;
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