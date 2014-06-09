#!/usr/bin/perl
# This is a one-shot script designed to take a list of Steve's fastq files and order them for a spreadsheet for the pipeline

use strict;
use Class::Struct;
use Forks::Super;

struct( read_sort => {
	'fq1' => '$',
	'fq2' => '$',
	'insert' => '$',
	'sd' => '$',
	'fnum' => '$',
	'fc' => '$',
	'segs' => '@',
	'sam' => '$',
});

my %key_conv = (
'GL68' => 'BTHO42',
'GL69' => 'BTHO44',
'GL70' => 'BTHO46',
'GL71' => 'BTHO47',
'GL72' => 'BTHO48',
'GL73' => 'BTHO49',
'GL74' => 'BTHO50',
'GL75' => 'BTHO51',
'GL76' => 'BTHO52',
'GL77' => 'BTHO53',
'GL78' => 'BTHO54',
'GL79' => 'BTHO56',
'GL80' => 'BIBR05',
'GL81' => 'BIBR07',
'GL82' => 'BTHO10',
'GL83' => 'BTHO57',
'GL84' => 'BTJE10',
'GL85' => 'BTLM04',
'GL86' => 'BTHO19',
'GL87' => 'BTHO20',
'GL88' => 'BTHO22',
'GL89' => 'BTHO25',
'GL90' => 'BTHO29',
'GL91' => 'BTHO32',
);

my %data_struct; # {an_name}->[]->read_sort
my $tmpdir = "/mnt/iscsi/vnx_gliu_7/tmp";

chomp(@ARGV);
open(IN, "< $ARGV[0]") || die "Could not open $ARGV[0]!\n";
while(my $line = <IN>){
	chomp $line;
	my @s = split(/\//, $line);
	my @r = split(/[_\.]/, $s[-1]);
	my $an;
	if(exists($key_conv{$r[0]})){
		print "Converting $r[0] to: " . $key_conv{$r[0]} . "\n";
		$an = $key_conv{$r[0]};
	}else {
		$an = $r[0];
	}
	my $found = 0;
	if($r[3] eq "R1"){
		if(exists $data_struct{$an}){
			foreach my $row (@{$data_struct{$an}}){
				if($row->fnum() eq $r[4] && $row->fc() eq $s[5]){
					$row->fq1($line);
					$found = 1;
				}
			}
		}
		unless($found){
			push(@{$data_struct{$an}}, read_sort->new(
				'fq1' => $line,
				'fnum' => $r[4],
				'segs' => \@r,
				'fc' => $s[5],
			));
		}
	}else{
		if(exists $data_struct{$an}){
			foreach my $row (@{$data_struct{$an}}){
				if($row->fnum() eq $r[4] && $row->fc() eq $s[5]){
					$row->fq2($line);
					$found = 1;
				}
			}
		}
		unless($found){
			push(@{$data_struct{$an}}, read_sort->new(
				'fq2' => $line,
				'fnum' => $r[4],
				'segs' => \@r,
				'fc' => $s[5],
			));

		}
	}
}
close IN;

open(OUT, "> temp_test_starter.txt");
foreach my $an (keys(%data_struct)){
	foreach my $row(@{$data_struct{$an}}){
		print OUT $row->fq1() . "\t" . $row->fq2 . "\t" . "\t" .  "\tp\t$an\n";
	}
}
close OUT;

print "Done loading files into spreadsheet\n";
print "Now testing insert sizes\n";
mkdir("$tmpdir");

my $counter = 0;
foreach my $an (keys(%data_struct)){
	print "Working on $an\n";
	
	foreach my $row (@{$data_struct{$an}}){
		$counter++;
		print "Working on fastq number: " . $row->fnum() . " ... ";
		my $fq1 = $row->fq1();
		my $fq2 = $row->fq2();
		my $exbase = $row->fc() . "_" . $an . "_" . $row->fnum() . "_$counter";
		
		my ($sam, $cmd) = cmd_str($fq1, $fq2, $exbase);
		$row->sam($sam);
		
		fork { cmd => $cmd, timeout => 3000, max_proc => 20};
		sleep(2);
		
		
	}
}
my $pidcount = waitall();
print "Waited for $pidcount processes\n";
foreach my $an (keys(%data_struct)){
	foreach my $row (@{$data_struct{$an}}){
		my @inserts;
		my $sam = $row->sam();
		open(IN, "< $sam") || die "Could not open $sam!\n";
		while(my $line = <IN>){
			chomp $line;
			if($line =~ /^\@/){next;}
			my @segs = split(/\t/, $line);
			push(@inserts, abs($segs[8]));
		}
		close IN;
		my $rev = remove_upper_lower_deci(\@inserts);
		my $avg = average($rev);
		my $sd = stdev($rev);
		
		$row->insert(int($avg));
		$row->sd($sd);
		print "done\t$an\r";
	}
}
print "\n";

open(OUT, "> spreadsheet_starter.txt");
foreach my $an (keys(%data_struct)){
	foreach my $row(@{$data_struct{$an}}){
		print OUT $row->fq1() . "\t" . $row->fq2 . "\t" . $row->insert() . "\t" . $row->sd() . "\tp\t$an\n";
	}
}
close OUT;

exit;

sub cmd_str{
	my ($fq1, $fq2, $exbase) = @_;
	my $tmpdir = "/mnt/iscsi/vnx_gliu_7/tmp";
	my $cmd = "gunzip -c $fq1 | head -n 400000 > $tmpdir/$exbase.fq1;";
	$cmd .= " gunzip -c $fq2 | head -n 400000 > $tmpdir/$exbase.fq2;";
	$cmd .= " bwa aln reference/bosTau6.fa $tmpdir/$exbase.fq1 > $tmpdir/$exbase.fq1.sai;";
	$cmd .= " bwa aln reference/bosTau6.fa $tmpdir/$exbase.fq2 > $tmpdir/$exbase.fq2.sai;";
	$cmd .= " bwa sampe reference/bosTau6.fa $tmpdir/$exbase.fq1.sai $tmpdir/$exbase.fq2.sai $tmpdir/$exbase.fq1 $tmpdir/$exbase.fq2 > $tmpdir/$exbase.aln.sam;";
	$cmd .= " rm $tmpdir/$exbase.fq1 $tmpdir/$exbase.fq2 $tmpdir/$exbase.fq1.sai $tmpdir/$exbase.fq2.sai";
	
	my $expect = "$tmpdir/$exbase.aln.sam";
	return($expect, $cmd);
}

sub remove_upper_lower_deci{
	my ($aref) = @_;
	my $ret;
	my @sorted = sort{$b <=> $a} @{$aref};
	my $quint = int(scalar(@sorted) * 0.1);
	for(my $x = $quint; $x < scalar(@sorted) - $quint; $x++){
		push(@{$ret}, $sorted[$x]);
	}
	return $ret;
}

sub stdev {
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
	if(scalar(@numb) == 0){
		return 0;
	}
	my $total3 = 0;
	foreach my $num1 (@numb) {
	$total3 += $num1;
	}
	my $mean3 = $total3 / (scalar @numb);
	return $mean3;
}