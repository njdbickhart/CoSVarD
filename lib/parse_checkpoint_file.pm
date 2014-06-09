#!/usr/bin/perl
# This is the module that parses a checkpoint file and retrieves information about the previous run
package parse_checkpoint_file;
use strict;
use Class::Struct;

Class::Struct::struct( CSingleFile =>{
	'mrsfast_sam' => '$',
	'bwa_bam' => '$',
});

Class::Struct::struct( CSRRPFileEntry =>{
	'anchor_file' => '$',
	'split_sam1' => '$',
	'split_sam2' => '$',
	'divet_file' => '$',
	'upper' => '$',
	'lower' => '$',
});

Class::Struct::struct( CSampleName => {
	'sample_name' => '$',
	'fastq_1' => '$',
	'fastq_2' => '$',
	'lib_type' => '$',
	'insert_size' => '$',
	'insert_stdev' => '$',
	'split_fq1_list' => '@',
	'split_fq2_list' => '@',
	'srrp_files' => '@', # contains list of SRRPFileEntry objects
	'single_files' => '@', # contains list of SingleFile objects
});

my $conf; # config file parser pointer
my $pwd; # current working directory

sub new {
	my $class = shift(@_);
	$conf = shift(@_);
	my $checkfile = shift(@_);
	my $fq = {}; #Empty hash to be stored with files later; 
	my $self = {
		'checkfile' => $checkfile,
		'fastqFiles' => $fq,
	};
	bless $self, $class;
	return $self;
}

sub parse_file{
	my ($self, $file) = @_;
	
	# Open file, check for parity
	open(IN, "< $file") || die "[Checkpoint Parser] Could not open $file for parsing\n";
	my $header = <IN>;
	my @wdsplit = split(/=/, $header);
	if($wdsplit[0] eq "dir"){
		$pwd = $wdsplit[1];
	}else{
		print "[Checkpoint Parser] Malformed checkpoint file! Failed at working directory header!\n";
	}
	
	# Run through file and extract information
	while(my $line = <IN>){
		chomp $line;
		
	}
}

1;