#!/usr/bin/perl
package parse_config_file;
use Class::Struct;

Class::Struct::struct( SingleFile =>{
	'mrsfast_sam' => '$',
	'bwa_bam' => '$',
});

Class::Struct::struct( sampleName => {
	'sample_name' => '$',
	'rgnum' => '$',
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

# The spreadsheet
Class::Struct::struct( fastqFiles => {
	'file_name' => '$',
	'filelist' => '@', # Contains an array of sampleName
	'has_samples' => '$', # Boolean for determining if align_org is filled
});

# Input genome files for use in the pipeline
Class::Struct::struct( genomeFiles => {
	'base_ref_genome' => '$',
	'alt_ref_genome' => '$',
	'gap_file' => '$',
	'fq_line_limit' => '$', # Default is 100000
	'align_org' => '%', # Hash of array filled in later with reorganized files based on the base name of the animal
});

# Necessary directories for running programs
Class::Struct::struct( utilityPaths => {
	'bin_dir' => '$',
	'java_path' => '$',
	'gatk_path' => '$',
	'mrsfast_path' => '$',
	'command_dir' => '$',
	'base_output_dir' => '$',
	'fastq_output_dir' => '$',
	'process_restart_log' => '$',
});

# Boolean flags that are designed to control aspects of the run
Class::Struct::struct( runFlags => { 
	'is_threaded_not_lsf' => '$',
	'gatk_call_snps_indels' => '$',
	'split_fastq' => '$',
	'align_mrsfast_splitread' => '$',
	'call_mrsfast_cnvs' => '$',
	'merge_cnvs' => '$',
	'clean_mrsfast_bam' => '$',
	'clean_bwa_bam' => '$',
});


sub new {
	my $class = shift(@_);
	my $fq = {}; #Empty hash to be stored with files later; Only one default value stored
	my $gf = genomeFiles->new('fq_line_limit' => 100000);
	my $up = utilityPaths->new();
	my $rf = runFlags->new();
	
	# Setting up default runFlags
	$rf->is_threaded_not_lsf(0);
	$rf->gatk_call_snps_indels(1);
	$rf->split_fastq(1);
	$rf->align_mrsfast_splitread(1);
	$rf->call_mrsfast_cnvs(1);
	$rf->merge_cnvs(1);
	$rf->clean_mrsfast_bam(0);
	$rf->clean_bwa_bam(0);
	
	my $self = {
		'fastqFiles' => $fq,
		'genomeFiles' => $gf,
		'utilityPaths' => $up,
		'runFlags' => $rf,
	};
	bless $self, $class;
	return $self;
}

sub open_file {
	my $self = shift;
	my $file = shift;
	
	my ($path_switch, $sample_switch, $output_switch, $flag_switch);
	$path_switch = 0; $sample_switch = 0; $output_switch = 0;
	open (IN, "< $file") || die "Could not open input configuration file!\n";
	while (my $line = <IN>){
		if ($line =~ /^\s+/){ next;} # Skipping comments and empty lines
		if ($line =~ /^#/){ next;}
		$line =~ s/\r//g;
		chomp $line;
		
		# Switch setting
		if ($line =~ /\[path\]/ && !($sample_switch) && !($output_switch) && !($flag_switch)){
			$path_switch = 1; next;
		}elsif($line =~ /\[files\]/ && !($path_switch) && !($output_switch) && !($flag_switch)){
			$sample_switch = 1; next;
		}elsif($line =~ /\[output\]/ && !($path_switch) && !($sample_switch) && !($flag_switch)){
			$output_switch = 1; next;
		}elsif($line =~ /\[flags\]/ && !($path_switch) && !($sample_switch) && !($output_switch)){
			$flag_switch = 1; next;
		}
		
		if($line =~ /\[\/path\]/){
			$path_switch = 0; next;
		}elsif($line =~ /\[\/files\]/){
			$sample_switch = 0; next;
		}elsif($line =~ /\[\/output\]/){
			$output_switch = 0; next;
		}elsif($line =~ /\[\/flags\]/){
			$flag_switch = 0; next;
		}
		
		# Line processing
		if($path_switch && !($sample_switch) && !($output_switch) && !($flag_switch)){
			_process_path_line($self, $line);
		}elsif($sample_switch && !($path_switch) && !($output_switch) && !($flag_switch)){
			_process_sample_line($self, $line);
		}elsif($output_switch && !($path_switch) && !($sample_switch) && !($flag_switch)){
			_process_output_line($self, $line);
		}elsif($flag_switch && !($sample_switch) && !($output_switch) && !($path_switch)){
			_process_flag_line($self, $line);
		}else{
			print STDERR "Encountered error with malformed configuration file\n";
			print STDERR "User claimed this line:\t$line\n";
			print STDERR "Was within the following fields:\n";
			if($path_switch){
				print STDERR "[path]\n";
			}
			if($sample_switch){
				print STDERR "[files]\n";
			}
			if($output_switch){
				print STDERR "[output]\n";
			}
			if($flag_switch){
				print STDERR "[flags]\n";
			}
			print STDERR "Exiting prematurely from pipeline\n";
			exit;
		}
	}
}

sub _process_path_line {
	my $self = shift;
	my $line = shift;
	
	my @segs = split(/=/, $line);
	if(scalar(@segs) != 2){
		print STDERR "Possible empty field in the [path] field in config file:\n";
		print STDERR "$line\n";
	}
	
	if($segs[0] eq "mrsfast_path"){ $self->{'utilityPaths'}->mrsfast_path($segs[1]);}
	elsif($segs[0] eq "base_ref_genome"){$self->{'genomeFiles'}->base_ref_genome($segs[1]);}
	elsif($segs[0] eq "alt_ref_genome"){$self->{'genomeFiles'}->alt_ref_genome($segs[1]);}
	elsif($segs[0] eq "gap_file"){$self->{'genomeFiles'}->gap_file($segs[1]);}
	elsif($segs[0] eq "gatk_path"){$self->{'utilityPaths'}->gatk_path($segs[1]);}
	elsif($segs[0] eq "java_path"){$self->{'utilityPaths'}->java_path($segs[1]);}
	elsif($segs[0] eq "bin_dir"){$self->{'utilityPaths'}->bin_dir($segs[1]);}
	else{
		print STDERR "Could not interpret this line in the [path] field in config file:\n";
		print STDERR "$line\n";
	}
}

sub _process_sample_line {
	my $self = shift;
	my $line = shift;
	
	my @segs = split(/[\t\s,]/, $line);
	if(scalar(@segs) != 6){
		print STDERR "Found " . scalar(@segs) . "segments in the [file] field in configuration file\n";
		print STDERR "This is the offending line:\n";
		print STDERR "$line\n";
	}
	my $temp;
	if(!exists($self->{'fastqFiles'}->{$segs[4]})){
		$self->{'fastqFiles'}->{$segs[4]} = fastqFiles->new('has_samples' => 1);
	}else{
		$temp = $self->{'fastqFiles'}->{$segs[4]}->filelist();
	}
	
	my $samp = sampleName->new();
	$samp->fastq_1($segs[0]);
	$samp->fastq_2($segs[1]);
	$samp->insert_size($segs[2]);
	$samp->insert_stdev($segs[3]);
	$samp->lib_type($segs[4]);
	$samp->sample_name($segs[5]);
		
	push(@{$temp}, $samp);
	$self->{'fastqFiles'}->{$segs[4]}->filelist($temp);
}

sub _process_output_line {
	my $self = shift;
	my $line = shift;
	
	my @segs = split(/=/, $line);
	if(scalar(@segs) != 2){
			print STDERR "Possible empty field in the [output] field in config file:\n";
			print STDERR "$line\n";
	}
	
	if($segs[0] eq "command_dir"){ $self->{'utilityPaths'}->command_dir($segs[1]);}
	elsif($segs[0] eq "base_output_dir"){ $self->{'utilityPaths'}->base_output_dir($segs[1]);}
	elsif($segs[0] eq "fastq_output_dir"){ $self->{'utilityPaths'}->fastq_output_dir($segs[1]);}
	elsif($segs[0] eq "process_restart_log"){ $self->{'utilityPaths'}->process_restart_log($segs[1]);}
	elsif($segs[0] eq "fastq_line_limit"){ $self->{'genomeFiles'}->fq_line_limit($segs[1]);}
	else{
		print STDERR "Could not interpret this line in the [output] field in config file\n";
		print STDERR "$line\n";
	}
}

sub _process_flag_line{
	my $self = shift;
	my $line = shift;
		
	my @segs = split(/=/, $line);
	if(scalar(@segs) != 2){
			print STDERR "Possible empty field in the [flags] field in config file:\n";
			print STDERR "$line\n";
	}
	
	if($segs[0] eq 'gatk_call_snps_indels'){ $self->{'runFlags'}->gatk_call_snps_indels(int($segs[1]));}
	elsif($segs[0] eq 'is_threaded_not_lsf'){$self->{'runFlags'}->is_threaded_not_lsf(int($segs[1]));}
	elsif($segs[0] eq 'split_fastq'){$self->{'runFlags'}->split_fastq(int($segs[1]));}
	elsif($segs[0] eq 'align_mrsfast_splitread'){$self->{'runFlags'}->align_mrsfast_splitread(int($segs[1]));}
	elsif($segs[0] eq 'call_mrsfast_cnvs'){$self->{'runFlags'}->call_mrsfast_cnvs(int($segs[1]));}
	elsif($segs[0] eq 'merge_cnvs'){$self->{'runFlags'}->merge_cnvs(int($segs[1]));}
	elsif($segs[0] eq 'clean_mrsfast_bam'){$self->{'runFlags'}->clean_mrsfast_bam(int($segs[1]));}
	elsif($segs[0] eq 'clean_bwa_bam'){$self->{'runFlags'}->clean_bwa_bam(int($segs[1]));}
	else{
		print STDERR "Could not interpret this line in the [flags] field in config file:\n";
		print STDERR "$line\n";
	}
}
1;