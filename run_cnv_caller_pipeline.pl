#!/usr/bin/perl
# This script makes use of the Forks::Super module to create and manage jobs on the cluster.
# It is a simple driving script that will produce STDOUT indicating the progress of the data
# I will take an input text file with instructions on each file (need to distinguish single end and mate pair)
# then run them through the cluster at about 79 processes at a time (master process takes up one slot)
# File types: "s" single end, "m" mate pair, "p" paired end
# Update 12/18/2012: Modified pipeline and reorganized everything to use my config-based setup

use strict;
use File::Basename;
use FileHandle;
use Cwd;
use lib dirname(__FILE__) . '/lib';
use lsf_pipeline_utils;
use thread_pipeline_utils;
use parse_checkpoint_file;
use Class::Struct;
use parse_config_file;
use Forks::Super;
use Getopt::Std;


Class::Struct::struct( SRRPFileEntry =>{
	'anchor_file' => '$',
	'split_sam1' => '$',
	'split_sam2' => '$',
	'divet_file' => '$',
	'upper' => '$',
	'lower' => '$',
});


$Forks::Super::ON_BUSY = "block";
my $usage = "$0 -c <config file> -p <max proc>\n
This script runs the cnv and snp calling pipeline on a series of files on an LSF server

	-c	Config file entry[Mandatory]
	-p	Maximum number of processes to run in the background [default = 20]
	-r	Name of checkpoint file to restart pipeline [Switch]\n";

my %opts;
getopt('cprh', \%opts);
my $svpackage = new lsf_pipeline_utils;
my $tsvpackage = new thread_pipeline_utils;

my $conf = parse_config_file->new();

my %filelist;	# Holds spreadsheet file information; {type}->[row]->[filename1, filename2, insert, SD]
my @job_array; 	# Holds jobs; each job is a command to run a subroutine with args
my @split_list; # Holds split file lists
my $initial_time = time();
my $remove = 1;	# Debugging variable to prevent premature deletion of test files


# Important structures
Class::Struct::struct (bwa_bams => {
	'list_of_bams' => '@',
	'rgstr' => '$',
	'directory' => '$',
	'sample_name' => '$',
});


#########################################
#	Parity checking			#
#########################################

unless((defined($opts{'c'}) ) && !defined($opts{'h'})){
	print $usage;
	exit;
}
my $max_p = 20;
if(defined($opts{'p'})){
	$max_p = $opts{'p'};
}
$conf->open_file($opts{'c'});

mkdir($conf->{'utilityPaths'}->base_output_dir()) || print $! . "\n";
mkdir($conf->{'utilityPaths'}->command_dir()) || print $! . "\n";
mkdir($conf->{'utilityPaths'}->fastq_output_dir()) || print $! . "\n";

# log file check and print beginning of run
my ($logfile_base, $logfile_dir, $logfile_suffix) = fileparse($opts{'c'}, qr/\.[^.]*/);
my $chk_time = file_timestamp();
my $logfile = $logfile_base . ".log";
my $checkfile = $logfile_base . "_$chk_time.checkpoint";
my $log_fh = FileHandle->new();
if(-e $logfile){
	# Place to put the restart log file parser
	print "Already have $logfile ... appending";
	$log_fh->open(">> $logfile");
	print $log_fh "End of Run\n";
	print $log_fh "---------------------------\n---------------------------\n";
	if(defined($opts{'r'})){
		print $log_fh "Restarting from checkpoint file: " . $opts{'r'} . "\n";
	}else{
		print $log_fh "New run, no checkpoint file selected\n";
	}
}else{
	$log_fh->open("> $logfile");
}
my $ts = timestamp_return();
print $log_fh "$ts : Beginning run of config file: $opts{c}";
foreach my $k (keys(%opts)){
	print $log_fh "Command option: $k\tvalue: $opts{$k}\n";
}

my $check_fh = FileHandle->new();
my $check_parse = parse_checkpoint_file->new();
if(-e $checkfile){
	# Terminate early, asking user if he/she wants to restart a run
	print "The checkpoint file for this configuration already exists: $checkfile\n";
	print "Did you want to restart this run? If so, use the -r option.\n";
	print "If you do not, then delete the checkpoint file before restarting. Exiting prematurely...\n";
	exit;
}elsif(defined($opts{'r'})){
	# TODO implement checkpoint parser
}else{
	print "Creating checkpoint file: $checkfile . Use this if you wish to restart this run\n";
	print $log_fh "$ts : Created checkpoint file: $checkfile. Use this if you wish to restart this run\n";
	$check_fh->open("> $checkfile");
	my $pwd = cwd();
	print $check_fh "dir=$pwd\n";
}
#########################################
#		Alignment		#
#########################################

# Work on regular paired end files

my @removalfiles;
my %output_files;
my %iterators;
my %readgroup;
my %bwabamstore; # Contains a hash of bwa_bams objects with the rgnum as the key
if(exists($conf->{'fastqFiles'}->{'p'})){
	if (scalar(@{$conf->{'fastqFiles'}->{'p'}->filelist()}) > 0){
		print "Working on paired end files\n";
		$ts = timestamp_return();
		print $log_fh "$ts : working on paired end files\n";
	}
	
	my $fit = 0;
	foreach my $f (@{$conf->{'fastqFiles'}->{'p'}->filelist()}){
		# SPLIT FQ Conditional
		if($conf->{'runFlags'}->split_fastq()){
			$readgroup{$f->sample_name()} += 1;
			$f->rgnum($readgroup{$f->sample_name()});
			# Iterator hash to keep fastq number iterators unique
			if(exists($iterators{$f->sample_name()})){
				$fit = $iterators{$f->sample_name()};
			}else{
				$fit = 0;
			}
			my ($rit, $filestore1, $filestore2) = $svpackage->split_pe_fastqs(
				$f->fastq_1(), 
				$f->fastq_2(), 
				$conf->{'utilityPaths'}->fastq_output_dir(), 
				$conf->{'genomeFiles'}->fq_line_limit(), 
				$fit
			);
			$f->split_fq1_list($filestore1);
			$f->split_fq2_list($filestore2);
			
			# Print to checkpoint file
			print $check_fh $f->sample_name() . "\tpfq1\t" . $readgroup{$f->sample_name()} . "\t" . join("\t", @{$filestore1}) . "\n";
			print $check_fh $f->sample_name() . "\tpfq2\t" . $readgroup{$f->sample_name()} . "\t" . join("\t", @{$filestore2}) . "\n";
			
			if($rit){
				$iterators{$f->sample_name()} = $rit;
			}else{
				print "Error on file iterator for " . $f->sample_name() . "! $rit " . scalar(@{$filestore1}) . "\n";
			}
		}else{
			# Calculate what files should have been created and check for them.
			# Fill parse config file with sampleName entries, including read group numbers and file lists
		}
		# Determine BWA read group information for later merger and parsing
		my $rgrstr = construct_rg_str($f->sample_name(), "ILLUMINA", $readgroup{$f->sample_name()});
		my $bwa_obj = bwa_bams->new(
			'rgstr' => $rgrstr,
			'directory' => $conf->{'utilityPaths'}->base_output_dir() . "/" . $f->sample_name(),
			'sample_name' => $f->sample_name(),
		);
		push(@{$bwabamstore{$readgroup{$f->sample_name()}}}, $bwa_obj);
		
		# VHSR file output
		my $vhsr_file = FileHandle->new();
		if($conf->{'runFlags'}->align_mrsfast_splitread()){
			my $vhsr_file_str = $conf->{'utilityPaths'}->base_output_dir() . "/" . $f->sample_name() . "/" . $f->sample_name() . "_vhsr.list";
			if(-e $vhsr_file_str){
				# Open existing file for appending
				$vhsr_file->open(">> $vhsr_file_str");
			}else{
				$vhsr_file->open("> $vhsr_file_str");
			}
			# Print to checkpoint file
			print $check_fh $f->sample_name() . "\tsrrp\t$vhsr_file_str\n";
		}
		
		my $i = $fit;
		for (my $k = 0; $k < scalar(@{$f->split_fq1_list()}); $k++){
			my $fq1 = $f->split_fq1_list()->[$k];
			my $fq2 = $f->split_fq2_list()->[$k];
			my ($command, $file_array, $outbase);
			
			# ALIGN split fastqs conditional
			if($conf->{'runFlags'}->align_mrsfast_splitread()){
				if($conf->{'runFlags'}->is_threaded_not_lsf()){
					($command, $file_array, $outbase) = $tsvpackage->t_generate_palign_command(
						$fq1, 
						$fq2, 
						$conf->{'genomeFiles'}->base_ref_genome(), 
						2, 
						$f->insert_size(), 
						$f->insert_stdev(), 
						$conf->{'utilityPaths'}->base_output_dir(), 
						$i, 
						$conf->{'genomeFiles'}->alt_ref_genome(), 
						$f->sample_name(), 
						$conf->{'utilityPaths'}->bin_dir(), 
						$conf->{'utilityPaths'}->command_dir(), 
						$conf->{'utilityPaths'}->java_path(), 
						$readgroup{$f->sample_name()}
					);
				}else{
					($command, $file_array, $outbase) = $svpackage->generate_palign_command(
						$fq1, 
						$fq2, 
						$conf->{'genomeFiles'}->base_ref_genome(), 
						2, 
						$f->insert_size(), 
						$f->insert_stdev(), 
						$conf->{'utilityPaths'}->base_output_dir(), 
						$i, 
						$conf->{'genomeFiles'}->alt_ref_genome(), 
						$f->sample_name(), 
						$conf->{'utilityPaths'}->bin_dir(), 
						$conf->{'utilityPaths'}->command_dir(), 
						$conf->{'utilityPaths'}->java_path(), 
						$readgroup{$f->sample_name()}
					);
				}				
				fork { timeout => 1800, cmd => $command, max_proc => $max_p };
				sleep(3);
				
				# Creating SRRP entry
				my $temp = SRRPFileEntry->new(
					'anchor_file' => $file_array->[3],
					'split_sam1' => $file_array->[0],
					'divet_file' => $file_array->[2],
					'upper' => $f->insert_size() + ($f->insert_stdev() * 3),
					'lower' => 50,
				);
				push(@{$f->srrp_files()}, $temp);
				
				# Print srrp files to flatfile list
				print $vhsr_file $file_array->[0] . "\t" . $file_array->[2] . "\t" . $file_array->[3] . "\t" . $f->insert_size() . "\t" . $f->insert_stdev() .  "\n";
				
				# Print files to checkpoint file
				print $check_fh $f->sample_name() . "\tdoc\t" . $f->rgnum() . "\t" . $file_array->[0] . "\t" . $file_array->[1] . "\n";
				print $check_fh $f->sample_name() . "\tsnp\t" . $f->rgnum() . "\t" . $file_array->[4] . "\n";
			}else{
				#Calculate what files should have been created and check for them
			}			
			$i++;
		}
		if(!exists($conf->{'genomeFiles'}->align_org()->{$f->sample_name()})){
			my %h;
			push(@{$h{$f->sample_name()}}, $f);
			$conf->{'genomeFiles'}->align_org(\%h);
		}else{
			my %h = %{$conf->{'genomeFiles'}->align_org()};
			push(@{$h{$f->sample_name()}}, $f);
			$conf->{'genomeFiles'}->align_org(\%h);
		}
		print "\n";
		$ts = timestamp_return();
		print $log_fh "$ts : Completed work on: " . $f->sample_name() . " " . $f->rgnum() . "\n";
		
		# Now to do cleanup work on the sam files and convert them to bams
		my $command;
		if($conf->{'runFlags'}->is_threaded_not_lsf()){
			($command) = t_format_mrsfast_bam(
				$conf->{'utilityPaths'}->base_output_dir() . "/" . $f->sample_name(),
				$conf->{'genomeFiles'}->base_ref_genome(),
				$conf->{'utilityPaths'}->bin_dir(),
				$conf->{'utilityPaths'}->command_dir(),
				$f->sample_name(),
				1
			);
		}else{
			($command) = format_mrsfast_bam(
				$conf->{'utilityPaths'}->base_output_dir() . "/" . $f->sample_name(),
				$conf->{'genomeFiles'}->base_ref_genome(),
				$conf->{'utilityPaths'}->bin_dir(),
				$conf->{'utilityPaths'}->command_dir(),
				$f->sample_name(),
				1
			);
		}
		fork { timeout => 5000, cmd => $command, max_proc => $max_p };
	}	
	#@removalfiles = ();
	my $endtime = time - $initial_time;
	print "Completed\t" .  "\tfiles in\t" . $endtime . "\n";
	$ts = timestamp_return();
	print $log_fh "$ts : Completed paired end files. Total time: $endtime\n";
}


# Work on mate pair files
if(exists($conf->{'fastqFiles'}->{'m'})){
	if (scalar(@{$conf->{'fastqFiles'}->{'m'}->filelist()}) > 0){
		print "Working on mate-pair files\n";
	}
	my $fit = 0;
	foreach my $f (@{$conf->{'fastqFiles'}->{'m'}->filelist()}){
		if($conf->{'runFlags'}->split_fastq()){
			$readgroup{$f->sample_name()} += 1;
			$f->rgnum($readgroup{$f->sample_name()});
			# Iterator hash to keep fastq number iterators unique
			if(exists($iterators{$f->sample_name()})){
				$fit = $iterators{$f->sample_name()};
			}else{
				$fit = 0;
			}
			my ($rit, $filestore1, $filestore2) = $svpackage->split_se_fastqs(
				$f->fastq_1(), 
				$conf->{'utilityPaths'}->fastq_output_dir(), 
				$conf->{'genomeFiles'}->fq_line_limit(), 
				$fit
			);
			$f->split_fq1_list($filestore1);
			if($rit){
				$iterators{$f->sample_name()} = $rit;
			}else{
				print "Error on file iterator for " . $f->sample_name() . "! $rit " . scalar(@{$filestore1}) . "\n";
			}
			my ($eit, $filestore1, $filestore2) = $svpackage->split_se_fastqs(
				$f->fastq_2(), 
				$conf->{'utilityPaths'}->fastq_output_dir(), 
				$conf->{'genomeFiles'}->fq_line_limit(), 
				$rit
			);
			if($eit){
				$iterators{$f->sample_name()} = $eit;
			}else{
				print "Error on file iterator for " . $f->sample_name() . "! $eit " . scalar(@{$filestore2}) . "\n";
			}
			$f->split_fq1_list($filestore2);
			
			# Print to checkpoint file
			print $check_fh $f->sample_name() . "\tmfq1\t" . $readgroup{$f->sample_name()} . "\t" . join("\t", @{$filestore1}) . "\n";
			print $check_fh $f->sample_name() . "\tmfq2\t" . $readgroup{$f->sample_name()} . "\t" . join("\t", @{$filestore2}) . "\n";
		}else{
			# Calculate what files should have been created and check for them.
		}
		# Determine BWA read group information for later merger and parsing
		my $rgrstr = construct_rg_str($f->sample_name(), "ILLUMINA", $readgroup{$f->sample_name()});
		my $bwa_obj = bwa_bams->new(
			'rgstr' => $rgrstr,
			'directory' => $conf->{'utilityPaths'}->base_output_dir() . "/" . $f->sample_name(),
			'sample_name' => $f->sample_name(),
		);
		push(@{$bwabamstore{$readgroup{$f->sample_name()}}}, $bwa_obj);
		
		my $i = $fit;
		for (my $k = 0; $k < scalar(@{$f->split_fq1_list()}); $k++){
			my ($command, $file_array, $outbase);
			if($conf->{'runFlags'}->align_mrsfast_splitread()){
				if($conf->{'runFlags'}->is_threaded_not_lsf()){
					($command, $file_array, $outbase) = $tsvpackage->t_generate_salign_command(
						$f->split_fq1_list()->[$k], 
						$conf->{'genomeFiles'}->base_ref_genome(), 
						2, 
						$conf->{'genomeFiles'}->base_output_dir(), 
						$i, 
						$f->sample_name(), 
						$conf->{'utilityPaths'}->bin_dir(), 
						$conf->{'utilityPaths'}->command_dir(),
						$conf->{'utilityPaths'}->java_path(),
						$readgroup{$f->sample_name()},
						$conf->{'genomeFlies'}->alt_ref_genome()
					);
				}else{
					($command, $file_array, $outbase) = $svpackage->generate_salign_command(
						$f->split_fq1_list()->[$k], 
						$conf->{'genomeFiles'}->base_ref_genome(), 
						2, 
						$conf->{'genomeFiles'}->base_output_dir(), 
						$i, 
						$f->sample_name(), 
						$conf->{'utilityPaths'}->bin_dir(), 
						$conf->{'utilityPaths'}->command_dir(),
						$conf->{'utilityPaths'}->java_path(),
						$readgroup{$f->sample_name()},
						$conf->{'genomeFlies'}->alt_ref_genome()
					);
				}
				fork { timeout => 1800, cmd => $command, max_proc => $max_p };
				sleep(3);
				my $temp = parse_config_file::SingleFile->new(
					'mrsfast_sam' => $file_array->[0],
					'bwa_bam' => $file_array->[1],
				);
				push(@{$f->single_files()}, $temp);
				push(@{$f->mrsfast_bam_list()}, $file_array->[0]);
				push(@{$f->bwa_bam_list()}, $file_array->[1]);
				# Print files to checkpoint file
				print $check_fh $f->sample_name() . "\tdoc\t" . $f->rgnum() . "\t" . $file_array->[0] . "\n";
				print $check_fh $f->sample_name() . "\tsnp\t" . $f->rgnum() . "\t" . $file_array->[1] . "\n";
			}else{
				# Calculate what files should have been created and check for them
			}			
			
			$i++;
		}
		print "\n";
		#push(@removalfiles, @filestore);
		if(!exists($conf->{'genomeFiles'}->align_org()->{$f->sample_name()})){
			my %h;
			push(@{$h{$f->sample_name()}}, $f);
			$conf->{'genomeFiles'}->align_org(\%h);
		}else{
			my %h = %{$conf->{'genomeFiles'}->align_org()};
			push(@{$h{$f->sample_name()}}, $f);
			$conf->{'genomeFiles'}->align_org(\%h);
		}
		
		my $command;
		if($conf->{'runFlags'}->is_threaded_not_lsf()){
			($command) = t_format_mrsfast_bam(
				$conf->{'utilityPaths'}->base_output_dir() . "/" . $f->sample_name(),
				$conf->{'genomeFiles'}->base_ref_genome(),
				$conf->{'utilityPaths'}->bin_dir(),
				$conf->{'utilityPaths'}->command_dir(),
				$f->sample_name(),
				0
			);
		}else{
			($command) = format_mrsfast_bam(
				$conf->{'utilityPaths'}->base_output_dir() . "/" . $f->sample_name(),
				$conf->{'genomeFiles'}->base_ref_genome(),
				$conf->{'utilityPaths'}->bin_dir(),
				$conf->{'utilityPaths'}->command_dir(),
				$f->sample_name(),
				0
			);
		}
		fork { timeout => 5000, cmd => $command, max_proc => $max_p };
	}

	
	#@removalfiles = ();
	my $endtime = time - $initial_time;
	print "Completed\t" . scalar(@{$filelist{'m'}}) . "\tfiles in\t" . $endtime . "\n";	
	
}


# Work on single end files
if(exists($conf->{'fastqFiles'}->{'s'})){
	if (scalar(@{$conf->{'fastqFiles'}->{'s'}->filelist()}) > 0){
		print "Working on single end files\n";
	}
	my $fit = 0;
	foreach my $f (@{$conf->{'fastqFiles'}->{'s'}->filelist()}){
		if($conf->{'runFlags'}->split_fastq()){
			$readgroup{$f->sample_name()} += 1;
			# Iterator hash to keep fastq number iterators unique
			if(exists($iterators{$f->sample_name()})){
				$fit = $iterators{$f->sample_name()};
			}else{
				$fit = 0;
			}
			my ($rit, $filestore1, $filestore2) = $svpackage->split_se_fastqs(
				$f->fastq_1(), 
				$conf->{'utilityPaths'}->fastq_output_dir(), 
				$conf->{'genomeFiles'}->fq_line_limit(), 
				$fit
			);
			$f->split_fq1_list($filestore1);
			if($rit){
				$iterators{$f->sample_name()} = $rit;
			}else{
				print "Error on file iterator for " . $f->sample_name() . "! $rit " . scalar(@{$filestore1}) . "\n";
			}
			
			# Print to checkpoint file
			print $check_fh $f->sample_name() . "\tsfq1\t" . $readgroup{$f->sample_name()} . "\t" . join("\t", @{$filestore1}) . "\n";
						
		}else{
			# Calculate what files should have been created and check for them.
		}
		# Determine BWA read group information for later merger and parsing
		my $rgrstr = construct_rg_str($f->sample_name(), "ILLUMINA", $readgroup{$f->sample_name()});
		my $bwa_obj = bwa_bams->new(
			'rgstr' => $rgrstr,
			'directory' => $conf->{'utilityPaths'}->base_output_dir() . "/" . $f->sample_name(),
			'sample_name' => $f->sample_name(),
		);
		push(@{$bwabamstore{$readgroup{$f->sample_name()}}}, $bwa_obj);
		my $i = $fit;
		for(my $k = 0; $k < scalar(@{$f->split_fq1_list()}); $k++){
			my ($command, $file_array, $outbase); 
			if($conf->{'runFlags'}->align_mrsfast_splitread()){
				if($conf->is_threaded_not_lsf()){
					($command, $file_array, $outbase) = $tsvpackage->t_generate_salign_command(
						$f->split_fq1_list()->[$k], 
						$conf->{'genomeFiles'}->base_ref_genome(), 
						2, 
						$conf->{'genomeFiles'}->base_output_dir(), 
						$i, 
						$f->sample_name(), 
						$conf->{'utilityPaths'}->bin_dir(), 
						$conf->{'utilityPaths'}->command_dir(),
						$conf->{'utilityPaths'}->java_path(),
						$readgroup{$f->sample_name()},
						$conf->{'genomeFlies'}->alt_ref_genome()
					);
				}else{
					($command, $file_array, $outbase) = $svpackage->generate_salign_command(
						$f->split_fq1_list()->[$k], 
						$conf->{'genomeFiles'}->base_ref_genome(), 
						2, 
						$conf->{'genomeFiles'}->base_output_dir(), 
						$i, 
						$f->sample_name(), 
						$conf->{'utilityPaths'}->bin_dir(), 
						$conf->{'utilityPaths'}->command_dir(),
						$conf->{'utilityPaths'}->java_path(),
						$readgroup{$f->sample_name()},
						$conf->{'genomeFlies'}->alt_ref_genome()
					);
				}
				fork { timeout => 1800, cmd => $command, max_proc => $max_p };
				sleep(3);
				my $temp = parse_config_file::SingleFile->new(
					'mrsfast_sam' => $file_array->[0],
					'bwa_bam' => $file_array->[1],
				);
				push(@{$f->single_files()}, $temp);
				push(@{$f->mrsfast_bam_list()}, $file_array->[0]);
				push(@{$f->bwa_bam_list()}, $file_array->[1]);
				# Print files to checkpoint file
				print $check_fh $f->sample_name() . "\tdoc\t" . $f->rgnum() . "\t" . $file_array->[0] . "\n";
				print $check_fh $f->sample_name() . "\tsnp\t" . $f->rgnum() . "\t" . $file_array->[1] . "\n";
			}else{
				# Calculate what files should have been created and check for them
			}
			
			$i++;
		}
		if(!exists($conf->{'genomeFiles'}->align_org()->{$f->sample_name()})){
			my %h;
			push(@{$h{$f->sample_name()}}, $f);
			$conf->{'genomeFiles'}->align_org(\%h);
		}else{
			my %h = %{$conf->{'genomeFiles'}->align_org()};
			push(@{$h{$f->sample_name()}}, $f);
			$conf->{'genomeFiles'}->align_org(\%h);
		}
		print "\n";
		
		my $command;
		if($conf->{'runFlags'}->is_threaded_not_lsf()){
			($command) = t_format_mrsfast_bam(
				$conf->{'utilityPaths'}->base_output_dir() . "/" . $f->sample_name(),
				$conf->{'genomeFiles'}->base_ref_genome(),
				$conf->{'utilityPaths'}->bin_dir(),
				$conf->{'utilityPaths'}->command_dir(),
				$f->sample_name(),
				0
			);
		}else{
			($command) = format_mrsfast_bam(
				$conf->{'utilityPaths'}->base_output_dir() . "/" . $f->sample_name(),
				$conf->{'genomeFiles'}->base_ref_genome(),
				$conf->{'utilityPaths'}->bin_dir(),
				$conf->{'utilityPaths'}->command_dir(),
				$f->sample_name(),
				0
			);
		}
		fork { timeout => 5000, cmd => $command, max_proc => $max_p };
	}

	
	#@removalfiles = ();
	my $endtime = time - $initial_time;
	print "Completed\t" . scalar(@{$filelist{'s'}}) . "\tfiles in\t" . $endtime . "\n";
}
my $pidwait = waitall();
if($pidwait > 0){
	# Catch any remaining jobs
	print "Had to wait for $pidwait child processes\n";
}
#########################################
#	   SNP and INDEL calling	#
#########################################
if($conf->{'runFlags'}->gatk_call_snps_indels()){
# All processed files should be in %bwabamstore, so iterate through them and do a merger
my $snp_dir = $conf->{'utilityPaths'}->base_output_dir() . "/snp";
mkdir("$snp_dir") || print "$!\n";

my @GATK_raw_bams;
foreach my $rgnum (keys(%bwabamstore)){
	foreach my $rows (@{$bwabamstore{$rgnum}}){
		my $dir = $rows->directory();
		my $samp_name = $rows->sample_name();
		my $bindir = $conf->{'utilityPaths'}->bin_dir();
		my $threads = $conf->{'runFlags'}->is_threaded_not_lsf();
		my $usage = "$0 -l \<number of files to merge\> -i \<input path to search\> -b \<bin dir\> -o \<output name\> -p \<max processors\> -f \<BOOLEAN: is threaded?\>";
		system("$bindir/merge_bams_sort_index.pl -i $dir/$samp_name.$rgnum -b $bindir -o $dir/$samp_name.$rgnum.full.sorted.merged.bam -p $max_p -f $threads");
		push(@GATK_raw_bams, "$dir/$samp_name.$rgnum.full.sorted.merged.bam");
		# print to checkpoint file
		print $check_fh "$samp_name\tsnp\t$rgnum\t$dir/$samp_name.$rgnum.full.sorted.merged.bam\n";
	}
}
# Bams are ready for GATK
# TODO implement GATK pipeline run
my $java = $conf->{'utilityPaths'}->java_path() . "/java";
my $gatk = $conf->{'utilityPaths'}->gatk_path() . "/GenomeAnalysisTK.jar";
my $ref = $conf->{'genomeFiles'}->alt_ref_genome();
my $jopt = "-Xmx9g";
my $indeltargetint = "$snp_dir/indeltarget.intervals";
my $realign_bam = "$snp_dir/indel_realigned_project.bam";
my $reduce_bam = "$snp_dir/indel_realigned_project_reduced.bam";
my $raw_snps = "$snp_dir/genotyper_raw_snp_indel.vcf";
my $filtered_snps = "$snp_dir/genotyper_filtered_snp_indel_calls.vcf";

my $input_bam_str;
foreach my $bam (@GATK_raw_bams){
	$input_bam_str .= "-I $bam ";
}

# Phase one, create analysis ready reads per sample
if($conf->{'runFlags'}->is_threaded_not_lsf()){
	system("$java $jopt -jar $gatk -R $ref -T RealignerTargetCreator -nt $max_p $input_bam_str -o $indeltargetint");
	system("$java $jopt -jar $gatk -R $ref -T IndelRealigner $input_bam_str -targetIntervals $indeltargetint -o $realign_bam");
	# Now to reduce the size of that massive bam file and then delete it
	system("$java $jopt -jar $gatk -R $ref -T ReduceReads -I $realign_bam -o $reduce_bam");
	#system("rm $realign_bam");
}else{
	
}
# Phase two, initial discovery and genotyping
if($conf->{'runFlags'}->is_threaded_not_lsf()){
	system("$java $jopt -jar $gatk -R $ref -nt $max_p -T UnifiedGenotyper -I $reduce_bam -o $raw_snps");
}else{
	
}

# Phase three, filter the variants
if($conf->{'runFlags'}->is_threaded_not_lsf()){
	# TODO add variant filter conditionals
	system("$java $jopt -jar $gatk -R $ref -T VariantFiltration -o $filtered_snps --variant $raw_snps ");
}else{
	
}

# Phase four, annotate remaining variants
# TODO implement this


}else{
	# Do nothing for now
}

#########################################
#		CNV calling		#
#########################################


# Generating chromosome lengths for variationhunter
if($conf->{'runFlags'}->call_mrsfast_cnvs()){
my $chrlens = $conf->{'utilityPaths'}->base_ref_genome() . ".lens";
if(-e $chrlens){
	print "Chrlen file exists, skipping creation\n";
}else{
	print "Calculating chromosome lengths before analysis loop\n";
	system($conf->{'utilityPaths'}->bin_dir() . "/calculate_chr_lengths_from_fasta.pl $opts{r}");
}
# Now, working on main routines for CNV calling
my %doc_files; # {animal}->{autodup/del|xdup/del} = file
my %split_files; # {animal} = events.final
my %vh_files; #{animal} = [final.SV]
foreach my $bnames (keys(%{$conf->{'genomeFiles'}->align_org()})){
	#TODO Update this with my new DOC and RPSR programs
	
	
}
my $endtime = time - $initial_time;
print "Completed main algorithm loop in $endtime seconds\n";
my $pidwait = waitall();
if($pidwait > 0){
	# Catch any remaining jobs
	print "Had to wait for $pidwait child processes\n";
}
}
# Premature exit while I work out the merger routine
exit;

#########################################
#		Merger 			#
#########################################

# Now, formating the files before merger

# doc bed: chr, start, end, gain/loss, copynumber
# vh bed: chr, insidestart, insideend, outsidestart, outsideend, ins/del/inv/, support, sumprob
# split bed: chr, start, end, del/inv, balanced support, unbalanced support, anchor edit, split edit
my %output_files;
my %finalmergefiles; #{name} = [doc, vh, split]
foreach my $bnames (keys(%output_files)){
	mkdir("$opts{o}/$bnames/merger") || print "$!\n";

	
}
my $pidwait = waitall();
if($pidwait > 0){
	# Catch any remaining jobs before merging all the files
	print "Had to wait for $pidwait child processes\n";
}

# Now, merge the files into expanded bed format 
# output1: chr outsidestart outsideend cnvtypestr score + insidestart insideend color
# output2: chr, insidestart, insideend, outsidestart, outsideend, cnvtypstr, copynumber, confidence, docconsis, splitprob, vhprob
# cnvtypestr = [gain/loss/ins/inv/del]_[dsv]_cnv#
foreach my $bnames (keys(%output_files)){
	
}

my $endtime = time - $initial_time;
print "Completed main algorithm loop in $endtime seconds\n";
my $pidwait = waitall();
if($pidwait > 0){
	# Catch any remaining jobs
	print "Had to wait for $pidwait child processes\n";
}


#########################################
#		Intersection		#
#########################################
# to be created

exit;	

sub progress_bar {
	my ( $got, $total, $width, $char ) = @_;
	$width ||= 25; $char ||= '=';
	my $num_width = length $total;
	printf ("|%-25s${width}| Queued %${num_width}s files of %s (%.2f%%)\n", $char x (($width-1)*$got/$total). '>', $got, $total, 100*$got/+$total);
}
sub construct_rg_str{
	my ( $bnames, $plat, $num ) = @_;
	my $sm = "$bnames.$num";
	my $rgstr = "\@RG\tID:$bnames\tPL:$plat\tSM:$sm";
	return $rgstr;
}
sub timestamp_return{
	my $now = localtime();
	return $now;
}
sub file_timestamp{
	my @now = localtime(time);
	return "$now[5]_$now[4]_$now[3]_$now[2]_$now[1]";
}