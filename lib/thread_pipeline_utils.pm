#!/usr/bin/perl
# These are subroutines separated from the original pipeline to unclutter it
# These subroutines are designed to fork out the actual command and not to submit an LSF job
package thread_pipeline_utils;
use warnings;
use strict;

sub new {
	my ($class,%arg) = @_;
	my $self={};
	bless($self,$class || ref($class));
	return $self;
}

sub t_samtools_merger {
	require strict;

	my $prefix = shift(@_);
	my $num = shift(@_);
	my @mergefiles = @_;
	chomp(@mergefiles);
	chomp $num;
	chomp $prefix;

	my $mergestr = join(" ", @mergefiles); 
	my $output =  "$prefix.$num.merge.sorted.bam";
	system("samtools merge $output $mergestr");
}
sub t_final_total_merger_command {
	my ($self, $doc, $vh, $sr, $out, $basename, $cnfile, $bin_dir, $command_dir) = @_;
	my $command = "$bin_dir/merge_final_cnv_calls.pl -d $doc -p $vh -s $sr -c $cnfile -o $out";
	return $command;
}

sub t_format_split_command{
	my ($self, $input, $output, $basename, $bin_dir, $command_dir) = @_;
	my $command = "$bin_dir/merge_split_read_output.pl $input $output";
	return $command;
}

sub t_format_doc_command {
	my ($self, $autodups, $autodels, $xdups, $xdels, $cnfile, $finaldocout, $basename, $bin_dir, $command_dir) = @_;
	my $command = "$bin_dir/merge_doc_output.pl -i $autodups ";
	$command .= "-d $autodels -x $xdups -y $xdels -c $cnfile -o $finaldocout";
	return $command; 
}

sub t_vh_merge_command {
	my ($self, $flistref, $output, $basename, $bin_dir, $command_dir) = @_;
	my $command = "$bin_dir/merge_vh_output.pl ";
	$command .= join(" ", @{$flistref});
	$command .= " $output";
	return $command; 
}

sub t_finalize_split_command {
	my ($self, $psplitfolder , $outfile, $basename, $bin_dir, $command_dir, $java_dir) = @_;
	my $command = "$bin_dir/finalize_split_read.pl -i $psplitfolder -o $outfile -b $basename -B $bin_dir -j $java_dir";
	return $command;	
}

sub t_generate_vh_command {
	my ($self, $pvhfolder, $insert_size, $insert_sd, $chr_length, $gap_file, $basename, $itname, $bin_dir, $command_dir) = @_;
	
	my $full_sv_events = "$pvhfolder/$basename.$itname.final.SV";
	my $max = $insert_size + (3 * $insert_sd);
	my $min = $insert_size - (2 * $insert_sd);
	if($min < 100){
		$min = 100;
	}
	
	my $command = "$bin_dir/run_variation_hunter_pipeline.pl -i $pvhfolder ";
	$command .= "-n $basename -m $min -M $max -p 0.001 -t 4 -c $chr_length -g $gap_file -I $itname -b $bin_dir";
	return ($command, $full_sv_events);
}


sub t_generate_wssd_command {
	my ($self, $folderpath, $outprefix, $reference, $basename, $bin_dir, $command_dir, $java_dir) = @_;

	# predicted files
	my $file1 = "$folderpath/doc_wins/$outprefix.file1.bed";
	my $file2 = "$folderpath/doc_wins/$outprefix.file2.bed";
	my $file3 = "$folderpath/doc_wins/$outprefix.file3.bed";
	my $autodups = "$folderpath/doc_wins/$outprefix.file1.bed.final.wssd";
	my $autodels = "$folderpath/doc_wins/$outprefix.file1.bed.final.deletions.tab";
	my $xdups = "$folderpath/doc_wins/$outprefix.file1.bed.X.final.wssd";
	my $xdels = "$folderpath/doc_wins/$outprefix.file1.bed.final.X.deletions.tab";
	my $cnfile = "$folderpath/doc_wins/$outprefix.file3.bed.gc.depth.normalized.CN";

	my $command = "$bin_dir/alkan_pipeline_no_controls_placeholder.pl ";
	$command .= "--File1 $file1 --File2 $file2 --File3 $file3 --bed $outprefix --hits $folderpath --ref $reference --java $java_dir --bin $bin_dir";
	
	return $command, $autodups, $autodels, $xdups, $xdels, $cnfile;
}

sub t_generate_palign_command {
	my ($self, $f1, $f2, $reference, $edit_dis, $insert, $sd, $out, $i, $alternate, $outbase, $bin_dir, $command_dir, $java_dir, $rgnum) = @_;
		
	my $bnames = $outbase;
	my $sr_out = "$out/$outbase/split_read";
	my $working = "$sr_out/$bnames.$i.sr";
	my $divet = "$working.divet.vh";
	my $single = "$working.single.txt";
	
	my $splitfq1 = "$out/$outbase/$bnames.$i.1.split.fq";
	my $splitfq2 = "$out/$outbase/$bnames.$i.2.split.fq";
	my $mrsfastsplitsam1 = "$working.mrsfast.bam";
	my $dup = "$out/$outbase/$bnames.$i.clean.sort.nodup.bam";

	# Organized list of return files
	my @return = ($mrsfastsplitsam1, "", $divet, $single, $dup);
			
	my $l = "50";
	my $m = "500000";
	my $u = $insert + (3 * $sd);
	my $x = "50";

	my $command = "perl $bin_dir/run_alignment_programs_pe.pl -1 $f1 -2 $f2 -o $out -r $reference -e $edit_dis ";
	$command .= "-l $l -m $m -u $u -n $i -v Xmx3G -a $alternate -c $outbase -b $bin_dir -j $java_dir -g $rgnum";
	
	open(OUT, "> $command_dir/P_$bnames\_align.$i.sh") || print "Could not create command output! $!\n";
	print OUT $command;
	close OUT;
	
	# Returning bsub command and list of anticipated files for merger
	return $command, \@return, $outbase;
}

sub t_generate_salign_command {
	my ($self, $file, $reference, $edit_dis, $out, $i, $outbase, $bin_dir, $command_dir, $java_dir, $rgnum, $alternate) = @_;
	my @basename = split(/\//, $file);
	my @split_name = split(/\./, $basename[-1]);

	my $bnames = $outbase;

	# List of anticipated file names
	my $mrsfastsam = "$out/$outbase/$bnames.$i.se.mrsfast.bam";
	my $dup = "$out/$outbase/$bnames.$i.clean.sort.nodup.bam";
	my $tmpdir = "$out/$outbase";

	# Organized list of return files
	my @return = ($mrsfastsam, $dup);
	
	my $command = "perl $bin_dir/run_alignment_programs_se.pl -i $file -o $out -r $reference -a $alternate -j $java_dir -v Xmx1G -e $edit_dis -n $i -c $outbase -g $rgnum -b $bin_dir";
	open(OUT, "> $command_dir/S_$bnames\_align.$i.sh") || print "Could not create command output! $!\n";
	print OUT $command;
	close OUT;
	return $command, \@return, $outbase;
}

sub t_format_mrsfast_bam {
	my($self, $folder, $reference, $bin, $command_dir, $bname, $split) = @_;
	my $fai = "$reference.fai";
	my $splits = "$folder/split_read";
	
	my $command = "perl $bin/convert_mrsfast_sams_to_bams.pl -i $folder -r $fai -b $bin;";
	if($split){
		$command .= "perl $bin/convert_mrsfast_sams_to_bams.pl -i $splits -r $fai -b $bin";
	}
	open(OUT, "> $command_dir/B_$bname\_samtools.sh") || print "Could not create command output! $!\n";
	print OUT $command;
	close OUT;
	return $command;
}
1;