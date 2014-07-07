#!/usr/bin/perl
# These are subroutines separated from the original pipeline to unclutter it
package lsf_pipeline_utils;
use warnings;
use strict;

sub new {
	my ($class,%arg) = @_;
	my $self={};
	bless($self,$class || ref($class));
	return $self;
}

sub samtools_merger {
	require strict;

	my $samtools = "/ibfs7/asg2/bickhartd/bin/samtools";

	my $prefix = shift(@_);
	my $num = shift(@_);
	my @mergefiles = @_;
	chomp(@mergefiles);
	chomp $num;
	chomp $prefix;

	my $mergestr = join(" ", @mergefiles); 
	my $output =  "$prefix.$num.merge.sorted.bam";
	system(qq{bsub -J merge$num -oo $prefix.$num.merge.OUT -K -R "rusage[mem=4000]" $samtools merge $output $mergestr });
}
sub final_total_merger_command {
	my ($self, $doc, $vh, $sr, $out, $basename, $cnfile) = @_;
	my $command_dir = "/ibfs7/asg0/bickhartd/commands";
	
	my $command = "bsub -J FM_$basename -K ";
	$command .= "-oo $command_dir/out/FM_$basename.out ";
	$command .= "-e $command_dir/err/FM_$basename.err ";
	$command .= qq{-R "rusage[mem=20000]" };
	$command .= "-q normal ";

	$command .= "/ibfs7/asg2/bickhartd/bin/merge_final_cnv_calls.pl -d $doc -p $vh -s $sr -c $cnfile -o $out";
	return $command;
}

sub format_split_command{
	my ($self, $input, $output, $basename) = @_;
	my $command_dir = "/ibfs7/asg0/bickhartd/commands";
	
	my $command = "bsub -J RS_$basename -K ";
	$command .= "-oo $command_dir/out/RS_$basename.out ";
	$command .= "-e $command_dir/err/RS_$basename.err ";
	$command .= qq{-R "rusage[mem=5000]" };
	$command .= "-q normal ";

	$command .= "/ibfs7/asg2/bickhartd/bin/merge_split_read_output.pl $input $output";
	return $command;
}

sub format_doc_command {
	my ($self, $autodups, $autodels, $xdups, $xdels, $cnfile, $finaldocout, $basename) = @_;
	my $command_dir = "/ibfs7/asg0/bickhartd/commands";
	
	my $command = "bsub -J FD_$basename -K ";
	$command .= "-oo $command_dir/out/FD_$basename.out ";
	$command .= "-e $command_dir/err/FD_$basename.err ";
	$command .= qq{-R "rusage[mem=5000]" };
	$command .= "-q normal ";

	$command .= "/ibfs7/asg2/bickhartd/bin/merge_doc_output.pl -i $autodups ";
	$command .= "-d $autodels -x $xdups -y $xdels -c $cnfile -o $finaldocout";
	return $command; 
}

sub vh_merge_command {
	my ($self, $flistref, $output, $basename) = @_;
	my $command_dir = "/ibfs7/asg0/bickhartd/commands";
	
	my $command = "bsub -J VM_$basename -K ";
	$command .= "-oo $command_dir/out/VM_$basename.out ";
	$command .= "-e $command_dir/err/VM_$basename.err ";
	$command .= qq{-R "rusage[mem=5000]" };
	$command .= "-q normal ";

	$command .= "/ibfs7/asg2/bickhartd/bin/merge_vh_output.pl ";
	$command .= join(" ", @{$flistref});
	$command .= " $output";
	return $command; 
}

sub finalize_split_command {
	my ($self, $psplitfolder , $outfile, $basename) = @_;
	my $command_dir = "/ibfs7/asg0/bickhartd/commands";
	
	my $command = "bsub -J FS_$basename -K ";
	$command .= "-oo $command_dir/out/FS_$basename.out ";
	$command .= "-e $command_dir/err/FS_$basename.err ";
	$command .= "-n 4 ";
	$command .= qq{-R "rusage[mem=9500] span[hosts=1]" };
	$command .= "-q normal ";
	
	$command .= "/ibfs7/asg2/bickhartd/bin/finalize_split_read.pl -i $psplitfolder -o $outfile -b $basename";
	return $command;	
}

sub generate_vh_command {
	my ($self, $pvhfolder, $insert_size, $insert_sd, $chr_length, $gap_file, $basename, $itname) = @_;
	my $command_dir = "/ibfs7/asg0/bickhartd/commands";
	
	my $command = "bsub -J VH_$basename -K ";
	$command .= "-oo $command_dir/out/VH_$basename.out ";
	$command .= "-e $command_dir/err/VH_$basename.err ";
	$command .= qq{-R "rusage[mem=5000]" };
	$command .= "-q normal ";
	
	my $full_sv_events = "$pvhfolder/$basename.$itname.final.SV";
	my $max = $insert_size + (3 * $insert_sd);
	my $min = $insert_size - (2 * $insert_sd);
	if($min < 100){
		$min = 100;
	}
	$command .= "'/ibfs7/asg2/bickhartd/bin/run_variation_hunter_pipeline.pl -i $pvhfolder ";
	$command .= "-n $basename -m $min -M $max -p 0.001 -t 4 -c $chr_length -g $gap_file -I $itname'";
	return ($command, $full_sv_events);
}


sub generate_wssd_command {
	my ($self, $folderpath, $outprefix, $reference, $basename) = @_;
	my $command_dir = "/ibfs7/asg0/bickhartd/commands";

	# predicted files
	my $file1 = "$folderpath/doc_wins/$outprefix.file1.bed";
	my $file2 = "$folderpath/doc_wins/$outprefix.file2.bed";
	my $file3 = "$folderpath/doc_wins/$outprefix.file3.bed";
	my $autodups = "$folderpath/doc_wins/$outprefix.file1.bed.final.wssd";
	my $autodels = "$folderpath/doc_wins/$outprefix.file1.bed.final.deletions.tab";
	my $xdups = "$folderpath/doc_wins/$outprefix.file1.bed.X.final.wssd";
	my $xdels = "$folderpath/doc_wins/$outprefix.file1.bed.final.X.deletions.tab";
	my $cnfile = "$folderpath/doc_wins/$outprefix.file3.bed.gc.depth.normalized.CN";

	my $command = "bsub -J W_$basename -K ";
	$command .= "-oo $command_dir/out/J_$basename.out ";
	$command .= "-e $command_dir/err/J_$basename.err ";
	$command .= qq{-R "rusage[mem=4000] span[hosts=1]" };
	$command .= "-n 8 ";
	$command .= "-q normal ";
	
	# Just because I'm daring, I want to incorporate the WSSD pipeline call in the same comand
	# This simplifies the calling statement and queue scheme that I have for processing all the data
	$command .= "\'/ibfs7/asg2/bickhartd/bin/alkan_pipeline_no_controls_placeholder.pl ";
	$command .= "--File1 $file1 --File2 $file2 --File3 $file3 --bed $outprefix --hits $folderpath --ref $reference\'";
	
	return $command, $autodups, $autodels, $xdups, $xdels, $cnfile;
}

sub split_pe_fastqs{
	my ($self, $fq1, $fq2, $out, $linelimit, $sit) = @_;
	my $f1 = FileHandle->new;
	my $f2 = FileHandle->new;
	if($fq1 =~ /\.gz$/){
		$f1->open("gunzip -c $fq1 |") || print "Gzip file open error! $!\n";
	}else{
		$f1->open("< $fq1") || print "File open error! $!\n";
	}
	if($fq2 =~ /\.gz$/){
		$f2->open("gunzip -c $fq2 |") || print "Gzip file open error! $!\n";
	}else{
		$f2->open("< $fq2") || print "File open error! $!\n";
	}
	my @f1base = split(/\//, $fq1);
	my @f2base = split(/\//, $fq2);

	open(O1, "> $out/$f1base[-1].ptmp_$sit.fq1") || print "File output error! $!\n";
	open(O2, "> $out/$f2base[-1].ptmp_$sit.fq2") || print "File output error! $1\n";


	my $fcount = $sit + 1;#changed $sit++ to $sit + 1
	my $linecount = 0;
	my (@filelist1, @filelist2); 	# list of returned chunk fastqs
	push(@filelist1, "$out/$f1base[-1].ptmp_$sit.fq1");
	push(@filelist2, "$out/$f2base[-1].ptmp_$sit.fq2");
	while (my $h1 = <$f1>){
		my $h2 = <$f2>;
		my $s1 = <$f1>;
		my $s2 = <$f2>;
		my $p1 = <$f1>;
		my $p2 = <$f2>;
		my $q1 = <$f1>;
		my $q2 = <$f2>;
		$h1 =~ s/[_\/]/-/g;
		$h2 =~ s/[_\/]/-/g;
		chomp($h1, $h2, $s1, $s2, $p1, $p2, $q1, $q2);
		
		# remove illumina filter tags. Check for parity
		my @h1seg = split(/\s+/, $h1);
		my @h2seg = split(/\s+/, $h2);
		if (scalar(@h1seg) > 1 && scalar(@h2seg) > 1){
			if(($h1seg[1] =~ m/.:Y:./) || ($h2seg[1] =~ m/.:Y:./) || length($s1) < 35 || length($s2) < 35){
				next;
			}
		}else{
			if(length($s1) < 35 || length($s2) < 35){next;}
			@h1seg = split(/\//, $h1);
			@h2seg = split(/\//, $h2);
		}

		# Now, split reads into 50 bp chunks 
		my (@s1ch, @q1ch, @s2ch, @q2ch);
		for (my $x = 0; $x < length($s1); $x += 50){
			my $tmps1 = substr($s1, $x, 50);
			my $qmps1 = substr($q1, $x, 50);
			if(length($tmps1) != 50){next;}
			elsif(length($qmps1) != 50){next;}
			push(@s1ch, $tmps1);
			push(@q1ch, $qmps1);
		}
		for (my $x = length($s2); $x - 50 > 0 ; $x -= 50){
			my $tmps2 = substr($s2, $x - 50, 50);
			my $qmps2 = substr($q2, $x - 50, 50);
			if(length($tmps2) != 50){next;}
			elsif(length($qmps2) != 50){next;}
			push(@s2ch, $tmps2);
			push(@q2ch, $qmps2);
		}

		for (my $x = 0; $x < scalar(@s2ch) && $x < scalar(@s1ch); $x++){
			print O1 "$h1seg[0]:$x\_1\n$s1ch[$x]\n$p1\n$q1ch[$x]\n";
			print O2 "$h2seg[0]:$x\_2\n$s2ch[$x]\n$p2\n$q2ch[$x]\n";
			$linecount++;
		}
		
		if($linecount != 0 && $linecount >= $linelimit){
			close O1;
			close O2;
			$linecount = 0;
			print "Finished with file:\t$fcount\r";
			if(-e "$out/$f1base[-1].ptmp_$fcount"){
			
			}
			open(O1, "> $out/$f1base[-1].ptmp_$fcount.fq1") || print "File output error! $!\n";
			open(O2, "> $out/$f2base[-1].ptmp_$fcount.fq2") || print "File output error! $1\n";
			push(@filelist1, "$out/$f1base[-1].ptmp_$fcount.fq1");
			push(@filelist2, "$out/$f2base[-1].ptmp_$fcount.fq2");
			$fcount++;
		}
	}
	print "\n";
	$f1->close();
	$f2->close();
	close O1;
	close O2;
	return $fcount, \@filelist1, \@filelist2;
}

sub _skip_pe_split_fq{

}

sub split_se_fastqs {
	my ($self, $f, $out, $linelimit, $sit) = @_;
	my $in = FileHandle->new;
	if($f =~ /\.gz$/){
		$in->open("gunzip -c $f |") || print "Gzip file open error! $!\n";
	}else{
		$in->open("< $f") || print "File open error! $!\n";
	}
	my @fbase = split(/\//, $f);
	open(OUT, "> $out/$fbase[-1].stmp_$sit.fq") || print "File output error! $!\n";
	
	my @filelist;
	push(@filelist, "$out/$fbase[-1].stmp_$sit.fq");
	my $fcount = $sit + 1;#changed $sit++ to $sit + 1
	my $linecount = 0;
	while(my $h = <$in>){
		my $s = <$in>;
		my $p = <$in>;
		my $q = <$in>;
		$h =~ s/[_\/]/-/g;
		chomp($h, $s, $p, $q);
		my @hseg = split(/\s+/, $h);
		if(scalar(@hseg) > 1){
			if ($hseg[1] =~ m/.:Y:./ || length($s) < 35){
				next;
			}
		}else{
			if(length($s) < 35){next;}
			@hseg = split(/\//, $h);
		}
		my (@s1ch, @q1ch, @s2ch, @q2ch);
		for (my $x = 0; $x < length($s); $x += 35){
			my $tmps1 = substr($s, $x, 35);
			my $tmps2 = substr($q, $x, 35);
			if(length($tmps1) != 35){next;}
			elsif(length($tmps2) != 35){next;}
			push(@s1ch, $tmps1);
			push(@q1ch, $tmps2);
		}
		for (my $x = 0; $x < scalar(@s1ch); $x++){
			print OUT "$hseg[0]:$x\_1\n$s1ch[$x]\n$p\n$q1ch[$x]\n";
			$linecount++;
		}
		
		if($linecount != 0 && $linecount >= 500000){
			$linecount = 0;
			close OUT;
			open(OUT, "> $out/$fbase[-1].stmp_$fcount.fq") || print "File output error! $!\n";
			push(@filelist, "$out/$fbase[-1].stmp_$fcount.fq");
			$fcount++;
		}
	}
	$in->close();
	close OUT;
	return $fcount, \@filelist;#changed @filelist to \@filelist
}
sub generate_palign_command {
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
	# Since the divets are flat files, and can be easily "cat" into large files manually, I will not include them
	my @return = ($mrsfastsplitsam1, "", $divet, $single, $dup);
			
	my $l = "50";
	my $m = "500000";
	my $u = $insert + (3 * $sd);
	my $x = "50";
	
	my $command = "bsub -J P_$bnames\_align.$i -K ";
	$command .= "-oo $command_dir/out/P_$bnames\_align.$i.out ";
	$command .= "-e $command_dir/err/P_$bnames\_align.$i.err ";
	$command .= qq{-R "rusage[mem=1000]" };
	$command .= "-q normal ";

	$command .= "\'perl $bin_dir/run_alignment_programs_pe.pl -1 $f1 -2 $f2 -o $out -r $reference -e $edit_dis ";
	$command .= "-l $l -m $m -u $u -n $i -a $alternate -c $outbase -b $bin_dir -j $java_dir -v Xmx1G -g $rgnum\'";
	
	open(OUT, "> $command_dir/P_$bnames\_align.$i.sh") || print "Could not create command output! $!\n";
	print OUT $command;
	close OUT;
	
	# Returning bsub command and list of anticipated files for merger
	return $command, \@return, $outbase;
}


sub generate_salign_command {
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

	my $command = "bsub -J S_$bnames\_align.$i -K ";
	$command .= "-o $command_dir/out/S_$bnames\_align.$i.out ";
	$command .= "-e $command_dir/err/S_$bnames\_align.$i.err ";
	$command .= qq{-R "rusage[mem=1000]" };
	$command .= "-q normal ";

	$command .= "\'perl $bin_dir/run_alignment_programs_se.pl -i $file -o $out -r $reference -a $alternate -j $java_dir -v Xmx1G -e $edit_dis -n $i -c $outbase -g $rgnum -b $bin_dir\'";
	open(OUT, "> $command_dir/S_$bnames\_align.$i.sh") || print "Could not create command output! $!\n";
	print OUT $command;
	close OUT;
	return $command, \@return, $outbase;
}

sub format_mrsfast_bam {
	my($self, $folder, $reference, $bin, $command_dir, $bname, $split) = @_;
	my $fai = "$reference.fai";
	my $splits = "$folder/split_read";
	
	my $command = "bsub -J B_$bname\_samtools -K ";
	$command .= "-o $command_dir/out/B_$bname\_samtools.out ";
	$command .= "-e $command_dir/err/B_$bname\_samtools.err ";
	$command .= qq{-R "rusage[mem=1000]" };
	$command .= "-q normal ";
	
	$command .= "perl $bin/convert_mrsfast_sams_to_bams.pl -i $folder -r $fai -b $bin;";
	if($split){
		$command .= "perl $bin/convert_mrsfast_sams_to_bams.pl -i $splits -r $fai -b $bin";
	}
	open(OUT, "> $command_dir/B_$bname\_samtools.sh") || print "Could not create command output! $!\n";
	print OUT $command;
	close OUT;
	return $command;
}

1;
