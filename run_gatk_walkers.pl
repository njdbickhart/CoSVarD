#!/usr/bin/perl
# This script is designed to run the GATK on a sample given an input merged bam file

use strict;
use Getopt::Std;

my %opts;

my $usage = "$0 -i input_merged.bam -o output_base_name -b bin_dir -p max_processors -r reference_genome.fa -f <BOOLEAN> is threaded?
-j java_executable -g gatk_jar -t training_set_files  -s java_opts -u temp_dir\n";

getopt('iobpfrjgtus', \%opts);

unless(defined($opts{i}) && defined($opts{o}) && defined($opts{b}) && defined($opts{p})){
	print $usage;
	exit;
}

my $recalflag = 0;
my $java = $opts{j};
my $gatk = $opts{g};
my $ref = $opts{r};
my $jopts = $opts{s};
my $proc = (defined($opts{p}))? $opts{p} : 1;
my $temp = (defined($opts{u}))? $opts{u} : "/tmp";

if(defined($opts{u})){
	mkdir("$opts{u}") || print "$!\n";
}

# Checking to see if opts{i} has comma separated values
my $input = $opts{i};
if($opts{i} =~ m/,/){
	my @temp = split(/,/, $opts{i});
	$input = join(" -I ", @temp); 
}

# GATK temp files
my $grp_file = "$opts{o}.recal.grp";
my $targ_file = "$opts{o}.targets.intervals";
my $realigned_bam = "$opts{o}.realigned.bam";
my $reduced_bam = "$opts{o}.realigned.reduced.bam";
my $raw_unfiltered = "$opts{o}.raw.unfiltered.vcf";
my $tranches = "$opts{o}.filt.tranches";
my $recal = "$opts{o}.filt.recal";
my $rplot = "$opts{o}.plots.R";
my $recalibrated_snps = "$opts{o}.recalibrated.unfiltered.vcf";
my $rfiltered_snps = "$opts{o}.recalibrated.filtered.vcf";
my $filtered_snps = "$opts{o}.raw.filtered.vcf";

# Creating base recalibration data 
#$ ~/jdk1.7.0/bin/java -jar ~/GenomeAnalysisTK-2.3-3-g4706074/GenomeAnalysisTK.jar -R ~/reference/umd3_kary_unmask_ngap.fa -T BaseRecalibrator -I BIBR07.clean.sort.nodup.merge.bam -knownSites ~/reference/dbsnp_filtered_umd3_snps.vcf -o BIBR07_recal_data.grp
my @recalentries;
if(defined($opts{t})){
	$recalflag = 1;
	my $recalsites;
	open(IN, "< $opts{t}") || die "Could not open $opts{t} training file!\n";
	while(my $line = <IN>){
		chomp $line;
		if($line =~ /-knownSites/){
			$recalsites = $line;
		}else{
			push(@recalentries, $line);
		}
	}
	close IN;
	print "$java $jopts -Djava.io.tmpdir=$temp -jar $gatk -T BaseRecalibrator -R $ref -I $input $recalsites -nct $proc -o $grp_file";
	system("$java $jopts -Djava.io.tmpdir=$temp -jar $gatk -T BaseRecalibrator -R $ref -I $input $recalsites -nct $proc -o $grp_file");
}

# Realigning indels
#$ ~/jdk1.7.0/bin/java -jar ~/GenomeAnalysisTK-2.3-3-g4706074/GenomeAnalysisTK.jar -R ~/reference/umd3_kary_unmask_ngap.fa -T RealignerTargetCreator -nt 8 -I BIBR07.clean.sort.nodup.merge.bam -o BIBR07_indels.targets.intervals
#$ ~/jdk1.7.0/bin/java -jar ~/GenomeAnalysisTK-2.3-3-g4706074/GenomeAnalysisTK.jar -R ~/reference/umd3_kary_unmask_ngap.fa -T IndelRealigner -I BIBR07.clean.sort.nodup.merge.bam -targetIntervals BIBR07_indels.targets.intervals -o BIBR07.clean.sort.nodup.merge.realign.bam
system("$java $jopts -Djava.io.tmpdir=$temp -jar $gatk -T RealignerTargetCreator -R $ref -nt $proc -I $input -o $targ_file");
system("$java $jopts -Djava.io.tmpdir=$temp -jar $gatk -T IndelRealigner -R $ref -I $input -targetIntervals $targ_file -o $realigned_bam");

# Read reduction to reduce size of bam
#$ ~/jdk1.7.0/bin/java -jar ~/GenomeAnalysisTK-2.3-3-g4706074/GenomeAnalysisTK.jar -R ~/reference/umd3_kary_unmask_ngap.fa -T ReduceReads -I BIBR07.clean.sort.nodup.merge.realign.bam -o BIBR07.clean.sort.nodup.merge.realign.reduced.bam
system("$java $jopts -Djava.io.tmpdir=$temp -jar $gatk -T ReduceReads -R $ref -I $realigned_bam -o $reduced_bam");
if(-s $reduced_bam){
	# Check for existence of reduced bam; if present, delete the realigned bam
	system("rm $realigned_bam");
}else{
	# Contingency: keep the realigned bam if it exists and use it for later
	if(-s $realigned_bam){
		$reduced_bam = $realigned_bam;
	}else{
		print "$0: ERROR problem with $realigned_bam !\n";
		exit;
	}
}
	
# Run the unified Genotyper		
#$ ~/jdk1.7.0/bin/java -jar ~/GenomeAnalysisTK-2.3-3-g4706074/GenomeAnalysisTK.jar -R ~/reference/umd3_kary_unmask_ngap.fa -T UnifiedGenotyper -nt 8 -I BIBR07.clean.sort.nodup.merge.realign.reduced.bam -o BIBR07.initial.raw.snps.vcf
system("$java $jopts -Djava.io.tmpdir=$temp -jar $gatk -T UnifiedGenotyper -R $ref -nt $proc -I $reduced_bam -o $raw_unfiltered");	
			
# if training sets are provided, do recalibration
#$ ~/jdk1.7.0/bin/java -jar ~/GenomeAnalysisTK-2.3-3-g4706074/GenomeAnalysisTK.jar -R ~/reference/umd3_kary_unmask_ngap.fa -T VariantRecalibrator -nt 8 -input BIBR07.initial.raw.snps.vcf -resource:dbsnp,known=true,training=true,truth=false,prior=6.0 ~/reference/dbsnp_filtered_umd3_snps.vcf -resource:ldchip,known=true,training=true,truth=true,prior=12.0 ~/reference/dbsnp_ld_filtered_sorted_known_snps.vcf -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ -mode SNP -recalFile BIBR07_variant.recal -tranchesFile BIBR07_variant.tranches -rscriptFile BIBR07_variant.plots.R
#$ ~/jdk1.7.0/bin/java -jar ~/GenomeAnalysisTK-2.3-3-g4706074/GenomeAnalysisTK.jar -R ~/reference/umd3_kary_unmask_ngap.fa -T ApplyRecalibration -nt 8 -input BIBR07.initial.raw.snps.vcf --ts_filter_level 99.0 -mode SNP -recalFile BIBR07_variant.recal -tranchesFile BIBR07_variant.tranches -o BIBR07.recalibrated.unfiltered.snps.vcf
if($recalflag){
	my $resource = join(' ', @recalentries);
	system("$java $jopts -Djava.io.tmpdir=$temp -jar $gatk -R $ref -T VariantRecalibrator -nt $proc -input $raw_unfiltered $resource -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ -mode SNP -recalFile $recal -tranchesFile $tranches -rscriptFile $rplot");
	system("$java $jopts -Djava.io.tmpdir=$temp -jar $gatk -R $ref -T ApplyRecalibration -nt $proc -input $raw_unfiltered --ts_filter_level 99.0 -mode SNP -recalFile $recal -tranchesFile $tranches -o $recalibrated_snps");
	
	system("$java $jopts -Djava.io.tmpdir=$temp -jar $gatk -R $ref -T VariantFiltration --variant $recalibrated_snps --filterExpression \'QD < 1.0 || MQ < 40.0 || FS > 60.0 || HaplotypeScore > 13.0 || ReadPosRankSum < -8.0\' --filterName \'GATKBest\' --filterExpression \'DP < 3\' --filterName \'DefaultRD\' --filterExpression \'DP > 6 && FS < 0.01\' --filterName \'StrandBias\' -o $rfiltered_snps");
}else{
	system("$java $jopts -Djava.io.tmpdir=$temp -jar $gatk -R $ref -T VariantFiltration --variant $raw_unfiltered --filterExpression \'QD < 1.0 || MQ < 40.0 || FS > 60.0 || HaplotypeScore > 13.0 || ReadPosRankSum < -8.0\' --filterName \'GATKBest\' --filterExpression \'DP < 3\' --filterName \'DefaultRD\' --filterExpression \'DP > 6 && FS < 0.01\' --filterName \'StrandBias\' -o $filtered_snps");
}

exit;