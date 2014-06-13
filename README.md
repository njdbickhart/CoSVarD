CoSVarD
=======

CoSVarD stands for "Complementary Sequence Variant Detection" and it is a pipline I have developed for processing sequence data for SNP, INDEL and CNV calling. The current pipeline is in Beta, and I am working on a replacement that is far more extensible. 

The current features include:
  - The ability to input simple spreadsheets to automate processing of your data
  - Integration of SNP and INDEL calling, as well as CNV calling from popular and custom methodologies
  - LSF and simple threading
  - Scalability through a "divide and conquer" algorithm

Please use the "sample_config_file.txt" configuration file in the "lib" folder as an example input for the program. Any line that begins with a '#' is a comment and is not interpretted by the configuration file parser.

To run the pipeline, you simply invoke it as follows:

$ perl run_cnv_caller_pipeline.pl -c [your configuration file] -p [The max number of threads you wish to use]


Here is a detailed list of instructions for installing and using CoSVarD:

1) First, setting up your environment
	a) Your first priority is to ensure that you have the proper Perl Modules and setup.
		- Use CPAN to install the following modules (For a brief tutorial on CPAN, please see: http://perl.about.com/od/packagesmodules/qt/perlcpan.htm)
		- Forks::Super (and Sys::CpuAffinity/Sys::CpuLoadX if you can!)
		- namespace::autoclean
		- FileHandle
		- File::Basename
		- Cwd
	b) Next, you need the following external programs and resources:
		- Samtools (version > 0.1.19)
		- BWA (version > 0.6.2)
		- PicardTools (version > 1.53)
		- The sam JDK toolkit (from the makers of PicardTools) (version > 1.83)
		- MrsFAST (version > 2.3.0.2)
		- The Java Development Kit (JDK) version 1.7 or 1.8 (gives java version 7 and 8, respectively)
		- My PairMatchMrsfastSam.jar program (found also on GitHub)
		- My BedUtils.jar library (found also on GitHub)
	c) Now, to piece it all all together into one coherent package, create one "bin" folder that you will use for the pipeline (it can be your ~/bin folder if need be!)
	d) Place the CoSVarD scripts into this bin folder directly, and ensure that the "lib" folder is within your bin folder!
	e) Place the following additional programs into your "bin" folder:
		- java
		- mrsfast
		- samtools
		- bwa
		- PairMatchMrsfastSam.jar
		- CleanSam.jar (Picardtools)
		- ReorderSam.jar (Picardtools)
		- MarkDuplicates.jar (Picardtools)
	f) Now, place the following additional libraries into your "bin/lib" folder:
		- BedUtils.jar
		- The "sam-1.83.jar" from the Sam JDK package
	e) Just to ensure that everything is setup correctly, run "perl run_cnv_caller_pipeline.pl" within your bin folder. If you get a "usage" text blurb, then you should have most of everything in place!
2) Now, you have to prepare your reference fasta files for alignment
	a) CoSVarD uses two reference fasta files in order to prepare files for CNV and SNP/INDEL calling
		- One file is a masked reference fasta file, generated either from repeatmasker output, or downloaded from the UCSC genome browser ("base_ref_genome" reference fasta in the config file)
		- The other file is the unmasked reference fasta ("alt_ref_genome" reference fasta in the config file)
	b) Repeat the following steps on both fasta files
	c) First, you must index your reference fasta using both the BWA and MrsFAST alignment programs:
		- bwa index reference.fa
		- mrsfast --index reference.fa
	d) Now, its time to create a fasta index using the samtools and Picard programs
		- samtools faidx reference.fa
		- java -jar picard-tools/CreateSequenceDictionary.jar REFERENCE=reference.fa OUTPUT=reference.dict
3) Next, setting up your configuration file
	a) In order to give you extensive control over your run, CoSVarD uses a configuration file system to parse your command options and input data
	b) Use the "bin/lib/sample_config_file.txt" example to help you start setting up your own config file
	c) **More Soon!**
4) Finally, the output data that you will receive:
	a) Within each animal's folder you will have the following files:
		- File type: *.mrsfast.bam is the MrsFAST alignment file template
		- File type: *.clean.sort.nodup.bam is the BWA alignment file template
		- Directory: split_read
			- Within this folder you should have:
			- File type: *.single.txt.gz is the one-end-anchor file template
			- File type: *.working.mrsfast.sam.gz is the split-read alignment file template
		- Directory: divet
			- Within this folder there should only be one file type
			- File type: *.divet.gz is the discordant read pair file template
	b) All other files are extraneous and can be deleted if space is at a premium.