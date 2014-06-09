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

$ perl run_cnv_caller_pipeline.pl -c [your configuration file] 
