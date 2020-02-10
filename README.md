# demultiplex
Demultiplexes Illumina reads with dual matched indexes. 

Example usage:

python3 demultiplex.py -r1 R1.fastq.gz -r2 R2.fastq.gz -r3 R3.fastq.gz -r4 R4.fastq.gz -q 30 -i index_file.csv -f 1.9

There are 6 required arguments: r1, r2, r3, r4, index_file, and format. 

The first four are the compressed fastq files. The code assumes that R1 and R4 are forward strand and reverse strand biological reads, respectively. Also assumes that R2 and R3 are your forward strand and reverse strand barcodes, respectively. 

index_file is a txt file which contains each valid index, separated by commas. 

f (--format) is the Illumina real-time analysis version. This is for quality score adjustment. Illumina 1.8+ is Phred+33, below that is Phred+64. 

There are 1 optional argument: -q (--qual_thresh.) By default, when any of the biological reads or barcode sequences have an average quality score less than 30, those records are automatically put in an "undetermined" file. 

Current speed is 1 million records (4 million lines)  in approximately 3.5 minutes. 

