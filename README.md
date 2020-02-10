# demultiplex
Demultiplexes Illumina reads (Illumina 1.8+) with dual matched indexes. 

Example usage:

python3 demultiplex.py -r1 R1.fastq.gz -r2 R2.fastq.gz -r3 R3.fastq.gz -r4 R4.fastq.gz -q 30

There are 4 required arguments: r1, r2, r3, and r4. 

This assumes that R1 and R4 are forward strand and reverse strand biological reads, respectively. Also assumes that R2 and R3 are your forward strand and reverse strand barcodes, respectively. 

There is 1 optional argument: -q (--qual_thresh.) By default, when any of the biological reads or barcode sequences have an average quality score less than 30, those records are automatically put in an "undetermined" file. 

Current speed is 1 million records (4 million lines)  in approximately 3.5 minutes. 

