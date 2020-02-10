import gzip
import argparse
from index import IndexRecord
from sequence import SequenceRecord

def write_file_dict(index_list):
    '''
    Creates and returns a dictionary which contains information on index and associated output file
    '''

    fw_file_list = []
    rv_file_list = []

    # make files for fw and rv
    for index in index_list:
        fw_file_list.append(open('fw_' + index + '.fastq','a'))
        rv_file_list.append(open('rv_' + index + '.fastq','a'))

    # creates tuple with (forward file wrapper, reverse file wrapper) for each fw/rv combination
    file_list_tup = tuple(zip(fw_file_list,rv_file_list))

    # creates dictionary where the key is the barcode, value is tuple with (forward file open wrapper, reverse file open wrapper)
    output_file_dict = dict(zip(index_list,file_list_tup))

    return(output_file_dict)

def process_files(filenames,threshold, index_list_fn):
    '''
    This function will process the FASTQ files.
    It will check for correct indexes, index hopping, and bad indexes
    It will also write to specific files based on the above
    '''

    seq_1 = filenames[0]
    seq_2 = filenames[3]
    index_1 = filenames[1]
    index_2 = filenames[2]

    # get valid indexes into set
    f = open(index_list_fn)
    index_list = f.readline().strip().split(',')
    f.close()

    # set counters to count number of occurrences for each situation
    index_hop_counter = 0
    correct_read_counter = 0
    und_read_counter = 0

    # get output file dictionary from function
    output_file_dict = write_file_dict(index_list)

    # open hopped files to write to
    fw_hop_output = open('fw_hopped.fastq','a')
    rv_hop_output = open('rv_hopped.fastq','a')
    # open undetermined files to write to
    fw_und_output = open('fw_undetermined.fastq','a')
    rv_und_output = open('rv_undetermined.fastq','a')

    # store amount of records read (to report later)
    record_counter = 0

    # open the four sequence files
    with gzip.open(seq_1,'r') as s1, gzip.open(seq_2,'r') as s2, gzip.open(index_1,'r') as i1, gzip.open(index_2,'r') as i2:
        print('Started demultiplexing...')
        while True:

            seq1_list = [next(s1).decode('UTF-8').strip() for x in range(4)]
            seq2_list = [next(s2).decode('UTF-8').strip() for x in range(4)]
            index1_list = [next(i1).decode('UTF-8').strip() for x in range(4)]
            index2_list = [next(i2).decode('UTF-8').strip() for x in range(4)]

            seq1 = SequenceRecord(seq1_list)
            seq2 = SequenceRecord(seq2_list)
            index1 = IndexRecord(index1_list)
            index2 = IndexRecord(index2_list)

            # if first line of seq1 is empty string, means at end of file
            # break out of loop
            if seq1.header == '':
                break

            record_counter += 1

            if record_counter % 1000000 == 0:
                print(record_counter/4., 'records processed.')

            # add that combined string to header of each sequence (1st line of seq1/seq2 += concatenated index string above)

            seq1.update_header(index1.index_seq, index2.index_seq)
            seq2.update_header(index1.index_seq, index2.index_seq)

            # if the indexes are not in the index list, or the quality score of any of the biological/index sequences are less than the threshold (default is 30), then store the record in an undetermined file.
            if (index1.index_seq not in index_list) or (index2.reverse_complement() not in index_list) or (seq1.average_quality() < threshold) or (seq2.average_quality() < threshold) or (index1.average_quality() < threshold) or (index2.average_quality() < threshold):
                und_read_counter += 1
                fw_und_output.write('\n'.join((seq1.header, seq1.sequence, seq1.optional_line, seq1.quality_line)) + '\n')
                rv_und_output.write('\n'.join((seq2.header, seq2.sequence, seq2.optional_line, seq2.quality_line)) + '\n')
                continue

            # if index_1 same as reverse_complement(index_2), this means NO index hopping
            # write record to barcode-specific output file
            if index1.index_seq == index2.reverse_complement():
                correct_read_counter += 1

                fw_output_file = output_file_dict[index1.index_seq][0]
                rv_output_file = output_file_dict[index2.reverse_complement()][1]

                # write R1 list to file, where filename is based on output file dictionary (key is index_1)
                # write R4 list to file, where filename is based on output file dictionary (key is index_1), except has reverse designation in filename
                fw_output_file.write('\n'.join((seq1.header, seq1.sequence, seq1.optional_line, seq1.quality_line)) + '\n')
                rv_output_file.write('\n'.join((seq2.header, seq2.sequence, seq2.optional_line, seq2.quality_line)) + '\n')
            # otherwise, index hopping has occurred
            else:

                index_hop_counter += 1

                # write R1 list to a index hopping designated file (forward strand)
                # write R4 list to index hopping designated file (reverse strand)
                fw_hop_output.write('\n'.join((seq1.header, seq1.sequence, seq1.optional_line, seq1.quality_line)) + '\n')
                rv_hop_output.write('\n'.join((seq2.header, seq2.sequence, seq2.optional_line, seq2.quality_line)) + '\n')

    fw_hop_output.close()
    rv_hop_output.close()
    fw_und_output.close()
    rv_und_output.close()

    # close all correctly paired index read files
    for key in output_file_dict:
        # gets both fw and rv file wrappers
        values = output_file_dict[key]
        for file_wrapper in values:
            # close both of them
            file_wrapper.close()

    # calculates percent of reads which were correct, undetermined, and index hopped.
    percent_correct = str(100*round(correct_read_counter/record_counter,3))
    percent_index_hop = str(100*round(index_hop_counter/record_counter,3))
    percent_und = str(100*round(und_read_counter/record_counter,3))


    # prints statistics about the demultiplexing process to a text file
    f = open('output_statistics.txt','w')
    print('Finished!')
    print('Total number of reads processed: ' + str(record_counter),file=f)
    print(percent_correct + '% of reads have correctly paired indexes and are above the quality threshold. (' + str(correct_read_counter) + ')',file=f)
    print(percent_index_hop + '% of reads are index hopped. (' + str(index_hop_counter) + ')',file=f)
    print(percent_und + '% of reads are undetermined or below the quality threshold. (' + str(und_read_counter) + ')',file=f)
    print('Statistics written to output_statistics.txt.')
    f.close()

# this is a common line in python scripts
# this evals to TRUE if the file is called directly, FALSE if it is imported
if __name__ == "__main__":
    # add arguments for argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-r1", "--r1", required=True, type = str, help="Enter R1 (read 1) filename.")
    parser.add_argument("-r2", "--r2", required=True, type = str, help="Enter R2 (index 1) filename.    ")
    parser.add_argument("-r3", "--r3", required=True, type = str, help="Enter R3 (index 2) filename.")
    parser.add_argument("-r4", "--r4", required=True, type = str, help="Enter R4 (read 2) filename.")
    parser.add_argument("-q", "--qual_thresh", type=int, default = 30, required=False, help="Enter average quality score threshold.")
    parser.add_argument("-i", "--index_file", type=str, required=True, help="Comma separated list of valid indexes.")
    args = parser.parse_args()

    # place user's arguments into list
    file_list = [args.r1, args.r2, args.r3, args.r4]

    # threshold refers to phred score
    process_files(file_list, args.qual_thresh, args.index_file)
