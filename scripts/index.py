class IndexRecord:

    def __init__(self, record):
        '''
        Each barcode record contains a header, sequence, optional line, and quality line (phred scores)
        '''
        self.record = record
        self.header = record[0]
        self.index_seq = record[1]
        self.optional_line = record[2]
        self.quality_line = record[3]
        self.complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}

    def reverse_complement(self):

        """
        Returns the reverse complement of the barcode
        """

        # [::-1] slice reverses a string
        reverse_barcode = self.index_seq[::-1]

        # convert string into list of chars
        reverse_barcode_list = list(reverse_barcode)

        # makes list of reverse complement, combines back into string using join function
        reverse_complement = ''.join([self.complement[i] for i in reverse_barcode_list])

        return(reverse_complement)

    def average_quality(self, phred_offset):
        '''
        Takes in a sequence of nucleotides, and returns the average quality score
        '''
        # for each character, get phred score using ord(), and adds it to cumulative sum variable
        #print(ord(char) - 33 for char in self.quality_line)
        sum_quality_score = sum(ord(char) - phred_offset for char in self.quality_line)

        #for char in self.quality_line:
            #phred_score = ord(char)-33
            #sum_quality_score += phred_score

        # gets avg by dividing sum by length of sequence
        avg_quality_score = sum_quality_score / len(self.quality_line)
        return(avg_quality_score)
