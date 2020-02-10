

class SequenceRecord:

    def __init__(self, record):
        '''
        Each sequence record contains a header, sequence, optional line, and quality line (phred scores)
        '''
        self.record = record
        self.header = record[0]
        self.sequence = record[1]
        self.optional_line = record[2]
        self.quality_line = record[3]


    def update_header(self, index1, index2):
        '''
        Updates header of the sequence record
        '''
        self.header = self.header +  ' ' + index1 + '-' + index2

    def average_quality(self, phred_offset):
        '''
        Takes in a sequence of nucleotides, and returns the average quality score
        '''
        # for each character, get phred score using ord(), and adds it to cumulative sum variable

        sum_quality_score = sum(ord(char) - phred_offset for char in self.quality_line)
        #for char in self.quality_line:
            #phred_score = ord(char)-33
            #sum_quality_score += phred_score

        # gets avg by dividing sum by length of sequence
        avg_quality_score = sum_quality_score / len(self.quality_line)
        return(avg_quality_score)
