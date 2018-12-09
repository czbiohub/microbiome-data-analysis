import argparse, pysam
import numpy as np

class bamProcessingTools:
    """
    Class Description:  A set of tools for processing bam files
                        -- Currently only handles BAM files.
                        -- Assumes bam is not sorted
    Revision History:   2017.09.28 Brian Yu Created
    """

    def __init__(self, bam_file):
        """
        Initialize by creating an AlignmentFile structure for the bam file.
        Assumes bam is NOT sorted. So all the name and coord are read in here.
        """
        # Created pysam AlignmentFile structure for bam
        self.input_bamfile = bam_file
        bamfile = pysam.AlignmentFile(self.input_bamfile,'rb')
        # Read in name and coord, do not sort here
        self.query_names = {}
        for read in bamfile.fetch(until_eof=True):
            if read.query_name in self.query_names.keys():
                self.query_names[read.query_name] = self.query_names[read.query_name] + 1
            else:
                self.query_names[read.query_name] = 1
        bamfile.close()
        print('Number of bam names with only one read is:', len([1 for x in self.query_names.keys() if self.query_names[x] == 1]), sep=' ')

        
    def remove_unpaired_alignment(self, outputfile):
        """
        After quality filtering of bam alignments, expecially with bwa, some
        paired alignments may have only 1 alignment left. This function deletes
        those alignments from the bam files.
        -- Assumes bam is not sorted
        -- Read the bam file and if 2 record for query then print
        """
        # create a bam read and write structure
        bamfile = pysam.AlignmentFile(self.input_bamfile, "rb")
        output_bamfile = pysam.AlignmentFile(outputfile, "wb", template=bamfile)
        count = 0
        # process bam reads
        for read in bamfile.fetch(until_eof=True):
            if self.query_names[read.query_name] == 2:
                output_bamfile.write(read)
                count += 1
        output_bamfile.close()
        bamfile.close()
        print('Number of bam records in output bamfile is:', count, sep=' ')

        
# When running the script from command line, the following lines are executed
if __name__ == "__main__":
    usage = "USAGE: python bamProcessingTools.py input_file output_file"

    # Making default argument list structures
    p = argparse.ArgumentParser(usage=usage)
    p.add_argument(dest='input_file_name', action='store', type=str)
    p.add_argument(dest='output_file_name', action='store', type=str)

    arguments = p.parse_args()

    try:
        f = bamProcessingTools(arguments.input_file_name)
        f.remove_unpaired_alignment(arguments.output_file_name)
    except ValueError as e:
        print("ERROR: ValueError:",e)
        print(usage)
    except TypeError as e:
        print("ERROR: TypeError:",e)
        print(usage)
    except IOError as e:
        print("ERROR: IOError %s",e)
        print(usage)

