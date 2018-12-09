import argparse, random

class threshold_scaffolds:
    """
    Class Description:
    User Notes:
    Revision History:   2015.04.13 Brian Yu Created
                        2015.05.28 Changed file open syntax to with open() as f
    """

    def __init__(self, inputscaffolds):
        """
        Reads in both fastq files, keeps index. This file could be used for both
        spades output and output of supercontigs so the contig labes might not be >NODE
        :param inputscaffolds: name of the input scaffold file in fasta format
        :return:
        """

        with open(inputscaffolds, 'r') as fid:

            self.label = []
            self.sequence = []

            for l in fid:
                # if line begins with '>NODE'
                if '>' in l:
                    # Assume label is on one line, this shouldn't remove spaces
                    self.label.append(l.strip('\n'))
                else:
                    # if length node > length seq then add to append to seq
                    # Basically this is the first line of sequence
                    if len(self.label) > len(self.sequence):
                        self.sequence.append(l.strip('\n'))
                    # else append to existing seq (this is for multi line fasta.
                    else:
                        self.sequence[-1] = self.sequence[-1] + l.strip('\n')

            assert(len(self.sequence) == len(self.label))
            input_scaffold_length = len(self.label)

        # print "Done Init"


    def pick_long_scaffolds(self, threshold, outputname):
        """
        Outputs a new scaffold file containing only scaffolds longer than threshold
        :param threshold: length threshold in bp
               outputname: name of the output scaffold file in fasta format
        :return:
        """

        with open(outputname, 'w') as fid:
        
            # if you add something to the beginning to the contig id then length is no longer #3
            for i in range(len(self.label)):
                # if int(self.label[i].split('_')[3]) >= threshold:
                if len(self.sequence[i]) >= threshold:
                    fid.write('%s\n%s\n' %(self.label[i], self.sequence[i]))


# When running the script from command line, the following lines are executed
if __name__ == "__main__":
    usage = "USAGE: python threshold_scaffolds threshold input_scaffolds output_scaffolds"

    # Making default argument list structures
    p = argparse.ArgumentParser(usage=usage)
    p.add_argument(dest='thresh', action='store', type=int)
    p.add_argument(dest='inputname', action='store', type=str)
    p.add_argument(dest='outputname', action='store', type=str)

    arguments = p.parse_args()

    try:
        f = threshold_scaffolds(arguments.inputname)
        f.pick_long_scaffolds(arguments.thresh, arguments.outputname)

    except ValueError, e:
        print "ERROR: ValueError %s" % e
        print usage
    except TypeError, e:
        print "ERROR: TypeError %s" % e
        print usage
    except IOError, e:
        print "ERROR: IOError %s" % e
        print usage

