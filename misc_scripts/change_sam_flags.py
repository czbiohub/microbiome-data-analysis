import argparse

class process_samfiles:
    """
    Class Description:  Clear or set specified flags in sam files
    User Notes:
    Revision History:   2017.03.31 Brian Yu Created
                        
    """

    
    def __init__(self, inputsam, outputsam, bit, setflag):
        """
        Checks sam flags in a sam file and perform the following:
        1. if a bit is set and needs to be cleared, clear it
        2. if a bit is not set and needs to be set, set it
        bit: 0-10
        setflag: True or False
        """
        if type(bit) is not int:
            raise TypeError('bit is not an integer')

        if type(setflag) is not bool:
            raise TypeError('setflag is not a boolean')

        with open(inputsam, 'r') as fin, open(outputsam, 'w') as fout:

            for l in fin:
                if l.strip() is not '':
                    # if line begins with '@' then it's a header
                    if l[0] == '@':
                        fout.write('%s' %(l))
                    else:
                        templine = l.split() # l should be tab delimited
                        if len(templine) < 9:
                            print(l)
                        samflag = int(templine[1]) # field 2 should be samflag
                        # make a mask
                        bitmask = 1 << bit
                        if setflag:
                            templine[1] = str(samflag | bitmask)
                        else:
                            templine[1] = str(samflag & ~bitmask)
                        fout.write('\t'.join(templine) + '\n')


# When running the script from command line, the following lines are executed
if __name__ == "__main__":
    usage = "USAGE: python change_sam_flags.py set/clear bitnumber input_samfile output_samfile"

    # Making default argument list structures
    p = argparse.ArgumentParser(usage=usage)
    p.add_argument(dest='setflag', action='store', type=str)
    p.add_argument(dest='bitnum', action='store', type=int)
    p.add_argument(dest='inputname', action='store', type=str)
    p.add_argument(dest='outputname', action='store', type=str)

    arguments = p.parse_args()

    try:
        if arguments.setflag == 'set' or arguments.setflag == 'Set' or arguments.setflag == 's':
            process_samfiles(arguments.inputname, arguments.outputname, arguments.bitnum, True)
        elif arguments.setflag == 'clear' or arguments.setflag == 'Clear' or arguments.setflag == 'c':
            process_samfiles(arguments.inputname, arguments.outputname, arguments.bitnum, False)
        else:
            raise TypeError('Input not correct')

    except ValueError, e:
        print "ERROR: ValueError %s" % e
        print usage
    except TypeError, e:
        print "ERROR: TypeError %s" % e
        print usage
    except IOError, e:
        print "ERROR: IOError %s" % e
        print usage

