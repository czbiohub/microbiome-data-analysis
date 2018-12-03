import argparse, os, pysam
import numpy as np
import pandas as pd
from subprocess import call
from multiprocessing import Process, Queue
from time import sleep


def process_one_reference(input_parameter, ref_name, sample_list): 
    """
    Tabulates all the wrong unique alignment on one genome over a certain threshold.
    """
    output_filename = input_parameter.folder + '/' + ref_name + '.excluded_locations.txt'
    # Create a dictionary to store positions. keys are contig, position tuples
    temp_contigs = {}
    # For each alignment depth file
    for aln_file in sample_list:
        # Find locations that are wrong but have over threshold depth
        with open(input_parameter.folder+'/'+aln_file+'.sortedByName.depth_unique.txt','r') as f:
            for l in f:
                l = l.split()
                if ref_name in l[0] and int(l[2]) >= input_parameter.threshold:
                    if (l[0], int(l[1])) in temp_contigs.keys():
                        temp_contigs[(l[0], int(l[1]))] += int(l[2])
                    else:
                        temp_contigs[(l[0], int(l[1]))] = int(l[2])
            
    # Organize dictionary into list of sorted positions
    ref_position_list = list(temp_contigs.keys())
    ref_position_list.sort() # in place sort
    # Store contig name, position, and coverage
    with open(output_filename, 'w') as outfile:
        for key in ref_position_list:
            # order is contig name \t position \t coverage
            t = outfile.write(key[0] + '\t' + str(key[1]) + '\t' + str(temp_contigs[key]) + '\n')
    print('Completed Processing of '+ref_name)

        
# When running the script from command line, the following lines are executed
if __name__ == "__main__":
    usage = "USAGE: python extract_masked_regions.py folder strain_file threshold"

    # Making default argument list structures
    p = argparse.ArgumentParser(usage=usage)
    p.add_argument(dest='folder', action='store', type=str)
    p.add_argument(dest='file', action='store', type=str)
    p.add_argument(dest='threshold', action='store', type=int)

    A = p.parse_args()

    try:
        print(A)
        coreNum = 4

        # Open file and read in all the strain names and files, should be a csv file
        with open(A.file,'r') as f:
            # get header, after the first comma
            sample_root = f.readline().rstrip().split(',')
            ref_grid = {};
            for l in f:
                l = l.rstrip().split(',')
                ref_grid[l[0]] = []
                for n in range(1,len(l)):
                    if l[n] != 'N':
                        # l[0] is the name of the genome
                        ref_grid[l[0]].append(sample_root[n])
        # Error checking
        for ref_name in list(ref_grid.keys()):
            if ref_name in ref_grid[ref_name]:
                print(ref_name + ' is wrong: ' + ref_grid[ref_name])

        # Setup the queue
        procQueue = Queue()
        for ref_name in list(ref_grid.keys()):
            procQueue.put(ref_name)

        # Divde all the files into groups corresponding to number of cores.
        coreCnt = 0
        proc = []
        while not procQueue.empty():
            if coreCnt < coreNum:
                ref_name = procQueue.get() # should be ref_name, which is a string
                p = Process(target=process_one_reference, args=(A, ref_name, ref_grid[ref_name]))
                p.start()
                proc.append(p)
                coreCnt += 1
                sleep(5) # put 10 sec in between calling new threads 
            else: # all 16 cores are used so wait until things finish
                for p in proc:
                    p.join()
                proc = []
                coreCnt = 0
        # Finish up the currently running processes
        for p in proc:
            p.join()


        # Print completion statement
        print('Script completed')

    except ValueError as e:
        print("ERROR: ValueError:",e)
        print(usage)
    except TypeError as e:
        print("ERROR: TypeError:",e)
        print(usage)
    except IOError as e:
        print("ERROR: IOError %s",e)
        print(usage)

