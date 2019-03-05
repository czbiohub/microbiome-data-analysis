import argparse, os

def bed_to_genome_length(bed_file_list, output_file):
    """
    take one file at a time and calculate the total length of all regions for each genome
    return to output
    """
    with open(bed_file_list,'r') as f:
        filelist = []
        for l in f:
            filelist.append(l.split()[0])
    
    arrays = [{} for i in range(len(filelist))]

    for m in range(len(filelist)):
        filename = filelist[m]
        with open(filename, 'r') as f:
            for l in f:
                l = l.split()
                genome_name = l[0].split('_')[0]
                region_size = int(l[2]) - int(l[1])
                if genome_name not in list(arrays[m].keys()):
                    arrays[m][genome_name] = region_size
                else:
                    arrays[m][genome_name] += region_size
    # output to file
    with open(output_file, 'w') as f:
        t = f.write(','+','.join(filelist)+'\n')
        for key in list(arrays[0].keys()):
            t = f.write(key)
            for i in range(len(filelist)):
                if key in list(arrays[i].keys()):
                    t = f.write(','+str(arrays[i][key]))
                else:
                    t = f.write(',0')
            t = f.write('\n')

        
# When running the script from command line, the following lines are executed
if __name__ == "__main__":
    usage = "USAGE: python depth_to_genome_length.py list_of_bed_file_names output_file"

    # Making default argument list structures
    p = argparse.ArgumentParser(usage=usage)
    p.add_argument(dest='bed_file', action='store', type=str)
    p.add_argument(dest='output', action='store', type=str)

    A = p.parse_args()

    try:        
        bed_to_genome_length(A.bed_file, A.output)
    except ValueError as e:
        print("ERROR: ValueError:",e)
        print(usage)
    except TypeError as e:
        print("ERROR: TypeError:",e)
        print(usage)
    except IOError as e:
        print("ERROR: IOError %s",e)
        print(usage)

