import argparse, os

def depth_to_bed(depth_file, bed_file, min_region_size, min_gap):
    """
    Create a basic 3 column bed file from a samtools depth output.
    The min_region_size and min_gap are used to adjust the granularity
    of the regions created if the depth file do not have continuous
    positions in all regions.
    min_region_size: If a regions is created, it needs to to this size
    min_gap: minimum distance in bp between two adjacent regions.
    """

    with open(depth_file, 'r') as input_file, open(bed_file, 'w') as output_file:
        # Initialize the variable place holders
        left = None
        right = None
        prev_contig_name = None
        for l in input_file:
            l = l.split()
            new_contig_name = l[0]
            pos = int(l[1]) # position on the contig
            if left is None: # new region has not been assigned at all
                left = pos - 1 # because in a bed file, left is 0 indexed
            elif right is None: # left is assigned but right is not
                # the current position is too far so a new region, or if pos on a new contig
                # This means that the region is only 1 bp. 
                # For this application I'm making such regions into min regions
                if pos - left >= min_gap or new_contig_name != prev_contig_name: 
                    right = left + int(min_region_size / 2) 
                    left = max(0, right - min_region_size)
                    t = output_file.write(prev_contig_name + '\t' + str(left) + '\t' + str(right) + '\n')
                    left = pos
                    right = None
                else: # current position is still within the min gap so still in the same region
                    right = pos
            else: # both left and right are assigned and it's on the same contig
                # Location of the current position is too far to the right. It belongs to a new region
                if pos - right >= min_gap or new_contig_name != prev_contig_name: 
                    # the region is already larger than min_region_size
                    if right - left >= min_region_size:
                        t = output_file.write(prev_contig_name + '\t' + str(left) + '\t' + str(right) + '\n')
                        left = pos
                        right = None
                    else: # not a large enough region, create an extended region centered on right - left
                        right = int((left + right) / 2 + (min_region_size / 2))
                        left = max(0, right - min_region_size)
                        t = output_file.write(prev_contig_name + '\t' + str(left) + '\t' + str(right) + '\n')
                        left = pos
                        right = None
                else: # current position is still in the same region
                    right = pos
            # Update contig name
            prev_contig_name = new_contig_name

        
# When running the script from command line, the following lines are executed
if __name__ == "__main__":
    usage = "USAGE: python depth_to_bed.py depth_file_name bed_file_name min_region_size min_gap"

    # Making default argument list structures
    p = argparse.ArgumentParser(usage=usage)
    p.add_argument(dest='depth_file', action='store', type=str)
    p.add_argument(dest='bed_file', action='store', type=str)
    p.add_argument(dest='min_region_size', action='store', type=int, default=100)
    p.add_argument(dest='min_gap', action='store', type=int, default=1000)

    A = p.parse_args()

    try:        
        depth_to_bed(A.depth_file, A.bed_file, A.min_region_size, A.min_gap)
    except ValueError as e:
        print("ERROR: ValueError:",e)
        print(usage)
    except TypeError as e:
        print("ERROR: TypeError:",e)
        print(usage)
    except IOError as e:
        print("ERROR: IOError %s",e)
        print(usage)

