import os, sys, argparse, subprocess, gzip

usage = "USAGE: python sample_sheet_process.py [options] -i input_sample_sheet.csv "

# Making default argument list structures
p = argparse.ArgumentParser(usage=usage)
p.add_argument('-i', '--input', dest='sample_sheet', action='store', type=str)
p.add_argument('-t', '--threads', dest='coreNum', action='store', type=int, default=8)
p.add_argument('-p', '--nanopore_path', dest='folder_path', action='store', type=str, default='s3://czb-seqbot/nanopore/nanopore-samplesheets')
p.add_argument('--trial', dest='test_sample_sheet', action='store_true', default=False)

# Parse arguments
args = p.parse_args()
correct_header = "Study_ID,Study_Description,BioSample_ID,BioSample_Description,Sample_ID,Sample_Name,Sample_Owner,Index".split(',')

# If trial run, copy the sample sheet over (assumption is that it exists)
# Check that all the appropriate lines exist and the number of libraries agree
if args.test_sample_sheet:
    # copy sample sheet from S3
    t = subprocess.call(['aws', 's3', 'cp', args.folder_path+'/'+args.sample_sheet, '.'])
    with open(args.sample_sheet, 'r') as f:
        c = []
        data_flag = False
        flowcell_flag = False
        kit_flag = False
        barcoded_flag = False
        barcodekit_flag = False
        library_num = False
        for l in f:
            if data_flag:
                c.append(l.rstrip().split(','))
            else:
                if '[Data]' in l:
                    data_flag = True
                if 'Flow Cell' in l and l.rstrip().split(',')[1]:
                    flowcell_flag = True
                if 'Kit' in l and l.rstrip().split(',')[1]:
                    kit_flag = True
                if 'Barcoded' in l and l.rstrip().split(',')[1]:
                    barcoded_flag = True
                if 'Barcode Kit' in l and l.rstrip().split(',')[1]:
                    barcodekit_flag = True
                if 'Library Number' in l and l.rstrip().split(',')[1]:
                    library_num = int(l.rstrip().split(',')[1])
    header = c[0]
    data = c[1:len(c)]

    # Could check to see that there are no spaces in Study_ID, Sample_ID and Sample_Name
    if len(data) != library_num:
        raise Exception("Library number and the number of [Data] lines in sample sheet do not agree.")
    elif not (flowcell_flag and kit_flag and barcoded_flag and barcodekit_flag):
        raise Exception("Some run parameters under [Heeader] are missing.")
    else:
        for x in correct_header:
            if x not in header:
                raise Exception(x + ' seems to be missing certain columns.')

# If this is not a trial run
else:

    with open(args.sample_sheet, 'r') as f:
        c = []
        data_flag = False
        library_num = False
        for l in f:
            if data_flag:
                c.append(l.rstrip().split(','))
            else:
                if '[Data]' in l:
                    data_flag = True
                if 'Library Number' in l and l.rstrip().split(',')[1]:
                    library_num = int(l.rstrip().split(',')[1])
    header = c[0]
    data = c[1:len(c)]
    # data is a dictionary
    data = {header[j]:[data[i][j] for i in range(len(data))] for j in range(len(header))}
    runID = args.sample_sheet.split('.')[0] # yyyymmdd_....
    # If there are different Study_IDs, split out by folder
    if len(list(set(data['Study_ID']))) > 1:
        split_by_studyID = True
    else:
        split_by_studyID = False

    # Perform porechop on all the Index that exist
    good_barcode = []
    for libnum in range(library_num):
        # check barcode exists is the same
        barcode_string = data['Index'][libnum]
        barcode_number = int(barcode_string.lower().split('barcode')[1])
        barcode_string = 'barcode{:02.0f}'.format(barcode_number) # to make barcode names consistent
        good_barcode.append(barcode_string)
        # If the barcode folder is there
        output_fastq_name = data['Sample_Name'][libnum]+'__'+barcode_string+'.fastq.gz';
        if os.path.isdir('guppy_fastqs/'+barcode_string):
            subprocess.call(['porechop','-i','guppy_fastqs/'+barcode_string,'-o',\
            '01_BASECALLED/'+output_fastq_name,'--discard_middle','--verbosity','1','--threads',str(args.coreNum)])
            if split_by_studyID and data['Sample_Name'][libnum] != data['Sample_ID'][libnum]:
                subprocess.call(['mv', '01_BASECALLED/'+output_fastq_name, '01_BASECALLED/'+data['Study_ID'][libnum]+'/'+data['Sample_ID'][libnum]+'/'+output_fastq_name])
            elif split_by_studyID and data['Sample_Name'][libnum] == data['Sample_ID'][libnum]:
                subprocess.call(['mv', '01_BASECALLED/'+output_fastq_name, '01_BASECALLED/'+data['Study_ID'][libnum]+'/'+output_fastq_name])
            elif not split_by_studyID and data['Sample_Name'][libnum] != data['Sample_ID'][libnum]:
                subprocess.call(['mv', '01_BASECALLED/'+output_fastq_name, '01_BASECALLED/'+data['Sample_ID'][libnum]+'/'+output_fastq_name])

    # for the rest of the barcodes, if the number of reads is high, then output it, else discard. This will include unclassified
    for root, dir, files in os.walk('guppy_fastqs'):
        for barcode in dir:
            if barcode not in good_barcode:
                # check how many reads total
                fastq_list = os.listdir('guppy_fastqs/'+barcode)
                fastq_list = [x for x in fastq_list if 'fastq.gz' in x]
                total_reads = 0
                for filename in fastq_list:
                    with gzip.open('guppy_fastqs/'+barcode+'/'+filename, 'rb') as f:
                        for i, l in enumerate(f):
                            total_reads += 1
                # If there are more than 1000 reads, trim and output the fastq
                if total_reads / 4 >= 1000:
                    subprocess.call(['porechop','-i','guppy_fastqs/'+barcode,'-o',\
                    '01_BASECALLED/'+barcode+'__'+runID+'.fastq.gz','--discard_middle','--verbosity','1','--threads',str(args.coreNum)])
