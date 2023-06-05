#!/usr/bin/env python2.7
# Chris Eisenhart 

"""

"""
import os
import sys
import argparse
import logging
import subprocess
import re

log = None # The global logger is defined, it will be initiated after the user commands are read

def parseArgs(args): 
    """
    Parse the command line arguments into useful python objects.  '--' variables are optional
    set_defaults only applies when the argument is not provided (it won't override)
    """
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument("pipeline",
                        help = " The pipeline which the samplesheet belongs too",
                        choices=["sarek","rnaseq"],
                        action = "store")
    parser.add_argument("s3_bucket",
                        help = " The input s3 bucket",
                        action = "store")
    parser.add_argument("outFile",
                        help = " The output file",
                        action = "store")
    parser.add_argument("--verbose",
                        help = " The verbosity level for stdout messages (default INFO)",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        action = "store")

    parser.set_defaults(verbose = "INFO")

    options = parser.parse_args()

    # Basic logging config and update the logging verbosity level
    logging.basicConfig(format='[%(filename)s %(funcName)s %(lineno)d %(levelname)s] %(message)s', level=options.verbose)

    # Define a file for the log object to write too
    fh = logging.FileHandler(__file__ + '.log')

    # Create the logger and add the file handler to it
    global log
    log = logging.getLogger(__name__)
    log.addHandler(fh)
    return options



def get_s3_fastq_files(s3_directory):
    """
    Takes an S3 directory and returns a list of the fastq files that are present.
    """
    # Execute the aws s3 ls --recursive command
    command = f"aws s3 ls s3:/{s3_directory} --recursive"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    # Check if the command was successful
    if result.returncode == 0:
        # Extract the file names with the .fastq.gz extension
        files = [line.split()[-1] for line in result.stdout.splitlines() if line.endswith(".fastq.gz")]
        return files
    else:
        # Handle the case when the command fails
        print("Failed to list S3 files.")
        print(result.stderr)
        return []


def parse_sample_name(s3_path):
    """
    Takes an S3 like path, ie; 
    p573/Fulgent/Batch1/SH8611/FT-SA146502D-FT-SPN03430_H52NCDSX5/SH8611_SA146502D_S63_L004_R1_001.fastq.gz
    and returns a sample name.
    SH8611_SA146502D_S63
    """
    def find_lane_index(elements):
        for i, element in enumerate(elements):
            if element.startswith('L') and len(element) == 4 and element[1:].isdigit():
                return i
        return None
    def find_orient_index(elements):
        for i, element in enumerate(elements):
            if element.startswith('R') and len(element) == 2 and element[1:].isdigit():
                return i
        return None
    # Split the S3 path by '/'
    components = s3_path.split('/')
    # Split the final element by _ and remove the extension
    sub_components = components[-1].replace(".fastq.gz","").split("_")
    # Identify the lane index
    lane_index = find_lane_index(sub_components)
    if lane_index:
        sample_name = "_".join(sub_components[:lane_index])
    else:
        # Identify the orientation index
        orient_index = find_orient_index(sub_components)
        sample_name = "_".join(sub_components[:orient_index])
    return sample_name


def pair_fastq_file_list(fastq_file_list):
    """
    Takes a list of fastq files and identifies if any are unpaired, if so the program exits.
    If no unpaired fastq are found, then a list of the paired fastq files is returned.
    """
    # Separate R1 and R2 files and flag any unpaired files
    paired_files = []
    unpaired_files = []
    for fastq_file in fastq_file_list:
        if "_R1" in fastq_file:
            sample_name = parse_sample_name(fastq_file)
            r1_file = fastq_file
            r2_file = fastq_file.replace("_R1", "_R2")
            if r2_file in fastq_file_list:
                paired_files.append((sample_name,r1_file, r2_file))
            else:
                unpaired_files.append(r1_file)
        elif "_R2" in fastq_file:
            # Skip R2 files as they will be paired with R1 files
            continue
        else:
            # Skip files that do not follow the R1 or R2 naming convention
            continue
    if unpaired_files:
        print("Unpaired fastq files found! Please handle manually.")
        exit(1)
    else:
        return paired_files


def generate_nf_core_rnaseq_sample_sheet(paired_files):
    """
    sample,fastq_1,fastq_2,strandedness
    SAMPLE_PE,SAMPLE_PE_RUN1_1.fastq.gz,SAMPLE_PE_RUN1_2.fastq.gz,forward
    SAMPLE_PE,SAMPLE_PE_RUN2_1.fastq.gz,SAMPLE_PE_RUN2_2.fastq.gz,forward
    SAMPLE_SE,SAMPLE_SE_RUN1_1.fastq.gz,,forward
    """
    results = ["sample,fastq_1,fastq_2,strandedness"]
    for data in paired_files:
        sample_name, r1_fastq, r2_fastq = data
        sample_sheet_entry = f'{sample_name},{r1_fastq},{r2_fastq},auto'
        results.append(sample_sheet_entry)
    return results


def main(args):
    """
    """ 
    options = parseArgs(args)
    if options.pipeline == "rnaseq":
        fastq_files = get_s3_fastq_files(options.s3_bucket)
        paired_files = pair_fastq_file_list(fastq_files)
        output_data = generate_nf_core_rnaseq_sample_sheet(paired_files)
    elif options.pipeline == "sarek":
        return 
    with open(options.outFile, 'w') as out_file:
        out_file.write("\n".join(output_data))

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
