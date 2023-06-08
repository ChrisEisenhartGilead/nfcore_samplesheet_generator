#!/usr/bin/env python2.7
# Chris Eisenhart
"""
Generate sample sheets for Gilead nf-core pipelines from an AWS s3 bucket.  It is assumed that
.fastq.gz files with a semi standard naming structure are in the provided S3 bucket.  The script
will attempt to identify the samples present, and create a samplesheet that can be passed into
the corresponding Gilead nf-core pipeline.
"""
import sys
import argparse
import logging
import subprocess

log = None  # The global logger is defined, it will be initiated after the user commands are read


def validate_options(options):
    """
    Use this to perform additional checks on parameters and to enable conditional mandatory arguments (ex Sarek)
    """
    if options.s3_bucket.startswith("s3"):
        print("The S3 bucket should not include 's3:/', it should be format /my_bucket/part1/part2/")
        exit(1)
    if options.pipeline == "sarek":
        if options.sarek_status is None:
            print("The --sarek_status option is required for the 'sarek' pipeline.")
            exit(1)
        if options.sarek_status == "paired":
            if options.sarek_mapping_file is None:
                print(
                    "The --sarek_mapping_file option is required for the 'sarek' 'paired' pipeline."
                )
                exit(1)


def parseArgs(args):
    """
    Parse the command line arguments into useful python objects.  '--' variables are optional
    set_defaults only applies when the argument is not provided (it won't override)
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "pipeline",
        help=" The pipeline which the samplesheet belongs too",
        choices=["sarek", "rnaseq"],
        action="store",
    )
    parser.add_argument("s3_bucket", help=" The input s3 bucket", action="store")
    parser.add_argument("outFile", help=" The output file", action="store")
    parser.add_argument(
        "--sarek_status",
        help=" The sarek run status (germline/tumor/paired)",
        choices=["tumor", "germline", "paired"],
        action="store",
    )
    parser.add_argument(
        "--sarek_mapping_file",
        help=" A two column tsv mapping sample names to 0 or 1 where 0 is germline and 1 is tumor",
        action="store",
    )
    parser.add_argument(
        "--verbose",
        help=" The verbosity level for stdout messages (default INFO)",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        action="store",
    )
    parser.set_defaults(verbose="INFO")
    options = parser.parse_args()
    validate_options(options)

    # Basic logging config and update the logging verbosity level
    logging.basicConfig(
        format="[%(filename)s %(funcName)s %(lineno)d %(levelname)s] %(message)s",
        level=options.verbose,
    )

    # Define a file for the log object to write too
    fh = logging.FileHandler(__file__ + ".log")

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
        files = [
            line.split()[-1]
            for line in result.stdout.splitlines()
            if line.endswith(".fastq.gz")
        ]
        return files
    else:
        # Handle the case when the command fails
        print("Failed to list S3 files.")
        print(result.stderr)
        return []


def parse_sample_data(s3_path):
    """
    Takes an S3 like path, ie;
    p573/Fulgent/Batch1/SH8611/FT-SA146502D-FT-SPN03430_H52NCDSX5/SH8611_SA146502D_S63_L004_R1_001.fastq.gz
    and returns a sample name.
    SH8611_SA146502D_S63
    """
    # Find the index of the lane element of the fastq file if it exists (ex; L001)
    def find_lane_index(elements):
        for i, element in enumerate(elements):
            if element.startswith("L") and len(element) == 4 and element[1:].isdigit():
                return i
        return None

    # Find the index of the orientation element of the fastq file (R1 or R2)
    def find_orient_index(elements):
        for i, element in enumerate(elements):
            if element.startswith("R") and len(element) == 2 and element[1:].isdigit():
                return i
        return None

    # Split the S3 path by '/'
    components = s3_path.split("/")
    # Split the final element by _ and remove the extension
    sub_components = components[-1].replace(".fastq.gz", "").split("_")
    # Identify the lane index
    lane_index = find_lane_index(sub_components)
    if lane_index:
        sample_name = "_".join(sub_components[:lane_index])
        lane = sub_components[lane_index]
    else:
        # Identify the orientation index
        orient_index = find_orient_index(sub_components)
        sample_name = "_".join(sub_components[:orient_index])
        lane = ""

    # Look for Fulgent specific patterns ex; FT-SA204671R-FT-SPN04508_H5VHMDSX7 or FT-SPN04508_H5VHMDSX7-FT-SA204671R
    for comp in components:
        if "FT-SA" in comp: # This identifies a Fulgent Sample accession
            acc_loc = comp.index("FT-SA") # Identify where in this component we will see the accession
            sub_components = comp.split("-")
            cur_loc = 0
            found_acc = False
            ft_sample_accession = ""
            for sub_comp in sub_components:
                if found_acc:
                    ft_sample_accession = sub_comp
                    assert ft_sample_accession.startswith("SA") # Verify at the correct component
                    break
                if cur_loc == acc_loc: # Reached the 'FT' prior to the Sample Accession in the string
                    found_acc = True
                cur_loc += len(sub_comp) + 1
            sample_name = sub_comp
            break

    return sample_name, lane


def pair_fastq_file_list(fastq_file_list):
    """
    Takes a list of fastq files and identifies if any are unpaired, if so the program exits. For
    each paired fastq, the sample name and lane are inferred and stored along the other data. If
    no unpaired fastq are found, then a list of the paired fastq files is returned.
    """
    # Separate R1 and R2 files and flag any unpaired files
    paired_files = []
    unpaired_files = []
    for fastq_file in fastq_file_list:
        if "_R1" in fastq_file:
            sample_name, lane = parse_sample_data(fastq_file)
            r1_file = fastq_file
            r2_file = fastq_file.replace("_R1", "_R2")
            if r2_file in fastq_file_list:
                paired_files.append((sample_name, lane, r1_file, r2_file))
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


def generate_nf_core_rnaseq_sample_sheet(paired_files, s3_bucket):
    """
    sample,fastq_1,fastq_2,strandedness
    SAMPLE_PE,SAMPLE_PE_RUN1_1.fastq.gz,SAMPLE_PE_RUN1_2.fastq.gz,forward
    SAMPLE_PE,SAMPLE_PE_RUN2_1.fastq.gz,SAMPLE_PE_RUN2_2.fastq.gz,forward
    SAMPLE_SE,SAMPLE_SE_RUN1_1.fastq.gz,,forward
    """
    results = ["sample,fastq_1,fastq_2,strandedness"]
    s3_base = s3_bucket.split("/")[1]
    for data in paired_files:
        sample_name, lane, r1_fastq, r2_fastq = data
        sample_sheet_entry = f"{sample_name},s3://{s3_base}/{r1_fastq},s3://{s3_base}/{r2_fastq},auto"
        results.append(sample_sheet_entry)
    return results


def get_sample_name_to_status(mapping_file, paired_files):
    """
    Parse user provided mapping file to generate a sample mapping dict which
    maps sample names to their status ex; {"SA204434": "germline"}
    """
    # Dict to map the sample name to the sample status (tumor/germline)
    sample_name_to_status = {}
    with open(mapping_file, "r") as file:
        count = 0
        for line in file:
            count += 1
            key, value = line.strip().split("\t")
            found_sample_name = None
            status = None
            if value == "tumor":
                status = 1
            elif value == "germline":
                status = 0
            else:
                print(
                    f'Mapping file {mapping_file} has value other than "tumor" or "germline" on line {count}'
                )
                exit(1)
            for data in paired_files:
                if (
                    key in data[0] or key in data[2]
                ):  # Check the sample name and the fastq
                    found_sample_name = data[0]
                if found_sample_name:
                    break
            if not found_sample_name:
                print(
                    f"Mapping file {mapping_file} has sample name {key} on line {count} which was \
                         not found in any fastq files in the S3 bucket indicated"
                )
                exit(1)
            sample_name_to_status[found_sample_name] = status
    return sample_name_to_status


def generate_nf_core_sarek_sample_sheet(paired_files, s3_bucket, in_status, mapping_file):
    """
    patient,sex,status,sample,lane,fastq_1,fastq_2
    patient1,XX,0,normal_sample,lane_1,test_L001_1.fastq.gz,test_L001_2.fastq.gz
    patient1,XX,0,normal_sample,lane_2,test_L002_1.fastq.gz,test_L002_2.fastq.gz
    patient2,NA,0,normal_sample,lane_1,test_L003_1.fastq.gz,test_L003_2.fastq.gz
    """
    results = ["patient,sex,status,sample,lane,fastq_1,fastq_2"]
    s3_base = s3_bucket.split("/")[1]
    if in_status == "paired":
        sample_name_to_status = get_sample_name_to_status(mapping_file, paired_files)
    status = 1 if in_status == "tumor" else 0
    for data in paired_files:
        sample_name, lane, r1_fastq, r2_fastq = data
        if in_status == "paired":
            status = sample_name_to_status[sample_name]
        sample_sheet_entry = (
            f"{sample_name},NA,{status},{sample_name},{lane},s3://{s3_base}/{r1_fastq},s3://{s3_base}/{r2_fastq}"
        )
        results.append(sample_sheet_entry)
    return results


def main(args):
    """
    Parse command line options and call the appropriate function.
    """
    options = parseArgs(args)
    if options.pipeline == "rnaseq":
        fastq_files = get_s3_fastq_files(options.s3_bucket)
        paired_files = pair_fastq_file_list(fastq_files)
        output_data = generate_nf_core_rnaseq_sample_sheet(paired_files, options.s3_bucket)
    elif options.pipeline == "sarek":
        fastq_files = get_s3_fastq_files(options.s3_bucket)
        paired_files = pair_fastq_file_list(fastq_files)
        output_data = generate_nf_core_sarek_sample_sheet(
            paired_files, options.s3_bucket, options.sarek_status, options.sarek_mapping_file
        )
    with open(options.outFile, "w") as out_file:
        out_file.write("\n".join(output_data))


if __name__ == "__main__":
    sys.exit(main(sys.argv))
