#!/usr/bin/env 
# Chris Eisenhart 

"""
Take in a nfcore sample sheet and determine whether the formatting is valid
"""
import sys
import argparse
import logging
import shutil
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
    parser.add_argument("sampleSheet",
                        help = " The input nfcore RNAseq sample sheet",
                        action = "store")
    parser.add_argument("outFile",
                        help = " The validated output sample sheet",
                        action = "store")
    parser.add_argument("--verbose",
                        help = " The verbosity level for stdout messages (default INFO)",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        action = "store")

    parser.set_defaults(verbose = "INFO")
    parser.set_defaults(isTsv = False)

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


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = f"ERROR: Please check samplesheet -> {error}"
    if context != "" and context_str != "":
        error_str = f"ERROR: Please check samplesheet -> {error}\n{context.strip()}: '{context_str.strip()}'"
    print(error_str)


def validate_rnaseq_samplesheet(file_in):
    """
    This function is based on code from https://github.com/nf-core/rnaseq/blob/master/bin/check_samplesheet.py

    This function checks that the samplesheet follows the following structure:

    sample,fastq_1,fastq_2,strandedness
    SAMPLE_PE,SAMPLE_PE_RUN1_1.fastq.gz,SAMPLE_PE_RUN1_2.fastq.gz,forward
    SAMPLE_PE,SAMPLE_PE_RUN2_1.fastq.gz,SAMPLE_PE_RUN2_2.fastq.gz,forward
    SAMPLE_SE,SAMPLE_SE_RUN1_1.fastq.gz,,forward

    For an example see:
    https://github.com/nf-core/test-datasets/blob/rnaseq/samplesheet/v3.1/samplesheet_test.csv
    """
    # Assume the file is valid to start, any criteria failing will invalidate it
    is_valid = True
    sample_mapping_dict = {}
    with open(file_in, "r", encoding="utf-8-sig") as fin:
        ## Check header
        MIN_COLS = 3
        HEADER = ["sample", "fastq_1", "fastq_2", "strandedness"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print(f"ERROR: Please check samplesheet header -> {','.join(header)} != {','.join(HEADER)}")
            is_valid = False

        ## Check sample entries
        for line in fin:
            if line.strip():
                lspl = [x.strip().strip('"') for x in line.strip().split(",")]

                ## Check valid number of columns per row
                if len(lspl) < len(HEADER):
                    print_error(
                        f"Invalid number of columns (minimum = {len(HEADER)})!",
                        "Line",
                        line,
                    )
                    is_valid = False

                num_cols = len([x for x in lspl[: len(HEADER)] if x])
                if num_cols < MIN_COLS:
                    print_error(
                        f"Invalid number of populated columns (minimum = {MIN_COLS})!",
                        "Line",
                        line,
                    )
                    is_valid = False

                ## Check sample name entries
                sample, fastq_1, fastq_2, strandedness = lspl[: len(HEADER)]
                if sample.find(" ") != -1:
                    print(f"WARNING: Spaces will be replaced by underscores for sample: {sample}")
                if not sample:
                    print_error("Sample entry has not been specified!", "Line", line)
                    is_valid = False

                ## Check FastQ file extension
                for fastq in [fastq_1, fastq_2]:
                    if fastq:
                        if fastq.find(" ") != -1:
                            print_error("FastQ file contains spaces!", "Line", line)
                        if not fastq.endswith(".fastq.gz") and not fastq.endswith(".fq.gz"):
                            print_error(
                                "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
                                "Line",
                                line,
                            )
                            is_valid = False

                ## Check strandedness
                strandednesses = ["unstranded", "forward", "reverse", "auto"]
                if strandedness:
                    if strandedness not in strandednesses:
                        print_error(
                            f"Strandedness must be one of '{', '.join(strandednesses)}'!",
                            "Line",
                            line,
                        )
                        is_valid = False
                else:
                    print_error(
                        f"Strandedness has not been specified! Must be one of {', '.join(strandednesses)}.",
                        "Line",
                        line,
                    )
                    is_valid = False

                ## Auto-detect paired-end/single-end
                sample_info = []  ## [single_end, fastq_1, fastq_2, strandedness]
                if sample and fastq_1 and fastq_2:  ## Paired-end short reads
                    sample_info = ["0", fastq_1, fastq_2, strandedness]
                elif sample and fastq_1 and not fastq_2:  ## Single-end short reads
                    sample_info = ["1", fastq_1, fastq_2, strandedness]
                else:
                    print_error("Invalid combination of columns provided!", "Line", line)
                    is_valid = False

                ## Create sample mapping dictionary = {sample: [[ single_end, fastq_1, fastq_2, strandedness ]]}
                sample_info = sample_info + lspl[len(HEADER) :]
                if sample not in sample_mapping_dict:
                    sample_mapping_dict[sample] = [sample_info]
                else:
                    if sample_info in sample_mapping_dict[sample]:
                        print_error("Samplesheet contains duplicate rows!", "Line", line)
                        is_valid = False
                    else:
                        sample_mapping_dict[sample].append(sample_info)

    ## Write validated samplesheet with appropriate columns
    if len(sample_mapping_dict) > 0:
        for sample in sorted(sample_mapping_dict.keys()):
            ## Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
            if not all(x[0] == sample_mapping_dict[sample][0][0] for x in sample_mapping_dict[sample]):
                print_error(
                    f"Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end!",
                    "Sample",
                    sample,
                )
                is_valid = False

            ## Check that multiple runs of the same sample are of the same strandedness
            if not all(x[3] == sample_mapping_dict[sample][0][3] for x in sample_mapping_dict[sample]):
                print_error(
                    f"Multiple runs of a sample must have the same strandedness!",
                    "Sample",
                    sample,
                )
                is_valid = False
    else:
        print_error(f"No entries to process!", "Samplesheet: {file_in}")
        is_valid = False
    return is_valid 


def validate_sarek_samplesheet(file_in):
    """
    This function checks that the samplesheet follows the following structure:

    patient,sex,status,sample,lane,fastq_1,fastq_2
    patient1,XX,0,normal_sample,lane_1,test_L001_1.fastq.gz,test_L001_2.fastq.gz
    patient1,XX,0,normal_sample,lane_2,test_L002_1.fastq.gz,test_L002_2.fastq.gz
    patient1,XX,0,normal_sample,lane_3,test_L003_1.fastq.gz,test_L003_2.fastq.gz
    patient1,XX,1,tumor_sample,lane_1,test2_L001_1.fastq.gz,test2_L001_2.fastq.gz
    patient1,XX,1,tumor_sample,lane_2,test2_L002_1.fastq.gz,test2_L002_2.fastq.gz
    patient1,XX,1,relapse_sample,lane_1,test3_L001_1.fastq.gz,test3_L001_2.fastq.gz
    """
    is_valid = True
    patient_mapping_dict = {}
    with open(file_in, "r", encoding="utf-8-sig") as fin:
        ## Check header
        MIN_COLS = 7
        HEADER = ["patient","sex","status","sample","lane","fastq_1","fastq_2"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print(f"ERROR: Please check samplesheet header -> {','.join(header)} != {','.join(HEADER)}")
            is_valid = False

        ## Check sample entries
        for line in fin:
            if line.strip():
                lspl = [x.strip().strip('"') for x in line.strip().split(",")]

                ## Check valid number of columns per row
                if len(lspl) < len(HEADER):
                    print_error(
                        f"Invalid number of columns (minimum = {len(HEADER)})!",
                        "Line",
                        line,
                    )
                    is_valid = False

                num_cols = len([x for x in lspl[: len(HEADER)] if x])
                if num_cols < MIN_COLS:
                    print_error(
                        f"Invalid number of populated columns (minimum = {MIN_COLS})!",
                        "Line",
                        line,
                    )
                    is_valid = False

                ## Check sample name entries
                patient, sex, status, sample, lane, fastq1, fastq2 = lspl[:len(HEADER)]
                if patient.find(" ") != -1:
                    print(f"WARNING: Spaces will be replaced by underscores for patient: {patient}")
                if not patient:
                    print_error("Patient entry has not been specified!", "Line", line)
                    is_valid = False

                ## Check FastQ file extension
                for fastq in [fastq1, fastq2]:
                    if fastq:
                        if fastq.find(" ") != -1:
                            print_error("FastQ file contains spaces!", "Line", line)
                        if not fastq.endswith(".fastq.gz") and not fastq.endswith(".fq.gz"):
                            print_error(
                                "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
                                "Line",
                                line,
                            )
                            is_valid = False

                ## Check sex
                allowable_sex_params = ["XX", "XY", "NA"]
                if sex:
                    if sex not in allowable_sex_params:
                        print_error(
                            f"Sex must be one of '{', '.join(allowable_sex_params)}'!",
                            "Line",
                            line,
                        )
                        is_valid = False
                else:
                    print_error(
                        f"Sex has not been specified! Must be one of {', '.join(allowable_sex_params)}.",
                        "Line",
                        line,
                    )
                    is_valid = False

                ## Check status (tumor/normal) 
                allowable_status_params = ["0", "1"]
                if status:
                    if status not in allowable_status_params:
                        print_error(
                            f"Sex must be one of '{', '.join(allowable_status_params)}'!",
                            "Line",
                            line,
                        )
                        is_valid = False
                else:
                    print_error(
                        f"Status has not been specified! Must be one of {', '.join(allowable_status_params)}.",
                        "Line",
                        line,
                    )
                    is_valid = False

                ## Create patient mapping dictionary = patient: [[patient, sex, status, sample, lane, fastq1, fastq2]]
                patient_info = lspl[:len(HEADER)]
                if patient not in patient_mapping_dict:
                    patient_mapping_dict[patient] = [patient_info]
                else:
                    is_duplicate = False
                    for prev_patient_info in patient_mapping_dict[patient]:
                        if patient_info == prev_patient_info:
                            print_error("Samplesheet contains duplicate rows!", "Line", line)
                            is_valid = False
                            is_duplicate = True
                    if not is_duplicate:
                        patient_mapping_dict[patient].append(patient_info)

    ## Write validated samplesheet with appropriate columns
    if len(patient_mapping_dict) > 0:
        for patient in sorted(patient_mapping_dict.keys()):
            ## Check that multiple runs of the same patient are of the same sexs
            if not all(x[1] == patient_mapping_dict[patient][0][1] for x in patient_mapping_dict[patient]):
                print_error(
                    f"Multiple runs of a patient must have the same sex!",
                    "Sample",
                    sample,
                )
                is_valid = False
    else:
        print_error(f"No entries to process!", "Samplesheet: {file_in}")
        is_valid = False
    return is_valid


def main(args):
    """
    """ 
    options = parseArgs(args)
    if options.pipeline == "rnaseq":
        is_valid = validate_rnaseq_samplesheet(options.sampleSheet)
    elif options.pipeline == "sarek":
        is_valid = validate_sarek_samplesheet(options.sampleSheet)

    if is_valid:
        print("The samplesheet provided passed validation")
        shutil.copyfile(options.sampleSheet, options.outFile)
    else:
        print("The samplesheet provided failed validation")

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
