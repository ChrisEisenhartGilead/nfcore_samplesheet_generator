#!/usr/bin/env nextflow
// Import the required Nextflow modules
// nextflow.enable.dsl=2

/*
 * pipeline input parameters
 */
params.s3_dir = ""
params.pipeline = ""
params.output_file = ""
params.sarek_status = ""
params.sarek_mapping_table = ""

log.info """\
    N F   P I P E L I N E   S A M P L E S H E E T S
    ===================================
    s3_dir       : ${params.transcriptome_file}
    pipeline     : ${params.reads}
    output_file       : ${params.outdir}
    """
    .stripIndent()


// Define the process for the first step
process generate_samplesheet {
  input:
  val s3Dir
 
  output:
  file 'output.csv'
  
  script:
  """
  # Assuming you have a Python script named 'generate_nf_core_sample_sheet.py' in the same directory
  python '${baseDir}/generate_nf_core_sample_sheet.py' rnaseq $s3Dir output.csv
  """
}

// Define the process for the second step
process validate_samplesheet {
  input:
  file csvFile
  val pipeline
  
  output:
  file 'validated_output.csv'
 
  script:
  """
  # Assuming you have a Python script named 'validate_nf_core_sample_sheet.py' in the same directory
  python '${baseDir}/validate_nf_core_sample_sheet.py' $pipeline $csvFile validated_output.csv
  """
}

// Define a process to upload the file to S3
process upload {
  input:
  file validatedCsvFile
  val outputFileName

  script:
  """
  # Assuming you have AWS CLI configured and installed
  aws s3 cp $validatedCsvFile s3://cb-multimodal-integration/work/runs/samplesheets/$outputFileName
  """
}

// Define the workflow execution order
workflow {
  // Define a channel to pass the S3 directory to the first step
  s3DirChannel = Channel.fromPath(params.s3Dir)

  // Define the first step and connect it to the input channel
  generate_samplesheet(s3DirChannel) 

  // Connect the output of the first step to the input of the second step
  validate_samplesheet(generate_samplesheet.out, params.pipeline)
  
  // Define the final step to upload the validated CSV file back to S3
  upload(validate_samplesheet.out, params.outputFileName)
}

