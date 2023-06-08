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
process generate_rnaseq_samplesheet {
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

process generate_sarek_samplesheet {
  input:
  val s3Dir
  val sarek_status
 
  output:
  file 'output.csv'
  
  script:
  """
  # Assuming you have a Python script named 'generate_nf_core_sample_sheet.py' in the same directory
  python '${baseDir}/generate_nf_core_sample_sheet.py' sarek $s3Dir output.csv --sarek_status $sarek_status
  """
}

process generate_sarek_samplesheet_paired {
  input:
  val s3Dir
  val sarek_status
  val sarek_mapping_table
 
  output:
  file 'output.csv'
  
  script:
  """
  # Assuming you have a Python script named 'generate_nf_core_sample_sheet.py' in the same directory
  python '${baseDir}/generate_nf_core_sample_sheet.py' sarek $s3Dir output.csv --sarek_status $sarek_status --sarek_mapping_file $sarek_mapping_table
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
  if (params.pipeline == "rnaseq") {
    out = generate_rnaseq_samplesheet(s3DirChannel) 
  } else if (params.pipeline == "sarek") {
    if (params.sarek_mapping_file) {
      out = generate_sarek_samplesheet_paired(s3DirChannel, params.sarek_status, params.sarek_mapping_file) 
    } else {
      out = generate_sarek_samplesheet(s3DirChannel, params.sarek_status) 
    }
  }

  // Connect the output of the first step to the input of the second step
  validate_samplesheet(out, params.pipeline)
  
  // Define the final step to upload the validated CSV file back to S3
  upload(validate_samplesheet.out, params.outputFileName)
}

