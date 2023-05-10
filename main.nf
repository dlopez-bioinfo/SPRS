#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2


/*
================================================================
  MODULES
================================================================
*/
include {VCF_INDEX} from "./modules/vcf_index"


/*
================================================================
  SUBWORKFLOWS
================================================================
*/
include { SUB_PREPARE_VCF } from "./subworkflows/prepare_vcf"
include { SUB_IMPUTE_VCF } from "./subworkflows/impute_vcf"
include { SUB_CREATE_SAMPLE_SHEET } from "./subworkflows/create_sample_sheet"


/*
================================================================
  FUNCTIONS
================================================================
*/
// Function which prints help message text
def helpMessage() {
    log.info"""
Usage:

nextflow run babelomics/csvs_prs <ARGUMENTS>

  Input Data:
    --input_folder        Folder containing input files (.vcf)
    --exome               False for WGS, True for WES or clinical exomes (default=false)
    --capture_bed         Capture bed file used for WES or clinical exomes. Ignored for WGS

  Reference:
    --shapeit4_map_files  Shapeit4 genome mapping files (.gmap.gz)
    --minimac4_ref_files  Minimac4 reference files (.bcf)
    --reference_genome    Reference genome (.fa)
    --chrom_mapping_file  Chr prefix to no chr prefix mapping file

  Output:
    --output_folder       Folder for output results. Default to ${launchDir}/results/

    """.stripIndent()
}


/*
================================================================
  MAIN WORKFLOW
================================================================
*/
workflow {
  // Show help message if the user specifies the --help flag at runtime
  // or if any required params are not provided
  if ( 
      params.help 
      || params.input_folder == false 
      || params.capture_bed == false
      || params.shapeit4_map_files == false 
      || params.minimac4_ref_files == false 
      || params.reference_genome == false
      || params.chrom_mapping_file == false){
    helpMessage()
    exit 1
  }

    
  // Make a channel with all of the files from the --input_folder
  // [ SAMPLE_ID, VCF_FILE ]
  Channel
    .fromFilePairs("${params.input_folder}/*.vcf.gz", size: 1, checkIfExists: true, flat: true)
    .dump(tag: '__INPUT_VCF__')      
    .set { vcf_ch }


  // Fix/clean VCFs
  SUB_PREPARE_VCF(vcf_ch, params.exome, params.capture_bed, params.chrom_mapping_file)

  // phase and impute missing values
  SUB_IMPUTE_VCF(SUB_PREPARE_VCF.out.prepared_vcf, params.shapeit4_map_files, params.minimac4_ref_files)

  // Create a sample sheet to feed pgsc_calc pipeline
  SUB_CREATE_SAMPLE_SHEET(SUB_IMPUTE_VCF.out)
}
