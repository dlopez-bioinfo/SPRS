#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process VCF_CONCAT {
  label 'parallel'
  container "${params.container__bcftools}"
  publishDir "${params.output_folder}/vcf/", mode: 'copy', overwrite: true


  input:
    tuple val(meta),
          path(vcf_with_index_list)

  output:
    tuple val(meta),
          path("${meta}_concat.vcf.gz") 

  script:
    """
      bcftools concat \
        --threads ${task.cpus} \
        -Oz -o ${meta}_concat.vcf.gz \
        *.vcf.gz 
    """
}
