#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process VCF_INDEX {
  label 'parallel'
  container "${params.container__bcftools}"
                  
  input:
    tuple val(meta), path(vcf)

  output:
    tuple val(meta), path("${vcf}"), path("${vcf}.csi") 

  script:   
    """
      bcftools index \
        --threads ${task.cpus} \
        ${vcf}
    """
}
