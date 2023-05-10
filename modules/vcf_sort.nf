#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process VCF_SORT {
  container "${params.container__bcftools}"
                  
  input:
    tuple val(sampleId), path(vcf)

  output:
    tuple val(sampleId), path("${sampleId}_sorted.vcf.gz") 

  script:   
    def out = "${sampleId}_sorted.vcf.gz"

    """
      bcftools sort \
        -Oz -o ${out} \
        ${vcf}
    """
}
