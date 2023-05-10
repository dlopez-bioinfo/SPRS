#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process VCF_REMOVE_AF {
  label 'parallel'
  container "${params.container__bcftools}"
                  
  input:
    tuple val(sampleId), path(vcf)

  output:
    tuple val(sampleId), path("${sampleId}_noAF.vcf.gz") 

  script:
    def out = "${sampleId}_noAF.vcf.gz"

    """
      bcftools annotate \
        --threads ${task.cpus} \
        -x 'INFO,FORMAT' \
        -Oz -o ${out} \
        ${vcf}
    """
}
