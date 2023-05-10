#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process VCF_REMOVE_CHR_PREFIX {
  label 'parallel'
  container "${params.container__bcftools}"
  errorStrategy 'retry'
  maxRetries 2
                  
  input:
    tuple val(sampleId), path(vcf)
    path rename_mapping_file

  output:
    tuple val(sampleId), path("${sampleId}_noChr.vcf.gz")

  script:   
    def out = "${sampleId}_noChr.vcf.gz"

    """
      if [[ "${task.attempt}" != "1" ]]
      then
        bcftools index --threads ${task.cpus} ${vcf}
      fi

      bcftools annotate \
        --threads ${task.cpus} \
        --rename-chrs ${rename_mapping_file} \
        -Oz -o ${out} \
        ${vcf}
    """
}
