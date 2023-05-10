#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process VCF_NORM {
  label 'parallel'
  container "${params.container__bcftools}"
                  
  input:
    tuple val(meta), path(vcf), path(ref)

  output:
    tuple val(meta), path("${meta}_norm.vcf.gz")

  script:
    def out = "${meta}_norm.vcf.gz"

    """
      bcftools norm \
        --threads ${task.cpus} \
        -m +any \
        ${vcf} \
        | bcftools norm \
            --threads ${task.cpus} \
            --rm-dup none \
            --fasta-ref ${ref} \
            -Oz -o ${out} \
    """
}
