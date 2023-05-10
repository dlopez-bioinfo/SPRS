#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process VCF_MERGE_BY_CHROM {
  label 'parallel'
  container "${params.container__bcftools}"
  publishDir "${params.output_folder}/vcf/", mode: 'copy', overwrite: true

  input:
    path(vcf_idx_list)
    val chrom

  output:
    tuple val(chrom), path("ALL_CHR${chrom}.vcf.gz") 

  script:
    def out = "ALL_CHR${chrom}.vcf.gz"

    """
      bcftools merge \
        --threads ${task.cpus} \
        --missing-to-ref \
        -r ${chrom} \
        -Oz -o ${out} \
        *.vcf.gz
    """
}
