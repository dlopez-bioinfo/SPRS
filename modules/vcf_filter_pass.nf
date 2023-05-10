#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process VCF_FILTER_PASS {
  label 'parallel'
  container "${params.container__bcftools}"
                  
  input:
    tuple val(sampleId), path(vcf)
    val exome
    path capture_bed

  output:
    tuple val(sampleId), path("${sampleId}_pass.vcf.gz")

  script:
    def arg_target = exome ? "-T " + capture_bed : "-t " + (1..22).join(",") + ",X,Y"
    def out = "${sampleId}_pass.vcf.gz"

    """
      bcftools view \
        --threads ${task.cpus} \
        ${arg_target} \
        ${vcf} \
      | bcftools view \
        --threads ${task.cpus} \
        --apply-filters PASS \
      | bcftools view \
        --threads ${task.cpus} \
        --trim-alt-alleles \
      | bcftools view \
        --threads ${task.cpus} \
        --min-ac 1 \
      | bcftools view \
        --threads ${task.cpus} \
        --exclude 'AN=1' \
        -Oz -o ${out} 
    """
}
