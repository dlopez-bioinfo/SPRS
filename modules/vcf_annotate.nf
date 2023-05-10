#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process VCF_ANNOTATE {
  label 'parallel'
  container "${params.container__bcftools}"
                  
  input:
    tuple path(vcf), path(index)
    tuple path(annotation_file), path(annotation_file_idx)
    val column_id

  output:
    path "*_annotated.vcf.gz" 

  script:
    def arg_output_vcf = "${vcf.name.replaceAll(/.vcf.gz/, '')}_annotated.vcf.gz"

    """
      bcftools annotate \
        --threads ${task.cpus} \
        --columns ${column_id} \
        --annotations ${annotation_file} \
        -Oz -o ${arg_output_vcf} \
        ${vcf} 
    """
}
