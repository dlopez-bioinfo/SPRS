#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process VCF_IMPUTE {
  label 'parallel'
  container "${params.container__minimac4}"
  publishDir "${params.output_folder}/imputation/", mode: 'copy', overwrite: true
  errorStrategy "${ workflow.profile =~ 'test' ? 'ignore' : 'terminate' }"
                  
  input:
    tuple val(chrom),
          path(reference),
          path(vcf)

  output:
    tuple val(chrom), path("*_minimac4.dose.vcf.gz")

  script:
    def out_prefix = "${vcf.name.replaceAll(/.vcf.gz/, '')}_minimac4"    

    """
      minimac4 \
        --refHaps ${reference} \
        --haps ${vcf} \
        --prefix ${out_prefix} \
        --ChunkLengthMb  300 \
        --ChunkOverlapMb 100 \
        --cpus ${task.cpus}
    """
}
