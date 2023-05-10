#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process VCF_PHASE {
  label 'parallel'
  container "${params.container__shapeit4}"
  publishDir "${params.output_folder}/vcf/", mode: 'copy', overwrite: true
                  
  input:
    tuple val(chrom), 
          path(vcf),      
          path(vcf_idx),    
          path(genome_map),
          val(is_PAR)

  output:
    tuple val(chrom), 
          path("${vcf.name.replaceAll(/.vcf.gz/, '')}_phased_${genome_map.name.replaceAll(/.gmap.gz/, '')}.vcf.gz"),
          val(is_PAR)

  script:
    def out = "${vcf.name.replaceAll(/.vcf.gz/, '')}_phased_${genome_map.name.replaceAll(/.gmap.gz/, '')}.vcf.gz"    

    """
      shapeit4 \
        --thread ${task.cpus} \
        --input ${vcf} \
        --map ${genome_map} \
        --sequencing \
        --region ${chrom} \
        --output ${out}          
    """
}
