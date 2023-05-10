#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2


/*
================================================================
  MODULES
================================================================
*/
include {VCF_CONCAT} from "../modules/vcf_concat"
include {VCF_SORT} from "../modules/vcf_sort"


/*
================================================================
  MAIN WORKFLOW
================================================================
*/
workflow SUB_CREATE_SAMPLE_SHEET {
  take:
    vcf_ch

  main:
    // Make a channel of autosomes and pseudoautosomes
    // [ CHROM, VCF]
    vcf_ch
    .branch { chrom, vcf -> 
      chrom_sex: chrom == 'X' 
      chrom_auto: chrom != 'X'}
    .set { vcf_final_ch }


    // Merge X (PAR1, PAR2 and non-pseudoautosomal region)
    VCF_CONCAT(vcf_final_ch.chrom_sex) \
      | VCF_SORT


    // write pgsc_calc pipeline sample sheet
    ss_path = "${params.output_folder}/psc_calc_samplesheet.txt"
    sample_sheet = file(ss_path)
    sample_sheet.text = "sampleset,vcf_path,bfile_path,pfile_path,chrom\n"
    vcf_final_ch.chrom_auto
      .mix(VCF_SORT.out)
      .map { chrom, vcf -> 
        sample_sheet.append("CSVS,${vcf},,,${chrom}\n")
      }
}
