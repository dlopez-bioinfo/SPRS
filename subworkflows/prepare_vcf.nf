#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2


/*
================================================================
  MODULES
================================================================
*/
include {VCF_REMOVE_CHR_PREFIX} from "../modules/vcf_remove_chr_prefix"
include {VCF_REMOVE_AF} from "../modules/vcf_remove_af"
include {VCF_FILTER_PASS} from "../modules/vcf_filter_pass"
include {VCF_SORT} from "../modules/vcf_sort"
include {
  VCF_INDEX as INDEX_FILTERED;
  VCF_INDEX as INDEX_MERGED} from "../modules/vcf_index"
include {VCF_MERGE_BY_CHROM} from "../modules/vcf_merge_by_chrom"
include {VCF_NORM} from "../modules/vcf_norm"


/*
================================================================
  MAIN WORKFLOW
================================================================
*/
workflow SUB_PREPARE_VCF {
  take:
    vcf_ch
    exome
    capture_bed
    mapping_file


  main:
    // Make a Channel with chromosome names
    // [ CHROM ]
    Channel
      .of(1..22, 'X')
      .map { it.toString()}
      .set { chrom_ch }


    // Remove "chr" prefix and clean some problematic fields
    VCF_REMOVE_CHR_PREFIX(vcf_ch, mapping_file) \
      | VCF_REMOVE_AF 
    
    
    // keep only pass variant
    VCF_FILTER_PASS (VCF_REMOVE_AF.out, exome, capture_bed) \
      | VCF_SORT \
      | INDEX_FILTERED
    return_ch = INDEX_FILTERED.out  


    // Prepare a channel with all the VCFs
    INDEX_FILTERED.out
      .map { sampleId, vcf, idx -> 
        [ vcf, idx] } 
      .collect() \
      .set { prepared_vcf_idx_list }


    // Merge VCFs by chromosome
    VCF_MERGE_BY_CHROM(prepared_vcf_idx_list, chrom_ch) \
      | map { chrom, vcf ->
          [ chrom, vcf, params.reference_genome ]} \
      | VCF_NORM \
      | INDEX_MERGED
    return_ch = INDEX_MERGED.out  

  emit:
    prepared_vcf = return_ch 
}
