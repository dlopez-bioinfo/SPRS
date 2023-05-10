#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2


/*
================================================================
  MODULES
================================================================
*/
include {VCF_PHASE} from "../modules/vcf_phase"
include {VCF_INDEX} from "../modules/vcf_index"
include {VCF_CONCAT} from "../modules/vcf_concat"
include {VCF_SORT} from "../modules/vcf_sort"
include {VCF_IMPUTE} from "../modules/vcf_impute"


/*
================================================================
  MAIN WORKFLOW
================================================================
*/
workflow SUB_IMPUTE_VCF {
  take:
    vcf_ch
    shapeit4_map
    minimac4_ref

  main:
    // Create a channel for shapeit4 genetic map
    // [ CHROM, MAP_FILE, IS_PAR ]
    channel
      .fromPath(shapeit4_map)
      .map { it -> 
        [ (it.name =~ /^chr([0-9X]+).*/)[0][1], it, (it =~ /X_par[12]/).getCount()>0 ]}
      .dump(tag: "__SHAPEIT4_MAP__")
      .set { shapeit4_map_ch }


    // Create a channel for shapeit4 genetic map
    // [ CHROM, MAP_FILE, IS_PAR ]
    channel
      .fromPath(minimac4_ref)
      .map { it ->
        [ (it.name =~ /^([0-9X]+).*/)[0][1], it, (it =~ /X.Pseudo.Auto/).getCount()>0 ]}
      .dump(tag: "__MINIMAC4_REF__")
      .set { minimac4_ref_ch }


    // Create a channel for phasing VCF
    // [ CHROM, VCF_FILE, VCF_INDEX, SHAPEIT4_MAP, is_PAR ]
    vcf_ch
      .cross(shapeit4_map_ch)
      .map { [ 
        it[0][0],  // Chromosome
        it[0][1],  // VCF file
        it[0][2],  // VCF index file
        it[1][1],  // Shapeit4 genomic map file
        it[1][2]]} // Is X_PAR?
      .dump(tag: "__PHASING_CH__")
      .set { phase_ch }


    // phase VCF
    VCF_PHASE(phase_ch) 


    // concat X pseudo-autosome regions (PAR1 and PAR2)
    // [ CHROM, VCF_FILE, is_PAR ]
    VCF_PHASE.out
      .filter { chrom, vcf, is_par -> 
        is_par == true }   // select X_PAR VCFs
      .map { chrom, phased_vcf, is_par -> 
        [ chrom, phased_vcf ]} \
      | VCF_INDEX \
      | groupTuple
      | map { chrom, par1, par2 -> [ 
          chrom,
          [ par1, par2 ].flatten() ]}   // list of VCFs and corresponding indexes  
      | VCF_CONCAT \
      | map { chrom, vcf -> 
        [ chrom, vcf, true ] } \
      | set { pseudo_autosomes_ch }


    // create channel for autosomes and pseudo-autosomes (PAR)
    minimac4_ref_ch                 //add autosomes
      .filter { chrom, ref, is_par -> 
        is_par == false } 
      .join(VCF_PHASE.out
        .filter { chrom, vcf, is_par -> 
          is_par == false })       
      .map { chrom, ref, is_par, vcf, _is_par -> 
        [ chrom, ref, vcf ]}
      .mix(minimac4_ref_ch          // add (merged) pseudo-autosomes
        .filter { chrom, ref, is_par -> 
          is_par == true } 
        .join(pseudo_autosomes_ch)  
        .map { chrom, ref, is_par, vcf, _is_par -> 
          [ chrom, ref, vcf ]})
      .set { impute_ch}

    // Impute VCFs
    VCF_IMPUTE(impute_ch)

  emit:
    VCF_IMPUTE.out //imputed VCFs

}
