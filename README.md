# Spanish Polygenic Risk Score reference distributions pipeline

## Table of contents
  * [Summary](#summary)
  * [Requisites](#requisites)
  * [Before you start](#before-you-start)
  * [Usage](#usage)
  * [Citation](#citation)

-------------------------------------------------------------

## Summary
*sprs* is a bioinformatics pipeline for calculating Polygenic Risk Scores (PRS) from [PDG catalog](https://www.pgscatalog.org) on a set of samples so that they can be compared to the [Spanish PRS reference distributions](http://csvs.clinbioinfosspa.es/?tab=prs). The pipeline was designed to fix common VCF malformations, impute missing values and eventually generate a sample sheet to feed the [pgsc_calc](https://github.com/PGScatalog/pgsc_calc) tool.

The pipeline is built using Nextflow, a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The Nextflow DSL2 implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies.


## Requisites
1. Software dependencies
    * [Nextflow](https://www.nextflow.io) >= 22.10.6
    * [Singularity](https://apptainer.org/docs/user/latest/) or [Docker](https://docs.docker.com/) containerization software 
2. Resources
    * **Reference genome**: we recommend [hs37d5](ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/), that includes data from GRCh37, the rCRS mitochondrial sequence, Human herpesvirus 4 type 1 and the concatenated decoy sequences. 
    * **Shapeit4 genetic maps**: maps for b37 can be downloaded directly from [shapeit4 github repository](https://github.com/odelaneau/shapeit4/blob/master/maps/genetic_maps.b37.tar.gz)
    * **Mimimac4 reference panel**: we recommend [1000 Genomes Phase 3 reference panel](ftp://share.sph.umich.edu/minimac3/G1K_P3_M3VCF_FILES_WITH_ESTIMATES.tar.gz)
3. Samples
    * For an accurate phasing prediction step, it is recommended to analyze 20 or more samples at the same time. 
    * VCFs are expected to be in hg19 or GRCh37. GRCh38 is not supported yet. If needed, [picard liftover tool](https://gatk.broadinstitute.org/hc/en-us/articles/360036831351-LiftoverVcf-Picard-) can be used to remap VCFs.


## Before you start
*nextflow.config* included in this workflow has profiles for *singularity* and *docker*, and can be executed *locally* or in a standard *Slurm* HPC. Other executor engines like *SGE*, *AWS*, *Google cloud* or *kubernetes* can easily be used by doing some minor changes to the *nextflow.config*. Although it is not mandatory, we recommend you to create a *profile* that fit the specific requirements of your institution. See [nextflow executors](https://www.nextflow.io/docs/latest/executor.html) for additional information.


## Usage
Firstly, the sample-sheet and the imputed VCFs files have to be generated. You can do this for a set of genomes in a slurm environment by running the following command:

```bash
$ nextflow run babelomics/sprs \
  -profile slurm,singularity \
  -resume \
  --input_folder input_folder \
  --exome 'false' \
  --shapeit4_map_files '/resources/shapeit4/chr*.b37.gmap.gz' \
  --minimac4_ref_files '/resources/minimac4/*.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz' \
  --reference_genome '/resources/ref/hs37d5.fa'
```

To see the full pipeline parameters list use:

```bash
$ nextflow run babelomics/sprs --help 
```

Note that *-profile* and *-resume* are nextflow parameters and therefore are preceded by a single hyphens. To run the pipeline locally on a docker environment use `-profile docker`. Custom configuration files can be included with `-config nextflow.config` option. To see additional nextflow parameters use `nextflow run -help`.


Once the sample sheet has been generated, you can calculate any polygenic risk score from [PDG catalog](https://www.pgscatalog.org) by running:

```bash
$ nextflow run pgscatalog/pgsc_calc \
  -profile docker \
  -resume \
  --input results/psc_calc_samplesheet.txt \
  --target_build GRCh37 \
  --pgs_id PGS000021
```

Your results, could then be compared with the [Spanish PRS reference distribution](http://csvs.clinbioinfosspa.es/?tab=prs).


## Citation
A manuscript describing the tool is in preparation. In the meantime if you use the tool we ask you to cite the repo and the paper describing the CSVS resource:
  * Spanish Polygenic Risk Score pipeline. https://github.com/babelomics/sprs
  * Peña-Chilet M et al. CSVS, a crowdsourcing database of the Spanish population genetic variability. Nucleic Acids Res. 2021 Jan 8;49(D1):D1130-D1137. doi: 10.1093/nar/gkaa794.


