# SNP-calling-GATK4

This repository contains an annotated bash script (SNP_calling_GATK4.sh) which outlines all the steps (see below) required for the bootstrapping workflow to conduct variant calling (.vcf) on contigs (here UCE).

## STEPS
1 Data preparation - get reference, index it and generate sequence dictionnary  
2 Map cleaned reads to reference, sort, mark and remove duplicates & index bams - SCRIPT 1_bwa_mapping_array.pbs (array job)  
3 Initial haplotype calling & consolidation - SCRIPT 2_haplotype_cohort.pbs (array job)  
4 Joint genotyping & snp extraction - SCRIPT 3_GVCF1.pbs  
5 Hard filtering for snps from haplotype caller callset - SCRIPT 4_HardFiltering.pbs  
6 1st recalibration on unrecalibrated data & before/after plots- SCRIPT 5_BQSR_recal1.pbs  
7 Haplotype calling on 1st recalibrated bam files & consolidation - SCRIPT 6_haplotype_cohort2.pbs (array job)  
8 Joint genotyping & snp extraction - SCRIPT 7_GVCF2.pbs  
9 Hard filtering - SCRIPT 8_HardFiltering.pbs  
10 2nd recalibration on unrecalibrated data & before/after plots - SCRIPT 6_BQSR_recal2.pbs  
######### CHECK FOR CONVERGENCE IN THE BEFORE/AFTER RECALIBRATION PLOTS ########  
If convergence - 11a "true" variant calling on recalibrated bam files - SCRIPT 10_final.pbs  
If no convergence - 11b Repeat steps from 7-10 until convergence (usually 2-4 times is enough)  

Each step is associated with a pbs script and was set to run on the JCU HPC machine. Note that some are best run as array jobs (e.g. mapping, haplotype calling and BQSR)

## REQUIRED CONTAINERS TO RUN THE WORKFLOW versions in parentheses were used in provided scripts
bwa (v.0.7.17)  
samtools (v.1.12)  
gatk4 (v.4.2.0.0)  
R (v.4.0.3)  

## CITATION
Baraf, L., Hung, J. and Cowman, P. (in press) Phylogenomics of marine angelfishes: diagnosing sources of systematic discordance for an iconic reef fish family (F: Pomacanthidae)."

## CONTACT
If you have any questions or issues running this worfkflow - feel free to reach out!
lauriane.baraf@my.jcu.edu.au

## BEST FISHES AND HAPPY CODING !



