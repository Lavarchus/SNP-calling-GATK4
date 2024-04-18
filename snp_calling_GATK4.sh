#! /bin/bash

## This script was modified from Erickson, K. L., Pentico, A., Quattrini, A. M., & McFadden, C. S. (2021). New approaches to species delimitation and population structure of anthozoans: Two case studies of octocorals using ultraconserved elements and exons. Molecular Ecology Resources, 21(1), 78-92
## It runs with GATK-4.2.* and aims to call variants from UCE data of non-model species that lacks "gold standard" snps or indels reference.
## Steps and required programs are listed below. Scripts for each step are available in the associated repository. 
# NOTE: Option for calling ONLY snps and removing those around or within indels by masking indels is outlined at hard filtering steps (5 & 9). This might be useful if you want to run a e.g. SNAPP analysis which requires an input free of indels. You can also remove indels with vcftools after running your GATK bootstrapping.

## Written by Lauriane Baraf, lauriane.baraf@my.jcu.edu.au on the 15th of September 2022.
# Modified last: 05/10/2022

## STEPS
# 1 Data preparation - get reference, index it and generate sequence dictionnary
# 2 Map cleaned reads to reference, sort, mark and remove duplicates & index bams - SCRIPT 1_bwa_mapping_array.pbs (array job)
# 3 Initial haplotype calling & consolidation - SCRIPT 2_haplotype_cohort.pbs (array job)
# 4 Joint genotyping & snp extraction - SCRIPT 3_GVCF1.pbs
# 5 Hard filtering for snps from haplotype caller callset - SCRIPT 4_HardFiltering.pbs
# 6 1st recalibration on unrecalibrated data & before/after plots- SCRIPT 5_BQSR_recal1.pbs
# 7 Haplotype calling on 1st recalibrated bam files & consolidation - SCRIPT 6_haplotype_cohort2.pbs (array job)
# 8 Joint genotyping & snp extraction - SCRIPT 7_GVCF2.pbs
# 9 Hard filtering - SCRIPT 8_HardFiltering.pbs
# 10 2nd recalibration on unrecalibrated data & before/after plots - SCRIPT 6_BQSR_recal2.pbs
######### CHECK FOR CONVERGENCE IN THE BEFORE/AFTER RECALIBRATION PLOTS ########
## IF CONVERGENCE
# 11a "true" variant calling on recalibrated bam files - SCRIPT 10_final.pbs
## IF NO CONVERGENCE
# 11b Repeat steps from 7-10 until convergence (usually 2-4 times is enough)

## Load containers ###
module load conda3
source $CONDA_PROF/conda.sh
shopt -s expand_aliases
conda activate bwa-0.7.17
conda activate --stack samtools-1.12
conda activate --stack gatk4-4.2.0.0
conda activate --stack R-4.0.3

## If you want to have your output directories ready - set for 2 runs, add more if more recal needed to reach convergence
mkdir -p snp_calling/{reference,logs,1_mapping_output,2_haplotype_cohort_output,3_GVCF,4_hardfiltering_output,5_BQSR_recal1_output,6_haplotype_cohort_output,7_GVCF_output,8_harfiltering_output,9_BQSR_recal2_output}

################################################ 1 DATA PREPARATION #############################################
# Run the summary stats command from Phyluce on the exploded fasta files (phyluce_assembly_get_fasta_lengths) to get SNP calling reference as the sample with the highest number of UCE contigs in the taxon set. In the output, number of contigs is the second column from the fasta summary file.
# Which sample has most contigs?
awk -F',' '{sub(/\.unaligned.fasta/, "", $1); print $1","$2}' /path/to/fasta_summary.csv | sort -t, -k2nr

# Store unaligned contigs in reference directory and set your reference variable for your sample
REFERENCE="path/to/reference/reference.unaligned.fasta"
# Index reference loci
bwa index ${REFERENCE} #OR bwa index -p prefix -a is /path/to/ref.fasta
# Generate fasta file index & a sequence dictionary for HaplotypeCaller
samtools faidx ${REFERENCE}
# Create sequence dictionary so header only contains sequence records
java -Xmx2g -jar /sw/picard/2.18.1/picard.jar CreateSequenceDictionary R=${REFERENCE} O=./reference.unaligned.dict

################### 2 MAP READS TO REFERENCE & MARK DUPLICATES - RUN 1 BWA MAPPING ARRAY PBS SCRIPT #################
# Might want to make a directory for log files
cd /path/to/snp_calling
mkdir -p logs/1_mapping_logs

# Run mapping array job 1_bwa_mapping_array.pbs - modify according to number of samples (e.g. #PBS -J 0-90 for array index)

# OPTIONAL AFTER MAPPING: Get statistics
    #Mapping statistics
METRICS="$PWD/1_mapping_output/*dedup_metricsfile"
# Create file with header
echo "LIBRARY UNPAIRED_READS_EXAMINED READ_PAIRS_EXAMINED UNMAPPED_READS	UNPAIRED_READ_DUPLICATES READ_PAIR_DUPLICATES READ_PAIR_OPTICAL_DUPLICATES PERCENT_DUPLICATION ESTIMATED_LIBRARY_SIZE" > metrics_stats.txt

for docs in $METRICS;
    do
    echo $docs;
    #create a file name based on document's name. Get path to folder and only keep last field
    FILE_NAME=$(echo $docs | cut -d/ -f7);
    echo $FILE_NAME
    #create sample name based on file's name
    SAMPLE_NAME=$(echo $FILE_NAME | cut -d. -f1);
    #print line with summary stats
    SUMMARY=$(sed -n '8p' $docs)
    echo "$SAMPLE_NAME$SUMMARY" >> metrics_stats.txt
done

    # SAMTOOLS statistics
STATS="$PWD/1_mapping_output/*stats.txt"
# Create file with header
echo "LIBRARY Total_QC_passed_reads duplicates MAPPED PAIRED read1 read2 properly_paired with_itself_mate_paired singletons	mate_mapped_different_chr mate_mapped_different_chrQ5" > samtools_stats.txt

for docs in $STATS;
    do
    echo $docs;
    #create file name based on document's name. Get path to folder and only keep last field
    FILE_NAME=$(echo $docs | cut -d/ -f7);
    echo $FILE_NAME
    #create sample name based on file's name
    SAMPLE_NAME=$(echo $FILE_NAME | cut -d. -f1);
    echo $SAMPLE_NAME
    #print line with summary stats
    S1=$(sed -n '1p' $docs)
    S2=$(sed -n '2p' $docs)
    S3=$(sed -n '3p' $docs)
    S4=$(sed -n '4p' $docs)
    S5=$(sed -n '5p' $docs)
    S6=$(sed -n '6p' $docs)
    S7=$(sed -n '7p' $docs)
    S8=$(sed -n '8p' $docs)
    S9=$(sed -n '9p' $docs)
    S10=$(sed -n '10p' $docs)
    S11=$(sed -n '11p' $docs)

    echo "$SAMPLE_NAME $S1 $S2 $S3 $S4 $S5 $S6 $S7 $S8 $S9 $S10 $S11" >> samtools_stats.txt

done

################## 3 FIRST ROUND OF HAPLOTYPE CALLING - CONSOLIDATION - RUN 2_HAPLOTYPE_COHORT PBS SCRIPT ##################
# Might want to make a directory for log files
mkdir -p logs/2_haplotype_cohort_logs

# Initial variant calling on mapped reads, probably the most time consumming step. All parameters were defined as per Erickson et al. (2021).
gatk HaplotypeCaller \
    -R ${REFERENCE} \
    -I ${SAMPLE}.ref_all_dedup.bam \
    -O 2_haplotype_cohort_output/${OUTPUT_NAME} \
    --emit-ref-confidence GVCF \
    --contamination-fraction-to-filter 0.0002 \
    --min-base-quality-score 20 \
    --phred-scaled-global-read-mismapping-rate 30 \
    --standard-min-confidence-threshold-for-calling 40.0

# Consolidate/Combine GVCFs files into one vcf file that can be used for joint genotyping (t=15mins)
ls -d -1 2_haplotype_cohort_output/*.g.vcf > gvcf.list

gatk CombineGVCFs \
    -R ${REFERENCE}
    --variant 2_haplotype_cohort_output/gvcf.list \
    -O 2_haplotype_cohort_output/combined_gvcf.vcf

############################## 4 JOINT GENOTYPING & EXTRACTION - RUN 3_GVCF1 PBS SCRIPT #############################
# Joint genotyping to determine the genotype of each individual at each site.
gatk GenotypeGVCFs \
    -R ${REFERENCE} \
    -V 2_haplotype_cohort_output/combined_gvcf.vcf \
    -O 3_GVCF_output/genotyped_X_samples.vcf

#Extract SNPSs from callset
gatk SelectVariants \
    -R ${REFERENCE} \
    -V 3_GVCF_output/genotyped_X_samples.vcf \
    -select-type SNP \
    -O 3_GVCF_output/genotyped_X_samples_snps.vcf

# Extract indels from callset 
gatk SelectVariants \
    -R ${REFERENCE} \
    -V 3_GVCF_output/genotyped_X_samples.vcf \
    -select-type INDEL \
    -O 3_GVCF_output/genotyped_X_samples_indels.vcf

############ OPTIONAL: GET TABLE FOR DENSITY PLOTS ############
# Make density plots for each annotations to define the best hardfiltering values
# Convert vcf file into a table that can be used in R
gatk VariantsToTable \
    -V 3_GVCF1_output/genotyped_X_samples_snps.vcf \
    -R ${REFERENCE} \
    -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
    --output gatk4_density_plot.table
# Use the outputted table to plot variant scores and get a better idea on how to define your thinning values for the variant hard filtration 
# Usually all next-gen seq runs have similar distribution patterns so the filtering values recommended by gatk are a good start, I'd always double check for every datasets, just so I can sleep at night.
################################################################

################## 5 HARD FILTERING - RUN 4_HARDFILTERING PBS SCRIPT ##################
# Haplotype caller algorithm is purposefully not too conservative so you want to pass the variants through a harder filter on variants to improve callset. This will help distinguish true variants from false variants that can emerge from e.g. sequencing errors as basecall noise. By default HaplotypeCaller and GenotypeGVCF do not emit variants with QUAL < 10.0.
# Filtering values were determine based on the results from the density plots. Not too far from "standard" gatk hardfiltering values so probably would be okay to just run with it on another next-gen seq data, although I'd always recommend to double check (reviewer 2 is always watching...). 
# Note: if working on exons or uces from invertebrates might be better to use the SOR filter, which is an alternative to estimating strand bias (like FS) because FS tend to penalize variants occuring at the ends of exons and uces are mainly exonic in inverts. Reads at the ends of exons tend to only be covered by reads in one direction and FS gives those variants a bad score so adding SOR will take into account the ratios of reads that cover both alleles. SOR is also better for SOR better for high coverage data.

# OPTION: Masking the indels from snps data will just flag them as such in the filtered callset, doing that will allow their removal when running the SelectVariant command --exclude-filtered which will only let the "PASSED" data filter through OR when preparing the input for e.g. SNAPP analysis. This is useful to identify SNPs that are contained within or around indels.
# NOTE: You could also run all the filtration on indels first and then use the hard filtered (bqsr) vcf file to mask the indels within the snp callset. Indel calling methods have high rates of false positives so I tried running it with the hard filtered indel vcf file to account for this and avoid loosing too much snps to false positive indels. Out of the 38,608 indels extracted from the initial callset, 3,385 were filtered out after hard filtering of which 34 were later identified as PASSED SNPs. There was no notable differences in the b/a recalibration plots so it's probably okay to mask with the "light" filtered indel file if only gaining a few SNPs.

gatk VariantFiltration \
    -R ${REFERENCE}  \
    -V 3_GVCF_output/genotyped_X_samples_snps.vcf \
    --mask 3_GVCF_output/genotyped_X_samples_indels.vcf \
    --mask-extension 5 \
    --mask-name InDel \
    --cluster-window-size 10 \
    --filter-expression "QUAL < 30.0" \
    --filter-name "LowQuality" \
    --filter-expression "QD < 4.0" \
    --filter-name "QualByDepth" \
    --filter-expression "MQ < 38.0" \
    --filter-name "RMSMappingQuality" \
    --filter-expression "FS > 60.0" \
    --filter-name "FisherStrand" \
    -O 4_hardfiltering_output/genotyped_X_samples_snps_filtered.vcf
# --filter-expression "SOR > 15.0" \
# --filter-name "StrandOddsRatio" \

# Filter variants according to filter expression: INDEL. Filter values were set according to GATK guidelines.
gatk VariantFiltration \
	-R ${REFERENCE} \
	-V 3_GVCF_output/genotyped_X_samples_indels.vcf \
	-filter-name "QD_filter" \
  -filter-expression "QD < 2.0" \
	-filter-name "FisherStrand_filter" \
  -filter-expression "FS > 200.0" \
	-filter-name "LowQual_filter" \
  -filter-expression "QUAL < 30.0" \
  -O 4_hardfiltering_output/genotyped_X_samples_indels_filtered.vcf
# -filter-name "SOR_filter" \
# -filter-expression "SOR > 10.0" \

# Check number of record pre and post filtration, should be the same because it's only labelling snp types, you can see how many variants will be left in the callset after the hard filtration by looking at the "PASS" ones in the next step
grep -v "^#" 3_GVCF_output/genotyped_X_samples_snps.vcf | wc -l 
grep -v "^#" 4_hardfiltering_output/genotyped_X_samples_snps_filtered_1.vcf | wc -l 
# And number of variants filtered out by checking the PASS column
grep -v "^#" 4_hardfiltering_output/genotyped_X_samples_snps_filtered_1.vcf | cut -f 7 | sort | uniq -c  

# Exclude filtered variants to get only the variants that PASSED the filtering as the input for the BaseRecalibration (BQSR) tool later on because, like many GATK tools, BaseRecalibration (BQSR) most likely doesn't ignore the "FAILED" status of variants. You can also only keep biallelic sites here by using --restrict-alleles-to BIALLELIC or do it after GATK using VCFtools.
gatk SelectVariants \
    --exclude-filtered true \
    -V 4_hardfiltering_output/genotyped_X_samples_snps_filtered.vcf \
    -O 4_hardfiltering_output/bqsr_snps.vcf

gatk SelectVariants \
	--exclude-filtered true \
	-V 4_hardfiltering_output/genotyped_X_samples_indels_filtered.vcf \
	-O 4_hardfiltering_output/bqsr_indels.vcf

################## 6 BASE RECALIBRATION - APPLYBSQR - A/F RECALIBRATION PLOTS LOOP ##################
## 1ST RECALIBRATION - RUN SCRIPT 5_BQSR_recal1.pbs 
# Run BaseRecalibrator on the orignial bam files to get the recal table with the original base scores. If did indel masking, comment out --known-sites bqsr_indels.vcf
gatk BaseRecalibrator \
    -R ${REFERENCE} \
    -I 1_mapping_output/sample.ref_all_dedup.bam \
    --known-sites 4_hardfiltering_output/bqsr_snps.vcf \
    --known-sites 4_hardfiltering_output/bqsr_indels.vcf \
    -O 5_BQSR_recal1_output/pre_recal1.table

# Run ApplyBQSR with the recal table from last step to create a new bam file with the adjusted base quality scores
gatk ApplyBQSR \
    -R ${REFERENCE} \
    -I 1_mapping_output/4{SAMPLE}.ref_all_dedup.bam \
    --bqsr-recal-file 5_BQSR_recal1_output/pre_recal1.table \
    -O 5_BQSR_recal1_output/${SAMPLE}_recal1.bams

# Run BaseRecalibrator again but on the adjusted bam files from last step to create a new recal table with the adjusted base quality scores and do a B/A plot. If did indel masking, comment out --known-sites bqsr_indels.vcf.
gatk BaseRecalibrator \
    -R ${REFERENCE} \
    -I 5_BQSR_recal1_output/recal1.bams \
    -O 5_BQSR_recal1_output/post_recal1.table \
    --known-sites bqsr_snps.vcf \
    --known-sites bqsr_indels.vcf

## BEFORE/AFTER RECALIBRATION PLOTS ## INCLUDED IN PREV SCRIPT 5_BQSR_RECAL1.PBS
# Before/After recalibration plots to see how the first recalibration worked out on the original data. You want to look for convergence. If it converges you can go ahead into a "true" variant calling run and call on the recalibrated bam files. If results aren't converging, you'll need to bootstrap your data by repeating all the steps from HaplotypeCaller on the unrecalibrated data to another recalibration on the unrecalibrated data and look at the plots for convergence (usually 2-4 times).

# R packages that needs to be installed for AnalyzeCovariates: gsalib & ggplot2
gatk AnalyzeCovariates \
    -before ${SAMPLE_NAME}_pre_recal1.table \
    -after ${SAMPLE_NAME}_post_recal1.table \
    -plots ${SAMPLE_NAME}_recalibration_plot_recal1.pdf

# You should see that the quality score accuracy is better after recalibration. To get convergence of your covariates you'll need at least a second run of recalibration.

############ 7 Call variant with HaplotypeCaller on recalibrated data - RUN SCRIPT 6_haplotype_cohort2.pbs ##########
gatk HaplotypeCaller \
    -R ${REFERENCE} \
    -I 5_BQSR_recal1_output/${SAMPLE}_recal1.bams \
    -O 6_haplotype_cohort_output/${OUTPUT_NAME} \
    --emit-ref-confidence GVCF \
    --contamination-fraction-to-filter 0.0002 \
    --min-base-quality-score 20 \
    --phred-scaled-global-read-mismapping-rate 30 \
    --standard-min-confidence-threshold-for-calling 40.0

# Consolidate/Combine GVCFs files to create the file for joint genotyping (t=15mins)
ls -d -1 6_haplotype_cohort_output/*.g.vcf > gvcf_recal1.list

gatk CombineGVCFs \
    -R ${REFERENCE} \
    -V gvcf_recal1.list \
    -O 6_haplotype_cohort_output/combined_gvcf_recal1.vcf

####### 8 Joint genotyping with cohort of GVCF variant files produced by HaplotypeCaller - RUN SCRIPT 7_GVCF.PBS ####### (t=10mins)
gatk GenotypeGVCFs \
    -R ${REFERENCE} \
    -V 6_haplotype_cohort_output/combined_gvcf_recal1.vcf \
    -O 7_GVCF_output/genotyped_X_samples_recal1.vcf

#Extract SNPSs from callset
gatk SelectVariants \
	-R ${REFERENCE} \
	-V 7_GVCF_output/genotyped_X_samples_recal1.vcf \
	-select-type SNP \
	-O 7_GVCF_output/genotyped_X_samples_snps_recal1.vcf &&

# Extract indels from callset
gatk SelectVariants \
	-R ${REFERENCE} \
	-V 7_GVCF_output/genotyped_X_samples_recal1.vcf \
	-select-type INDEL \
	-O 7_GVCF_output/genotyped_X_samples_indels_recal1.vcf

########################## 9 Hard filtering - RUN SCRIPT 8_HardFiltering.PBS ########################## (t=1sec)
# OPTION: Masking the indels from snps data will just flag them as such in the filtered callset, doing that will allow their removal when running the SelectVariant command --exclude-filtered which will only let the "PASSED" data filter through OR when preparing the input for e.g. SNAPP analysis (see above step 5 for more information)
gatk VariantFiltration \
    -R ${REFERENCE}  \
    -V 7_GVCF2_output/genotyped_X_samples_snps_recal1.vcf \
    -O 8_hardfiltering_output/genotyped_X_samples_snps_filtered_recal1.vcf \
    --mask 7_GVCF_output/genotyped_X_samples_indels_recal1.vcf \
    --mask-extension 5 \
    --mask-name InDel \
    --cluster-window-size 10 \
    --filter-expression "QUAL < 30.0" \
    --filter-name "LowQuality" \
    --filter-expression "QD < 5.0" \
    --filter-name "LowCoverage" \
    --filter-expression "MQ < 35.0" \
    --filter-name "RMSMappingQuality" \
    --filter-expression "FS > 60.0" \
    --filter-name "FisherStrand" 
# --filter-expression "SOR > 4.0" \
# --filter-name "StrandOddsRatio" 

# Filter variants according to filter expression: INDEL.
gatk VariantFiltration \
	-R ${REFERENCE} \
	-V 7_GVCF2_output/genotyped_X_samples_indels_recal1.vcf \
  -O 8_hardfiltering_output/genotyped_X_samples_indels_filtered_recal1.vcf \
	-filter-name "QD_filter" \
  -filter-expression "QD < 2.0" \
	-filter-name "FisherStrand_filter" \
  -filter-expression "FS > 200.0" \
	-filter-name "LowQual_filter" \
  -filter-expression "QUAL < 30.0" 
# -filter-name "SOR_filter" \
# -filter-expression "SOR > 10.0"

# Check number of record post filtration, should be the same
grep -v "^#" 8_hardfiltering_output/genotyped_X_samples_snps_filtered_recal1.vcf | wc -l 
# # And number of variants filtered out by checking the PASS column
grep -v "^#" 8_hardfiltering_output/genotyped_X_samples_snps_filtered_recal1.vcf | cut -f7 | sort | uniq -c 

# Exclude filtered variants to get only the variants that PASSED the filtering as the input for the BaseRecalibration (BQSR) tool.
gatk SelectVariants \
    --exclude-filtered true \
    -V 8_hardfiltering_output/genotyped_X_samples_snps_filtered_recal1.vcf \
    -O 8_hardfiltering_output/bqsr_snps_recal1.vcf &&

gatk SelectVariants \
    --exclude-filtered true \
    -V 8_hardfiltering_output/genotyped_X_samples_indels_filtered_recal1.vcf \
    -O 8_hardfiltering_output/bqsr_indels_recal1.vcf

################# 10 2ND RECALIBRATION ON UNCALIBRATED BAMS - RUN SCRIPT 9_BQSR_recal2.PBS #################
# Get recalibration table after 2nd pass on unrecalibrated data. If did indel masking, comment out --known-sites bqsr_indels_recal1.vcf
gatk BaseRecalibrator \
    -R ${REFERENCE} \
    -I ${SAMPLE} \
    --known-sites 8_hardfiltering_output/bqsr_snps_recal1.vcf \
    --known-sites 8_hardfiltering_output/bqsr_indels_recal1.vcf \
    -O 9_BQSR_recal2_output/${SAMPLE_NAME}_recal2.table

# Run ApplyBQSR with the recal table from last step to create a new bam file with the adjusted base quality scores
gatk ApplyBQSR \
    -R ${REFERENCE} \
    -I ${SAMPLE} \
    --bqsr-recal-file 9_BQSR_recal2_output/${SAMPLE_NAME}_recal2.table \
    -O 9_BQSR_recal2_output/${SAMPLE_NAME}_recal2.bams

# Run BaseRecalibrator on the recalibrated bam files to get the table with corrected base scores after 2nd recalibration to plot before/after. If did indel masking, comment out --known-sites bqsr_indels_recal1.vcf
gatk BaseRecalibrator \
    -R ${REFERENCE} \
    -I 9_BQSR_recal2_output/${SAMPLE_NAME}_recal2.bams \
    --known-sites 8_hardfiltering_output/bqsr_snps_recal1.vcf \
    --known-sites 8_hardfiltering_output/bqsr_indels_recal1.vcf \
    -O ${SAMPLE_NAME}_post_recal2.table

# Before/After recalibration plots to see how the first recalibration worked out on the original data. You want to look for convergence. If it converges you can go ahead into a "true" variant calling run and call on the recalibrated bam files. If results aren't converging, you'll need to bootstrap your data by repeating all the steps from HaplotypeCaller on the unrecalibrated data until the plots show convergence (usually 2-4 times).

# Before/After plot comparing 1st and 2nd recalibrations
gatk AnalyzeCovariates \
    -before 5_BQSR_recal1_output/${SAMPLE_NAME}_post_recal1.table \
    -after ${SAMPLE_NAME}_post_recal2.table \
    -plots ${SAMPLE_NAME}_recalibration_plot_recal2.pdf

# if your covariates haven't converged yet, no need to break the duck's third leg (weird french expression to say, don't panic! odd I know..) just keep repeating from steps 7-10 until convergence (usually 2-4 times is enough). May the covariate gods be with you!

################### 11 Final Haplotype Calling on recalibrated data | JOINT GENOTYPING | HARD FILTERING ###################
# RUN SCRIPT 10_final.PBS

# Check number of record pre and post filtration, should be the same
grep -v "^#" $PWD/10b_GVCF_output/genotyped_X_samples_snps_recal2.vcf | wc -l 
grep -v "^#" $PWD/10c_hardfiltering_output/genotyped_X_samples_snps_filtered_recal2.vcf| wc -l 
# And number of variants filtered out by checking the PASS column
grep -v "^#" $PWD/10c_hardfiltering_output/genotyped_X_samples_snps_filtered_recal2.vcf | cut -f 7 | sort | uniq -c 

# Compile statistics with BCFtools
conda activate bcftools-1.12
# Get info on final callset
# Get number of sites
bcftools view -H $PWD/10c_hardfiltering_output/genotyped_PASS_snps_raw_recal2.vcf | wc -l
# or grep -v "^#" $PWD/10c_hardfiltering_output/genotyped_PASS_snps_raw_recal2.vcf | wc -l

# Get stats
bcftools stats $PWD/10c_hardfiltering_output/PMCO_genotyped_PASS_snps_recal2.vcf > stats_snps_recal2.csv
# A single report file is generated with summary statistics for each library processed containing the following pieces of information: Reads, Aligned Reads (% Aligned), Aligned Bases (Read Length, % Paired, % Duplicate, Mean Insert Size), SNPs, Filtered SNPs, SNPs after BQSR, Filtered SNPs after BQSR (Average Coverage)

##############!# VCF FILES READY FOR FORMATTING FOR DOWNSTREAM ANALYSIS #!#################

# Filter and reformat SNPs using vcftools for downstream analysis following parameters set by Erickson et al. (2021)
# --min-alleles 2 \ #to only include biallelic sites
# --max-alleles 2 \ #to only include biallelic sites
# --thin 1000 \ #to ensure that no two snps are within 1000 bp of one another (a proxy filter to get 1 snp per locus)
# --max-missing 0.75 \ #to set a threshold for taxon completeness (e.g. 0.75 only includes SNPs present for at least 75% of individuals, thus less than 25% of individuals have missing data)
# --max-non-ref-af 0.99 \ #includes only sites with all Non-Reference (ALT) Allele Frequencies (af) within the range specified, to filter out any SNPs that are the same across all samples
# --recode \ #to create a new vcf of filtered SNPs

module load conda3
source $CONDA_PROF/conda.sh
conda activate vcftools-0.1.16

vcftools --vcf 10c_hardfiltering_output/PMCO_genotyped_PASS_snps_recal2.vcf \
    --min-alleles 2 \
    --max-alleles 2 \
    --thin 1000 \
    --max-missing 0.75 \
    --recode \
    --out 11_formatting/PMCO_filtered_snps_vcf75
# --max-non-ref-af 0.99 implied when only keeping biallelic sites

## DAPC - Reformat with Plink to use in DAPC and structure analysis which can be use for species delimitation analysis
plink --vcf filtered_snps_vcf75.recode.vcf \
    --const-fid 0 \ #converts sample IDs to within-family IDs while setting all family IDs to a single value (default '0'). Basically it tells Plink to use the same name as those set in the vcf file
    --allow-extra-chr \ #allows new contig name format e.g uce-...
    --recode structure \
    --out structure_snps_vcf75

## iPYRAD vcf files for ipyrad-analysis will have to be pre-filtered to remove indels and low-quality SNPs (bcftools) and converted into HDF5 format (see ipyrad.pbs script)

## Format for SNAPP (BEAST2) species tree - needs input SNP file to be in a nexus format ##
# use perl script vcf2nex.pl provided by BEAST - I'll be running SNAPP on my laptop so ending the script here. Find SNAPP walkthrough in the SNAPP script.

## If you end up having to add one or more sample to your final callset, instead of running all the samples again and wait long enough to see teeth growing on chickens (another unnecessary french expression! you're very welcome), you can run the new samples down this pipeline separately MAPPED TO THE SAME REFERENCE SEQUENCE OBVISOULY, and then merge the final vcf files.

java -Xmx2g -jar /sw/picard/2.18.1/picard.jar MergeVcfs I= $PWD/10c_hardfiltering_output/PMCO_genotyped_PASS_snps_recal2.vcf I= $PWD/extra_sample_genotyped_PASS_snps_recal2.vcf O= updated_snp.vcf



