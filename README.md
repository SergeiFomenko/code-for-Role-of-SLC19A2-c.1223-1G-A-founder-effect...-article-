# Haplotype Analysis for SLC19A2 Founder Effect Study

This repository contains code used in the publication:

**"Role of SLC19A2 c.1223+1G\>A founder effect in forming a new endemic
region of patients with thiamine-responsive megaloblastic anemia"**

------------------------------------------------------------------------

## Overview

This script performs haplotype reconstruction and mismatch rate analysis
for a candidate founder variant using whole-genome sequencing data. 
Main step of analysis:

- Identifies homozygous and heterozygous carriers of a target
variant
- Constructs a consensus haplotype from homozygous patients
- Computes genotype concordance and mismatch rates
- Generates summary statistics and plot for region visualisation

------------------------------------------------------------------------

## Input data

Required inputs:

-   --multi_vcf : Multisample VCF (joint-called, e.g. GLnexus)\
-   --variant : Target variant (CHROM_POS_REF_ALT)\
-   --region : Region of interest (CHROM_START_END)\
-   --nokin_subset : List of unrelated individuals (one sample ID per
    line)\
-   --genes : Gene annotation (BED format)\
-   --results_path : Output Excel file\
-   --plot_path : Output figure

------------------------------------------------------------------------

## Methodology

### Haplotype construction

The consensus haplotype is defined using variants that are: - uniformly
homozygous (0/0 or 1/1) across all homozygous patients\
- heterozygous sites are excluded

### Concordance in patients

-   patient_heterorate: fraction of patients with 0/1 genotype\
-   patients_match_percent: maximal fraction sharing same homozygous
    genotype

### Mismatch rate

Mismatch is defined as discordant homozygous genotype relative to
patient haplotype: - 0/0 vs 1/1\
- 1/1 vs 0/0

------------------------------------------------------------------------

## Output

### Excel file
Sheet: snp_concordance_patients\
- patient_heterorate - fraction of homozygous variant carriers with 0/1 genotype in this position
- homoref_rate_patients - fraction of homozygous variant carriers with 0/0 genotype in this position
- homoalt_rate_patients - fraction of homozygous variant carriers with 1/1 genotype in this position
- patients_match_percent - patients genotype concordance, maximal fraction of patient that have either 0/0 or 1/1 genotype in this position 
Sheet: mismatch_rate_population\
-mismatch_carriers - fraction of mismatches in heterozygous variant carriers compared to concensus homozygous carriers haplotype (0/0 vs 1/1 or 1/1 vs 0/0) in this position
-mismatch_non-carriers - fraction of mismatches in non-carriers compared to concensus homozygous carriers haplotype (0/0 vs 1/1 or 1/1 vs 0/0) in this position

------------------------------------------------------------------------

## Usage example

python analyse_haplotype.py\
--multi_vcf cohort.vcf.gz\
--variant chr1_169549811_G_A\
--region chr1_168981142_171250829\
--nokin_subset unrelated.txt\
--genes annotation/gene.bed\
--results_path results.xlsx\
--plot_path plot.png

------------------------------------------------------------------------

## Installation

conda env create -f environment.yml

------------------------------------------------------------------------

## Notes

-   VCF must contain biallelic SNPs\
-   Genotype encoding follows cyvcf2 convention\
-   Sample IDs must match between inputs

