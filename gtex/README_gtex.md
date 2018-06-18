
# README for the gtex directory
# YoSon Park



On hpc pmacs, I have made a directory where others in the chrbrolab project space may access gtex data for their own project(s). Most data available on pmacs are officially released data unless otherwise explained here. I tried to organize the subdirectories such that I don't have to add a lot of descriptions here.



## Index

<!--ts-->
* [citations for data usage](#citations-for-data-usage)

* [data subdirectories](#data)
  * [subdir for GTEx v4 data](#gtex-v4-data)
    * [data/GTEx_v4_auxiliary](#datagtex_v4_auxiliary)
    * [data/GTEx_v4_metadata](#datagtex_v4_metadata)
    * [data/GTEx_Analysis_v4_eQTL_raw_expression_matrices](#datagtex_analysis_v4_eqtl_raw_expression_matrices)
  * [subdir for GTEx v6 data](#gtex-v6-data)
    * [data/GTEx_v6p_allelecount](#datagtex_v6p_allelecount)
    * [data/GTEx_v6p_auxiliary](#datagtex_v6p_auxiliary)
    * [data/GTEx_v6p_genotype](#datagtex_v6p_genotype)
    * [data/GTEx_v6p_genotype_dosage_adjusted_by_covariates](#datagtex_v6p_genotype_dosage_adjusted_by_covariates)
    * [data/GTEx_v6p_genotype_original](#datagtex_v6p_genotype_original)
    * [data/GTEx_v6p_genotype_phased](#datagtex_v6p_genotype_phased)
    * [data/GTEx_v6p_metadata](#datagtex_v6p_metadata)
    * [data/GTEx_Analysis_v6p_eQTL_expression_matrices_internal](#datagtex_analysis_v6p_eqtl_expression_matrices_internal)
    * [data/GTEx_Analysis_v6p_eQTL_expression_matrices_tpm](#datagtex_analysis_v6p_eqtl_expression_matrices_tpm)
  * [subdir for GTEx v7 data](#gtex-v7-data)
    * [data/GTEx_v7_genotype](#datagtex_v7_genotype)
    * [data/GTEx_v7_metadata](#datagtex_v7_metadata)
    * [data/GTEx_Analysis_v7_eQTL_raw_expression_matrices](#datagtex_analysis_v7_eqtl_raw_expression_matrices)
    * [data/GTEx_Analysis_v7_eQTL_normalized_expression_matrices](#datagtex_analysis_v7_eqtl_normalized_expression_matrices)
  * [subdir for GTEx v8 data](#gtex-v8-data)
    * [data/GTEx_v8_genotype](#datagtex_v8_genotype)
    * [data/GTEx_v8_genotype_dosage_adjusted_by_covariates](#datagtex_v8_genotype_dosage_adjusted_by_covariates)
    * [data/GTEx_v8_metadata](#datagtex_v8_metadata)
    * [data/GTEx_Analysis_v8_eQTL_normalized_expression_matrices](#datagtex_analysis_v8_eqtl_normalized_expression_matrices)
    * [data/GTEx_Analysis_v8_eQTL_normalized_expression_matrices_adjusted_by_covariates](#datagtex_analysis_v8_eqtl_normalized_matrices_adjusted_by_covariates)


* [results subdirectories](#results)
  * [subdir for GTEx v4 results](#gtex-v4-results)
    * [results/GTEx_Analysis_v4_eQTL](#resultsgtex_analysis_v4_eqtl)
  * [subdir for GTEx v6 results](#gtex-v6-results)
    * [results/GTEx_Analysis_v6_eQTL](#resultsgtex_analysis_v6_eqtl)
    * [results/GTEx_Analysis_v6_metasoft](#resultsgtex_analysis_v6_metasoft)
    * [results/GTEx_Analysis_v6p_eQTL](#resultsgtex_analysis_v6p_eqtl)
    * [results/GTEx_Analysis_v6p_eQTL_ld_proxies](#resultsgtex_analysis_v6p_eqtl_ld_proxies)
    * [results/GTEx_Analysis_v6p_eQTL_ld_proxies_100kbps_sentinel](#resultsgtex_analysis_v6p_eqtl_ld_proxies_100kbps_sentinel)
    * [results/GTEx_Analysis_v6p_metasoft](#resultsgtex_analysis_v6p_metasoft)
    * [results/GTEx_Analysis_v6p_metatissue](#resultsgtex_analysis_v6p_metatissue)
  * [subdir for GTEx v7 results](#gtex-v7-results)
    * [results/GTEx_Analysis_v7_eQTL](#resultsgtex_analysis_v7_eqtl)
  * [subdir for GTEx v8 results](#gtex-v8-results)
    * [results/GTEx_Analysis_v8_eQTL](#resultsgtex_analysis_v8_eqtl)
    * [results/GTEx_Analysis_v8_metasoft](#resultsgtex_analysis_v8_metasoft)
  * [subdir for other files](#other-files)
    * [results/downloaded_files](#resultsdownloaded_files)
<!--te-->



# citations for data usage

There are four different GTEx datasets on hpc pmacs. There are several publications specific for each application. You can find many papers from the consortium members [here](https://gtexportal.org/home/publicationPage). Here, I list the main (_flagship_ so to speak) papers you can use for citations.

[GTEx v3 or v4 data](http://www.sciencemag.org/content/348/6235/648.short)

[GTEx v6 or v6p or v7 data](http://dx.doi.org/10.1038/nature24277)



# data

GTEx main consortium provides genotypes, whole exome and genome sequences and RNA-seq data. One could also find [histological image](https://gtexportal.org/home/histologyPage) data or apply for [biobank samples](https://gtexportal.org/home/samplesPage) on the portal. For candidate gene studies, analyses using summary statistics / data with no privacy concerns, and/or visualizations, etc., check out the [GTEx portal](https://gtexportal.org/home/) for the latest version. 

Also available locally are the final release of the GTEx data, v8, which I use for analyses and contribute to the consortium publication(s). **Preliminary analyses results I have either downloaded or generated on pmacs may change at the time of official release so always check for the latest versions and documentations before using them.**



## GTEx v4 data

The GTEx v3 pilot or v4 data publication is [here](http://science.sciencemag.org/content/348/6235/648.full). The v3 pilot data or subsequent data freeze v4 are seldom used as of data freeze v6 release. From dbGaP, one may check full documentation and raw data using its accession number, `phs000424.v4.p1`.


### data/GTEx_v4_auxiliary

For GTEx v4, the reference was a gene-level model based on transcript models of [gencode](https://www.gencodegenes.org/) v18. For GTEx v3, gencode v12 was used. References for v4 analyses are downloaded here. I also downloaded the list of variants include in the v4 analyses. The genotypes used for v4 were imputed based on the 1000 Genomes phase 1 reference.


### data/GTEx_v4_metadata

Sample information and other data. Also includes covariates used for eQTL mapping.


### data/GTEx_Analysis_v4_eQTL_raw_expression_matrices

Expression matrices for the following are available:

```
fraction of intron that is covered by reads
intron read count
junction read count
transcript read count
transcript rpkm
exon read count
gene read count
gene rpkm
```




## GTEx v6 data

For v6, the gene-level annotations were updated mid-release, and all subsequent analyses refer to this updated data as v6p for _v6 patched_. For all official uses, v6p replaces v6 or are interchangeably used for data descriptions in papers, etc.


### data/GTEx_v6p_allelecount

I generated allele counts using star 2-pass alignment, wasp qc processing steps and samtools. Example scripts are here: [run_wasp](https://github.com/ypar/chrbrolab/tree/master/gtex/data_processing/run_wasp.py) and [run_star2pass.py](https://github.com/ypar/chrbrolab/tree/master/gtex/data_processing/run_star2pass.py). Allele counts generated by LDACC using tophat v1.4 aligned bam files are available on dbGaP.



### data/GTEx_v6p_auxiliary



### data/GTEx_v6p_genotype

### data/GTEx_v6p_genotype_dosage_adjusted_by_covariates

### data/GTEx_v6p_genotype_original

### data/GTEx_v6p_genotype_phased

### data/GTEx_v6p_metadata

De-identified, open access versions of the sample and subject information and corresponding dictionaries.

### data/GTEx_Analysis_v6p_eQTL_expression_matrices_internal



### data/GTEx_Analysis_v6p_eQTL_expression_matrices_tpm

I generated tpm expression matrices using the gtex provided rpkm matrices. I include the script [here](https://github.com/ypar/chrbrolab/tree/master/gtex/data_processing/tpms_from_rpkms.py).




## GTEx v7 data

The latest publicly available data freeze now is v7, which contains

|V7 Release	          |# Tissues   |# Donors   | # Samples   |
|---------------------|------------|-----------|-------------|
|Total	              | 53	       |  714      | 11688       |
|With Genotype        | 53	       |  635      | 10361       |
|Has eQTL Analysis\*  | 48         |  620      | 10294       |

\* Number of samples with genotype >= 70

Sample count per tissue:
![sample_count](https://github.com/ypar/chrbrolab/blob/master/gtex/misc_plots/tdsNumSamplesByTissuesBarChart.svg)

^ above sample count plot is downloaded from [the official website](https://gtexportal.org/home/tissueSummaryPage)


### data/GTEx_v7_genotype

In this directory, I downloaded the official release of GTEx Whole Genome Sequences (WGS) for dbGaP accession number: phs000424.v7.p1 (release v7, GTEx_Analysis_2016-01-15). Detailed information is provided within the README of the directory. 


### data/GTEx_v7_metadata

### data/GTEx_Analysis_v7_eQTL_raw_expression_matrices

### data/GTEx_Analysis_v7_eQTL_normalized_expression_matrices





## GTEx v8 data

### data/GTEx_v8_genotype

### data/GTEx_v8_genotype_dosage_adjusted_by_covariates

### data/GTEx_v8_metadata

### data/GTEx_Analysis_v8_eQTL_normalized_expression_matrices

### data/GTEx_Analysis_v8_eQTL_normalized_expression_matrices_adjusted_by_covariates






# results


## GTEx v4 results

### results/GTEx_Analysis_v4_eQTL

(All) cis-eQTL association results from 14 tissues are available:
```
Adipose_Subcutaneous  
Cells_Transformed_fibroblasts  
Muscle_Skeletal		
Stomach
Artery_Aorta	  
Esophagus_Mucosa		     
Heart_Left_Ventricle			
Nerve_Tibial		
Thyroid
Artery_Tibial	  
Esophagus_Muscularis	     
Lung					
Skin_Sun_Exposed_Lower_leg	
Whole_Blood
```




## GTEx v6 results

### results/GTEx_Analysis_v6_eQTL

### results/GTEx_Analysis_v6_metasoft

### results/GTEx_Analysis_v6p_eQTL

(All) cis-eQTL association results from 44 tissues are available:

```
Adipose_Subcutaneous
Adipose_Visceral_Omentum
Adrenal_Gland
Artery_Aorta
Artery_Coronary
Artery_Tibial
Brain_Anterior_cingulate_cortex_BA24
Brain_Caudate_basal_ganglia
Brain_Cerebellar_Hemisphere
Brain_Cerebellum
Brain_Cortex
Brain_Frontal_Cortex_BA9
Brain_Hippocampus
Brain_Hypothalamus
Brain_Nucleus_accumbens_basal_ganglia
Brain_Putamen_basal_ganglia
Breast_Mammary_Tissue
Cells_EBV-transformed_lymphocytes
Cells_Transformed_fibroblasts
Colon_Sigmoid
Colon_Transverse
Esophagus_Gastroesophageal_Junction
Esophagus_Mucosa
Esophagus_Muscularis
Heart_Atrial_Appendage
Heart_Left_Ventricle
Liver
Lung
Muscle_Skeletal
Nerve_Tibial
Ovary
Pancreas
Pituitary
Prostate
Skin_Not_Sun_Exposed_Suprapubic
Skin_Sun_Exposed_Lower_leg
Small_Intestine_Terminal_Ileum
Spleen
Stomach
Testis
Thyroid
Uterus
Vagina
Whole_Blood
```


### results/GTEx_Analysis_v6p_eQTL_ld_proxies

I used plink v 1.09 to estimate ld and find ld proxy variants from each of the genome-wide significant eQTLs (q<0.05). an example script is [here](https://github.com/ypar/chrbrolab/blob/master/gtex/finemapping/ld_gtex450_allsigeqtl.py).



### results/GTEx_Analysis_v6p_eQTL_ld_proxies_100kbps_sentinel

For this directory, I generated ld proxies for sentinel SNPs only. I have also generated annotated tables for downstream analyses.  

the following columns are included per file  

within header, the suffix \_variant indicates the variant of interest and \_sentinel indicates the variant that has been tagged as the sentinel variant
in one or more eGenes in the given tissue for the eQTL scan by GTEx  

* rs_id_dbSNP142_GRCh37p13_variant: rsID of the variant in the genome-wide significant eqtl file
* SNP_variant: gtex common variant ID
* is_genomewide_eqtl_variant: 1 if variant is included in the genome-wide eqtl file
* MAF_variant: within-data maf of the bvariant
* SNP_sentinel: sentinel SNP
* MAF_sentinel: within data maf of the sentinel SNP
* R2: r2 estimated from 450 inds in gtex genotype between the variant and the sentinel
* gene_variant: gene assigned to the variant in the genome-wide eqtl results file
* gene_name_variant: common gene name for gene_variant
* orientation_variant: strand orientation of the gene_variant
* tss_distance_variant: distance to tss from the variant
* gene_sentinel: gene assigned to the sentinel in the genome-wide eqtl results file
* gene_name_sentinel: common gene name for the sentinel
* rs_id_dbSNP142_GRCh37p13_sentinel: rsID of the sentinel
* orientation_sentinel: strand orientation of the gene_sentinel
* tss_distance_sentinel: distance to tss from the sentinel





### results/GTEx_Analysis_v6p_metasoft

Multi-tissue association results from analyses using [metasoft](http://genetics.cs.ucla.edu/meta/).

### results/GTEx_Analysis_v6p_metatissue

I generated tissue specificity analysis results using metatissue and metasoft. 



## GTEx v7 results

### results/GTEx_Analysis_v7_eQTL





## GTEx v8 results

### results/GTEx_Analysis_v8_eQTL

### results/GTEx_Analysis_v8_metasoft

Multi-tissue association results from analyses using [metasoft](http://genetics.cs.ucla.edu/meta/).





## other files

### results/downloaded_files

Here I have backup copies of the v6(p) eQTL and metasoft results. In case of corrupted / overwritten summary stats, etc., extract these tar archives to respective directories.






