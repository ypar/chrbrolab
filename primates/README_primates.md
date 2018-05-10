

# README for the primates project directories
# YoSon Park

### last modified 05/10/2018


The primates project space is used by mtriz and me. The paper has been [published](https://www.ncbi.nlm.nih.gov/pubmed/28855262). Some of the relevant scripts have already been published as a separate [git repo](https://github.com/ypar/cre_evo_primates). Here I only list most of the workspace primarily used by me.



## Index 

<!--ts-->


* [ChIP-seq data and analysis results](#chipseq)
  * [aligned bam files](#aln)
  * [allele-specific expression analysis results](#ase)
  * [ensembl epo mapping data and analysis results](#ensembl_epo)
  * [fastq sequencing files](#fastq)
  * [_de novo_ and known motif search and enrichment analyses results](#motif)
  * [peak identification and relevant analyses results](#peaks)
  * [human-specific te flanking region data](#te_flanking)

* [RNA-seq data and analysis results](#rnaseq)
  *
  *
*
*
*


<!--te-->




# chipseq

## aln

## ase

## ensembl_epo

each dir w prefix countmatrix contain count matrices combining all species  
each dir w prefix readcounts contain per species read counts  

file name examples  

E.g.  
homo_sapiens_cons_FDR1_peak_e39_coords_1kbps_no_composite_count_matrix.txt  
homo_sapiens_cons_FDR1_peak_e39_coords_1kbps_no_composite_count_matrix_summary.txt  
homo_sapiens_cons_FDR1_peak_e39_coords_1kbps_no_composite_count_matrix_noNA.txt   
homo_sapiens_cons_FDR1_peak_e39_coords_1kbps_no_composite_count_matrix_noNA_summary.txt   


homo_sapiens: <prefix for reference species>   
cons: consensus peak among all included samples   
FDR1: consensus peaks are filtered down to fdr 1% only   
Peak: chipseq peaks   
e39 coords: eutharian 39 mammals database in ensembl msa coordinates are used  
1kbps: 1kbps +/- from the called peak are included for broader assessments before pinning down candidates in detail  
No composite: all regions with composite references are excluded    
Count matrix: read counts from chipseq experiments  
summary: descriptive summary statistics such as min, max. mean, etc are calculated from the count matrix  
noNA: any rows with NAs are removed  

following subdirectories are present.

### bedfiles_500bps  
### countmatrix_1kbps  
### countmatrix_500bps	
### countmatrix_fdr5_q10_w_500bps  
### readcounts_1kbps  
### readcounts_500bps  
### readcounts_fdr5_q10_w_500bps  
### sub_rate



## fastq

## motif

## peaks

## te_flanking





# rnaseq





