

# README for the primates project directories
# YoSon Park

### last modified 05/16/2018


The primates project space is used by mtriz and me. The paper has been [published](https://www.ncbi.nlm.nih.gov/pubmed/28855262). Some of the relevant scripts have already been published as a separate [git repo](https://github.com/ypar/cre_evo_primates). Here I only list most of the workspace primarily used by me. for any unlisted directories, contact mtriz.



## Index 

<!--ts-->

* [miscellaneous documents](#documents)


* [ChIP-seq data and analysis results](#chipseq)
  * [aligned bam files](#aln)
  * [allele-specific expression analysis results](#ase)
  * [ensembl epo mapping data and analysis results](#ensembl_epo)
  * [fastq sequencing files](#fastq)
  * [_de novo_ and known motif search and enrichment analyses results](#motif)
  * [peak identification and relevant analyses results](#peaks)
  * [human-specific te flanking region data](#te_flanking)

* [RNA-seq data and analysis results](#rnaseq)
  * [aligned bam files](#aln)
  * [fastq sequencing files](#fastq)
  * [kallisto quatification files](#kallisto_outputs)
  * [kallisto index files](#kallisto_index)



<!--te-->



# documents

there are some readmes generated for the working directories in the primates project space.



# chipseq

## aln

aligned chipseq data are in this directory.

## ensembl_epo

I put resulting files from homologous region analysis. most file names contain minimal information regarding the content of the file. 


file name examples:  

E.g.  
homo_sapiens_cons_FDR1_peak_e39_coords_1kbps_no_composite_count_matrix.txt  
homo_sapiens_cons_FDR1_peak_e39_coords_1kbps_no_composite_count_matrix_summary.txt  
homo_sapiens_cons_FDR1_peak_e39_coords_1kbps_no_composite_count_matrix_noNA.txt   
homo_sapiens_cons_FDR1_peak_e39_coords_1kbps_no_composite_count_matrix_noNA_summary.txt   


homo_sapiens: <prefix for reference species>   
cons: consensus peak among all included samples   
FDR1: consensus peaks are filtered down to fdr 10% only   
peak: chipseq peaks   
e39 coords: eutharian 39 mammals database in ensembl msa coordinates are used  
1kbps: 1kbps +/- from the called peak are included for broader assessments before pinning down candidates in detail  
no composite: all regions with composite references are excluded    
count matrix: read counts from chipseq experiments  
summary: descriptive summary statistics such as min, max. mean, etc are calculated from the count matrix  
noNA: any rows with NAs are removed  



following subdirectories are currently available:

### bedfiles_500bps  

bed files containing peak boundaries +/- 500 bps used for ChIP-seq analysis. following information is included:  


> head -2 bedfiles_500bps/callithrix_jacchus_cons_H3K27Ac_FDR1_CSpeaks_peak_e39_coords_500bps_no_composite_homo_sapiens.bed  

>> 17	80990879	80992910	MS_cons_H3K27Ac_FDR1_CSpeaks_peak_10000	0	1     
>> 17	81340999	81342614	MS_cons_H3K27Ac_FDR1_CSpeaks_peak_10007	0	1  

column 1: chromosome id  
column 2: starting position of the peak - 500 bps  
column 3: ending position of the peak + 500 bps  
column 4: name of the peak in the format of two-letter-species-name + conserved + histone modification type + multiple testing correction threshold + consensus across six species + peak + index\#  
column 5: score used for ucsc genome browser. not used for [bedtools](http://bedtools.readthedocs.io/en/latest/content/general-usage.html) and other tools used for our analysis. dummy variable 0 filled in.   
column 6: strand indicated as +1 and -1.  

For details regarding bed files, check out the [ucsc genome browser faq](http://bedtools.readthedocs.io/en/latest/content/general-usage.html) or [bedtools description of bed files](http://bedtools.readthedocs.io/en/latest/content/general-usage.html).


### countmatrix_1kbps  

I use the same naming convention I described above. an example is here:  

`otolemur_garnettii_cons_FDR1_peak_e39_coords_1kbps_no_composite_count_matrix_noNA.txt`
* this file includes count matrix for Otolemur garnettii consensus peaks at FDR<0.1, generated using ensembl epo 39 mammalians multi-species alignment mapping. reads within +/- 1kbps of called ChIP-seq peaks are considered. only non-composite regions are presented and any rows with N/As are removed. 



for these files, following information is included:  

> head -2 countmatrix_1kbps/otolemur_garnettii_cons_FDR1_peak_e39_coords_1kbps_no_composite_count_matrix_noNA.txt  

>> peakname	HS_B16_H3K27Ac	HS_B16_H3K4me1	HS_B16_INPUT	HS_B2_H3K27Ac	HS_B2_H3K4me1	HS_B2_INPUT	HS_B48_H3K27Ac	HS_B48_H3K4me1	HS_B48_INPUT	CH_40191_H3K27Ac	CH_40191_H3K4me1	CH_40191_INPUT	CH_4X0523_H3K27Ac	CH_4X0523_H3K4me1	CH_4X0523_INPUT	RH_27415_H3K27Ac	RH_27415_H3K4me1	RH_27415_INPUT	RH_27425_H3K27Ac	RH_27425_H3K4me1	RH_27425_INPUT	ML_7022f_Basil_H3K27Ac	ML_7022f_Basil_H3K4me1	ML_7022f_Basil_INPUT	ML_ALFA_ALFA_H3K27Ac	ML_ALFA_ALFA_H3K4me1	ML_ALFA_ALFA_INPUT	MS_17086_H3K27Ac	MS_17086_H3K4me1	MS_17086_INPUT	MS_31480_H3K27Ac	MS_31480_H3K4me1	MS_31480_INPUT	MS_32842_H3K27Ac	MS_32842_H3K4me1	MS_32842_INPUT	BB_1413m_Tuff_H3K27Ac	BB_1413m_Tuff_H3K4me1	BB_1413m_Tuff_INPUT	BB_8012_H3K27Ac	BB_8012_H3K4me1	BB_8012_INPUT  
>> BB_cons_H3K27Ac_FDR1_CSpeaks_peak_10011	694	31	165	235	121	56	39	43	51	137.0	40.0	26.0	50.0	18.0	31.0	3.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	163.0	56.0	10.0	24.0	59.0	43.0	295.0	58.0	67.0	463	105	45	168	53	28  

wherein each column represents the read counts for the peak id (rownames) in each sample. prefixes are species and sample IDs of each sample. the suffix indicates whether they are derived from H3K27acetylation and H2K4monomethylation histone modifications or baseline input ChIP-seq data.



### countmatrix_500bps	
see above for counmatrix file descriptions

### countmatrix_fdr5_q10_w_500bps  
see above for counmatrix file descriptions

### readcounts_1kbps  


I use the same naming convention as described above for files. an example is provided below:

`microcebus_murinus_cons_FDR1_peak_e39_coords_1kbps_no_composite_homo_sapiens_readcounts_w_header.txt`
* this file includes read counts for peaks called in Microcebus murinus, filtered to only included consensus regions among included samples at FDR<0.1, mapped to Homo sapiens ChiP-seq by using ensembl-epo-39-mammalians multi-species alignment with +/- 1kbps boundaries. only non-composite sites are presented with a header for column names.


for these files, following information is included with varying number of columns based on samples presented:

> head -2 readcounts_1kbps/microcebus_murinus_cons_FDR1_peak_e39_coords_1kbps_no_composite_homo_sapiens_readcounts_w_header.txt  
>> chrom	start	end	peakname	score	strand	HS_B16_H3K27Ac	HS_B16_H3K4me1	HS_B16_INPUT	HS_B2_H3K27Ac	HS_B2_H3K4me1	HS_B2_INPUT	HS_B48_H3K27Ac	HS_B48_H3K4me1	HS_B48_INPUT  
>> KI270736.1	11566	12411	ML_cons_H3K4me1_FDR1_CSpeaks_peak_5231	0	1	1	2	20	0	2	3	17	0	1  

chrom: chromosome id  
start: starting position of the peak - boundary (in this case 1kbps)  
end: ending position of the peak + boundary (in this case 1kbps)  
peakname: name of the peak for identification purposes  
score: score. not used and are filled with 0.  
strand: strand indicated as +1 or -1.  
HS_B16_H3K27Ac: readcounts for homo sapiens sample ID B16, ChIP-seq measuring histone modification H3K27acetylation  
HS_B16_H3K4me1: readcounts for homo sapiens sample ID B16, ChIP-seq measuring histone modification H3K4monomethylation  
HS_B16_INPUT: readcounts for homo sapiens sample ID B16, ChIP-seq measuring baseline without modification  
HS_B2_H3K27Ac: readcounts for homo sapiens sample ID B2, ChIP-seq measuring histone modification H3K27acetylation  
HS_B2_H3K4me1: readcounts for homo sapiens sample ID B2, ChIP-seq measuring histone modification H3K4monomethylation  
HS_B2_INPUT: readcounts for homo sapiens sample ID B2, ChIP-seq measuring baseline without modification  
HS_B48_H3K27Ac: readcounts for homo sapiens sample ID B48, ChIP-seq measuring histone modification H3K27acetylation  
HS_B48_H3K4me1: readcounts for homo sapiens sample ID B48, ChIP-seq measuring histone modification H3K4monomethylation  
HS_B48_INPUT: readcounts for homo sapiens sample ID B48, ChIP-seq measuring baseline without modification  



### readcounts_500bps  
see above for readcounts matrix descriptions

### readcounts_fdr5_q10_w_500bps  
see above for readcounts matrix descriptions


### sub_rate

I placed substitution rate analysis results here. 


## fastq

fastq files for the primates project ChIP-seq data. 

reminder: All raw sequence data from this study have been submitted to the [NCBI BioProject database](https://www.ncbi.nlm.nih.gov/bioproject/) under accession numbers [PRJNA349047 (RNA-seq)](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA349047) and [PRJNA349046 (ChIP-seq)](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA349046).


## motif

I ran motif search and enrichment analyses using the [MEME-ChIP suite](http://meme-suite.org/tools/meme-chip). here's a paper describing it:
[Timothy L. Bailey, Mikael Bodén, Fabian A. Buske, Martin Frith, Charles E. Grant, Luca Clementi, Jingyuan Ren, Wilfred W. Li, William S. Noble, "MEME SUITE: tools for motif discovery and searching", Nucleic Acids Research, 37:W202-W208, 2009.](http://nar.oxfordjournals.org/content/37/suppl_2/W202.full)



## peaks


### peaks/consensus_peaks

consensus peaks refer to peaks that are identified across all studied samples of the species. 


### peaks/plots

I placed some plots used to assess ChIP-seq analysis results.



## te_flanking






# rnaseq


## aln

the latest, star-aligned and q10 filtered bam files are here.


## fastq

fastq files for the primates project RNA-seq data. 

reminder: All raw sequence data from this study have been submitted to the [NCBI BioProject database](https://www.ncbi.nlm.nih.gov/bioproject/) under accession numbers [PRJNA349047 (RNA-seq)](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA349047) and [PRJNA349046 (ChIP-seq)](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA349046).


## kallisto_outputs

kallisto quantifications. for details regarding kallisto or any updated software releases, check out the official [website](https://pachterlab.github.io/kallisto/about).


## kallisto_index

kallisto index files used for quantifications.











