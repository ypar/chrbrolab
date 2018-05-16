
# README for the reference directory
# YoSon Park

### last modified 04/26/2018

On hpc pmacs environment, I have set up directories and downloaded reference files for several useful resources of general use. For future projects, I also list websites where one may find additional or updated information on relevant files. While I stopped updating these directories, hopefully other users of the chrbrolab project space find the references and citations useful.



## Index

<!--ts-->
* [1kgp_phase3](#1kgp_phase3)
  * 1kgp_phase3 vcf files
  * [1kgp_phase3/haplotype_reference_for_impute2](#1kgp_phase3haplotype_reference_for_impute2)
* [dbsnp](#dbsnp)
* [ensembl](#ensembl)
  * [ensembl/chromosome_sizes](#ensemblchromosome_sizes)
  * [ensembl/epo_39_eutherian](#ensemblepo_39_eutherian)
  * [ensembl/gtf_files](#ensemblgtf_files)
  * [ensembl/transcriptomes](#ensembltranscriptomes)
* [exac](#exac)
  * exac/release0.1
  * exac/release0.2
  * exac/release0.3
  * exac/release1
* [fusion](#fusion)
* [gatk](#gatk)
* [gencode](#gencode)
* [geuvadis](#geuvadis)
* [hapmap](#hapmap)
* [homologene](#homologene)
* [liftover](#liftover)
* [omni](#omni)
* [roadmap](#roadmap)
* [rsem](#rsem)
* [star](#star)
<!--te-->







# 1kgp_phase3

This directory contains vcf files of the 1000 Genomes Project (Phase 3) named as `ALL.chr${chromosome_number}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz`. I downloaded these files from the official ftp site here (`ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/`). Relevant documents are in the directory and additional information or other files can be fetched from the official ftp site.

## 1kgp_phase3/haplotype_reference_for_impute2/

This directory contains haplotype, legend files, etc. used for [shapeit2](http://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html), [impute2](http://mathgen.stats.ox.ac.uk/impute/impute_v2.html) or related software.



# dbsnp

I downloaded a list of SNP coordinates for the hg19 dbSNP version 142 here.



# ensembl


## ensembl/chromosome_sizes

Files contain information of chromosome sizes for seven species studied for primate evolution of regulatory elements. I downloaded these reference files indicating chromosome sizes for our paper published [here](https://genome.cshlp.org/content/early/2017/09/08/gr.218149.116). I renamed files to represent an abbreviated version of the species' full names, which may not be identical to their genome assembly names officially used by the ensembl database.

* ensembl/chromosome_sizes/OLD_UCSC  
For the paper, we used ensembl references. All previously downaloded UCSC reference-based chromosome sizes are relocated into this directory by MT.

* ensembl/chromosome_sizes/calJac  
Chromosome size information for _Callithrix jacchus_ (Marmoset) ensembl reference version calJac3. Latest version should be checked on the official [database](https://useast.ensembl.org/Callithrix_jacchus/Info/Index) before analysis. 

* ensembl/chromosome_sizes/macMul  
Chromosome size information for _Macaca mulatta_ (Rhesus macaque) ensembl reference version rheMac2, also previously prefixed Mmul and macMul in different contexts. Latest version should be checked on the official [database](https://useast.ensembl.org/Macaca_mulatta/Info/Index) before analysis. 

* ensembl/chromosome_sizes/otoGar  
Chromosome size information for _Otolemur garnettii_ (Bushbaby) ensembl reference version otoGar3. Latest version should be checked on the official [database](https://useast.ensembl.org/Otolemur_garnettii/Info/Index) before analysis.  

* ensembl/chromosome_sizes/tupBel
Chromosome size information for _Tupaia belangeri_ (Tree Shrew) ensembl reference version tupBel1. Latest version should be checked on the official [database](https://useast.ensembl.org/Tupaia_belangeri/Info/Index) before analysis. 

* ensembl/chromosome_sizes/homSap  
Chromosome size information for _Homo sapiens_ (Human) ensembl reference version GRCh38. Latest version should be checked on the official [database](https://useast.ensembl.org/Homo_sapiens/Info/Index) before analysis. 

* ensembl/chromosome_sizes/micMur  
Chromosome size information for _Microcebus murinus_ (Mouse lemur) ensembl reference version micMur3. The ensembl reference is named Mmur. Latest version should be checked on the official [database](https://useast.ensembl.org/Macaca_mulatta/Info/Index) before analysis. 

* ensembl/chromosome_sizes/panTro  
Chromosome size information for _Pan troglodytes_ (Chimpanzee) ensembl reference version panTro3. On ensembl, version names may be prefixed Pan_tro or CHIMP. Latest version should be checked on the official [database](https://useast.ensembl.org/Pan_troglodytes/Info/Index) before analysis. 

The files look like this:

```
> head -2 ensembl/chromosome_sizes/otoGar/otoGar3.chrom.sizes
GL873520	73491278	/gbdb/otoGar3/otoGar3.2bit
GL873521	55144334	/gbdb/otoGar3/otoGar3.2bit
```


## ensembl/epo_39_eutherian

Ensembl comparative genomics resources were described [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4761110/). There now exists an extensive documentation for using ensembl-compara [here](https://media.readthedocs.org/pdf/ensembl-compara/master/ensembl-compara.pdf), if you are interested in using more of this database. I have used EPO_39_Eutherian dataset and some supplemental files are downdloaded here. All queries for [my paper](https://genome.cshlp.org/content/early/2017/09/08/gr.218149.116) were done using the database API directly and did not use any offline reference files.


## ensembl/gtf_files

GTF files for the ensembl Homo sapiens versions GRCh37 (ensembl 65), GRCh37 (ensembl 75) and GRCh38 (ensembl 79).



## ensembl/transcriptomes

These files are primarily downloaded in 2015. I used the following commands to download reference transcriptome files from ensembl for [our primates evolution of regulatory elements paper](https://genome.cshlp.org/content/early/2017/09/08/gr.218149.116). 

```
wget -nd --mirror ftp://ftp.ensembl.org/pub/release-81/fasta/pan_troglodytes/cdna/
#for hg19, use release 75 or gencode19
wget -nd --mirror ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cdna/
#for hg38, use release 81, or gencode 22/23
wget -nd --mirror ftp://ftp.ensembl.org/pub/release-81/fasta/homo_sapiens/cdna/
wget -nd --mirror ftp://ftp.ensembl.org/pub/release-81/fasta/microcebus_murinus/cdna/
wget -nd --mirror ftp://ftp.ensembl.org/pub/release-81/fasta/tupaia_belangeri/cdna/
wget -nd --mirror ftp://ftp.ensembl.org/pub/release-81/fasta/macaca_mulatta/cdna/
wget -nd --mirror ftp://ftp.ensembl.org/pub/release-81/fasta/callithrix_jacchus/cdna/
```


* transcriptomes/calJac  
Transcriptome reference for _Callithrix jacchus_ (Marmoset) ensembl reference version calJac3. Latest version should be checked on the official [database](https://useast.ensembl.org/Callithrix_jacchus/Info/Index) before analysis. 

* transcriptomes/macMul  
Transcriptome reference for _Macaca mulatta_ (Rhesus macaque) ensembl reference version rheMac2, also previously prefixed Mmul and macMul in different contexts. Latest version should be checked on the official [database](https://useast.ensembl.org/Macaca_mulatta/Info/Index) before analysis. 

* transcriptomes/otoGar  
Transcriptome reference for _Otolemur garnettii_ (Bushbaby) ensembl reference version otoGar3. Latest version should be checked on the official [database](https://useast.ensembl.org/Otolemur_garnettii/Info/Index) before analysis.  

* transcriptomes/tupBel
Transcriptome reference for _Tupaia belangeri_ (Tree Shrew) ensembl reference version tupBel1. Latest version should be checked on the official [database](https://useast.ensembl.org/Tupaia_belangeri/Info/Index) before analysis. 

* transcriptomes/homSap  
Transcriptome reference for _Homo sapiens_ (Human) ensembl reference version GRCh38. Latest version should be checked on the official [database](https://useast.ensembl.org/Homo_sapiens/Info/Index) before analysis. 

* transcriptomes/micMur  
Transcriptome reference for _Microcebus murinus_ (Mouse lemur) ensembl reference version micMur3. The ensembl reference is named Mmur. Latest version should be checked on the official [database](https://useast.ensembl.org/Macaca_mulatta/Info/Index) before analysis. 

* transcriptomes/panTro  
Transcriptome reference for _Pan troglodytes_ (Chimpanzee) ensembl reference version panTro3. On ensembl, version names may be prefixed Pan_tro or CHIMP. Latest version should be checked on the official [database](https://useast.ensembl.org/Pan_troglodytes/Info/Index) before analysis. 



# exac

I downloaded reference files for analysis from the ExAC [(Exome Aggregation Consortium)](http://exac.broadinstitute.org/) ftp site here `ftp://ftp.broadinstitute.org/pub/ExAC_release`. Detailed descriptions can be found in relevant directory READMEs or on the official website [here](http://exac.broadinstitute.org/downloads). Subdirectory names refer to the data freeze / release version numbers.



# fusion

I downloaded files from the [fusion website](http://gusevlab.org/projects/fusion/) for integrative analysis of transcriptome-wide association data. The reference paper is [here](https://www.ncbi.nlm.nih.gov/pubmed/26854917). This is an excerpt from the official website:

> Functional Summary-based Imputation
>> FUSION is a suite of tools for performing a transcriptome-wide (or any other ome-wide) association study by predicting functional/molecular phenotypes into GWAS using only summary statistics. The goal is to identify associations between a GWAS phenotype and a functional phenotype that was only measured in reference data. We provide precomputed functional weights (primarily gene expression) from multiple studies to facilitate this analysis.



# gatk

HG19/GRCh37 genome reference files I generated for analyses using [GATK](https://software.broadinstitute.org/gatk/).



# gencode

The gencode annotation file downloaded from the [official website](https://www.gencodegenes.org/). Currently, gencode v19 gene annotations for GRCh37 and gencode v28 gene annotation for GRCh38 are in respectively named subdirectories.



# geuvadis

I downloaded some RNA-seq expression files and metadata from the [geuvadis project](https://www.ebi.ac.uk/Tools/geuvadis-das/). Geuvadis RNA-sequencing project uses mRNA and small RNA sequencing data from the same lymphoblastoid cell lines of 465 individuals included in the 1000 genomes project. Files were downloaded from ftp site (example fastq file: `ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR188021/ERR188021_1.fastq.gz`).



# hapmap

The HapMap phase3 reference files are here. Plink and hapmap formats are available. Specifically, plink_format contains data from the draft release 2 for genome-wide SNP genotyping in DNA samples from 11 HapMap 3 populations.



# homologene

I downloaded HomoloGene from here `ftp://ftp.ncbi.nih.gov/pub/HomoloGene/build68/`. Detailed file descriptions are here `ftp://ftp.ncbi.nih.gov/pub/HomoloGene/HomoloGene_Field_Description.txt`. Here is the explanation of the HomoloGene Build Procedure from the [website](https://www.ncbi.nlm.nih.gov/homologene/build-procedure/):

> HomoloGene Build Procedure
>> The input for HomoloGene processing consists of the proteins from the input organisms. These sequences are compared to one another (using blastp) and then are matched up and put into groups, using a tree built from sequence similarity to guide the process, where closer related organisms are matched up first, and then further organisms are added as the tree is traversed toward the root. The protein alignments are mapped back to their corresponding DNA sequences, where distance metrics can be calculated (e.g. molecular distance, Ka/Ks ratio). Sequences are matched using synteny when applicable. Remaining sequences are matched up by using an algorithm for maximizing the score globally, rather than locally, in a bipartite matching. Cutoffs on bits per position and Ks values are set to prevent unlikely "orthologs" from being grouped together. These cutoffs are calculated based on the respective score distribution for the given groups of organisms. Paralogs are identified by finding sequences that are closer within species than other species.


# liftover

I downloaded the UCSC liftover references for [this paper](https://genome.cshlp.org/content/early/2017/09/08/gr.218149.116). 




# omni

Reference files to process Illumina HumanOmni BeadChips.



# roadmap

I downloaded liver-related data from the [epigenomics roadmap project](http://www.roadmapepigenomics.org/). 



# rsem

HG19 transcriptome reference files I generated for analyses using [RSEM](https://deweylab.github.io/RSEM/).




# star

HG19 transcriptome reference files I generated for RNA-seq data processing using [STAR](https://github.com/alexdobin/STAR).




