
# README for the ipsc project directories
# YoSon Park



For sample ID conversion related documents, see [this README](https://github.com/ypar/chrbrolab/blob/master/ipsc/README_master_fastq_list_for_dbgap.md).

additional information may be updated via public dbGaP profile under the project's full name [The National Heart, Lung, and Blood Institute (NHLBI)-funded Next Generation Genetic Association Studies (NextGen) Consortium: Phenotyping Lipid traits in iPS derived hepatocytes Study (PhLiPS Study)](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001341.v1.p1).


## Index

<!--ts-->

* [citations for data usage](#citations-for-data-usage)

* [chrbrolab directories](#chrbrolab-directories)
  * [subdir for cis-eqtl analysis results](#ciseqtl-results) 
  * [subdir for replication data](#replication-data)
  * [subdir for splicing analysis data](#splicing-analysis-data)
  * [subdir for meta-analysis results](#metaanalysis-results)

* [chrb_rader_lab data subdirectories](#chrb_rader_lab-directories)

* [raderlab_subdirectories](#raderlab-directories)


<!--te-->


# citations for data usage

if you are using the ipsc project data, please cite this paper: [Large, Diverse Population Cohorts of hiPSCs and Derived Hepatocyte-like Cells Reveal Functional Genetic Variation at Blood Lipid-Associated Loci.](https://www.ncbi.nlm.nih.gov/pubmed/28388432)


# chrbrolab directories

## ciseqtl results


## replication data

I used embrionic stem cell and induced pluripotent stem cell gene expression matrices from choi et al (dataset accession number; GSE73211) from the following paper: [A comparison of genetically matched cell lines reveals the equivalence of human iPSCs and ESCs](https://www.nature.com/articles/nbt.3388).



## splicing analysis data

I used leafcutter to quantify reads a splicing junctions. following files are available:

'${tissue}_${sample}_bam2junc.bed.gz'

head -2 view here:
> zcat bam2junc/primaryheps/primaryheps_HU8110_bam2junc.bed.gz| head -2
>> chr1	14716	14981	HWI-D00712:1297:C7MJLANXX:3:1309:10584:68289	0	-	14716	14981	0	2	113,12	0,253
>> chr1	14720	14985	HWI-D00712:1315:C89UCANXX:3:1112:10160:40095	0	-	14720	14985	0	2	109,16	0,249

'${tissue}_${sample}_bam2junc.junc'

> head -2 bam2junc/primaryheps/primaryheps_HU8110_bam2junc.junc
>> chr17	48222595	48236670	.	15	-
>> chr3	195245982	195250486	.	8	-


here is an an example script for leafcutter bam2junc runs

> ${tissue}_${sample}_bam2junc_bsub.sh

```
#!/bin/bash

leafcutterdir=bin/src/leafcutter
samtoolsdir=bin/local/

wkdir=ipsc/splicing/bam2junc/${tissue}/
mkdir -p \${wkdir}

bedfile=${tissue}_${sample}_bam2junc.bed
juncfile=${tissue}_${sample}_bam2junc.junc
echo "converting $bamfile to \$juncfile"

$samtoolsdir/samtools view $bamfile | python2.7 $leafcutterdir/scripts/filter_cs.py | $leafcutterdir/scripts/sam2bed.pl --use-RNA-strand - \${wkdir}/\$bedfile

$leafcutterdir/scripts/bed2junc.pl \${wkdir}/\$bedfile \${wkdir}/\$juncfile

```

## metaanalysis results

I used metatissue to generated gene expression and genotype matrices adjusted for correlated data structure implementing mixed model estimators.


# chrb_rader_lab directories


# raderlab directories


