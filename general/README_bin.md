
# README for the lab-wide software and resource directories
# YoSon Park

### last modified 04/27/2018


For general use by chrbrolab project space users, I have, in the past, downloaded, installed and/or organized several commonly used tools. I no longer maintain these directories but I highly encourage others to keep following the directory organizing convention. I do not install or maintain user-specific resources for other users, e.g. R packages on user-specific environments. Also hpc pmacs common environment now maintains much of the hard-to-install software or modules. I list these subdirectories here so one may look up relevant websites and references for their own applications.

**Always check the latest on pmacs before considering (re)installation in the home environment. Always check the official websites or github repositories for updates and patches for each tool before use.**



## Index

<!--ts-->

* [bin/local](#binlocal)
* [bin/src](#binsrc)
  * [bin/src/ensembl](#binsrcensembl)
  * [bin/src/ensembl_bioperl](#binsrcensembl_bioperl)
  * [bin/src/ensembl-compara](#binsrcensembl-compara)
  * [bin/src/ensembl-funcgen](#binsrcensembl-funcgen)
  * [bin/src/ensembl-git-tools](#binsrcensembl-git-tools)
  * [bin/src/ensembl-io](#binsrcensemblio)
  * [bin/src/ensembl-tools-release-83](#binsrcensembl-tools-release-83)
  * [bin/src/ensembl-variation](#binsrcensembl-variation)

  * [bin/src/FastQC](#binsrcfastqc)
  * [bin/src/slate](#binsrcslate)
  * [bin/src/TrimGalore](#binsrctrimgalore)
  * [bin/src/Q](#binsrcq)
  * [bin/src/germline](#binsrcgermline)
  * [bin/src/spp](#binsrcspp)
  * [bin/src/ifluatex](#binsrcifluatex)
  * [bin/src/impute_v2.3.2_x86_64_static](#binsrcimpute_v232_x86_64_static)

  * [bin/src/liftover](#binsrcliftover)
  * [bin/src/rsem-1.2.20](#binsrcrsem-1220)
  * [bin/src/STAR-STAR_2.4.2a](#binsrcstar-star_242a)
  * [bin/src/aseq-v1.1.8](#binsrcaseq-v118)
  * [bin/src/bolt-lmm](#binsrcbolt-lmm)
  * [bin/src/EXTREME-2.0.0](#binsrcextyreme-200)
  * [bin/src/pandoc](#binsrcpandoc)

  * [bin/src/squid](#binsrcsquid)
  * [bin/src/rasqual](#binsrcrasqual)
  * [bin/src/BioPerl-1.6.1](#binsrcbioperl-161)
  * [bin/src/metaseq](#binsrcmetaseq)

  * [bin/src/phantompeakqualtools-read-only](#binsrcphantompeakqualtools-read-only)
  * [bin/src/phantompeakqualtools_spp1.13](#binsrcphantompeakqualtools_spp113)

  * [bin/src/bioperl-live](#binsrcbioperl-live)
  * [bin/src/plink](#binsrcplink)
  * [bin/src/kallisto](#binsrckallisto)
  * [bin/src/kallisto-0.42.2.1](#binsrckallisto-04221)
  * [bin/src/leafcutter](#binsrcleafcutter)
  * [bin/src/fusion_twas-master](#binsrcfusion_twas-master)
  * [bin/src/vcftools](#binsrcvcftools)
  * [bin/src/CrossMap-0.2.5](#binsrccrossmap-025)

  * [bin/src/stack](#binsrcstack)
  * [bin/src/picard](#binsrcpicard)
  * [bin/src/texlive](#binsrctexlive)

  * [bin/src/subread-1.4.6-p5-source](#binsrcsubread-146-p5-source)
  * [bin/src/bioperl](#binsrcbioperl)
  * [bin/src/WASP](#binsrcwasp)
  * [bin/src/idr](#binsrcidr)
  * [bin/src/gcta_1.24.7](#binsrcgcta_1247)
  * [bin/src/shapeit.v2](#binsrcshapeitv2)
  * [bin/src/subread-1.4.6-p5-Linux-x86_64](#binsrcsubread-146-p5-linux-x86_64)
  * [bin/src/leafcutter-master](#binsrcleafcutter-master)
  * [bin/src/eqtlbma](#binsrceqtlbma)
  * [bin/src/GEMMA](#binsrcgemma)
  * [bin/src/sleuth](#binsrcsleuth)
  * [bin/src/FastQTL_2.814](#binsrcfastqtl_2814)
  * [bin/src/bcftools](#binsrcbcftools)
  * [bin/src/peer](#binsrcpeer)
  * [bin/src/gatk](#binsrcgatk)
  * [bin/src/bedtools2](#binsrcbedtools2)
  * [bin/src/cutadapt-1.3](#binsrccutadapt-13)


<!--te-->




# bin/local

This directory contains many of the binaries or executable scripts. One may use them without path indication by including this directory in their $PATH variable. 


# bin/src

Source codes, etc. for the executables in the local directory are located here. I list reference websites, github repositories, etc. so one many look up the latest.


## bin/src/ensembl

## bin/src/ensembl_bioperl

## bin/src/ensembl-compara

## bin/src/ensembl-funcgen

## bin/src/ensembl-git-tools

## bin/src/ensembl-io

## bin/src/ensembl-tools-release-83

## bin/src/ensembl-variation


## bin/src/FastQC

## bin/src/slate

## bin/src/TrimGalore

## bin/src/Q

## bin/src/germline

## bin/src/spp

## bin/src/ifluatex

## bin/src/impute_v2.3.2_x86_64_static



## bin/src/liftover

## bin/src/rsem-1.2.20

## bin/src/STAR-STAR_2.4.2a

I primarily locally installed star for consistency in large-scale data processing (i.e. to prevent potential discrepancies between different subversions installed without explicit binary change). pmacs now maintains several versions of star as modules, which one can use by, e.g.

```
module load STAR-2.5.2a
```


## bin/src/aseq-v1.1.8

## bin/src/bolt-lmm

## bin/src/EXTREME-2.0.0

## bin/src/pandoc


## bin/src/squid

## bin/src/rasqual

## bin/src/BioPerl-1.6.1

## bin/src/metaseq



## bin/src/phantompeakqualtools-read-only

## bin/src/phantompeakqualtools_spp1.13



## bin/src/bioperl-live

## bin/src/plink

## bin/src/kallisto

## bin/src/kallisto-0.42.2.1

## bin/src/leafcutter

## bin/src/fusion_twas-master

## bin/src/vcftools

## bin/src/CrossMap-0.2.5



## bin/src/stack

## bin/src/picard

## bin/src/texlive



## bin/src/subread-1.4.6-p5-source

## bin/src/bioperl

## bin/src/WASP

## bin/src/idr

## bin/src/gcta_1.24.7

## bin/src/shapeit.v2

## bin/src/subread-1.4.6-p5-Linux-x86_64

## bin/src/leafcutter-master

## bin/src/eqtlbma

No longer locally installed. One may check pmacs for the latest and use relevant modules, e.g. by 

```
module load eQtlBma-1.3.1
````

## bin/src/GEMMA

## bin/src/sleuth

## bin/src/FastQTL_2.814

## bin/src/bcftools

## bin/src/peer

## bin/src/gatk

## bin/src/bedtools2

## bin/src/cutadapt-1.3




