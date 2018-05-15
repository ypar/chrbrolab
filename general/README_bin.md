
# README for the lab-wide software and resource directories
# YoSon Park

### last modified 05/14/2018


For general use by chrbrolab project space users, I have, in the past, downloaded, installed and organized several commonly used tools for many lab personnel and collaborators. I no longer maintain these directories but I highly encourage others to keep following the directory organizing convention. I do not install or maintain user-specific resources for other users, e.g. R packages on user-specific environments. Also hpc pmacs common environment now (2017 and later) maintains much of the hard-to-install software or modules. I list these subdirectories here so one may look up relevant websites and references for their own applications. Index is provided for quick lookups but one may search for specific software as well.

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

  * [bin/src/cutadapt-1.3](#binsrccutadapt-13)
  * [bin/src/FastQC](#binsrcfastqc)
  * [bin/src/TrimGalore](#binsrctrimgalore)

  * [bin/src/Q](#binsrcq)
  * [bin/src/germline](#binsrcgermline)

  * [bin/src/slate](#binsrcslate)
  * [bin/src/spp](#binsrcspp)

  * [bin/src/ifluatex](#binsrcifluatex)

  * [bin/src/impute_v2.3.2_x86_64_static](#binsrcimpute_v232_x86_64_static)
  * [bin/src/shapeit.v2](#binsrcshapeitv2)
  

  * [bin/src/liftover](#binsrcliftover)
  * [bin/src/rsem-1.2.20](#binsrcrsem-1220)
  * [bin/src/STAR-STAR_2.4.2a](#binsrcstar-star_242a)
  * [bin/src/aseq-v1.1.8](#binsrcaseq-v118)
  * [bin/src/bolt-lmm](#binsrcbolt-lmm)
  * [bin/src/EXTREME-2.0.0](#binsrcextyreme-200)

  * [bin/src/pandoc](#binsrcpandoc)

  * [bin/src/vcftools](#binsrcvcftools)

  * [bin/src/squid](#binsrcsquid)
  * [bin/src/rasqual](#binsrcrasqual)

  * [bin/src/bioperl-live](#binsrcbioperl-live)
  * [bin/src/bioperl](#binsrcbioperl)  
  * [bin/src/BioPerl-1.6.1](#binsrcbioperl-161)
  
  * [bin/src/metaseq](#binsrcmetaseq)

  * [bin/src/phantompeakqualtools-read-only](#binsrcphantompeakqualtools-read-only)
  * [bin/src/phantompeakqualtools_spp1.13](#binsrcphantompeakqualtools_spp113)

  * [bin/src/plink](#binsrcplink)

  * [bin/src/kallisto](#binsrckallisto)
  * [bin/src/kallisto-0.42.2.1](#binsrckallisto-04221)
  * [bin/src/leafcutter](#binsrcleafcutter)
  * [bin/src/fusion_twas-master](#binsrcfusion_twas-master)

  * [bin/src/CrossMap-0.2.5](#binsrccrossmap-025)

  * [bin/src/stack](#binsrcstack)
  * [bin/src/picard](#binsrcpicard)
  * [bin/src/texlive](#binsrctexlive)

  * [bin/src/subread-1.4.6-p5-source](#binsrcsubread-146-p5-source)

  * [bin/src/WASP](#binsrcwasp)
  * [bin/src/idr](#binsrcidr)
  * [bin/src/gcta_1.24.7](#binsrcgcta_1247)
  
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


<!--te-->




# bin/local

This directory contains many of the binaries or executable scripts. One may use them without path indication by including this directory in their $PATH variable. 


# bin/src

Source codes, etc. for the executables in the local directory are located here. I list reference websites, github repositories, etc. so one many look up the latest.


## bin/src/ensembl

Most ensembl tools were used to test for the [primates project](https://github.com/ypar/chrbrolab/blob/master/primates/README_primates.md). However, 

## bin/src/ensembl_bioperl

## bin/src/ensembl-compara

## bin/src/ensembl-funcgen

## bin/src/ensembl-git-tools

## bin/src/ensembl-io

## bin/src/ensembl-tools-release-83

## bin/src/ensembl-variation




## bin/src/cutadapt-1.3



## bin/src/FastQC

FastQC is a quality control tool for high throughput sequence data and any updates can be checked [here](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). It requires a suitable java runtime enviroment and [picard](http://broadinstitute.github.io/picard/). The following is a description of what FastQC can do from the official website:

> The main functions of FastQC are

>> Import of data from BAM, SAM or FastQ files (any variant)
>> Providing a quick overview to tell you in which areas there may be problems
>> Summary graphs and tables to quickly assess your data
>> Export of results to an HTML based permanent report
>> Offline operation to allow automated generation of reports without running the interactive application





## bin/src/TrimGalore

TrimGalore! is a perl script wrapper for [Cutadapt](https://github.com/marcelm/cutadapt) and [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/). Any updates and detailed documentations may be found [here](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) and also in its [github repo](https://github.com/FelixKrueger/TrimGalore). This is an excerpt from the official website:

> Trim Galore! is a wrapper script to automate quality and adapter trimming as well as quality control, with some added functionality to remove biased methylation positions for RRBS sequence files (for directional, non-directional (or paired-end) sequencing). It's main features are:

> Trim Galore is now also available from GitHub. You are invited to leave comments, feature request or bug reports over there!

>> For adapter trimming, Trim Galore! uses the first 13 bp of Illumina standard adapters ('AGATCGGAAGAGC') by default (suitable for both ends of paired-end libraries), but accepts other adapter sequence, too
>> For MspI-digested RRBS libraries, Trim Galore! performs quality and adapter trimming in two subsequent steps. This allows it to remove 2 additional bases that contain a cytosine which was artificially introduced in the end-repair step during the library preparation
>> For any kind of FastQ file other than MspI-digested RRBS, Trim Galore! can perform single-pass adapter- and quality trimming
>> The Phred quality of basecalls and the stringency for adapter removal can be specified individually
>> Trim Galore! can remove sequences if they become too short during the trimming process. For paired-end files Trim Galore! removes entire sequence pairs if one (or both) of the two reads became shorter than the set length cutoff. Reads of a read-pair that are longer than a given threshold but for which the partner read has become too short can optionally be written out to single-end files. This ensures that the information of a read pair is not lost entirely if only one read is of good quality
>> Trim Galore! can trim paired-end files by 1 additional bp from the 3' end of all reads to avoid problems with invalid alignments with Bowtie 1
>> Trim Galore! accepts and produces standard or gzip compressed FastQ files
>> FastQC can be run on the resulting output files once trimming has completed (optional)




## bin/src/Q

## bin/src/germline



## bin/src/slate




## bin/src/spp

I do not maintain installation of spp for the lab. One may install it for their favorite R version by using devtools.


```
require(devtools)
devtools::install_github('hms-dbmi/spp', build_vignettes = FALSE)
```

Following features are listed for spp on the [official website](http://compbio.med.harvard.edu/Supplements/ChIP-seq/):

> Features
>> Assess overall DNA-binding signals in the data and select appropriate quality of tag alignment.
>> Discard or restrict positions with abnormally high number of tags.
>> Calculate genome-wide profiles of smoothed tag density and save them in WIG files for viewing in other browsers.
>> Calculate genome-wide profiles providing conservative statistical estimates of fold enrichment ratios along the genome. These can be exported for browser viewing, or thresholded to determine regions of significant enrichment/depletion.
>> Determine statistically significant point binding positions
>> Assess whether the set of point binding positions detected at a current sequencing depth meets saturation criteria, and if does not, estimate what sequencing depth would be required to do so.


The following paper provides the rationale for the methods implemented in spp:
[Kharchenko PK, Tolstorukov MY, Park PJ "Design and analysis of ChIP-seq experiments for DNA-binding proteins" Nat. Biotech. doi:10.1038/nbt.1508](https://www.nature.com/articles/nbt.1508)


## bin/src/ifluatex




## bin/src/impute_v2.3.2_x86_64_static



## bin/src/shapeit.v2





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

Pandoc is a swiss knife among document converter tools. Detailed documentations can be found [here](https://pandoc.org/MANUAL.html). I highly recommend [anaconda](https://www.anaconda.com/download/#linux) or similarly packaged versions of python and set of modules useful for computing needs. The latest anaconda includes pandoc among hundreds of other useful packages so for most people using python, the standalone installation of pandoc is no longer necessary. The full list of packages installed by anaconda distributed version of python is [here](https://docs.anaconda.com/anaconda/packages/py3.6_linux-64).

Following is an excerpt from the official website:

> If you need to convert files from one markup format into another, pandoc is your swiss-army knife. Pandoc can convert documents in (several dialects of) Markdown, reStructuredText, textile, HTML, DocBook, LaTeX, MediaWiki markup, TWiki markup, TikiWiki markup, Creole 1.0, Vimwiki markup, OPML, Emacs Org-Mode, Emacs Muse, txt2tags, Microsoft Word docx, LibreOffice ODT, EPUB, or Haddock markup to

> HTML formats
>> XHTML, HTML5, and HTML slide shows using Slidy, reveal.js, Slideous, S5, or DZSlides

> Word processor formats
>> Microsoft Word docx, OpenOffice/LibreOffice ODT, OpenDocument XML, Microsoft PowerPoint.

> Ebooks
>> EPUB version 2 or 3, FictionBook2

> Documentation formats
>> DocBook version 4 or 5, TEI Simple, GNU TexInfo, Groff man, Groff ms, Haddock markup

> Archival formats
>> JATS

> Page layout formats
>> InDesign ICML

> Outline formats
>> OPML

> TeX formats
>> LaTeX, ConTeXt, LaTeX Beamer slides

> PDF
>> via pdflatex, xelatex, lualatex, pdfroff, wkhtml2pdf, prince, or weasyprint.

> Lightweight markup formats
>> Markdown (including CommonMark and GitHub-flavored Markdown), reStructuredText, AsciiDoc, Emacs Org-Mode, Emacs Muse, Textile, txt2tags, MediaWiki markup, DokuWiki markup, TikiWiki markup, TWiki markup, Vimwiki markup, and ZimWiki markup.

> Custom formats
>> custom writers can be written in lua.



## bin/src/squid


## bin/src/rasqual


## bin/src/bioperl

## bin/src/BioPerl-1.6.1

## bin/src/bioperl-live



## bin/src/metaseq



## bin/src/phantompeakqualtools-read-only

## bin/src/phantompeakqualtools_spp1.13




## bin/src/plink

## bin/src/kallisto

## bin/src/kallisto-0.42.2.1

## bin/src/leafcutter

Leafcutter paper is now [published](https://www.nature.com/articles/s41588-017-0004-9).


## bin/src/fusion_twas-master

## bin/src/vcftools




## bin/src/CrossMap-0.2.5



## bin/src/stack

## bin/src/picard

## bin/src/texlive



## bin/src/subread-1.4.6-p5-source


## bin/src/WASP

## bin/src/idr

## bin/src/gcta_1.24.7





## bin/src/subread-1.4.6-p5-Linux-x86_64

## bin/src/leafcutter-master


## bin/src/FastQTL_2.814



## bin/src/eqtlbma

I no longer maintain the installation of eqtlbma for the lab. One may check pmacs for the latest and use relevant modules, e.g. by 

```
module load eQtlBma-1.3.1

```


## bin/src/GEMMA



## bin/src/sleuth


## bin/src/bcftools

## bin/src/peer

## bin/src/gatk

## bin/src/bedtools2




