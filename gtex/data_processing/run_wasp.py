#!/usr/bin/env python3

# run_wasp.py
# YoSon
# 01/02/2016

from sys import argv
import os, sys, subprocess
import pandas as pd

usage = """
usage:
python run_wasp.py time sample_file scriptsuffix rootdir
"""

if len(argv) < 4:
    print usage
else:
    time = argv[1]
    samplefile = argv[2]
    mastername = str(argv[3])
    rootdir = str(argv[4])

master_script = rootdir + "/ase/low_level_script/star_2pass_hg19_wasp_" + str(mastername) + ".sh"
master_handle=open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")
workdir = rootdir + '/data/phenotype/expression/mapped_rna_seq_reads'

samples = pd.read_csv(samplefile, sep='\t', header=False, dtype=object, names=['samples'])

for i, row in samples.iterrows():
    sample = row['samples']
    bsubbatchfile = rootdir + "/ase/low_level_script/star_2pass_hg19_wasp_" + str(sample) + "_bsub"
    bsubbatchhandle=open(bsubbatchfile, 'w')
    # create the batch script
    bsubbatchlist = []
    header=r"""#!/bin/bash
#BSUB -J star_2pass_hg19_wasp_bsub_%s
#BSUB -o low_level_log/star_2pass_hg19_wasp_bsub_%s.o
#BSUB -e low_level_log/star_2pass_hg19_wasp_bsub_%s.e
#BSUB -M 65000
#BSUB -n 6
#BSUB -R 'span[ptile=6]'
#BSUB -t %s
"""%(sample, sample, sample, time)
    bsubbatchlist.append(header)
    bsubbatch="""
rootdir='%s'
sample='%s'

wkdir="${rootdir}/data/phenotype/expression/mapped_rna_seq_reads/star_2pass_hg19_wasp/${sample}/"
bamdir="${rootdir}/data/phenotype/expression/mapped_rna_seq_reads/star_2pass_hg19_genomealigned/${sample}/"

mkdir -p -m 770 ${wkdir}
cd ${bamdir}

echo "job begun" >> ${wkdir}/${sample}_wasp_preproc.log

module load STAR-2.5.2a

#wasp ref file
waspdir="${rootdir}/ase/snpdir"
# genome dir
genomedir="${rootdir}/data/auxiliary/alignment/star_2pass_hg19"

echo "sample is $sample" >> ${wkdir}/${sample}_wasp_preproc.log

#bam file from star2pass
starbam="${rootdir}/data/phenotype/expression/mapped_rna_seq_reads/star_2pass_hg19_genomealigned/${sample}/Aligned.q10.bam"

if [ ! -s ${wkdir}/Aligned.q10.wasp.rmdup.bam.bai ]; then

  if [[ -s ${starbam} && ! -s ${starbam}.bai ]]; then
    ~/bin/src/anaconda3/bin/samtools index ${starbam}
  fi

  if [[ -s ${starbam} && -s ${starbam}.bai ]]; then

    rm -f ${wkdir}/${sample}_wasp_preproc.log
    echo "${starbam} and ${starbam}.bai exist" >> ${wkdir}/${sample}_wasp_preproc.log
    echo "${sample} run_wasp job is initiated" >> ${wkdir}/${sample}_wasp_preproc.log
    # rm -f Aligned.toGenome.bam.flagstat
    # # collect stats on STAR bam file
    # ~/src/anaconda3/bin/samtools flagstat ${starbam} > Aligned.out.bam.flagstat
    if [ ! -s ${starbam} ]; then
      #perform post star processing qc
      # q10 filter step
      ~/bin/src/anaconda3/bin/samtools view -@ 6 -q 10 -b ${bamdir}/Aligned.toGenome.bam > ${starbam}
    fi

    # if [ -s Aligned.q10.bam ]; then
    #   echo "STAR q10 filter run for ${sample} done" >> ${wkdir}/${sample}_wasp_preproc.log
    # fi
    #
    # if [[ -s Aligned.q10.bam && ! -s Aligned.q10.sort.bam ]];then
    #   ~/bin/src/anaconda3/bin/samtools flagstat ${bamdir}/Aligned.q10.bam > ${bamdir}/Aligned.q10.bam.flagstat
    #   ~/bin/src/anaconda3/bin/samtools sort -@ 6 -m 10G -T ${bamdir}/Aligned.q10.sort.temp -o ${bamdir}/Aligned.q10.sort.bam ${bamdir}/Aligned.q10.bam
    #   ~/bin/src/anaconda3/bin/samtools flagstat ${bamdir}/Aligned.q10.sort.bam > ${bamdir}/Aligned.q10.sort.bam.flagstat
    #   mv ${bamdir}/Aligned.q10.sort.bam ${bamdir}/Aligned.q10.bam
    # fi

    if [[ ! -s Aligned.q10.remap.fq1.gz && ! -s Aligned.q10.to.remap.num.gz ]]; then
      # STEP 2: identify reads that need to be remapped
      python2 ~/bin/src/WASP/mapping/find_intersecting_snps.py -s -p ${bamdir}/Aligned.q10.bam $waspdir
    fi

    if [[ -s Aligned.q10.remap.fq1.gz && -s Aligned.q10.remap.fq2.gz ]]; then
      if [ ! -s Aligned.q10.remap.bam ]; then

        # STEP 3: Remap suspect reads
        # star zip deconpressing step is quite slow. unzip for processing
        gzip -dc Aligned.q10.remap.fq1.gz > Aligned.q10.remap.fq1
        gzip -dc Aligned.q10.remap.fq2.gz > Aligned.q10.remap.fq2

        /opt/software/STAR/2.5.2a/bin/STAR \
        --runThreadN 6 \
        --genomeLoad LoadAndKeep \
        --genomeDir $genomedir \
        --readFilesIn Aligned.q10.remap.fq1 Aligned.q10.remap.fq2 \
        --outSAMattributes NH HI AS NM MD \
        --alignSJDBoverhangMin 1 \
        --outFilterIntronMotifs RemoveNoncanonical \
        --outFileNamePrefix Aligned.q10.remap.

        ~/bin/src/anaconda3/bin/samtools view -@ 6 -Sb -q 10 Aligned.q10.remap.Aligned.out.sam > Aligned.q10.remap.bam

        if [ -s Aligned.q10.remap.bam ]; then
          echo "remap for ${sample} done" >> ${wkdir}/${sample}_wasp_preproc.log
        fi

      fi

      if [ ! -s Aligned.q10.wasp.rmdup.bam ]; then

        if [ ! -s Aligned.q10.wasp.bam ]; then
          python2 ~/bin/src/WASP/mapping/filter_remapped_reads.py -p Aligned.q10.to.remap.bam Aligned.q10.remap.bam Aligned.q10.remap.keep.bam Aligned.q10.to.remap.num.gz
          # Merge output
          ~/bin/src/anaconda3/bin/samtools merge -f -@ 6 Aligned.q10.wasp.bam Aligned.q10.remap.keep.bam Aligned.q10.keep.bam
        fi

        if [[ -s Aligned.q10.wasp.bam && ! -s Aligned.q10.remap.keep.bam.flagstat ]]; then
          echo "WASP and Merging for ${sample} done" >> ${wkdir}/${sample}_wasp_preproc.log
          rm -f Aligned.q10.wasp.bam.flagstat
          ~/bin/src/anaconda3/bin/samtools flagstat Aligned.q10.wasp.bam > Aligned.q10.wasp.bam.flagstat
          rm -f Aligned.q10.remap.out.sam

          if [[ -s Aligned.q10.remap.fq1.gz && -s Aligned.q10.remap.fq2.gz ]]; then
            rm -f Aligned.q10.remap.fq1 Aligned.q10.remap.fq2
          fi

          ~/bin/src/anaconda3/bin/samtools flagstat Aligned.q10.remap.bam > Aligned.q10.remap.bam.flagstat
          ~/bin/src/anaconda3/bin/samtools flagstat Aligned.q10.to.remap.bam > Aligned.q10.to.remap.bam.flagstat
          ~/bin/src/anaconda3/bin/samtools flagstat Aligned.q10.keep.bam > Aligned.q10.keep.bam.flagstat
          ~/bin/src/anaconda3/bin/samtools flagstat Aligned.q10.remap.keep.bam > Aligned.q10.remap.keep.bam.flagstat
        fi

        if [ ! -s Aligned.q10.wasp.sort.bam ]; then
          # Sort BAM
          ~/bin/src/anaconda3/bin/samtools sort -@ 6 Aligned.q10.wasp.bam Aligned.q10.wasp.sort
        fi

        if [ -s Aligned.q10.wasp.sort.bam ]; then
          echo "Sort for ${sample} done" >> ${wkdir}/${sample}_wasp_preproc.log
          ~/bin/src/anaconda3/bin/samtools flagstat Aligned.q10.wasp.sort.bam > Aligned.q10.wasp.sort.bam.flagstat
          mv Aligned.q10.wasp.sort.bam Aligned.q10.wasp.bam
        fi

        if [ ! -s Aligned.q10.wasp.rmdup.bam ]; then
          # remove dups
          python2 ~/bin/src/WASP/mapping/rmdup_pe.py Aligned.q10.wasp.bam Aligned.q10.wasp.rmdup.bam
        fi

        if [ -s Aligned.q10.wasp.rmdup.bam ]; then
          echo "Duplicate removal for ${sample} done" >> ${wkdir}/${sample}_wasp_preproc.log

          if [ ! -s Aligned.q10.wasp.rmdup.bam.bai ]; then
            # ~/bin/src/anaconda3/bin/samtools sort -@ 6 Aligned.q10.wasp.rmdup.bam Aligned.q10.wasp.rmdup.sort
            ~/bin/src/anaconda3/bin/samtools sort -@ 6 -m 10G -T Aligned.q10.wasp.rmdup.sort.temp -o Aligned.q10.wasp.rmdup.sort.bam Aligned.q10.wasp.rmdup.bam
            mv Aligned.q10.wasp.rmdup.sort.bam Aligned.q10.wasp.rmdup.bam
            ~/bin/src/anaconda3/bin/samtools index Aligned.q10.wasp.rmdup.bam
          fi
          rm -f ${wkdir}/Aligned.q10.wasp.rmdup.bam.flagstat Aligned.q10.wasp.rmdup.bam.flagstat
          ~/bin/src/anaconda3/bin/samtools flagstat Aligned.q10.wasp.rmdup.bam > Aligned.q10.wasp.rmdup.bam.flagstat

        fi

      fi
    fi
  fi
fi

# if [[ -s ${wkdir}/Aligned.q10.wasp.rmdup.bam && -s ${wkdir}/Aligned.q10.wasp.rmdup.bam.bai ]]; then
#   echo "Cleaning up..." >> ${wkdir}/${sample}_wasp_preproc.log
#   rm -f Aligned.q10.wasp.sort.bam Aligned.q10.bam Aligned.q10.keep.bam Aligned.q10.remap.bam Aligned.q10.remap.keep.bam
#   rm -f Aligned.q10.remap.fq1 Aligned.q10.remap.fq2
#   mv * ${wkdir}
#   chmod -f -R 770 ${wkdir}
# fi

"""%(rootdir, sample)
    bsubbatchlist.append(bsubbatch)
    for x in bsubbatchlist:
        bsubbatchhandle.write(x)
    bsubbatchhandle.close()
    samplefile = workdir + '/star_2pass_hg19_genomealigned/' + sample + '/Aligned.toGenome.bam'
    outfile = workdir + '/star_2pass_hg19_wasp/' + sample + '/Aligned.q10.wasp.rmdup.bam'
    if os.path.isfile(samplefile) and not os.path.isfile(outfile):
        master_handle.write("bsub < " + bsubbatchfile  + " \n")

