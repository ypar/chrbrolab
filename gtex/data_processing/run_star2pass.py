#!/usr/bin/env python3

# run_star2pass.py
# YoSon
# 01/02/2016

from sys import argv
import os, sys, subprocess
import pandas as pd

usage = """
usage:
python run_star2pass.py time sample_file scriptsuffix rootdir
"""

if len(argv) < 5:
    print usage
else:
    time = argv[1]
    samplefile = argv[2]
    mastername = str(argv[3])
    rootdir = str(argv[4])

master_script = rootdir + "/ase/low_level_script/star_2pass_hg19_genomealigned_" + str(mastername) + ".sh"
master_handle=open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")
workdir = rootdir + '/data/phenotype/expression/mapped_rna_seq_reads'

samples = pd.read_csv(samplefile, sep='\t', header=False, dtype=object, names=['samples'])

for i, row in samples.iterrows():
    sample = row['samples']
    bsubbatchfile = rootdir + "/ase/low_level_script/star_2pass_hg19_genomealigned_" + str(sample) + "_bsub"
    bsubbatchhandle=open(bsubbatchfile, 'w')
    # create the batch script
    bsubbatchlist = []
    header=r"""#!/bin/bash
#BSUB -J star_2pass_genomealigned_%s
#BSUB -o low_level_log/star_2pass_hg19_genomealigned_%s.o
#BSUB -e low_level_log/star_2pass_hg19_genomealigned_%s.e
#BSUB -M 65000
#BSUB -n 6
#BSUB -R 'span[ptile=6]'
#BSUB -t %s
"""%(sample, sample, sample, time)
    bsubbatchlist.append(header)
    #for sample_name in samples:
    bsubbatch="""
rootdir='%s'
sample='%s'

module load STAR-2.5.2a

# genome dir
genomedir="${rootdir}/data/auxiliary/alignment/star_2pass_hg19"

sampledir="${rootdir}/data/phenotype/expression/trimmed_rna_seq_reads"
wkdir="${rootdir}/data/phenotype/expression/mapped_rna_seq_reads/star_2pass_hg19_genomealigned/${sample}/"
mkdir -p -m 770 ${wkdir}
cd ${wkdir}

#bam file from star2pass
gstarbam="Aligned.toGenome.bam"

echo "sample is $sample" >> ${wkdir}/${sample}_genomealigned.log

if [[ ! -s ${sampledir}/${sample}_1.trimmed.P.fastq || ! -s ${sampledir}/${sample}_2.trimmed.P.fastq ]]; then

# star decompression is slow. unzip first before processing.
  bzip2 -dc ${sampledir}/${sample}_1.trimmed.P.fastq.bz2 > ${sampledir}/${sample}_1.trimmed.P.fastq
  bzip2 -dc ${sampledir}/${sample}_2.trimmed.P.fastq.bz2 > ${sampledir}/${sample}_2.trimmed.P.fastq

fi

if [[ ! -s Aligned.toGenome.Aligned.out.sam && ! -s Aligned.toGenome.bam ]]; then
# for featurecount quants

/opt/software/STAR/2.5.2a/bin/STAR  \
--genomeDir $genomedir \
--runThreadN 6 \
--genomeLoad LoadAndKeep \
--readFilesIn \
${sampledir}/${sample}_1.trimmed.P.fastq \
${sampledir}/${sample}_2.trimmed.P.fastq \
--outSAMattributes NH HI AS NM MD \
--outFilterIntronMotifs RemoveNoncanonical \
--outFileNamePrefix Aligned.toGenome.
fi


if [[ -s Aligned.toGenome.Aligned.out.sam && ! -s Aligned.toGenome.bam ]]; then
  ~/bin/src/anaconda3/bin/samtools view -bS Aligned.toGenome.Aligned.out.sam > Aligned.toGenome.bam
fi

if [[ -s Algined.toGenome.bam && ! -s Aligned.toGenome.bam.flagstat ]]; then
  ~/bin/src/anaconda3/bin/samtools flagstat Algined.toGenome.bam > Aligned.toGenome.bam.flagstat
fi

if [[ -s Aligned.toGenome.Aligned.out.sam && -s Aligned.toGenome.bam ]]; then
  rm Aligned.toGenome.Aligned.out.sam
fi

if [ ! -s Aligned.q10.bam ]; then
  #perform post star processing qc
  # q10 filter step
  ~/bin/src/anaconda3/bin/samtools view -@ 6 -q 10 -b ${gstarbam} > Aligned.q10.bam
fi

if [ -s Aligned.q10.bam ]; then
  echo "STAR q10 filter run for ${sample} done" >> ${wkdir}/${sample}_genomealigned.log
fi

if [[ -s Aligned.q10.bam && ! -s Aligned.q10.sort.bam.flagstat ]];then
  ~/bin/src/anaconda3/bin/samtools flagstat Aligned.q10.bam > Aligned.q10.bam.flagstat
  ~/bin/src/anaconda3/bin/samtools sort -@ 6 -m 10G -T Aligned.q10.sort.temp -o Aligned.q10.sort.bam Aligned.q10.bam
  ~/bin/src/anaconda3/bin/samtools flagstat Aligned.q10.sort.bam > Aligned.q10.sort.bam.flagstat
  mv Aligned.q10.sort.bam Aligned.q10.bam
fi


"""%(rootdir, sample)
    bsubbatchlist.append(bsubbatch)
    for x in bsubbatchlist:
        bsubbatchhandle.write(x)
    bsubbatchhandle.close()
    samplefile = workdir + '/star_2pass_hg19_genomealigned/' + sample + '/Aligned.out.bam'
    outfile = workdir + '/star_2pass_hg19_genomealigned/' + sample + '/Aligned.q10.bam'
    if not os.path.isfile(outfile):
        master_handle.write("bsub < " + bsubbatchfile  + " \n")



