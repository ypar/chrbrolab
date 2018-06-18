#!/usr/bin/env python3

# YoSon
# making tpm matrices from ldacc provided rpkm matrices

import os, sys
import pandas as pd
import logging

logging.basicConfig(filename='tpms_from_rpkms.log',level=logging.DEBUG)

# calculate tpms from rpkms

# test data
# tissuefile = 'v6p_fastQTL_FOR_QC_ONLY/Liver_Analysis.v6p.FOR_QC_ONLY.rpkm.gct'
tissuefile = sys.argv[1]

logging.info('tissue file is '+tissuefile)

tissuedf = pd.read_csv(tissuefile, sep='\t', header=0, skiprows=2)

if isinstance(tissuedf, pd.DataFrame):
  
  logging.info('finished loading '+tissuefile)
  
  # first change sample IDs to subject IDs
  for i, x in enumerate(tissuedf.columns):
    if x.startswith('GTE'):
      subid = x.rstrip().split('-')[0]+'-'+x.rstrip().split('-')[1]
      tissuedf.rename(columns={x: subid}, inplace=True)
  
  tissuedf.index = tissuedf.Name
  tissuedf.drop(['Name', 'Description'], axis=1, inplace=True)
  sums = tissuedf.sum(axis=0)
  
  tpms = tissuedf.divide(sums, axis='columns') * 1000000
  
  # assuming file name is in the format of the test file above
  outfile = tissuefile.replace('v6p_fastQTL_FOR_QC_ONLY', 'tpms_from_rpkms')
  outfile = outfile.replace('FOR_QC_ONLY.rpkm.gct', 'tpm.txt')
  tpms.to_csv(outfile, sep='\t', index=True, na_rep='NA')
  
  logging.info('finished writing '+outfile)
  
  gmeans = tpms.mean(axis=1)
  
  gmeansqc = gmeans[gmeans>1]
  
  outfile = outfile.replace('tpm.txt', 'tpm_min1.txt')
  tpmsqc = tpms.ix[gmeansqc.index]
  tpmsqc.to_csv(outfile, sep='\t', index=True, na_rep='NA')
  logging.info('finished writing '+outfile)
  
  gmeansqc = gmeans[gmeans>10]
  
  outfile = outfile.replace('tpm_min1.txt', 'tpm_min10.txt')
  tpmsqc = tpms.ix[gmeansqc.index]
  tpmsqc.to_csv(outfile, sep='\t', index=True, na_rep='NA')
  logging.info('finished writing '+outfile)
  



