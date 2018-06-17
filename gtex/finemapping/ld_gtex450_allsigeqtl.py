#!/usr/bin/env python3

###
# YoSon
# 03/11/2016
# ld estimate for all significant eqtl snps + ld 0.000001 in gtex vcf
###

import pandas as pd
import subprocess
import os, sys, csv
from sys import argv

# usage: python3 ld_gtex450_allsigeqtl.py ${eqtl_file}

# example eqtl file

# head -2 Brain_Hippocampus_Analysis.snpgenes
# snp	gene	beta	t_stat	se	p_value	nom_thresh	min_p	gene_emp_p	k	n	gene_q_value	beta_noNorm	snp_chrom	snp_pos	minor_allele_samples	minor_allele_count	maf	ref_factor	ref	alt	snp_id_1kg_project_phaseI_v3	rs_id_dbSNP142_GRCh37p13	num_alt_per_site	has_best_p	is_chosen_snp	gene_name	gene_source	gene_chr	gene_start	gene_stop	orientation	tss_position	gene_type	gencode_attributes	tss_distance
# 3_132424172_C_G_b37	ENSG00000113971.14	-0.63439639754645	-7.81684582255339	0.08115759	1.11705996506884e-10	1.8524778E-6	5.0791446E-12	9.999E-5	0	10000	0.003375747	-1.33222274393573	3	132424172	22	23	0.14197531	1	C	G	rs10935029	rs10935029	1	0	0	NPHP3	HAVANA	3	132339976	132441303	-	132441303protein_coding	gene_id "ENSG00000113971.14"; transcript_id "ENSG00000113971.14"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "NPHP3"; transcript_type "nonsense_mediated_decay"; transcript_status "KNOWN"; transcript_name "NPHP3-001"; level 2; havana_gene "OTTHUMG00000159713.3"; havana_transcript "OTTHUMT00000357015.1";	17131


_wkdir = '/project/chrbrolab/gtex/results/GTEx_Analysis_v6p_eQTL_ld_proxies/'
_genodir = '/project/chrbrolab/gtex/data/GTEx_v6p_genotype_phased/'

_eqtl = argv[1]

# assuming file name follows the standard I use as an example above
# extract tissue name 
_tissue = _eqtl.split('/')[-1]
_tissue = _tissue.replace('_Analysis.snpgenes', '')
_masterfile = _wkdir + _tissue + 'ld_gtex_eqtl_master_bsub'

if not os.path.isdir(_wkdir):
  os.makedirs(_wkdir)

if not os.path.isdir(_tissue):
  os.makedirs(_tissue)

if not os.path.isdir(_tissue+'/low_level_script/'):
  os.makedirs(_tissue+'/low_level_script/')

if not os.path.isdir(_tissue+'/low_level_log/'):
  os.makedirs(_tissue+'/low_level_log/')

_eqtldf = pd.read_csv(_eqtl, sep='\t', header=0)
_snps = pd.unique(_eqtldf.snp.ravel())

for _snp in _snps:
  
  _chrom = _eqtldf.loc[_eqtldf['snp']==_snp].iloc[0]['snp_chrom']
  _vcffile = _genodir + 'GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_chr' + str(_chrom) + '_phased_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz'
  
  _ldfile = _wkdir + str(_tissue) + '/' +  str(_tissue) + '_chr' + str(_snp) + '_gtex450_r2_1e-6'
  _bsubfile = str(_ldfile) + '.bsub'
  _bsubfile = str(_bsubfile).replace(_tissue+'/', _tissue+'/low_level_script/')
  
  with open(_bsubfile, 'w') as _bsub:
    
    _stderr = _bsubfile + '.e'
    _stderr = _stderr.replace('script', 'log')
    _stdout = _bsubfile + '.o'
    _stdout = _stdout.replace('script', 'log')
    _jobid = _bsubfile.split('/')[-1]
    
    _header = """#!/bin/bash
#BSUB -e {_stderr}
#BSUB -o {_stdout}
#BSUB -J {_jobid}
#BSUB -M 5500

"""
    _params = {
      "_stderr" : _stderr,
      "_stdout" : _stdout,
      "_jobid" : _jobid
    }
    _cmd = 'plink --vcf ' + str(_vcffile) + ' --memory 4500 --r2 with-freqs --ld-window-r2 0.000001 --ld-snp ' + str(_snp) + ' --out ' + str(_ldfile)
    
    _bsub.write(_header.format(**_params))
    _bsub.write(_cmd)
  
  _call = 'sh ' + _bsubfile
  
  with open(_masterfile, 'a') as _master:
    _master.write(_call)
  
  
  #subprocess.call(_call, shell=True)
  #subprocess.Popen(_bsubfile)



