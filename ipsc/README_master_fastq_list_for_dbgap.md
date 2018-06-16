
###
# README_master_fastq_list_for_dbgap
# YoSon Park
###

As of 2017 Dec, I have re-organized the ipsc project space on pmacs. raderlab users can find most raw files in raderlab project space. other collaborators may find relevant files in chrb_rader_lab project space. many downstream files such as summary statistics from Pashos, Park et al (2017) will be also copied over to chrbrolab project space. 


So far (Dec 2017), I have made the following changes to make this master list from merging final iPSC subjects list with three StudyInfo files provided by the NGS core.



# StudyInfo file condition (cell type) mislabelling or multi-labelling

There are five cell/tissue types used for the Pashos, Park et al (2017) CSC paper: iPSC, HLC, human primary hepatocytes, human liver (chicago) and human liver (gtex). GTEx samples were processed from fastq files I transferred from dbGaP. The other four were sequenced and processed at Penn. There were three different 'project' names these samples were assigned to internally: DREP, DRSE, and DRRB. As far as samples and all downstream processing and analyses are concerned, separating samples into these three projects does not seem necessary. As I was applying corrected sample IDs, I have also collapsed three project spaces to make shorter, more comprehensible directory structures.

* COND_name = ‘Human hepatocytes iPSC’ are changed to ‘HLC’ 
* COND_name = ‘hepatocyte’ are changed to ‘HLC’ (per Pashos, these were simply different labels used for HLC at the earlier phase of the project)
* COND_name = ‘primary’ are changed to ‘PrimaryHepatocyte’
* COND_name = ‘Human Liver’ are changed to ‘Liver’ but dropped from the dbGaP submission list 



# error in the StudyInfo file format

* New column ‘subjectid' is created by merging missing data in ‘SAMP_sid’ (lists subject IDs for most samples) with ‘SAMP_name’ (lists subject IDs for some samples) 



# error in sample submission (retrospective discovery from gene expression and genotype QC)

* Subject 781, COND_name= ‘hepatocyte’ in batch DREP-iPSHepat (but not in DRSE-EviPilot or DRRB-Liver, if any) are changed to ‘iPSC’ — said to be swapped at NGSC submission stage
* Subject 781, COND_name = ‘iPSC’ in batch DREP-iPSHepat (but not in DRSE-Evipilot or DRRB-Liver, if any) are changed to ‘HLC’ — said to be swapped at NGSC submission stage



# typoes in Study Info

* Subject id ‘126-1050’ is an incorrectly labelled ID of ‘126-1650’. The final list will replace all ‘126-1050’ to ‘126-1650’



# samples with multiple subject IDs

* Subjects with prefix ‘hiPS’ now have prefixes removed to be consistent with the remaining samples.
  * List of subject IDs from three StudyInfo files combined that do not match any sample dictionary subjects: 'hiPS 113-1138', 'hiPS 161', 'hiPS 225’, 'hiPS 245', 'hiPS 254', 'hiPS 286', 'hiPS 312', 'hiPS 316’, 'hiPS 321', 'hiPS 334', 'hiPS 340', 'hiPS 44-136', 'hiPS 486’,




# samples with multiple clone IDs, arbitrary NGSC core-assigned suffixes, possible typoes and excel errors.

* ‘M1’ is labelled as ‘M1-H10-1304’ in the genotypes. the final list will be matched to ‘M1’
* ‘MR43’ is labelled as ‘SV20’ in the genotypes. There are also other SV# and MR# in other records said to be the same subject. The final list will be matched to ‘MR43’

* There are several samples from the supplemental table 1 of the CSC paper that do not match directly to any StudyInfo subjects list. The following IDs are merged to the simpler (first column) IDs of the 2+ different IDs used from StudyInfo files to create my version of the master sample list:
  * '131-4128', '131-4128-1’, 
  * '134', '134-2’, 
  * '145', '145-2’, 
  * '16-1174', '16-1174-2’, 
  * '20-138', '20-138-2’, 
  * '25-1158', '25-1158-1’, 
  * '44', '44-136’, '44-1336'
  * '486', '486-1460-080’, 
  * '553', '553-1631’, 
  * '57-073', '57-073-2’, '57-1673’ <— this particular sample is an enigma with no relevant records. Assuming 0 turned into 16 somehow by a typo?
  * '60-099', '60-099-2’, 
  * '61-1559', '61-1559-2’, '61-1559-060'
  * '71-44', '71-44-1’, '71-044'
  * '75-1681', '75-1681-2’, 
  * '76-1544’, '76-1544-127’   



# formatting changes for convenience

* Hyphens from all raw fastq files are removed for raderlab directory storage. E.g. FGC1297_s_4_1_TCCGCGAATATAGCCT.fastq.gz rather than FGC1297_s_4_1_TCCGCGAA-TATAGCCT.fastq.gz
* All raw fastq files subdirectories as well as dbGaP submission will be based on dbGaP IDs rather than internal IDs. E.g. 126-1650, iPSC, would be referred by dbGaP ID Penn126i-120-1-iPSC. This was done due to cleaner/newer nature of the dbGaP IDs created after most other iterations/versions of the internal IDs have taken place.






