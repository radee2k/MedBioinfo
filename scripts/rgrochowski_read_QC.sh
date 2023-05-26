#!/bin/bash

# 2023-05, Radoslaw Grochowski
# This is a script documenting my workflow for assignment 10 of MedBioInfo Apllied Bioinformatics course.

date

# A minimal set of commands to reprodcue my analyses:

# check what patients were assigned to rgrochowski
# sqlite3 -batch /shared/ifbstor1/projects/2314_medbioinfo/pascal/central_database/sample_collab.db 'select run_accession from sample_annot left join sample2bioinformatician using(patient_code) where username is "rgrochowski";'

# run_accession
# -------------
# ERR6913112   
# ERR6913113   
# ERR6913320   
# ERR6913209   
# ERR6913122
# ERR6913208   
# ERR6913218   
# ERR6913349   

# Write the output to a file used to download assigned files
# sqlite3 -batch -noheader -csv /shared/ifbstor1/projects/2314_medbioinfo/pascal/central_database/sample_collab.db 'select run_accession from sample_annot left join sample2bioinformatician using(patient_code) where username is "rgrochowski";' > ../analyses/rgrochowski_run_accessions.txt

# load the api for sequence download from NCBI
# module load sra-tools

# configure the sra-toolkit
# vdb-config --interactive

# fastq-dump --gzip --split-3 --readids --disable-multithreading -X 10 --outdir /shared/home/rgrochowski/medbioinfo_folder/radoslaw/MedBioinfo/data/sra_fastq/ ERR6913113

# cat rgrochowski_run_accessions.txt | xargs fastq-dump --gzip --split-3 --readids --disable-multithreading -X 10 --outdir /shared/home/rgrochowski/medbioinfo_folder/radoslaw/MedBioinfo/data/sra_fastq/

cd ~/medbioinfo_folder/radoslaw/MedBioinfo/data/sra_fastq; for file in ./* ; do ( echo $file ; zcat $file | grep ^@ | wc -l ) >> ../../analyses/seq_count.txt ; done


 srun -t 30 -c 1 seqkit stats -j 1 ./* >> ../../analyses/seq_count_seqkit.txt

# the reads seem to be trimmed and quality filtered as nearly all the bases have the same, good quality score and inspected reads had no low quality ends and adapters present and about standard size (150bp)

srun -t 30 -c 1 ls | grep gz$ | xargs -I {} seqkit rmdup -j 1 {} -D ../../analyses/seq_count_dup.txt -o dedup/{}

# output
# [INFO] 0 duplicated records removed 16 times and the seq_count_dup.txt is empty

# reads weren't duplicated

srun -t 30 -c 4 seqkit locate -j 4 ./*.gz -p AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA >> ../../analyses/trimming.txt

# no adaptors were found among the sequences

# more complicated but a nice way of joining stdout 
# srun -c 4 -t 30 ls ./*.gz | xargs -i sh -c ' echo {} ; seqkit locate -j 4 {} -p GAGCACACGTCTGAACTCCAGTCA,GAGCGTCGTGTAGGGAAAGAGTGT'  >> ../../analyses/trimming_locate2.txt

# looking for shorter adaptors
srun -t 30 -c 4 seqkit locate -j 4 ./*.gz -p GAGCACACGTCTGAACTCCAGTCA,GAGCGTCGTGTAGGGAAAGAGTGT >> ../../analyses/trimming2.txt




 srun -c 4 -t 30 ls ./*.gz | xargs -I {} seqkit locate -j 4 {} -p AACTCCAGTCA,GGAAAGAGTGT >> ../../analyses/trimming_locate3.txt

# MAKE COMMENTS

 module fastqc 
 module seqkit

 srun -t 30 -c 4 fastqc --noextract -o analyses/fastqc/ -t 4 data/sra_fastq/ERR6913122_1.fastq.gz data/sra_fastq/ERR6913122_2.fastq.gz
/n srun -t 30 -c 4 xargs -a analyses/rgrochowski_run_accessions.txt -I {} fastqc --noextract -o analyses/fastqc/ -t 4 data/sra_fastq/{}_1.fastq.gz data/sra_fastq/{}_2.fastq.gz

module load flash2

srun -t 30 -c 2 xargs -a analyses/rgrochowski_run_accessions.txt -I {} flash2 -t 2 -d data/merged_pairs/ -M 151  -o {}.flash data/sra_fastq/{}_1.fastq.gz data/sra_fastq/{}_2.fastq.gz 2>&1 | tee -a analyses/rgrochowski_flash2.log

srun -t 30 -c 2 flash2 -t 2 -d data/merged_pairs/ -M 90  -o ERR6913112.flash data/sra_fastq/ERR6913112_1.fastq.gz data/sra_fastq/ERR6913112_2.fastq.gz 2>&1 | tee -a analyses/rgrochowski_flash2.log

srun -t 30 -c 2 flash2 -d data/merged_pairs/ -o ERR6913112.flash data/sra_fastq/ERR6913112_1.fastq.gz data/sra_fastq/ERR6913112_2.fastq.gz 2>&1 | tee -a analyses/rgrochowski_flash2.log

# this time 88% of pairs were merged which indicated, that the library size wasn't appropriate for paired-end sequencing. 

# around 13% of reads had an overlap lower than the maximum set as the default, it was hence changed with -M parameter
# less data/merged_pairs/ERR6913112.flash.histogram - histogram of the reads reinforces this statement

seqkit stat data/merged_pairs/ERR6913112.flash.extendedFrags.fastq

# setup bowtie2 alignment and download the reference for PhiX phage libraries were spiked with


# setup bowtie2 alignment and download the reference for PhiX phage libraries were spiked

efetch -db nuccore -id NC_001422 -format fasta > data/reference_seqs/PhiX_NC_001422.fna

# create an bowtie2 index database

srun bowtie2-build -f data/reference_seqs/PhiX_NC_001422.fna data/bowtie2_DBs/PhiX_bowtie2_DB

# run the alignment of merged reads with PhiX spike in

srun -t 20 -c 10 bowtie2 -p 10 --no-unal -S analyses/bowtie/rgrochowski_merged2PhiX.sam -x data/bowtie2_DBs/PhiX_bowtie2_DB -U data/merged_pairs/*extendedFrags.fastq  2>&1 | tee -a analyses/bowtie/rgrochowski_bowtie2_PhiX.log
