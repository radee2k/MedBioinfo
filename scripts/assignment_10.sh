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
module load sra-tools

# configure the sra-toolkit
vdb-config --interactive

# go to the project directory
cd /shared/home/rgrochowski/medbioinfo_folder/radoslaw/MedBioinfo/

# download the database
# fastq-dump --gzip --split-3 --readids --disable-multithreading -X 10 --outdir /shared/home/rgrochowski/medbioinfo_folder/radoslaw/MedBioinfo/data/sra_fastq/ ERR6913113
cat rgrochowski_run_accessions.txt | xargs fastq-dump --gzip --split-3 --readids --disable-multithreading -X 10 --outdir /shared/home/rgrochowski/medbioinfo_folder/radoslaw/MedBioinfo/data/sra_fastq/

# count the number of reads per file
cd ~/medbioinfo_folder/radoslaw/MedBioinfo/data/sra_fastq; for file in ./* ; do ( echo $file ; zcat $file | grep ^@ | wc -l ) >> ../../analyses/seq_count.txt ; done

# get back to the project directory
cd /shared/home/rgrochowski/medbioinfo_folder/radoslaw/MedBioinfo/
# use seqkit stats to inspect the reads, -j - threads
module load seqkit
srun -t 30 -c 1 seqkit stats -j 1 ./* >> /analyses/seq_count_seqkit.txt

# the reads seem to be trimmed and quality filtered as nearly all the bases have the same, good quality score and inspected reads had no low quality ends and adapters present and about standard size (150bp)

# check if the reads were deduplicated/the dataset contains unwanted read enrichment
srun -t 30 -c 1 ls | grep gz$ | xargs -I {} seqkit rmdup -j 1 {} -D analyses/seq_count_dup.txt -o dedup/{}

# output
# [INFO] 0 duplicated records removed 16 times and the seq_count_dup.txt is empty

# hence, reads weren't duplicated

# use seqkit locate to look for adaptor sequences
srun -t 30 -c 4 seqkit locate -j 4 ./data/sra_fastq/*.gz -p AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA >> analyses/trimming.txt

# no adaptors were found among the sequences

# more complicated but a nice way of joining stdout 
# srun -c 4 -t 30 ls ./data/sra_fastq/*.gz | xargs -i sh -c ' echo {} ; seqkit locate -j 4 {} -p GAGCACACGTCTGAACTCCAGTCA,GAGCGTCGTGTAGGGAAAGAGTGT'  >> ../../analyses/trimming_locate2.txt

# looking for shorter adaptors
srun -t 30 -c 4 seqkit locate -j 4 ./data/sra_fastq/*.gz -p GAGCACACGTCTGAACTCCAGTCA,GAGCGTCGTGTAGGGAAAGAGTGT >> analyses/trimming2.txt

# after shortening of queries, some adaptor parts could be found in a few reads - reads were at least partially trimmed

# use fastqc for quality control
module load fastqc 

# srun -t 30 -c 4 fastqc --noextract -o analyses/fastqc/ -t 4 data/sra_fastq/ERR6913122_1.fastq.gz data/sra_fastq/ERR6913122_2.fastq.gz
srun -t 30 -c 4 xargs -a analyses/rgrochowski_run_accessions.txt -I {} fastqc --noextract -o analyses/fastqc/ -t 4 data/sra_fastq/{}_1.fastq.gz data/sra_fastq/{}_2.fastq.gz


# use flash2 to merge paired end reads with overlaps
module load flash2
srun -t 30 -c 2 xargs -a analyses/rgrochowski_run_accessions.txt -I {} flash2 -t 2 -d data/merged_pairs/ -M 151  -o {}.flash data/sra_fastq/{}_1.fastq.gz data/sra_fastq/{}_2.fastq.gz 2>&1 | tee -a analyses/rgrochowski_flash2.log

# around 13% of reads had an overlap lower than the maximum set as the default, it was hence changed with -M parameter
srun -t 30 -c 2 flash2 -t 2 -d data/merged_pairs/ -M 90  -o ERR6913112.flash data/sra_fastq/ERR6913112_1.fastq.gz data/sra_fastq/ERR6913112_2.fastq.gz 2>&1 | tee -a analyses/rgrochowski_flash2.log

# most reads were still above the threshold
srun -t 30 -c 2 flash2 -t 2 -d data/merged_pairs/ -M 151  -o ERR6913112.flash data/sra_fastq/ERR6913112_1.fastq.gz data/sra_fastq/ERR6913112_2.fastq.gz 2>&1 | tee -a analyses/rgrochowski_flash2.log

# this time 88% of pairs were merged which indicated, that the library size wasn't appropriate for paired-end sequencing. 

# inspect the merged reads
# less data/merged_pairs/ERR6913112.flash.histogram - histogram of the reads reinforces this statement
# seqkit stat data/merged_pairs/ERR6913112.flash.extendedFrags.fastq

# download the reference for PhiX phage libraries were spiked with and setup bowtie2 alignment 
sh -c "$(wget -q ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
efetch -db nuccore -id NC_001422 -format fasta > data/reference_seqs/PhiX_NC_001422.fna

# create an bowtie2 index database
module bowtie2
srun bowtie2-build -f data/reference_seqs/PhiX_NC_001422.fna data/bowtie2_DBs/PhiX_bowtie2_DB

# run the alignment of merged reads with PhiX spike in
srun -t 20 -c 10 bowtie2 -p 10 --no-unal -S analyses/bowtie/rgrochowski_merged2PhiX.sam -x data/bowtie2_DBs/PhiX_bowtie2_DB -U data/merged_pairs/*extendedFrags.fastq  2>&1 | tee -a analyses/bowtie/rgrochowski_bowtie2_PhiX.log

# No reads were found to align to PhiX phage's genome.

# download the reference for SARS-CoV-2
efetch -db nuccore -id NC_045512 -format fasta > data/reference_seqs/sc2_NC_045512.fna

# index it
srun bowtie2-build -f data/reference_seqs/sc2_NC_045512.fna data/bowtie2_DBs/sc2_bowtie2_DB

# Align patient's reads to the SARS-CoV-2 reference
srun -t 20 -c 10 bowtie2 -p 10 --no-unal -S analyses/bowtie/rgrochowski_merged2sc2.sam -x data/bowtie2_DBs/sc2_bowtie2_DB -U data/merged_pairs/*extendedFrags.fastq  2>&1 | tee -a analyses/bowtie/rgrochowski_bowtie2_sc2.log

# 2284 (0.08%) sequences aligned exactly 1 time indicating presence of SC2 genome in the sample.
samtools view -bh analyses/bowtie/rgrochowski_merged2sc2.sam > analyses/bowtie/rgrochowski_merged2sc2.bam

# samtools view --no-header bowtie/rgrochowski_merged2sc2.sam | cut -d' ' -f 1 | cut -d'.' -f 1 | sort | uniq | wc -l
# output: 5 - the virus was present in 5 samples (DNA or cDNA)

# prepare the file to be viewed with tview
# Upon inspection with tview, the reads are appear well aligned and across all the viral genome.
samtools sort analyses/bowtie/rgrochowski_merged2sc2.bam >| analyses/bowtie/rgrochowski_merged2sc2_sort.bam
samtools tview analyses/bowtie/rgrochowski_merged2sc2_sort_index.bam
# samtools tview analyses/bowtie/rgrochowski_merged2sc2_sort.bam data/reference_seqs/sc2_NC_045512.fna

samtools view  -b analyses/bowtie/rgrochowski_merged2sc2.sam analyses/bowtie/rgrochowski_merged2sc2.bam

# aggregate all the logs and other results with multiqc
module load multiqc
cd analyses/
srun multiqc --force --title "rgrochowski sample sub-set" ../data/merged_pairs/ ./fastqc/ ./rgrochowski_flash2.log ./bowtie/
