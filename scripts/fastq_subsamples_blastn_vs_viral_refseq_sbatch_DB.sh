#!/bin/bash
#
#SBATCH --partition=fast             # long, fast, etc.
#SBATCH --ntasks=1                   # nb of *tasks* to be run in // (usually 1), this task can be multithreaded (see cpus-per-task)
#SBATCH --nodes=1                    # nb of nodes to reserve for each task (usually 1)
#SBATCH --cpus-per-task=1            # nb of cpu (in fact cores) to reserve for each task /!\ job killed if commands below use more cores
#SBATCH --mem=70GB                  # amount of RAM to reserve for the tasks /!\ job killed if commands below use more RAM
#SBATCH --time=02:00:00               # maximal wall clock duration (D-HH:MM) /!\ job killed if commands below take more time than reservation
#SBATCH -o ./outputs_array_db/slurm.%A.%a.out   # standard output (STDOUT) redirected to these files (with Job ID and array ID in file names) HAS TO EXIST BEFORE BEING SUBMITTED!!!!!!!!!
#SBATCH -e ./outputs_array_db/slurm.%A.%a.err   # standard error  (STDERR) redirected to these files (with Job ID and array ID in file names)
# /!\ Note that the ./outputs/ dir above needs to exist in the dir where script is submitted **prior** to submitting this script
#SBATCH --array=1-8                # 1-N: clone this script in an array of N tasks: $SLURM_ARRAY_TASK_ID will take the value of 1,2,...,N
#SBATCH --job-name=MedBioinfo        # name of the task as displayed in squeue & sacc, also encouraged as srun optional parameter
#SBATCH --mail-type END              # when to send an email notiification (END = when the whole sbatch array is finished)
#SBATCH --mail-user radoslaw.grochowski@ki.se

#################################################################
# Preparing work (cd to working dir, get hold of input data, convert/un-compress input data when needed etc.)
workdir="/shared/ifbstor1/projects/2314_medbioinfo/radoslaw/MedBioinfo/analyses/blastn/array_db"
datadir="/shared/ifbstor1/projects/2314_medbioinfo/radoslaw/MedBioinfo/data/merged_pairs"
input_file1="/shared/ifbstor1/projects/2314_medbioinfo/radoslaw/MedBioinfo/data/merged_pairs/samples_for_blast/ERR6913112_sample_100.fna"
accnum_file="/shared/ifbstor1/projects/2314_medbioinfo/radoslaw/MedBioinfo/analyses/rgrochowski_run_accessions.txt"
output_file=${workdir}"analyses/blastn/db/fastq_subsample_blastn_vs_NT"
db="/shared/bank/nt/current/blast/nt"


echo START: `date`

module load seqkit blast #as required

mkdir -p ${workdir}      # -p because it creates all required dir levels **and** doesn't throw an error if the dir exists :)
cd ${workdir}

echo ${accnum_file}

accnum=$( sed -n  "$SLURM_ARRAY_TASK_ID"p ${accnum_file} )
input_file="${datadir}/${accnum}.flash.extendedFrags"

echo ${accnum}
echo ${input_file}


seqkit range -r 1:10 ${input_file}.fasta > ${input_file}_10.fasta


output_file="${workdir}/${accnum}_blastn_vs_viral.out"

# Start work
srun --job-name=${accnum} blastn -num_threads ${SLURM_CPUS_PER_TASK} -query ${input_file}_10.fasta -db ${db} -evalue 1E-10 -outfmt 6 -perc_identity 75 -out ${output_file}_10.out -max_target_seqs 5


echo END: `date`




# this extracts the item number $SLURM_ARRAY_TASK_ID from the file of accnums
# accnum=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${accnum_file})
# input_file="${datadir}/${accnum}.fastq"
# alternatively, just extract the input file as the item number $SLURM_ARRAY_TASK_ID in the data dir listing
# this alternative is less handy since we don't get hold of the isolated "accnum", which is very handy to name the srun step below :)
# input_file=$(ls "${datadir}/*.fastq.gz" | sed -n ${SLURM_ARRAY_TASK_ID}p)

# if the command below can't cope with compressed input
# srun gunzip "${input_file}.gz"

# because there are mutliple jobs running in // each output file needs to be made unique by post-fixing with $SLURM_ARRAY_TASK_ID and/or $accnum
# output_file="${workdir}/ABCjob.${SLURM_ARRAY_TASK_ID}.${accnum}.out"

#################################################################
# Start work
# srun --job-name=${accnum} some_abc_software --threads ${SLURM_CPUS_PER_TASK} --in ${input_file} --out ${output_file}

#################################################################
# Clean up (eg delete temp files, compress output, recompress input etc)
# srun gzip ${input_file}
#srun gzip ${output_file}

