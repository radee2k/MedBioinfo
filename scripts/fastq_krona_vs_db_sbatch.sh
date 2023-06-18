#!/bin/bash
#
#SBATCH --partition=fast             # long, fast, etc.
#SBATCH --ntasks=1                   # nb of *tasks* to be run in // (usually 1), this task can be multithreaded (see cpus-per-task)
#SBATCH --nodes=1                    # nb of nodes to reserve for each task (usually 1)
#SBATCH --cpus-per-task=16            # nb of cpu (in fact cores) to reserve for each task /!\ job killed if commands below use more cores
#SBATCH --mem=80GB                  # amount of RAM to reserve for the tasks /!\ job killed if commands below use more RAM
#SBATCH --time=00:15:00               # maximal wall clock duration (D-HH:MM) /!\ job killed if commands below take more time than reservation
#SBATCH -o /shared/ifbstor1/projects/2314_medbioinfo/radoslaw/MedBioinfo/analyses/kraken2/all/outputs_kraken/slurm.%A.%a.out   # standard output (STDOUT) redirected to these files (with Job ID and array ID in file names) HAS TO EXIST BEFORE BEING SUBMITTED!!!!!!!!!
#SBATCH -e /shared/ifbstor1/projects/2314_medbioinfo/radoslaw/MedBioinfo/analyses/kraken2/all/outputs_kraken/slurm.%A.%a.err   # standard error  (STDERR) redirected to these files (with Job ID and array ID in file names)
# /!\ Note that the ./outputs/ dir above needs to exist in the dir where script is submitted **prior** to submitting this script
#SBATCH --array=1-8                # 1-N: clone this script in an array of N tasks: $SLURM_ARRAY_TASK_ID will take the value of 1,2,...,N
## SBATCH --job-name=rg_kraken2        # name of the task as displayed in squeue & sacc, also encouraged as srun optional parameter
#SBATCH --mail-type END              # when to send an email notiification (END = when the whole sbatch array is finished)
#SBATCH --mail-user radoslaw.grochowski@ki.se

#################################################################
# Preparing work (cd to working dir, get hold of input data, convert/un-compress input data when needed etc.)
workdir="/shared/ifbstor1/projects/2314_medbioinfo/radoslaw/MedBioinfo/"
#datadir="/shared/ifbstor1/projects/2314_medbioinfo/radoslaw/MedBioinfo/data/merged_pairs"
input_path="/shared/ifbstor1/projects/2314_medbioinfo/radoslaw/MedBioinfo/data/sra_fastq/"
accnum_file="/shared/ifbstor1/projects/2314_medbioinfo/radoslaw/MedBioinfo/analyses/rgrochowski_run_accessions.txt"
db="/shared/projects/2314_medbioinfo/kraken2/arch_bact_vir_hum_protoz_fung/"


echo START: `date`

module load kraken2 bracken #as required

cd ${workdir}


accnum=$( sed -n  "$SLURM_ARRAY_TASK_ID"p ${accnum_file} )
output_file=${workdir}"analyses/kraken2/all/${accnum}_kraken2"
input_file=${input_path}${accnum}
#input_file1=${input_file}_1.fastq.gzl
#input_file2=${input_file}_2.fastq.gz

echo ${accnum}
echo ${input_file}
echo ${output_file}

# Start work
srun  --job-name=${accnum} kraken2 --paired --threads ${SLURM_CPUS_PER_TASK} --db ${db}  --gzip-compressed --output ${output_file}.out --report ${output_file}_report.tsv  ${input_file}_1.fastq.gz ${input_file}_2.fastq.gz

#srun kraken2 --paired --threads ${SLURM_CPUS_PER_TASK} --db ${db}  --gzip-compressed --output ${output_file}.out --report ${output_file}_report.tsv  ${input_file1} ${input_file2}

srun --job-name=${accnum}_rg_bracken bracken -d ${db} -i ${output_file}_report.tsv -o ${output_file}_bracken.tsv -w  ${output_file}_bracken_report.tsv


# run krona on all the reports
# cd /shared/ifbstor1/projects/2314_medbioinfo/radoslaw/MedBioinfo/analyses/krona

output_file2="${workdir}analyses/krona/${accnum}_kraken2"

srun --job-name=${accnum}_rg_krona1 /shared/projects/2314_medbioinfo/kraken2/KrakenTools/kreport2krona.py -r ${output_file}_bracken_report.tsv -o ${output_file2}_krona.out

sed -Ei "s/.{1}__//g" ${output_file2}_krona.out # -E for easy Perl 

srun --job-name=${accnum}_rg_krona2 /shared/projects/2314_medbioinfo/kraken2/bin/ktImportText -o ${output_file2}_krona.html ${output_file2}_krona.out



echo END: `date`