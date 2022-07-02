#! /bin/bash -l
#SBATCH --job-name=count
#SBATCH --output=logs/counts_%A_%a.log
#SBATCH --time=168:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=8
#SBATCH --partition=panda

samplefile=$1

sample=$(cat  $samplefile | cut -f1 -d' ' |sed -n "${SLURM_ARRAY_TASK_ID}p")

conda activate sambcfenv

samtools sort -@ 8 -o "$sample".sorted $sample
samtools index -@ 8 "$sample".sorte
coverage=$(bcftools mpileup --threads 8 --no-reference "$sample".sorted | awk -v X="3" '$4>=X' | wc -l)
echo -e "$sample\t$coverage" > "$sample"_covinfo
