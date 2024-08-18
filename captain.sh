#!/usr/bin/env bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=32g
#SBATCH --time=10-00:00:00
#SBATCH --parsable
#SBATCH -J "[[PIPENICKNAME]]"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output "[[WORKDIR]]/Pipe_runtime/[[SNAPSHOT]]/logs/Pipe_runtime.%j.out"
#SBATCH --error  "[[WORKDIR]]/Pipe_runtime/[[SNAPSHOT]]/logs/Pipe_runtime.%j.err"
set -euo pipefail

current_hostname=$(hostname)

# Check hostname to avoid running snakemake in the biowulf login node
if [ "$current_hostname" == "biowulf.nih.gov" ]; then
    echo "Current hostname is biowulf.nih.gov. Stopping the script."
    exit 1
else
    module load snakemake singularity
fi

uid=$(uuidgen|cut -d '-' -f1)
echo 'Snakemake captain uid:' $uid
echo $uid > [[WORKDIR]]/Pipe_runtime/current_uid 
echo [[SNAPSHOT]] > [[WORKDIR]]/Pipe_runtime/current_snapshot

# Main process of pipeline
snakemake --snakefile Snakefile -d "[[WORKDIR]]" \
  --cores 8 \
  --local-cores 8 \
  --jobs 2000 \
  --latency-wait 120 all \
  --max-jobs-per-second 1 \
  --max-status-checks-per-second 0.01 \
  --use-singularity \
  --singularity-args "-B [[BINDPATH]]" \
  --configfile="[[WORKDIR]]/Pipe_runtime/[[SNAPSHOT]]/config.yaml" \
  --printshellcmds \
  --cluster-config "[[WORKDIR]]/Pipe_runtime/[[SNAPSHOT]]/cluster.yaml" \
  --cluster "sbatch --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} {cluster.extra} --out [[WORKDIR]]/Pipe_runtime/[[SNAPSHOT]]/logs/slurm/{rule}.%j.out" \
  --jobname "{name}.{jobid}.uid_"$uid".sh" \
  --stats "[[WORKDIR]]/Pipe_runtime/[[SNAPSHOT]]/logs/runtime_statistics.json" \
  --keep-going \
  --rerun-incomplete \
  --keep-remote 2>&1

# Create summary report
snakemake -d "[[WORKDIR]]/Pipe_runtime/[[SNAPSHOT]]" --report "[[WORKDIR]]/Pipe_runtime/[[SNAPSHOT]]/Snakemake_Report.html"

