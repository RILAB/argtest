#!/bin/bash
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --job-name=smk-controller
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=14-00:00:00
#SBATCH --output=logs/slurm/controller-%j.out
# NOTE: --error is intentionally unset so stderr merges into the --output file.
#
# Snakemake CONTROLLER driver job.
# ---------------------------------
# Runs the long-lived snakemake controller on a COMPUTE node instead of the
# head node. Its per-job polling threads and status-script forks then count
# against this allocation's cgroup, not the head node's per-user 512-task cap
# (user.slice pids.max) that throttled us before ("fork: Resource temporarily
# unavailable").
#
# Partition is `high` ON PURPOSE: it is non-preemptable (150-day max walltime).
# The controller must NOT run on `low` -- a preempted controller would take the
# whole run down with it. The per-rule jobs it submits still go to whatever
# partition the configfile/Snakefile requests (typically `low`).
#
# Resources: the controller is poll/IO-bound, not compute-heavy. 4 CPUs + 8 GB
# comfortably hosts the ~jobs(=250) polling threads. Walltime is generous
# because step5 alone runs ~1-2 days per replicate.
#
# Usage (from the repo root):
#   mkdir -p logs/slurm
#   sbatch profiles/slurm/run-controller.sbatch [CONFIGFILE] [EXTRA snakemake args...]
# Defaults to config/snakemake.yaml. --rerun-incomplete is always applied.
# Examples:
#   sbatch profiles/slurm/run-controller.sbatch
#   sbatch profiles/slurm/run-controller.sbatch config/snakemake.yaml -n   # dry-run
#
# Watch progress:  tail -f logs/slurm/controller-<jobid>.out

set -euo pipefail

# First positional arg is the configfile, but only if it isn't itself a flag
# (so `sbatch ... -n` still works without a configfile). Everything left over
# is forwarded to snakemake verbatim.
configfile="config/snakemake.yaml"
if [[ "${1:-}" && "${1:-}" != -* ]]; then
    configfile="$1"
    shift
fi

# Activate the conda env named in environment.yml (detect, don't hardcode --
# the dataset/env can vary across configs).
source /etc/profile.d/modules.sh 2>/dev/null || true
module load conda
env_name="$(awk -F': *' '/^name:/{print $2; exit}' environment.yml)"
conda activate "$env_name"

echo "=== snakemake controller ==="
echo "host:       $(hostname)"
echo "slurm job:  ${SLURM_JOB_ID:-<none>}"
echo "configfile: $configfile"
echo "conda env:  $env_name"
echo "extra args: $*"
echo "============================"

exec snakemake \
    --profile profiles/slurm \
    --configfile "$configfile" \
    --rerun-incomplete \
    "$@"
