#!/bin/bash
# Job-status command for the snakemake cluster-generic executor.
#
# The executor invokes this with one argument — the external job ID returned
# by sbatch --parsable — and expects exactly one of these tokens on stdout:
#   success   -> job done, output files are ready
#   running   -> job still in flight; keep polling
#   failed    -> job will never finish; advance the DAG with a failure
#
# Without this, the executor falls back to "wait for output file to appear"
# mode and hangs forever when SLURM kills a job before it writes output
# (e.g. TIMEOUT, OUT_OF_MEMORY, NODE_FAIL).

set -euo pipefail

jobid="$1"

# Prefer sacct: it sees both running and terminal states for jobs in the
# accounting DB. Fall back to squeue while sacct is briefly empty for a
# just-submitted job.
state=$(sacct -j "$jobid" --format=State --noheader --parsable2 2>/dev/null \
        | head -n 1 | awk '{print $1}')
if [ -z "$state" ]; then
    state=$(squeue -h -j "$jobid" -o '%T' 2>/dev/null || true)
fi

case "$state" in
    COMPLETED)
        echo success
        ;;
    PENDING|RUNNING|COMPLETING|CONFIGURING|REQUEUED|RESIZING|SUSPENDED)
        echo running
        ;;
    BOOT_FAIL|CANCELLED|CANCELLED+|DEADLINE|FAILED|NODE_FAIL|OUT_OF_MEMORY|PREEMPTED|TIMEOUT)
        echo failed
        ;;
    *)
        # Unknown / empty state: keep the executor polling instead of failing
        # the job — sacct sometimes lags behind submission by a few seconds.
        echo running
        ;;
esac
