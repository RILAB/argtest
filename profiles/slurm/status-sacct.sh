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
#
# Keep forks to a minimum: this script runs once per in-flight job on every
# poll, and on a head node with a per-user cgroup PID cap a fork-heavy poll
# (sacct | head | awk + squeue) trips "fork: Resource temporarily
# unavailable". We do the first-line/first-field parse with bash builtins
# (read + parameter expansion) so the only fork is the sacct call itself, and
# we only fork squeue when sacct came back empty.
state=""
while IFS='|' read -r first _; do
    state="$first"
    break
done < <(sacct -j "$jobid" --format=State --noheader --parsable2 2>/dev/null)
state="${state%% *}"

if [ -z "$state" ]; then
    state=$(squeue -h -j "$jobid" -o '%T' 2>/dev/null || true)
fi

case "$state" in
    COMPLETED)
        echo success
        ;;
    # PREEMPTED is "running" here, NOT failed: the `low` partition is
    # PreemptMode=REQUEUE (and the cluster default is JobRequeue=1), so a
    # preempted job is automatically requeued and runs again to completion.
    # If we reported it as failed, Snakemake would give up on a job SLURM is
    # about to finish — and without --keep-going that halts the whole DAG.
    # (On a CANCEL-mode partition this mapping would be wrong; revisit if the
    # submit partition changes.)
    PENDING|RUNNING|COMPLETING|CONFIGURING|REQUEUED|RESIZING|SUSPENDED|PREEMPTED)
        echo running
        ;;
    BOOT_FAIL|CANCELLED|CANCELLED+|DEADLINE|FAILED|NODE_FAIL|OUT_OF_MEMORY|TIMEOUT)
        echo failed
        ;;
    *)
        # Unknown / empty state: keep the executor polling instead of failing
        # the job — sacct sometimes lags behind submission by a few seconds.
        echo running
        ;;
esac
