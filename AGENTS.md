1. Git repo check
- If the current directory is not a git repository, ask the user whether they want to create one before making changes.

2. Data-file permission (always ask for modifications)
- Always request permission before modifying an existing data file.
- Permission is not required for creating a new data file.
- User instructions such as “modify,” “remove,” “edit,” or similar do **not** count as permission for data-file modifications. You must still ask.

3. Data-file definition
- A data file is any non-code project file used as biological input, output, or reference data.
- Examples include: `*.txt`, `*.csv`, `*.tsv`, `*.fastq`, `*.fastq.gz`, `*.sam`, `*.bam`, `*.maf`, `*.map`, `*.bed`, `*.vcf`, `*.gvcf`, `*.fai`, and similar biological data or tabular reference files.
- If unsure whether a file is a data file, treat it as a data file and ask permission first.

4. Uncommitted changes check (target file only)
- Before modifying a file, check whether the **target file** has uncommitted changes.
- If the target file has uncommitted changes, ask whether the user wants to commit those changes first.
- This rule is target-file-only (not repo-wide).

5. One permission can cover multiple files
- A single permission request is sufficient if the user clearly authorizes modifying multiple specific data files in the same task.

6. Backups for data files
- When modifying a data file (with permission), create a backup copy first using the `.bak` suffix (for example, `file.txt.bak`).

7. Symlink write policy
- If a path to be modified is a symbolic link, never modify the symlink target.
- If modification is required (with permission), create a regular-file copy in the current working directory, modify that copy, and replace the symlink path with the modified regular file.
- Prefer `path.tmp` + atomic rename (`mv path.tmp path`) so the symlink is replaced by a regular file.
- Do not use in-place editors (`sed -i`, `perl -pi`, etc.) on symlink paths.
- If a backup is required, back up the symlink path as it exists before replacement.

8. Write scope restriction
- Do not create, modify, or delete files outside the current working directory unless the user explicitly requests it.
- Exception: temporary files may be created or modified in `/tmp` and system temporary directories (for example, macOS `/var/folders/...`) when needed for task execution.
- Files written to `/tmp` should be treated as temporary working files, not final outputs, unless the user explicitly requests otherwise.
- This rule does not permit modifying symlink targets outside the current working directory; symlink paths must follow the symlink write policy.
- Tool-generated temporary files, caches, and logs are allowed in the current working directory or approved temporary directories when required to complete the task.


9. Data-file permission exception (temporary paths)
- Permission is not required for creating or modifying data files under `/tmp` or system temporary directories (for example macOS `/var/folders/...`) when they are temporary working files used to run, test, or validate the project.
- This exception does not apply to files in the repository working directory (including `tests/`, `results/`, `example_data/`, `config.yaml`, `README.md`, `*.md`, etc.) or any other non-temporary location.
- Temporary files created under `/tmp` remain subject to the symlink write policy if the target path is a symlink.
- The agent should prefer `/tmp/<project>-<purpose>/...` paths for temporary data outputs to make scope explicit.

10. Bootstrap AGENTS.md into repos

If AGENTS.md is missing in the current git repository root, copy the default AGENTS.md into that repository before doing other work.
Do not overwrite an existing AGENTS.md.
If the current directory is not a git repository, follow the existing git-repo check rule first.

11. Conda environment bootstrap
- Before running any Python, pytest, pip, or project CLI command, first initialize conda in the shell with `module load conda`.
- If the current working directory contains `environment.yml`, detect the environment name from its `name:` field and run `conda activate <that-name>` before continuing.
- If both `environment.yml` and an already-active conda environment are present, prefer the environment named in `environment.yml`.
- Run Python-related commands in a login bash shell so `module` and `conda activate` work correctly.
- If activation fails, stop and report the error before running project commands in another Python environment.

12. Offload heavy commands to SLURM (this is a head node)
- Do not run heavy work directly on the head node. Any command expected to use more than 1 CPU, or substantially more than ~1 GB of RAM, must be launched via the helper `~/.claude/bin/hpc_run 'COMMAND'`.
- `hpc_run` submits with `srun`, preferring partition `low` and falling back to `high` if `low` cannot allocate a node within 60s. Account defaults to `jrigrp`.
- Set resources per task with env vars, estimating generously: `HPC_CPUS`, `HPC_MEM` (for example `64G`), `HPC_TIME` (for example `04:00:00`). Defaults are 1 CPU / 8G / 02:00:00.
- The helper preserves cwd, points caches at node-local scratch, and activates the conda env from `./environment.yml` automatically, so do not also wrap the command in `module load conda` or `conda activate`; that would double-activate.
- Keep genuinely light commands on the head node: git, ls, grep, file reads, `squeue`/`sinfo`/`sacct`, `--version` checks, and small/quick scripts.
- When unsure whether a command is heavy, prefer `hpc_run`.
- Example: `HPC_MEM=64G ~/.claude/bin/hpc_run 'python measure_merge_mem.py 0'`

13. Standard options prompt for data-file decisions
- When a task requires user confirmation about committing target-file changes and/or modifying existing data files, do not ask an open-ended question.
- Use this exact options prompt and ask the user to choose one option:
  1) commit, backup, and modify
  2) don't commit, but backup and modify
  3) do not modify
  4) other
- Apply the user's choice to the relevant target files for that task.

--- project-doc ---

## Repo Rules

- Do not ask whether to create a git commit or push after making changes.
- When the user says "commit and push", immediately run the necessary git commands (including `git add`, `git status`, `git commit`, `git push`) without any further confirmation or questions.
- Never ask whether to run `git commit` or `git push` after the user has already said "commit and push".
- For Snakemake runs, use a longer command timeout to avoid premature termination.
