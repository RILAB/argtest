#!/usr/bin/env python3
"""Warn about proposed ARG individual/ploidy contract violations."""

import argparse
from pathlib import Path

from argtest_common import audit_individual_contract, load_ts


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ts", nargs="+", type=Path, help="Tree sequences to audit")
    parser.add_argument("--name-substring-to-remove", default="")
    return parser.parse_args()


def main():
    args = parse_args()
    for path in args.ts:
        findings = audit_individual_contract(
            load_ts(path),
            name_substring_to_remove=args.name_substring_to_remove,
        )
        if findings:
            for finding in findings:
                print(f"WARNING: {path}: {finding}")
        else:
            print(f"OK: {path}")


if __name__ == "__main__":
    main()
