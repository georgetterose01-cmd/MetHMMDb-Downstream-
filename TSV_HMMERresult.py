#!/usr/bin/env python3

import csv
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

columns = [
    "target_name",
    "target_accession",
    "target_length",
    "query_name",
    "query_accession",
    "query_length",
    "full_Evalue",
    "full_score",
    "full_bias",
    "domain_number",
    "domain_count",
    "conditional_Evalue",
    "independent_Evalue",
    "domain_score",
    "domain_bias",
    "hmm_from",
    "hmm_to",
    "ali_from",
    "ali_to",
    "env_from",
    "env_to",
    "accuracy",
    "description"
]

with open(input_file, "r") as infile, \
     open(output_file, "w", newline="") as outfile:

    writer = csv.writer(outfile, delimiter="\t")
    writer.writerow(columns)

    for line in infile:
        if line.startswith("#") or not line.strip():
            continue

        # HMMER separates the final description using extra spaces.
        fields = line.strip().split(maxsplit=22)

        if len(fields) == 23:
            writer.writerow(fields)
        else:
            print(
                f"Warning: unexpected number of fields ({len(fields)}): "
                f"{line.strip()}",
                file=sys.stderr
            )

print(f"Output written to: {output_file}")
