
#!/usr/bin/env python3
# Analyze HMMER results obtained with MetHMMDB.
# Uses only the Python 3 standard library.

import csv
import sys
from pathlib import Path
from collections import defaultdict

# ---------- Select input and output paths ----------
input_file = Path(
    input("Path to input TSV file: ").strip().strip('"').strip("'")
).expanduser()

output_dir = Path(
    input("Path to output directory: ").strip().strip('"').strip("'")
).expanduser()

if not input_file.is_file():
    sys.exit(f"ERROR: Input file not found: {input_file}")

output_dir.mkdir(parents=True, exist_ok=True)

# ---------- Read the TSV file ----------
with input_file.open("r", encoding="utf-8-sig", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    if not reader.fieldnames:
        sys.exit("ERROR: The input file has no header.")

    reader.fieldnames = [name.strip() for name in reader.fieldnames]
    required = [
        "target_name", "query_name", "full_Evalue",
        "full_score", "domain_score", "query_length",
        "hmm_from", "hmm_to", "description"
    ]

    missing = [col for col in required if col not in reader.fieldnames]
    if missing:
        sys.exit("ERROR: Missing columns: " + ", ".join(missing))

    rows = list(reader)

def number(row, column):
    """Convert a field to float; return None if unavailable."""
    try:
        return float(row[column])
    except (ValueError, TypeError, KeyError):
        return None

def write_tsv(filename, fieldnames, data):
    """Write a tab-separated output file."""
    with (output_dir / filename).open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames,
            delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(data)

# ---------- Calculate HMM profile coverage ----------
for row in rows:
    qlen = number(row, "query_length")
    start = number(row, "hmm_from")
    end = number(row, "hmm_to")

    if qlen and start is not None and end is not None and qlen > 0:
        row["profile_coverage_percent"] = round(
            (end - start + 1) / qlen * 100, 2
        )
    else:
        row["profile_coverage_percent"] = ""

# ---------- Apply the published full-sequence E-value cutoff ----------
# The paper reports HMMER E-value <= 1e-10 for MetHMMDB.
cutoff = 1e-10

filtered = [
    row for row in rows
    if number(row, "full_Evalue") is not None
    and number(row, "full_Evalue") <= cutoff
]

# Preserve all qualifying rows, including repeated domain hits.
write_tsv(
    "01_Coincidences.tsv",
    list(rows[0].keys()) if rows else required,
    filtered
)

# ---------- Unique protein-profile pairs ----------
# If a pair occurs more than once, retain the row with the best full score.
best_pair = {}

for row in filtered:
    key = (row["target_name"], row["query_name"])
    previous = best_pair.get(key)

    score = number(row, "full_score")
    previous_score = number(previous, "full_score") if previous else None

    if previous is None or (
        score is not None and
        (previous_score is None or score > previous_score)
    ):
        best_pair[key] = row

pairs = list(best_pair.values())

write_tsv(
    "02_protein_profiles.tsv",
    list(rows[0].keys()) + (
        ["profile_coverage_percent"]
        if rows and "profile_coverage_percent" not in rows[0]
        else []
    ) if rows else required + ["profile_coverage_percent"],
    pairs
)

# ---------- Summarize by unique protein ----------
by_protein = defaultdict(list)
for row in pairs:
    by_protein[row["target_name"]].append(row)

protein_summary = []

for protein, hits in sorted(by_protein.items()):
    evals = [number(r, "full_Evalue") for r in hits]
    scores = [number(r, "full_score") for r in hits]
    coverages = [
        number(r, "profile_coverage_percent") for r in hits
    ]

    protein_summary.append({
        "target_name": protein,
        "number_of_profiles": len(set(r["query_name"] for r in hits)),
        "profiles": "; ".join(sorted(set(r["query_name"] for r in hits))),
        "annotations": "; ".join(sorted(set(
            r["description"] for r in hits if r.get("description")
        ))),
        "best_full_Evalue": min(x for x in evals if x is not None)
            if any(x is not None for x in evals) else "",
        "best_full_score": max(x for x in scores if x is not None)
            if any(x is not None for x in scores) else "",
        "maximum_profile_coverage": max(
            (x for x in coverages if x is not None), default=""
        )
    })

write_tsv(
    "03_unique_proteins.tsv",
    [
        "target_name", "number_of_profiles", "profiles",
        "annotations", "best_full_Evalue", "best_full_score",
        "maximum_profile_coverage"
    ],
    protein_summary
)

# ---------- Summarize by HMM profile ----------
by_profile = defaultdict(list)
for row in pairs:
    by_profile[row["query_name"]].append(row)

profile_summary = []

for profile, hits in by_profile.items():
    evals = [number(r, "full_Evalue") for r in hits]
    scores = [number(r, "full_score") for r in hits]
    coverages = [
        number(r, "profile_coverage_percent") for r in hits
    ]

    profile_summary.append({
        "query_name": profile,
        "unique_proteins": len(set(r["target_name"] for r in hits)),
        "best_Evalue": min(x for x in evals if x is not None)
            if any(x is not None for x in evals) else "",
        "best_score": max(x for x in scores if x is not None)
            if any(x is not None for x in scores) else "",
        "mean_profile_coverage": round(
            sum(x for x in coverages if x is not None) /
            sum(x is not None for x in coverages), 2
        ) if any(x is not None for x in coverages) else ""
    })

profile_summary.sort(
    key=lambda r: r["unique_proteins"], reverse=True
)

write_tsv(
    "04_profiles.tsv",
    [
        "query_name", "unique_proteins", "best_Evalue",
        "best_score", "mean_profile_coverage"
    ],
    profile_summary
)

# ---------- Summary by annotation ----------
by_annotation = defaultdict(list)
for row in pairs:
    by_annotation[row.get("description", "")].append(row)

annotation_summary = []

for annotation, hits in by_annotation.items():
    annotation_summary.append({
        "description": annotation,
        "unique_proteins": len(set(r["target_name"] for r in hits)),
        "profiles": "; ".join(sorted(set(r["query_name"] for r in hits)))
    })

annotation_summary.sort(
    key=lambda r: r["unique_proteins"], reverse=True
)

write_tsv(
    "05_annotation.tsv",
    ["description", "unique_proteins", "profiles"],
    annotation_summary
)

# ---------- Report ----------
print("\nAnalysis completed.")
print(f"Input rows: {len(rows)}")
print(f"Rows passing full-sequence E-value <= {cutoff}: {len(filtered)}")
print(f"Unique protein-profile pairs: {len(pairs)}")
print(f"Unique proteins: {len(by_protein)}")
print(f"Unique HMM profiles: {len(by_profile)}")
print(f"Output directory: {output_dir.resolve()}")
