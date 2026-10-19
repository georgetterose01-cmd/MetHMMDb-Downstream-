
#!/usr/bin/env python3
# Search a bacterial proteome against MetHMMDB FASTA references.
# Uses Python standard library and local NCBI BLAST+ executables.

import csv
import subprocess
from pathlib import Path
from collections import defaultdict

# Initial exploratory thresholds; adjust after inspecting results.
MAX_EVALUE = 1e-10
MIN_IDENTITY = 40.0
MIN_QUERY_COVERAGE = 40.0
MIN_REFERENCE_COVERAGE = 40.0
THREADS = 4

FASTA_EXTENSIONS = {
    ".fa", ".fas", ".fasta", ".fsa", ".faa", ".pep"
}


def read_fasta(path):
    """Yield (header, sequence) records from a FASTA file."""
    header = None
    sequence = []

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header = line[1:].strip()
                sequence = []
            else:
                sequence.append(line)

    if header is not None:
        yield header, "".join(sequence)


def write_tsv(path, rows, fields):
    """Write dictionaries to a tab-separated file."""
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def interval_coverage(intervals, sequence_length):
    """Calculate coverage from the union of alignment intervals."""
    if sequence_length <= 0 or not intervals:
        return 0.0

    intervals = sorted(
        (min(a, b), max(a, b)) for a, b in intervals
    )
    start, end = intervals[0]
    covered = 0

    for next_start, next_end in intervals[1:]:
        if next_start <= end + 1:
            end = max(end, next_end)
        else:
            covered += end - start + 1
            start, end = next_start, next_end

    covered += end - start + 1
    return min(100.0, covered / sequence_length * 100.0)


def main():
    proteome_text = input("Path to genome protein FASTA (.faa): ").strip()
    refs_text = input("Folder containing the 254 reference FASTA files: ").strip()
    output_text = input("Output folder: ").strip()

    proteome = Path(proteome_text).expanduser().resolve()
    refs_dir = Path(refs_text).expanduser().resolve()
    outdir = Path(output_text).expanduser().resolve()

    if not proteome.is_file():
        raise SystemExit(f"ERROR: Proteome not found: {proteome}")

    if not refs_dir.is_dir():
        raise SystemExit(f"ERROR: Reference folder not found: {refs_dir}")

    if not proteome.suffix.lower() in FASTA_EXTENSIONS:
        print("WARNING: The proteome extension is unusual; proceeding anyway.")

    for executable in ("blastp", "makeblastdb"):
        try:
            subprocess.run(
                [executable, "-version"],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )
        except (FileNotFoundError, subprocess.CalledProcessError):
            raise SystemExit(
                f"ERROR: {executable} is unavailable. "
                "Install BLAST+ with: sudo apt install ncbi-blast+"
            )

    outdir.mkdir(parents=True, exist_ok=True)

    ref_files = sorted(
        p for p in refs_dir.iterdir()
        if p.is_file() and p.suffix.lower() in FASTA_EXTENSIONS
    )

    if not ref_files:
        raise SystemExit("ERROR: No FASTA files found in the reference folder.")

    catalog = []
    combined = outdir / "MetHMMDB_references_combined.faa"
    ref_number = 0

    # Create unique reference IDs while retaining their original provenance.
    with combined.open("w", encoding="utf-8") as output:
        for ref_file in ref_files:
            for header, sequence in read_fasta(ref_file):
                if not sequence:
                    continue

                ref_number += 1
                ref_id = f"REF{ref_number:08d}"

                output.write(f">{ref_id}\n")
                for i in range(0, len(sequence), 80):
                    output.write(sequence[i:i + 80] + "\n")

                catalog.append({
                    "reference_id": ref_id,
                    "source_file": ref_file.name,
                    "original_header": header,
                    "reference_length": len(sequence),
                    "file_stem": ref_file.stem
                })

    if ref_number == 0:
        raise SystemExit("ERROR: No protein sequences found in the FASTA files.")

    catalog_by_id = {r["reference_id"]: r for r in catalog}

    write_tsv(
        outdir / "00_reference_catalog.tsv",
        catalog,
        ["reference_id", "source_file", "original_header",
         "reference_length", "file_stem"]
    )

    print(f"Reference files: {len(ref_files)}")
    print(f"Reference protein sequences: {ref_number}")
    print("Building local BLAST protein database...")

    db_prefix = outdir / "MetHMMDB_protein_db"

    subprocess.run([
        "makeblastdb",
        "-in", str(combined),
        "-dbtype", "prot",
        "-parse_seqids",
        "-out", str(db_prefix)
    ], check=True)

    raw_file = outdir / "01_blastp_alignments.tsv"

    # Query = genome proteins; subject = MetHMMDB reference proteins.
    fields = (
        "qseqid sseqid pident length qlen slen "
        "qstart qend sstart send evalue bitscore"
    )

    print("Running BLASTp. This may take some time...")

    with raw_file.open("w", encoding="utf-8") as raw:
        subprocess.run([
            "blastp",
            "-query", str(proteome),
            "-db", str(db_prefix),
            "-evalue", str(MAX_EVALUE),
            "-num_threads", str(THREADS),
            "-max_target_seqs", "10000",
            "-outfmt", f"6 {fields}"
        ], check=True, stdout=raw)

    # Group HSP alignments by genome protein and reference protein.
    pairs = defaultdict(list)

    with raw_file.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 12:
                continue

            (qid, sid, pid, length, qlen, slen,
             qs, qe, ss, se, ev, bits) = parts

            pairs[(qid, sid)].append({
                "pident": float(pid),
                "length": int(length),
                "qlen": int(qlen),
                "slen": int(slen),
                "qstart": int(qs),
                "qend": int(qe),
                "sstart": int(ss),
                "send": int(se),
                "evalue": float(ev),
                "bitscore": float(bits)
            })

    pair_rows = []

    for (qid, sid), hsps in pairs.items():
        ref = catalog_by_id.get(sid)
        if ref is None:
            continue

        aligned_length = sum(h["length"] for h in hsps)
        weighted_identity = (
            sum(h["pident"] * h["length"] for h in hsps)
            / aligned_length
            if aligned_length else 0.0
        )

        query_cov = interval_coverage(
            [(h["qstart"], h["qend"]) for h in hsps],
            hsps[0]["qlen"]
        )
        reference_cov = interval_coverage(
            [(h["sstart"], h["send"]) for h in hsps],
            hsps[0]["slen"]
        )

        best_evalue = min(h["evalue"] for h in hsps)
        best_bitscore = max(h["bitscore"] for h in hsps)

        pair_rows.append({
            "genome_protein": qid,
            "reference_id": sid,
            "source_file": ref["source_file"],
            "reference_header": ref["original_header"],
            "reference_family_file": ref["file_stem"],
            "identity_percent": round(weighted_identity, 3),
            "alignment_length_sum": aligned_length,
            "genome_protein_length": hsps[0]["qlen"],
            "reference_length": hsps[0]["slen"],
            "query_coverage_percent": round(query_cov, 3),
            "reference_coverage_percent": round(reference_cov, 3),
            "evalue": best_evalue,
            "bitscore_best_HSP": best_bitscore,
            "number_of_HSPs": len(hsps)
        })

    pair_fields = [
        "genome_protein", "reference_id", "source_file",
        "reference_header", "reference_family_file",
        "identity_percent", "alignment_length_sum",
        "genome_protein_length", "reference_length",
        "query_coverage_percent", "reference_coverage_percent",
        "evalue", "bitscore_best_HSP", "number_of_HSPs"
    ]

    pair_rows.sort(key=lambda r: (r["evalue"], -r["bitscore_best_HSP"]))

    write_tsv(
        outdir / "02_all_protein_reference_pairs.tsv",
        pair_rows, pair_fields
    )

    # Filter candidate protein-reference pairs.
    filtered = [
        r for r in pair_rows
        if r["evalue"] <= MAX_EVALUE
        and r["identity_percent"] >= MIN_IDENTITY
        and r["query_coverage_percent"] >= MIN_QUERY_COVERAGE
        and r["reference_coverage_percent"] >= MIN_REFERENCE_COVERAGE
    ]

    write_tsv(
        outdir / "03_filtered_candidate_pairs.tsv",
        filtered, pair_fields
    )

    # One best reference assignment per genome protein.
    # All alternative matches remain available in the previous table.
    best_by_protein = {}

    for row in filtered:
        qid = row["genome_protein"]
        if qid not in best_by_protein:
            best_by_protein[qid] = row

    best_rows = list(best_by_protein.values())
    best_rows.sort(key=lambda r: r["genome_protein"])

    write_tsv(
        outdir / "04_best_reference_per_protein.tsv",
        best_rows, pair_fields
    )

    # Unique candidate proteins, regardless of how many references match.
    unique_rows = []
    for qid, row in sorted(best_by_protein.items()):
        unique_rows.append({
            "genome_protein": qid,
            "genome_protein_length": row["genome_protein_length"],
            "best_reference_id": row["reference_id"],
            "best_reference_file": row["source_file"],
            "best_reference_header": row["reference_header"],
            "identity_percent": row["identity_percent"],
            "query_coverage_percent": row["query_coverage_percent"],
            "reference_coverage_percent": row["reference_coverage_percent"],
            "evalue": row["evalue"],
            "bitscore_best_HSP": row["bitscore_best_HSP"]
        })

    unique_fields = [
        "genome_protein", "genome_protein_length",
        "best_reference_id", "best_reference_file",
        "best_reference_header", "identity_percent",
        "query_coverage_percent", "reference_coverage_percent",
        "evalue", "bitscore_best_HSP"
    ]

    write_tsv(
        outdir / "05_unique_candidate_proteins.tsv",
        unique_rows, unique_fields
    )

    # Summary by original reference FASTA file.
    summary = []

    for source in sorted({r["source_file"] for r in filtered}):
        rows = [r for r in filtered if r["source_file"] == source]
        summary.append({
            "source_file": source,
            "unique_genome_proteins": len(
                {r["genome_protein"] for r in rows}
            ),
            "unique_reference_sequences": len(
                {r["reference_id"] for r in rows}
            ),
            "protein_reference_pairs": len(rows)
        })

    summary.sort(
        key=lambda r: (-r["unique_genome_proteins"], r["source_file"])
    )

    write_tsv(
        outdir / "06_summary_by_reference_file.tsv",
        summary,
        ["source_file", "unique_genome_proteins",
         "unique_reference_sequences", "protein_reference_pairs"]
    )

    print("\nAnalysis complete.")
    print(f"Reference sequences: {ref_number}")
    print(f"Protein-reference pairs before filtering: {len(pair_rows)}")
    print(f"Filtered protein-reference pairs: {len(filtered)}")
    print(f"Unique candidate genome proteins: {len(unique_rows)}")
    print(f"Results folder: {outdir}")


if __name__ == "__main__":
    main()
