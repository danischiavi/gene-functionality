#!/usr/bin/env python3
"""
Standalone MaxEntScan.

This pipeline is fully independent of the SpliceAI evaluation pipeline.
It runs on ALL regions in the input CSVs (no subsampling) and produces
MaxEntScan splice-site scores in a single TSV per (dataset, set) pair.

============================================================================
STAGE 2 — Flanked sequence extraction (no subsampling)
============================================================================
For every region in each of the 10 input CSVs:

  - Extract the exon plus 100 nt of flanking intronic sequence on each side
    from the GRCh38 reference FASTA (using the chromosome/start/end coords
    from the CSV, treating coords as 1-based inclusive).
  - For positive CSVs (which have a Strand column): if Strand == "-",
    reverse-complement the entire extracted sequence so it is in transcribed
    (5'→3') orientation. This is necessary because the reference FASTA
    returns the forward strand regardless of which strand the gene is on,
    and MaxEntScan can only score splice motifs in transcribed orientation.
  - For negative-control CSVs (no Strand column): intergenic regions have
    no biological transcribed direction, so they are extracted on the
    forward strand as-is. No reverse-complementation.

Output: one .fa + one _metadata.csv per (dataset, set) in:
  data/datasets/maxentscan_eval/02_flanked_fasta_full/

============================================================================
STAGE 3 — MaxEntScan scoring
============================================================================
For every extracted region:

  - Slide MaxEntScan's acceptor (score3, 23-nt window) and donor (score5,
    9-nt window) scoring windows across the entire exon body plus 50 nt of
    intronic flank on each side. Both motif types use the same symmetric
    scan range — no distinction between canonical and cryptic positions.
  - Record the max acceptor score, max donor score, the positions at which
    each maximum occurred, and a combined score = max_acceptor + max_donor.

Output: one .tsv per (dataset, set) in:
  data/datasets/maxentscan_eval/scores_full/

Columns:
  region_id              -- matches the ID column in the input CSV
  max_acceptor           -- best score3 across the scan range
  max_acceptor_position  -- 0-based position of the AG anchor at the max
                            (in coordinates of the extracted flanked sequence,
                             where the exon's 5' edge sits at position 100)
  max_donor              -- best score5 across the scan range
  max_donor_position     -- 0-based position of the GT anchor at the max
                            (same coordinate convention as max_acceptor_position)
  splice_combined        -- max_acceptor + max_donor
                            (NaN if either component is NaN)

============================================================================
Safety
============================================================================
At startup, the script checks both output directories. If either is
non-empty, its contents are copied to a timestamped sibling folder before
any new files are written. So re-running the script never silently loses
previous results — you can recover them from the backup folder.

============================================================================
Environment
============================================================================
This script uses the Python 3.13 module with maxentpy and pyfaidx installed
to ~/.local/lib/python3.13/site-packages/.

If pyfaidx is not installed there yet, install with:
  module load python/3.13.0-e2wrr3t
  pip install --user pyfaidx

(maxentpy should already be installed from the previous MaxEntScan work.)
"""

from __future__ import annotations

import shutil
import sys
import time
from datetime import datetime
from pathlib import Path

import pandas as pd

try:
    from pyfaidx import Fasta
except ImportError:
    print(
        "FATAL: pyfaidx is not installed in this Python environment.\n"
        "       Install it with: pip install --user pyfaidx\n"
        "       (using the Python 3.13 module: module load python/3.13.0-e2wrr3t)",
        file=sys.stderr,
    )
    sys.exit(1)

from maxentpy import maxent
from maxentpy.maxent import load_matrix3, load_matrix5


# =========================================================================
# CONFIGURATION
# =========================================================================

# Project root, derived from this file's location (scripts/ -> project root),
# so every path below is relative to the repository rather than absolute.
PROJECT_DIR = Path(__file__).resolve().parent.parent

REFERENCE_FASTA = PROJECT_DIR / "data" / "datasets" / "hg38.fa"
INPUT_DIR = PROJECT_DIR / "data" / "datasets" / "maxentscan_eval" / "input"

# Dataset → input filenames (update if the filenames on Aoraki differ)
DATASETS = {
    "protein_exon2": {
        "positives": "functional-protein-exon2-dataset_with_strand.csv",
        "controls":  "protein-exon2-negative-control.csv",
    },
    "protein_exon3": {
        "positives": "functional-protein-exon3-dataset_with_strand.csv",
        "controls":  "protein-exon3-negative-control.csv",
    },
    "lncrna_exon1": {
        "positives": "functional-lncrna-exon1-dataset_with_strand.csv",
        "controls":  "lncrna-exon1-negative-control.csv",
    },
    "lncrna_exon2": {
        "positives": "functional-lncrna-exon2-dataset_with_strand.csv",
        "controls":  "lncrna-exon2-negative-control.csv",
    },
    "sncrna": {
        "positives": "functional-short-ncrna-dataset_with_strand.csv",
        "controls":  "short-ncrna-negative-control-dataset.csv",
    },
}

# Output folders (under maxentscan_eval/, fully separate from spliceai_eval/)
MAXENT_DIR = PROJECT_DIR / "data" / "datasets" / "maxentscan_eval"
STAGE2_OUT = MAXENT_DIR / "02_flanked_fasta_full"
STAGE3_OUT = MAXENT_DIR / "scores_full"

# Flank sizes
FLANK_BP_EXTRACT = 100   # stage 2: nt of intronic flank extracted on each side
FLANK_BP_SCAN    = 50    # stage 3: nt of flank into which motif anchors may slide

# MaxEntScan motif window structure (Yeo & Burge 2004)
SCORE3_WINDOW = 23   # 20 intronic + AG + 3 exonic
SCORE5_WINDOW = 9    # 3 exonic + GT + 4 intronic
ACCEPTOR_WINDOW_LEFT_OF_ANCHOR = 20   # nt to the left of AG anchor in score3 window
DONOR_WINDOW_LEFT_OF_ANCHOR    = 3    # nt to the left of GT anchor in score5 window


# =========================================================================
# HELPERS
# =========================================================================

# Lookup table for reverse complement (uppercase + lowercase + N)
_COMPLEMENT_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    """Reverse-complement a DNA sequence, preserving case for N's."""
    return seq.translate(_COMPLEMENT_TABLE)[::-1]


def backup_existing_output(output_dir: Path) -> Path | None:
    """If output_dir exists and contains files, copy it to a timestamped
    sibling folder.

    Naming: <output_dir>_backup_YYYYMMDD_HHMMSS

    Returns the backup path on success, or None if no backup was needed.
    Each backup carries its own timestamp so repeated re-runs leave
    distinct snapshots that don't overwrite each other.
    """
    if not output_dir.exists():
        return None
    if not any(output_dir.iterdir()):
        return None
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_dir = output_dir.parent / f"{output_dir.name}_backup_{timestamp}"
    if backup_dir.exists():
        raise FileExistsError(
            f"Backup target {backup_dir} already exists; refusing to overwrite."
        )
    shutil.copytree(output_dir, backup_dir)
    return backup_dir


# =========================================================================
# STAGE 2 — extract flanked sequences from hg38
# =========================================================================

def stage2_extract(
    dataset: str,
    set_name: str,
    csv_path: Path,
    fasta: Fasta,
) -> int:
    """Extract every region in csv_path with FLANK_BP_EXTRACT nt of flank on
    each side. Reverse-complements minus-strand entries (positives only).

    Returns the number of regions successfully written.

    Writes:
      STAGE2_OUT / {dataset}_{set_name}.fa
      STAGE2_OUT / {dataset}_{set_name}_metadata.csv
    """
    df = pd.read_csv(csv_path)
    has_strand = "Strand" in df.columns
    if not has_strand:
        print(
            f"  [stage2] {set_name}: no Strand column found "
            f"(intergenic control) — using forward strand for all entries"
        )

    fa_out = STAGE2_OUT / f"{dataset}_{set_name}.fa"
    meta_out = STAGE2_OUT / f"{dataset}_{set_name}_metadata.csv"
    fa_out.parent.mkdir(parents=True, exist_ok=True)

    meta_rows: list[dict] = []
    n_written = 0
    n_skipped_bounds = 0
    n_skipped_unknown_chrom = 0
    n_minus_revcomp = 0
    t0 = time.time()

    with fa_out.open("w") as fh_fa:
        for _, row in df.iterrows():
            rid = str(row["ID"])
            chrom = str(row["Chromosome"])
            start_1based = int(row["Start"])
            end_1based = int(row["End"])
            strand = str(row["Strand"]).strip() if has_strand else "+"

            # CSV coords are 1-based inclusive [Start, End].
            # Convert to 0-based half-open [start_0, end_0) for pyfaidx.
            exon_start_0 = start_1based - 1
            exon_end_0 = end_1based
            exon_length = exon_end_0 - exon_start_0  # end - start (both 0-based)

            # Add flank on each side
            extract_start = exon_start_0 - FLANK_BP_EXTRACT
            extract_end = exon_end_0 + FLANK_BP_EXTRACT

            # Bounds checks
            try:
                chrom_len = len(fasta[chrom])
            except KeyError:
                # Chromosome not in the reference (e.g., chrM, chrUn_*)
                n_skipped_unknown_chrom += 1
                continue
            if extract_start < 0 or extract_end > chrom_len:
                # Exon too close to the chromosome end to take a full flank
                n_skipped_bounds += 1
                continue

            # Extract the flanked region from the forward strand
            seq = str(fasta[chrom][extract_start:extract_end].seq).upper()

            # For minus-strand positives, reverse-complement so we end up in
            # transcribed (5'→3') orientation. Note: after RC the exon body
            # still occupies the middle of the sequence, exactly
            # FLANK_BP_EXTRACT nt from either end, so the boundary positions
            # below remain correct in either case.
            if strand == "-":
                seq = reverse_complement(seq)
                n_minus_revcomp += 1

            # Boundaries inside the extracted sequence
            boundary_5 = FLANK_BP_EXTRACT          # 0-based first exonic base
            boundary_3 = FLANK_BP_EXTRACT + exon_length  # 1-based last exonic base + 1
            #   → 0-based last exonic base = boundary_3 - 1

            # FASTA record. Header layout matches what SpliceAI pipeline
            # produced, so downstream stage 3 reader can be reused as-is.
            header = (
                f">{set_name.upper()}|{rid}|{chrom}:{start_1based}-{end_1based}"
                f"|strand:{strand}|exon:{boundary_5}-{boundary_3}"
            )
            fh_fa.write(header + "\n")
            # 80-col line wrapping
            for i in range(0, len(seq), 80):
                fh_fa.write(seq[i:i + 80] + "\n")

            meta_rows.append({
                "ID": rid,
                "Chromosome": chrom,
                "Start": start_1based,
                "End": end_1based,
                "Strand": strand,
                "exon_length": exon_length,
                "boundary_5_in_seq": boundary_5,
                "boundary_3_in_seq": boundary_3,
                "extracted_length": len(seq),
            })
            n_written += 1

    pd.DataFrame(meta_rows).to_csv(meta_out, index=False)

    dt = time.time() - t0
    msg = f"  [stage2] {set_name}: extracted {n_written} regions in {dt:.1f}s"
    if n_minus_revcomp:
        msg += f" ({n_minus_revcomp} reverse-complemented)"
    if n_skipped_bounds:
        msg += f"  [{n_skipped_bounds} skipped: near chromosome end]"
    if n_skipped_unknown_chrom:
        msg += f"  [{n_skipped_unknown_chrom} skipped: unknown chromosome]"
    print(msg)
    return n_written


# =========================================================================
# STAGE 3 — score MaxEntScan
# =========================================================================

def scan_motif(
    seq: str,
    scan_start: int,
    scan_end: int,
    window_size: int,
    window_offset_left: int,
    scorer,
    matrix,
) -> tuple[float, int]:
    """Slide a motif anchor across [scan_start, scan_end] inclusive and
    return (max_score, anchor_position_at_max).

    `window_offset_left` is the number of nucleotides of the scoring window
    that sit to the LEFT of the anchor:
      acceptor (score3): AG anchor at right edge of 20-nt intron chunk → 20
      donor    (score5): GT anchor at right edge of 3-nt exon chunk    → 3

    Returns (nan, -1) if no valid window could be scored (every window
    contained N's, or the scan range fell outside the sequence).
    """
    best_score = float("-inf")
    best_anchor = -1
    L = len(seq)
    for anchor in range(scan_start, scan_end + 1):
        start = anchor - window_offset_left
        end = start + window_size
        if start < 0 or end > L:
            continue
        window = seq[start:end]
        if "N" in window:
            continue
        try:
            score = float(scorer(window, matrix=matrix))
        except Exception:
            continue
        if score > best_score:
            best_score = score
            best_anchor = anchor
    if best_score == float("-inf"):
        return float("nan"), -1
    return best_score, best_anchor


def read_flanked_fasta(path: Path) -> dict[str, str]:
    """Read a FASTA produced by stage 2, keyed by region ID."""
    records: dict[str, str] = {}
    current_key: str | None = None
    current_chunks: list[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if current_key is not None:
                    records[current_key] = "".join(current_chunks)
                parts = line[1:].split("|")
                current_key = parts[1] if len(parts) >= 2 else line[1:].split()[0]
                current_chunks = []
            else:
                current_chunks.append(line.upper())
        if current_key is not None:
            records[current_key] = "".join(current_chunks)
    return records


def stage3_score(dataset: str, set_name: str, m3, m5) -> None:
    """Read stage 2 output for (dataset, set_name), score MaxEntScan, write TSV."""
    fa = STAGE2_OUT / f"{dataset}_{set_name}.fa"
    meta = STAGE2_OUT / f"{dataset}_{set_name}_metadata.csv"
    if not fa.exists() or not meta.exists():
        print(f"  [stage3] {set_name}: stage 2 outputs missing — skipping")
        return

    seqs = read_flanked_fasta(fa)
    df = pd.read_csv(meta)

    rows: list[dict] = []
    n_missing_seq = 0
    n_skipped_bounds = 0
    t0 = time.time()

    for _, m in df.iterrows():
        rid = str(m["ID"])
        seq = seqs.get(rid)
        if seq is None:
            n_missing_seq += 1
            continue

        five_edge = int(m["boundary_5_in_seq"])         # 0-based first exonic base
        three_edge = int(m["boundary_3_in_seq"]) - 1     # 0-based last exonic base

        # Scan range: exon body + FLANK_BP_SCAN on each side
        scan_start = five_edge - FLANK_BP_SCAN
        scan_end = three_edge + FLANK_BP_SCAN

        if scan_start < 0 or scan_end >= len(seq):
            n_skipped_bounds += 1
            continue

        acc_score, acc_pos = scan_motif(
            seq, scan_start, scan_end,
            SCORE3_WINDOW, ACCEPTOR_WINDOW_LEFT_OF_ANCHOR,
            maxent.score3, m3,
        )
        don_score, don_pos = scan_motif(
            seq, scan_start, scan_end,
            SCORE5_WINDOW, DONOR_WINDOW_LEFT_OF_ANCHOR,
            maxent.score5, m5,
        )

        # Combined score: NaN propagates (nan + x = nan in Python float math)
        combined = acc_score + don_score

        rows.append({
            "region_id":             rid,
            "max_acceptor":          acc_score,
            "max_acceptor_position": acc_pos,
            "max_donor":             don_score,
            "max_donor_position":    don_pos,
            "splice_combined":       combined,
        })

    out = STAGE3_OUT / f"{dataset}_{set_name}.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)

    dt = time.time() - t0
    msg = f"  [stage3] {set_name}: scored {len(rows)} regions in {dt:.1f}s  -> {out.name}"
    if n_missing_seq:
        msg += f"  [{n_missing_seq} missing seq]"
    if n_skipped_bounds:
        msg += f"  [{n_skipped_bounds} scan range out of sequence]"
    print(msg)


# =========================================================================
# DRIVER
# =========================================================================

def main() -> None:
    # Sanity checks on prerequisites
    if not REFERENCE_FASTA.exists():
        print(f"FATAL: reference FASTA not found at {REFERENCE_FASTA}", file=sys.stderr)
        sys.exit(1)
    if not INPUT_DIR.exists():
        print(f"FATAL: input dir not found at {INPUT_DIR}", file=sys.stderr)
        sys.exit(1)

    print(f"Reference:   {REFERENCE_FASTA}")
    print(f"Input dir:   {INPUT_DIR}")
    print(f"Stage 2 out: {STAGE2_OUT}")
    print(f"Stage 3 out: {STAGE3_OUT}")
    print(f"Flank (extract): {FLANK_BP_EXTRACT} nt on each side of exon")
    print(f"Flank (scan):    {FLANK_BP_SCAN} nt on each side of exon")
    print()

    # Auto-backup of any existing outputs
    for d in (STAGE2_OUT, STAGE3_OUT):
        b = backup_existing_output(d)
        if b is not None:
            print(f"Backed up {d.name}/ -> {b.name}/")
        else:
            print(f"No prior output to back up in {d.name}/")
    print()

    # Open the reference FASTA once
    print("Opening reference FASTA (may take a moment if .fai needs to be built)...")
    fasta = Fasta(str(REFERENCE_FASTA))
    print(f"Reference loaded ({len(fasta.keys())} contigs)")

    # ---- STAGE 2 -------------------------------------------------------
    print("\n========== STAGE 2: extract flanked sequences ==========")
    stage2_t0 = time.time()
    for dataset, files in DATASETS.items():
        print(f"\n=== {dataset} ===")
        for set_name, fname in files.items():
            csv_path = INPUT_DIR / fname
            if not csv_path.exists():
                print(f"  [stage2] {fname} not found — skipping")
                continue
            stage2_extract(dataset, set_name, csv_path, fasta)
    print(f"\nStage 2 complete in {time.time() - stage2_t0:.1f}s")

    # ---- STAGE 3 -------------------------------------------------------
    print("\n========== STAGE 3: MaxEntScan scoring ==========")
    print("Loading MaxEntScan matrices...")
    m3 = load_matrix3()
    m5 = load_matrix5()

    stage3_t0 = time.time()
    for dataset in DATASETS:
        print(f"\n=== {dataset} ===")
        for set_name in ("positives", "controls"):
            stage3_score(dataset, set_name, m3, m5)
    print(f"\nStage 3 complete in {time.time() - stage3_t0:.1f}s")

    print("\nAll done.")


if __name__ == "__main__":
    main()
