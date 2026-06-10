"""
Compute per-sequence mononucleotide frequencies for all positive and control
sequences across the 5 datasets (protein exon 2/3, lncRNA exon 1/2, sncRNA).

For each sequence:
    - Percent T, A, C, G (excluding N and any non-ACGT bases from the
      denominator). N's and non-ACGT characters do not contribute.

This is a streamlined version of `compute_nucleotide_composition.py`: it
keeps the mononucleotide-frequency calculation exactly as before but no
longer computes GC content in the first or last half of each sequence.
Because `base_percentages` is unchanged, the resulting T/A/C/G values are
byte-identical to the corresponding columns in the previous composition
output.

Runs on the FULL input CSVs (all 1,000 positives + all ~9,600 controls
per dataset).

Output per dataset x kind (positives/controls):
    <dataset>_<kind>_nucleotide_frequencies.csv
    Columns: ID, T, A, C, G

Plus a combined master file with two extra leading columns (dataset, kind):
    ALL_nucleotide_frequencies.csv

If the previous composition output (`06_composition/`) is still present
alongside the new output, the script automatically verifies that the new
T, A, C, G values match the previous run row-for-row and reports the
result at the end.

Usage:
    python compute_nucleotide_frequencies.py [input_dir] [output_dir]

Defaults:
    input_dir:  ../../data/datasets/spliceai_eval/input
    output_dir: ../../data/datasets/spliceai_eval/output/07_nucleotide_frequencies
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

# Match the dataset structure of the SpliceAI pipeline
DATASETS = {
    "protein_exon2": {
        "positive_csv": "functional-protein-exon2-dataset__1_.csv",
        "control_csv": "protein-exon2-negative-control.csv",
    },
    "protein_exon3": {
        "positive_csv": "functional-protein-exon3-dataset__1_.csv",
        "control_csv": "protein-exon3-negative-control.csv",
    },
    "lncrna_exon1": {
        "positive_csv": "functional-lncrna-exon1-dataset__2_.csv",
        "control_csv": "lncrna-exon1-negative-control.csv",
    },
    "lncrna_exon2": {
        "positive_csv": "functional-lncrna-exon2-dataset__2_.csv",
        "control_csv": "lncrna-exon2-negative-control.csv",
    },
    "sncrna": {
        "positive_csv": "functional-short-ncrna-dataset__2_.csv",
        "control_csv": "short-ncrna-negative-control-dataset__5_.csv",
    },
}


def base_percentages(sequence: str) -> dict:
    """Return a dict of T, A, C, G percentages.

    N and any non-ACGT bases are excluded from the denominator.
    Returns NaN for all four if the sequence has no valid bases.

    This function is unchanged from `compute_nucleotide_composition.py`,
    which is why the new T/A/C/G output is byte-identical to the previous
    run.
    """
    seq = sequence.upper()
    counts = {"A": 0, "C": 0, "G": 0, "T": 0}
    for base in seq:
        if base in counts:
            counts[base] += 1
    total_valid = sum(counts.values())
    if total_valid == 0:
        return {b: float("nan") for b in "ACGT"}
    return {b: counts[b] / total_valid * 100 for b in "ACGT"}


def compute_frequencies(input_csv: Path) -> pd.DataFrame:
    """Compute mononucleotide frequencies for every sequence in an input CSV.

    Returns a DataFrame with columns: ID, T, A, C, G.
    """
    df = pd.read_csv(input_csv)
    rows = []
    for _, row in df.iterrows():
        seq_id = row["ID"]
        seq = str(row["Sequence"])
        pct = base_percentages(seq)
        rows.append({
            "ID": seq_id,
            "T": pct["T"],
            "A": pct["A"],
            "C": pct["C"],
            "G": pct["G"],
        })
    return pd.DataFrame(rows)


# --- Verification --------------------------------------------------------

def verify_against_old(new_output_dir: Path, old_output_dir: Path) -> None:
    """Compare the new per-dataset outputs against the previous
    `_composition.csv` files. Reports any T/A/C/G differences.

    Silently skips if the old output dir doesn't exist (e.g., it was
    moved or deleted) -- verification is opportunistic, not required.
    """
    print(f"\n[verify] Comparing against previous composition output at:")
    print(f"         {old_output_dir}")
    if not old_output_dir.exists():
        print(f"[verify] Old output dir not found -- skipping verification.")
        return

    cols = ["T", "A", "C", "G"]
    all_ok = True
    n_compared = 0

    for new_path in sorted(new_output_dir.glob("*_nucleotide_frequencies.csv")):
        # The combined master file uses a different name; skip in this loop
        if new_path.name.startswith("ALL_"):
            continue
        base = new_path.stem.replace("_nucleotide_frequencies", "")
        old_path = old_output_dir / f"{base}_composition.csv"
        if not old_path.exists():
            print(f"  [skip] {base}: no matching old file ({old_path.name})")
            continue

        old_df = pd.read_csv(old_path)
        new_df = pd.read_csv(new_path)

        # Align by ID to be robust to any row ordering differences
        if "ID" in old_df.columns and "ID" in new_df.columns:
            merged = old_df[["ID"] + cols].merge(
                new_df[["ID"] + cols],
                on="ID",
                suffixes=("_old", "_new"),
                how="outer",
                indicator=True,
            )
            if (merged["_merge"] != "both").any():
                missing = (merged["_merge"] != "both").sum()
                print(f"  [DIFF] {base}: {missing} IDs present in one file but not the other")
                all_ok = False
                continue
            differences = 0
            for c in cols:
                differences += (merged[f"{c}_old"] != merged[f"{c}_new"]).sum()
            if differences == 0:
                print(f"  [ok]   {base}: T, A, C, G match exactly ({len(merged):,} rows)")
            else:
                print(f"  [DIFF] {base}: {differences} cell-level differences across T/A/C/G")
                all_ok = False
        else:
            # Fallback: compare row-by-row (assumes same ordering)
            if old_df[cols].equals(new_df[cols]):
                print(f"  [ok]   {base}: T, A, C, G match exactly ({len(new_df):,} rows)")
            else:
                diff_rows = (old_df[cols] != new_df[cols]).any(axis=1).sum()
                print(f"  [DIFF] {base}: {diff_rows} rows differ in T/A/C/G")
                all_ok = False
        n_compared += 1

    if n_compared == 0:
        print("[verify] No old files found to compare -- skipping.")
    elif all_ok:
        print(f"\n[verify] All {n_compared} datasets match the previous output exactly for T, A, C, G.")
    else:
        print(f"\n[verify] WARNING: at least one dataset's frequencies differ from the previous output. Investigate before using.")


# --- Main ---------------------------------------------------------------

def main() -> None:
    # Allow overriding input/output paths from the command line
    if len(sys.argv) >= 2:
        input_dir = Path(sys.argv[1])
    else:
        # Default: same input dir used by the SpliceAI pipeline
        script_dir = Path(__file__).resolve().parent
        project_dir = script_dir.parent  # scripts/ -> project root
        input_dir = project_dir / "data" / "datasets" / "spliceai_eval" / "input"

    if len(sys.argv) >= 3:
        output_dir = Path(sys.argv[2])
    else:
        script_dir = Path(__file__).resolve().parent
        project_dir = script_dir.parent  # scripts/ -> project root
        output_dir = (
            project_dir / "data" / "datasets" / "spliceai_eval" /
            "output" / "07_nucleotide_frequencies"
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Input:  {input_dir}")
    print(f"Output: {output_dir}\n")

    all_rows = []
    for ds_name, ds_cfg in DATASETS.items():
        print(f"=== {ds_name} ===")
        for kind, csv_key in [("positives", "positive_csv"),
                              ("controls", "control_csv")]:
            csv_path = input_dir / ds_cfg[csv_key]
            if not csv_path.exists():
                print(f"  SKIPPED {kind}: file not found ({csv_path.name})")
                continue
            df_out = compute_frequencies(csv_path)
            out_path = output_dir / f"{ds_name}_{kind}_nucleotide_frequencies.csv"
            df_out.to_csv(out_path, index=False)
            print(f"  {kind:9s}: {len(df_out)} sequences -> {out_path.name}")
            # Tag rows for the combined master file
            df_out_tagged = df_out.copy()
            df_out_tagged.insert(0, "kind", kind)
            df_out_tagged.insert(0, "dataset", ds_name)
            all_rows.append(df_out_tagged)
        print()

    if all_rows:
        master = pd.concat(all_rows, ignore_index=True)
        master_path = output_dir / "ALL_nucleotide_frequencies.csv"
        master.to_csv(master_path, index=False)
        print(f"Combined master: {len(master)} sequences -> ALL_nucleotide_frequencies.csv")

    # Automatic verification against the previous composition output
    old_output_dir = output_dir.parent / "06_composition"
    verify_against_old(output_dir, old_output_dir)


if __name__ == "__main__":
    main()
