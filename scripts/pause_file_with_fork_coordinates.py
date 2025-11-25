#!/usr/bin/env python3
"""
Script to add fork coordinates to a pause file based on left/right fork information.

Usage:
    python pause_file_with_fork_coordinates.py \
        --pause_file /path/to/test_pause_file.txt \
        --left_fork_file /path/to/test_left_fork.txt \
        --right_fork_file /path/to/test_right_fork.txt \
        --output_file /path/to/test_updated_pause_file.txt
"""

import pandas as pd
import argparse


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Update a pause file with fork information.")
    parser.add_argument("--pause_file", required=True, help="Path to the pause file")
    parser.add_argument("--left_fork_file", required=True, help="Path to the left fork file")
    parser.add_argument("--right_fork_file", required=True, help="Path to the right fork file")
    parser.add_argument("--output_file", required=True, help="Path to save the updated pause file")
    return parser.parse_args()


def load_fork_file(fork_file, fork_type):
    """Load fork file and return a DataFrame with relevant columns."""
    cols = ["contig", f"{fork_type}_start", f"{fork_type}_end", "readID", "contig2", "x1", "x2", "strand"]
    df = pd.read_csv(fork_file, sep=r'\s+', header=None, names=cols)
    
    # Convert only the start/end columns to int
    df[f"{fork_type}_start"] = df[f"{fork_type}_start"].astype(int)
    df[f"{fork_type}_end"] = df[f"{fork_type}_end"].astype(int)
    
    return df[["readID", f"{fork_type}_start", f"{fork_type}_end"]]


def assign_fork(row, lf_df, rf_df):
    """Assign left or right fork coordinates to a row of the pause file."""
    ps = row["pauseSite"]
    readID = row["readID"]

    match_L = lf_df[(lf_df["readID"] == readID) & (lf_df["lf_start"] <= ps) & (lf_df["lf_end"] >= ps)]
    if not match_L.empty:
        row["direction"] = "L"
        row["left_fork_start"] = match_L.iloc[0]["lf_start"]
        row["left_fork_end"] = match_L.iloc[0]["lf_end"]
        return row

    match_R = rf_df[(rf_df["readID"] == readID) & (rf_df["rf_start"] <= ps) & (rf_df["rf_end"] >= ps)]
    if not match_R.empty:
        row["direction"] = "R"
        row["right_fork_start"] = match_R.iloc[0]["rf_start"]
        row["right_fork_end"] = match_R.iloc[0]["rf_end"]
        return row

    return row


def update_detect_index(row):
    """Update the detectIndex column based on fork direction and coordinates."""
    parts = row["detectIndex"].split("_")
    parts[5] = row["direction"]  # update dirn
    if row["direction"] == "L":
        parts[6] = str(row["left_fork_start"])
        parts[7] = str(row["left_fork_end"])
    else:  # R
        parts[6] = str(row["right_fork_start"])
        parts[7] = str(row["right_fork_end"])
    return "_".join(parts)


def main():
    # Parse arguments
    args = parse_args()

    # Load fork files
    lf = load_fork_file(args.left_fork_file, "lf")
    rf = load_fork_file(args.right_fork_file, "rf")

    # Load pause file, skipping comment lines
    pause = pd.read_csv(args.pause_file, comment="#", delim_whitespace=True)
    pause["pauseSite"] = pause["pauseSite"].astype(int)

    # Keep track of original number of rows
    original_rows = pause.shape[0]

    # Extract readID from detectIndex
    pause["readID"] = pause["detectIndex"].apply(lambda x: x.split("_")[0])

    # Initialize fork columns
    pause["left_fork_start"] = "NA"
    pause["left_fork_end"] = "NA"
    pause["right_fork_start"] = "NA"
    pause["right_fork_end"] = "NA"
    pause["direction"] = "NA"

    # Assign forks
    pause = pause.apply(assign_fork, axis=1, lf_df=lf, rf_df=rf)

    # Keep only rows with a fork match
    pause = pause[pause["direction"] != "NA"].copy()

    # Print number of skipped rows
    skipped_rows = original_rows - pause.shape[0]
    print(f"Number of lines skipped (excluding comment lines): {skipped_rows}")

    # Update detectIndex
    pause["detectIndex"] = pause.apply(update_detect_index, axis=1)

    # Drop temporary readID column
    pause.drop(columns=["readID"], inplace=True)

    # Save updated pause file
    pause.to_csv(args.output_file, sep="\t", index=False)
    print(f"Updated pause file saved to: {args.output_file}")


if __name__ == "__main__":
    main()
