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
import os


def parse_args():
    parser = argparse.ArgumentParser(description="Update a pause file with fork information.")
    parser.add_argument("--pause_file", required=True)
    parser.add_argument("--left_fork_file", required=True)
    parser.add_argument("--right_fork_file", required=True)
    parser.add_argument("--output_file", required=True)
    return parser.parse_args()


def load_fork_file(fork_file, fork_type):
    """
    Fork file columns (no header): contig fork_start fork_end read_id contig start_read end_read strand
    """
    cols = ["contig", f"{fork_type}_start", f"{fork_type}_end", "readID",
            "contig2", "start_read", "end_read", "strand"]
    df = pd.read_csv(fork_file, sep=r"\s+", header=None, names=cols)
    df[f"{fork_type}_start"] = df[f"{fork_type}_start"].astype(int)
    df[f"{fork_type}_end"] = df[f"{fork_type}_end"].astype(int)

    # Map strand
    df["strand"] = df["strand"].map({"fwd": "+", "rev": "-"})

    # Set direction
    df["direction"] = "L" if fork_type == "lf" else "R"

    return df


def assign_fork(row, lf_df, rf_df):
    ps = row["pauseSite"]
    readID = row["readID"]

    # Check left forks
    match_L = lf_df[(lf_df["readID"] == readID) &
                    (lf_df["lf_start"] <= ps) &
                    (lf_df["lf_end"] >= ps)]
    if not match_L.empty:
        r = match_L.iloc[0]
        row["direction"] = "L"
        row["left_fork_start"] = r["lf_start"]
        row["left_fork_end"] = r["lf_end"]
        row["strand"] = r["strand"]
        row["alignLen"] = r["end_read"]
        row["contig"] = r["contig"]
        row["start_read"] = r["start_read"]
        row["end_read"] = r["end_read"]
        return row

    # Check right forks
    match_R = rf_df[(rf_df["readID"] == readID) &
                    (rf_df["rf_start"] <= ps) &
                    (rf_df["rf_end"] >= ps)]
    if not match_R.empty:
        r = match_R.iloc[0]
        row["direction"] = "R"
        row["right_fork_start"] = r["rf_start"]
        row["right_fork_end"] = r["rf_end"]
        row["strand"] = r["strand"]
        row["alignLen"] = r["end_read"]
        row["contig"] = r["contig"]
        row["start_read"] = r["start_read"]
        row["end_read"] = r["end_read"]
        return row

    return row


def construct_detect_index(row):
    """
    Construct detectIndex:
    <readID>_<contig>_<start_read>_<end_read>_<strand>_<direction>_<fork_start>_<fork_end>
    """
    fork_start = row["left_fork_start"] if row["direction"] == "L" else row["right_fork_start"]
    fork_end = row["left_fork_end"] if row["direction"] == "L" else row["right_fork_end"]

    return f"{row['readID']}_{row['contig']}_{row['start_read']}_{row['end_read']}_{row['strand']}_{row['direction']}_{fork_start}_{fork_end}"


def add_non_paused_forks(lf, rf, used_forks, template_columns):
    rows = []

    for df, start_col, end_col, direction in [
        (lf, "lf_start", "lf_end", "L"),
        (rf, "rf_start", "rf_end", "R")
    ]:
        for _, r in df.iterrows():
            key = (r["readID"], r[start_col], r[end_col], direction)
            if key in used_forks:
                continue

            # Start row with template columns
            row = {col: "NA" for col in template_columns}

            # Fill actual fork values
            row["contig"] = r["contig"]
            row["strand"] = r["strand"]
            row["alignLen"] = r["end_read"]
            row["start_read"] = r["start_read"]
            row["end_read"] = r["end_read"]
            if direction == "L":
                row["left_fork_start"] = r["lf_start"]
                row["left_fork_end"] = r["lf_end"]
            else:
                row["right_fork_start"] = r["rf_start"]
                row["right_fork_end"] = r["rf_end"]

            row["readID"] = r["readID"]
            row["direction"] = direction
            row["keep"] = False

            # Construct detectIndex after all values are filled
            row["detectIndex"] = construct_detect_index(row)

            rows.append(row)

    df_non_paused = pd.DataFrame(rows)

    # Remove temporary readID column
    df_non_paused.drop(columns=["readID"], inplace=True)

    return df_non_paused


def main():
    args = parse_args()

    pause = pd.read_csv(args.pause_file, comment="#", delim_whitespace=True)
    pause["pauseSite"] = pause["pauseSite"].astype(int)

    lf = load_fork_file(args.left_fork_file, "lf")
    rf = load_fork_file(args.right_fork_file, "rf")

    # Extract readID from detectIndex
    pause["readID"] = pause["detectIndex"].apply(lambda x: x.split("_")[0])

    # Initialize fork columns
    pause["left_fork_start"] = "NA"
    pause["left_fork_end"] = "NA"
    pause["right_fork_start"] = "NA"
    pause["right_fork_end"] = "NA"
    pause["direction"] = "NA"
    pause["strand"] = "NA"
    pause["alignLen"] = "NA"
    pause["start_read"] = "NA"
    pause["end_read"] = "NA"

    # Assign forks to paused rows
    pause = pause.apply(assign_fork, axis=1, lf_df=lf, rf_df=rf)

    paused = pause[pause["direction"] != "NA"].copy()
    paused["detectIndex"] = paused.apply(construct_detect_index, axis=1)
    paused["keep"] = True

    # Track used forks
    used_forks = set()
    for _, r in paused.iterrows():
        if r["direction"] == "L":
            used_forks.add((r["readID"], r["left_fork_start"], r["left_fork_end"], "L"))
        else:
            used_forks.add((r["readID"], r["right_fork_start"], r["right_fork_end"], "R"))

    paused.drop(columns=["readID"], inplace=True)

    # Add non-paused forks
    non_paused = add_non_paused_forks(lf, rf, used_forks, paused.columns)

    output = pd.concat([paused, non_paused], ignore_index=True)

    # Save to output
    output.to_csv(args.output_file, sep="\t", index=False)
    print(f"Updated pause file saved to: {os.path.abspath(args.output_file)}")


if __name__ == "__main__":
    main()
