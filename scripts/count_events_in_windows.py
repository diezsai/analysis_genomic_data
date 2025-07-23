import sys
import argparse

def parse_fai(fai_file):
    """
    Parse the FASTA index (.fai) file to get contig sizes.
    Returns:
        dict: A dictionary mapping contig names to their lengths.
    """
    contig_sizes = {}
    with open(fai_file, 'r') as f:
        for line in f:
            if not line.strip():
                continue  # Skip empty lines

            parts = line.strip().split()
            if len(parts) < 2:
                print(f"Warning: Skipping malformed line in FAI file: {line.strip()}")
                continue

            try:
                contig_name = parts[0]
                contig_length = int(parts[1])
                contig_sizes[contig_name] = contig_length
            except ValueError:
                print(f"Error: Invalid contig size on line: {line.strip()}")
                continue

    return contig_sizes

def parse_events(events_file):
    """
    Parse the file containing genomic events (e.g., origin positions).
    Returns:
        dict: A dictionary mapping contigs to lists of event positions.
    """
    events = {}
    with open(events_file, 'r') as f:
        header = True
        for line in f:
            line = line.strip()
            if not line:
                continue  # Skip empty lines

            if header:
                header = False  # Skip header line
                continue

            parts = line.split("\t")
            if len(parts) < 2:
                print(f"Warning: Skipping malformed line in events file: {line}")
                continue

            contig = parts[0]
            try:
                position = int(parts[1])
            except ValueError:
                print(f"Warning: Invalid position value: {parts[1]} on line: {line}")
                continue

            if contig not in events:
                events[contig] = []
            events[contig].append(position)
    return events

def count_events_in_window(contig, start, end, event_positions):
    """
    Count how many events fall within a specific genomic window.
    """
    count = 0
    for event in event_positions:
        if start <= event < end:
            count += 1
    return count

def sliding_window(contig_sizes, events, window_size, slide_size, output_file):
    """
    Perform sliding window analysis across the genome.
    """
    with open(output_file, 'w') as out:
        # Write header
        out.write("contig\tstart_window\tend_window\tevent_count\n")

        for contig, size in contig_sizes.items():
            print(f"Processing contig: {contig} of size {size}")

            event_positions = events.get(contig, [])
            print(f"Found {len(event_positions)} events for contig {contig}")

            for start in range(0, size, slide_size):
                end = min(start + window_size, size)
                event_count = count_events_in_window(contig, start, end, event_positions)

                print(f"Window: {start}-{end}, Event count: {event_count}")

                out.write(f"{contig}\t{start}\t{end}\t{event_count}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
        Count the number of genomic events (e.g., origins of replication) 
        that fall within specified genomic windows using a sliding window approach.
        
        Input:
            - FASTA .fai index file for contig sizes
            - Event positions file (e.g., origin midpoints)
        
        Output:
            - Tab-separated file listing event counts per window
    """)
    parser.add_argument("fai_file", help="Input FASTA index (.fai) file with contig sizes")
    parser.add_argument("events_file", help="Input file with genomic events (tab-separated: contig, position)")
    parser.add_argument("window_size", type=int, help="Window size to use for counting (in bases)")
    parser.add_argument("--slide_size", type=int, default=None,
                        help="Step size for sliding window (in bases). Defaults to window_size (no overlap).")
    parser.add_argument("output_file", help="Output file for writing the results")

    args = parser.parse_args()

    # Default slide_size to window_size if not provided
    if args.slide_size is None:
        args.slide_size = args.window_size
        print(f"No slide_size provided. Using window_size ({args.window_size}) as slide_size.")

    # Sanity checks
    if args.slide_size <= 0 or args.window_size <= 0:
        sys.exit("Error: Both window_size and slide_size must be positive integers.")

    if args.slide_size > args.window_size:
        print("Warning: slide_size is greater than window_size. Windows will be non-overlapping with gaps.")

    # Run the analysis
    contig_sizes = parse_fai(args.fai_file)
    events = parse_events(args.events_file)
    sliding_window(contig_sizes, events, args.window_size, args.slide_size, args.output_file)
