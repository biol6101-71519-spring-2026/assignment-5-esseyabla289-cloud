#!/usr/bin/env python3

#!/usr/bin/env python3
import argparse

def summarize(input_files, output_file):
    total_lines = 0
    total_chars = 0

    for input_file in input_files:       # loop over each file
        with open(input_file, 'r') as f:
            lines = f.readlines()
        total_lines += len(lines)
        total_chars += sum(len(line) for line in lines)

    with open(output_file, 'w') as out:
        out.write("Summary Report\n")
        out.write(f"Total lines: {total_lines}\n")
        out.write(f"Total characters: {total_chars}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarize alignment stats files")
    parser.add_argument("--input", required=True, nargs='+', help="Input files")
    parser.add_argument("--output", required=True, help="Output summary file")
    args = parser.parse_args()
    summarize(args.input, args.output)
