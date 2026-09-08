#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
@Author: Zijie Jiang
@Contact: jzjlab@163.com
@File: filter_coords.py
@Time: 2025/3/26 10:32
@Function: filter non-overlapping alignments from nucmer coords file
"""
import sys


def parse_nucmer_coords(coords_file):
    """
        parse nucmer coordinates file
    Args:
        coords_file: nucmer coords file path

    Returns:
        list of alignment dicts
    """

    alignments = []
    with open(coords_file, 'r') as f:
        for _ in range(5):
            next(f)

        for line in f:
            parts = [part.strip() for part in line.split('|')]

            fields = []
            for part in parts:
                subfields = [field for field in part.split() if field]
                fields.extend(subfields)

            if len(fields) < 13:
                continue

            ref_start = int(fields[0])
            ref_end = int(fields[1])
            query_start = int(fields[2])
            query_end = int(fields[3])
            ref_len = int(fields[4])
            query_len = int(fields[5])
            identity = float(fields[6])
            ref_name = fields[11]
            query_name = fields[12]

            direction = 'Forward' if query_start < query_end else 'Reverse'

            alignments.append({
                'ref_start': ref_start,
                'ref_end': ref_end,
                'query_start': query_start,
                'query_end': query_end,
                'ref_len': ref_len,
                'query_len': query_len,
                'identity': identity,
                'ref_name': ref_name,
                'query_name': query_name,
                'align_len': abs(ref_end - ref_start) + 1,
                'direction': direction
            })
    return alignments


def filter_non_overlapping(alignments):
    """
        filter non-overlapping alignments
    Args:
        alignments: list of alignment dicts

    Returns:
        filtered list of alignment dicts
    """

    ref_groups = {}
    for align in alignments:
        ref_name = align['ref_name']
        if ref_name not in ref_groups:
            ref_groups[ref_name] = []
        ref_groups[ref_name].append(align)

    results = []
    for ref_name, aligns in ref_groups.items():
        aligns.sort(key=lambda x: x['align_len'], reverse=True)

        covered_regions = []
        selected_aligns = []
        for align in aligns:
            current_start = align['ref_start']
            current_end = align['ref_end']
            is_overlapping = False

            for start, end in covered_regions:
                if not (current_end < start or current_start > end):
                    is_overlapping = True
                    break

            if not is_overlapping:
                covered_regions.append((current_start, current_end))
                selected_aligns.append(align)

        selected_aligns.sort(key=lambda x: x['ref_start'])
        results.extend(selected_aligns)

    return results


def write_output(results, output_file):
    with open(output_file, 'w') as f:
        f.write(
            "Ref_start\tRef_end\tQuery_start\tQuery_end\tRef_len\tQuery_len\tRef_name\tQuery_name\tDirection\tAlign_len\tIdentity\n")
        for align in results:
            line = f"{align['ref_start']}\t{align['ref_end']}\t{align['query_start']}\t{align['query_end']}\t" \
                   f"{align['ref_len']}\t{align['query_len']}\t{align['ref_name']}\t" \
                   f"{align['query_name']}\t{align['direction']}\t{align['align_len']}\t{align['identity']}\n"
            f.write(line)


def filter_coords(input_file, output_file):
    """
        filter non-overlapping alignments from nucmer coords file
    Args:
        input_file: nucmer coords file path
        output_file: output file path

    Returns:
        None
    """

    # parse coords file
    alignments = parse_nucmer_coords(input_file)

    # filter non-overlapping alignments
    filtered_alignments = filter_non_overlapping(alignments)

    # write output
    write_output(filtered_alignments, output_file)


def main():
    if len(sys.argv) != 3:
        print("Usage: python filter_coords.py <input_coords_file> <output_file>")
        sys.exit(1)
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    filter_coords(input_file, output_file)


if __name__ == "__main__":
    main()
