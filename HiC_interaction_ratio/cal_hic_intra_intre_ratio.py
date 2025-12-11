#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
@Author: Zijie Jiang
@Contact: jzjlab@163.com
@File: cal_hic_intra_intre_ratio.py
@Time: 2025/8/4 21:47
@Function: calculate intra- and inter-scaffold contact ratios from Hi-C contact matrix
"""

import sys

import numpy as np
import pandas as pd


# read chromosome/scaffold lengths
def load_scaffold_lengths(file_path):
    """
        Read scaffold lengths from a file.
    Args:
        file_path: Path to the scaffold lengths file.

    Returns:
        A dictionary mapping scaffold names to their lengths.
    """

    scaffolds = {}
    with open(file_path, 'r') as f:
        for line in f:
            name, length = line.strip().split()
            scaffolds[name] = int(length)
    return scaffolds


def calculate_inter_scaffold_ratio(matrix_file, chrom_lengths_file, bin_size=10000):
    """
        Calculate intra- and inter-scaffold contact counts from a Hi-C contact matrix.
    Args:
        matrix_file: Path to the Hi-C contact matrix file.
        chrom_lengths_file: Path to the scaffold lengths file.
        bin_size: Size of each bin in the contact matrix.

    Returns:
        A tuple containing intra-scaffold contact count, inter-scaffold contact count, and total contact count.
    """

    scaffolds = load_scaffold_lengths(chrom_lengths_file)

    data = pd.read_csv(matrix_file, sep='\s+', names=['bin1', 'bin2', 'count'])

    bin_to_scaffold = {}
    current_bin = 0
    for scaffold, length in scaffolds.items():
        num_bins = (length // bin_size) + 1
        for i in range(current_bin, current_bin + num_bins):
            bin_to_scaffold[i] = scaffold
        current_bin += num_bins

    intra_contact = 0
    inter_contact = 0
    total_contact = 0

    for _, row in data.iterrows():
        bin1, bin2, count = int(row['bin1'] // bin_size), int(row['bin2'] // bin_size), row['count']
        if np.isnan(count):
            continue
        total_contact += count
        if bin1 in bin_to_scaffold and bin2 in bin_to_scaffold:
            if bin_to_scaffold[bin1] == bin_to_scaffold[bin2]:
                intra_contact += count
            else:
                inter_contact += count

    return intra_contact, inter_contact, total_contact


def main():
    if len(sys.argv) != 4:
        print("Usage: python cal_hic_intra_intre_ratio.py <matrix_file> <chrom_lengths_file> <bin_size>\n")
        print("Recommended bin_size: 500000")
        sys.exit(1)

    matrix_file = sys.argv[1]
    print(f"Processing file: {matrix_file}\n")

    chrom_lengths_file = sys.argv[2]
    print(f"Using chromosome lengths file: {chrom_lengths_file}\n")

    bin_size = int(sys.argv[3])
    print(f"Using bin size: {bin_size}\n")

    intra_contact, inter_contact, total_contact = calculate_inter_scaffold_ratio(matrix_file, chrom_lengths_file,
                                                                                 bin_size)
    inter_ratio = inter_contact / total_contact if total_contact > 0 else 0
    intra_ratio = intra_contact / total_contact if total_contact > 0 else 0
    intra_inter_ratio = intra_contact / inter_contact if inter_contact > 0 else 0
    print(
        f"Intra-scaffold Contacts\t\tInter-scaffold Contacts\tTotal Contacts\tIntra Ratio\tInter Ratio\tIntra/Inter Ratio")

    print(
        f"{int(intra_contact)}\t{int(inter_contact)}\t{int(total_contact)}\t{intra_ratio:.2f}\t{inter_ratio:.2f}\t{intra_inter_ratio:.2f}\n")

    print(f"Finished processing file: {matrix_file}\n")


if __name__ == "__main__":
    main()
