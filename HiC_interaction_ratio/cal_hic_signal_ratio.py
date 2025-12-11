#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
@Author: Zijie Jiang
@Contact: jzjlab@163.com
@File: cal_hic_signal_ratio.py
@Time: 2025/8/4 21:18
@Function: calculate Hi-C signal-to-noise ratio from contact matrix
"""

import sys

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix


def load_contact_matrix_sparse(file_path, bin_size=10000):
    """
        Load Hi-C contact matrix from a dumped file into a sparse matrix.
    Args:
        file_path: Path to the dumped matrix file.
        bin_size: Size of each bin in the contact matrix.

    Returns:
        A scipy sparse matrix representing the Hi-C contact matrix.
    """

    # Read the dumped matrix file
    data = pd.read_csv(file_path, sep='\s+', names=['bin1', 'bin2', 'count'])

    # Convert counts to integers, replacing NaN with 0
    data['count'] = data['count'].fillna(0).astype(np.int32)

    max_bin = max(data['bin1'].max(), data['bin2'].max())

    matrix_size = (max_bin // bin_size) + 1

    row = data['bin1'] // bin_size
    col = data['bin2'] // bin_size

    values = data['count']

    # Use int32 for sparse matrix
    matrix = coo_matrix((values, (row, col)), shape=(matrix_size, matrix_size), dtype=np.int32)
    return matrix


def calculate_diagonal_intensity_sparse(matrix, k=5):
    """
        Calculate the average intensity of the diagonal region within k bins.
    Args:
        matrix: Hi-C contact sparse matrix.
        k: Number of bins around the diagonal to consider.

    Returns:
        Average intensity of the diagonal region.
    """

    diagonal_sum = 0
    diagonal_count = 0
    for i, j, v in zip(matrix.row, matrix.col, matrix.data):
        if abs(i - j) <= k:
            diagonal_sum += v
            diagonal_count += 1
    return diagonal_sum / diagonal_count if diagonal_count > 0 else 0


def calculate_background_intensity_sparse(matrix, k=5):
    """
        Calculate the average intensity of the background region outside k bins from the diagonal.
    Args:
        matrix: Hi-C contact sparse matrix.
        k: Number of bins around the diagonal to exclude.

    Returns:
        Average intensity of the background region.
    """

    background_sum = 0
    background_count = 0
    for i, j, v in zip(matrix.row, matrix.col, matrix.data):
        if abs(i - j) > k:
            background_sum += v
            background_count += 1
    return background_sum / background_count if background_count > 0 else 0


def main():
    if len(sys.argv) != 4:
        print("Usage: python cal_hic_signal_ratio.py <matrix_file> <bin_size> <bin_windows>\n")
        print("Recommended bin_size: 500000, bin_windows: 50")
        sys.exit(1)
    matrix_file = sys.argv[1]
    print(f"Processing file: {matrix_file}\n")

    bin_size = int(sys.argv[2])
    print(f"Using bin size: {bin_size}\n")

    bin_windows = int(sys.argv[3])
    print(f"Using bin windows: {bin_windows}\n")

    contact_matrix = load_contact_matrix_sparse(matrix_file, bin_size)

    diag_intensity = int(calculate_diagonal_intensity_sparse(contact_matrix, bin_windows))

    bg_intensity = int(calculate_background_intensity_sparse(contact_matrix, bin_windows))

    snr = diag_intensity / bg_intensity if bg_intensity > 0 else float('inf')

    print(f"Diagonal Intensity\tBackground Intensity\tSNR")
    print(f"{diag_intensity}\t{bg_intensity}\t{snr:.2f}\n")
    print(f"Finished processing file: {matrix_file}")


if __name__ == "__main__":
    main()
