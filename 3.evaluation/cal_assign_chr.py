#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
@Author: Zijie Jiang
@Contact: jzjlab@163.com
@File: cal_assign_chr.py
@Time: 2025/6/9 17:34
@Function: calculate contig to chromosome mapping and assign sub-indices for contigs mapping to the same chromosome
"""
import sys

import numpy as np
import pandas as pd


def get_contig_chr_mapping(cluster_file, contig_scaffold_file, output_file):
    cluster_df = pd.read_csv(cluster_file, sep="\t", usecols=[0, 1], header=0, names=["chr", "scaffold", ])

    chr_info = cluster_df.set_index('scaffold')['chr'].to_dict()

    contig_df = pd.read_csv(contig_scaffold_file, sep="\t", usecols=[6, 7], header=0, names=["scaffold", "contig"])

    contig_df['assigned_chr'] = contig_df['scaffold'].map(chr_info).fillna('unassigned')

    contig_df['actual_chr'] = contig_df['contig'].str.split("_").str[0]

    contig_df['actual_chr_number'] = contig_df['actual_chr'].str.strip("chr").astype(int)

    def assign_sub_index(group):
        n = len(group)
        if n == 1:
            group['sub_index'] = group['actual_chr_number']
        else:

            step = 1.0 / (n + 1)
            sub_indices = np.linspace(-0.5, 0.5, n, endpoint=False) + step / 2
            group['sub_index'] = group['actual_chr_number'] + sub_indices
        return group

    contig_df = contig_df.groupby('actual_chr_number').apply(assign_sub_index).reset_index(drop=True)

    contig_df.to_excel(output_file, index=False)


def main():
    input = sys.argv[1]

    contig_scaffold_file = sys.argv[2]

    output_file = sys.argv[3]

    get_contig_chr_mapping(input, contig_scaffold_file, output_file)


if __name__ == '__main__':
    main()
