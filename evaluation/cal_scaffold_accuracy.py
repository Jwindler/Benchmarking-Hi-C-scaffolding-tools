#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
@Author: Zijie Jiang
@Contact: jzjlab@163.com
@File: cal_scaffold_accuracy.py
@Time: 2025/3/26 14:58
@Function: calculate scaffold accuracy
"""
import os.path
import sys
from collections import defaultdict

import pandas as pd
from Bio import SeqIO


def read_genome(genome):
    genome_info = {}
    for record in SeqIO.parse(genome, "fasta"):
        genome_info[record.id] = {
            "length": len(record.seq)
        }

    return genome_info


def cal_cluster_rate(genome, scaffold_info, contig2scaffold, output_file):
    genome_info = read_genome(genome)

    cluster_info = get_cluster_info(scaffold_info, contig2scaffold)

    for chrom in genome_info:
        cluster_info[chrom]["cluster_rate"] = cluster_info[chrom]["cluster_length"] / genome_info[chrom]["length"]

    with open(output_file, 'w') as f:
        f.write(
            "Chr\tScaffold\tCluster_rate\tCluster_count\tCluster_len\tChr_len\n")

        for chrom in genome_info:
            temp_rate = cluster_info[chrom]["cluster_length"] / genome_info[chrom]["length"]
            cluster_info[chrom]["cluster_rate"] = temp_rate

            line = f'{chrom}\t{cluster_info[chrom]["scaffold"]}\t{temp_rate}\t{cluster_info[chrom]["cluster_count"]}\t{cluster_info[chrom]["cluster_length"]}\t{genome_info[chrom]["length"]}\n'
            f.write(line)

    return cluster_info


def get_cluster_info(scaffold_info, contig2scaffold):
    cluster_info = {}

    contig2scaffold_df = pd.read_csv(contig2scaffold, sep='\t')  # 制表符分隔
    contig2scaffold_df.columns = ["ref_start", "ref_end", "query_start", "query_end", "ref_len", "query_len",
                                  "ref_name", "query_name", "direction", "align_len", "identity"]

    for key in scaffold_info:
        max_query_name = scaffold_info[key]['max_query_name']

        temp_ref = contig2scaffold_df[contig2scaffold_df['ref_name'] == max_query_name]

        unique_rows = len(temp_ref.drop_duplicates(subset='query_name'))

        mask = temp_ref['query_name'].str.startswith(key + '_')

        temp_len = temp_ref.loc[mask, 'align_len'].sum()

        cluster_info[key] = {
            "cluster_length": temp_len,
            "cluster_count": unique_rows,
            "scaffold": max_query_name,
            "chromosome": key
        }

    return cluster_info


def cal_order_rate(genome, scaffold_info, scaffold2ref, output_file):
    scaffold2ref_df = pd.read_csv(scaffold2ref, sep='\t')  # 制表符分隔
    scaffold2ref_df.columns = ["ref_start", "ref_end", "query_start", "query_end", "ref_len", "query_len",
                               "ref_name", "query_name", "direction", "align_len", "identity"]

    order_info = {}
    genome_info = read_genome(genome)

    with open(output_file, 'w') as f:
        f.write(
            "Chr\tScaffold\tOrder_rate\tOrder_len\tChr_len\n")
        for key in genome_info:
            if key == 'chr1':
                print("DEBUG")
            max_query_name = scaffold_info[key]['max_query_name']

            temp_ref = scaffold2ref_df[
                (scaffold2ref_df['query_name'] == max_query_name) & (scaffold2ref_df['ref_name'] == key)]

            temp_len = temp_ref['align_len'].sum()

            temp_rate = temp_len / genome_info[key]['length']

            line = f'{key}\t{max_query_name}\t{temp_rate}\t{temp_len}\t{genome_info[key]["length"]}\n'
            order_info[key] = {
                "order_rate": temp_rate,
                "order_len": temp_len,
                "scaffold": max_query_name,
                "chromosome": key,
                "chromosome_len": genome_info[key]["length"]

            }
            f.write(line)

    return order_info


def cal_orientation_rate(genome, scaffold_info, scaffold2ref, output_file):
    scaffold2ref_df = pd.read_csv(scaffold2ref, sep='\t')  # 制表符分隔
    scaffold2ref_df.columns = ["ref_start", "ref_end", "query_start", "query_end", "ref_len", "query_len",
                               "ref_name", "query_name", "direction", "align_len", "identity"]

    orientation_info = {}
    genome_info = read_genome(genome)

    with open(output_file, 'w') as f:
        f.write(
            "Chr\tScaffold\tOrientation_rate\tOrientation\tOrientation_len\tForward_count\tForward_len\tReverse_count\tReverse_len\tChr_len\n")
        for key in genome_info:
            max_query_name = scaffold_info[key]['max_query_name']

            temp_ref = scaffold2ref_df[
                (scaffold2ref_df['query_name'] == max_query_name) & (scaffold2ref_df['ref_name'] == key)]

            temp_forward_df = temp_ref[temp_ref['direction'] == 'Forward']
            temp_reverse_df = temp_ref[temp_ref['direction'] == 'Reverse']

            temp_forward_df_len = temp_forward_df['align_len'].sum()
            temp_reverse_df_len = temp_reverse_df['align_len'].sum()

            if temp_forward_df_len >= temp_reverse_df_len:
                temp_rate = temp_forward_df_len / genome_info[key]['length']
                orientation = 'Forward'
                orientation_len = temp_forward_df_len
            else:
                temp_rate = temp_reverse_df_len / genome_info[key]['length']
                orientation = 'Reverse'
                orientation_len = temp_reverse_df_len

            line = f'{key}\t{max_query_name}\t{temp_rate}\t{orientation}\t{orientation_len}\t{len(temp_forward_df)}\t{temp_forward_df_len}\t{len(temp_reverse_df)}\t{temp_reverse_df_len}\t{genome_info[key]["length"]}\n'
            orientation_info[key] = {
                "orientation_rate": temp_rate,
                "scaffold": max_query_name,
                "chromosome": key,
                "orientation": orientation,
                "orientation_len": orientation_len,
                "chromosome_len": genome_info[key]["length"]

            }
            f.write(line)


def parse_filtered_coords(scaffold2ref):
    """
        parse filtered coords file to get scaffold info
    Args:
        scaffold2ref: filtered scaffold to reference coords file path

    Returns:

    """

    ref_stats = defaultdict(lambda: {
        'count': 0,
        'total_length': 0,
        'query_lengths': defaultdict(int)
    })

    try:
        with open(scaffold2ref, 'r') as f:

            lines = f.readlines()

            for line in lines[1:]:
                fields = line.split('\t')
                ref_name = fields[6]
                query_name = fields[7]
                align_len = int(fields[9])

                # 更新统计数据
                ref_stats[ref_name]['count'] += 1
                ref_stats[ref_name]['total_length'] += align_len
                ref_stats[ref_name]['query_lengths'][query_name] += align_len

    except FileNotFoundError:
        print(f"错误：文件 {scaffold2ref} 未找到")
        return None
    except Exception as e:
        print(f"解析文件时出错：{str(e)}")
        return None

    scaffold_info = {}

    for ref_name, stats in sorted(ref_stats.items()):
        max_query = max(stats['query_lengths'].items(), key=lambda x: x[1])
        max_query_name = max_query[0]
        max_query_length = max_query[1]

        # record scaffold info
        scaffold_info[ref_name] = {
            'count': stats['count'],
            'total_length': stats['total_length'],
            'max_query_name': max_query_name,
            'max_query_length': max_query_length
        }

    return scaffold_info


def main():
    genome = sys.argv[1]  # genome fasta file

    contig_scaffold = sys.argv[2]  # contig compare scaffold file

    scaffold_ref = sys.argv[3]  # scaffold compare reference file

    output_dir = sys.argv[4]  # output directory

    cluster_output_file = os.path.join(output_dir, "cluster.txt")
    order_output_file = os.path.join(output_dir, "order.txt")
    orientation_output_file = os.path.join(output_dir, "orientation.txt")

    scaffold_info = parse_filtered_coords(scaffold_ref)

    cal_cluster_rate(genome, scaffold_info, contig_scaffold, cluster_output_file)

    cal_order_rate(genome, scaffold_info, scaffold_ref, order_output_file)

    cal_orientation_rate(genome, scaffold_info, scaffold_ref, orientation_output_file)


if __name__ == '__main__':
    main()
