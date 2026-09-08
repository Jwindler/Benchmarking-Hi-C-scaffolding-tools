#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
@Author: Zijie Jiang
@Contact: jzjlab@163.com
@File: parse_contig_coords.py
@Time: 2025/3/26 10:32
@Function: main program entry
"""

import pysam


def parse_fasta(fasta):
    fa_dict = dict()

    with open(fasta) as f:

        for line in f:
            if line.startswith('>'):
                temp_id = line.split()[0][1:]
                fa_dict[temp_id] = [0, 0, 0]

    return fa_dict


def parse_bam(bam, fa_dict, threads):
    with pysam.AlignmentFile(bam, threads=threads, format_options=[b'filter=flag.read1 && refid != mrefid']) as f:

        for aln in f:
            reference_name, next_reference_name = aln.reference_name, aln.next_reference_name
            if next_reference_name == reference_name:
                fa_dict[reference_name][0] += 1
            elif reference_name.split('_')[0] == next_reference_name.split('_')[0]:
                fa_dict[reference_name][1] += 1
                fa_dict[next_reference_name][1] += 1
            else:
                fa_dict[reference_name][2] += 1
                fa_dict[next_reference_name][2] += 1


def output_statistics(fa_dict, tag):
    with open('{}_HiC_links.txt'.format(tag), 'w') as f:
        for ctg, link_list in fa_dict.items():
            f.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                ctg, link_list[0], link_list[1], link_list[2], sum(link_list), tag))


def main():
    bam = "sample.bam"
    fasta = "sample.fa"
    tag = "sample"
    threads = 4

    fa_dict = parse_fasta(fasta)

    parse_bam(bam, fa_dict, threads)

    output_statistics(fa_dict, tag)


if __name__ == '__main__':
    main()
