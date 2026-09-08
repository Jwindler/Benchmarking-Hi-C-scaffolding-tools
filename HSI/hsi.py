#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
@Author: Zijie Jiang
@Contact: jzjlab@163.com
@File: hsi.py
@Time: 2026/8/18 10:02
@Function:
"""
import os

import pandas as pd


def get_len(file_path):
	len_dict = {}
	with open(file_path, 'r') as f:
		for line in f:
			if line.startswith('#') or not line.strip():
				continue
			chrom = line.strip().split()[0]
			length = line.strip().split()[1]
			len_dict[chrom] = int(length)
	return len_dict


def cal_index(matrix_file, len_file, bin_size, index_output, inter_index_output,
              inner_index_output):
	scaffolds_len = get_len(len_file)

	genome_len = sum(scaffolds_len.values())

	matrix: pd.DataFrame = pd.read_csv(matrix_file, sep=r'\s+',  # type: ignore
	                                   names=['bin1', 'bin2', 'count'])

	matrix['count'] = pd.to_numeric(matrix['count'], errors='coerce').fillna(0)

	matrix_max_len = matrix.iat[-1, 1]

	len_ratio = genome_len // matrix_max_len

	print(f"Genome size: {genome_len}, Matrix max length: {matrix_max_len}, ratio: {len_ratio}\n")

	cols = ['bin1', 'bin2']
	scaled_matrix = matrix.copy()
	scaled_matrix[cols] = (scaled_matrix[cols] * len_ratio).round().astype(int)

	scaffolds_index = {}

	inter_list = []
	inner_list = []

	last_len = 0
	for scaffold, length in scaffolds_len.items():
		start = last_len
		end = last_len + length

		last_len += length

		if length < bin_size:  # TODO: filter out scaffolds smaller than bin size
			continue

		inner = scaled_matrix.loc[
			(scaled_matrix['bin1'] >= start) & (scaled_matrix['bin1'] < end) & (
					scaled_matrix['bin2'] >= start) & (
					scaled_matrix['bin2'] < end), 'count'].sum()

		inter = scaled_matrix.loc[
			(scaled_matrix['bin1'] >= start) & (scaled_matrix['bin1'] < end) &
			(scaled_matrix['bin2'] >= end), 'count'].sum()

		index_val = inner / inter if inter != 0 else 0.0

		scaffolds_index[scaffold] = {'length': length, 'inner': round(inner),
		                             'inter': round(inter), 'index': index_val}

		# scaffolds inter index
		inter_matrix = scaled_matrix.loc[
			(scaled_matrix['bin1'] >= start) & (scaled_matrix['bin1'] < end) &
			(scaled_matrix['bin2'] >= end)]

		inter_sum = inter_matrix.groupby('bin2', as_index=False)['count'].sum()

		inter_sum['Scaffold'] = scaffold
		inter_sum['index'] = inter_sum['count'] / inner

		inter_list.append(inter_sum)

		# scaffolds inner index
		cond = (
				(scaled_matrix['bin1'] >= start) & (scaled_matrix['bin1'] < end) &
				(scaled_matrix['bin2'] >= start) & (scaled_matrix['bin2'] < end)
		)
		sub_df = scaled_matrix[cond]

		sum_bin1 = sub_df.groupby('bin1')['count'].sum()
		sum_bin2 = sub_df.groupby('bin2')['count'].sum()
		diag = sub_df[sub_df['bin1'] == sub_df['bin2']].set_index('bin1')['count']

		inner_contact = sum_bin1.add(sum_bin2, fill_value=0).sub(diag, fill_value=0)

		col_sums = pd.DataFrame({
			'region': inner_contact.index,
			'inner-contact': inner_contact.values,
			'Scaffold': scaffold,
			'index': inter / inner_contact.values
			# FIXME: RuntimeWarning: divide by zero encountered in divide
		})

		inner_list.append(col_sums)

	total_inner = sum(item['inner'] for item in scaffolds_index.values())
	print(f"Sum: {total_inner}\n")

	total_inter = sum(item['inter'] for item in scaffolds_index.values())
	print(f"Sum: {total_inter}\n")

	index = total_inner / total_inter if total_inner != 0 else float('inf')

	print(f"Index (inner/inter): {index}\n")

	# Index output
	with open(index_output, 'w') as f:
		f.write("Scaffold\tLength\tInner\tInter\tIndex\n")
		f.write(f"Genome\t{genome_len}\t{total_inner}\t{total_inter}\t{index}\n")
		for scaffold, value in scaffolds_index.items():
			f.write(
				f"{scaffold}\t{value['length']}\t{value['inner']}\t{value['inter']}\t{value['index']}\n")

	# Region index output
	each_inter_contact = pd.concat(inter_list, ignore_index=True)

	each_inter_contact = each_inter_contact.rename(
		columns={'count': 'inter-contact', 'bin2': 'region'})

	cols_order = ['Scaffold', 'region', 'inter-contact', 'index']

	each_inter_contact.to_csv(
		inter_index_output,
		sep='\t',
		columns=cols_order,
		index=False,
		header=True
	)

	# Inner region index output
	each_inner_contact = pd.concat(inner_list, ignore_index=True)

	cols_order = ['Scaffold', 'region', 'inner-contact', 'index']

	each_inner_contact.to_csv(
		inner_index_output,
		sep='\t',
		columns=cols_order,
		index=False,
		header=True
	)


def main():
	if len(sys.argv) != 4:
		print(
			"Usage: python hsi.py <matrix_file> <chrom_lengths_file> <bin_size> <outdir>\n")
		print("Recommended bin_size: 500000")
		sys.exit(1)

	matrix_file = sys.argv[1]
	print(f"Processing file: {matrix_file}\n")

	chrom_lengths_file = sys.argv[2]
	print(f"Using chromosome lengths file: {chrom_lengths_file}\n")

	bin_size = int(sys.argv[3])
	print(f"Using bin size: {bin_size}\n")

	outdir = sys.argv[4]
	print(f"Using output directory: {outdir}\n")

	index_output = os.path.join(outdir, f"{sample}-index.txt")

	inter_index_output = os.path.join(outdir, f"{sample}-inter-index.txt")

	inner_index_output = os.path.join(outdir, f"{sample}-inner-index.txt")

	cal_index(matrix_file, chrom_lengths_file, bin_size, index_output, inter_index_output,
	          inner_index_output)

	print(f"Finished processing file: {matrix_file}\n")


if __name__ == "__main__":
	main()
