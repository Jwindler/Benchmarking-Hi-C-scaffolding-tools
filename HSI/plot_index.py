#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
@Author: Zijie Jiang
@Contact: jzjlab@163.com
@File: plot_index.py
@Time: 2026/8/24 15:53
@Function:
"""
import os

import matplotlib.pyplot as plt
import pandas as pd


def plot_index(input, outdir, inner=False):
	os.makedirs(outdir, exist_ok=True)

	df = pd.read_csv(input, sep=r"\s+")

	for scaffold, group_data in df.groupby("Scaffold"):
		group_data = group_data.sort_values(by="region")

		region_mb = group_data["region"] / 1e6

		if inner:
			region_mb = region_mb - region_mb.iat[0]

		plt.figure(figsize=(8, 5))

		plt.plot(
			region_mb,
			group_data["index"],
			marker="o",
			linestyle="-",
			linewidth=1.5,
			markersize=.7,
			label=scaffold,
		)

		plt.title(f"{scaffold}", fontsize=12)
		plt.xlabel("Region (Mb)", fontsize=10)
		plt.ylabel("Scaffold inter index", fontsize=10)
		if inner:
			plt.ylabel("Interaction score", fontsize=10)

		plt.grid(True, linestyle="--", alpha=0.6)

		plt.savefig(os.path.join(outdir, f"{scaffold}.png"), dpi=300, bbox_inches="tight")

		plt.close()

	print(f"Processing complete! All plots have been saved to: {outdir}\n")


def main():
	input_file = "/home/jzj/downloads/Fig6-example/Ok/500k/Pin_hic-inter-index.txt"

	output_dir = "/home/jzj/downloads/output_plots"

	plot_index(input_file, output_dir)


if __name__ == '__main__':
	main()
