#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
@Author: Zijie Jiang
@Contact: jzjlab@163.com
@File: asy2chr_len.py
@Time: 2025/8/17 15:32
@Function: convert .assembly file to chr len file
"""


def main():
    """
        convert .assembly file to chr len file
    Returns:
        None
    """

    asy_file = "sample.assembly"
    output_file = "sample_chr_len.txt"

    scaffold_info = {}
    scaffold_index = {}
    scaffold_index_counter = 1

    with open(asy_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                scaffold_name, index, scaffold_len = line.strip().split(" ")
                scaffold_info[index] = {
                    "name": scaffold_name[1:],  # Remove '>' from the scaffold name
                    "length": int(scaffold_len)
                }
            else:
                scaffold_index[scaffold_index_counter] = line.strip()
                scaffold_index_counter += 1
    scaffold_info_new = {}
    counter = 1
    for key, value in scaffold_index.items():
        temp_value = value.split(" ")

        for i in temp_value:
            temp_i = i.replace("-", "")

            print(temp_i)
            scaffold_info_new[counter] = scaffold_info[temp_i]["length"]

            counter += 1

    with open(output_file, 'w') as f:
        for key, value in scaffold_info_new.items():
            line = f"{key}\t{value}\n"
            f.write(line)


if __name__ == '__main__':
    main()
