# Hi-C-based scaffolding

The configuration files for building the scaffolder runtime environment are already in the current directory. The preparation process for all scaffolder input files and the commands for running them are as follows:



## Content

- [Hi-C-based scaffolding](#hi-c-based-scaffolding)
  - [Content](#content)
  - [Preprocessing](#preprocessing)
    - [Juicer](#juicer)
    - [Bam](#bam)
  - [3D-DNA](#3d-dna)
  - [HapHiC](#haphic)
  - [Pin\_hic](#pin_hic)
  - [SALSA2](#salsa2)
  - [YaHS](#yahs)




## Preprocessing



### Juicer

The input file for `3D-DNA` is <merged_nodups.txt>, which is constructed using the `Juicer` workflow, as follows:

```sh
# Align Hi-C data to the assembly to get <merged_nodups.txt>
bash juicer.sh -z <contig.fa> -p <contig.size> -y <restriction.site.file> -s <RE> -d output -D <juicer_dir> -S early 

```



### Bam

For other scaffolders, use `<hic.filtered.bam>` or other similar format, and the build command is as follows:

```sh
bwa index <contig.fa>

# Align Hi-C data to the assembly, remove PCR duplicates and filter out secondary and supplementary alignments
bwa mem -5SP <contig.fa> $fq1 $fq2 | samblaster | samtools view - -S -h -b -F 3340 -o hic.bam

# Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
filter_bam hic.bam 1 --nm 3 | samtools view - -b -o hic.filtered.bam
```





## 3D-DNA

```sh
bash run-asm-pipeline.sh <contig.fa> <merged_nodups.txt>
```



## HapHiC

```sh
haphic pipeline <contig.fa> <hic.filtered.bam> <nchrs> --correct_nrounds 2
```



## Pin_hic

```sh
pin_hic_it -x < contig.fa.fai> -r < contig.fa> -O output <hic.filtered.bam>
```



## SALSA2

```sh
# convert bam format to bed format
bamToBed -i <hic.filtered.bam> > <hic.bed>

# sort bed
sort -k 4 <hic.bed> > tmp && mv tmp <hic.bed>

# run SALSA2
python run_pipeline.py -a < contig.fa> -l < contig.fa.fai> -b <hic.bed> -e <RE> -m yes -o output
```



## YaHS

```sh
yahs <contig.fa> <hic.filtered.bam> -e <RE>
```
