# Simulation with various factors

To investigate the factors that may influence scaffolder performance, we conducted a series of simulations varying contig assembly contig N50, coefficient of variation in contig length (CV), proportion of chimeric contigs, and Hi-C sequencing depth.

​	 

​	   

## Content

- [Simulation with various factors](#simulation-with-various-factors)
  - [N50](#n50)
  - [CV](#cv)
  - [Chimeric](#chimeric)
  - [Hi-C depth](#hi-c-depth)

​	   

​	  

## N50

For simulating varying contig N50 values (CV: 0.3), we fragmented the reference genome to produce assemblies with target N50 values of approximately 50 kb, 100 kb, 300 kb, 500 kb, 1 mb, and 3 mb.

```sh
smi_contigs=sim_contigs.py

genome=sample.fa

for i in 50K 100K 300K 500K 1M 3M; do
  case $i in
    *K)
      value=$(($(echo $i | sed 's/K//') * 1000))
      ;;
    *M)
      value=$(($(echo $i | sed 's/M//') * 1000 * 1000))
      ;;
  esac
  
  python $smi_contigs $genome $value 0.3 --seed 888
done
```

​	 



## CV

For simulating varying CV values, we employed the genome assembly with contig N50 of 300 kb and simulated CV values of 0.1, 0.3, 0.5, 1, 2, and 3.

```sh
smi_contigs=sim_contigs.py

genome=sample.fa

for i in 0.1 0.3 0.5 1 2 3; do
  python $smi_contigs $genome 300000 $i --seed 888  
done
```

​	 



## Chimeric

Using the assembly with a contig N50 of 300 kb and CV of 0 as the baseline, we generated a series of contig assemblies with chimeric proportions of 0.03, 0.05, 0.07, 0.1, 0.2 and 0.3. 

```sh
sim_chimeric_contigs=sim_chimeric_contigs.py

genome=sample.fa

for i in 0.1 0.3 0.5 1 2 3; do
  python $sim_chimeric_contigs $genome --chimeric_ratio $i --seed 888
done
```

​	 



## Hi-C depth

For simulating Hi-C sequencing depth, we first determined the effective Hi-C coverage (after removal of invalid pairs). Then, we random sampled the original alignments with SAMtools (version: 1.22.1) to generate datasets at 1X, 6 X, 12 X, 18 X, and 24 X effective depth.

- Calculate effective Hi-C coverage.

```sh
input=sample.bam

threads=20

samtools flagstat $input -@ $threads

# Number of read bases
samtools view $input | awk '{sum += length($10)} END {print sum}'

# sequence depth
samtools sort $input -o result.bam -@ $threads

samtools depth -@ $threads -a result.bam > depth.txt

awk '{sum+=$3} END {print sum/NR}' depth.txt
```



- Downsampling

```sh
# bam
input=sample.bam

sampled_bam=sampled.bam

sample_ratio=

samtools view $input -s $sample_ratio -o $sampled_bam --threads $threads

# 
sample_mnd=sample_mnd.py

input=merged_nodups.txt

sampled_mnd=sampled_merged_nodups.txt

read_pairs=

sample_ratio=

python3 $sample_mnd $input $read_pairs $sample_ratio --seed 888 > $sampled_bam
```



**Note**: `sim_contigs.py`, `sim_chimeric_contigs.py`, and `ample_mnd.py` are modified versions of the scripts from [HapHiC](https://github.com/zengxiaofei/HapHiC) repository.  Please refer to that repository if needed.
