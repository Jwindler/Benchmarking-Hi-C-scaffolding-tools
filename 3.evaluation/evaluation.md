# Performance evaluation

 	

​	 		

## Content


- [Performance evaluation](#performance-evaluation)
  - [Contiguity](#contiguity)
  - [Accuracy](#accuracy)
    - [Preprocessing](#preprocessing)
    - [filter](#filter)
    - [calculate accuracy](#calculate-accuracy)
  - [Chromosome assignment](#chromosome-assignment)

​	 

​	   

## Contiguity

The contiguity metrics (scaffold number, total assembly size, NG50, LG50, and mis-assembly count) mentioned in this study were evaluated using `QUAST` (version: 5.3.0).

```sh
scaffolds=sample.fa
# scaffold-level genome

reference=T2T-CHM13.fa
# reference T2T genome

threads=40

output=results

quast.py $scaffolds -r $reference --threads $threads -o $output
```

​	 

 	

## Accuracy

The process for calculating the accuracy of scaffolding using simulation data is as follows:

Briefly, scaffolds were aligned to both the T2T reference genome and the original contig-level assembly using nucmer (MUMmer4, version: 4.0.1). The resulting .delta files were converted to tabular coordinates with show-coords. Then filtered alignments and computed statistics. 

​	 

​	 

### Preprocessing

```sh
threads=40

contig=sample_contig.fa

scaffold=sample_scaffold.fa

reference=T2T-CHM13.fa

# conitg compare scaffold
nucmer --threads $threads -p conitg_scaffold $scaffold $contig

show-coords -rcl conitg_scaffold.delta > conitg_scaffold.coords

# scaffold compare reference
nucmer --threads $threads -p scaffold_ref $reference $scaffold

show-coords -rcl scaffold_ref.delta > scaffold_ref.coords
```

​	 

​	 

### filter

```sh
coords_file=conitg_scaffold.coords

filtered_coords_file=conitg_scaffold_filtered.txt

python3 filter_coords.py $coords_file $filtered_coords_file
```

​	 

​	 

### calculate accuracy

```sh
 # genome fasta file
genome=sample_scaffold.fa

# contig compare scaffold file
contig_scaffold=conitg_scaffold_filtered.txt

# scaffold compare reference file
scaffold_ref=scaffold_reference_filtered.txt

# output directory
output_dir=./

python3 cal_scaffold_accuracy.py $reference_genome $contig_scaffold $scaffold_ref $output_dir
```

- output：

    - cluster.txt: the accuracy rate of the successful cluster in the scaffold.
    - order.txt: the accuracy rate of the successful order in the scaffold.
    - orientation.txt: the accuracy rate of the successful orientation in the scaffold.


​	 

​	   

## Chromosome assignment

![](https://s2.loli.net/2025/12/11/3jPQe5brTS8dosD.png)

Each dot represents a scaffold sequence, and the scaffolds with different colours indicate their original chromosome assignment in the genome. Unanchored scaffolds are placed above each figure. A dot off-diagonal represents incorrectly anchored scaffolds.

```sh
input=cluster.txt

contig_scaffold_file=conitg_scaffold_filtered.txt

output_file=sample.xlsx

python3 cal_assign_chr.py $input $contig_scaffold_file $output_file
```

- output：The chromosome number of each scaffold sequence corresponds to the actual chromosome number.
