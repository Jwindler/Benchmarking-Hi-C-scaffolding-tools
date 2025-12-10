

## Contiguity







## Accuracy

### Preprocessing

```sh

threads=40

contig=/lustre/home/acct-clswcc/jiangzijie/2.projects/3.ScaffoldBenchmark/1.data/3.Fragaria_ananassa/Fa_300k_0.3.fa

scaffold=/lustre/home/acct-clswcc/jiangzijie/2.projects/3.ScaffoldBenchmark/4.eva_scaffold/2.Fa/1X/1.3D-DNA/Fa_1X.3d-dna.fasta

ref=/lustre/home/acct-clswcc/jiangzijie/2.projects/3.ScaffoldBenchmark/1.data/3.Fragaria_ananassa/Fa_complete_rename.fa

# conitg_scaffold
nucmer --threads $threads -p conitg_scaffold $scaffold $contig

show-coords -rcl conitg_scaffold.delta > conitg_scaffold.coords

# scaffold_ref
nucmer --threads $threads -p scaffold_ref $ref $scaffold

show-coords -rcl scaffold_ref.delta > scaffold_ref.coords

```



### filter



```sh

coords_file=

filtered_coords_file=

python filter_coords.py $coords_file $filtered_coords_file

```



### cal

