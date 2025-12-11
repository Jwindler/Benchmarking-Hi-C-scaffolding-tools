## Hi-C Interaction Ratio

This document describes the **Hi-C interaction ratio** calculation method and the calculation process used in this paper.



![Schematic diagram illustrating the calculation of the Hi-C interaction ratio.](https://s2.loli.net/2025/12/11/kmlBxdTDGyMfz7L.png)

The Hi-C interaction ratio is computed in three main steps. First, the whole genome interaction matrix is extracted from the input file, which can be provided in either dense (.hic) or sparse (.bed) format. This matrix serves as the basis for subsequent calculations. In the second step, the whole genome matrix is divided into submatrices according to the user-specified resolution and window size to facilitate localized counting. In the final step, two metrics are calculated using predefined formulas: the Hi-C signal ratio and the Hi-C intra/inter ratio.



### 1. Preprocessing

First, the whole-genome interaction matrix is extracted from the Hi-C results based on the specified resolution.

```sh
juicer_tool=/paht/juicer/scripts/common/juicer_tools.jar

hic_path=sample_scaffolder.hic

resolution=500000

output=sample_scaffolderr_500K.txt

java -Xmx100g -jar $juicer_tool dump observed KR $hic_path assembly assembly BP $resolution $output
```

Nots:

- `juicer_tools.jar` is from https://github.com/aidenlab/juicer

- The `.hic` file is generated from the scaffolding results; please refer to the respective scaffolders' tutorials for more information.

- The `resolution` used in this article were 100kb, 500kb, and 2.5Mb. Higher resolutions require more time to extract the matrix and calculate the Hi-C interaction ratio.  It is recommended to calculate the Hi-C interaction ratio at multiple resolutions (here use 500kb for example).

    

### 2. Hi-C signal ratio

The Hi-C signal ratio is defined as the ratio of interaction strength within the diagonal regions of the whole genome matrix to the background (off-diagonal) interaction strength. A high ratio indicates strong self-interaction of correctly assembled scaffolds. Mis-assemblies cause spurious long-range interactions in the background region, thereby lowering the ratio.

```sh
matrix=sample_scaffolderr_500K.txt

bin_size=500000

bin_windows=50

python3 cal_hic_signal_ratio.py $matrix $bin_size $bin_windows
```

Nots:

- `matrix` from the preceding preprocessing step.
- The `bin_size` used in this article were 100kb, 500kb, and 2.5Mb. It is recommended to calculate at multiple resolutions (here use 500kb for example).
- `bin_windows` is used to specify the number of bin extensions on both sides of the matrix (here use 50 for example).



### 3. Hi-C intra/inter ratio

The Hi-C intra/inter ratio is the ratio of total intra-chromosomal interactions to inter-chromosomal interactions in the whole genome matrix. Because Hi-C signals reflect spatial proximity, properly scaffolded chromosomes exhibit markedly higher intra-chromosomal than inter-chromosomal contact frequencies. Consequently, higher accuracy in clustering, ordering, and orientation of scaffolds results in a higher intra/inter ratio.

```sh
matrix=sample_scaffolderr_500K.txt

chrom_lengths_file=

bin_size=500000

python3 cal_hic_intra_intre_ratio.py $matrix $chrom_lengths_file $bin_size
```

Nots:

- `matrix` from the preceding preprocessing step.
- `chrom_lengths_file` contains the lengths of each sequence in the scaffolds (format: name /t length).
- The `bin_size` used in this article were 100kb, 500kb, and 2.5Mb. It is recommended to calculate at multiple resolutions (here use 500kb for example).



The `cal_hic_signal_ratio.py` and `cal_hic_intra_intre_ratio.py` calculation scripts are available in this directory.
