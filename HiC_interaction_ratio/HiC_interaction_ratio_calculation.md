

## Preprocessing

```sh

juicer_tool=

hic_path=

resolution=

output=

java -Xmx100g -jar $juicer_tool dump observed KR $hic_path assembly assembly BP $resolution $output

```





## Hi-C signal ratio



```sh

matrix=

bin_size=

bin_windows=

python3 cal_hic_signal_ratio.py $matrix $bin_size $bin_windows
```





## Hi-C intra/inter ratio



```sh

matrix=

chrom_lengths_file=

bin_size=

python3 cal_hic_intra_intre_ratio.py $matrix $chrom_lengths_file $bin_size
```

