# Summary statistics from the uSFS data

```shell script
ls ../sfs/*_sfs_data.txt | python regional_stats.py -sel > summary_stats_usfs.csv
Rscript sum_stats_comp.R
```

![](summary_stats.png)