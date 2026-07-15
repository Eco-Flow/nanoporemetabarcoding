Second test FASTQ generated using:

```
fastQpick -f 0.8 -z test.fastq.gz test2.old.fastq.gz
mv fastQpick_output/test.fastq.gz test2.fastq.gz
```

Looks like fastQpick doesn't work as expected and uses the same output name for both the sampled and subsampled.
