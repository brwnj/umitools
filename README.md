#umitools

Tools to handle reads sequenced with unique molecular identifiers (UMIs).

## Trim the UMI

Incorporate the UMI into the read name in order to later identify while
processing mapped reads.

```
python umitools.py trim --end 5 unprocessed_fastq NNNNNV > out.fq
```

## Remove Duplicates

For any given start site, save only one read per UMI.

```
python umitools.py rmdup unprocessed.bam out.bam NNNNNV
```

##Requires

+ pysam
+ toolshed
