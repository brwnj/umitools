#umitools

Tools to handle reads sequenced with unique molecular identifiers (UMIs).

Right now, this only handles UMIs at the 3' end ([READ][UMI]). Support will be 
added soon for 5' UMIs.

## Trim the UMI

Incorporate the UMI into the read sequence in order to later identify among mapped reads.
```
umitools trim unprocessed_fastq NNNNNV > out.fq
```

## Remove Duplicates

For any given start site, save only one read per UMI.
```
umitools rmdup unprocess.bam out.bam
```

#Requirements

+ pysam
+ toolshed