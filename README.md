#umitools

Tools to handle reads sequenced with unique molecular identifiers (UMIs).

Right now, this only handles UMIs at the 5' end ([READ][UMI]).

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