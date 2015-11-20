#umitools

[![DOI](https://zenodo.org/badge/19296/brwnj/umitools.svg)](https://zenodo.org/badge/latestdoi/19296/brwnj/umitools)

Tools to handle reads sequenced with unique molecular identifiers (UMIs).

## Trim the UMI

Incorporate the UMI into the read name in order to later identify while
processing mapped reads.

```
umitools trim --end 5 unprocessed_fastq NNNNNV > out.fq
```

## Remove Duplicates

For any given start site, save only one read per UMI. Writes bed3+ to stdout
with before and after counts per start.

```
umitools rmdup unprocessed.bam out.bam > before_after.bed
```

Specifying `--mismatches` will, for a given start site, merge all UMIs within that
edit distance into a single unique hit. For example, if a new UMI is within a single
mismatch of any existing observed UMIs for a start position, it will be merged and
considered a duplicate. The mismatch can occur at any position, regardless of the
IUPAC sequence you're using.

##Requires

```
pip install pysam editdist
```
