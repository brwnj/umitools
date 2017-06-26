# UMI Tools

[![DOI](https://zenodo.org/badge/19296/brwnj/umitools.svg)](https://zenodo.org/badge/latestdoi/19296/brwnj/umitools)

Tools to handle reads sequenced with unique molecular identifiers (UMIs).

## Trim the UMI

Incorporate the UMI into the read name in order to later identify while
processing mapped reads.

```
umitools trim --end 5 unprocessed_fastq NNNNNV > out.fq
```

If you want to save reads with invalid UMI sequences, you can specify `--invalid`.

```
umitools trim --end 5 --invalid bad_umi.fq unprocessed_fastq NNNNNV > out.fq
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

## Installation

umitools has two requirements: [pysam][] and [editdist][].
Use pip to install [pysam].

```
pip install pysam
```

[editdist] has to be downloaded and installed from source ([Downloads page][editdist-download]).

```
wget https://py-editdist.googlecode.com/files/py-editdist-0.3.tar.gz
tar xzf py-editdist-0.3.tar.gz
cd py-editdist-0.3/
python setup.py install
```

Finally download and install umitools from source.

```
wget -O umitools-master.zip https://github.com/brwnj/umitools/archive/master.zip
unzip umitools-master.zip
cd umitools-master
python setup.py install
```

[pysam]: https://pypi.python.org/pypi/pysam
[editdist]: https://pypi.python.org/pypi/editdist/0.1
[editdist-download]: https://code.google.com/p/py-editdist/downloads/list
