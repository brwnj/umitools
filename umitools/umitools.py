#!/usr/bin/env python
# encoding: utf-8
"""
Tools to handle reads sequenced with unique molecular identifiers (UMIs).
"""
from __future__ import print_function

import editdistance
import gzip
import os
import re
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import Counter, defaultdict
from itertools import islice, izip
from pysam import Samfile

try:
    from itertools import izip as zip
except ImportError:
    pass

from ._version import __version__


IUPAC = {
    "A": "A",
    "T": "T",
    "C": "C",
    "G": "G",
    "R": "GA",
    "Y": "TC",
    "M": "AC",
    "K": "GT",
    "S": "GC",
    "W": "AT",
    "H": "ACT",
    "B": "GTC",
    "V": "GCA",
    "D": "GAT",
    "N": "GATC"
}
UMI_REGEX = re.compile(r'UMI_([\w]*)')
gzopen = lambda f: gzip.open(f) if f.endswith(".gz") else open(f)


class Fastq(object):
    """FASTQ record. Holds record of name, sequence, and quality scores.
    """
    def __init__(self, args):
        self.name = args[0][1:]
        self.seq = args[1]
        self.qual = args[3]
        assert len(self.seq) == len(self.qual)

    def __repr__(self):
        return "Fastq(%s)" % self.name

    def __str__(self):
        return "@%s\n%s\n+\n%s" % (self.name, self.seq, self.qual)


def umi_from_name(name):
    """Extracts the UMI sequence from the read name.

    Args:
        name (str): Name of the sequence

    Returns:
        str: UMI sequence

    >>> umi_from_name("cluster_1017333:UMI_GCCGCA")
    'GCCGCA'
    """
    return UMI_REGEX.findall(name)[0].strip()


def passing_distances(query, targets, n):
    """Tests target set of sequences to the query.

    Args:
        query (str): query sequence
        targets (set): unique sequences
        n (int): allowable mismatches when comparing a query to a given sequence of the targets

    Returns:
        bool

    >>> s = "ACTGA"
    >>> ts_1 = {"ACTGG"}
    >>> ts_2 = {"ACTCC", "ACTGG"}
    >>> ts_3 = {"ACTCC", "ACTTT"}
    >>> n = 1
    >>> passing_distances(s, ts_1, n)
    True
    >>> passing_distances(s, ts_2, n)
    True
    >>> passing_distances(s, ts_3, n)
    False
    """
    for target in targets:
        if editdistance.distance(target, query) <= n:
            return True
    return False


def process_bam(abam, bbam, mismatches=0):
    """Removes duplicate reads characterized by their UMI at any given start location.

    Args:
        abam (str): Input bam with potential duplicate UMIs
        bbam (str): Output bam after removing duplicate UMIs
        mismatches (Optional[int]): Allowable edit distance between UMIs
    """
    with Samfile(args.abam, 'rb') as in_bam, Samfile(args.bbam, 'wb', template=in_bam) as out_bam:

        for chrom in in_bam.references:
            print("processing chromosome", chrom, file=sys.stderr)

            umi_idx = defaultdict(set)
            read_counts = Counter()

            for read in in_bam.fetch(chrom):
                if read.is_unmapped:
                    continue

                # get the iupac umi sequence
                umi = umi_from_name(read.qname)

                # get actual read start
                # read.pos accounts for 5' soft clipping
                if read.is_reverse:
                    # read.alen alignment length accounting for 3' soft clipping
                    # UMIs are then compared to reads with the same start
                    read_start = read.pos + read.alen
                else:
                    read_start = read.pos

                # add count for this start; counts all reads
                read_counts[read_start] += 1

                # check if UMI seen
                if umi in umi_idx[read_start]:
                    continue
                elif args.mismatches > 0 and passing_distances(umi, umi_idx[read_start], args.mismatches):
                    # add UMI as unique hit
                    umi_idx[read_start].add(umi)
                    continue

                # keep track of unique UMIs - set eliminates duplicates
                umi_idx[read_start].add(umi)

                out_bam.write(read)

            # process before and after counts over chrom
            for start, before_count in sorted(read_counts.items()):
                print(chrom, start, start + 1, before_count, len(umi_idx[start]), sep="\t")


def readfq(filehandle):
    """Fastq iterator.

    Args:
        filehandle (file): open file handle

    Yields:
        Fastq
    """
    fqclean = (x.strip("\r\n") for x in filehandle if x.strip())
    while True:
        rd = [x for x in islice(fqclean, 4)]
        if not rd:
            raise StopIteration
        assert all(rd) and len(rd) == 4
        yield Fastq(rd)


def valid_umi(iupac, umi):
    """Parse UMI sequence to validate against IUPAC sequence.

    Args:
        iupac (str): IUPAC sequence
        umi (str): observed sequence

    Returns:
        bool

    >>> valid_umi("NNNV", "ACGT")
    False
    >>> valid_umi("NNNV", "ACGG")
    True
    """
    for code, base in zip(iupac, umi):
        try:
            if base not in IUPAC[code]:
                return False
        except KeyError:
            return False
    return True


def clip_umi(record, iupac_umi, n, end):
    """Removed UMI sequence from read, trims respective length from qual, then appends UMI onto read name.

    Args:
        record (Fastq): `Fastq` record
        iupac_umi (str): IUPAC sequence of the UMI
        n (int): Length of the UMI
        end (int): The end of the read on which the UMI resides

    Returns:
        Fastq else str: The record or the failed UMI sequence

    >>> fq = Fastq(["@cluster_455 2","GGGGGAGCCACGAGGTGTGTTTTATTTTCATTATTC","+","C===>=B=@:<;4A;8=9?6EEC0?DDA72B@3EB4"])
    >>> clip_umi(fq, "NNNNNV", 6, 5)
    Fastq(cluster_455:UMI_GGGGGA 2)
    >>> fq = Fastq(["@cluster_455 2","GGXXGAGCCACGAGGTGTGTTTTATTTTCATTATTC","+","C===>=B=@:<;4A;8=9?6EEC0?DDA72B@3EB4"])
    >>> clip_umi(fq, "NNNNNV", 6, 5)
    'GGXXGA'
    """
    if end == 5:
        umi = record.seq[:n]
        record.seq = record.seq[n:]
        record.qual = record.qual[n:]
    else:
        umi = record.seq[-n:]
        record.seq = record.seq[:-n]
        record.qual = record.qual[:-n]
    if not valid_umi(iupac_umi, umi):
        return umi
    try:
        name, pair = record.name.split(" ", 1)
        record.name = "{name}:UMI_{umi} {pair}".format(name=name, umi=umi, pair=pair)
    except ValueError:
        record.name = "{name}:UMI_{umi}".format(name=record.name, umi=umi)
    return record


def process_fastq(fastq, umi, end=5, verbose=False, top=10):
    """For every valid umi, trim while incorporating UMI into read name.

    Args:
        fastq (str): file path to unprocessed FASTQ file
        umi (str): IUPAC sequence of UMI
        end (Optional[int]): 5 or 3, which ever end you're UMI is located on
        verbose (Optional[bool]): True prints basic stats on observed UMIs
        top (Optional[int]): Number of the the top invalid UMIs to print out
    """
    umi_stats = Counter()
    umi = umi.upper()
    u_leng = len(umi)
    with gzopen(fastq) as fq:
        for read in readfq(fq):
            read = clip_umi(read, umi, u_leng, end)
            if type(read) is Fastq:
                print(read)
            else:
                umi_stats.update([read])
    if verbose:
        print("Invalid UMI Total:", sum(umi_stats.values()), file=sys.stderr)
        print("Unique UMIs Removed:", len(list(umi_stats)), file=sys.stderr)
        print("Top", top, "Invalid UMIs:", file=sys.stderr)
        for umi, val in umi_stats.most_common(top):
            print(umi, val, sep="\t", file=sys.stderr)


def main():

    def _file_exists(parser, arg):
        if not os.path.exists(arg):
            parser.error("The file %s does not exist" % arg)
        if not os.path.isfile(arg):
            parser.error("Expected file, not folder (%s)" % arg)
        return arg

    p = ArgumentParser(description=__doc__)
    p.add_argument('--version', action='version', version='%(prog)s {version}'.format(version=__version__))
    subp = p.add_subparsers(help='commands', dest='command')

    # fastq processing
    fastq = subp.add_parser('trim', description=("Trims the UMI sequence from the read, incorporating the unique "
                                                 "sequence in the read name facilitating filtering of the alignments."),
                            formatter_class=ArgumentDefaultsHelpFormatter,
                            help="trim UMI and incorporate sequence into read name")
    fastq.add_argument('fastq', metavar='FASTQ', type=lambda x: _file_exists(p, x),
                       help='reads with untrimmed UMI')
    fastq.add_argument('umi', metavar='UMI',
                       help='IUPAC UMI sequence, e.g. NNNNNV')
    fastq.add_argument('--end', choices=['5', '3'], default=5, type=int,
                       help="UMI location on the read")
    fastq.add_argument('--verbose', action='store_true',
                       help="print UMI stats to stderr")
    fastq.add_argument('--top', type=int, default=10,
                       help="when verbose, print this many of the top filtered UMI sequences")

    # bam processing
    bam = subp.add_parser('rmdup', description=("Removes duplicate reads, that were previously characterized by "
                                                "their UMI, at any given start location. Coverage differences before "
                                                "and after are written to STDOUT as BED3+."),
                          formatter_class=ArgumentDefaultsHelpFormatter,
                          help="remove duplicate UMI entries from all start positions")
    bam.add_argument('abam', metavar='INPUT_BAM', type=lambda x: _file_exists(p, x),
                     help='bam with UMI in read name')
    bam.add_argument('bbam', metavar='OUTPUT_BAM',
                     help='non-duplicate UMIs at any given start position')
    bam.add_argument('-m', '--mismatches', default=0, type=int,
                     help="allowable mismatches when comparing UMIs at any given start location")

    args = p.parse_args()
    if args.command == 'trim':
        process_fastq(args.fastq, args.umi, args.end, args.verbose, args.top)
    elif args.command == 'rmdup':
        process_bam(args.abam, args.bbam, args.mismatches)


if __name__ == "__main__":
    main()
