#!/usr/bin/env python
# encoding: utf-8
"""
Tools to handle reads sequenced with unique molecular identifiers (UMIs).
"""
import os
import sys
import doctest
from re import findall
from pysam import Samfile
from toolshed import nopen
from collections import Counter
from itertools import islice, izip
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

os.environ['TERM'] = 'linux'

__version__ = "version 0.1.5"

IUPAC = {"A":"A","T":"T","C":"C","G":"G","R":"GA","Y":"TC",
         "M":"AC","K":"GT","S":"GC","W":"AT","H":"ACT",
         "B":"GTC","V":"GCA","D":"GAT","N":"GATC"}

class Fastq(object):
    def __init__(self, args):
        self.name = args[0][1:]
        self.seq = args[1]
        self.qual = args[3]
        assert len(self.seq) == len(self.qual)

    def __repr__(self):
        return "Fastq({name})".format(name=self.name)

    def __str__(self):
        return "@{name}\n{seq}\n+\n{qual}".format(name=self.name,
                seq=self.seq, qual=self.qual)

def umi_from_name(name):
    """
    extract the UMI sequence from the read name.

    >>> umi_from_name("cluster_1017333:UMI_GCCGCA")
    'GCCGCA'
    """
    return findall(r'UMI_([\w]*)', name)[0].strip()

def add_strand(umi, is_reverse):
    """
    add strand info onto UMI, allowing pos and neg strand to saturate UMI.

    >>> add_strand("GCCGCA", True)
    'GCCGCAneg'
    >>> add_strand("GCCGCA", False)
    'GCCGCApos'
    """
    return "%sneg" % umi if is_reverse else "%spos" % umi

def get_chromosomes(sam_header):
    """
    parse sam header to return SN values within SQ lines.
    """
    chromosomes = []

    for line in sam_header.split("\n"):
        if not line.startswith("@SQ"): continue
        line = line.strip().split("\t")

        for token in line:
            if not token.startswith("SN"): continue
            chromosome = token.split(":", 1)[1]
            chromosomes.append(chromosome)

    return chromosomes

def process_bam(args):
    """
    removes duplicate reads characterized by their UMI at any given start
    location.
    """
    with Samfile(args.abam, 'rb') as in_bam, Samfile(args.bbam, 'wb', template=in_bam) as out_bam:
        chromosomes = get_chromosomes(in_bam.text)

        for chrom in chromosomes:
            print >>sys.stderr, "processing chromosome", chrom
            umi_idx = {}

            for read in in_bam.fetch(chrom):
                if read.is_unmapped: continue
                # get the iupac umi sequence
                umi = umi_from_name(read.qname)
                # add strand onto umi before adding to index
                umi = add_strand(umi, read.is_reverse)
                # get actual read start
                # read.pos accounts for 5' soft clipping
                if read.is_reverse:
                    # read.alen alignment length accounting for 3' soft clipping
                    # UMIs are then compared to reads with the same start
                    read_start = read.pos + read.alen
                else:
                    read_start = read.pos
                # check for duplicate UMI
                try:
                    if umi in umi_idx[read_start]:
                        continue
                    umi_idx[read_start].add(umi)
                except KeyError:
                    umi_idx[read_start] = {umi}

                out_bam.write(read)

def readfq(fq):
    with nopen(fq) as fh:
        fqclean = (x.strip("\r\n") for x in fh if x.strip())
        while True:
            rd = [x for x in islice(fqclean, 4)]
            if not rd: raise StopIteration
            assert all(rd) and len(rd) == 4
            yield Fastq(rd)

def valid_umi(iupac, umi):
    """
    parse UMI sequence to validate against IUPAC sequence.

    >>> valid_umi("NNNV", "ACGT")
    False
    >>> valid_umi("NNNV", "ACGG")
    True
    """
    for code, base in izip(iupac, umi):
        try:
            if not base in IUPAC[code]:
                return False
        except KeyError:
            return False
    return True

def clip_umi(record, iupac_umi, n, end):
    """
    >>> fq = Fastq(["@cluster_455 2",\
                        "GGGGGAGCCACGAGGTGTGTTTTATTTTCATTATTC",\
                        "+",\
                        "C===>=B=@:<;4A;8=9?6EEC0?DDA72B@3EB4"])
    >>> clip_umi(fq, "NNNNNV", 6, "5")
    Fastq(cluster_455:UMI_GGGGGA 2)
    >>> fq = Fastq(["@cluster_455 2",\
                        "GGXXGAGCCACGAGGTGTGTTTTATTTTCATTATTC",\
                        "+",\
                        "C===>=B=@:<;4A;8=9?6EEC0?DDA72B@3EB4"])
    >>> clip_umi(fq, "NNNNNV", 6, "5")
    'GGXXGA'
    """
    if end == "5":
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
        record.name = "{name}:UMI_{umi} {pair}".format(name=name,
                                                        umi=umi,
                                                        pair=pair)
    except ValueError:
        record.name = "{name}:UMI_{umi}".format(name=record.name, umi=umi)
    return record

def process_fastq(args):
    """
    for every valid umi, trim while incorporating into read name.
    """
    umi_stats = Counter()
    iupac = args.umi
    u_leng = len(args.umi)
    end = args.end
    for r in readfq(args.fastq):
        r = clip_umi(r, iupac, u_leng, end)
        if type(r) is Fastq:
            print r
        else:
            umi_stats.update([r])
    if args.verbose:
        print >>sys.stderr, "Invalid UMI Total:   {count}".format(count=sum(umi_stats.values()))
        print >>sys.stderr, "Unique UMIs Removed: {count}".format(count=len(list(umi_stats)))
        print >>sys.stderr, "Top {count} Invalid UMIs:".format(count=args.top)
        for umi, val in umi_stats.most_common(args.top):
            print >>sys.stderr, "\t".join([umi, str(val)])


if __name__ == "__main__":

    p = ArgumentParser(description=__doc__, version=__version__)
    subp = p.add_subparsers(help='commands')

    # fastq processing
    fastq = subp.add_parser('trim', description="Trims the UMI sequence from \
            the read, incorporating the unique sequence in the read name \
            facilitating filtering of the alignments.",
            formatter_class=ArgumentDefaultsHelpFormatter,
            help="trim UMI and incorporate sequence into read name")
    fastq.add_argument('fastq', metavar='FASTQ',
            help='reads with untrimmed UMI')
    fastq.add_argument('umi', metavar='UMI',
            help='IUPAC UMI sequence, e.g. NNNNNV')
    fastq.add_argument('--end', choices=['5', '3'], default="5",
            help="UMI location on the read")
    fastq.add_argument('--verbose', action='store_true',
            help="print UMI stats to stderr")
    fastq.add_argument('--top', type=int, default=10,
            help="when verbose, print this many of the top filtered \
            UMI sequences")
    fastq.set_defaults(func=process_fastq)

    # bam processing
    bam = subp.add_parser('rmdup', description="Removes duplicate reads, that \
            were previously characterized by their UMI, at any given start \
            location.",
            formatter_class=ArgumentDefaultsHelpFormatter,
            help="remove duplicate UMI entries from all start positions")
    bam.add_argument('abam', metavar='INPUT_BAM',
            help='bam with UMI in read name')
    bam.add_argument('bbam', metavar='OUTPUT_BAM',
            help='non-duplicate UMIs at any given start position')
    bam.add_argument('umi', metavar='UMI',
            help='IUPAC sequence of the UMI, e.g. NNNNNV')
    bam.set_defaults(func=process_bam)

    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        args = p.parse_args()
        args.func(args)
