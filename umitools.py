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
from itertools import islice, izip, groupby
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

os.environ['TERM'] = 'linux'

__version__ = "version 0.1.4"

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

def decode(x):
    """not accounting for other scoring
    
    >>> decode("4")
    19
    >>> decode(";")
    26
    """
    return ord(x) - 33

def average(quals):
    """
    >>> average("4578;;77;;H@GFF>DFFFFFEEE")
    30.84...
    >>> average("/==96996<FGCHHHGGGFFE=EDFFFEEB")
    32.53...
    """
    vals = map(decode, quals)
    return sum(vals)/float(len(quals))

def write_reads(bam, reads, cap):
    """writes the read into bam file and verifies UMI saturation cap."""
    pos = 0
    neg = 0
    for umi, read in reads.iteritems():
        if umi.endswith("neg"): neg += 1
        if umi.endswith("pos"): pos += 1
        bam.write(read)

    assert pos <= cap
    assert neg <= cap

def total_possible(umi):
    """total possible UMI combinations.
    >>> total_possible("NNNNNV")
    3072
    >>> total_possible("NNNNNR")
    2048
    """
    return reduce(lambda x, y: x * y, [len(IUPAC[n]) for n in umi])

def umi_from_name(name):
    """extract the UMI sequence from the read name.
    >>> umi_from_name("cluster_1017333:UMI_GCCGCA")
    'GCCGCA'
    """
    return findall(r'UMI_([\w]*)', name)[0].strip()

def add_strand(umi, reverse):
    """add strand info onto UMI, allowing pos and neg strand to saturate UMI.
    >>> add_strand("GCCGCA", True)
    'GCCGCAneg'
    >>> add_strand("GCCGCA", False)
    'GCCGCApos'
    """
    return "%sneg" % umi if reverse else "%spos" % umi

def group_starts(fh):
    for key, grp in groupby(fh, key=lambda read: read.pos):
        yield grp

def process_bam(args):
    """removes duplicate reads characterized by their UMI at any given start
    location.
    """
    possible = total_possible(args.umi)
    filtered, seen = 0, 0
    with Samfile(args.abam, 'rb') as sam, Samfile(args.bbam, 'wb', template=sam) as bam:
        for reads in group_starts(sam.fetch()):
            unique_reads = {}
            for read in reads:
                seen += 1
                umi = add_strand(umi_from_name(read.qname), read.is_reverse)
                try:
                    if average(unique_reads[umi].qqual) < average(read.qqual):
                        unique_reads[umi] = read
                    filtered += 1
                except KeyError:
                    unique_reads[umi] = read
            write_reads(bam, unique_reads, possible)
    if args.verbose:
        print >>sys.stderr, \
            "Input Sequences:     {seen}".format(seen=seen)
        print >>sys.stderr, \
            "Sequences Removed:   {filtered}".format(filtered=filtered)
        print >>sys.stderr, \
            "Remaining Sequences: {difference}".format(difference=seen - filtered)

def readfq(fq):    
    with nopen(fq) as fh:
        fqclean = (x.strip("\r\n") for x in fh if x.strip())
        while True:
            rd = [x for x in islice(fqclean, 4)]
            if not rd: raise StopIteration
            assert all(rd) and len(rd) == 4
            yield Fastq(rd)

def valid_umi(iupac, umi):
    """parse UMI sequence to validate against IUPAC sequence.
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
    """for every valid umi, trim while incorporating into read name."""
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

def main(args):
    args.func(args)

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
    bam.add_argument('--verbose', action='store_true',
            help="print rmdup stats")
    bam.set_defaults(func=process_bam)

    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        args = p.parse_args()
        main(args)
