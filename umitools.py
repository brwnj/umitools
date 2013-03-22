#!/usr/bin/env python
# encoding: utf-8
"""
Tools to handle reads sequenced with unique molecular identifiers (UMIs).
"""
import sys
from pysam import Samfile
from toolshed import nopen
from itertools import islice, izip

__version__ = "v0.1.1"

def get_umi(name):
    """parses read name and returns umi string"""
    toks = name.split(":")
    for tok in toks:
        if tok.startswith("UMI"):
            return tok.lstrip("UMI_")

def process_bam(args):
    """removes duplicate reads characterized by their UMI at any given start 
    location.
    """
    samfile = Samfile(args.abam, 'rb')
    with Samfile(args.bbam, 'wb', template=samfile) as bamfile:
        reads = []
        umis = []
        pos = 0
        for read in samfile.fetch():
            umi = get_umi(read.qname)
            # same start coord
            if read.pos == pos:
                # skip non-unique UMIs
                if not umi in umis:
                    umis.append(umi)
                    reads.append(read)
            else:
                if len(reads) > 0:
                    for alignedread in reads:
                        bamfile.write(alignedread)
                reads = [read]
                umis = [umi]
                pos = read.pos
        if len(reads) > 0:
            for alignedread in reads:
                bamfile.write(alignedread)

def read_fastq(fh):
    """fastq parser that returns name, seq, and qual."""
    values = list(islice(fh, 4))
    if len(values) == 4:
        id1, seq, id2, qual = values
    elif len(values) == 0:
        raise StopIteration
    else:
        raise EOFError("unexpected end of file")
    assert id1.startswith('@'),\
            ">> Fastq out of sync at read:\n%s\n" % id1
    assert id2.startswith('+'),\
            ">> Fastq out of sync at read:\n%s\n" % id1
    assert len(seq) == len(qual),\
            ">> Sequence and Quality are not the same length \
            for read:\n%s\n" % id1
    yield id1[1:-1], seq[:-1], qual[:-1]

def valid_umi(iupac, umi):
    """parse UMI sequence to validate against IUPAC sequence."""
    IUPAC_definitions = {"A":"A","T":"T","C":"C","G":"G","R":"GA","Y":"TC",
                            "M":"AC","K":"GT","S":"GC","W":"AT","H":"ACT",
                            "B":"GTC","V":"GCA","D":"GAT","N":"GATC"}
    for code, base in izip(iupac, umi):
        try:
            if not base in IUPAC_definitions[code]:
                return False
        except KeyError:
            return False
    return True

def process_fastq(args):
    """reads in fastq with UMI still in the sequence:
        + verifies the UMI as valid
        + trims the UMI
        + incorporate UMI sequence into read name
    """
    umi_length = len(args.umi)
    with nopen(fastq) as fq:
        for label, seq, qual in read_fastq(fq):
            # TODO: add 3' support
            umi = seq[:umi_length]
            if valid_umi(args.umi, umi):
                print "@%s:UMI_%s\n%s\n+\n%s" % \
                    (label, umi, seq[umi_length:], qual[umi_length:])
            else:
                print >> sys.stderr, ">> Read %s ignored. Invalid UMI (%s)\n" % \
                    (label, umi)

def main(args):
    args.func(args)

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--version', action="version", version="%(prog)s " + __version__)
    subp = p.add_subparsers(help='commands')
    # bam processing
    bam = subp.add_parser('rmdup', help="remove duplicate UMI \
            entries at originating from all start positions")
    bam.add_argument('abam', metavar='INPUT_BAM', help='bam with UMI in \
            read name')
    bam.add_argument('bbam', metavar='OUTPUT_BAM', help="bam with no \
            duplicate UMIs at any given 5' location")
    bam.set_defaults(func=process_bam)
    # fastq processing
    fastq = subp.add_parser('trim', help="trim 5' UMI and \
            incorporate sequence into read name")
    fastq.add_argument('fastq', metavar='FASTQ', help='reads with \
            untrimmed UMI')
    fastq.add_argument('umi', metavar='UMI', help='IUPAC UMI sequence, \
            e.g. NNNNNV')
    fastq.set_defaults(func=process_fastq)
    main(p.parse_args())