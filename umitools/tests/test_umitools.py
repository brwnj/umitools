import os
import sys
from contextlib import contextmanager
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
from nose.tools import raises

from umitools.umitools import is_indexed, process_bam, process_fastq, UMINotFound


DATA = os.path.join(os.path.dirname(__file__), "data")


@contextmanager
def captured_output():
    new_out, new_err = StringIO(), StringIO()
    old_out, old_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = new_out, new_err
        yield sys.stdout, sys.stderr
    finally:
        sys.stdout, sys.stderr = old_out, old_err


def test_process_fastq():
    umi = "NNNNNGGG"
    end = 5
    verbose = True
    top = 6
    with captured_output() as (out, err):
        process_fastq(os.path.join(DATA, "t.fastq"), umi, end=end, verbose=verbose, top=top)
    it = iter(out.getvalue().split("\n"))
    assert it.next().strip() == "@HWI-700819F:306:C72JMACXX:2:1101:1095:2117:UMI_ATTTAGGG 1:N:0:ACGAGTCT"
    assert it.next().strip() == "TCTTCCAAGGTGACAAATTATATAATGAAAAAGCTGTTACCAGAAACTTTCAGCAGACATCTTATTGATAATATTTAATCAGCATTCTCATT"
    assert it.next().strip() == "+"
    assert it.next().strip() == "FFFFFIFFFBBBFFFIIIFFIIIIIIFIIFIIIFFFFFIIIFFFFFFFIIIIIFFIIIIIIIIIFIIBFFFFFBFBFFFFFFFFBFFFFFFF"
    assert it.next().strip() == "@HWI-700819F:306:C72JMACXX:2:1101:1059:2161:UMI_AGATAGGG 1:N:0:ACGAGTCT"
    assert it.next().strip() == "GGTAGCTCACTTCCACTATGTCCTATCAATAGGAGCTGTATTTGCCATCATAGGAGGCTTCATTCACTGATTTCCCCTATTCTCAGGCTACA"
    assert it.next().strip() == "+"
    assert it.next().strip() == "FFBFFIFFIIBFFFBFF<FFBFFFFIFFBFFFIIIBFF<BFFFFI7BFB<77B<FFFFFFFF<07BBFBBBBFFBBBBBBBBFFFB<B7<B<"
    assert it.next().strip() == "@HWI-700819F:306:C72JMACXX:2:1101:1153:2164:UMI_TGATAGGG 1:N:0:ACGAGTCT"
    assert it.next().strip() == "ATCTAATGGTAAATTGATTACCTAATTAGCTGTCTCTTATACACATCTGACGCACGAGTCTTCGTATGCCGTCTTCTGCTTGAAAAAAAAAA"
    assert it.next().strip() == "+"
    assert it.next().strip() == "FFFFFIIIIFFFIIIIFIIFFIIIIIIBFIIIFFIIIIFFBBFFFFBFFIFIIIIIIIFBBBFFBFFFFFFBFFBBFB<<B<7<BBB#####"
    it = iter(err.getvalue().split("\n"))
    assert it.next().strip() == "Invalid UMI Total: 1"
    assert it.next().strip() == "Unique UMIs Removed: 1"
    assert it.next().strip() == "Top 6 Invalid UMIs:"
    assert it.next().strip() == "NGCAAGGN	1"


def test_indexed_bam():
$ umitools rmdup unsorted.bam tmp.bam
[E::hts_idx_push] unsorted positions
index not found for unsorted.bam and indexing failed

def test_process_bam():
$ umitools rmdup --mismatches 1 ordered_umi.bam tmp.bam
processing chromosome 1
1	9	10	4	2
1	11	12	2	1
1	29	30	2	1

$ samtools view tmp.bam
read8:UMI_ATTCAGGG	16	1	5	255	25M	=	142618765	25	CGACCCACTCCGCCATTTTCATCCG	IIGIIIHIGIIFIIIIIIIGIGIII	NM:i:0
read1:UMI_AAAAAGGG	1	1	10	255	25M	=	142618765	25	CGACCCACTCCGCCATTTTCATCCG	IIGIIIHIGIIFIIIIIIIGIGIII	NM:i:0
read4:UMI_AAAGGGGG	1	1	10	255	25M	=	142618765	25	CGACCCACTCCGCCATTTTCATCCG	IIGIIIHIGIIFIIIIIIIGIGIII	NM:i:0
read5:UMI_ATTTAGGG	1	1	12	255	25M	=	142618765	25	CGACCCACTCCGCCATTTTCATCCG	IIGIIIHIGIIFIIIIIIIGIGIII	NM:i:0

for line in Samfile("tmp.bam"):
    print line.pos, line.qname
   ....:
4 read8:UMI_ATTCAGGG
9 read1:UMI_AAAAAGGG
9 read4:UMI_AAAGGGGG
11 read5:UMI_ATTTAGGG


@raises(UMINotFound)
def test_process_bam_raises():
$ umitools rmdup --mismatches 1 no_umi.bam tmp.bam
processing chromosome 1
You may be processing alignments that haven't been annotated with UMIs!
Traceback (most recent call last):
  File "/usr/local/bin/umitools", line 9, in <module>
    load_entry_point('umitools==1.0.3', 'console_scripts', 'umitools')()
  File "/Users/brownj/devel/umitools/umitools/umitools.py", line 368, in main
    process_bam(args.abam, args.bbam, mismatches=args.mismatches)
  File "/Users/brownj/devel/umitools/umitools/umitools.py", line 169, in process_bam
    umi = umi_from_name(read.qname)
  File "/Users/brownj/devel/umitools/umitools/umitools.py", line 111, in umi_from_name
    raise UMINotFound(name)
umitools.umitools.UMINotFound: read8
