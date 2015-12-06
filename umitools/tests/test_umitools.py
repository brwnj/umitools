import os
import sys
from contextlib import contextmanager
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
from nose.tools import raises
from pysam import Samfile

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


def test_process_fastq_without_invalid():
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


def test_process_fastq_with_invalid():
    umi = "NNNNNGGG"
    end = 5
    verbose = True
    top = 6
    invalid_fastq = os.path.join(DATA, "invalid.fastq")
    with captured_output() as (out, err):
        process_fastq(os.path.join(DATA, "t.fastq"), umi, end=end, invalid=invalid_fastq, verbose=verbose, top=top)
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
    with open(invalid_fastq) as fh:
        it = iter(fh)
        assert it.next().strip() == "@HWI-700819F:306:C72JMACXX:2:1101:1902:2060:UMI_NGCAAGGN 1:N:0:ACGAGTCT"
        assert it.next().strip() == "ATGACCCACCAATCGCATGCCTATCATATAGTAAAACCCAGCCCATGACCCCTAACAGGGGCCCTCTCAGCCCTCCTAATGACCTCCGGCCT"
        assert it.next().strip() == "+"
        assert it.next().strip() == "FFFFFIIIIIIIIIIIIIIIIIIIIIIIIIIIIIBFIIIIIIIIIIIFIIIIIFFFFFFFFFFFFFFBBBBFFFFFFFFFFBFFFFB<BFFF"


@raises(SystemExit)
def test_indexed_unsorted_bam():
    bam = os.path.join(DATA, "unsorted.bam")
    is_indexed(bam)


def test_indexed_sorted_bam():
    bam = os.path.join(DATA, "ordered_umi.bam")
    bai = os.path.join(DATA, "ordered_umi.bam.bai")
    if os.path.exists(bai):
        os.remove(bai)
    is_indexed(bam)
    assert os.path.exists(bai)
    os.remove(bai)


def test_process_bam_no_mismatches():
    tbam = os.path.join(DATA, "tmp.bam")
    bam = os.path.join(DATA, "ordered_umi.bam")
    if os.path.exists(tbam):
        os.remove(tbam)
    with captured_output() as (out, err):
        process_bam(bam, tbam, mismatches=0)
    assert os.path.exists(tbam)
    it = iter(out.getvalue().split("\n"))
    assert it.next().strip() == "1\t9\t10\t4\t4"
    assert it.next().strip() == "1\t11\t12\t2\t1"
    assert it.next().strip() == "1\t29\t30\t2\t2"

    bam_reader = Samfile(tbam)
    it = iter(bam_reader)
    r = it.next()
    assert r.pos == 4
    assert r.qname == "read8:UMI_ATTCAGGG"
    r = it.next()
    assert r.pos == 9
    assert r.qname == "read1:UMI_AAAAAGGG"
    r = it.next()
    assert r.pos == 9
    assert r.qname == "read2:UMI_AAAATGGG"
    r = it.next()
    assert r.pos == 9
    assert r.qname == "read3:UMI_AAAACGGG"
    r = it.next()
    assert r.pos == 9
    assert r.qname == "read4:UMI_AAAGGGGG"
    r = it.next()
    assert r.pos == 11
    assert r.qname == "read5:UMI_ATTTAGGG"
    r = it.next()
    assert r.pos == 29
    assert r.qname == "read7:UMI_ATTTAGGG"
    bam_reader.close()
    os.remove(tbam)


def test_process_bam_mismatches():
    tbam = os.path.join(DATA, "tmp.bam")
    bam = os.path.join(DATA, "ordered_umi.bam")
    if os.path.exists(tbam):
        os.remove(tbam)
    with captured_output() as (out, err):
        process_bam(bam, tbam, mismatches=1)
    assert os.path.exists(tbam)
    it = iter(out.getvalue().split("\n"))
    assert it.next().strip() == "1\t9\t10\t4\t2"
    assert it.next().strip() == "1\t11\t12\t2\t1"
    assert it.next().strip() == "1\t29\t30\t2\t1"

    bam_reader = Samfile(tbam)
    it = iter(bam_reader)
    r = it.next()
    assert r.pos == 4
    assert r.qname == "read8:UMI_ATTCAGGG"
    r = it.next()
    assert r.pos == 9
    assert r.qname == "read1:UMI_AAAAAGGG"
    r = it.next()
    assert r.pos == 9
    assert r.qname == "read4:UMI_AAAGGGGG"
    r = it.next()
    assert r.pos == 11
    assert r.qname == "read5:UMI_ATTTAGGG"
    bam_reader.close()
    os.remove(tbam)


@raises(UMINotFound)
def test_process_bam_raises():
    bam = os.path.join(DATA, "no_umi.bam")
    tbam = os.path.join(DATA, "tmp.bam")
    with captured_output() as (o, e):
        process_bam(bam, tbam)
    if os.path.exists(tbam):
        os.remove(tbam)
