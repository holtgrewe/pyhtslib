#!/usr/bin/env python
"""Tests for reading BCF/VCF files sequentially or through indices"""

import pyhtslib.bcf as bcf

from tests.bcf_fixtures import *  # NOQA

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'

QNAMES = ['I', 'II.14978392', 'III', 'IV', 'V', 'VI']
END_POS = [102, 102, 102, 102, 102, 100101]
FLAGS = [16, 16, 16, 16, 16, 2048]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def check_file(f):
    records = []
    for i, record in enumerate(f):
        assert record.r_id == 0
        assert record.ref == 'CHROMOSOME_I'
        assert record.qname == QNAMES[i]
        assert record.mapq == 1
        assert record.begin_pos == 1
        assert record.end_pos == END_POS[i]
        assert record.flag == FLAGS[i]
        records.append(record.detach())

    assert len(records) == 6

    assert records[5].is_supplementary
    assert records[0].is_reversed
    assert all([records[i].flag == 16 for i in range(5)])
    assert all([records[i].is_reversed for i in range(5)])

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_six_records_read_sequential_vcf(six_records_vcf):
    with bcf.BCFFile(str(six_records_vcf)) as f:
        check_file(f)


def test_six_records_read_sequential_vcf_gz(six_records_vcf_gz):
    with bcf.BCFFile(str(six_records_vcf_gz)) as f:
        check_file(f)


def test_six_records_read_sequential_bcf(six_records_bcf):
    with bcf.BCFFile(str(six_records_bcf)) as f:
        check_file(f)


def test_two_hundred_read_sequential_vcf(two_hundred_vcf):
    with bcf.BCFFile(str(two_hundred_vcf)) as f:
        num = len(list(f))
        assert num == 200


def test_two_hundred_read_sequential_vcf_gz(two_hundred_vcf_gz):
    with bcf.BCFFile(str(two_hundred_vcf_gz)) as f:
        num = len(list(f))
        assert num == 200


def test_two_hundred_read_sequential_bcf(two_hundred_bcf):
    with bcf.BCFFile(str(two_hundred_bcf)) as f:
        num = len(list(f))
        assert num == 200


def test_two_hundread_through_index_vcf_gz(
        two_hundred_vcf_gz, two_hundred_tbi):
    with bcf.BCFIndex(str(two_hundred_vcf_gz)) as idx:
        assert len(list(idx.query('chr17:10,000,000-11,000,000'))) == 2


def test_two_hundread_through_index_bcf(two_hundred_bcf, two_hundred_bai):
    with bcf.BCFIndex(str(two_hundred_bcf)) as idx:
        assert len(list(idx.query('chr17:10,000,000-11,000,000'))) == 2
        assert len(list(idx.query('chr17:10,000,000-15,000,000'))) == 12
