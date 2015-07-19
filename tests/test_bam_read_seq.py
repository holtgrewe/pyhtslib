#!/usr/bin/env python
"""Tests for reading BAM files sequentially"""

# TODO(holtgrewe): tests for CRAM

import pyhtslib.bam as bam

from tests.bam_fixtures import *  # NOQA

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


def test_read_sequential_sam(six_records_sam):
    with bam.BAMFile(str(six_records_sam)) as f:
        check_file(f)


def test_read_sequential_sam_gz(six_records_sam_gz):
    with bam.BAMFile(str(six_records_sam_gz)) as f:
        check_file(f)

def test_read_sequential_bam(six_records_bam):
    with bam.BAMFile(str(six_records_bam)) as f:
        check_file(f)
