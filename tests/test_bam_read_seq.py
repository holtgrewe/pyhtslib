#!/usr/bin/env python
"""Tests for reading BAM files sequentially

"""

# TODO(holtgrewe): tests for CRAM

import pyhtslib.bam as bam

from tests.bam_fixtures import *  # NOQA

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'


def test_read_sequential_sam(six_records_sam):
    with bam.BAMFile(str(six_records_sam)) as f:
        count = 0
        for record in f:
            count += 1
        assert count == 6


def test_read_sequential_sam_gz(six_records_sam_gz):
    with bam.BAMFile(str(six_records_sam_gz)) as f:
        count = 0
        for record in f:
            count += 1
        assert count == 6


def test_read_sequential_bam(six_records_bam):
    with bam.BAMFile(str(six_records_bam)) as f:
        count = 0
        for record in f:
            count += 1
        assert count == 6
