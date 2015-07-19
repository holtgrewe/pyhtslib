#!/usr/bin/env python
"""Tests for module pyhtslib.bam

This file contains tests for the header-related code.
"""

# TODO(holtgrewe): tests for CRAM

from collections import OrderedDict  # for brevity

import pyhtslib.bam as bam

from tests.bam_fixtures import *  # NOQA

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def check_header_of_header_only_sam(header):
    assert header.struct_ptr
    assert header.struct
    assert header.target_infos == [
        bam.BAMHeaderTargetInfo('CHROMOSOME_I', 1009800),
        bam.BAMHeaderTargetInfo('CHROMOSOME_II', 5000),
        bam.BAMHeaderTargetInfo('CHROMOSOME_III', 5000),
        bam.BAMHeaderTargetInfo('CHROMOSOME_IV', 5000),
        bam.BAMHeaderTargetInfo('CHROMOSOME_V', 5000),
        bam.BAMHeaderTargetInfo('CHROMOSOME_X', 5000),
        bam.BAMHeaderTargetInfo('CHROMOSOME_MtDNA', 5000)]
    assert header.header_records == [
        bam.BAMHeaderRecord('HD', OrderedDict(
            [('VN', '1.0'), ('SO', 'unsorted')])),
        bam.BAMHeaderRecord('SQ', OrderedDict(
            [('SN', 'CHROMOSOME_I'), ('LN', '1009800')])),
        bam.BAMHeaderRecord('SQ', OrderedDict(
            [('SN', 'CHROMOSOME_II'), ('LN', '5000')])),
        bam.BAMHeaderRecord('SQ', OrderedDict(
            [('SN', 'CHROMOSOME_III'), ('LN', '5000')])),
        bam.BAMHeaderRecord('SQ', OrderedDict(
            [('SN', 'CHROMOSOME_IV'), ('LN', '5000')])),
        bam.BAMHeaderRecord('SQ', OrderedDict(
            [('SN', 'CHROMOSOME_V'), ('LN', '5000')])),
        bam.BAMHeaderRecord('SQ', OrderedDict(
            [('SN', 'CHROMOSOME_X'), ('LN', '5000')])),
        bam.BAMHeaderRecord('SQ', OrderedDict(
            [('SN', 'CHROMOSOME_MtDNA'), ('LN', '5000')])),
        bam.BAMHeaderRecord('RG', OrderedDict(
            [('ID', 'UNKNOWN'), ('SM', 'UNKNOWN')])),
        bam.BAMHeaderRecord('PG', OrderedDict(
            [('ID', 'bowtie2'), ('PN', 'bowtie2'), ('VN', '2.0.0-beta5')]))]

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_bam_header_record():
    rec = bam.BAMHeaderRecord('HD', [('VN', '1.0'), ('SO', 'unsorted')])
    assert rec.for_sam_header() == '@HD\tVN:1.0\tSO:unsorted'


def test_bam_header_comment():
    co = bam.BAMHeaderComment('foo bar baz!')
    assert co.for_sam_header() == '@CO\tfoo bar baz!'


def test_read_header_from_sam(header_only_sam):
    with bam.BAMFile(str(header_only_sam)) as f:
        check_header_of_header_only_sam(f.header)


def test_read_header_from_sam_gz(header_only_sam_gz):
    with bam.BAMFile(str(header_only_sam_gz)) as f:
        check_header_of_header_only_sam(f.header)


def test_read_header_from_bam(header_only_bam):
    with bam.BAMFile(str(header_only_bam)) as f:
        check_header_of_header_only_sam(f.header)
