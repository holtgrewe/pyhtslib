#!/usr/bin/env python
"""Tests for module pyhtslib.bcf

This file contains tests for the header-related code.
"""

from collections import OrderedDict  # for brevity

import pyhtslib.bcf as bcf

from tests.bcf_fixtures import *  # NOQA

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def check_header_of_header_only_vcf(header):
    assert header.struct_ptr
    assert header.struct
    assert header.target_infos == [
        bcf.BCFHeaderTargetInfo('CHROMOSOME_I', 1009800),
        bcf.BCFHeaderTargetInfo('CHROMOSOME_II', 5000),
        bcf.BCFHeaderTargetInfo('CHROMOSOME_III', 5000),
        bcf.BCFHeaderTargetInfo('CHROMOSOME_IV', 5000),
        bcf.BCFHeaderTargetInfo('CHROMOSOME_V', 5000),
        bcf.BCFHeaderTargetInfo('CHROMOSOME_X', 5000),
        bcf.BCFHeaderTargetInfo('CHROMOSOME_MtDNA', 5000)]
    assert header.header_records == [
        bcf.BCFHeaderRecord('HD', OrderedDict(
            [('VN', '1.0'), ('SO', 'unsorted')])),
        bcf.BCFHeaderRecord('SQ', OrderedDict(
            [('SN', 'CHROMOSOME_I'), ('LN', '1009800')])),
        bcf.BCFHeaderRecord('SQ', OrderedDict(
            [('SN', 'CHROMOSOME_II'), ('LN', '5000')])),
        bcf.BCFHeaderRecord('SQ', OrderedDict(
            [('SN', 'CHROMOSOME_III'), ('LN', '5000')])),
        bcf.BCFHeaderRecord('SQ', OrderedDict(
            [('SN', 'CHROMOSOME_IV'), ('LN', '5000')])),
        bcf.BCFHeaderRecord('SQ', OrderedDict(
            [('SN', 'CHROMOSOME_V'), ('LN', '5000')])),
        bcf.BCFHeaderRecord('SQ', OrderedDict(
            [('SN', 'CHROMOSOME_X'), ('LN', '5000')])),
        bcf.BCFHeaderRecord('SQ', OrderedDict(
            [('SN', 'CHROMOSOME_MtDNA'), ('LN', '5000')])),
        bcf.BCFHeaderRecord('RG', OrderedDict(
            [('ID', 'UNKNOWN'), ('SM', 'UNKNOWN')])),
        bcf.BCFHeaderRecord('PG', OrderedDict(
            [('ID', 'bowtie2'), ('PN', 'bowtie2'), ('VN', '2.0.0-beta5')]))]

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_bcf_header_record():
    rec = bcf.BCFHeaderRecord('HD', [('VN', '1.0'), ('SO', 'unsorted')])
    assert rec.for_vcf_header() == '@HD\tVN:1.0\tSO:unsorted'


def test_bcf_header_comment():
    co = bcf.BCFHeaderComment('foo bar baz!')
    assert co.for_vcf_header() == '@CO\tfoo bar baz!'


def test_read_header_from_vcf(header_only_vcf):
    with bcf.BCFFile(str(header_only_vcf)) as f:
        check_header_of_header_only_vcf(f.header)


def test_read_header_from_vcf_gz(header_only_vcf_gz):
    with bcf.BCFFile(str(header_only_vcf_gz)) as f:
        check_header_of_header_only_vcf(f.header)


def test_read_header_from_bcf(header_only_bcf):
    with bcf.BCFFile(str(header_only_bcf)) as f:
        check_header_of_header_only_vcf(f.header)
