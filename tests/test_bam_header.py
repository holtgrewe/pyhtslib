#!/usr/bin/env python
"""Tests for module pyhtslib.bam

This file contains tests for the header-related code.
"""

# TODO(holtgrewe): tests for CRAM

import os
import py
import pytest
import sys  # NOQA  TODO(holtgrew): remove again

from collections import OrderedDict  # for brevity

import pyhtslib.bam as bam
import pyhtslib.bam_internal as bam_internal
import pyhtslib.hts_internal as hts_internal

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.yield_fixture
def header_only_sam(tmpdir):
    """Copy the header_only.sam file to temporary directory."""
    src = py.path.local(os.path.dirname(__file__)).join(
        'files', 'header_only.sam')
    dst = tmpdir.join('header_only.sam')
    src.copy(dst)
    yield dst
    dst.remove()


@pytest.yield_fixture
def header_only_sam_header(header_only_sam):
    hts_file = hts_internal._hts_open(
        str(header_only_sam).encode('utf-8'), 'r')
    hdr = bam_internal._sam_hdr_read(hts_file)
    yield hdr
    bam_internal._bam_hdr_destroy(hdr)
    hts_internal._hts_close(hts_file)


@pytest.yield_fixture
def header_only_sam_gz(tmpdir):
    """Copy the header_only.sam.gz file to temporary directory."""
    src = py.path.local(os.path.dirname(__file__)).join(
        'files', 'header_only.sam.gz')
    dst = tmpdir.join('header_only.sam.gz')
    src.copy(dst)
    yield dst
    dst.remove()


@pytest.yield_fixture
def header_only_sam_gz_header(header_only_sam_gz):
    hts_file = hts_internal._hts_open(
        str(header_only_sam_gz).encode('utf-8'), 'r')
    hdr = bam_internal._sam_hdr_read(hts_file)
    yield hdr
    bam_internal._bam_hdr_destroy(hdr)
    hts_internal._hts_close(hts_file)


@pytest.yield_fixture
def header_only_tbi(tmpdir):
    """Copy the header_only.sam.gz.tbi file to temporary directory."""
    src = py.path.local(os.path.dirname(__file__)).join(
        'files', 'header_only.sam.gz.tbi')
    dst = tmpdir.join('header_only.sam.gz.tbi')
    src.copy(dst)
    yield dst
    dst.remove()


@pytest.yield_fixture
def header_only_bam(tmpdir):
    """Copy the header_only.bam file to temporary directory."""
    src = py.path.local(os.path.dirname(__file__)).join(
        'files', 'header_only.bam')
    dst = tmpdir.join('header_only.bam')
    src.copy(dst)
    yield dst
    dst.remove()


@pytest.yield_fixture
def header_only_bam_header(header_only_bam):
    hts_file = hts_internal._hts_open(
        str(header_only_bam).encode('utf-8'), 'r')
    hdr = bam_internal._sam_hdr_read(hts_file)
    yield hdr
    bam_internal._bam_hdr_destroy(hdr)
    hts_internal._hts_close(hts_file)


@pytest.yield_fixture
def header_only_bai(tmpdir):
    """Copy the header_only.bam.bai file to temporary directory."""
    src = py.path.local(os.path.dirname(__file__)).join(
        'files', 'header_only.bam.bai')
    dst = tmpdir.join('header_only.bam.bai')
    src.copy(dst)
    yield dst
    dst.remove()

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_bam_header_record():
    rec = bam.BAMHeaderRecord('HD', [('VN', '1.0'), ('SO', 'unsorted')])
    assert rec.for_sam_header() == '@HD\tVN:1.0\tSO:unsorted'


def test_bam_header_comment():
    co = bam.BAMHeaderComment('foo bar baz!')
    assert co.for_sam_header() == '@CO\tfoo bar baz!'


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


def test_read_header_from_sam(header_only_sam_header):
    header = bam.BAMHeader(header_only_sam_header)
    check_header_of_header_only_sam(header)


def test_read_header_from_sam_gz(header_only_sam_gz_header):
    header = bam.BAMHeader(header_only_sam_gz_header)
    check_header_of_header_only_sam(header)


def test_read_header_from_bam(header_only_bam_header):
    header = bam.BAMHeader(header_only_bam_header)
    check_header_of_header_only_sam(header)
