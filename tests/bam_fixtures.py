#!/usr/bin/env python
"""Fixture files for the pyhtslib.bam tests"""

import os
import py
import pytest

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
def header_only_bam(tmpdir):
    """Copy the header_only.bam file to temporary directory."""
    src = py.path.local(os.path.dirname(__file__)).join(
        'files', 'header_only.bam')
    dst = tmpdir.join('header_only.bam')
    src.copy(dst)
    yield dst
    dst.remove()


@pytest.yield_fixture
def header_only_bai(tmpdir):
    """Copy the header_only.bam.bai file to temporary directory."""
    src = py.path.local(os.path.dirname(__file__)).join(
        'files', 'header_only.bam.bai')
    dst = tmpdir.join('header_only.bam.bai')
    src.copy(dst)
    yield dst
    dst.remove()
