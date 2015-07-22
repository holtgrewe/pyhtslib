#!/usr/bin/env python
"""Fixture files for the pyhtslib.bcf tests"""

import os
import py
import pytest

import pyhtslib.bcf_internal as bcf_internal
import pyhtslib.hts_internal as hts_internal

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.yield_fixture
def header_only_vcf(tmpdir):
    """Copy the header_only.vcf file to temporary directory."""
    src = py.path.local(os.path.dirname(__file__)).join(
        'files', 'header_only.vcf')
    dst = tmpdir.join('header_only.vcf')
    src.copy(dst)
    yield dst
    dst.remove()


@pytest.yield_fixture
def header_only_vcf_header(header_only_vcf):
    hts_file = hts_internal._hts_open(
        str(header_only_vcf).encode('utf-8'), 'r')
    hdr = bcf_internal._vcf_hdr_read(hts_file)
    yield hdr
    bcf_internal._bcf_hdr_destroy(hdr)
    hts_internal._hts_close(hts_file)


@pytest.yield_fixture
def header_only_vcf_gz(tmpdir):
    """Copy the header_only.vcf.gz file to temporary directory."""
    src = py.path.local(os.path.dirname(__file__)).join(
        'files', 'header_only.vcf.gz')
    dst = tmpdir.join('header_only.vcf.gz')
    src.copy(dst)
    yield dst
    dst.remove()


@pytest.yield_fixture
def header_only_vcf_gz_header(header_only_vcf_gz):
    hts_file = hts_internal._hts_open(
        str(header_only_vcf_gz).encode('utf-8'), 'r')
    hdr = bcf_internal._vcf_hdr_read(hts_file)
    yield hdr
    bcf_internal._bcf_hdr_destroy(hdr)
    hts_internal._hts_close(hts_file)


@pytest.yield_fixture
def header_only_bcf(tmpdir):
    """Copy the header_only.bcf file to temporary directory."""
    src = py.path.local(os.path.dirname(__file__)).join(
        'files', 'header_only.bcf')
    dst = tmpdir.join('header_only.bcf')
    src.copy(dst)
    yield dst
    dst.remove()


@pytest.yield_fixture
def header_only_csi(tmpdir):
    """Copy the header_only.bcf.csi file to temporary directory."""
    src = py.path.local(os.path.dirname(__file__)).join(
        'files', 'header_only.bcf.csi')
    dst = tmpdir.join('header_only.bcf.csi')
    src.copy(dst)
    yield dst
    dst.remove()


