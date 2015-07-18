#!/usr/bin/env python3
"""Tests for the module pyhtslib.tabix"""

import os
import py
import pytest

import pyhtslib.tabix as tabix

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.yield_fixture
def reduced_pg_vcf(tmpdir):
    """Copy the reduced_pg.vcf.gz file to temporary directory

    This file contains 112 variants from the
    Illumina platinum genome for NA12877.
    """
    src = py.path.local(os.path.dirname(__file__)).join(
        'files', 'reduced_pg.vcf.gz')
    dst = tmpdir.join('reduced_pg.vcf.gz')
    src.copy(dst)
    yield dst
    dst.remove()


@pytest.yield_fixture
def reduced_pg_tbi(tmpdir):
    """Copy the reduced_pg.vcf.gz.tbi file to temporary directory

    This is the index for reduced_pg.vcf.gz.
    """
    src = py.path.local(os.path.dirname(__file__)).join(
        'files', 'reduced_pg.vcf.gz.tbi')
    dst = tmpdir.join('reduced_pg.vf.gz.tbi')
    src.copy(dst)
    yield dst
    dst.remove()

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_vcf_tabix_index_auto_build(reduced_pg_vcf):
    pass


def test_vcf_tabix_index_auto_load(reduced_pg_vcf):
    pass


def test_vcf_tabix_index_require_index_fails(reduced_pg_vcf):
    pass


def test_vcf_tabix_index_build_index(reduced_pg_vcf):
    pass


def test_vcf_tabix_load_whole_file(reduced_pg_vcf):
    pass


def test_vcf_tabix_load_chr3(reduced_pg_vcf):
    pass


def test_vcf_tabix_load_chr3_region(reduced_pg_vcf):
    pass
