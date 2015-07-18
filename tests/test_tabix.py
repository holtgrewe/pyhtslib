#!/usr/bin/env python3
"""Tests for the module pyhtslib.tabix"""

import os
import py
import pytest
import sys

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
    dst = tmpdir.join('reduced_pg.vcf.gz.tbi')
    src.copy(dst)
    yield dst
    dst.remove()

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def XXtest_vcf_tabix_index_auto_build(reduced_pg_vcf):
    idx = tabix.TabixIndex(str(reduced_pg_vcf))
    assert idx.path == reduced_pg_vcf
    assert idx.tbi_path == str(reduced_pg_vcf) + '.tbi'
    assert idx.file
    assert not idx.require_index
    assert idx.auto_load
    assert idx.auto_build
    assert idx.struct_ptr
    assert idx.struct


def XXtest_vcf_tabix_index_auto_load(reduced_pg_vcf, reduced_pg_tbi):
    idx = tabix.TabixIndex(str(reduced_pg_vcf), auto_build=False)
    assert idx.path == reduced_pg_vcf
    assert idx.tbi_path == str(reduced_pg_vcf) + '.tbi'
    assert idx.file
    assert not idx.require_index
    assert idx.auto_load
    assert not idx.auto_build
    assert idx.struct_ptr
    assert idx.struct


def XXtest_vcf_tabix_index_require_index_fails(reduced_pg_vcf):
    with pytest.raises(tabix.TabixIndexException):
        tabix.TabixIndex(str(reduced_pg_vcf), auto_build=False,
                         auto_load=False, require_index=True)


def XXtest_vcf_tabix_index_build_index(reduced_pg_vcf):
    tabix.TabixIndex.build(str(reduced_pg_vcf))
    assert os.path.exists(str(reduced_pg_vcf) + '.tbi')
    idx = tabix.TabixIndex(str(reduced_pg_vcf), auto_build=False,
                           auto_load=True, require_index=True)
    assert idx.struct_ptr
    assert idx.struct


def test_vcf_tabix_load_whole_file(reduced_pg_vcf, reduced_pg_tbi):
    with tabix.TabixIndex(str(reduced_pg_vcf), require_index=True) as t:
        count = 0
        for l in t.from_start():
            count += 1
    assert count == 112


def test_vcf_tabix_load_chr3(reduced_pg_vcf, reduced_pg_tbi):
    with tabix.TabixIndex(str(reduced_pg_vcf), require_index=True) as t:
        count = 0
        for l in t.query('chr3'):
            count += 1
    assert count == 7


def test_vcf_tabix_load_chr3_region(reduced_pg_vcf, reduced_pg_tbi):
    with tabix.TabixIndex(str(reduced_pg_vcf), require_index=True) as t:
        count = 0
        for l in t.query('chr3:45,000,000-150,000,000'):
            count += 1
    assert count == 3


def test_vcf_tabix_get_header(reduced_pg_vcf, reduced_pg_tbi):
    with tabix.TabixIndex(str(reduced_pg_vcf), require_index=True) as t:
        header = t.get_header()
        assert header.startswith('##fileformat=VCFv4.1')
        assert len(header) == 3598
