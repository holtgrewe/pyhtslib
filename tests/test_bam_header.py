#!/usr/bin/env python
"""Tests for module pyhtslib.bam

This file contains tests for the header-related code.
"""

import os
import py
import pytest

import pyhtslib.bam as bam

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
def header_only_sam_gz(tmpdir):
    """Copy the header_only.sam.gz file to temporary directory."""
    src = py.path.local(os.path.dirname(__file__)).join(
        'files', 'header_only.sam.gz')
    dst = tmpdir.join('header_only.sam.gz')
    src.copy(dst)
    yield dst
    dst.remove()


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


def test_read_header_from_sam(header_only_sam):
    pass
