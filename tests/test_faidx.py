#!/usr/bin/env python
"""Tests for module pyhtslib.faidx."""

import os
import py
import pytest

import pyhtslib.faidx as faidx

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.yield_fixture
def two_genes_fasta(tmpdir):
    """Copy the two_genes.fa file to temporary directory."""
    src = py.path.local(os.path.dirname(__file__)).join(
        'files', 'two_genes.fa')
    dst = tmpdir.join('two_genes.fa')
    src.copy(dst)
    yield dst
    dst.remove()


@pytest.yield_fixture
def two_genes_fai(tmpdir):
    """Copy the two_genes.fa.fai file to temporary directory."""
    src = py.path.local(os.path.dirname(__file__)).join(
        'files', 'two_genes.fa.fai')
    dst = tmpdir.join('two_genes.fa.fai')
    src.copy(dst)
    yield dst
    dst.remove()

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_fasta_index_auto_build(two_genes_fasta):
    with faidx.FASTAIndex(str(two_genes_fasta)) as idx:
        assert idx.fasta_path == two_genes_fasta
        assert idx.struct_ptr
        assert len(idx.seq_dict) == 2
        assert idx.seq_dict['HSBGPG'].length == 1231
        assert idx.seq_dict['HSGLTH1'].length == 1020
        idx.close()  # otherwise there are leaks


def test_fasta_index_auto_load(two_genes_fasta, two_genes_fai):
    with faidx.FASTAIndex(str(two_genes_fasta)) as idx:
        assert idx.fasta_path == two_genes_fasta
        assert idx.struct_ptr
        assert len(idx.seq_dict) == 2
        assert idx.seq_dict['HSBGPG'].length == 1231
        assert idx.seq_dict['HSGLTH1'].length == 1020


def test_fasta_index_require_index_fails(two_genes_fasta):
    with pytest.raises(faidx.FASTAIndexException):
        faidx.FASTAIndex(str(two_genes_fasta), require_index=True,
                         auto_build=False, auto_load=False)

def test_fasta_index_load(two_genes_fasta, two_genes_fai):
    with faidx.FASTAIndex(str(two_genes_fasta)) as idx:
        s = idx.fetch('HSBGPG:5-100')
