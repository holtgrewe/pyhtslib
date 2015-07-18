#!/usr/bin/env python3
"""Tests for the module pyhtslib"""

import pyhtslib

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_genome_interval():
    gitv = pyhtslib.GenomeInterval('chr1', 1000, 2000)
    assert str(gitv) == 'chr1:1,001-2,000'
