#!/usr/bin/env python
"""Tests for reading BCF/VCF files sequentially or through indices"""

import pyhtslib.bcf as bcf

from tests.bcf_fixtures import *  # NOQA

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'

IDS = [['rs17497801'], [], ['rs4149436'], [], ['rs4815270'], ['rs41279999']]
CHROMS = ['2', '2', '2', '19', '20', '22']
RIDS = [11, 11, 11, 10, 12, 14]
BEGIN_POS = [32961690, 44559528, 108999785, 15473140, 2474080, 36694953]
END_POS = [32961691, 44559529, 108999786, 15473141, 2474081, 36694954]
REFS = ['C', 'C', 'T', 'G', 'T', 'C']
ALTS = [['T'], ['T'], ['C'], ['T'], ['C'], ['G']]
QUALS = [89.59, 234.97, 73.25, 242.29, 150.25, 471.66]
GENOTYPES = ['C/T', 'C/T', 'T/C', 'G/T', 'T/C', 'C/G']

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def check_file(f):
    records = []
    for i, record in enumerate(f):
        assert record.r_id == RIDS[i]
        assert record.chrom == CHROMS[i]
        assert record.begin_pos == BEGIN_POS[i]
        assert record.end_pos == END_POS[i]
        assert record.ids == IDS[i]
        assert record.ref == REFS[i]
        assert record.alts == ALTS[i]
        assert abs(record.qual - QUALS[i]) < 0.01
        assert record.filters == ['PASS']
        assert len(record.info.keys()) >= 15
        KEYS = ['AC', 'AF', 'AN', 'BaseQRankSum', 'DP', 'Dels', 'FS',
                'HRun', 'HaplotypeScore', 'MQ0', 'MQ', 'MQRankSum', 'QD',
                'ReadPosRankSum', 'set']
        for key in KEYS:
            assert key in record.info.keys()
        assert record.format == ['GT', 'AD', 'DP', 'GQ', 'PL']
        assert len(record.genotypes) == 1
        assert str(record.genotypes[0].gt) == \
            'GenotypeCall({})'.format(repr(GENOTYPES[i]))
        records.append(record.detach())

    assert len(records) == 6

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_six_records_read_sequential_vcf(six_records_vcf):
    with bcf.BCFFile(str(six_records_vcf)) as f:
        check_file(f)


def test_six_records_read_sequential_vcf_gz(six_records_vcf_gz):
    with bcf.BCFFile(str(six_records_vcf_gz)) as f:
        check_file(f)


def test_six_records_read_sequential_bcf(six_records_bcf):
    with bcf.BCFFile(str(six_records_bcf)) as f:
        check_file(f)


def test_two_hundred_read_sequential_vcf(two_hundred_vcf):
    with bcf.BCFFile(str(two_hundred_vcf)) as f:
        num = len(list(f))
        assert num == 200


def test_two_hundred_read_sequential_vcf_gz(two_hundred_vcf_gz):
    with bcf.BCFFile(str(two_hundred_vcf_gz)) as f:
        num = len(list(f))
        assert num == 200


def test_two_hundred_read_sequential_bcf(two_hundred_bcf):
    with bcf.BCFFile(str(two_hundred_bcf)) as f:
        num = len(list(f))
        assert num == 200


def test_two_hundread_through_index_vcf_gz(
        two_hundred_vcf_gz, two_hundred_tbi):
    with bcf.BCFIndex(str(two_hundred_vcf_gz)) as idx:
        assert len(list(idx.query('17:10,000,000-11,000,000'))) == 2


def test_two_hundread_through_index_bcf(two_hundred_bcf, two_hundred_csi):
    with bcf.BCFIndex(str(two_hundred_bcf)) as idx:
        assert len(list(idx.query('17:10,000,000-11,000,000'))) == 2
        assert len(list(idx.query('17:10,000,000-15,000,000'))) == 3
