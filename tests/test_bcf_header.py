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


def check_header_of_header_only_vcf(header, line_add=0):
    # conversion to BCF appends two header lines
    assert header.struct_ptr
    assert header.struct
    assert header.vcf_version == 'VCFv4.2'
    assert header.sample_names == ['sample_normal_wes', 'sample_tumor_wes']
    assert len(header.target_infos) == 85
    assert header.target_infos[0] == bcf.BCFHeaderTargetInfo('1', 249250621)
    assert header.target_infos[-1] == \
        bcf.BCFHeaderTargetInfo('NC_007605', 171823)
    assert header.key_ids == \
        {'AD': 2, 'BQ': 3, 'DB': 10, 'DP': 4, 'FA': 5, 'GQ': 6, 'GT': 7,
         'MQ0': 11, 'PASS': 0, 'PL': 8, 'REJECT': 1, 'SOMATIC': 12,
         'SS': 9, 'VT': 13}
    assert len(header.header_records) == 101 + line_add
    assert set(header.id_to_filter_record.keys()) == set(['PASS', 'REJECT'])
    assert set(header.id_to_info_record.keys()) == \
        set(['DB', 'MQ0', 'SOMATIC', 'VT'])
    assert set(header.id_to_format_record) == \
        set(['AD', 'BQ', 'DP', 'FA', 'GQ', 'GT', 'PL', 'SS'])
    assert len(header.id_to_contig_record) == 85


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_bcf_header_target_info():
    info = bcf.BCFHeaderTargetInfo('name', 33)
    assert info.name == 'name'
    assert info.length == 33
    assert info == bcf.BCFHeaderTargetInfo('name', 33)


def test_bcf_generic_header_record():
    info = bcf.BCFGenericHeaderRecord('name', 'text')
    assert info.tag == 'name'
    assert info.value == 'text'
    assert info == bcf.BCFGenericHeaderRecord('name', 'text')
    assert info.to_header_line() == '##name=text'


def test_bcf_structured_header_record():
    info = bcf.BCFStructuredHeaderRecord(
        'name', [('key', 'value'), ('key2', 'value2')], ['key'])
    assert info.tag == 'name'
    assert info.entries == OrderedDict([('key', 'value'), ('key2', 'value2')])
    assert info.repr_entries == set(['key'])
    assert info == bcf.BCFStructuredHeaderRecord(
        'name', [('key', 'value'), ('key2', 'value2')], ['key'])
    assert info.to_header_line() == \
        '##name=<key="value",key2=value2>'


def test_bcf_filter_header_record():
    info = bcf.BCFFilterHeaderRecord(
        'FILTER', 'q50', 'quality', [('key', 'value'), ('key2', 'value2')],
        ['key'])
    assert info.tag == 'FILTER'
    assert info.entries == OrderedDict([
        ('id', 'q50'), ('description', 'quality'),
        ('key', 'value'), ('key2', 'value2')])
    assert info.repr_entries == set(['description', 'key'])
    assert info == bcf.BCFFilterHeaderRecord(
        'FILTER', 'q50', 'quality',
        [('key', 'value'), ('key2', 'value2')], ['key'])
    assert info.to_header_line() == \
        '##FILTER=<id=q50,description="quality",key="value",key2=value2>'


def test_bcf_info_header_record():
    info = bcf.BCFInfoHeaderRecord(
        'INFO', 'Foo', 1, 'String', 'My description',
        [('key', 'value'), ('key2', 'value2')], ['key'])
    assert info.tag == 'INFO'
    assert info.entries == OrderedDict([
        ('id', 'Foo'), ('number', 1), ('type', 'String'),
        ('description', 'My description'), ('key', 'value'),
        ('key2', 'value2')])
    assert info.repr_entries == set(['description', 'key'])
    assert info == bcf.BCFInfoHeaderRecord(
        'INFO', 'Foo', 1, 'String', 'My description',
        [('key', 'value'), ('key2', 'value2')], ['key'])
    assert info.to_header_line() == \
        ('##INFO=<id=Foo,number=1,type=String,'
         'description="My description",key="value",key2=value2>')


def test_bcf_format_header_record():
    info = bcf.BCFFormatHeaderRecord(
        'FORMAT', 'Foo', 1, 'String', 'My description',
        [('key', 'value'), ('key2', 'value2')], ['key'])
    assert info.tag == 'FORMAT'
    assert info.entries == OrderedDict([
        ('id', 'Foo'), ('number', 1), ('type', 'String'),
        ('description', 'My description'), ('key', 'value'),
        ('key2', 'value2')])
    assert info.repr_entries == set(['description', 'key'])
    assert info == bcf.BCFFormatHeaderRecord(
        'FORMAT', 'Foo', 1, 'String', 'My description',
        [('key', 'value'), ('key2', 'value2')], ['key'])
    assert info.to_header_line() == \
        ('##FORMAT=<id=Foo,number=1,type=String,'
         'description="My description",key="value",key2=value2>')


def test_bcf_contig_header_record():
    info = bcf.BCFContigHeaderRecord(
        'contig', 'MT', 12345, [('key', 'value'), ('key2', 'value2')],
        ['key'])
    assert info.tag == 'contig'
    assert info.entries == OrderedDict([('id', 'MT'), ('length', 12345),
                                        ('key', 'value'), ('key2', 'value2')])
    assert info.repr_entries == set(['key'])
    assert info == bcf.BCFContigHeaderRecord(
        'contig', 'MT', 12345, [('key', 'value'), ('key2', 'value2')],
        ['key'])
    assert info.to_header_line() == \
        '##contig=<id=MT,length=12345,key="value",key2=value2>'


def test_bcf_record_construct_minimal():
    """Construt a minimal bcf.BCFHeader from scratch"""
    header = bcf.BCFHeader()
    assert header.to_vcf_header() == '##fileformat=VCFv4.2\n'


def test_bcf_record_construct_small():
    """Construt a small bcf.BCFHeader from scratch"""
    header = bcf.BCFHeader()
    header.add_header_record(bcf.BCFFilterHeaderRecord(
        'FILTER', 'q50', 'Quality'))
    header.add_header_record(bcf.BCFInfoHeaderRecord(
        'INFO', 'Bar', 1, 'String', 'Some info'))
    header.add_header_record(bcf.BCFFormatHeaderRecord(
        'FORMAT', 'Foo', 1, 'String', 'Some format'))
    header.add_header_record(bcf.BCFContigHeaderRecord('contig', 'MT', 12345))
    assert header.to_vcf_header() == \
        ('##fileformat=VCFv4.2\n'
         '##FILTER=<id=q50,description="Quality">\n'
         '##INFO=<id=Bar,number=1,type=String,description="Some info">\n'
         '##FORMAT=<id=Foo,number=1,type=String,description="Some format">\n'
         '##contig=<id=MT,length=12345>\n')


def test_read_header_from_vcf(header_only_vcf):
    with bcf.BCFFile(str(header_only_vcf)) as f:
        check_header_of_header_only_vcf(f.header)


def test_read_header_from_vcf_gz(header_only_vcf_gz):
    with bcf.BCFFile(str(header_only_vcf_gz)) as f:
        check_header_of_header_only_vcf(f.header)


def test_read_header_from_bcf(header_only_bcf):
    with bcf.BCFFile(str(header_only_bcf)) as f:
        check_header_of_header_only_vcf(f.header, 2)
