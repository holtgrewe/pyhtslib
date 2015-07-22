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


def check_header_of_header_only_vcf(header):
    assert header.struct_ptr
    assert header.struct
    assert header.target_infos == [
        bcf.BCFHeaderTargetInfo('CHROMOSOME_I', 1009800),
        bcf.BCFHeaderTargetInfo('CHROMOSOME_II', 5000),
        bcf.BCFHeaderTargetInfo('CHROMOSOME_III', 5000),
        bcf.BCFHeaderTargetInfo('CHROMOSOME_IV', 5000),
        bcf.BCFHeaderTargetInfo('CHROMOSOME_V', 5000),
        bcf.BCFHeaderTargetInfo('CHROMOSOME_X', 5000),
        bcf.BCFHeaderTargetInfo('CHROMOSOME_MtDNA', 5000)]
    assert header.header_records == []

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
    return
    with bcf.BCFFile(str(header_only_vcf)) as f:
        check_header_of_header_only_vcf(f.header)


def test_read_header_from_vcf_gz(header_only_vcf_gz):
    return
    with bcf.BCFFile(str(header_only_vcf_gz)) as f:
        check_header_of_header_only_vcf(f.header)


def test_read_header_from_bcf(header_only_bcf):
    return
    with bcf.BCFFile(str(header_only_bcf)) as f:
        check_header_of_header_only_vcf(f.header)
