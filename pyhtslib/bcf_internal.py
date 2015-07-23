#!/usr/bin/env python3
"""Implementation details for VCF/BCF access"""

import ctypes

import pyhtslib.load_dll as pl
from pyhtslib.hts_internal import *  # NOQA

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'

# export everything from this submodule manually, including the code that
# starts with an underscore, importing modules will not import the latter
__all__ = [
    # constants

    '_BCF_HL_FLT',
    '_BCF_HL_INFO',
    '_BCF_HL_FMT',
    '_BCF_HL_CTG',
    '_BCF_HL_STR',
    '_BCF_HL_GEN',

    '_BCF_HT_FLAG',
    '_BCF_HT_INT',
    '_BCF_HT_REAL',
    '_BCF_HT_STR',

    '_BCF_VL_FIXED',
    '_BCF_VL_VAR',
    '_BCF_VL_A',
    '_BCF_VL_G',
    '_BCF_VL_R',

    '_BCF_DT_ID',
    '_BCF_DT_CTG',
    '_BCF_DT_SAMPLE',

    '_UNKNOWN_CATEGORY',
    '_SEQUENCE_DATA',
    '_VARIANT_DATA',
    '_INDEX_FILE',
    '_REGION_LIST',
    '_CATEGORY_MAXIMUM',

    '_UNKNOWN_FORMAT',
    '_BINARY_FORMAT',
    '_TEXT_FORMAT',
    '_SAM',
    '_BAM',
    '_BAI',
    '_CRAM',
    '_CRAI',
    '_VCF',
    '_BCF',
    '_CSI',
    '_GZI',
    '_TBI',
    '_BED',
    '_FORMAT_MAXIMUM',

    '_NO_COMPRESSION',
    '_GZIP',
    '_BGZF',
    '_CUSTOM',
    '_COMPRESSION_MAXIMUM',

    '_BCF_BT_NULL',
    '_BCF_BT_INT8',
    '_BCF_BT_INT16',
    '_BCF_BT_INT32',
    '_BCF_BT_FLOAT',
    '_BCF_BT_CHAR',

    '_VCF_REF',
    '_VCF_SNP',
    '_VCF_MNP',
    '_VCF_INDEL',
    '_VCF_OTHER',

    '_BCF1_DIRTY_ID',
    '_BCF1_DIRTY_ALS',
    '_BCF1_DIRTY_FLT',
    '_BCF1_DIRTY_INF',

    '_BCF_ERR_CTG_UNDEF',
    '_BCF_ERR_TAG_UNDEF',
    '_BCF_ERR_NCOLS',
    '_BCF_ERR_LIMITS',

    '_BCF_UN_STR',
    '_BCF_UN_FLT',
    '_BCF_UN_INFO',
    '_BCF_UN_SHR',
    '_BCF_UN_FMT',
    '_BCF_UN_IND',
    '_BCF_UN_ALL',

    # structures

    '_bcf_hrec_t',
    '_bcf_idinfo_t',
    '_bcf_hdr_t',
    '_variant_t',

    # routines

    '_bcf_fmt_t',
    '_bcf_info_t_v1_member',
    '_bcf_info_t',
    '_bcf_dec_t',
    '_bcf1_t',
    '_bcf_hdr_name2id',
    '_bcf_readrec',
    '_bcf_hdr_read',
    '_bcf_alt_hdr_read',
    '_bcf_read',
    '_bcf_read1',
    '_bcf_hdr_get_version',
    '_bcf_itr_next',
    '_bcf_itr_destroy',
    '_bcf_init',
    '_bcf_init1',
    '_bcf_destroy1',
    '_bcf_hdr_seqnames',
    '_bcf_hdr_id2int',
    '_bcf_hdr_id2length',
    '_bcf_hdr_id2number',
    '_bcf_hdr_id2type',
    '_bcf_hdr_id2coltype',
    '_bcf_hdr_info_exists',
    '_bcf_hdr_id2hrec',
    '_bcf_unpack',
    '_bcf_get_info_values',
    '_bcf_get_info_int32',
    '_bcf_get_info_float',
    '_bcf_get_info_string',
    '_bcf_get_info_flag',
    '_bcf_hdr_int2id',
    '_bcf_hdr_nids',
    '_bcf_hdr_nsequences',
    '_bcf_hdr_nsamples',
    '_bcf_hdr_get_sample_name',
    '_bcf_hdr_destroy',
    '_bcf_rec_nalleles',
    '_bcf_itr_querys',
    '_bcf_index_load',
    '_bcf_index_load2',

    '_vcf_read1',
    '_vcf_read',
    '_vcf_parse1',
    '_vcf_parse',

    '_bcf_get_fmt',
    '_bcf_get_info',
    '_bcf_get_fmt_id',
    '_bcf_get_info_id',

    '_bcf_get_format_int32',
    '_bcf_get_format_float',
    '_bcf_get_format_char',
    '_bcf_get_genotypes',
    '_bcf_get_format_string',
    '_bcf_get_format_values',

    '_bcf_gt_phased',
    '_bcf_gt_unphased',
    '_bcf_gt_missing',
    '_bcf_gt_is_missing',
    '_bcf_gt_is_phased',
    '_bcf_gt_allele',
]

# ----------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------

_BCF_HL_FLT = 0  # header line
_BCF_HL_INFO = 1
_BCF_HL_FMT = 2
_BCF_HL_CTG = 3
_BCF_HL_STR = 4  # structured header line TAG=<A=..,B=..>
_BCF_HL_GEN = 5  # generic header line

_BCF_HT_FLAG = 0  # header type
_BCF_HT_INT = 1
_BCF_HT_REAL = 2
_BCF_HT_STR = 3

_BCF_VL_FIXED = 0  # variable length
_BCF_VL_VAR = 1
_BCF_VL_A = 2
_BCF_VL_G = 3
_BCF_VL_R = 4

_BCF_DT_ID = 0  # dictionary type
_BCF_DT_CTG = 1
_BCF_DT_SAMPLE = 2

_UNKNOWN_CATEGORY = 0
_SEQUENCE_DATA = 1
_VARIANT_DATA = 2
_INDEX_FILE = 3
_REGION_LIST = 4
_CATEGORY_MAXIMUM = 32767

_UNKNOWN_FORMAT = 0
_BINARY_FORMAT = 1
_TEXT_FORMAT = 2
_SAM = 3
_BAM = 4
_BAI = 5
_CRAM = 6
_CRAI = 7
_VCF = 8
_BCF = 9
_CSI = 10
_GZI = 11
_TBI = 12
_BED = 13
_FORMAT_MAXIMUM = 32767

_NO_COMPRESSION = 0
_GZIP = 1
_BGZF = 2
_CUSTOM = 3
_COMPRESSION_MAXIMUM = 32767

_BCF_BT_NULL = 0
_BCF_BT_INT8 = 1
_BCF_BT_INT16 = 2
_BCF_BT_INT32 = 3
_BCF_BT_FLOAT = 5
_BCF_BT_CHAR = 7

_VCF_REF = 0
_VCF_SNP = 1
_VCF_MNP = 2
_VCF_INDEL = 4
_VCF_OTHER = 8

_BCF1_DIRTY_ID = 1
_BCF1_DIRTY_ALS = 2
_BCF1_DIRTY_FLT = 4
_BCF1_DIRTY_INF = 8

_BCF_ERR_CTG_UNDEF = 1
_BCF_ERR_TAG_UNDEF = 2
_BCF_ERR_NCOLS = 4
_BCF_ERR_LIMITS = 8

_BCF_UN_STR = 1  # up to ALT inclusive
_BCF_UN_FLT = 2  # up to FILTER
_BCF_UN_INFO = 4  # up to INFO
_BCF_UN_SHR = (_BCF_UN_STR | _BCF_UN_FLT | _BCF_UN_INFO)  # all shared
_BCF_UN_FMT = 8  # unpack format and each sample
_BCF_UN_IND = _BCF_UN_FMT  # a synonymo of _BCF_UN_FMT
_BCF_UN_ALL = (_BCF_UN_SHR | _BCF_UN_FMT)  # everything

# ----------------------------------------------------------------------------
# Structures
# ----------------------------------------------------------------------------


class _bcf_hrec_t(ctypes.Structure):

    _fields_ = [('type', ctypes.c_int),
                ('key', ctypes.c_char_p),
                ('value', ctypes.c_char_p),
                ('nkeys', ctypes.c_int),
                ('keys', ctypes.POINTER(ctypes.c_char_p)),
                ('vals', ctypes.POINTER(ctypes.c_char_p))]


class _bcf_idinfo_t(ctypes.Structure):

    _fields_ = [('info', ctypes.c_uint32 * 3),
                ('hrec', _bcf_hrec_t * 3),
                ('id', ctypes.c_int)]


class _bcf_idpair_t(ctypes.Structure):

    _fields_ = [('key', ctypes.c_char_p),
                ('val', ctypes.POINTER(_bcf_idinfo_t))]


class _bcf_hdr_t(ctypes.Structure):

    _fields_ = [('n', ctypes.c_int32 * 3),
                ('id', ctypes.POINTER(_bcf_idpair_t) * 3),
                ('dict', ctypes.c_void_p * 3),
                ('samples', ctypes.POINTER(ctypes.c_char_p)),
                ('hrec', ctypes.POINTER(ctypes.POINTER(_bcf_hrec_t))),
                ('nhrec', ctypes.c_int),
                ('dirty', ctypes.c_int),
                ('ntransl', ctypes.c_int),
                ('transl', ctypes.POINTER(ctypes.c_int) * 2),
                ('nsamples_ori', ctypes.c_int),
                ('keep_samples', ctypes.POINTER(ctypes.c_uint8)),
                ('mem', _kstring_t)]


class _variant_t(ctypes.Structure):

    _fields_ = [('type', ctypes.c_int),
                ('n', ctypes.c_int)]


class _bcf_fmt_t(ctypes.Structure):

    _fields_ = [('id', ctypes.c_int),
                ('n', ctypes.c_int),
                ('size', ctypes.c_int),
                ('type', ctypes.c_int),
                ('p', ctypes.POINTER(ctypes.c_uint8)),
                ('p_len', ctypes.c_uint32),
                ('p_off', ctypes.c_uint32, 31),
                ('p_free', ctypes.c_uint32, 1)]


class _bcf_info_t_v1_member(ctypes.Union):

    _fields_ = [('i', ctypes.c_int32),
                ('f', ctypes.c_float)]


class _bcf_info_t(ctypes.Structure):

    _fields_ = [('key', ctypes.c_int),
                ('type', ctypes.c_int),
                ('len', ctypes.c_int),
                ('v1', _bcf_info_t_v1_member),
                ('vptr', ctypes.POINTER(ctypes.c_uint8)),
                ('vptr_len', ctypes.c_uint32),
                ('vptr_off', ctypes.c_uint32, 31),
                ('vptr_free', ctypes.c_uint32, 1)]


class _bcf_dec_t(ctypes.Structure):

    _fields_ = [('m_fmt', ctypes.c_int),
                ('m_info', ctypes.c_int),
                ('m_id', ctypes.c_int),
                ('m_als', ctypes.c_int),
                ('m_allele', ctypes.c_int),
                ('m_flt', ctypes.c_int),
                ('n_flt', ctypes.c_int),
                ('flt', ctypes.POINTER(ctypes.c_int)),
                ('id', ctypes.c_char_p),
                ('als', ctypes.c_char_p),
                ('allele', ctypes.POINTER(ctypes.c_char_p)),
                ('info', ctypes.POINTER(_bcf_info_t)),
                ('fmt', ctypes.POINTER(_bcf_fmt_t)),
                ('var', ctypes.POINTER(_variant_t)),
                ('n_var', ctypes.c_int),
                ('n_var_type', ctypes.c_int),
                ('n_shared_dirty', ctypes.c_int),
                ('n_indiv_dirty', ctypes.c_int)]


class _bcf1_t(ctypes.Structure):

    _fields_ = [('rid', ctypes.c_int32),
                ('pos', ctypes.c_int32),
                ('rlen', ctypes.c_int32),
                ('qual', ctypes.c_float),
                ('n_info', ctypes.c_uint32, 16),
                ('n_allele', ctypes.c_uint32, 16),
                ('n_fmt', ctypes.c_uint32, 8),
                ('n_sample', ctypes.c_uint32, 24),
                ('shared', _kstring_t),
                ('indiv', _kstring_t),
                ('d', _bcf_dec_t),
                ('max_unpack', ctypes.c_int),
                ('unpacked', ctypes.c_int),
                ('unpack_size', ctypes.c_int * 3),
                ('errcode', ctypes.c_int)]

# ----------------------------------------------------------------------------
# C functions and their return types
# ----------------------------------------------------------------------------

htslib = pl.load_htslib()

# put into .so by us, was static inline
_bcf_hdr_name2id = htslib.pyhtslib_bcf_hdr_name2id
_bcf_hdr_name2id.restype = ctypes.c_int

_bcf_readrec = htslib.bcf_readrec
_bcf_readrec.restype = ctypes.c_int

_bcf_hdr_read = htslib.bcf_hdr_read
_bcf_hdr_read.restype = ctypes.POINTER(_bcf_hdr_t)

# our extension, from vt
_bcf_alt_hdr_read = htslib.bcf_alt_hdr_read
_bcf_alt_hdr_read.restype = ctypes.POINTER(_bcf_hdr_t)

_bcf_read = htslib.bcf_read
_bcf_read.restype = ctypes.c_int


def _bcf_read1(fp, h, v):
    return _bcf_read(fp, h, v)

_bcf_hdr_get_version = htslib.bcf_hdr_get_version
_bcf_hdr_get_version.restype = ctypes.c_char_p


def _bcf_itr_next(htsfp, itr, r):
    """Replacement for C macro ``bcf_itr_next``."""
    return _hts_itr_next(htsfp[0].fp.bgzf, itr, r, 0)


def _bcf_itr_destroy(iter):
    return _hts_itr_destroy(iter)

_bcf_init = htslib.bcf_init
_bcf_init.restype = ctypes.POINTER(_bcf1_t)

_bcf_destroy = htslib.bcf_destroy
_bcf_destroy.restype = None

_bcf_hdr_seqnames = htslib.bcf_hdr_seqnames
_bcf_hdr_seqnames.restype = ctypes.POINTER(ctypes.c_char_p)

_bcf_hdr_id2int = htslib.bcf_hdr_id2int
_bcf_hdr_id2int.restype = ctypes.c_int

_bcf_hdr_destroy = htslib.bcf_hdr_destroy
_bcf_hdr_destroy.restype = None


def _bcf_hdr_id2length(hdr, type_, int_id):
    return (hdr[0].id[_BCF_DT_ID][int_id].val[0].info[type_] >> 8) & 0xf


def _bcf_hdr_id2number(hdr, type_, int_id):
    return (hdr[0].id[_BCF_DT_ID][int_id].val[0].info[type_] >> 12)


def _bcf_hdr_id2type(hdr, type_, int_id):
    return (hdr[0].id[_BCF_DT_ID][int_id].val[0].info[type_] >> 4) & 0xf


def _bcf_hdr_id2coltype(hdr, type_, int_id):
    return hdr[0].id[_BCF_DT_ID][int_id].val[0].info[type_] & 0xf


def _bcf_hdr_info_exists(hdr, type_, int_id):
    return not (int_id < 0 or (_bcf_hdr_id2coltype(hdr, type_, int_id) == 0xf))


def _bcf_hdr_id2hrec(hdr, dict_type, col_type, int_id):
    return hdr[0].id[_BCF_DT_CTG if (dict_type == _BCF_DT_CTG)
                     else _BCF_DT_ID][0].hrec[
                         0 if (dict_type == _BCF_DT_CTG) else col_type]


_bcf_unpack = htslib.bcf_unpack
_bcf_unpack.restype = ctypes.c_int

_bcf_get_info_values = htslib.bcf_get_info_values
_bcf_get_info_values.restype = ctypes.c_int


def _bcf_get_info_int32(hdr, line, tag, dst, ndst):
    return _bcf_get_info_values(hdr, line, tag, dst, ndst, _BCF_HT_INT)


def _bcf_get_info_float(hdr, line, tag, dst, ndst):
    return _bcf_get_info_values(hdr, line, tag, dst, ndst, _BCF_HT_REAL)


def _bcf_get_info_string(hdr, line, tag, dst, ndst):
    return _bcf_get_info_values(hdr, line, tag, dst, ndst, _BCF_HT_STR)


def _bcf_get_info_flag(hdr, line, tag, dst, ndst):
    return _bcf_get_info_values(hdr, line, tag, dst, ndst, _BCF_HT_FLAG)


def _bcf_hdr_int2id(hdr, type_, int_id):
    return hdr[0].id[type_][int_id].key


def _bcf_hdr_nids(hdr):
    """Helper function that returns the number of IDs"""
    return hdr[0].n[_BCF_DT_ID]


def _bcf_hdr_nsequences(hdr):
    """Helper function that returns the number of sequences"""
    return hdr[0].n[_BCF_DT_CTG]


def _bcf_hdr_nsamples(hdr):
    """Helper function that returns the number of sequences
    Replacement for C macro ``bcf_hdr_nsamples``.
    """
    return hdr[0].n[_BCF_DT_SAMPLE]


def _bcf_hdr_get_sample_name(hdr, i):
    """Replacement for C macro ``bcf_hdr_get_sample_name`` from vt."""
    return hdr[0].samples[i]


def _bcf_rec_nalleles(rec_ptr):
    """Helper function for obtaining number of alleles from a bcf1_t record."""
    return rec_ptr[0].n_allele


def _bcf_init1():
    """Replacement for C macro ``bcf_init``."""
    return _bcf_init()


def _bcf_destroy1(ptr):
    """Replacement for C macro ``bcf_destroy1``."""
    return _bcf_destroy(ptr)


def _bcf_itr_querys(idx, hdr, s):
    """Replacement for the C macro ``tbx_itr_querys()``"""
    return _hts_itr_querys(idx, s, _bcf_hdr_name2id, hdr, _hts_itr_query,
                           _bcf_readrec)


def _bcf_index_load(fn):
    return _hts_idx_load(fn, HTS_FMT_CSI)


_bcf_index_load2 = htslib.bcf_index_load2
_bcf_index_load2.restype = ctypes.POINTER(_hts_idx_t)

_vcf_read = htslib.vcf_read
_vcf_read.restype = ctypes.c_int


def _vcf_read1(p, h, v):
    return _vcf_read(p, h, v)

_vcf_parse = htslib.vcf_parse
_vcf_parse.restype = ctypes.c_int


def _vcf_parse1(p, h, v):
    return _vcf_parse(p, h, v)

_bcf_get_fmt = htslib.bcf_get_fmt
_bcf_get_fmt.restype = ctypes.POINTER(_bcf_fmt_t)

_bcf_get_info = htslib.bcf_get_info
_bcf_get_info.restype = ctypes.POINTER(_bcf_info_t)

_bcf_get_fmt_id = htslib.bcf_get_fmt_id
_bcf_get_fmt_id.restype = ctypes.POINTER(_bcf_info_t)

_bcf_get_info_id = htslib.bcf_get_info_id
_bcf_get_info_id.restype = ctypes.POINTER(_bcf_info_t)

_bcf_get_format_string = htslib.bcf_get_format_string
_bcf_get_format_string.restype = ctypes.c_int

_bcf_get_format_values = htslib.bcf_get_format_values
_bcf_get_format_values.restype = ctypes.c_int


def _bcf_get_format_int32(hdr, line, tag, dst, ndst):
    return _bcf_get_format_values(hdr, line, tag, dst, ndst, _BCF_HT_INT)


def _bcf_get_format_float(hdr, line, tag, dst, ndst):
    return _bcf_get_format_values(hdr, line, tag, dst, ndst, _BCF_HT_REAL)


def _bcf_get_format_char(hdr, line, tag, dst, ndst):
    return _bcf_get_format_values(hdr, line, tag, dst, ndst, _BCF_HT_STR)


def _bcf_get_genotypes(hdr, line, dst, ndst):
    return _bcf_get_format_values(hdr, line, "GT".encode('utf-8'),
                                  dst, ndst, _BCF_HT_INT)


def _bcf_gt_phased(idx):
    return ((idx + 1) << 1 | 1)


def _bcf_gt_unphased(idx):
    return ((idx + 1) << 1)


def _bcf_gt_missing():
    return 0


def _bcf_gt_is_missing(val):
    return (0 if (val >> 1) else 1)


def _bcf_gt_is_phased(idx):
    return (idx & 1)


def _bcf_gt_allele(val):
    return (((val) >> 1) - 1)
