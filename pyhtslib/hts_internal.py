#!/usr/bin/env python3
"""Shared htslib structs and functions"""

import ctypes

import pyhtslib.load_dll as pl

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'


# export everything from this submodule manually, including the code that
# starts with an underscore, importing modules will not import the latter
__all__ = [
    # handle to htslib
    'htslib',
    '_libc',
    # constants (through ``#define``)
    '_HTS_IDX_NOCOOR',
    '_HTS_IDX_START',
    '_HTS_IDX_REST',
    '_HTS_IDX_NONE',
    # htslib types
    '_BGZF',
    '_cram_fd',
    '_hFILE',
    '_htsFormat',
    '_htsFile',
    '_hts_idx_t',
    '_hts_itr_t',
    # klib types
    '_kstring_t',
    '_KS_SEP_SPACE',
    '_KS_SEP_TAB',
    '_KS_SEP_LINE',
    '_KS_SEP_MAX',
    # htslib functions
    '_bgzf_is_bgzf',
    '_hts_open',
    '_hts_close',
    '_hts_getline',
    '_hts_itr_destroy',
    '_hts_itr_next',
    '_hts_itr_query',
    '_hts_itr_querys',
    '_tbx_readrec',
    # wrapper Types
    'HTSFormatCategory',
    'HTSExactFormat',
    'HTSCompression',
]


htslib = pl.load_htslib()
_libc = pl.load_libc()

_bgzf_is_bgzf = htslib.bgzf_is_bgzf
_bgzf_is_bgzf.restype = ctypes.c_int


_HTS_IDX_NOCOOR = -2
_HTS_IDX_START = -3
_HTS_IDX_REST = -4
_HTS_IDX_NONE = -5


class _BGZF(ctypes.Structure):
    """Type for representing a bgzip-compressed file"""


class _cram_fd(ctypes.Structure):
    pass


class _hFILE(ctypes.Structure):
    pass


_KS_SEP_SPACE = 0  # isspace(): \t, \n, \v, \f, \r
_KS_SEP_TAB = 1  # isspace() && !' '
_KS_SEP_LINE = 2  # line separator: "\n" (Unix) or "\r\n" (Windows)
_KS_SEP_MAX = 2


class _kstring_t(ctypes.Structure):
    """Type from klib for string representation"""

    # we use ``ctypes.c_void_p`` for the type of s here as we do not
    # want to have ctype's behaviour of casting to ``bytes`` and losing
    # the original pointer
    _fields_ = [('l', ctypes.c_size_t),
                ('m', ctypes.c_size_t),
                ('p', ctypes.c_void_p)]

    @property
    def s(self):
        return ctypes.cast(self.p, ctypes.c_char_p).value

    def free_p(self):
        _libc.free(self.p)


class HTSFormatCategory:
    """Wrapper for the htslib type `htsFormatCategory`."""

    UNKNOWN_CATEGORY = 0
    SEQUENCE_DATA = 1
    VARIANT_DATA = 2
    INDEX_FILE = 3
    REGION_LIST = 4
    CATEGORY_MAXIMUM = 32767

    @staticmethod
    def to_str(value):
        DICT = {
            0: 'UNKNOWN_CATEGORY',
            1: 'SEQUENCE_DATA',
            2: 'VARIANT_DATA',
            3: 'INDEX_FILE',
            4: 'REGION_LIST',
            32767: 'CATEGORY_MAXIMUM',
        }
        return DICT[value]


class HTSExactFormat:
    """Wrapper for the htslib type `htsExactFormat`."""

    UNKNOWN_FORMAT = 0
    BINARY_FORMAT = 1
    TEXT_FORMAT = 2
    SAM = 3
    BAM = 4
    BAI = 5
    CRAM = 6
    CRAI = 7
    VCF = 8
    BCF = 9
    CSI = 10
    GZI = 11
    TBI = 12
    BED = 13
    FORMAT_MAXIMUM = 32767

    @staticmethod
    def to_str(value):
        DICT = {
            0: 'UNKNOWN_FORMAT',
            1: 'BINARY_FORMAT',
            2: 'TEXT_FORMAT',
            3: 'SAM',
            4: 'BAM',
            5: 'BAI',
            6: 'CRAM',
            7: 'CRAI',
            8: 'VCF',
            9: 'BCF',
            10: 'CSI',
            11: 'GZI',
            12: 'TBI',
            13: 'BED',
            32767: 'FORMAT_MAXIMUM',
        }
        return DICT[value]


class HTSCompression:
    """Wrapper for the htslib type htsExactFormat."""

    NO_COMPRESSION = 0
    GZIP = 1
    BGZF = 2
    CUSTOM = 3
    COMPRESSION_MAXIMUM = 32767

    @staticmethod
    def to_str(value):
        DICT = {
            0: 'NO_COMPRESSION',
            1: 'GZIP',
            2: 'BGZF',
            3: 'CUSTOM',
            3276: 'COMPRESSION_MAXIUM',
        }
        return DICT[value]


class _htsFormat_version_t(ctypes.Structure):
    """Wrapper for htsFormat's anonymous struct for version type."""

    _fields_ = [('major', ctypes.c_short),
                ('minor', ctypes.c_short)]


class _htsFormat(ctypes.Structure):
    """Wrapper for the htslib type htsFormat."""

    _fields_ = [('category', ctypes.c_uint),
                ('format', ctypes.c_uint),
                ('version', _htsFormat_version_t),
                ('compression', ctypes.c_uint),
                ('compression_level', ctypes.c_short),
                ('specific', ctypes.c_void_p)]


class _htsFile_fp_member(ctypes.Union):
    """Union type for the "fp" member of HTSFileStrut"""

    _fields_ = [('bgzf', ctypes.POINTER(_BGZF)),
                ('cram', ctypes.POINTER(_cram_fd)),
                ('hfile', ctypes.POINTER(_hFILE)),
                ('voidp', ctypes.c_void_p)]


class _htsFile(ctypes.Structure):
    """Wrapper for htslib type htsFile."""

    _fields_ = [('is_bin', ctypes.c_uint32, 1),
                ('is_write', ctypes.c_uint32, 1),
                ('is_be', ctypes.c_uint32, 1),
                ('is_cram', ctypes.c_uint32, 1),
                ('dummy', ctypes.c_uint32, 28),
                ('lineno', ctypes.c_int64),
                ('line', _kstring_t),
                ('fn', ctypes.c_char_p),
                ('fn_aux', ctypes.c_char_p),
                ('fp', _htsFile_fp_member),
                ('ftype', _htsFormat)]


class _hts_idx_t(ctypes.Structure):
    """Generic type for htslib index types"""


class _hts_itr_t(ctypes.Structure):
    """Generic type for htslib iterator types"""


_hts_open = htslib.hts_open
_hts_open.restype = ctypes.POINTER(_htsFile)

_hts_close = htslib.hts_close
_hts_close.restype = ctypes.c_int

_hts_getline = htslib.hts_getline
_hts_getline.restype = ctypes.c_int

_hts_itr_destroy = htslib.hts_itr_destroy
_hts_itr_destroy.restype = None

_hts_itr_next = htslib.hts_itr_next
_hts_itr_next.restype = ctypes.c_int

_hts_itr_query = htslib.hts_itr_query
_hts_itr_query.restype = ctypes.POINTER(_hts_itr_t)

_hts_itr_querys = htslib.hts_itr_querys
_hts_itr_querys.restype = ctypes.POINTER(_hts_itr_t)


_tbx_name2id = htslib.tbx_name2id
_tbx_name2id.restype = ctypes.c_int


_tbx_readrec = htslib.tbx_readrec
_tbx_readrec.restype = ctypes.c_int
