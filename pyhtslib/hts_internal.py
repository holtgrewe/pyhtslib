#!/usr/bin/env python3
"""Shared htslib structs and functions"""

import ctypes

import pyhtslib.load_dll as pl

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'


# export everything from this submodule manually, including the code that
# starts with an underscore, importing modules will not import the latter
__all__ = [
    # constants (through ``#define``)
    '_HTS_IDX_NOCOOR',
    '_HTS_IDX_START',
    '_HTS_IDX_REST',
    '_HTS_IDX_NONE',
    # htslib types
    '_hts_idx_t',
    '_BGZF',
    # htslib functions
    '_hts_itr_destroy',
    '_hts_itr_next',
    '_hts_itr_query',
    '_hts_itr_queryi',
    '_hts_itr_querys',
    '_tbx_readrec',
]


_HTS_IDX_NOCOOR = -2
_HTS_IDX_START = -3
_HTS_IDX_REST = -4
_HTS_IDX_NONE = -5


class _hts_idx_t(ctypes.Structure):
    """Generic type for htslib index types"""


class _hts_itr_t(ctypes.Structure):
    """Generic type for htslib iterator types"""


class _BGZF(ctypes.Structure):
    """Type for representing a bgzip-compressed file"""


# TODO(holtgrew): return type of the following functions
htslib = pl.load_htslib()

_hts_itr_destroy = htslib.itr_destroy
_hts_itr_destroy.restype = None

_hts_itr_next = htslib.hts_itr_next
_hts_itr_next.restype = ctypes.c_int

_hts_itr_query = htslib.itr_queryi
_hts_itr_query.restype = ctypes.POINTER(_hts_itr_t)


def _hts_itr_queryi(tbx, tid, beg, end):
    return _hts_itr_query(tbx[0].idx, tid, beg, end, _tbx_readrec)


def _hts_itr_querys(tbx, s):
    return _hts_itr_query(tbx[0].idx, s, _tbx_name2id, tbx, _hts_itr_query,
                          _tbx_readrec)


_tbx_name2id = htslib.tbx_name2id
_tbx_name2id.restype = ctypes.c_int


_tbx_readrec = htslib.tbx_readrec
_tbx_readrec.restype = ctypes.c_int
