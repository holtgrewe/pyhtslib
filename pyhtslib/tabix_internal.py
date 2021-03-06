#!/usr/bin/env python3
"""Implementation details for the Tabix format"""

import ctypes

import pyhtslib.load_dll as pl
from pyhtslib.hts_internal import *  # NOQA

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'

# export everything from this submodule manually, including the code that
# starts with an underscore, importing modules will not import the latter
__all__ = [
    # htslib types
    '_tbx_conf_t',
    '_tbx_t',
    # htslib functions
    '_hts_get_bgzfp',
    '_tbx_destroy',
    '_tbx_index_build',
    '_tbx_index_build2',
    '_tbx_index_load',
    '_tbx_index_load2',
    '_tbx_itr_destroy',
    '_tbx_itr_next',
    '_tbx_itr_queryi',
    '_tbx_itr_querys',
    '_tbx_name2id',
    '_tbx_readrec',
    '_tbx_seqnames',
]

# ----------------------------------------------------------------------------
# Structures
# ----------------------------------------------------------------------------


class _tbx_conf_t(ctypes.Structure):
    """C structure for the tabix index configuration"""

    _fields_ = [('preset', ctypes.c_int32),
                ('sc', ctypes.c_int32),
                ('bc', ctypes.c_int32),
                ('ec', ctypes.c_int32),
                ('meta_char', ctypes.c_int32),
                ('line_skip', ctypes.c_int32)]


class _tbx_t(ctypes.Structure):
    """C structure for a tabix index"""

    _fields_ = [('conf', _tbx_conf_t),
                ('idx', ctypes.POINTER(_hts_idx_t)),
                ('dict', ctypes.c_void_p)]

# ----------------------------------------------------------------------------
# C functions and their return types
# ----------------------------------------------------------------------------

htslib = pl.load_htslib()

_hts_get_bgzfp = htslib.hts_get_bgzfp
_hts_get_bgzfp.restype = ctypes.POINTER(_BGZF)

_tbx_destroy = htslib.tbx_destroy
_tbx_destroy.restype = None

_tbx_index_build = htslib.tbx_index_build
_tbx_index_build.restype = ctypes.c_int

_tbx_index_build2 = htslib.tbx_index_build2
_tbx_index_build2.restype = ctypes.c_int

_tbx_index_load = htslib.tbx_index_load
_tbx_index_load.restype = ctypes.POINTER(_tbx_t)

_tbx_index_load2 = htslib.tbx_index_load
_tbx_index_load2.restype = ctypes.POINTER(_tbx_t)

_tbx_name2id = htslib.tbx_name2id
_tbx_name2id.restype = ctypes.c_int

_tbx_readrec = htslib.tbx_readrec
_tbx_readrec.restype = ctypes.c_int

_tbx_seqnames = htslib.tbx_seqnames
_tbx_seqnames.restype = ctypes.POINTER(ctypes.c_char_p)


def _tbx_itr_destroy(itr):
    return _hts_itr_destroy(itr)


def _tbx_itr_next(htsfp, tbx, itr, r):
    return _hts_itr_next(_hts_get_bgzfp(htsfp), itr, r, tbx)


def _tbx_itr_queryi(tbx, tid, beg, end):
    return _hts_itr_queryi(tbx, tid, beg, end)


def _tbx_itr_querys(tbx, s):
    return _hts_itr_querys(tbx[0].idx, s, _tbx_name2id, tbx,
                           _hts_itr_query, _tbx_readrec)


def _tbx_bgzf_itr_next(bgzfp, tbx, itr, r):
    return _hts_itr_next(bgzfp, itr, r, tbx)
