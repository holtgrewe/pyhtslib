#!/usr/bin/env python3
"""Implementation details for the Tabix format."""

import ctypes

import pyhtslib.load_dll as pl
import pyhtslib.hts_internal as ph

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'

# export everything from this submodule manually, including the code that
# starts with an underscore, importing modules will not import the latter
__all__ = [
    # htslib types
    '_tbx_conf_t', '_tbx_t',
    # htslib functions
    '_hts_get_bgzfp', '_tbx_destroy', '_tbx_index_build', '_tbx_index_load',
    '_tbx_index_load2', '_tbx_name2id', '_tbx_readrec', '_tbx_seqname',
]

# ----------------------------------------------------------------------------
# Structures
# ----------------------------------------------------------------------------


class _tbx_conf_t(ctypes.struct):
    """C structure for the tabix index configuration"""

    _fields_ = [('preset', ctypes.int32),
                ('sc', ctypes.int32),
                ('bc', ctypes.int32),
                ('ec', ctypes.int32),
                ('meta_char', ctypes.int32),
                ('line_skip', ctypes.int32)]


class _tbx_t(ctypes.struct):
    """C structure for a tabix index"""

    _fields_ = [('conf', _tbx_conf_t),
                ('idx', ctypes.POINTER(ph._hts_idx_t)),
                ('dict', ctypes.c_void_pt)]

# ----------------------------------------------------------------------------
# C functions and return types
# ----------------------------------------------------------------------------

htslib = pl.load_htslib()

_hts_get_bgzfp = htslib.hts_get_bgzfp
_hts_get_bgzfp.restype = ctypes.POINTER(ph._BGZF)

_tbx_destroy = htslib.tbx_destroy
_tbx_destroy.restype = None

_tbx_index_build = htslib.tbx_index_build
_tbx_index_build.restype = ctypes.c_int

_tbx_index_load = htslib.tbx_index_load
_tbx_index_load.restype = ctypes.POINTER(_tbx_t)

_tbx_index_load2 = htslib.tbx_index_load
_tbx_index_load2.restype = ctypes.POINTER(_tbx_t)

_tbx_name2id = htslib.tbx_name2id
_tbx_name2id.restype = ctypes.c_int

_tbx_readrec = htslib.hts_readrec
_tbx_readrec.restype = ctypes.c_int

_tbx_seqname = htslib.tbx_seqname
_tbx_seqname.restype = ctypes.c_char_p


def _tbx_itr_destroy(itr):
    return ph._hts_itr_destroy(itr)


def _tbx_itr_queryi(tbx, tid, beg, end):
    return ph._hts_itr_queryi(tbx, tid, beg, end)


def _tbx_itr_querys(tbx, s):
    return ph._hts_itr_querys(tbx[0].idx, s, _tbx_name2id, tbx,
                              ph._hts_itr_query, _tbx_readrec)


def _tbx_itr_next(htsfp, tbx, itr, r):
    return ph._hts_itr_next(_hts_get_bgzfp(htsfp), itr, r, tbx)


def _tbx_bgzf_itr_next(bgzfp, tbx, itr, r):
    return ph._hts_itr_next(bgzfp, itr, r, tbx)
