#!/usr/bin/env python3
"""C wrapper code."""

import ctypes

import pyhtslib.load_dll as pl

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'

htslib = pl.load_htslib()

# export everything from this submodule manually, including the code that
# starts with an underscore, importing modules will not import the latter
__all__ = [
    'FAIDXStruct',
    '_fai_build',
    '_fai_destroy',
    '_fai_load',
    '_fai_fetch',
    '_faidx_nseq',
    '_faidx_iseq',
    '_faidx_seq_len',
]


class FAIDXStruct(ctypes.Structure):
    """Opaque structure from htslib."""

_fai_build = htslib.fai_build
_fai_build.restype = ctypes.c_int

_fai_destroy = htslib.fai_destroy
_fai_destroy.restype = None

_fai_load = htslib.fai_load
_fai_load.restype = ctypes.POINTER(FAIDXStruct)

_fai_fetch = htslib.fai_fetch
_fai_fetch.restype = ctypes.c_char_p

_faidx_nseq = htslib.faidx_nseq
_faidx_nseq.restype = ctypes.c_int

_faidx_iseq = htslib.faidx_iseq
_faidx_iseq.restype = ctypes.c_char_p

_faidx_seq_len = htslib.faidx_seq_len
_faidx_seq_len.restype = ctypes.c_int
