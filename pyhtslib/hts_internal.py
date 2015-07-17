#!/usr/bin/env python3
"""Shared htslib structs and functions"""

import ctypes

import pyhtslib.load_dll as pl

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'


__all__ = [
    # exported types
    '_hts_idx_t', '_BGZF',
    # exported functions
    '_hts_itr_destroy', '_hts_itr_queryi', '_hts_itr_querys',
    '_hts_itr_next',
]


class _hts_idx_t(ctypes.Structure):
    """Generic type for htslib index types"""


class _BGZF(ctypes.Structure):
    """Type for representing a bgzip-compressed file"""


# TODO(holtgrew): return type of the following functions
htslib = pl.load_htslib()

_hts_itr_destroy = htslib.itr_destroy
_hts_itr_destroy.restype = None

_hts_itr_queryi = htslib.itr_queryi
_hts_itr_queryi.restype = None

_hts_itr_querys = htslib.itr_querys
_hts_itr_querys.restype = None

_hts_itr_next = htslib.hts_itr_next
_hts_itr_next.restype = None
