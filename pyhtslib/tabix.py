#!/usr/bin/env python3
"""Wrapper for accessing tabix-indexed files."""

import os

from tabix_internal import *  # NOQA

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'


class TabixIndexException(Exception):
    """Raised when there is a problem with a TabixIndex file."""


class TabixConfig:
    """Configuration for tabix-indxed files"""

    @staticmethod
    def from_c_struct(self, name):
        res = _tbx_conf_t.in_dll(htslib, name)
        return TabixConfig(res.preset, res.seq_col, res.begin_col,
                           res.end_col, res.meta_char, res.line_skip)

    def __init__(self, preset, seq_col, begin_col, end_col, meta_char,
                 line_skip):
        self.preset = preset
        self.seq_col = seq_col
        self.begin_col = begin_col
        self.end_col = end_col
        self.meta_char = meta_char

    def to_c_struct(self):
        """Convert to C struct"""
        res = _tbx_conf_t()
        res.preset = self.preset
        res.sc = self.sc
        res.bc = self.bc
        res.ec = self.ec
        res.meta_char = self.meta_char


TBX_CONF_GFF = TabixConfig.from_c_struct('tbx_conf_gff')
TBX_CONF_BED = TabixConfig.from_c_struct('tbx_conf_bed')
TBX_CONF_PLSTBL = TabixConfig.from_c_struct('tbx_conf_plstbl')
TBX_CONF_SAM = TabixConfig.from_c_struct('tbx_conf_sam')
TBX_CONF_VCF = TabixConfig.from_c_struct('tbx_conf_vcf')


class TabixFileIter:
    """Allows iteration over tabix files"""

    def __init__(self, index, struct_ptr):
        self.index = index
        self.struct_ptr = struct_ptr
        self.struct = self.struct_pr[0]

    def __iter__(self):
        return self

    def __next__(self):
        if _tbx_itr_next(htsfp, self.index.struct_ptr, r, data):
            pass
        raise StopIteration()


class TabixIndex:
    """Tabix index"""

    @staticmethod
    def build(path, min_shift, config, tbi_path=None):
        """Build tabix index for file at path

        If ``tbi_path`` is given then the index is built at the given path.
        Otherwise, it is built at ``${path}.tbi``.
        """
        if tbi_path:
            _tbx_index_build2(path, tbi_path, min_shift,
                              ctypes.byref(config.to_c_struct_ptr))
        else:
            _tbx_index_build2(path, min_shift,
                              ctypes.byref(config.to_c_struct_ptr))

    def __init__(self, path, tbi_path=None, require_index=False,
                 auto_load=True, auto_build=True):
        #: path to indexed file
        self.path = path
        #: path to index file
        self.tbi_path = tbi_path or path + '.tbi'
        #: whether or not to require existence of FAI index file
        self.require_index = require_index
        #: whether or not to load index if it exists
        self.auto_load = auto_load
        #: whether or not to automatically build index if it does not
        #: exist yet, overrides ``require_index``
        self.auto_build = auto_build

        #: wrapped C struct
        self.struct = None
        #: pointer to C struct
        self.struct_ptr = None

        self._check_auto_build()
        self._check_auto_load()
        self._check_require_index()

    def _check_auto_build(self):
        if not self.auto_build:
            return
        if not os.path.exists(self.tai_path):
            TabixIndex.build(self.path)

    def _check_auto_load(self):
        if not self.auto_load:
            return
        if os.path.exists(self.tai_path):
            self.load()
        else:
            tpl = 'Tabix index required for {} for loading but not found.'
            raise TabixIndexException(tpl.format(self.path))

    def _check_require_index(self):
        if not self.require_index:
            return
        if not os.path.exists(self.tai_path):
            tpl = 'Tabix index required for {} but not found.'
            raise TabixIndexException(tpl.format(self.path))

    def query(self, region_str):
        return TabixFileIter(self, _tbx_itr_querys(self.struct_ptr,
                                                   region_str))

    def query_region(self, seq, begin, end):
        return TabixFileIter(self, _tbx_itr_queryi(
            self.struct_ptr, seq, begin, end))

    def from_start(self):
        return TabixFilteIter(self, _tbx_itr_query(self.struct_ptr))

    def from_current(self):
        return TabixFilteIter(self, _hts_itr_query(
            self.struct_ptr, _HTS_IDX_REST, 0, 0, _tbx_readrec))

    def load(self):
        self.close()
        if self.tai_path:
            self.struct_ptr = _tbx_index_load2(self.path, self.tbi_path)
            if not self.struct_ptr:
                tpl = 'Could not load tabix index for {}'
                raise TabixIndexException(tpl.format(self.path))
        else:
            self.struct_ptr = _tbx_index_lod(self.path)
            if not self.struct_ptr:
                tpl = 'Could not load tabix index {} for {}'
                raise TabixIndexException(tpl.format(self.tbi_path,
                                                     self.path))
        self.struct = self.struct_ptr[0]

    def close(self):
        if self.struct_ptr:
            _tbx_destroy(self.struct_ptr)
        self.struct = None
        self.struct_ptr = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
