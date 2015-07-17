#!/usr/bin/env python3
"""Wrapper for accessing tabix-indexed files."""

import os

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'


class TabixIndexException(Exception):
    """Raised when there is a problem with a TabixIndex file."""


class TabixConfig:
    """Configuration for tabix-indxed files"""

    def __init__(self, preset, seq_col, begin_col, end_col, meta_char,
                 line_skip):
        self.preset = preset
        self.seq_col = seq_col
        self.begin_col = begin_col
        self.end_col = end_col
        self.meta_char = meta_char


class TabixGFFConfig(TabixConfig):
    """Tabix presets for GFF"""

    def __init__(self):
        pass


class TabixBEDConfig(TabixConfig):
    """Tabix presets for BED"""

    def __init__(self):
        pass


class TabixPSLTBLConfig(TabixConfig):
    """Tabix presets for PSLTBL

    This is the output of Jim Kent's BLAT, for example.
    """

    def __init__(self):
        pass


class TabixSAMConfig(TabixConfig):
    """Tabix presets for SAM"""

    def __init__(self):
        pass


class TabixVCFConfig(TabixConfig):
    """Tabix presets for VCF"""

    def __init__(self):
        pass


class TabixIndex:
    """Tabix index"""

    @staticmethod
    def build(path, min_shift, config, tbi_path=None):
        pass

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

        self._check_auto_build()
        self._check_auto_load()
        self._check_require_index()

    def _check_auto_build(self):
        if not self.auto_build:
            return
        if not os.path.exists(self.fai_path):
            TabixIndex.build(self.fasta_path)

    def _check_auto_load(self):
        if not self.auto_load:
            return
        if os.path.exists(self.fai_path):
            self.load()
        else:
            tpl = 'Tabix index required for {} for loading but not found.'
            raise TabixIndexException(tpl.format(self.fasta_path))

    def _check_require_index(self):
        if not self.require_index:
            return
        if not os.path.exists(self.fai_path):
            tpl = 'Tabix index required for {} but not found.'
            raise TabixIndexException(tpl.format(self.fasta_path))

    def query(self, region_str):
        pass

    def query_region(self, seq, begin, end):
        pass

    def load(self):
        pass

    def close(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
