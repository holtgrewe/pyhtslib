#!/usr/bin/env python3
"""Random access to FASTA files."""

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'

import collections
import ctypes
import logging
import os.path
import sys

from faidx_internal import *


# TODO(holtgrewe): test this stuff!


class FASTAIndexException(Exception):
    """Raised when there is a problem with a FASTAIndex file."""


class GenomeInterval:
    """Zero-based genome interval."""

    def __init__(self, seq, begin_pos, end_pos):
        self.seq = seq
        self.begin_pos = begin_pos
        self.end_pos = end_pos

    def __str__(self):
        return '{}:{:}-{:}'.format(self.seq, self.begin_pos + 1, self.end_pos)


class FASTAIndexSequenceRecord:
    """FAI entry."""

    def __init__(self, name, length):
        #: name of the sequence
        self.name = name
        #: length of the seqeuence
        self.length = length

    def __repr__(self):
        tpl = '<FASTAIndexSequenceRecord({}, {})>'
        return tpl.format(*map(repr, [self.name, self.length]))


class FASTAIndex:
    """Index for a FASTA file, allows random access to FASTA file."""

    def __init__(self, fasta_path, require_index=False, auto_load=True,
                 auto_build=True):
        #: path to FASTA file with the index, index has path
        #: "${fasta_path}.fai"
        self.fasta_path = fasta_path
        #: whether or not to require existence of FAI index file
        self.require_index = require_index
        #: whether or not to load index if it exists
        self.auto_load = auto_load
        #: whether or not to automatically build index if it does not
        #: exist yet, overrides ``require_index``
        self.auto_build = True
        #: the pointer to the ``FAIDXStruct``
        self.struct_ptr = None
        #: the sequence dictionary
        self.seq_dict = None

        if self.auto_build:
            self._check_auto_build()
        if self.auto_load:
            self._check_auto_load()
        if self.require_index:
            self._check_require_index()

    def _check_auto_build(self):
        if not os.path.exists(self.fai_path):
            FASTAIndex.build(self.fasta_path)

    def _check_auto_load(self):
        if os.path.exists(self.fai_path):
            self.load()

    def _check_require_index(self):
        if not os.path.exists(self.fai_path):
            tpl = 'FASTA index required for {} but not found.'
            raise FastaIndexException(tpl.format(self.fasta_path))

    def load(self):
        """Load FAI index."""
        self.struct_ptr = _fai_load(self.fasta_path.encode('utf-8'))
        if not self.struct_ptr:
            tpl = 'Failed to load FASTA index for FASTA file {}'
            raise FastaIndexException(tpl.format(self.fasta_path))
        self.seq_dict = self._build_seq_dict()

    def _build_seq_dict(self):
        """Build sequence dictionary."""
        result = collections.OrderedDict()
        num_seqs = _faidx_nseq(self.struct_ptr)
        for i in range(num_seqs):
            seq = _faidx_iseq(self.struct_ptr, i)
            length = _faidx_seq_len(self.struct_ptr, seq)
            seq = seq.decode('utf-8')
            result[seq] = FASTAIndexSequenceRecord(seq, length)
        return result

    def fetch(self, region):
        if type(region) is GenomeInterval:
            region = str(region)
        region = region.encode('utf-8')
        res_len = ctypes.c_int()
        res = _fai_fetch(self.struct_ptr, region, ctypes.byref(res_len))
        return res.decode('utf-8')

    def close(self):
        """Free memory for self.struct_ptr if any."""
        if not self.struct_ptr:
            return
        logging.debug('Freeing FAI for %s', self.fasta_path)
        _fai_destroy(self.struct_ptr)
        self.struct_ptr = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    @property
    def fai_path(self):
        """Path to FAI file for ``self.fasta_path``."""
        return '{}.fai'.format(self.fasta_path)

    @staticmethod
    def build(fasta_path):
        """Build index for the file at the given path."""
        logging.debug('Building FAI file for %s', fasta_path)
        _fai_build(fasta_path.encode('utf-8'))

if __name__ == '__main__':
    logging.debug('starting...')
    logging.basicConfig(level=logging.DEBUG)
    with FASTAIndex('test.fa') as faidx:
        print(faidx.seq_dict, file=sys.stderr)
        print('in with statement', file=sys.stderr)
        print(faidx.fetch('HSBGPG:2-30'), file=sys.stderr)
    logging.debug('done!')
