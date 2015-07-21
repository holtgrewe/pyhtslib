#!/usr/bin/env python3
"""Wrapper for accessing tabix-indexed files"""

import ctypes
import logging
import os
import os.path

from pyhtslib.hts_internal import *  # NOQA
from pyhtslib.tabix_internal import *  # NOQA

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'


class TabixIndexException(Exception):
    """Raised when there is a problem with a TabixIndex file"""


class TabixFileException(Exception):
    """Raised when there is a problem with a TabixFile file"""


class TabixConfig:
    """Configuration for tabix-indxed files"""

    @staticmethod
    def create_empty():
        """Create TabixConfig with zero entries only.

        This is used for BCF and BAM format.
        """
        return TabixConfig(0, 0, 0, 0, 0, 0, 0)

    @staticmethod
    def from_c_struct(name, min_shift=None):
        res = _tbx_conf_t.in_dll(htslib, name)
        return TabixConfig(res.preset, res.sc, res.bc, res.ec,
                           res.meta_char, res.line_skip, min_shift)

    @staticmethod
    def from_extension(path, min_shift=None):
        """Return TabixConfig for the given extension"""
        if path.endswith('.gff.gz'):
            return TabixConfig.from_c_struct('tbx_conf_gff', min_shift)
        elif path.endswith('.bed.gz'):
            return TabixConfig.from_c_struct('tbx_conf_bed', min_shift)
        elif path.endswith('.sam.gz'):
            return TabixConfig.from_c_struct('tbx_conf_sam', min_shift)
        elif path.endswith('.vcf.gz'):
            return TabixConfig.from_c_struct('tbx_conf_vcf', min_shift)
        elif path.endswith('.bcf') or path.endswith('.bam'):
            # BCF and BAM are auto-detected, no need for preset
            return TabixConfig.create_empty()
        else:
            tpl = 'Unknown file type for {}'
            raise TabixIndexException(tpl.format(path))

    def __init__(self, preset, seq_col, begin_col, end_col, meta_char,
                 line_skip, min_shift=None):
        self.preset = preset
        self.seq_col = seq_col
        self.begin_col = begin_col
        self.end_col = end_col
        self.meta_char = meta_char
        self.min_shift = min_shift

    def to_c_struct(self):
        """Convert to C struct"""
        res = _tbx_conf_t()
        res.preset = self.preset
        res.sc = self.seq_col
        res.bc = self.begin_col
        res.ec = self.end_col
        res.meta_char = self.meta_char
        return res


TBX_CONF_GFF = TabixConfig.from_c_struct('tbx_conf_gff')
TBX_CONF_BED = TabixConfig.from_c_struct('tbx_conf_bed')
TBX_CONF_PSLTBL = TabixConfig.from_c_struct('tbx_conf_psltbl')
TBX_CONF_SAM = TabixConfig.from_c_struct('tbx_conf_sam')
TBX_CONF_VCF = TabixConfig.from_c_struct('tbx_conf_vcf')


class NormalTabixFileIter:
    """Allows iteration over tabix files after querying"""

    def __init__(self, index, struct_ptr):
        self.index = index
        self.struct_ptr = struct_ptr
        self.struct = self.struct_ptr[0]
        self._buffer = _kstring_t(0, 0, None)

    def __iter__(self):
        return self

    def __next__(self):
        if _tbx_itr_next(self.index.file.struct_ptr, self.index.struct_ptr,
                         self.struct_ptr, ctypes.byref(self._buffer)) >= 0:
            return self._buffer.s.decode('utf-8')
        else:
            self.close()
            if self._buffer:
                self._buffer.free_p()
            raise StopIteration()

    def close(self):
        """Free all associated resources

        This function is idempotent.
        """
        if self.struct_ptr:
            _hts_itr_destroy(self.struct_ptr)
            self.struct_ptr = None
            self.struct = None
        if self._buffer:
            self._buffer.free_p()
            self._buffer = None


class AllTabixFileIter:
    """Allows iteration over the whole tabix file"""

    def __init__(self, index):
        self.index = index
        self._buffer = _kstring_t(0, 0, None)
        self.current_chrom = iter(self._fetch_chroms())
        seq = next(self.current_chrom)
        self.struct_ptr = _tbx_itr_querys(self.index.struct_ptr, seq)
        self.struct = self.struct_ptr[0]

    def __iter__(self):
        return self

    def _fetch_chroms(self):
        nseq = ctypes.c_int()
        seq = _tbx_seqnames(self.index.struct_ptr, ctypes.byref(nseq))
        result = [seq[i] for i in range(nseq.value)]
        _libc.free(seq)
        return result

    def __next__(self):
        if not self.current_chrom:
            raise StopIteration()
        if _tbx_itr_next(self.index.file.struct_ptr, self.index.struct_ptr,
                         self.struct_ptr, ctypes.byref(self._buffer)) >= 0:
            return self._buffer.s.decode('utf-8')
        else:
            while True:
                try:
                    seq = next(self.current_chrom)
                    _tbx_itr_destroy(self.struct_ptr)
                    self.struct_ptr = _tbx_itr_querys(
                        self.index.struct_ptr, seq)
                    self.struct = self.struct_ptr[0]
                    if _tbx_itr_next(self.index.file.struct_ptr,
                                     self.index.struct_ptr,
                                     self.struct_ptr,
                                     ctypes.byref(self._buffer)) >= 0:
                        return self._buffer.s.decode('utf-8')
                except StopIteration:
                    self.close()
                    raise StopIteration()

    def close(self):
        """Free resources associated with the iterator

        This function is idempotent.
        """
        self.current_chrom = None
        if self.struct_ptr:
            _tbx_itr_destroy(self.struct_ptr)
            self.struct_ptr = None
            self.struct = None
        if self._buffer:
            self._buffer.free_p()
            self._buffer = None


class TabixFile:
    """Tabix file"""

    def __init__(self, path):
        #: path to the indexed file
        self.path = path
        #: wrapped C struct
        self.struct = None
        #: pointer to C struct
        self.struct_ptr = None

    def open(self):
        # check that the file was opened using bgzip
        if _bgzf_is_bgzf(self.path.encode('utf-8')) != 1:
            tpl = '{} was not compressed using bgzip'
            raise TabixFileException(tpl.format(self.path))
        self.struct_ptr = _hts_open(self.path.encode('utf-8'), 'r')
        if not self.struct_ptr:
            tpl = 'Opening tabix file {} failed'
            raise TabixFileException(tpl.format(self.path))
        self.struct = self.struct_ptr[0]

    def close(self):
        """Free all associated resources

        This function is idempotent.
        """
        if self.struct_ptr:
            if _hts_close(self.struct_ptr) != 0:
                tpl = 'Problem closing tabix file {}'
                raise TabixFileException(tpl.format(self.path))
            self.struct_ptr = None
            self.struct = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


class TabixIndex:
    """Tabix index"""

    @staticmethod
    def build(path, config=None, tbi_path=None):
        """Build tabix index for file at path

        If ``tbi_path`` is given then the index is built at the given path.
        Otherwise, it is built at ``${path}.tbi``.

        If ``min_shift`` and ``config`` are both not given (``None``) then
        we attempt to recognize it from the file format.
        """
        # auto-build configuration if necessary
        config = config or TabixConfig.from_extension(path)
        if not config:
            raise TabixIndexException('Cannot build tabix index as we could '
                                      'not recognize the file type')
        # build tabix index
        if tbi_path:
            _tbx_index_build2(path.encode('utf-8'), tbi_path.encode('utf-8'),
                              config.min_shift,
                              ctypes.byref(config.to_c_struct()))
        else:
            _tbx_index_build(path.encode('utf-8'), config.min_shift,
                             ctypes.byref(config.to_c_struct()))

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

        #: the ``TabixFile`` to use for reading
        self.file = TabixFile(self.path)
        self.file.open()

        #: wrapped C struct
        self.struct = None
        #: pointer to C struct
        self.struct_ptr = None

        # collection of iterators, we will call close() on all of them
        # in our own close to ensure that all memory is freed
        self.iterators = []

        self._check_auto_build()
        self._check_auto_load()
        self._check_require_index()

    def _check_file_ages(self):
        mtime_file = os.path.getmtime(self.path)
        mtime_index = os.path.getmtime(self.tbi_path)
        if mtime_file > mtime_index:
            tpl = 'The file {} is older than the index file {}'
            raise TabixIndexException(tpl.format(self.path, self.tbi_path))

    def _check_auto_build(self):
        if not self.auto_build:
            return
        logging.debug('auto-building index for {}'.format(self.path))
        if not os.path.exists(self.tbi_path):
            TabixIndex.build(self.path, tbi_path=self.tbi_path)

    def _check_auto_load(self):
        if not self.auto_load:
            return
        if os.path.exists(self.tbi_path):
            # check that the index is not older than the file, this is a
            # common source of errors
            self._check_file_ages()
            # load index
            self.load()
        else:
            tpl = 'Tabix index required for {} for loading but not found.'
            raise TabixIndexException(tpl.format(self.path))

    def _check_require_index(self):
        if not self.require_index:
            return
        if not os.path.exists(self.tbi_path):
            tpl = 'Tabix index required for {} but not found.'
            raise TabixIndexException(tpl.format(self.path))

    def get_header(self):
        """Return string with header lines resets iteration to file"""
        result = []
        buf = _kstring_t(0, 0, None)
        while _hts_getline(self.file.struct_ptr, _KS_SEP_LINE,
                           ctypes.byref(buf)) >= 0:
            if not buf.l or buf.s[0] != self.struct.conf.meta_char:
                break
            else:
                result.append(buf.s.decode('utf-8'))
                result.append('\n')
        buf.free_p()
        return ''.join(result)

    # TODO(holtgrewe): fix exception display if not region_string
    def query(self, region_str=None, seq=None, begin=None, end=None):
        if (region_str is None and
                (seq is None or begin is None or end is None)):
            raise TabixIndexException(
                'You have to either give region_str or seq/begin/end')
        if region_str:
            ptr = _tbx_itr_querys(self.struct_ptr,
                                  region_str.encode('utf-8'))
        else:
            ptr = _tbx_itr_queryi(self.struct_ptr, seq, begin, end)
        if not ptr:
            tpl = 'Could not jump to {}'
            raise TabixIndexException(tpl.format(region_str))
        self.iterators.append(NormalTabixFileIter(self, ptr))
        return self.iterators[-1]

    def from_start(self):
        self.iterators.append(AllTabixFileIter(self))
        return self.iterators[-1]

    def __iter__(self):
        return iter(self.from_start())

    def load(self):
        self.close(close_file=False)
        if self.tbi_path:
            self.struct_ptr = _tbx_index_load2(self.path.encode('utf-8'),
                                               self.tbi_path.encode('utf-8'))
            if not self.struct_ptr:
                tpl = 'Could not load tabix index for {}'
                raise TabixIndexException(tpl.format(self.path))
        else:
            self.struct_ptr = _tbx_index_load(self.path.encode('utf-8'))
            if not self.struct_ptr:
                tpl = 'Could not load tabix index {} for {}'
                raise TabixIndexException(tpl.format(self.tbi_path,
                                                     self.path))
        self.struct = self.struct_ptr[0]

    def close(self, close_file=True):
        if self.struct_ptr:
            _tbx_destroy(self.struct_ptr)
            self.struct = None
            self.struct_ptr = None
        if close_file:
            self.file.close()
        for it in self.iterators:
            it.close()
        self.iterators = []

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close(close_file=True)
