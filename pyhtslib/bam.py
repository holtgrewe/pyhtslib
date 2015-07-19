#!/usr/bin/env python3
"""Access to SAM and BAM files through htslib"""

import collections
import ctypes  # NOQA TODO(holtgrew): remove "NOQA"
import logging
import os
import os.path
import sys  # NOQA TODO(holtgrew): remove?

from pyhtslib.hts_internal import *  # NOQA
from pyhtslib.bam_internal import *  # NOQA
from pyhtslib.tabix_internal import *  # NOQA

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'


class BAMIndexException(Exception):
    """Raised when there is a problem with a BAMIndex file"""


class BAMFileException(Exception):
    """Raised when there is a problem with a BAMFile file"""


class BAMHeaderTargetInfo:
    """Information (name, length) for the reference/target sequence"""

    def __init__(self, name, length):
        #: name of the target/reference
        self.name = name
        #: length of the target/reference
        self.length = length

    def __eq__(self, other):
        return (self.name == other.name and self.length == other.length)

    def __repr__(self):
        return 'BAMHeaderTargetInfo({}, {})'.format(
            repr(self.name), self.length)


class BAMHeaderRecord:
    """Header record"""

    def __init__(self, key, tags):
        #: type of the header record
        self.key = key
        #: tags of the header record, ordered dict of ``str`` to ``str``
        #: mappings
        self.tags = collections.OrderedDict(tags)

    def for_sam_header(self):
        """Return SAM header line representation, excluding line ending"""
        def f(item):
            return '{}:{}'.format(*item)
        return '@{}\t{}'.format(self.key,
                                '\t'.join(map(f, self.tags.items())))

    def from_sam_header(line):
        arr = line.strip().split('\t')
        key = arr[0][1:]
        items = [entry.split(':', 1) for entry in arr[1:]]
        return BAMHeaderRecord(key, items)

    def __eq__(self, other):
        return self.key == other.key and self.tags == other.tags

    def __repr__(self):
        return 'BAMHeaderRecord({}, {})'.format(repr(self.key), self.tags)


class BAMHeaderComment:
    """Comments that can also appear in the header"""

    def __init__(self, text):
        #: comment text
        self.text = text

    def for_sam_header(self):
        """Return SAM header line representation, excluding line ending"""
        return '@CO\t' + self.text

    @staticmethod
    def from_sam_header(line):
        return BAMHeaderComment(line.strip().split('\t', 1)[1])

    def __eq__(self, other):
        return self.text == other.text

    def __repr__(self):
        return 'BAMHeaderComment.format({})'.format(self.text)


class BAMHeader:
    """Information stored in the BAM header"""

    @staticmethod
    def read_from_file(file_ptr):
        return BAMHeader(_sam_hdr_read(file_ptr))

    def __init__(self, struct_ptr=None):
        #: pointer to internal representation
        self.struct_ptr = struct_ptr
        #: internal representation
        self.struct = struct_ptr[0] if struct_ptr else None
        #: target/reference informations [(name, length)]
        self.target_infos = []
        #: list of ``BAMHeaderRecord and ``BAMComment``
        self.header_records = []

        self._fill_from_struct()

    def _fill_from_struct(self):
        """Fill BAM header from self.struct if given"""
        if not self.struct:
            return  # nothing given
        s = self.struct
        # read the list of (name, length) pairs
        for i in range(self.struct.n_targets):
            self.target_infos.append(
                BAMHeaderTargetInfo(s.target_name[i].decode('utf-8'),
                                    s.target_len[i]))
        # parse headers from the text
        recs = self.header_records
        for line in s.text.decode('utf-8').splitlines(True):
            if line.startswith('@CO'):
                recs.append(BAMHeaderComment.from_sam_header(line))
            else:
                recs.append(BAMHeaderRecord.from_sam_header(line))

    def free(self):
        """Free structs associated with the header

        If you use ``BAMFile`` as a context manager (using ``with``) then
        this function is called automatically when the context is left.
        """
        if not self.struct_ptr:
            return  # not set
        _bam_hdr_destroy(self.struct_ptr)
        self.struct_ptr = None
        self.struct = None

    def to_sam_header(self):
        """Return SAM header representation"""
        def f(entry):
            return entry.to_sam_header()
        return '\n'.join(map(f, self.header_record))


class BAMFile:
    """Wrapper for SAM/BAM/CRAM access

    It's strongly recommended to use as a context manager
    """

    def __init__(self, path):
        #: path to BAM file
        self.path = path
        #: wrapped C struct
        self.struct = None
        #: pointer to C struct
        self.struct_ptr = None
        #: representation of BAM header
        self.header = None

        self.open()

    def open(self):
        """Open file and read header"""
        if self.struct_ptr:
            return  # already open
        # open file and store handles
        self.struct_ptr = _hts_open(self.path.encode('utf-8'), 'r')
        self.struct = self.struct_ptr[0]
        # read header
        self.header = BAMHeader.read_from_file(self.struct_ptr)

    def close(self):
        """Close file again and free header and other data structures"""
        self.header.free()
        _hts_close(self.struct_ptr)
        self.struct_ptr = None
        self.struct = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


class BAMIndex:
    """Random-access access to BAM file"""

    @staticmethod
    def build(path, bai_path=None):
        assert False, 'Implement me!'

    def __init__(self, path, bai_path=None, require_index=True,
                 auto_load=True, auto_build=False):
        #: path to BAM file
        self.path = path
        #: path to BAI (BAM index) file
        self.bai_path = bai_path or path + '.bai'

        #: the ``BAMFile`` to use for reading
        self.file = BAMFile(self.path)
        self.file.open()

        #: wrapped C struct
        self.struct = None
        #: pointer to C struct
        self.struct_ptr = None

    def _check_file_ages(self):
        mtime_file = os.path.getmtime(self.path)
        mtime_index = os.path.getmtime(self.bai_path)
        if mtime_file > mtime_index:
            tpl = 'The file {} is older than the index file {}'
            raise BAMIndexException(tpl.format(self.path, self.bai_path))

    def _check_auto_build(self):
        if not self.auto_build:
            return
        logging.debug('auto-building index for {}'.format(self.path))
        if not os.path.exists(self.bai_path):
            BAMIndex.build(self.path, bai_path=self.bai_path)

    def _check_auto_load(self):
        if not self.auto_load:
            return
        if os.path.exists(self.bai_path):
            # check that the index is not older than the file, this is a
            # common source of errors
            self._check_file_ages()
            # load index
            self.load()
        else:
            tpl = 'BAM index required for {} for loading but not found.'
            raise BAMIndexException(tpl.format(self.path))

    def _check_require_index(self):
        if not self.require_index:
            return
        if not os.path.exists(self.bai_path):
            tpl = 'BAM index required for {} but not found.'
            raise BAMIndexException(tpl.format(self.path))

    def close(self, close_file=True):
        if self.struct_ptr:
            _hts_idx_destroy(self.struct_ptr)
        self.struct = None
        self.struct_ptr = None
        if close_file:
            self.file.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close(close_file=True)
