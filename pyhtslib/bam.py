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


class CIGARElement:
    """CIGAR element"""

    def __init__(self, count, operation):
        #: CIGAR count
        self.count = count
        #: CIGAR operation
        self.operation = operation

    def __repr__(self):
        return 'CIGARElement({}, {})'.format(
            self.count, repr(self.operation))

    def __str__(self):
        return '{}{}'.format(self.count, self.operation)


class BAMAuxTagParser:
    """Helper for parsing aux field tags from BAM records"""

    def __init__(self, aux, aux_len):
        self.aux = aux
        self.aux_len = aux_len
        self.keys = self.get_keys()
        print('keys={}'.format(self.keys))

    def get_keys(self):
        """Get keys"""
        def aux_type2size(t):
            DIR = {'A': 1, 'c': 1, 'C': 1, 's': 2, 'S': 2, 'i': 4, 'I': 4,
                   'f': 4, 'd': 8, 'Z': 'Z', 'H': 'H', 'B': 'B'}
            return DIR.get(t, 0)

        def skip_aux(s):
            # skip the type
            size = aux_type2size(s[0])
            s += 1
            if size in ['Z', 'H']:
                while s[0]:
                    s += 1
                return s + 1
            elif size == 'B':
                size = aux_type2size(s[0])
                s += 1
                n = ctypes.c_uint32(0)
                ctypes.memmove(ctypes.byref(n), s, 4)
                return s + size * n.value
            elif size == 0:
                raise BAMFileException('Problem parsing auxiliary fields!')
            else:
                return s + size

        result = []
        s = self.aux
        while s < aux + aux_len:
            result.append('{}{}'.format(chr(s[0]), char(s[1])))
            s = skip_aux(s)
        return result

    def __call__(self):
        result = collections.OrderedDict()
        for key in self.keys:
            # TODO(holtgrewe): fix/finish this
            result[key] = None
        return result


class BAMRecordImpl:
    """Information extracted from C internals of ``BAMRecord``"""

    @staticmethod
    def from_struct(ptr, header):
        """Return ``BAMRecordImpl`` from internal C structure"""
        return BAMRecordImpl(
            str(_bam_get_qname(ptr)),
            ptr[0].core.flag,
            ptr[0].core.tid,
            header.target_infos[ptr[0].core.tid].name,
            ptr[0].core.pos,
            _bam_endpos(ptr),
            ptr[0].core.qual,
            BAMRecordImpl._cigar(ptr),
            ptr[0].core.mtid,
            header.target_infos[ptr[0].core.mtid].name,
            ptr[0].core.mpos,
            ptr[0].core.isize,
            str(_bam_get_seq(ptr)),
            str(_bam_get_qual(ptr)),
            BAMAuxTagParser(_bam_get_aux(ptr), _bam_get_l_aux(ptr))(),
            )

    @staticmethod
    def _cigar(ptr):
        def to_ce(cigar):
            return CIGARElement(_bam_cigar_oplen(cigar),
                                _bam_cigar_opchr(cigar))

        cigar_arr = _bam_get_cigar(ptr)
        return [to_ce(cigar_arr[i]) for i in range(ptr[0].core.l_qseq)]

    def __init__(self, qname, flag, r_id, ref, begin_pos, end_pos, mapq,
                 cigar, r_id_next, ref_next, pos_next, tlen, seq, qual,
                 tags):
        #: read name (QNAME)
        self.qname = qname
        #: numeric flag (FLAG)
        self.flag = flag
        #: reference id of fragment's alignment (RID)
        self.r_id = r_id
        #: reference name of the alignment
        self.ref = ref
        #: begin position of the alignment (POS)
        self.begin_pos = begin_pos
        #: 0-base end position of alignment (inferred from ``begin_pos``
        #: and ``cigar``
        self.end_pos = end_pos
        #: mapping quality (MAPQ)
        self.mapq = mapq
        #: CIGAR string
        self.cigar = cigar
        #: reference id of next fragment's alignment
        self.r_id_next = r_id_next
        #: reference name of next fragment's alignment
        self.ref_next = ref_next
        #: 0-based position of next fragment's alignment
        self.pos_next = pos_next
        #: template length (TLEN)
        self.tlen = tlen
        #: sequence string (SEQ)
        self.seq = seq
        #: quality string (QUAL)
        self.qual = qual
        #: tags, as ``OrderedDict``
        self.tags = tags

    @property
    def is_paired(self):
        return (self.flag & _BAM_FPAIRED) != 0

    @property
    def is_proper_pair(self):
        return (self.flag & _BAM_FPROPER_PAIR) != 0

    @property
    def is_unmapped(self):
        return (self.flag & _BAM_UNMAP) != 0

    @property
    def is_mate_unmapped(self):
        return (self.flag & _BAM_MUNMAP) != 0

    @property
    def is_reversed(self):
        return (self.flag & _BAM_FREVERSE) != 0

    @property
    def is_mate_reversed(self):
        return (self.flag & _BAM_FMREVERSE) != 0

    @property
    def is_first(self):
        return (self.flag & _BAM_FREAD1) != 0

    @property
    def is_last(self):
        return (self.flag & _BAM_FREAD2) != 0

    @property
    def is_secondary(self):
        return (self.flag & _BAM_FSECONDARY) != 0

    @property
    def is_qc_fail(self):
        return (self.flag & _BAM_FQC_FAIL) != 0

    @property
    def is_duplicate(self):
        return (self.flag & _BAM_FDUPLICATE) != 0

    @property
    def is_supplementary(self):
        return (self.flag & _BAM_FSUPPLEMENTARY) != 0


class BAMRecord:
    """Record from a BAM file"""

    def __init__(self, struct_ptr, header):
        #: pointer to wrapped C struct
        self.struct_ptr = struct_ptr
        #: wrapped C struct
        self.struct = self.struct_ptr[0]
        #: ``BAMHeader`` for references
        self.header = header
        #: ``BAMRecordImpl`` instance used for the representation
        self.impl = None

    def detach(self):
        """Detach from wrapped C struct

        The only way of obtaining ``BAMRecords`` is through iterating
        ``BAMFile`` objects, either through an index or not.  For efficiency,
        the reading reuses the same buffer for reading.  If you want to keep a
        ``BAMRecord`` around for longer than the current iteration then you
        have to obtain a copy that is independent of the current buffer
        through the use of ``detach()``.
        """
        if not self.impl:
            self.impl = BAMRecordImpl.from_struct(self.struct_ptr)
        self.struct_ptr = None
        self.struct = None

    def _reset(self):
        """Reset Python side, as if freshly constructed"""
        self.impl = None

    def __getattr__(self, name):
        """Delegation to self.impl if set, auto-set from struct"""
        if not self.impl and self.struct:
            self.impl = BAMRecordImpl.from_struct(self.struct)
        elif not self.impl and not self.struct:
            raise AttributeError('self.impl is None and cannot rebuild from '
                                 'None self.struct')
        else:
            return getattr(self.impl, name)


class BAMFileIter:
    """Iterate over a ``BAMFile``

    Do not use directly but by iterating over ``BAMFile``.  Iteration must
    be completed or ``close()`` must be called to prevent resource leaks.
    """

    def __init__(self, bam_file):
        #: the ``BAMFile`` to iterate through
        self.bam_file = bam_file
        #: pointer to buffer for reading in the file record by record
        self.struct_ptr = _bam_init1()
        #: buffer for readin in the file itself
        self.struct = self.struct_ptr[0]
        #: ``BAMRecord`` meant for consumption by the user
        self.record = BAMRecord(self.struct_ptr, self.bam_file.header)

    def __next__(self):
        r = _sam_read1(self.bam_file.struct_ptr,
                       self.bam_file.header.struct_ptr,
                       self.struct_ptr)
        if r >= 0:
            # successfully read record from file
            self.record._reset()
            return self.record
        else:
            # end of file or something went wrong
            self.close()
            if r < -1:
                tpl = 'truncated file {}'
                raise BAMFileException(tpl.format(self.bam_file.path))
            else:
                self.close()
                raise StopIteration

    def close(self):
        if not self.struct_ptr:
            return
        _bam_destroy1(self.struct_ptr)
        self.struct_ptr = None
        self.struct = None


class BAMFile:
    """Wrapper for SAM/BAM/CRAM access

    It's strongly recommended to use as a context manager or through
    ``BAMIndex``.
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

    def __iter__(self):
        return BAMFileIter(self)

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
