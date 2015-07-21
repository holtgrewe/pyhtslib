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
    # TODO(holtgrewe): This probably works better using ctypes.string_at and
    # then performing all parsing in native Python.

    def __init__(self, struct_ptr, aux, aux_len):
        self.struct_ptr = struct_ptr
        self.aux = aux
        self.aux_len = aux_len
        self.keys = self.get_keys()

    def get_keys(self):
        """Get keys"""
        def aux_type2size(t):
            DIR = {'A': 1, 'c': 1, 'C': 1, 's': 2, 'S': 2, 'i': 4, 'I': 4,
                   'f': 4, 'd': 8, 'Z': 'Z', 'H': 'H', 'B': 'B'}
            return DIR.get(t, 0)

        def skip_aux(s):
            # skip the type
            def deref(s):
                return ctypes.c_char.from_address(s).value
            size = aux_type2size(deref(s).decode('utf-8'))
            # print('c={}, size={}'.format(deref(s), size), file=sys.stderr)
            s += 1
            if size in ['Z', 'H']:
                while deref(s) != b'\x00':
                    s += 1
                return s + 1
            elif size == 'B':
                size = aux_type2size(deref(s).decode('utf-8'))
                s += 1
                n = ctypes.c_uint32(0)
                ctypes.memmove(ctypes.byref(n), s, 4)
                return s + size * n.value
            elif size == 0:
                raise BAMFileException('Problem parsing auxiliary fields!')
            else:
                return s + size

        result = []
        s = ctypes.addressof(self.aux.contents)
        while s < ctypes.addressof(self.aux.contents) + self.aux_len:
            c1 = ctypes.c_char.from_address(s).value.decode('utf-8')
            c2 = ctypes.c_char.from_address(s + 1).value.decode('utf-8')
            result.append('{}{}'.format(c1, c2))
            # print('result={}'.format(result), file=sys.stderr)
            s = skip_aux(s + 2)
        return result

    def __call__(self):
        result = collections.OrderedDict()
        for key in self.keys:
            # print('key={}'.format(key), file=sys.stderr)
            ptr = _bam_aux_get(self.struct_ptr, key.encode('utf-8'))
            # print('ptr={}'.format(ptr), file=sys.stderr)
            t = chr(ptr[0])
            # print('t={}'.format(t), file=sys.stderr)
            if t in ['c', 'C', 's', 'S', 'i', 'I']:
                result[key] = _bam_aux2i(ptr)
            elif t in ['d', 'f']:
                result[key] = _bam_aux2f(ptr)
            elif t in ['A']:
                result[key] = _bam_aux2A(ptr)
            elif t in ['Z', 'H']:
                result[key] = ctypes.string_at(
                    ctypes.addressof(_bam_aux2Z(ptr).contents)).decode('utf-')
            else:
                # TODO(holtgrewe): finish for arrays
                print('WARNING: unknown type {}'.format(t), file=sys.stderr)
                result[key] = None
        return result


class BAMRecordImpl:
    """Information extracted from C internals of ``BAMRecord``"""

    @staticmethod
    def from_struct(ptr, header):
        """Return ``BAMRecordImpl`` from internal C structure"""
        qname = _bam_get_qname(ptr).value.decode('utf-8')
        # print('QNAME={} ({})'.format(qname, type(qname)))
        flag = ptr[0].core.flag
        # print('FLAG={} ({})'.format(flag, type(flag)))
        tid = ptr[0].core.tid
        # print('TID={} ({})'.format(tid, type(tid)))
        ref = header.target_infos[ptr[0].core.tid].name
        # print('REF={} ({})'.format(ref, type(ref)))
        pos = ptr[0].core.pos
        # print('POS={} ({})'.format(pos, type(pos)))
        end_pos = _bam_endpos(ptr)
        # print('END_POS={} ({})'.format(end_pos, type(end_pos)))
        mapq = ptr[0].core.qual
        # print('QUAL={} ({})'.format(qual, type(qual)))
        cigar = BAMRecordImpl._cigar(ptr)
        # print('CIGAR={} ({})'.format(cigar, type(cigar)))
        mtid = ptr[0].core.mtid
        # print('MATE TID={} ({})'.format(mtid, type(mtid)))
        mref = header.target_infos[ptr[0].core.mtid].name
        # print('MATE REF={} ({})'.format(mref, type(mref)))
        mpos = ptr[0].core.mpos
        # print('MPOS={} ({})'.format(mpos, type(mpos)))
        isize = ptr[0].core.isize
        # print('ISIZE={} ({})'.format(isize, type(isize)))
        _l_qseq = ptr[0].core.l_qseq
        _seq_ptr = ctypes.cast(_bam_get_seq(ptr),
                               ctypes.POINTER(ctypes.c_uint8))
        seq = ''.join([_BAM_SEQ_STR[_bam_seqi(_seq_ptr, i)]
                       for i in range(_l_qseq)])
        # print('SEQ={} ({})'.format(seq, type(seq)))
        _qual_ptr = _bam_get_qual(ptr)
        qual = ''.join([chr(ord('!') + _qual_ptr[i]) for i in range(_l_qseq)])
        # print('QUAL={} ({})'.format(qual, type(qual)), file=sys.stderr)
        # TODO(holtgrewe): make parsing of tags lazy
        tags = BAMAuxTagParser(ptr, _bam_get_aux(ptr), _bam_get_l_aux(ptr))()
        # print('TAGS={} ({})'.format(tags, type(tags)))
        return BAMRecordImpl(qname, flag, tid, ref, pos, end_pos, mapq, cigar,
                             mtid, mref, mpos, isize, seq, qual, tags)

    @staticmethod
    def _cigar(ptr):
        def to_ce(cigar):
            return CIGARElement(_bam_cigar_oplen(cigar),
                                _bam_cigar_opchr(cigar))

        cigar_arr = _bam_get_cigar(ptr)
        return [to_ce(cigar_arr[i]) for i in range(ptr[0].core.n_cigar)]

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

    def __init__(self, struct_ptr=None, header=None, impl=None):
        #: pointer to wrapped C struct
        self.struct_ptr = struct_ptr
        #: wrapped C struct
        self.struct = None
        if self.struct_ptr:
            self.struct = self.struct_ptr[0]
        #: ``BAMHeader`` for references
        self.header = header
        #: ``BAMRecordImpl`` instance used for the representation
        self.impl = impl

    def detach(self):
        """Return copy that is detached from the underlying C object

        The only way of obtaining ``BAMRecords`` is through iterating
        ``BAMFile`` objects, either through an index or not.  For efficiency,
        the reading reuses the same buffer for reading.  If you want to keep a
        ``BAMRecord`` around for longer than the current iteration then you
        have to obtain a copy that is independent of the current buffer
        through the use of ``detach()``.
        """
        impl = self.impl
        if not impl:
            impl = BAMRecordImpl.from_struct(self.struct_ptr, self.header)
        self.impl = None
        return BAMRecord(impl=impl)

    def _reset(self):
        """Reset Python side, as if freshly constructed"""
        self.impl = None

    def __getattr__(self, name):
        """Delegation to self.impl if set, auto-set from struct"""
        if not self.impl and not self.struct:
            raise AttributeError('self.impl is None and cannot rebuild from '
                                 'None self.struct')
        elif not self.impl and self.struct:
            self.impl = BAMRecordImpl.from_struct(self.struct_ptr, self.header)
        return getattr(self.impl, name)


class BAMFileIter:
    """Iterate over a ``BAMFile``

    Do not use directly but by iterating over ``BAMFile``.  Iteration must
    be completed or ``close()`` must be called to prevent resource leaks.
    """

    def __init__(self, bam_file):
        #: the ``BAMFile`` to iterate through
        self.bam_file = bam_file
        #: buffer for readin in the file itself
        self.struct_ptr = _bam_init1()
        #: pointer to buffer for reading in the file record by record
        self.struct = self.struct_ptr[0]
        #: ``BAMRecord`` meant for consumption by the user
        self.record = BAMRecord(self.struct_ptr, self.bam_file.header)

    def __next__(self):
        r = _sam_read1(self.bam_file.struct_ptr,
                       self.bam_file.header.struct_ptr,
                       ctypes.byref(self.struct))
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


# TODO(holtgrewe): we probably want to differentiate BAM/CRAM and SAM.gz with
# two classes
class BAMIndexIter:
    """Iterate over query results from a ``BAMIndex``

    Do not use directly but by iterating over query result of``BAMIndex``.
    Iteration must be completed or ``close()`` must be called to prevent
    resource leaks.
    """

    def __init__(self, bam_index, itr):
        #: the ``BAMIndex`` to iterate through
        self.bam_index = bam_index
        #: the ``BAMFile`` used
        self.bam_file = self.bam_index.bam_file
        #: buffer for reading in the file itself
        self.struct_ptr = _bam_init1()
        #: pointer to buffer for reading in the file record by record
        self.struct = self.struct_ptr[0]
        #: pointer to iterator struct to for iteration
        self.itr_ptr = itr
        #: iterator struct to use for iteration
        self.itr = self.itr_ptr[0]
        #: ``BAMRecord`` meant for consumption by the user
        self.record = BAMRecord(self.struct_ptr, self.bam_file.header)
        # buffer to use in case of SAM.gz
        self._buffer = None
        if not self.bam_index.is_bam_or_cram:
            self._buffer = _kstring_t(0, 0, None)

    def __iter__(self):
        return self

    def __next__(self):
        if self.bam_index.is_bam_or_cram:
            r = _sam_itr_next(self.bam_file.struct_ptr,
                              self.itr_ptr,
                              ctypes.byref(self.struct))
        else:
            r = _tbx_itr_next(self.bam_file.struct_ptr,
                              self.bam_index.struct_ptr,
                              self.itr_ptr,
                              ctypes.byref(self._buffer))
        if r >= 0:
            # attempt to parse SAM, if SAM
            if not self.bam_index.is_bam_or_cram:
                _sam_parse1(ctypes.byref(self._buffer),
                            self.bam_file.header.struct_ptr,
                            ctypes.byref(self.struct))
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
        if self.struct_ptr:
            _bam_destroy1(self.struct_ptr)
            self.struct_ptr = None
            self.struct = None
        if self._buffer:
            self._buffer.free_p()
            self._buffer = None
        if self.itr_ptr and self.bam_index.is_bam_or_cram:
            _sam_itr_destroy(self.itr_ptr)
            self.itr_ptr = None
            self.itr = None
        elif self.itr_ptr and not self.bam_index.is_bam_or_cram:
            _tbx_itr_destroy(self.itr_ptr)
            self.itr_ptr = None
            self.itr = None


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

    @staticmethod
    def _get_index_ext(path):
        if path.endswith('.sam.gz'):
            return '.tbi'
        elif path.endswith('.bam'):
            return '.bai'
        else:
            tpl = 'Not a valid alignment file extension: {}'
            raise BAMIndexException(tpl.format(path))

    def __init__(self, path, bai_path=None, require_index=True,
                 auto_load=True, auto_build=False):
        #: path to BAM file
        self.path = path
        #: path to BAI (BAM index) file
        self.bai_path = bai_path or path + BAMIndex._get_index_ext(path)
        #: whether or not an index is required on construction
        self.require_index = require_index
        #: whether or not to automatically load the index
        self.auto_load = auto_load
        #: whether or not to automatically build the index
        self.auto_build = auto_build
        #: whether or not is BAM/CRAM (alternative is SAM+tabix)
        self.is_bam_or_cram = not self.path.endswith('.sam.gz')

        #: the ``BAMFile`` to use for reading
        self.bam_file = BAMFile(self.path)
        self.bam_file.open()

        #: wrapped C struct
        self.struct = None
        #: pointer to C struct
        self.struct_ptr = None

        self._check_auto_build()
        self._check_auto_load()
        self._check_require_index()

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
            tpl = 'Index {} required for {} for loading but not found.'
            raise BAMIndexException(tpl.format(self.bai_path, self.path))

    def _check_require_index(self):
        if not self.require_index:
            return
        if not os.path.exists(self.bai_path):
            tpl = 'BAM index required for {} but not found.'
            raise BAMIndexException(tpl.format(self.path))

    # TODO(holtgrewe): fix exception display if not region_string
    def query(self, region_str=None, seq=None, begin=None, end=None):
        if (region_str is None and
                (seq is None or begin is None or end is None)):
            raise BAMIndexException(
                'You have to either give region_str or seq/begin/end')
        if not self.is_bam_or_cram:
            if region_str:
                ptr = _tbx_itr_querys(self.struct_ptr,
                                      region_str.encode('utf-8'))
            else:
                ptr = _tbx_itr_queryi(self.struct_ptr, seq, begin, end)
        else:
            if region_str:
                ptr = _sam_itr_querys(self.struct_ptr,
                                      self.bam_file.header.struct_ptr,
                                      region_str.encode('utf-8'))
            else:
                ptr = _sam_itr_queryi(self.struct_ptr, seq, begin, end)
        if not ptr:
            tpl = 'Could not jump to {}'
            raise BAMIndexException(tpl.format(region_str))
        return BAMIndexIter(self, ptr)

    def load(self):
        self.close(close_file=False)
        if not self.is_bam_or_cram and self.bai_path:
            self.struct_ptr = _tbx_index_load2(self.path.encode('utf-8'),
                                               self.bai_path.encode('utf-8'))
            if not self.struct_ptr:
                tpl = 'Could not load tabix index {} for {}'
                raise BAMIndexException(tpl.format(self.bai_path, self.path))
        elif not self.is_bam_or_cram and not self.bai_path:
            self.struct_ptr = _tbx_index_load(self.path.encode('utf-8'))
            if not self.struct_ptr:
                tpl = 'Could not load tabix index for {}'
                raise BAMIndexException(tpl.format(self.path))
        elif self.bai_path:
            self.struct_ptr = _sam_index_load2(self.bam_file.struct_ptr,
                                               self.path.encode('utf-8'),
                                               self.bai_path.encode('utf-8'))
            if not self.struct_ptr:
                tpl = 'Could not load BAM/CRAM index {} for {}'
                raise BAMIndexException(tpl.format(self.bai_path, self.path))
        else:
            self.struct_ptr = _sam_index_load(self.bam_file.struct_ptr,
                                              self.path.encode('utf-8'))
            if not self.struct_ptr:
                tpl = 'Could not load BAM/CRAM index for {}'
                raise BAMIndexException(tpl.format(self.path))
        self.struct = self.struct_ptr[0]

    def close(self, close_file=True):
        if self.struct_ptr:
            if self.path.endswith('.sam.gz'):
                _tbx_destroy(self.struct_ptr)
            else:
                _hts_idx_destroy(self.struct_ptr)
        self.struct = None
        self.struct_ptr = None
        if close_file:
            self.bam_file.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close(close_file=True)
