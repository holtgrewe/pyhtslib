#!/usr/bin/env python3
"""Access to BCF and BCF files"""

import collections
# import ctypes  # TODO(holtgrew): remove?
# import logging
# import os
# import os.path

from pyhtslib.hts_internal import *  # NOQA
from pyhtslib.bcf_internal import *  # NOQA

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'


class BCFIndexException(Exception):
    """Raised when there is a problem with a BCFIndex file"""


class BCFFileException(Exception):
    """Raised when there is a problem with a BCFFile file"""


class BCFInfoFieldException(Exception):
    """Raised on problems in parsing ``INFO`` fields from ``BCFRecord``s"""

    def __init__(self, res, key, type_):
        MSGS = {
            -1: 'no such INFO tag defined in the header',
            -2: ('clash between types defined n the header and encountered '
                 'in the BCF record'),
            -3: 'tag is not present in the BCF record',
        }
        msg = MSGS.get(res, 'Unknown problem with BCF info field.')
        tpl = 'Problem with {} INFO field {}: {}'
        Exception.__init__(self, tpl.format(type_, key, msg))


class BCFHeaderTargetInfo:
    """Information (name, length) for the reference/target sequence"""

    def __init__(self, name, length):
        #: name of the target/reference
        self.name = name
        #: length of the target/reference
        self.length = int(length)

    def __eq__(self, other):
        return (self.name == other.name and self.length == other.length)

    def __repr__(self):
        return 'BCFHeaderTargetInfo({}, {})'.format(
            repr(self.name), self.length)


class BCFGenericHeaderRecord:
    """Generic BCF header record with any content"""

    def __init__(self, tag, value):
        #: tag of the generic record
        self.tag = key
        #: value of the generic record
        self.value = value

    def to_header_line(self):
        return '##{}={}'.format(self.key, self.value)


class BCFStructuredHeaderRecord(BCFGenericHeaderRecord):
    """Structured BCF header record (``TAG=<A=...,B=...>``)"""

    def __init__(self, tag, entries=[], repr_entries=[]):
        #: ``str`` with tag name
        self.tag = tag
        #: ``collections.OrderedDict`` with ``str``-to-``str`` mapping
        self.entries = collections.OrderedDict(entries)
        # set of keys in self.entries where we want to use repr but rather
        # plain string conversion
        self.repr_entries = set(repr_entries)

    def to_header_line(self):
        def hdr_repr(key, val):
            """Representation of ``val`` in BCF header line"""
            if key not in self.repr_entries:
                return str(val)
            elif hasattr(val, 'encode'):
                return val.encode('string_escape').replace('"', r'\"')
            else:
                return repr(val)

        tuple_ = ','.join(['{}={}'.format(k, hdr_repr(key, v))
                           for k, v in self.entries.items() if v])
        return '##{tag}=<{tuple_}>'.format(tag=self.tag, tuple_=tuple_)


class BCFFilterHeaderRecord(BCFStructuredHeaderRecord):
    """BCF header record describing FILTER values"""

    def __init__(self, id_, description):
        BCFStructuredHeaderRecord.__init__(
            self, 'FILTER', [('ID', id_), ('Description', description)],
            ['Description'])


class BCFInfoHeaderRecord(BCFStructuredHeaderRecord):
    """BCF header record describing INFO values"""

    def __init__(self, id_, number, type_, description):
        BCFStructuredHeaderRecord.__init__(
            self, 'INFO',
            [('ID', id_), ('Number', number), ('Type', type_),
             ('Description', description)],
            ['Description'])


class BCFFormatHeaderRecord(BCFStructuredHeaderRecord):
    """BCF header record describing FORMAT values"""

    def __init__(self, id_, number, type_, description):
        BCFStructuredHeaderRecord.__init__(
            self, 'FORMAT',
            [('ID', id_), ('Number', number), ('Type', type_),
             ('Description', description)],
            ['Description'])


class BCFContigHeaderRecord(BCFStructuredHeaderRecord):
    """BCF header record describing a contig line"""

    def __init__(self, id_, length, other_entries=[], escape_args=[]):
        BCFStructuredHeaderRecord.__init__(
            self, 'contig',
            [('ID', id_), ('Length', int(length))] + other_entries,
            escape_args)


class BCFHeader:
    """The information stored in the header of a BCF file

    All information is completely extracted from the information in the BCF
    file and available as native Python objects.
    """

    @staticmethod
    def _read_from_file(file_ptr):
        """Read header from BCF file handle"""

    def __init__(self, struct_ptr=None):
        """Initialize object from pointer to htslib BCF header type

        When no argument is given, an empty header is constructed.
        """

        #: pointer to htslib BCF representation, if any
        self.struct_ptr = struct_ptr
        #: internal htslib BCF representation, if any
        self.struct = None if not self.struct_ptr else self.struct_ptr[0]

        #: VCF version, defaulting to ``"VCFv4.2"``
        self.vcf_version = "VCFv4.2"
        #: list of strings with sample IDs, defaulting to ``[]``
        self.sample_ids = []
        #: list of ``BCFHeaderTargetInfo``s, defaulting to ``[]``
        self.target_infos = []
        #: mapping of string key to numeric id, defaulting to ``{}``
        self.key_ids = {}
        #: list of BCFHeaderRecord objects, defaulting to ``[]``
        self.header_records = []
        #: dict mapping string ID keys to FILTER BCFHeaderRecord objects
        self.id_to_filter_record = {}
        #: dict mapping string ID keys to INFO BCFHeaderRecord objects
        self.id_to_info_record = {}
        #: dict mapping string ID keys to FORMAT BCFHeaderRecord objects
        self.id_to_format_record = {}
        #: dict mapping string ID keys to CONTIG BCFHeaderRecord objects
        self.id_to_contig_record = {}

        self._fill_from_struct()

    def add_header_record(self, record):
        """Add a BCFHeaderRecord to ``self.header_records``

        Also updates ``self.id_to_*`` fields appropriately.
        """
        if record.key in ['FILTER', 'INFO', 'FORMAT']:
            self.key_ids.setdefault(self.record.tag, len(self.header.key_ids))
        elif self.header.key == 'contig':
            self.target_infos.append(BCFHeaderTargetInfo(
                record.entries['ID'], record.entries['Length']))
        self.header_records.append(self.record)

    def free(self):
        """Free all memory associated with this object"""
        if self.struct_ptr:
            _bcf_hdr_destroy(self.struct_ptr)
            self.struct_ptr = None
            self.struct = None

    def _fill_from_struct(self):
        """Fill object with the information from ``self.struct_ptr``, if any
        """
        if not self.struct_ptr:
            return
        raise Exception('Implement me!')


class BAMRecordImpl:
    """Information extracted from C internals of ``BAMRecord``"""

    @staticmethod
    def from_struct(ptr, header):
        # TODO(holtgrewe): write me!
        return BAMRecord(ptr, header)

    def __init__(self, *args, **kwargs):
        # TODO(holtgrewe): write me!
        pass


class BCFRecord:
    """Record from a BCF file"""

    def __init__(self, struct_ptr=None, header=None, impl=None):
        #: pointer to wrapped C struct
        self.struct_ptr = struct_ptr
        #: wrapped C struct
        self.struct = None
        if self.struct_ptr:
            self.struct = self.struct_ptr[0]
        #: ``BCFHeader`` for references
        self.header = header
        #: ``BCFRecordImpl`` instance used for the representation
        self.impl = impl

    def detach(self):
        """Return copy that is detached from the underlying C object

        The only way of obtaining ``BCFRecord``s is through iterating
        ``BCFFile`` objects, either through an index or not.  For efficiency,
        the reading reuses the same buffer for reading.  If you want to keep a
        ``BCFRecord`` around for longer than the current iteration then you
        have to obtain a copy that is independent of the current buffer
        through the use of ``detach()``.
        """
        impl = self.impl
        if not impl:
            impl = BCFRecordImpl.from_struct(self.struct_ptr, self.header)
        self.impl = None
        return BCFRecord(impl=impl)

    def _reset(self):
        """Reset Python side, as if freshly constructed"""
        self.impl = None

    def __getattr__(self, name):
        """Delegation to self.impl if set, auto-set from struct"""
        if not self.impl and not self.struct:
            raise AttributeError('self.impl is None and cannot rebuild from '
                                 'None self.struct')
        elif not self.impl and self.struct:
            self.impl = BCFRecordImpl.from_struct(self.struct_ptr, self.header)
        return getattr(self.impl, name)


class BCFFileIter:
    """Iterate over a ``BCFFile``

    Do not use directly but by iterating over ``BCFFile``.  Iteration must
    be completed or ``close()`` must be called to prevent resource leaks.
    """

    def __init__(self, bcf_file):
        #: the ``BCFFile`` to iterate through
        self.bcf_file = bcf_file
        #: buffer for readin in the file itself
        self.struct_ptr = _bcf_init1()
        #: pointer to buffer for reading in the file record by record
        self.struct = self.struct_ptr[0]
        #: ``BCFRecord`` meant for consumption by the user
        self.record = BCFRecord(self.struct_ptr, self.bcf_file.header)
        # whether or not iterating over BCF file
        self.is_bcf = (self.bcf_file.format == 'BCF')

    def __next__(self):
        read = _bcf_read1 if self.is_bcf else _vcf_read1
        r = read(self.bcf_file.struct_ptr, self.bcf_file.header.struct_ptr,
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
                raise BCFFileException(tpl.format(self.bcf_file.path))
            else:
                self.close()
                raise StopIteration

    def close(self):
        if self.struct_ptr:
            _bcf_destroy1(self.struct_ptr)
            self.struct_ptr = None
            self.struct = None


class BCFFile:
    """Representation of a VCF/BCF file"""

    def __init__(self, path, mode='r', worker_threads=1):
        #: path to BCF file
        self.path = path
        #: mode to open file with
        self.mode = mode
        #: number of worker threads
        self.worker_threads = worker_threads
        #: wrapped C struct
        self.struct = None
        #: pointer to C struct
        self.struct_ptr = None
        #: representation of BAM header
        self.header = None

        # collection of iterators, we will call close() on all of them
        # in our own close to ensure that all memory is freed
        self.iterators = []

    @property
    def file_format(self):
        """Return file format, either 'VCF', 'BCF', or ``None``"""
        if not self.struct:
            return None
        DICT = {_HTSExactFormat.VCF: 'VCF', _HTSExactFormat.BCF: 'BCF'}
        return DICT.get(self.struct.ftype.format)

    def open(self):
        """Open file and read header"""
        if self.struct_ptr:
            return  # already open
        # open file and store handles
        self.struct_ptr = _hts_open(self.path.encode('utf-8'), self.mode)
        self.struct = self.struct_ptr[0]
        # check file format
        if self.file_format not in ['VCF', 'BCF']:
            self.close()
            raise BCFFileException('Not a VCF/BCF file: {}'.format(self.path))
        # read header
        self.header = BCFHeader._read_from_file(self.struct_ptr)

    def close(self):
        """Close file again and free header and other data structures

        This function is idempotent
        """
        if self.header:
            self.header.free()
            self.header = None
        if self.struct_ptr:
            _hts_close(self.struct_ptr)
            self.struct_ptr = None
            self.struct = None
        for it in self.iterators:
            it.close()
        self.iterators = []

    def __iter__(self):
        self.iterators.append(BCFFileIter(self))
        return self.iterators[-1]

    def __enter__(self):
        self.open()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
