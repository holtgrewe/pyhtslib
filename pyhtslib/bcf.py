#!/usr/bin/env python3
"""Access to BCF and BCF files"""

import collections
import ctypes  # TODO(holtgrew): remove?

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

    @staticmethod
    def _from_struct(struct):
        """Construct BCFGenericHeaderRecord from ``_bcf_hrec_t``"""
        return BCFGenericHeaderRecord(struct.key.decode('utf-8'),
                                      struct.value.decode('utf-8'))

    def __init__(self, tag, value):
        #: tag of the generic record
        self.tag = tag
        #: value of the generic record
        self.value = value

    def __eq__(self, other):
        return (self.tag == other.tag and self.value == other.value)

    def to_header_line(self):
        return '##{}={}'.format(self.tag, self.value)


class BCFStructuredHeaderRecord(BCFGenericHeaderRecord):
    """Structured BCF header record (``TAG=<A=...,B=...>``)"""

    @staticmethod
    def _from_struct_impl(klass, required, struct, repr_entries):
        """Construct klass object from ``_bcf_hrec_t``

        Used by subtypes only
        """
        tag = struct.key.decode('utf-8')
        entries = dict([(struct.keys[i].decode('utf-8'),
                         struct.vals[i].decode('utf-8'))
                        for i in range(struct.nkeys)])
        keys = dict([(k.lower(), k) for k in entries.keys()])
        for key in required:
            if not keys.get(key):
                tpl = ('FILTER record does not have entry for {} column (case '
                       'insensitive lookup) ')
                raise BCFFileException(tpl.format(key.upper()))
        required_list = [entries[keys[k]] for k in required]
        other_list = [(struct.keys[i].decode('utf-8'),
                       struct.vals[i].decode('utf-8'))
                      for i in range(struct.nkeys)
                      if struct.keys[i].decode('utf-').lower()
                      not in required]
        tpl = 'tag={}, required_list={}, other_entries={}, repr_entries={}'
        return klass(tag, *required_list, other_entries=other_list,
                     repr_entries=repr_entries)

    @staticmethod
    def _from_struct(struct):
        """Construct BCFStructuredHeaderRecord from ``_bcf_hrec_t``"""
        tag = struct.key
        entries = [(struct.keys[i], struct.vals[i])
                   for i in range(struct.nkeys)]
        return BCFStructuredHeaderRecord(tag, entries)

    def __init__(self, tag, entries=[], repr_entries=[]):
        #: ``str`` with tag name
        self.tag = tag
        #: ``collections.OrderedDict`` with ``str``-to-``str`` mapping
        self.entries = collections.OrderedDict(entries)
        # set of keys in self.entries where we want to use repr but rather
        # plain string conversion
        self.repr_entries = set(repr_entries)

    def __eq__(self, other):
        return (self.tag == other.tag and self.entries == other.entries and
                self.repr_entries == other.repr_entries)

    def to_header_line(self):
        def hdr_repr(key, val):
            """Representation of ``val`` in BCF header line"""
            if key not in self.repr_entries:
                return str(val)
            elif hasattr(val, 'encode'):
                return '"{}"'.format(val.replace('"', r'\"'))
            else:
                return repr(val)

        tuple_ = ','.join(['{}={}'.format(k, hdr_repr(k, v))
                           for k, v in self.entries.items() if v])
        return '##{tag}=<{tuple_}>'.format(tag=self.tag, tuple_=tuple_)


class BCFFilterHeaderRecord(BCFStructuredHeaderRecord):
    """BCF header record describing FILTER values"""

    @staticmethod
    def _from_struct(struct):
        return BCFFilterHeaderRecord._from_struct_impl(
            BCFFilterHeaderRecord, ['id', 'description'], struct,
            ['description'])

    def __init__(self, tag, id_, description, other_entries=[],
                 repr_entries=[]):
        BCFStructuredHeaderRecord.__init__(
            self, 'FILTER',
            [('id', id_), ('description', description)] + list(other_entries),
            ['description'] + repr_entries)


class BCFInfoHeaderRecord(BCFStructuredHeaderRecord):
    """BCF header record describing INFO values"""

    @staticmethod
    def _from_struct(struct):
        return BCFInfoHeaderRecord._from_struct_impl(
            BCFInfoHeaderRecord, ['id', 'number', 'type', 'description'],
            struct, ['description'])

    def __init__(self, tag, id_, number, type_, description, other_entries=[],
                 repr_entries=[]):
        BCFStructuredHeaderRecord.__init__(
            self, 'INFO',
            [('id', id_), ('number', number), ('type', type_),
             ('description', description)] + list(other_entries),
            ['description'] + repr_entries)


class BCFFormatHeaderRecord(BCFStructuredHeaderRecord):
    """BCF header record describing FORMAT values"""

    @staticmethod
    def _from_struct(struct):
        return BCFFormatHeaderRecord._from_struct_impl(
            BCFFormatHeaderRecord, ['id', 'number', 'type', 'description'],
            struct, ['description'])

    def __init__(self, tag, id_, number, type_, description, other_entries=[],
                 repr_entries=[]):
        BCFStructuredHeaderRecord.__init__(
            self, 'FORMAT',
            [('id', id_), ('number', number), ('type', type_),
             ('description', description)] + list(other_entries),
            ['description'] + repr_entries)


class BCFContigHeaderRecord(BCFStructuredHeaderRecord):
    """BCF header record describing a contig line"""

    @staticmethod
    def _from_struct(struct):
        return BCFContigHeaderRecord._from_struct_impl(
            BCFContigHeaderRecord, ['id', 'length'],
            struct, ['description'])

    def __init__(self, tag, id_, length, other_entries=[], repr_entries=[]):
        BCFStructuredHeaderRecord.__init__(
            self, 'contig',
            [('id', id_), ('length', int(length))] + list(other_entries),
            repr_entries)


class BCFHeader:
    """The information stored in the header of a BCF file

    All information is completely extracted from the information in the BCF
    file and available as native Python objects.
    """

    @staticmethod
    def _read_from_file(file_ptr):
        """Read header from BCF file handle"""
        ptr = _bcf_hdr_read(file_ptr)
        if not ptr:
            raise BCFFileException('Could not load BCF header from file!')
        return BCFHeader(ptr)

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
        #: list of strings with sample names, defaulting to ``[]``
        self.sample_names = []
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
        if record.tag in ['FILTER', 'INFO', 'FORMAT']:
            self.key_ids.setdefault(record.tag, len(self.key_ids))
            DICT = {
                'FILTER': self.id_to_filter_record,
                'INFO': self.id_to_info_record,
                'FORMAT': self.id_to_format_record,
            }
            DICT[record.tag][record.entries['id']] = record
        elif record.tag == 'contig':
            self.target_infos.append(BCFHeaderTargetInfo(
                record.entries['id'], record.entries['length']))
            self.id_to_contig_record[record.entries['id']] = record
        self.header_records.append(record)

    def free(self):
        """Free all memory associated with this object"""
        if self.struct_ptr:
            _bcf_hdr_destroy(self.struct_ptr)
            self.struct_ptr = None
            self.struct = None

    def to_vcf_header(self):
        """Return multi-line string with VCF header representation"""
        lines = [BCFGenericHeaderRecord('fileformat', self.vcf_version)]
        for r in self.header_records:
            if r.tag == 'fileformat':
                lines = []
        for r in self.header_records:
            lines.append(r)
        return '\n'.join([l.to_header_line() for l in lines] + [''])

    def _fill_from_struct(self):
        """Fill object with the information from ``self.struct_ptr``, if any
        """
        if not self.struct_ptr:
            return
        # get sequence names, lengths will be parsed from the headers below
        num_seqs = ctypes.c_int(0)
        arr_ptr = _bcf_hdr_seqnames(self.struct_ptr, ctypes.byref(num_seqs))
        num_seqs = num_seqs.value
        seq_names = [arr_ptr[i].decode('utf-8') for i in range(num_seqs)]
        _libc.free(arr_ptr)
        # get sample names
        num_samples = _bcf_hdr_nsamples(self.struct_ptr)
        self.sample_names = [
            _bcf_hdr_get_sample_name(self.struct_ptr, i).decode('utf-8')
            for i in range(num_samples)]
        # get mapping from FILTER/INFO/FORMAT key to numeric id
        num_ids = _bcf_hdr_nids(self.struct_ptr)
        key_to_id = [
            _bcf_hdr_int2id(self.struct_ptr, _BCF_DT_ID, i).decode('utf-8')
            for i in range(num_ids)]
        self.key_ids = dict([(key, i) for i, key in enumerate(key_to_id)])
        # parse out headers
        self._parse_headers_from_struct()
        # create self.target_infos
        for seq_name in seq_names:
            ctg = self.id_to_contig_record.get(seq_name)
            if not ctg:
                self.target_infos.append(BCFHeaderTargetInfo(seq_name, 0))

    def _parse_headers_from_struct(self):
        """Fill the header member fields from ``self.struct_ptr``

        ``self.struct_ptr`` must not be ``None``.
        """
        # mapping of bcf header record id to BCFGenericHeaderRecord subclass
        DIR = {
            _BCF_HL_FLT: BCFFilterHeaderRecord,
            _BCF_HL_INFO: BCFInfoHeaderRecord,
            _BCF_HL_FMT: BCFFormatHeaderRecord,
            _BCF_HL_CTG: BCFContigHeaderRecord,
            _BCF_HL_STR: BCFStructuredHeaderRecord,
            _BCF_HL_GEN: BCFGenericHeaderRecord,
        }
        for i in range(self.struct.nhrec):
            rec = self.struct.hrec[i][0]  # bcf_hrec_t
            self.add_header_record(DIR[rec.type]._from_struct(rec))


class BCFRecordImpl:
    """Information extracted from C internals of ``BCFRecord``"""

    @staticmethod
    def from_struct(ptr, header):
        # TODO(holtgrewe): write me!
        return BCFRecord(ptr, header)

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
