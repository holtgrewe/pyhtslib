#!/usr/bin/env python3
"""Access to BCF and BCF files"""

import collections
import ctypes  # TODO(holtgrew): remove?
import sys  # NOQA  # TODO(holtgrew): remove?

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
        #: list of ids
        self.ids = []
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
        self.ids = [
            _bcf_hdr_int2id(self.struct_ptr, _BCF_DT_ID, i).decode('utf-8')
            for i in range(num_ids)]
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


class BCFGenotype:
    """Representation of a genotype in a BCF/VCF record"""


class _BCFTypedInfoFieldConverter:
    # TODO(holtgrew): interpret as scalar if count is 1

    def __init__(self, header, info, struct, key):
        self.header = header
        self.info = info
        self.struct = struct
        self.key = key

    def __call__(self):
        # ignore any values if the header says that there are none; otherwise
        # we can get problems with things like "SNP=true"
        num = _bcf_hdr_id2int(self.header.struct_ptr, _BCF_DT_ID,
                              self.key.encode('utf-8'))
        length = _bcf_hdr_id2length(self.header.struct_ptr, _BCF_HL_INFO, num)
        number = _bcf_hdr_id2number(self.header.struct_ptr, _BCF_HL_INFO, num)
        type_ = _bcf_hdr_id2type(self.header.struct_ptr, _BCF_HL_INFO, num)
        if length == _BCF_VL_FIXED and number == 0:
            return self.convert_flag()

        # branch based on the type of the info field
        BRANCH = {
            _BCF_BT_NULL: self.convert_flag,
            _BCF_BT_INT8: self.convert_int,
            _BCF_BT_INT16: self.convert_int,
            _BCF_BT_INT32: self.convert_int,
            _BCF_BT_FLOAT: self.convert_float,
            _BCF_BT_CHAR: self.convert_char,
        }
        result = BRANCH[self.info.type]()

        # extract single value in case of scalars and split strings at comma
        # in case of vectors
        is_scalar = (length == _BCF_VL_FIXED and number == 1)
        if type_ == _BCF_HT_STR:
            return result if is_scalar else result.split(',')
        elif is_scalar:
            assert len(result) == 1, 'len(result) == {}'.format(len(result))
            return collections.OrderedDict(result[0])
        else:
            return collections.OrderedDict(result)

    def convert_flag(self):
        ndst = ctypes.c_int(0)
        dst = ctypes.c_char_p()
        res = _bcf_get_info_flag(
            self.header.struct_ptr, ctypes.byref(self.struct),
            self.key.encode('utf-8'), ctypes.byref(dst), ctypes.byref(ndst))
        return (res == 1)

    def convert_int(self):
        ndst = ctypes.c_int(0)
        dst = ctypes.POINTER(ctypes.c_int32)()
        res = _bcf_get_info_int32(
            self.header.struct_ptr, ctypes.byref(self.struct),
            self.key.encode('utf-8'), ctypes.byref(dst), ctypes.byref(ndst))
        if res < 0:
            raise BCFInfoFieldException(res, self.key, 'int32')
        return [dst[i] for i in range(ndst.value)]

    def convert_float(self):
        ndst = ctypes.c_int(0)
        dst = ctypes.POINTER(ctypes.c_float)()
        res = _bcf_get_info_float(
            self.header.struct_ptr, ctypes.byref(self.struct),
            self.key.encode('utf-8'), ctypes.byref(dst), ctypes.byref(ndst))
        if res < 0:
            raise BCFInfoFieldException(res, self.key, 'float')
        return [dst[i] for i in range(ndst.value)]

    def convert_char(self):
        ndst = ctypes.c_int(0)
        dst = ctypes.c_char_p()
        res = _bcf_get_info_string(
            self.header.struct_ptr, ctypes.byref(self.struct),
            self.key.encode('utf-8'), ctypes.byref(dst), ctypes.byref(ndst))
        if res < 0:
            raise BCFInfoFieldException(res, self.key, 'string')
        return dst.value.decode('utf-8')


class BCFRecordImpl:
    """Information extracted from C internals of ``BCFRecord``"""

    @staticmethod
    def from_struct(struct, header):
        """Return ``BCFRecordImpl`` from internal C structure"""
        r_id = struct.rid
        print('r_id={} ({})'.format(r_id, type(r_id)))
        chrom = header.target_infos[r_id].name
        print('chrom={} ({})'.format(chrom, type(chrom)))
        begin_pos = struct.pos
        print('begin_pos={} ({})'.format(begin_pos, type(begin_pos)))
        end_pos = struct.pos + struct.rlen
        print('end_pos={} ({})'.format(end_pos, type(end_pos)))
        tmp_id = struct.d.id
        print('tmp_id={} ({})'.format(tmp_id, type(tmp_id)))
        ids = struct.d.id.decode('utf-8').split(';')
        print('ids={} ({})'.format(ids, type(ids)))
        ref = struct.d.allele[0].decode('utf-8')
        print('ref={} ({})'.format(ref, type(ref)))
        alts = [struct.d.allele[i].decode('utf-8')
                for i in range(1, struct.n_allele)]
        print('alts={} ({})'.format(alts, type(alts)))
        qual = struct.qual
        print('qual={} ({})'.format(qual, type(qual)))
        filters = [header.ids[i] for i in range(struct.d.n_flt)]
        print('filters={} ({})'.format(filters, type(filters)))
        info = BCFRecordImpl._build_info_field(struct, header)
        print('info={} ({})'.format(info, type(info)))
        format_ = None
        print('format_={} ({})'.format(format_, type(format_)))
        genotypes = None
        print('genotypes={} ({})'.format(genotypes, type(genotypes)))

        return BCFRecordImpl(r_id, chrom, begin_pos, end_pos, ids, ref,
                             alts, qual, filters, info, format_, genotypes)

    def _build_info_field(struct, header):
        """Return INFO column entries."""
        result = []
        for i in range(struct.n_info):
            info = struct.d.info[i]
            key_str = header.ids[info.key]
            converter = _BCFTypedInfoFieldConverter(
                header, info, struct, key_str)
            result.append((key_str, converter()))
        return result

    def __init__(self, r_id, chrom, begin_pos, end_pos, ids, ref, alts,
                 qual, filters, info, format_, genotypes):
        #: target sequence id of variant (cmp. CHROM)
        self.r_id = r_id
        #: target sequence name (cmp. CHROM)
        self.chrom = chrom
        #: begin position of the variant in the reference (POS)
        self.begin_pos = begin_pos
        #: end position of the variant in the reference (POS)
        self.end_pos = end_pos
        #: list of uniq ids of the (e.g. dbSNP or COSMIC, cmp. ID), split at
        #: ``';'``
        self.ids = list(ids)
        #: reference sequence (REF)
        self.ref = ref
        #: alternative alleles, split at ``','``
        self.alts = list(alts)
        #: list of all alleles
        self.alleles = [self.ref] + self.alts
        #: alignment quality
        self.qual = qual
        #: filter entries, split at ``';'``
        self.filters = list(filters)
        #: ``OrdereDictionary`` with values from the INFO column
        self.info = info
        #: list of ids for FORMAT column
        self.format = format_
        #: list of ``BCFGenotype``, one for each genotype
        self.genotypes = genotypes


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
            self.impl = BCFRecordImpl.from_struct(
                self.struct_ptr[0], self.header)
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
        self.is_bcf = (self.bcf_file.file_format == 'BCF')

    def __next__(self):
        read = _bcf_read1 if self.is_bcf else _vcf_read1
        r = read(self.bcf_file.struct_ptr, self.bcf_file.header.struct_ptr,
                 self.struct_ptr)
        if r >= 0:
            r = _bcf_unpack(self.struct_ptr, _BCF_UN_ALL)
            # TODO(holtgrewe): check r from unpack?
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
