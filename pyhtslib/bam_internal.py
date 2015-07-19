#!/usr/bin/env python3
"""C structs and routines for SAM/BAM access"""

import ctypes

import pyhtslib.load_dll as pl
from pyhtslib.hts_internal import *  # NOQA

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'

# export everything from this submodule manually, including the code that
# starts with an underscore, importing modules will not import the latter

__all__ = [
    # constants
    '_BAM_CMATCH',
    '_BAM_CINS',
    '_BAM_CDEL',
    '_BAM_CREF_SKIP',
    '_BAM_CSOFT_CLIP',
    '_BAM_CHARD_CLIP',
    '_BAM_CPAD',
    '_BAM_CEQUAL',
    '_BAM_CDIFF',
    '_BAM_CBACK',

    '_BAM_CIGAR_STR',
    '_BAM_CIGAR_SHIFT',
    '_BAM_CIGAR_MASK',
    '_BAM_CIGAR_TYPE',

    '_BAM_FPAIRED',
    '_BAM_FPROPER_PAIR',
    '_BAM_FUNMAP',
    '_BAM_FMUNMAP',
    '_BAM_FREVERSE',
    '_BAM_FMREVERSE',
    '_BAM_FREAD1',
    '_BAM_FREAD2',
    '_BAM_FSECONDARY',
    '_BAM_FQCFAIL',
    '_BAM_FDUP',
    '_BAM_FSUPPLEMENTARY',

    # macros
    '_bam_cigar_op',
    '_bam_cigar_oplen',
    '_bam_cigar_opchr',
    '_bam_cigar_gen',
    '_bam_cigar_type',

    '_bam_is_rev',
    '_bam_is_mrev',
    '_bam_get_qname',
    '_bam_get_cigar',
    '_bam_get_seq',
    '_bam_get_qual',
    '_bam_get_aux',
    '_bam_get_l_aux',
    '_bam_seqi',
    '_bam_itr_destroy',
    '_bam_itr_queryi',
    '_bam_itr_querys',
    '_bam_itr_next',
    '_bam_index_load',
    '_bam_index_build',
    '_sam_itr_destroy',
    '_sam_itr_next',
    '_sam_open',
    '_sam_close',

    # structures
    '_bam_hdr_t',
    '_bam1_core_t',
    '_bam1_t',
    '_samFile',

    # functions
    '_bam_hdr_init',
    '_bam_hdr_read',
    '_bam_hdr_write',
    '_bam_hdr_destroy',
    '_bam_name2id',
    '_bam_hdr_dup',
    '_bam_init1',
    '_bam_destroy1',
    '_bam_read1',
    '_bam_write1',
    '_bam_copy1',
    '_bam_dup1',
    '_bam_endpos',
    '_bam_str2flag',
    '_bam_flag2str',
    '_bam_flag2str',
    '_sam_index_load',
    '_sam_index_load2',
    '_sam_index_build',
    '_sam_index_build2',
    '_sam_itr_queryi',
    '_sam_itr_querys',
    '_sam_hdr_parse',
    '_sam_hdr_read',
    '_sam_hdr_write',
    '_sam_parse1',
    '_sam_format1',
    '_sam_read1',
    '_sam_write1',
    '_bam_aux_get',
    '_bam_aux2i',
    '_bam_aux2f',
    '_bam_aux2A',
    '_bam_aux2Z',
    '_bam_aux_append',
    '_bam_aux_del',
]

# ----------------------------------------------------------------------------
# Constants and Macros
# ----------------------------------------------------------------------------

# CIGAR-related constants

_BAM_CMATCH = 0
_BAM_CINS = 1
_BAM_CDEL = 2
_BAM_CREF_SKIP = 3
_BAM_CSOFT_CLIP = 4
_BAM_CHARD_CLIP = 5
_BAM_CPAD = 6
_BAM_CEQUAL = 7
_BAM_CDIFF = 8
_BAM_CBACK = 9

_BAM_CIGAR_STR = "MIDNSHP=XB"
_BAM_CIGAR_SHIFT = 4
_BAM_CIGAR_MASK = 0xf
_BAM_CIGAR_TYPE = 0x3C1A7


def _bam_cigar_op(c):
    return c & _BAM_CIGAR_MASK


def _bam_cigar_oplen(c):
    return c >> _BAM_CIGAR_SHIFT


def _bam_cigar_opchr(c):
    return _BAM_CIGAR_STR[_bam_cigar_op(c)]


def _bam_cigar_gen(l, o):
    return (l << _BAM_CIGAR_SHIFT) | o


def _bam_cigar_type(o):
    # bit 1: consume query; bit 2: consume reference
    return (_BAM_CIGAR_TYPE >> (o << 1)) & 3

# the read is paired in sequencing, no matter whether it is mapped in a pair
_BAM_FPAIRED = 1
# the read is mapped in a proper pair
_BAM_FPROPER_PAIR = 2
# the read itself is unmapped; conflictive with BAM_FPROPER_PAIR
_BAM_FUNMAP = 4
# the mate is unmapped
_BAM_FMUNMAP = 8
# the read is mapped to the reverse strand
_BAM_FREVERSE = 16
# the mate is mapped to the reverse strand
_BAM_FMREVERSE = 32
# this is read1
_BAM_FREAD1 = 64
# this is read2
_BAM_FREAD2 = 128
# not primary alignment
_BAM_FSECONDARY = 256
# QC failure
_BAM_FQCFAIL = 512
# optical or PCR duplicate
_BAM_FDUP = 1024
# supplementary alignment
_BAM_FSUPPLEMENTARY = 2048


def _bam_is_rev(b):
    return (b[0].core.flag & _BAM_FREVERSE) != 0


def _bam_is_mrev(b):
    return (b[0].core.flag & _BAM_FMREVERSE) != 0


def _bam_get_qname(b):
    return ctypes.cast(b[0].data, ctypes.c_char_ptr)


def _bam_get_cigar(b):
    return ctypes.cast(b[0].data + b[0].core.l_qname,
                       ctypes.POINTER(ctypes.c_uint32))


def _bam_get_seq(b):
    return b[0].data + (b[0].core.n_cigar << 2) + b[0].core.l_qname


def _bam_get_qual(b):
    return _bam_get_seq + ((b[0].core.l_qseq) >> 1)


def _bam_get_aux(b):
    return _bam_get_qual + b[0].core.l_seq


def _bam_get_l_aux(b):
    return (b[0].l_data - (b[0].core.n_cigar << 2) - b[0].core.l_qname -
            b[0].core.l_qseq - ((b[0].core.l_qseq + 1) >> 1))


def _bam_seqi(s, i):
    return s[i >> 1] >> ((((~i) & 1) << 2) & 0xf)


def _bam_itr_destroy(it):
    return _hts_itr_destroy(it)


def _bam_itr_queryi(idx, tid, beg, end):
    return _sam_itr_queryi(idx, tid, beg, end)


def _bam_itr_querys(idx, hdr, region):
    return _sam_itr_querys(idx, hdr, region)


def _bam_itr_next(htsfp, itr, r):
    return _hts_itr_next(htsfp[0].fp.bgzf, itr, r, 0)


def _bam_index_load(fn):
    return _hts_idx_load(f, HTS_FMT_BAI)


def _bam_index_build(fn, min_shift):
    return _sam_index_build(fn, min_shift)


def _sam_itr_destroy(itr):
    return _hts_ir_destroy(itr)


def _sam_itr_next(htsfp, itr, r):
    return _hts_itr_next(htsfp[0].fp.bgzf, itr, r, htsfp)


def _sam_open(fn, mode):
    return _hts_open(fn, mode)


def _sam_close(fp):
    return hts_close(fp)

# ----------------------------------------------------------------------------
# Structures
# ----------------------------------------------------------------------------


class _bam_hdr_t(ctypes.Structure):
    """C structure for htslib type bam_hdr_t"""

    _fields_ = [('n_targets', ctypes.c_int32),
                ('ignore_sam_err', ctypes.c_int32),
                ('l_text', ctypes.c_uint32),
                ('target_len', ctypes.POINTER(ctypes.c_uint32)),
                ('cigar_tab', ctypes.POINTER(ctypes.c_int8)),
                ('target_name', ctypes.POINTER(ctypes.c_char_p)),
                ('text', ctypes.c_char_p),
                ('sdic', ctypes.c_void_p)]


class _bam1_core_t(ctypes.Structure):
    """C structure for htslib type bam1_core_t"""

    _fields_ = [('tid', ctypes.c_int32),
                ('pos', ctypes.c_int32),
                ('bin', ctypes.c_uint32, 16),
                ('qual', ctypes.c_uint32, 8),
                ('l_qname', ctypes.c_uint32, 8),
                ('l_qseq', ctypes.c_int32),
                ('mtid', ctypes.c_int32),
                ('mpos', ctypes.c_int32),
                ('isize', ctypes.c_int32)]


class _bam1_t(ctypes.Structure):
    """C structure for htslib type bam1_t"""

    _fields_ = [('core', _bam1_core_t),
                ('l_data', ctypes.c_int),
                ('m_data', ctypes.c_int),
                ('data', ctypes.POINTER(ctypes.c_uint8)),
                ('id', ctypes.c_uint64)]


_samFile = _htsFile

# ----------------------------------------------------------------------------
# C functions and their return types
# ----------------------------------------------------------------------------

htslib = pl.load_htslib()

_bam_hdr_init = htslib.bam_hdr_init
_bam_hdr_init.restype = ctypes.POINTER(_bam_hdr_t)

_bam_hdr_read = htslib.bam_hdr_read
_bam_hdr_read.restype = ctypes.POINTER(_bam_hdr_t)

_bam_hdr_write = htslib.bam_hdr_write
_bam_hdr_write.restype = ctypes.c_int

_bam_hdr_destroy = htslib.bam_hdr_destroy
_bam_hdr_destroy.restype = None

_bam_name2id = htslib.bam_name2id
_bam_name2id.restype = ctypes.c_int

_bam_hdr_dup = htslib.bam_hdr_dup
_bam_hdr_dup.restype = ctypes.POINTER(_bam_hdr_t)

_bam_init1 = htslib.bam_init1
_bam_init1.restype = ctypes.POINTER(_bam1_t)

_bam_destroy1 = htslib.bam_destroy1
_bam_destroy1.restype = None

_bam_read1 = htslib.bam_read1
_bam_read1.restype = ctypes.c_int

_bam_write1 = htslib.bam_write1
_bam_write1.restype = ctypes.c_int

_bam_copy1 = htslib.bam_copy1
_bam_copy1.restype = ctypes.POINTER(_bam1_t)

_bam_dup1 = htslib.bam_dup1
_bam_dup1.restype = ctypes.POINTER(_bam1_t)

_bam_endpos = htslib.bam_endpos
_bam_endpos.restype = ctypes.c_int32

_bam_str2flag = htslib.bam_str2flag
_bam_str2flag.restype = ctypes.c_int

_bam_flag2str = htslib.bam_flag2str
_bam_flag2str = ctypes.c_void_p  # actually a string that must be freed

_sam_index_load = htslib.sam_index_load
_sam_index_load.restype = ctypes.POINTER(_hts_idx_t)

_sam_index_load2 = htslib.sam_index_load2
_sam_index_load2.restype = ctypes.POINTER(_hts_idx_t)

_sam_index_build = htslib.sam_index_build
_sam_index_build.restype = ctypes.c_int

_sam_index_build2 = htslib.sam_index_build2
_sam_index_build2.restype = ctypes.c_int

_sam_itr_queryi = htslib.sam_itr_queryi
_sam_itr_queryi.restype = ctypes.POINTER(_hts_itr_t)

_sam_itr_querys = htslib.sam_itr_querys
_sam_itr_querys.restype = ctypes.POINTER(_hts_itr_t)

_sam_hdr_parse = htslib.sam_hdr_parse
_sam_hdr_parse.restype = ctypes.POINTER(_bam_hdr_t)

_sam_hdr_read = htslib.sam_hdr_read
_sam_hdr_read.restype = ctypes.POINTER(_bam_hdr_t)

_sam_hdr_write = htslib.sam_hdr_write
_sam_hdr_write.restype = ctypes.c_int

_sam_parse1 = htslib.sam_parse1
_sam_parse1.restype = ctypes.c_int

_sam_format1 = htslib.sam_format1
_sam_format1.restype = ctypes.c_int

_sam_read1 = htslib.sam_read1
_sam_read1.restype = ctypes.c_int

_sam_write1 = htslib.sam_write1
_sam_write1.restype = ctypes.c_int

_bam_aux_get = htslib.bam_aux_get
_bam_aux_get.restype = ctypes.c_void_p  # actually a string, user frees

_bam_aux2i = htslib.bam_aux2i
_bam_aux2i.restype = ctypes.c_int32

_bam_aux2f = htslib.bam_aux2f
_bam_aux2f.restype = ctypes.c_double

_bam_aux2A = htslib.bam_aux2A
_bam_aux2A.restype = ctypes.c_char

_bam_aux2Z = htslib.bam_aux2Z
_bam_aux2Z.restype = ctypes.c_void_p  # actually a string, freed by user

_bam_aux_append = htslib.bam_aux_append
_bam_aux_append.restype = None

_bam_aux_del = htslib.bam_aux_del
_bam_aux_del.restype = ctypes.c_int
