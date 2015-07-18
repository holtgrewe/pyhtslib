#!/usr/bin/env python
"""Helper code for loading the htslib dynamic library."""

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'

import ctypes
import ctypes.util
import logging
import os


def load_libc():
    """Load libc so we can free the memory later on."""
    libc = ctypes.CDLL(ctypes.util.find_library('c'))
    if not libc:
        raise Exception('Could not load libc.')
    return libc


def load_htslib():
    """Attempt to load htslib dynamic library.

    Try to get path to library file through environment variable
    ``HTSLIB_PATH``.  If this fails, attempt to find through
    ``ctypes.util.find_library``.  If this fails, raise an Exception.
    """
    # obtain path to libhts.so library (or similar)
    logging.debug('attempting to obtain htslib path from environment '
                  'HTSLIB_PATH')
    path = os.environ.get('HTSLIB_PATH')
    if not path:
        logging.debug('attempting to obtain htslib path through '
                      'find_library()')
        path = ctypes.util.find_library('hts')
    if not path:
        raise Exception('Could not load obtain path to htslib, tried '
                        'environment for HTSLIB_PATH and find_library().')
    # actually load the library
    logging.debug('loading libhts from %s', path)
    htslib = ctypes.CDLL(path)
    if not htslib:
        raise Exception('Could not load htslib, tried from {}'.format(path))
    return htslib
