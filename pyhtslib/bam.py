#!/usr/bin/env python3
"""Access to SAM and BAM files through htslib."""

import ctypes
import logging
import os
import os.path
import sys

from pyhtslib.hts_internal import *  # NOQA
from pyhtslib.tabix_internal import *  # NOQA

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'


