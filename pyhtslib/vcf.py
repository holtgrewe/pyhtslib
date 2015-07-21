#!/usr/bin/env python3
"""Access to VCF and BCF files"""

import ctypes
import logging
import os
import os.path

from pyhtslib.hts_internal import *  # NOQA
from pyhtslib.vcf_internal import *  # NOQA

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'
