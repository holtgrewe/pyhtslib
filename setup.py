#!/usr/bin/env python3
"""Setup script for the package pyhtslib"""

import os.path
import setuptools

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'


def read_version():
    """Read version from VERSION file relative to this file"""
    with open(os.path.join(os.path.dirname(__file__), 'VERSION')) as f:
        return f.read().strip()


# actually define the setuptools package
setuptools.setup(
    # name, version, and description of the package
    name='pyhtslib',
    description='Python wrapper for htslib',
    version=read_version(),
    url='https://bihealth.github.io/pyhtslib',

    # author information
    author='Manuel Holtgrewe',
    author_email='manuel.holtgrewe@bihealth.de',

    # contained packages
    packages=setuptools.find_packages(
        exclude=['_cffi_build', '_cffi_build.*'],
    ),
    ext_package='pyhtslib',

    # dependencies
    setup_requires=['cffi>=1.1'],
    cffi_modules=['pyhtslib/build_pyhtslib.py:ffi'],
    install_requires=['cffi>=1.1'],
)
