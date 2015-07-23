# pyhtslib -- HTS file access from Python

The Python library pyhtslib is a wrapper around [htslib](http://www.htslib.org/), a C library for accessing files used for high-throughput sequencing (HTS) file formats.
The aim of pyhtslib is to provide the I/O functionality from htslib through an easy-to-use and well-tested Python interface.

## Status

Read access to the following formats works:

- SAM/BAM -- sequential and indexed (BAM through `.bai` files, SAM through tabix indices)
- VCF/BCF -- sequential and indexed (BCF through `.csi` files, VCF through tabix indices)
- tabix -- reading of arbitrary TSV files
- FAI -- indexed FASTA

What is missing:

- writing of files [v0.6]
- CRAM support [v0.5]
- sequential FASTA and FASTQ through the `kseq.h` library [v0.4]

Other things on the roadmap for a v1.0:

- support for Python 2 (a lot of Python Bioinformatics tools have not moved to Python 3 yet) [v0.2]
- comprehensive documentation [v0.2]
- switching to Cython instead of using `ctypes` [v0.3]

The plan is then to support the more advanced features of htslib:

- support for `bcf\_synced\_reader` [v1.2]
- more lazy loading of VCF properties [v1.1]
- pileup functionality [v1.3]
- reading of multiple sorted BAM files at the same time [v1.4]


## Contributors

- Manuel Holtgrewe, Berlin Institute of Health/Charite University Medicine Berlin
