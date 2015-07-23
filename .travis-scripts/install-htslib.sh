#!/bin/sh
set -ex

test -f "${HOME}/htslib-1.2.1/libhts.so" && exit 0

wget https://github.com/samtools/htslib/releases/download/1.2.1/htslib-1.2.1.tar.bz2 \
    -O /tmp/htslib-1.2.1.tar.bz2
tar -xjf /tmp/htslib-1.2.1.tar.bz2

cat <<EOF >>htslib-1.2.1/vcf.c
int _bcf_hdr_name2id(const bcf_hdr_t *hdr, const char *id) { return bcf_hdr_id2int(hdr, BCF_DT_CTG, id); }
EOF

cd htslib-1.2.1
perl -p -i -e "s/^prefix.*/prefix=\/home\/travis\/htslib-1.2.1/g" Makefile
make -j 2
make install
