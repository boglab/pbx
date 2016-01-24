#!/usr/bin/env bash

source /opt/smrtanalysis/current/etc/setup.sh

cd ${1}
for f in *.bas.h5; do fprefix=`echo "$f" | rev | cut -c 8- | rev`; bash5tools.py ${f} --outFilePref ${fprefix} --outType fastq --minReadScore 0.75 --minLength 4000; done;
cat *.fastq > ${2}
rm *.fastq