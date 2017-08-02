#!/usr/bin/env bash
# given the path to a PDF file, this returns a count of science buzzwords found.

set -euo pipefail
IFS=$'\n\t'

if [ -z "$1" ] || [ -z "$2" ]
then
    echo "mandatory arguments: pdf, category"
    exit 1
fi

BUZZWORDS=$(pdftotext "$1" - | grep -E -i "novel|CRISPR|Bayes|leverage|machine learning|robust|paradigm|pipeline|modulate|automated|next-generation|throughput|dynamic|omics|niche|unique|tome |elucida" | wc -l)
echo "This document has ${BUZZWORDS} buzzwords in it!"
