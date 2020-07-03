#!/bin/bash

RVAL=0
set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color



../src/bam2prof  orig.bam

../src/bam2prof  3p.bam
../src/bam2prof  5p3p.bam
../src/bam2prof  5p.bam
../src/bam2prof  mutmiddle.bam


../src/bam2prof  lowqual.bam
../src/bam2prof -minq 10 lowqual.bam
