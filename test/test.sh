#!/bin/bash

RVAL=0
set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color


echo -n "Testing no mutation:"
../src/bam2prof  orig.bam  > orig.out

if diff orig.out orig.out_ > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi


echo -n "Testing mutation at the 5p end:"

../src/bam2prof  5p.bam   > 5p.out

if diff 5p.out 5p.out_ > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi

echo -n "Testing mutation at the 3p end:"

../src/bam2prof  3p.bam  > 3p.out

if diff 3p.out 3p.out_ > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi

echo -n "Testing mutation at the 5p and 3p end:"

../src/bam2prof  5p3p.bam  > 5p3p.out


if diff 5p3p.out 5p3p.out_ > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi



echo -n "Testing mutation in the middle:"

../src/bam2prof  -length 40  mutmiddle.bam >   mutmiddle.out

if diff mutmiddle.out mutmiddle.out_ > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi

echo -n "Testing base quality filter test#1:"

../src/bam2prof          lowqual.bam >  lowqualq0.out

if diff lowqualq0.out lowqualq0.out_ > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi

echo -n "Testing base quality filter test#2:"

../src/bam2prof -minq 10 lowqual.bam >  lowqualq10.out

if diff lowqualq10.out lowqualq10.out_ > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi

echo -n "Testing base quality filter test#3:"

../src/bam2prof -minq 11 lowqual.bam >  lowqualq11.out

if diff lowqualq11.out lowqualq11.out_ > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi
