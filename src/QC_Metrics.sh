#!/bin/bash

ARG1=$1

UTR=$(samtools view $ARG1 | column -t | grep -c "XF:Z:UTR")
CODING=$(samtools view $ARG1 |column -t | grep -c "XF:Z:CODING")
INTRONIC=$(samtools view $ARG1 |column -t | grep -c "XF:Z:INTRONIC")
INTERGENIC=$(samtools view $ARG1 |column -t | grep -c "XF:Z:INTERGENIC")
asUTR=$(samtools view $ARG1 |grep -v GE:Z |grep -c "XF:Z:UTR")
asCODING=$(samtools view $ARG1 |grep -v GE:Z |grep -c "XF:Z:CODING")

echo number of UTR
echo $UTR
echo number of CODING
echo $CODING
echo number of INTRONIC
echo $INTRONIC
echo number of INTERGENIC
echo $INTERGENIC
echo number of as.UTR
echo $asUTR
echo number of as.CODING
echo $asCODING

echo ----------------------------------------------------
samtools stats $ARG1 | head -30



