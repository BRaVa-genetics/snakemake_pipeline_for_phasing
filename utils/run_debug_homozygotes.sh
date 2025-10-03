#!/bin/bash

# this script will compare a resulting set of hom:alt genotypes (e.g. as those we used in BRaVa)
# with the corresponding genotype in the raw VCF, to assess the extend of imputation performed by Shapeit5

# requirements: pysam

# for LSF submit as
# bsub -J assess_homz.[1-22] -o assess_homz.%I.%J -q normal -n 1 bash run_debug_homozygotes.sh

if [[ -z $1 ]]; then
    C=$LSB_JOBINDEX
else
    C=$1
fi

wd_pref="/PATH/TO/WD"
to_check=$wd_pref/tmp/gt_to_assess.$C.txt
output="$wd_pref/assess_homz.report.$C.txt"

# INPUT
rawVCF="PATH/TO/QC_BUT_UNPHASED/chr${C}_hard_filters.tidy.vep.vcf.gz"
outVCF="PATH/TO/PHASED/GNH_49k.notrios.phased_all.chr$C.bcf"
genotypes="PATH/TO/OUTPUT_FROM_call_chets/chr$C.pLOF.txt.gz"

mkdir -p $wd_pref/tmp/
# prepare a list of (sample, variant) to assess
# note: in case of "double" homozygotes, we only work with the first variant
zcat $genotypes | grep hom | grep -v '_' | awk '{split($NF, a, "|"); print $1, a[1]}' > $to_check
python3 $wd_pref/debug_homozygotes.py -i $rawVCF -o $outVCF --focals $to_check --hom alt > $output

# rm tmp.gt.chr$C

