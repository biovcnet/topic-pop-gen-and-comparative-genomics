#!/bin/sh

PYSCRIPT=$1
FASTA=$2
PREFIX=$3

#You need to have pandas and pyfaidx installed in python2 for this to work
#python2 -m pip install pandas
#python2 -m pip install pyfaidx

#Run script from iMKT to process your alignment
#Your alignment MUST have 5 seqs (4 for polymorphism, 1 outgroup which should be the last sequence)
python2 $PYSCRIPT --multiFasta $FASTA --daf $PREFIX.daf.tmp --div $PREFIX.div.tmp --codonTable standard

#The iMKT helper script provided outputs columns in the wrong order (I suspect they changed it and forgot to update)
awk '{print $1"\t"$3"\t"$2}' $PREFIX.daf.tmp > $PREFIX.daf
awk '{print $3"\t"$2"\t"$4"\t"$1}' $PREFIX.div.tmp > $PREFIX.div

#Cleanup tmp files
rm $PREFIX.daf.tmp
rm $PREFIX.div.tmp
