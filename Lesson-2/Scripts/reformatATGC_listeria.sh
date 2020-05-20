#!/bin/bash

# JLW - 2020
# This script unpacks the Listeria monocytogenes ATGC cluster into individual gene alignments

cd Data/

OUTGROUP="272626"

../Scripts/fasta2tbl.sh listeria_monocytogenes.fa > listeria_monocytogenes.tbl

awk '{print $1}' listeria_monocytogenes.tbl | awk -F"|" '{print $4}' | sort | uniq | grep -v "$OUTGROUP" > tax_ids.txt

awk -F"|" '{print $1" "$2}' listeria_monocytogenes.tbl | sort | uniq -c | grep "12 synCogId" | awk '{print $3}' > universal_cogs.txt
awk -F"|" 'NR==FNR {id[$1]; next} $2 in id' universal_cogs.txt listeria_monocytogenes.tbl > listeria_monocytogenes_uc.tbl

grep -v "taxonId|272626" listeria_monocytogenes_uc.tbl > listeria_monocytogenes_uc_notout.tbl
grep "taxonId|272626" listeria_monocytogenes_uc.tbl > listeria_monocytogenes_uc_out.tbl

mkdir COGs
while read COG; do
	grep "synCogId|$COG" listeria_monocytogenes_uc_notout.tbl > COGs/synCOG$COG.tbl
        grep "synCogId|$COG" listeria_monocytogenes_uc_notout.tbl >> COGs/synCOG$COG.tbl
        ../Scripts/tbl2fasta.sh COGs/synCOG$COG.tbl > COGs/synCOG$COG.fasta
        rm COGs/synCOG$COG.tbl
done < universal_cogs.txt



