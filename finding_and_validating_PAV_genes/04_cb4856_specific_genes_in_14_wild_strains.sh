#!/bin/bash

## Checking the states of PD1074-specific genes in the 14 wild strains

export workD="/mnt/3main/lotusbo/isoseq.analysis/cb4856/sqanti3/cbspecific"
export dbD="/mnt/3main/lotusbo/isoseq.analysis/resources/PRJNA692613/ncbi-genomes-2021-04-26/blast_db"

mkdir ${workD}/otherstrains
mkdir ${workD}/otherstrains/seq
mkdir ${workD}/otherstrains/blast
mkdir ${workD}/otherstrains/notblasted/

# 01. BLAST search of the CB4856-specific genes in the 14 wild strains
bioawk -c fastx '{print ">"$name; print $seq}' ${workD}/seq/sqanti3transcriptgenome.large.fasta >> ${workD}/otherstrains/seq/sqanti3transcriptgenome.large.fasta
awk '{print $1}' ${workD}/cbspecific_transcript/s.l.notblasted.txt >> ${workD}/otherstrains/seq/cb_pavs.txt

list=`cat ${workD}/otherstrains/seq/cb_pavs.txt`
for i in $list; do
    j=`echo $i | sed 's/\./\\\./g' | sed 's/$/$/g'`
   grep -A 1 $j ${workD}/otherstrains/seq/sqanti3transcriptgenome.large.fasta >> ${workD}/otherstrains/seq/cb_pavs_sqanti3transcriptgenome.large.fasta
done

ls -l /mnt/3main/lotusbo/isoseq.analysis/resources/PRJNA692613/ncbi-genomes-2021-04-26/blast_db/ | awk '{print $9}' >> ${workD}/otherstrains/seq/dblist.txt
sed -i 's/_genome//g' ${workD}/otherstrains/seq/dblist.txt

list=`cat ${workD}/otherstrains/seq/dblist.txt`
for i in $list; do
    j=`echo $i | sed 's/$/_genome/g'`
    blastn -query ${workD}/otherstrains/seq/cb_pavs_sqanti3transcriptgenome.large.fasta -db ${dbD}/${j}_db/$j \
        -outfmt  "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \
        -num_threads 20 -out ${workD}/otherstrains/blast/cb_pavs_sqanti3transcriptgenome.large.blast.$j.out
done

list=`cat ${workD}/otherstrains/seq/dblist.txt`
for i in $list; do
    j=`echo $i | sed 's/$/_genome/g'`
    mkdir ${workD}/otherstrains/notblasted/$j
    bioawk -c fastx '{print $name}' ${workD}/otherstrains/seq/cb_pavs_sqanti3transcriptgenome.large.fasta >> ${workD}/otherstrains/notblasted/$j/s.l.notblasted.txt
    listk=`awk '{print $1}' ${workD}/otherstrains/blast/cb_pavs_sqanti3transcriptgenome.large.blast.$j.out | sort -u | sed "s/\./\\\\\\\\./g"`
    for k in $listk; do
        sed -i "/$k$/d" ${workD}/otherstrains/notblasted/$j/s.l.notblasted.txt
    done
done

list=`cat ${workD}/otherstrains/seq/dblist.txt`
for i in $list; do
    j=`echo $i | sed 's/$/_genome/g'`
    blastn -query ${workD}/otherstrains/seq/cb_pavs_sqanti3transcriptgenome.large.fasta -db ${dbD}/${j}_db/$j \
        -outfmt  "0" \
        -num_threads 20 -out ${workD}/otherstrains/blast/cb_pavs_sqanti3transcriptgenome.large.blast.outfmt0.$j.out
done

# 02. Extracting the CB4856-specific genes absent in the 14 wild strains
export workD="/mnt/3main/lotusbo/isoseq.analysis/cb4856/sqanti3/cbspecific"
export dbD="/mnt/3main/lotusbo/isoseq.analysis/resources/PRJNA692613/ncbi-genomes-2021-04-26/blast_db"

list=`cat ${workD}/otherstrains/seq/dblist.txt`
for i in $list; do
    j=`echo $i | sed 's/$/_genome/g'`
    listk=`bioawk -c fastx '{print $name}' ${workD}/otherstrains/seq/cb_pavs_sqanti3transcriptgenome.large.fasta`
    for k in $listk; do
        l=`echo $k | sed "s/\./\\\\\\\\./g"`
        m=`grep -c $l ${workD}/otherstrains/notblasted/$j/s.l.notblasted.txt`
        if [ $m == 0 ] ; then
            n=`echo Y`
        else
            n=`echo N`
        fi
    echo -e "$k\t$n" >> ${workD}/otherstrains/notblasted/cb_pavs_otherstrain.txt
    done
done

# 03. Extracting the CB4856-specific genes present (>=90% identity) in the 14 wild strains
export workD="/mnt/3main/lotusbo/isoseq.analysis/cb4856/sqanti3/cbspecific"

mkdir ${workD}/otherstrains/over90
list=`cat ${workD}/otherstrains/seq/dblist.txt`
for i in $list; do
    awk '$15 >= 90 {print $1}' ${workD}/otherstrains/blast/cb_pavs_sqanti3transcriptgenome.large.blast.${i}_genome.out | awk '!x[$0]++ {print $0}' \
        >> ${workD}/otherstrains/over90/over90_cb_pavs_s.l.${i}_genome.txt
done

list=`cat ${workD}/otherstrains/seq/dblist.txt`
for i in $list; do
    j=`echo $i | sed 's/$/_genome/g'`
    listk=`bioawk -c fastx '{print $name}' ${workD}/otherstrains/seq/cb_pavs_sqanti3transcriptgenome.large.fasta`
    for k in $listk; do
        l=`echo $k | sed "s/\./\\\\\\\\./g"`
        m=`grep -c $l ${workD}/otherstrains/over90/over90_cb_pavs_s.l.${i}_genome.txt`
        if [ $m -gt 0 ] ; then
            n=`echo 90`
        else
            n=`echo N`
        fi
    echo -e "$i\t$k\t$n" >> ${workD}/otherstrains/over90/cb_pavs_over90_otherstrain.txt
    done
done

