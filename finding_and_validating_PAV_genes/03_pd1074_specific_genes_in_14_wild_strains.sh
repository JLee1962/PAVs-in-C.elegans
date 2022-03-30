#!/bin/bash

## Checking the states of PD1074-specific genes in the 14 wild strains

export workD="/mnt/3main/lotusbo/isoseq.analysis/pd1074/sqanti3/pdspecific_2"
export dbD="/mnt/3main/lotusbo/isoseq.analysis/resources/PRJNA692613/ncbi-genomes-2021-04-26/blast_db"

mkdir ${workD}/otherstrains
mkdir ${workD}/otherstrains/seq
mkdir ${workD}/otherstrains/blast
mkdir ${workD}/otherstrains/notblasted/

# 01. BLAST search of the PD1074-specific genes in the 14 wild strains
bioawk -c fastx '{print ">"$name; print $seq}' ${workD}/seq/unmappedtranscriptgenome.large.fasta >> ${workD}/otherstrains/seq/unmappedtranscriptgenome.large.fasta
awk '{print $1}' ${workD}/pdspecific_2_transcript/pdspecific_2_transcript.notblasted.txt >> ${workD}/otherstrains/seq/pd_pavs_unmapped.large.txt
 modify manually ${workD}/otherstrains/seq/pd_pavs_unmapped.large.txt
list=`cat ${workD}/otherstrains/seq/pd_pavs_unmapped.large.txt`
for i in $list; do
    j=`echo $i | sed 's/\./\\\./g' | sed 's/$/$/g'`
   grep -A 1 $j ${workD}/otherstrains/seq/unmappedtranscriptgenome.large.fasta >> ${workD}/otherstrains/seq/pd_pavs_unmappedtranscriptgenome.large.fasta
done

ls -l /mnt/3main/lotusbo/isoseq.analysis/resources/PRJNA692613/ncbi-genomes-2021-04-26/blast_db/ | awk '{print $9}' >> ${workD}/otherstrains/seq/dblist.txt
sed -i 's/_genome//g' ${workD}/otherstrains/seq/dblist.txt

list=`cat ${workD}/otherstrains/seq/dblist.txt`
for i in $list; do
    j=`echo $i | sed 's/$/_genome/g'`
    blastn -query ${workD}/otherstrains/seq/pd_pavs_sqanti3transcriptgenome.large.fasta -db ${dbD}/${j}_db/$j \
        -outfmt  "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \
        -num_threads 20 -out ${workD}/otherstrains/blast/pd_pavs_sqanti3transcriptgenome.large.blast.$j.out
    blastn -task blastn-short -query ${workD}/otherstrains/seq/pd_pavs_unmappedtranscriptgenome.small.fasta -db ${dbD}/${j}_db/$j \
        -outfmt  "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \
        -num_threads 20 -out ${workD}/otherstrains/blast/pd_pavs_unmappedtranscriptgenome.small.blast.$j.out
    blastn -query ${workD}/otherstrains/seq/pd_pavs_unmappedtranscriptgenome.large.fasta -db ${dbD}/${j}_db/$j \
        -outfmt  "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \
        -num_threads 20 -out ${workD}/otherstrains/blast/pd_pavs_unmappedtranscriptgenome.large.blast.$j.out
done

list=`cat ${workD}/otherstrains/seq/dblist.txt`
for i in $list; do
    j=`echo $i | sed 's/$/_genome/g'`
    mkdir ${workD}/otherstrains/notblasted/$j
    bioawk -c fastx '{print $name}' ${workD}/otherstrains/seq/pd_pavs_sqanti3transcriptgenome.large.fasta >> ${workD}/otherstrains/notblasted/$j/s.l.notblasted.txt
    bioawk -c fastx '{print $name}' ${workD}/otherstrains/seq/pd_pavs_unmappedtranscriptgenome.small.fasta >> ${workD}/otherstrains/notblasted/$j/u.s.notblasted.txt
    bioawk -c fastx '{print $name}' ${workD}/otherstrains/seq/pd_pavs_unmappedtranscriptgenome.large.fasta >> ${workD}/otherstrains/notblasted/$j/u.l.notblasted.txt
    listk=`awk '{print $1}' ${workD}/otherstrains/blast/pd_pavs_sqanti3transcriptgenome.large.blast.$j.out | sort -u | sed "s/\./\\\\\\\\./g"`
    for k in $listk; do
        sed -i "/$k$/d" ${workD}/otherstrains/notblasted/$j/s.l.notblasted.txt
    done
    listl=`awk '{print $1}' ${workD}/otherstrains/blast/pd_pavs_unmappedtranscriptgenome.small.blast.$j.out | sort -u | sed "s/\./\\\\\\\\./g"`
    for l in $listl; do
        sed -i "/$l$/d" ${workD}/otherstrains/notblasted/$j/u.s.notblasted.txt
    done
    listm=`awk '{print $1}' ${workD}/otherstrains/blast/pd_pavs_unmappedtranscriptgenome.large.blast.$j.out | sort -u | sed "s/\./\\\\\\\\./g"`
    for m in $listm; do
        sed -i "/$m$/d" ${workD}/otherstrains/notblasted/$j/u.l.notblasted.txt
    done
done

# 02. Extracting the PD1074-specific genes absent in the 14 wild strains
export workD="/mnt/3main/lotusbo/isoseq.analysis/pd1074/sqanti3/pdspecific_2"
export dbD="/mnt/3main/lotusbo/isoseq.analysis/resources/PRJNA692613/ncbi-genomes-2021-04-26/blast_db"

list=`cat ${workD}/otherstrains/seq/dblist.txt`
for i in $list; do
    j=`echo $i | sed 's/$/_genome/g'`
    listk=`bioawk -c fastx '{print $name}' ${workD}/otherstrains/seq/pd_pavs_sqanti3transcriptgenome.large.fasta`
    for k in $listk; do
        l=`echo $k | sed "s/\./\\\\\\\\./g"`
        m=`grep -c $l ${workD}/otherstrains/notblasted/$j/s.l.notblasted.txt`
        if [ $m == 0 ] ; then
            n=`echo Y`
        else
            n=`echo N`
        fi
    echo -e "$i\t$k\t$n" >> ${workD}/otherstrains/notblasted/pd_pavs_otherstrain.txt
    done
done

list=`cat ${workD}/otherstrains/seq/dblist.txt`
for i in $list; do
    j=`echo $i | sed 's/$/_genome/g'`
    listk=`bioawk -c fastx '{print $name}' ${workD}/otherstrains/seq/pd_pavs_unmappedtranscriptgenome.large.fasta`
    for k in $listk; do
        l=`echo $k | sed "s/\./\\\\\\\\./g"`
        m=`grep -c $l ${workD}/otherstrains/notblasted/$j/u.l.notblasted.txt`
        if [ $m == 0 ] ; then
            n=`echo Y`
        else
            n=`echo N`
        fi
    echo -e "$i\t$k\t$n" >> ${workD}/otherstrains/notblasted/pd_pavs_otherstrain.txt
    done
done

list=`cat ${workD}/otherstrains/seq/dblist.txt`
for i in $list; do
    j=`echo $i | sed 's/$/_genome/g'`
    listk=`bioawk -c fastx '{print $name}' ${workD}/otherstrains/seq/pd_pavs_unmappedtranscriptgenome.small.fasta`
    for k in $listk; do
        l=`echo $k | sed "s/\./\\\\\\\\./g"`
        m=`grep -c $l ${workD}/otherstrains/notblasted/$j/u.s.notblasted.txt`
        if [ $m == 0 ] ; then
            n=`echo N`
        else
            n=`echo Y`
        fi
    echo -e "$i\t$k\t$n" >> ${workD}/otherstrains/notblasted/pd_pavs_otherstrain.txt
    done
done

# 03. Extracting the PD1074-specific genes present (>=90% identity) in the 14 wild strains
export workD="/mnt/3main/lotusbo/isoseq.analysis/pd1074/sqanti3/pdspecific_2"

mkdir ${workD}/otherstrains/over90
list=`cat ${workD}/otherstrains/seq/dblist.txt`
for i in $list; do
    awk '$15 >= 90 {print $1}' ${workD}/otherstrains/blast/pd_pavs_sqanti3transcriptgenome.large.blast.${i}_genome.out | awk '!x[$0]++ {print $0}' \
        >> ${workD}/otherstrains/over90/over90_pd_pavs_s.l.${i}_genome.txt
    awk '$15 >= 90 {print $1}' ${workD}/otherstrains/blast/pd_pavs_unmappedtranscriptgenome.large.blast.${i}_genome.out | awk '!x[$0]++ {print $0}' \
        >> ${workD}/otherstrains/over90/over90_pd_pavs_u.l.${i}_genome.txt
done

list=`cat ${workD}/otherstrains/seq/dblist.txt`
for i in $list; do
    j=`echo $i | sed 's/$/_genome/g'`
    listk=`bioawk -c fastx '{print $name}' ${workD}/otherstrains/seq/pd_pavs_sqanti3transcriptgenome.large.fasta`
    for k in $listk; do
        l=`echo $k | sed "s/\./\\\\\\\\./g"`
        m=`grep -c $l ${workD}/otherstrains/over90/over90_pd_pavs_s.l.${i}_genome.txt`
        if [ $m -gt 0 ] ; then
            n=`echo 90`
        else
            n=`echo N`
        fi
    echo -e "$i\t$k\t$n" >> ${workD}/otherstrains/over90/pd_pavs_over90_otherstrain.txt
    done
done
list=`cat ${workD}/otherstrains/seq/dblist.txt`
for i in $list; do
    j=`echo $i | sed 's/$/_genome/g'`
    listk=`bioawk -c fastx '{print $name}' ${workD}/otherstrains/seq/pd_pavs_unmappedtranscriptgenome.large.fasta`
    for k in $listk; do
        l=`echo $k | sed "s/\./\\\\\\\\./g"`
        m=`grep -c $l ${workD}/otherstrains/over90/over90_pd_pavs_u.l.${i}_genome.txt`
        if [ $m -gt 0 ] ; then
            n=`echo 90`
        else
            n=`echo N`
        fi
    echo -e "$i\t$k\t$n" >> ${workD}/otherstrains/over90/pd_pavs_over90_otherstrain.txt
    done
done

