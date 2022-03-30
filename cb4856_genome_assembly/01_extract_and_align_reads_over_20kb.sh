#!/bin/bash

# 01. Extract reads over 20kb
mkdir over20kb
bioawk -cfastx 'length($seq) > 20000 {print "@"$name"\n"$seq"\n+\n"$qual}' \
    /mnt/3main/lotusbo/cb4856.genome.assembly/raw_data/cb4856.ont.raw_read.fastq \
    > /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/cb4856.over20kb.ont.fastq

# 02. Align reads over 20kb on the cb4856 genome
minimap2 -ax map-ont -t 23 /mnt/3main/lotusbo/cb4856.genome.assembly/raw_data/cb4856.scaffold.gap.mt.fasta \
    /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/cb4856.over20kb.ont.fastq \
    | samtools sort -@ 20 -O sam -o /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/cb4856.over20kb.ont.minimap2.sam

# 03. Sort and re-align primary reads with >2 mapq
samtools view -f 0x10 -q 2 /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/cb4856.over20kb.ont.minimap2.sam \
    | awk '{print $1}' >> /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/primary.list.txt
samtools view -q 2 /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/cb4856.over20kb.ont.minimap2.sam \
    | awk '$2 == 0 {print $1}' >> /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/primary.list.txt
sort -u /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/primary.list.txt \
    > /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/primary.u.list.txt

list=`cat /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/primary.u.list.txt`
for i in $list ; do
    grep -A 3 $i /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/cb4856.over20kb.ont.fastq \
    >> /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/primary.cb4856.over20kb.ont.fastq
done

minimap2 -ax map-ont -t 23 /mnt/3main/lotusbo/cb4856.genome.assembly/raw_data/cb4856.scaffold.gap.mt.fasta \
    /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/primary.cb4856.over20kb.ont.fastq \
    | samtools sort -@ 20 -O bam -o /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/primary.cb4856.over20kb.ont.bam
samtools index -b -@ 20 /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/primary.cb4856.over20kb.ont.bam

# I calculated the gap location of cb4856.scaffold.gap.mt.fasta and saved it to /mnt/3main/lotusbo/cb4856.genome.assembly/raw_data/gap_location.txt
# I could get gap location by using this command: seqtk gap.
# Using IGV, I checked the reads which located near gap region, and saved it to /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/read_for_gaps/gap_read.txt
# I modified gap_read.txt for next step, and saved it to /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/read_for_gaps/gap_read_list.txt

# 04. Extract the genome sequence of reference region
mkdir /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/reference
cat /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/read_for_gaps/gap_read_list.txt \
    | awk '{print $3"\t"$4}' \
    | sort -u > /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/read_for_gaps/rregion.txt

cat /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/read_for_gaps/rregion.txt | while read rname rregion ; do
    samtools faidx /mnt/3main/lotusbo/cb4856.genome.assembly/raw_data/cb4856.scaffold.gap.mt.fasta $rregion \
    > /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/reference/$rname.fasta
done

# 05. Extract the query reads (the reads of # 03) and align to the reference with nucmer
list=`cat /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/read_for_gaps/gap_read_list.txt | awk '{print $3}' | sort -u`
for i in $list; do
    mkdir /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/$i
done

cat /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/read_for_gaps/gap_read_list.txt | while read -r qname ori rname rregion
do
    mkdir /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/$rname
    if [ $ori == "+" ]; then
        grep -A 3 $qname /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/cb4856.over20kb.ont.fastq \
            | bioawk -c fastx '{print">"$name; print $seq}' \
            > /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/$rname/$qname.fasta
        /mnt/main/tools/mummer-4.0.0beta2/nucmer -t 10 -p /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/$rname/$qname \
            /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/reference/$rname.fasta \
            /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/$rname/$qname.fasta
        /mnt/main/tools/mummer-4.0.0beta2/mummerplot -p /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/$rname/$qname \
            --png /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/$rname/$qname.delta
        sed 's/png tiny size 800,800/pdf/' /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/$rname/$qname.gp \
            | sed 's/png/pdf/' | sed 's/line 1  lt 1 lw 3 pt 6 ps 1/line 1  lt rgb "#FF0000" lw 2 pt 2 ps 0.1/' \
            | sed 's/line 2  lt 3 lw 3 pt 6 ps 1/line 2  lt rgb "#0066CC" lw 2 pt 2 ps 0.1/' \
            > /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/$rname/$qname.pdf.gp
    elif [ $ori == "-" ]; then
        mkdir /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/$rname
        grep -A 3 $qname /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/cb4856.over20kb.ont.fastq \
            | bioawk -c fastx '{print">"$name; print revcomp($seq)}' \
            > /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/$rname/$qname.fasta
        /mnt/main/tools/mummer-4.0.0beta2/nucmer -t 10 -p /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/$rname/$qname \
            /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/reference/$rname.fasta \
            /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/$rname/$qname.fasta
        /mnt/main/tools/mummer-4.0.0beta2/mummerplot -p /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/$rname/$qname \
            --png /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/$rname/$qname.delta
        sed 's/png tiny size 800,800/pdf/' /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/$rname/$qname.gp \
            | sed 's/png/pdf/' | sed 's/line 1  lt 1 lw 3 pt 6 ps 1/line 1  lt rgb "#FF0000" lw 2 pt 2 ps 0.1/' \
            | sed 's/line 2  lt 3 lw 3 pt 6 ps 1/line 2  lt rgb "#0066CC" lw 2 pt 2 ps 0.1/' \
            > /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/$rname/$qname.pdf.gp
    else
        echo "none"
    fi
done

cat /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/read_for_gaps/gap_read_list.txt | while read -r qname ori rname rregion
do
    gnuplot /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/$rname/$qname.pdf.gp
done
