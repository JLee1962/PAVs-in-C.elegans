#!/bin/bash

mkdir /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling

# After checking dot plot using nucmer and gnuplot, I picked reads for gap filling.
# When there were two or more reads for one gap region, I picked the longest one.
# Next, I checked aligned region between reference and reads using show-coords utility of MUMmer 3.

# 01. Remove the region which have accordance with gap filling region and split reference
# I saved information for spliting reference in /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/ref_for_gap_filling.txt

mkdir /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split
cat /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/ref_for_gap_filling.txt | while read rname chr rs re rregion ; do
    samtools faidx /mnt/3main/lotusbo/cb4856.genome.assembly/raw_data/cb4856.scaffold.gap.mt.fasta $rregion > /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/$rname.fasta
    sed -i '1d' /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/$rname.fasta
done

# 02. Extract gap filling sequence from the nanopore reads

cat /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/reads_for_gap_filling.txt | while read qname refname readname qs qe qregion ; do
    samtools faidx /mnt/3main/lotusbo/cb4856.genome.assembly/over20kb/$refname/$readname.fasta $qregion > /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/$qname.fasta
    sed -i '1d' /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/$qname.fasta
done

# 03. Make the gap-filled CB4856 genome

echo ">I" > /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/chrI.tmp.fasta
echo ">II" > /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/chrII.tmp.fasta
echo ">III" > /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/chrIII.tmp.fasta
echo ">IV" > /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/chrIV.tmp.fasta
echo ">V" > /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/chrV.tmp.fasta
echo ">X" > /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/chrX.tmp.fasta

list=`cat /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/I_list.txt`
for i in $list; do
    cat /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/$i.fasta >> /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/chrI.tmp.fasta
done

bioawk -c fastx '{print ">"$name; print $seq}' /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/chrI.tmp.fasta > /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/chrI.fasta

list=`cat /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/II_list.txt`
for j in $list; do
    cat /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/$j.fasta >> /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/chrII.tmp.fasta
done

bioawk -c fastx '{print ">"$name; print $seq}' /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/chrII.tmp.fasta > /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/chrII.fasta

list=`cat /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/III_list.txt`
for k in $list; do
    cat /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/$k.fasta >> /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/chrIII.tmp.fasta
done

bioawk -c fastx '{print ">"$name; print $seq}' /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/chrIII.tmp.fasta > /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/chrIII.fasta

list=`cat /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/IV_list.txt`
for l in $list; do
    cat /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/$l.fasta >> /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/chrIV.tmp.fasta
done

bioawk -c fastx '{print ">"$name; print $seq}' /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/chrIV.tmp.fasta > /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/chrIV.fasta

list=`cat /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/V_list.txt`
for m in $list; do
    cat /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/$m.fasta >> /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/chrV.tmp.fasta
done

bioawk -c fastx '{print ">"$name; print $seq}' /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/chrV.tmp.fasta > /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/chrV.fasta

list=`cat /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/X_list.txt`
for n in $list; do
    cat /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/$n.fasta >> /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/chrX.tmp.fasta
done

bioawk -c fastx '{print ">"$name; print $seq}' /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/chrX.tmp.fasta > /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/chrX.fasta

grep -A 233 "MtDNA" /mnt/3main/lotusbo/cb4856.genome.assembly/raw_data/cb4856.scaffold.gap.mt.fasta > /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/MtDNA.fasta

cat /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/chrI.fasta /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/chrII.fasta /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/chrIII.fasta /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/chrIV.fasta /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/MtDNA.fasta /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/chrV.fasta /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/split/chrX.fasta > /mnt/3main/lotusbo/cb4856.genome.assembly/gap_filling/cb4856.scaffold.gapfilling.mt.fasta
