#!/bin/bash
mkdir /mnt/3main/lotusbo/isoseq.analysis/cb4856/sqanti3/cbspecific
mkdir /mnt/3main/lotusbo/isoseq.analysis/cb4856/sqanti3/cbspecific/seq
export comD="/mnt/3main/lotusbo/isoseq.analysis/cb4856"
export workD="/mnt/3main/lotusbo/isoseq.analysis/cb4856/sqanti3/cbspecific"

## Finding pd1074-specific genes
# 01. Extracting the genomic sequence of transcripts of SQANTI3 analysis
awk '{ if ($3=="transcript") {print $0}}' ${comD}/sqanti3/sqanti3_qc/cb4856.collapsed.filtered.rep_corrected.gtf.cds.gff \
    | awk '{ split ($12, arr, "\""); print $1"\t"$4"\t"$5"\t"arr[2] }' | while read chr start stop name
do
    echo \>$name >> ${workD}/seq/sqanti3transcriptgenome.fasta
    samtools faidx /mnt/3main/lotusbo/isoseq.analysis/resources/cb4856/cb4856.scaffold.gapfilling.mt.wrapped.fasta \
    $chr:$start-$stop >> ${workD}/seq/sqanti3transcriptgenome.fasta
done
sed -i '/chr/d' ${workD}/seq/sqanti3transcriptgenome.fasta
bioawk -c fastx 'length($seq) <= 50 {print ">"$name"\n"$seq}' ${workD}/seq/sqanti3transcriptgenome.fasta >> ${workD}/seq/sqanti3transcriptgenome.small.fasta 
bioawk -c fastx 'length($seq) > 50 {print ">"$name"\n"$seq}' ${workD}/seq/sqanti3transcriptgenome.fasta >> ${workD}/seq/sqanti3transcriptgenome.large.fasta


# 02. Making a list of genes which have not similarity with the PD1074 genome
mkdir ${workD}/blast \
blastn -task blastn-short -query ${workD}/seq/sqanti3transcriptgenome.small.fasta -db /mnt/3main/lotusbo/isoseq.analysis/resources/blast_db/pd1074_genome/pd1074_genome \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \
    -num_threads 20 -out ${workD}/blast/sqanti3transcriptgenome.small.blast.pd1074_genome.out
blastn -query ${workD}/seq/sqanti3transcriptgenome.large.fasta -db /mnt/3main/lotusbo/isoseq.analysis/resources/blast_db/pd1074_genome/pd1074_genome \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \
    -num_threads 20 -out ${workD}/blast/sqanti3transcriptgenome.large.blast.pd1074_genome.out

mkdir ${workD}/notblasted
bioawk -c fastx '{print $name}' ${workD}/seq/sqanti3transcriptgenome.small.fasta >> ${workD}/notblasted/s.s.notblasted.txt
bioawk -c fastx '{print $name}' ${workD}/seq/sqanti3transcriptgenome.large.fasta >> ${workD}/notblasted/s.l.notblasted.txt

list=`awk '{print $1}' ${workD}/blast/sqanti3transcriptgenome.large.blast.pd1074_genome.out | sort -u | sed "s/\./\\\\\\\\./g"`
for i in $list; do
    sed -i "/$i$/d" ${workD}/notblasted/s.l.notblasted.txt
done

# 03. Extracting protein sequences of the cb4856-specific genes for BLASTp and Batch CD search

list=`cat ${workD}/cbspecific_transcript/s.l.notblasted.txt`
for i in $list; do
    j=`echo $i | sed 's/\./\\\./g'`
    echo ">$i" >> ${workD}/cbspecific_transcript/s.notblasted.faa
    grep -A 1 $j ${comD}/sqanti3//sqanti3_qc/cb4856.collapsed.filtered.rep_corrected.faa |sed '1d' >> ${workD}/cbspecific_transcript/s.notblasted.faa
done

