#!/bin/bash

mkdir /mnt/3main/lotusbo/isoseq.analysis/pd1074/sqanti3/pdspecific_2
mkdir /mnt/3main/lotusbo/isoseq.analysis/pd1074/sqanti3/pdspecific_2/seq
export comD="/mnt/3main/lotusbo/isoseq.analysis/pd1074"
export workD="/mnt/3main/lotusbo/isoseq.analysis/pd1074/sqanti3/pdspecific_2"

## Finding pd1074-specific genes
# 01. Extracting genomic sequence of transcripts of SQANTI3 analysis
awk '{ if ($3=="transcript") {print $0}}' ${comD}/sqanti3/sqanti3_qc/pd1074.collapsed.filtered.rep_corrected.gtf.cds.gff \
    | awk '{ split ($12, arr, "\""); print $1"\t"$4"\t"$5"\t"arr[2] }' | while read chr start stop name
do
    echo \>$name >> ${workD}/seq/sqanti3transcriptgenome.fasta
    samtools faidx /mnt/3main/lotusbo/isoseq.analysis/resources/pd1074/c_elegans.PRJEB28388.WS274.genomic.rename.fa \
    $chr:$start-$stop >> ${workD}/seq/sqanti3transcriptgenome.fasta
done
sed -i '/chr/d' ${workD}/seq/sqanti3transcriptgenome.fasta
bioawk -c fastx 'length($seq) <= 50 {print ">"$name"\n"$seq}' ${workD}/seq/sqanti3transcriptgenome.fasta >> ${workD}/seq/sqanti3transcriptgenome.small.fasta
bioawk -c fastx 'length($seq) > 50 {print ">"$name"\n"$seq}' ${workD}/seq/sqanti3transcriptgenome.fasta >> ${workD}/seq/sqanti3transcriptgenome.large.fasta


# 02. Extracting genomic sequence of PD1074 genes which is not transferred to CB4856 genome
sed '/#/d' /mnt/3main/lotusbo/isoseq.analysis/cb4856/isoseq3.one/liftOver/pd1074.reorganized_annotation.ptcunmapped.gff3 | \
    awk '$3=="mRNA"' >> ${workD}/seq/unmapped_transcript.gff3

awk '$3=="exon" {split ($9, arr, "="); print arr[3] }' /mnt/3main/lotusbo/isoseq.analysis/cb4856/isoseq3.one/liftOver/pd1074.reorganized_annotation.ptcliftOver.gff3 \
    >> ${workD}/seq/deletion_list_1

awk '$3=="gene" {split ($9, arr, "=|;"); print arr[2] }' /mnt/3main/lotusbo/isoseq.analysis/cb4856/isoseq3.one/liftOver/pd1074.reorganized_annotation.ptcliftOver.gff3 \
    >> ${workD}/seq/deletion_list_2

list=`sort -u ${workD}/seq/deletion_list_1`
for i in $list; do
    j="ID=$i;"
    awk -v j="$j" '$9!~j' ${workD}/seq/unmapped_transcript.gff3 >> ${workD}/seq/tmp.unmapped_transcript.gff3
    rm ${workD}/seq/unmapped_transcript.gff3
    mv ${workD}/seq/tmp.unmapped_transcript.gff3 ${workD}/seq/unmapped_transcript.gff3
done

list=`sort -u ${workD}/seq/deletion_list_2`
for i in $list; do
    j="Parent=$i$"
    awk -v j="$j" '$9!~j' ${workD}/seq/unmapped_transcript.gff3 >> ${workD}/seq/tmp.unmapped_transcript.gff3
    rm ${workD}/seq/unmapped_transcript.gff3
    mv ${workD}/seq/tmp.unmapped_transcript.gff3 ${workD}/seq/unmapped_transcript.gff3
done

awk '{ split ($9, arr, "=|;"); print $1"\t"$4"\t"$5"\t"arr[2] }' ${workD}/seq/unmapped_transcript.gff3 | while read chr start stop name
do
    echo \>$name >> ${workD}/seq/unmappedtranscriptgenome.fasta
    samtools faidx /mnt/3main/lotusbo/isoseq.analysis/resources/pd1074/c_elegans.PRJEB28388.WS274.genomic.rename.fa \
    $chr:$start-$stop >> ${workD}/seq/unmappedtranscriptgenome.fasta
done
sed -i '/>chr/d' ${workD}/seq/unmappedtranscriptgenome.fasta

bioawk -c fastx 'length($seq) <= 50 {print ">"$name"\n"$seq}' ${workD}/seq/unmappedtranscriptgenome.fasta >> ${workD}/seq/unmappedtranscriptgenome.small.fasta
bioawk -c fastx 'length($seq) > 50 {print ">"$name"\n"$seq}' ${workD}/seq/unmappedtranscriptgenome.fasta >> ${workD}/seq/unmappedtranscriptgenome.large.fasta

# 03. Making a list of genes which have not similarity with the CB4856 genome
mkdir ${workD}/blast
blastn -task blastn-short -query ${workD}/seq/sqanti3transcriptgenome.small.fasta -db /mnt/3main/lotusbo/isoseq.analysis/resources/blast_db/cb4856_genome/cb4856_genome \
    -outfmt  "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \
    -num_threads 20 -out ${workD}/blast/sqanti3transcriptgenome.small.blast.cb4856_genome.out
blastn -query ${workD}/seq/sqanti3transcriptgenome.large.fasta -db /mnt/3main/lotusbo/isoseq.analysis/resources/blast_db/cb4856_genome/cb4856_genome \
    -outfmt  "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \
    -num_threads 20 -out ${workD}/blast/sqanti3transcriptgenome.large.blast.cb4856_genome.out
blastn -task blastn-short -query ${workD}/seq/unmappedtranscriptgenome.small.fasta -db /mnt/3main/lotusbo/isoseq.analysis/resources/blast_db/cb4856_genome/cb4856_genome \
    -outfmt  "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \
    -num_threads 20 -out ${workD}/blast/unmappedtranscriptgenome.small.blast.cb4856_genome.out
blastn -query ${workD}/seq/unmappedtranscriptgenome.large.fasta -db /mnt/3main/lotusbo/isoseq.analysis/resources/blast_db/cb4856_genome/cb4856_genome \
    -outfmt  "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \
    -num_threads 20 -out ${workD}/blast/unmappedtranscriptgenome.large.blast.cb4856_genome.out

mkdir ${workD}/notblasted
bioawk -c fastx '{print $name}' ${workD}/seq/sqanti3transcriptgenome.small.fasta >> ${workD}/notblasted/s.s.notblasted.txt
bioawk -c fastx '{print $name}' ${workD}/seq/sqanti3transcriptgenome.large.fasta >> ${workD}/notblasted/s.l.notblasted.txt
bioawk -c fastx '{print $name}' ${workD}/seq/unmappedtranscriptgenome.small.fasta >> ${workD}/notblasted/u.s.notblasted.txt
bioawk -c fastx '{print $name}' ${workD}/seq/unmappedtranscriptgenome.large.fasta >> ${workD}/notblasted/u.l.notblasted.txt

list=`awk '{print $1}' ${workD}/blast/sqanti3transcriptgenome.small.blast.cb4856_genome.out | sort -u | sed "s/\./\\\\\\\\./g"`
for i in $list; do
    sed -i "/$i$/d" ${workD}/notblasted/s.s.notblasted.txt
done
list=`awk '{print $1}' ${workD}/blast/sqanti3transcriptgenome.large.blast.cb4856_genome.out | sort -u | sed "s/\./\\\\\\\\./g"`
for i in $list; do
    sed -i "/$i$/d" ${workD}/notblasted/s.l.notblasted.txt
done
list=`awk '{print $1}' ${workD}/blast/unmappedtranscriptgenome.small.blast.cb4856_genome.out | sort -u | sed "s/\./\\\\\\\\./g"`
for i in $list; do
    sed -i "/$i$/d" ${workD}/notblasted/u.s.notblasted.txt
done
list=`awk '{print $1}' ${workD}/blast/unmappedtranscriptgenome.large.blast.cb4856_genome.out | sort -u | sed "s/\./\\\\\\\\./g"`
for i in $list; do
    sed -i "/$i$/d" ${workD}/notblasted/u.l.notblasted.txt
done


# 04. Extracting annotations of pd1074-specific genes
echo -e "name\tlocation\tlength\tnumberofexon\tcds\tgeneID" >> ${workD}/pdspecific_2_transcript/pdspecific_2_transcript.notblasted.txt
list=`cat ${workD}/notblasted/s.l.notblasted.txt`
for i in $list; do
    j=\"$i\"\;
    k=`awk '$3=="transcript"' ${comD}/sqanti3/sqanti3_qc/pd1074.collapsed.filtered.rep_corrected.gtf.cds.gff | awk -v j="$j" '$12==j {print $1":"$4"-"$5}'`
    l=`awk -v i="$i" '$1==i {print $2}' ${workD}/seq/sqanti3transcript.fasta.length`
    m=`awk '$3=="exon"' ${comD}/sqanti3/sqanti3_qc/pd1074.collapsed.filtered.rep_corrected.gtf.cds.gff | grep -c $j`
    p=`awk '$3=="CDS"' ${comD}/sqanti3/sqanti3_qc/pd1074.collapsed.filtered.rep_corrected.gtf.cds.gff | grep -c $j`
    if [ $p -ge 1 ]; then
        o=`echo Coding`
    else
        o=`echo Noncoding`
    fi
    q=`awk '$3=="transcript"' ${comD}/sqanti3/sqanti3_qc/pd1074.collapsed.filtered.rep_corrected.gtf.cds.gff | awk -v j="$j" '$12==j {print $10}'`
    r=`echo $q | sed 's/"//g' | sed 's/;//g'`
    echo -e "$i\t$k\t$l\t$m\t$o\t$r" >> ${workD}/pdspecific_2_transcript/pdspecific_2_transcript.notblasted.txt
done

list=`cat ${workD}/notblasted/u.l.notblasted.txt`
for i in $list; do
    j=$i\;
    k=`awk '$3=="mRNA"' /mnt/3main/lotusbo/isoseq.analysis/resources/pd1074/pd1074.reorganized_annotation.gff3 | awk -v j="$j" '$9~j {print $1":"$4"-"$5}'`
    l=`awk -v i="$i" '$1==i {print $2}' ${workD}/seq/pd1074.transcript.fasta.length`
    n=$i$
    m=`awk '$3=="exon"' /mnt/3main/lotusbo/isoseq.analysis/resources/pd1074/pd1074.reorganized_annotation.gff3 | grep -c $n`
    p=`cat ${workD}/seq/n2_CDS_WormBase.gff3 | grep -c $i`
    if [ $p -ge 1 ]; then
        o=`echo Coding`
    else
        o=`echo Noncoding`
    fi
    r=`awk '$3=="mRNA"' /mnt/3main/lotusbo/isoseq.analysis/resources/pd1074/pd1074.reorganized_annotation.gff3 | awk -v j="$j" '$9~j { split($9, arr, "="); print arr[3]}'`
    echo -e "$i\t$k\t$l\t$m\t$o\t$r" >> ${workD}/pdspecific_2_transcript/pdspecific_2_transcript.notblasted.txt
done

list=`cat ${workD}/notblasted/u.s.notblasted.txt`
for i in $list; do
    j=$i\;
    k=`awk '$3=="mRNA"' /mnt/3main/lotusbo/isoseq.analysis/resources/pd1074/pd1074.reorganized_annotation.gff3  | awk -v j="$j" '$9~j {print $1":"$4"-"$5}'`
    l=`awk -v i="$i" '$1==i {print $2}' ${workD}/seq/pd1074.transcript.fasta.length`
    n=$i$
    m=`awk '$3=="exon"' /mnt/3main/lotusbo/isoseq.analysis/resources/pd1074/pd1074.reorganized_annotation.gff3  | grep -c $n`
    p=`cat ${workD}/seq/n2_CDS_WormBase.gff3 | grep -c $i`
    if [ $p -ge 1 ]; then
        o=`echo Coding`
    else
        o=`echo Noncoding`
    fi
    r=`awk '$3=="mRNA"' /mnt/3main/lotusbo/isoseq.analysis/resources/pd1074/pd1074.reorganized_annotation.gff3 | awk -v j="$j" '$9~j { split($9, arr, "="); print arr[3]}'`
    echo -e "$i\t$k\t$l\t$m\t$o\t$r" >> ${workD}/pdspecific_2_transcript/pdspecific_2_transcript.notblasted.txt
done

# 05. Extracting protein sequences of the pd1074-specific genes for BLASTp and Batch CD search
list=`cat ${workD}/notblasted/s.l.notblasted.txt`
for i in $list; do
    j=`echo $i | sed 's/\./\\\./g'`
    echo ">$i" >> ${workD}/pdspecific_2_transcript/s.notblasted.faa
    grep -A 1 $j ${comD}/sqanti3//sqanti3_qc/pd1074.collapsed.filtered.rep_corrected.faa |sed '1d' >> ${workD}/pdspecific_2_transcript/s.notblasted.faa
done

list=`cat ${workD}/notblasted/u.l.notblasted.txt`
for i in $list; do
    j=`echo $i | sed 's/\./\\\./g' | sed 's/$/$/g'`
    grep -A 1 "$j" /mnt/3main/lotusbo/isoseq.analysis/resources/n2/c_elegans.PRJNA13758.WS264.protein_rename.fa >> u.l.notblasted.faa
done

list=`cat ${workD}/notblasted/u.s.notblasted.txt`
for i in $list; do
    j=`echo $i | sed 's/\./\\\./g' | sed 's/$/$/g'`
    grep -A 1 "$j" /mnt/3main/lotusbo/isoseq.analysis/resources/n2/c_elegans.PRJNA13758.WS264.protein_rename.fa >> u.s.notblasted.faa
done
