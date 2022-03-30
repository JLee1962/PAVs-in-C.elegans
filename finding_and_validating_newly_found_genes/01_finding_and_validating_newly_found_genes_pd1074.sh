#!/bin/bash

## Finding newly found genes

# 01. Making a list of newly found gene candidate
mkdir /mnt/3main/lotusbo/isoseq.analysis/pd1074/sqanti3/novel_gene
export comD="/mnt/3main/lotusbo/isoseq.analysis/pd1074"
export workD="/mnt/3main/lotusbo/isoseq.analysis/pd1074/sqanti3/novel_gene"

cat ${comD}/sqanti3/sqanti3_qc/pd1074.collapsed.filtered.rep_classification.txt | awk '{ if (index($7,"novelGene")!=0) print $0}' \
    >> ${workD}/novel_gene_candidate_list

list=`cat ${workD}/novel_gene_candidate_list | awk '{print $1}'`
for i in $list; do
    cat ${comD}/isoseq3.one/collapse/pd1074.collapsed.read_stat.txt | grep "$i$" >> ${workD}/novel_gene_read_list
done

list=`cat ${workD}/novel_gene_candidate_list | awk '{print $1}' | sed 's/$/";/g'`
for i in $list; do
    cat ${comD}/sqanti3/sqanti3_qc/pd1074.collapsed.filtered.rep_corrected.gtf.cds.gff | grep $i >> ${workD}/novel_gene_candidate.gff
done


# 02_Comparing newly found gene candidates with N2 transcripts and known PD1074 transcripts
samtools view ${comD}/isoseq3.one/ccs/pd1074.isoseq.ccs.bam >> ${comD}/isoseq3.one/ccs/pd1074.isoseq.ccs.sam

list=`cat ${workD}/novel_gene_read_list | awk '{print $1}'`
for i in $list; do
    echo $i | awk -v i="$i" '{ if ($1==i) {print ">"$1; print $10}}' ${comD}/isoseq3.one/ccs/pd1074.isoseq.ccs.sam >> ${workD}/novel_gene_ccs.fa
done

minimap2 -ax splice -t 4 /mnt/3main/lotusbo/isoseq.analysis/resources/pd1074/c_elegans.PRJEB28388.WS274.genomic.rename.fa ${workD}/novel_gene_ccs.fa | samtools sort -o ${workD}/novel_gene.bam
samtools index ${workD}/novel_gene.bam

mkdir ${workD}/blast
blastn -query ${workD}/novel_gene_ccs.fa -db /mnt/3main/lotusbo/isoseq.analysis/resources/blast_db/pd1074_transcript/pd1074_transcript \
       -outfmt  "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \
       -num_threads 20 -out ${workD}/blast/novel_gene_ccs.blast.pd1074_transcript.out

blastn -query ${workD}/novel_gene_ccs.fa -db /mnt/3main/lotusbo/isoseq.analysis/resources/blast_db/n2_mRNA_transcript/n2_mRNA_transcript \
       -outfmt  "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \
       -num_threads 20 -out ${workD}/blast/novel_gene_ccs.blast.n2_mRNA_transcript.out

blastn -query ${workD}/novel_gene_ccs.fa -db /mnt/3main/lotusbo/isoseq.analysis/resources/blast_db/n2_ncRNA_transcript/n2_ncRNA_transcript \
       -outfmt  "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \
       -num_threads 20 -out ${workD}/blast/novel_gene_ccs.blast.n2_ncRNA_transcript.out

blastn -query ${workD}/novel_gene_ccs.fa -db /mnt/3main/lotusbo/isoseq.analysis/resources/blast_db/n2_pseudogenic_transcript/n2_pseudogenic_transcript \
       -outfmt  "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \
       -num_threads 20 -out ${workD}/blast/novel_gene_ccs.blast.n2_pseudogenic_transcript.out

blastn -query ${workD}/novel_gene_ccs.fa -db /mnt/3main/lotusbo/isoseq.analysis/resources/blast_db/n2_transposon_transcript/n2_transposon_transcript \
       -outfmt  "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \
       -num_threads 20 -out ${workD}/blast/novel_gene_ccs.blast.n2_transposon_transcript.out

# 03. Making a list of newly found gene candidates that do not have similarities in the above analysis
list=`awk '{print $1}' ${workD}/blast/novel_gene_ccs.blast.pd1074_transcript.out | sort -u`
for i in $list ; do
    awk -v i="$i" '$1==i {print $5}' ${workD}/novel_gene_read_list >> ${workD}/blasted.txt
done

list=`awk '{print $1}' ${workD}/blast/novel_gene_ccs.blast.n2_mRNA_transcript.out | sort -u`
for i in $list ; do
    awk -v i="$i" '$1==i {print $5}' ${workD}/novel_gene_read_list >> ${workD}/blasted.txt
done

list=`awk '{print $1}' ${workD}/blast/novel_gene_ccs.blast.n2_ncRNA_transcript.out | sort -u`
for i in $list ; do
    awk -v i="$i" '$1==i {print $5}' ${workD}/novel_gene_read_list >> ${workD}/blasted.txt
done

list=`awk '{print $1}' ${workD}/blast/novel_gene_ccs.blast.n2_pseudogenic_transcript.out | sort -u`
for i in $list ; do
    awk -v i="$i" '$1==i {print $5}' ${workD}/novel_gene_read_list >> ${workD}/blasted.txt
done

list=`awk '{print $1}' ${workD}/blast/novel_gene_ccs.blast.n2_transposon_transcript.out | sort -u`
for i in $list ; do
    awk -v i="$i" '$1==i {print $5}' ${workD}/novel_gene_read_list >> ${workD}/blasted.txt
done

cat ${workD}/blasted.txt | sort -u >> ${workD}/blasted.tmp.txt
rm ${workD}/blasted.txt
mv ${workD}/blasted.tmp.txt ${workD}/blasted.txt

cp ${workD}/novel_gene_candidate_list ${workD}/notblasted.novel_gene_candidate_list.txt
list=`cat ${workD}/blasted.txt`
for i in $list ; do
    awk -v  i="$i" '$1!=i' ${workD}/notblasted.novel_gene_candidate_list.txt >> ${workD}/temp.txt
    rm ${workD}/notblasted.novel_gene_candidate_list.txt
    mv ${workD}/temp.txt ${workD}/notblasted.novel_gene_candidate_list.txt
done
cat ${workD}/notblasted.novel_gene_candidate_list.txt | sort -u >> ${workD}/test.tmp.txt
mv ${workD}/test.tmp.txt ${workD}/notblasted.novel_gene_candidate_list.txt

# 04. After finding newly found gene candidates, I listed newly found genes by manually removing the transcript of PD1074 found in the Yoshimura et al. (2019) paper.





## Validating newly found genes

mkdir ${workD}/novel.spec

# 01. Extracting gff of newly found genes
list=`awk '{print $1}' ${workD}/notblasted.novel_gene_candidate_list.txt | sed 's/PB/"PB/g' | sed 's/$/";/g'`
for i in $list; do
    awk -v i="$i" '$12==i' ${workD}/novel_gene_candidate.gff >> ${workD}/novel.spec/notblasted.novel_gene_candidate.gff
done

# 02. Finding the exon number, starting site, ending site, and abundance of newly found genes
bioawk -c fastx '{ print $name"\t"length($seq) }' ${comD}/sqanti3/sqanti3_qc/pd1074.collapsed.filtered.rep_corrected.fasta \
    > ${comD}/sqanti3/sqanti3_qc/pd1074.collapsed.filtered.rep_corrected.fasta.length

echo -e "name\tlength\tnumberofexon\tchr\tstart\tend\tabun_num\tabun_rate" >> ${workD}/novel.spec/notblasted.novel_gene_candidate.txt
list=`awk '{print $1}' ${workD}/notblasted.novel_gene_candidate_list.txt`
for i in $list; do
    namelen=`awk -v i="$i" '$1==i' ${comD}/sqanti3/sqanti3_qc/pd1074.collapsed.filtered.rep_corrected.fasta.length`
    j=$i\"\;
    numofexon=`awk '$3=="exon"' ${workD}/novel.spec/notblasted.novel_gene_candidate.gff | grep -c $j`
    chrstartend=`awk '$3=="transcript"' ${workD}/novel.spec/notblasted.novel_gene_candidate.gff | grep $j | awk '{print $1"\t"$4"\t"$5}'`
    abun=`awk -v i="$i" '$1==i {print $2"\t"$3}' ${comD}/isoseq3.one/collapse/pd1074.collapsed.abundance.txt`
    echo -e "$namelen\t$numofexon\t$chrstartend\t$abun" >> ${workD}/novel.spec/notblasted.novel_gene_candidate.txt
done

# 03. Extracting the protein sequence of newly found genes for blastp search and Batch CD-search
gffread -y ${workD}/novel.spec/notblasted.novel_gene_candidate.faa -g /mnt/3main/lotusbo/isoseq.analysis/resources/pd1074/c_elegans.PRJEB28388.WS274.genomic.rename.fa \
    ${workD}/novel.spec/notblasted.novel_gene_candidate.gff
