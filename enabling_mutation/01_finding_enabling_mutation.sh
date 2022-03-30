#!/bin/bash

### Finding enabling mutations

mkdir 25_enabling_mutation
mkdir 25_enabling_mutation/pd1074
mkdir 25_enabling_mutation/cb4856

## using non-coding transcripts of PD1074 (I used only Iso-Seq results)

# 01. Extract sequences of the non-coding transcripts of PD1074
awk '$30=="non_coding" {print $1}' pd1074/sqanti3/sqanti3_qc/pd1074.collapsed.filtered.rep_classification.txt \
    >> 25_enabling_mutation/pd1074/pd1074_noncoding_list.txt
for i in `cat 25_enabling_mutation/pd1074/pd1074_noncoding_list.txt`
do
    bioawk -c fastx -v var="$i" '{if(index($name, var)!=0) print ">"$name"\n"$seq}' pd1074/isoseq3.one/collapse/pd1074.collapsed.filtered.rep.fa \
        >> 25_enabling_mutation/pd1074/pd1074_noncoding.fa
done

# 02. BLAST search of the non-coding transcripts of PD1074 to the CB4856 genome
blastn -query 25_enabling_mutation/pd1074/pd1074_noncoding.fa -db /mnt/3main/lotusbo/isoseq.analysis/resources/blast_db/cb4856_genome/cb4856_genome \
    -outfmt 6 -num_threads 30 -out 25_enabling_mutation/pd1074/pd1074_noncoding_cb4856genome.out

# 03. Identifying non-coding genes which have similarity with coding-genes of CB4856
# Among the above BLAST results, I extracted the results with >95% identity and >90% coverage with CB4856 genes and manually checked whether the CB4856 gene is a coding or non-coding gene.



## using non-coding transcripts of CB4856 (I used only Iso-Seq results)

# 01. Extract sequences of the  non-coding transcripts of CB4856
awk '$30=="non_coding" {print $1}' cb4856/sqanti3/sqanti3_qc/cb4856.collapsed.filtered.rep_classification.txt \
    >> 25_enabling_mutation/cb4856/cb4856_noncoding_list.txt
for i in `cat 25_enabling_mutation/cb4856/cb4856_noncoding_list.txt`
do
    j=`echo $i | sed 's/$/|/g'`
    bioawk -c fastx -v var="$j" '{if(index($name, var)!=0) print ">"$name"\n"$seq}' cb4856/isoseq3.one/collapse/cb4856.collapsed.filtered.rep.fa \
        >> 25_enabling_mutation/cb4856/cb4856_noncoding.fa
done

# 02. BLAST search of the non-coding transcripts of CB4856 to the PD1074 genome
blastn -query 25_enabling_mutation/cb4856/cb4856_noncoding.fa -db /mnt/3main/lotusbo/isoseq.analysis/resources/blast_db/pd1074_genome/pd1074_genome \
    -outfmt 6 -num_threads 30 -out 25_enabling_mutation/cb4856/cb4856_noncoding_pd1074genome.out

# 03. Identifying non-coding genes which have similarity with coding-genes of PD1074
Among the above BLAST results, I extracted the results with >95% identity and >90% coverage with PD1074 genes and manually checked whether the PD1074 gene is a coding or non-coding gene.
