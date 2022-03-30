#!/bin/bash

mkdir /mnt/3main/lotusbo/isoseq.analysis/pd1074/sqanti3
mkdir /mnt/3main/lotusbo/isoseq.analysis/pd1074/sqanti3/sqanti3_qc


# Running SQANTI3 quality control
source activate SQANTI3.env
export PYTHONPATH=$PYTHONPATH:/mnt/3main/lotusbo/tools/cDNA_Cupcake/sequence/
export PYTHONPATH=$PYTHONPATH:/mnt/3main/lotusbo/tools/cDNA_Cupcake/

gffread /mnt/3main/lotusbo/isoseq.analysis/resources/pd1074/pd1074.reorganized_annotation.gff3 \
 -T -o /mnt/3main/lotusbo/isoseq.analysis/resources/pd1074/pd1074.reorganized_annotation.gtf

cd /mnt/3main/lotusbo/isoseq.analysis/pd1074/sqanti3/sqanti3_qc
python /mnt/3main/lotusbo/tools/SQANTI3/sqanti3_qc.py /mnt/3main/lotusbo/isoseq.analysis/pd1074/isoseq3.one/collapse/pd1074.collapsed.filtered.rep.fa \
    /mnt/3main/lotusbo/isoseq.analysis/resources/pd1074/pd1074.reorganized_annotation.gtf \
    /mnt/3main/lotusbo/isoseq.analysis/resources/pd1074/c_elegans.PRJEB28388.WS274.genomic.rename.fa \
    --aligner_choice=minimap2 \
    --fl_count /mnt/3main/lotusbo/isoseq.analysis/pd1074/isoseq3.one/collapse/pd1074.collapsed.filtered.abundance.txt
