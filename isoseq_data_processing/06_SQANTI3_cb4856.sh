
#!/bin/bash

mkdir /mnt/3main/lotusbo/isoseq.analysis/cb4856/isoseq3.one/liftOver
liftOver -gff /mnt/3main/lotusbo/isoseq.analysis/resources/pd1074/pd1074.reorganized_annotation.gff3 \
    /mnt/3main/lotusbo/isoseq.analysis/liftover/pdtocb/pd1074Tocb4856.over.chain.gz \
    /mnt/3main/lotusbo/isoseq.analysis/cb4856/isoseq3.one/liftOver/pd1074.reorganized_annotation.ptcliftOver.gff3 \
    /mnt/3main/lotusbo/isoseq.analysis/cb4856/isoseq3.one/liftOver/pd1074.reorganized_annotation.ptcunmapped.gff3
gffread /mnt/3main/lotusbo/isoseq.analysis/cb4856/isoseq3.one/liftOver/pd1074.reorganized_annotation.ptcliftOver.gff3 \
 -T -o /mnt/3main/lotusbo/isoseq.analysis/cb4856/isoseq3.one/liftOver/pd1074.reorganized_annotation.ptcliftOver.gtf

mkdir /mnt/3main/lotusbo/isoseq.analysis/cb4856/sqanti3
mkdir /mnt/3main/lotusbo/isoseq.analysis/cb4856/sqanti3/sqanti3_qc

source activate SQANTI3.env
export PYTHONPATH=$PYTHONPATH:/mnt/3main/lotusbo/tools/cDNA_Cupcake/sequence/
export PYTHONPATH=$PYTHONPATH:/mnt/3main/lotusbo/tools/cDNA_Cupcake/

cd /mnt/3main/lotusbo/isoseq.analysis/cb4856/sqanti3/sqanti3_qc
python /mnt/3main/lotusbo/tools/SQANTI3/sqanti3_qc.py /mnt/3main/lotusbo/isoseq.analysis/cb4856/isoseq3.one/collapse/cb4856.collapsed.filtered.rep.fa \
    /mnt/3main/lotusbo/isoseq.analysis/cb4856/isoseq3.one/liftOver/pd1074.reorganized_annotation.ptcliftOver.gtf \
    /mnt/3main/lotusbo/isoseq.analysis/resources/cb4856/cb4856.scaffold.gapfilling.mt.wrapped.rename.fasta \
    --aligner_choice=minimap2 \
    --fl_count /mnt/3main/lotusbo/isoseq.analysis/cb4856/isoseq3.one/collapse/cb4856.collapsed.filtered.abundance.txt
