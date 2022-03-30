#!/bin/bash

export comD="/mnt/3main/lotusbo/isoseq.analysis/resources"
export workD="/mnt/3main/lotusbo/isoseq.analysis/resources/blast_db"

export sname="pd1074_genome"
mkdir ${workD}/${sname}
dustmasker -in ${comD}/pd1074/c_elegans.PRJEB28388.WS274.genomic.rename.fa -infmt fasta -parse_seqids -outfmt maskinfo_asn1_bin -out ${workD}/${sname}/${sname}_dust.asnb
makeblastdb -in ${comD}/pd1074/c_elegans.PRJEB28388.WS274.genomic.rename.fa -input_type fasta -dbtype nucl -parse_seqids -mask_data ${workD}/${sname}/${sname}_dust.asnb -out ${workD}/${sname}/${sname} -title "Caenorhabditis elegans, PD1074 strain"

export sname="cb4856_genome"
mkdir ${workD}/${sname}
dustmasker -in ${comD}/cb4856/cb4856.scaffold.gapfilling.mt.wrapped.rename.fasta -infmt fasta -parse_seqids -outfmt maskinfo_asn1_bin -out ${workD}/${sname}/${sname}_dust.asnb
makeblastdb -in ${comD}/cb4856/cb4856.scaffold.gapfilling.mt.wrapped.rename.fasta -input_type fasta -dbtype nucl -parse_seqids -mask_data ${workD}/${sname}/${sname}_dust.asnb -out ${workD}/${sname}/${sname} -title "Caenorhabditis elegans, CB4856 strain"

mkdir /mnt/3main/lotusbo/isoseq.analysis/resources/PRJNA692613/ncbi-genomes-2021-04-26/blast_db
export comD="/mnt/3main/lotusbo/isoseq.analysis/resources/PRJNA692613/ncbi-genomes-2021-04-26"
export workD="/mnt/3main/lotusbo/isoseq.analysis/resources/PRJNA692613/ncbi-genomes-2021-04-26/blast_db"

list=`ls -l /mnt/3main/lotusbo/isoseq.analysis/resources/PRJNA692613/ncbi-genomes-2021-04-26 | awk '$9~"GCA" {print $9}'`
for i in $list; do
    sname=`echo $i | awk '{split ($1, arr, "_"); print arr[3]}' | sed 's/$/_genome/g'`
    echo $i
    echo ${sname}
    mkdir ${workD}/${sname}
    makeblastdb -in ${comD}/$i  -input_type fasta -dbtype nucl -parse_seqids -out ${workD}/${sname}_db/${sname} -title "Caenorhabditis elegans, $j genome"
done

makeblastdb -in ${comD}/pd_intermediate_seq.fasta -input_type fasta -dbtype nucl -parse_seqids -out ${workD}/pd_intermediate_seq/pd_intermediate_seq -title "pd_intermediate_seq"

makeblastdb -in ${comD}/cb_intermediate_seq.fasta -input_type fasta -dbtype nucl -parse_seqids -out ${workD}/cb_intermediate_seq/cb_intermediate_seq -title "cb_intermediate_seq"

makeblastdb -in ${comD}/pd_intermediate_3.fasta -input_type fasta -dbtype nucl -parse_seqids -out ${workD}/pd_intermediate_3/pd_intermediate_3 -title "pd_intermediate_seq"

