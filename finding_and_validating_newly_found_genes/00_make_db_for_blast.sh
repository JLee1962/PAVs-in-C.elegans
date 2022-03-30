#!/bin/bash

mkdir /mnt/3main/lotusbo/isoseq.analysis/resources/blast_db
export comD="/mnt/3main/lotusbo/isoseq.analysis/resources"
export workD="/mnt/3main/lotusbo/isoseq.analysis/resources/blast_db"

export sname="pd1074_transcript"
mkdir ${workD}/${sname}
gffread -w ${workD}/${sname}/pd1074.reorganized_annotation.all_transcripts.fa -g ${comD}/pd1074/c_elegans.PRJEB28388.WS274.genomic.rename.fa ${comD}/pd1074/pd1074.reorganized_annotation.gff3
dustmasker -in ${workD}/${sname}/pd1074.reorganized_annotation.all_transcripts.fa -infmt fasta -parse_seqids -outfmt maskinfo_asn1_bin -out ${workD}/${sname}/${sname}_dust.asnb
makeblastdb -in ${workD}/${sname}/pd1074.reorganized_annotation.all_transcripts.fa -input_type fasta -dbtype nucl -parse_seqids -mask_data ${workD}/${sname}/${sname}_dust.asnb -out ${workD}/${sname}/${sname} -title "Caenorhabditis elegans, pd1074 strain, transcript"

export sname="pd1074_isoseq3"
mkdir ${workD}/${sname}
bioawk -c fastx '{ split ($name, arr, "|"); print ">"arr[1] ; print $seq}'  /mnt/3main/lotusbo/isoseq.analysis/pd1074/isoseq3.one/collapse/pd1074.collapsed.rep.fa > ${workD}/${sname}/pd1074.collapsed.rep.rename.fa
dustmasker -in ${workD}/${sname}/pd1074.collapsed.rep.rename.fa -infmt fasta -parse_seqids -outfmt maskinfo_asn1_bin -out ${workD}/${sname}/${sname}_dust.asnb
makeblastdb -in ${workD}/${sname}/pd1074.collapsed.rep.rename.fa -input_type fasta -dbtype nucl -parse_seqids -mask_data ${workD}/${sname}/${sname}_dust.asnb -out ${workD}/${sname}/${sname} -title "Caenorhabditis elegans, pd1074 strain, transcript_lab"

export sname="n2_mRNA_transcript"
mkdir ${workD}/${sname}
dustmasker -in ${comD}/n2/c_elegans.PRJNA13758.WS274.mRNA_transcripts.fa -infmt fasta -parse_seqids -outfmt maskinfo_asn1_bin -out ${workD}/${sname}/${sname}_dust.asnb
makeblastdb -in ${comD}/n2/c_elegans.PRJNA13758.WS274.mRNA_transcripts.fa -input_type fasta -dbtype nucl -parse_seqids -mask_data ${workD}/${sname}/${sname}_dust.asnb -out ${workD}/${sname}/${sname} -title "Caenorhabditis elegans, N2 mRNA_transcripts"

export sname="n2_ncRNA_transcript"
mkdir ${workD}/${sname}
dustmasker -in ${comD}/n2/c_elegans.PRJNA13758.WS274.ncRNA_transcripts.fa -infmt fasta -parse_seqids -outfmt maskinfo_asn1_bin -out ${workD}/${sname}/${sname}_dust.asnb
makeblastdb -in ${comD}/n2/c_elegans.PRJNA13758.WS274.ncRNA_transcripts.fa -input_type fasta -dbtype nucl -parse_seqids -mask_data ${workD}/${sname}/${sname}_dust.asnb -out ${workD}/${sname}/${sname} -title "Caenorhabditis elegans, N2 ncRNA_transcripts"

export sname="n2_pseudogenic_transcript"
mkdir ${workD}/${sname}
dustmasker -in ${comD}/n2/c_elegans.PRJNA13758.WS274.pseudogenic_transcripts.fa -infmt fasta -parse_seqids -outfmt maskinfo_asn1_bin -out ${workD}/${sname}/${sname}_dust.asnb
makeblastdb -in ${comD}/n2/c_elegans.PRJNA13758.WS274.pseudogenic_transcripts.fa -input_type fasta -dbtype nucl -parse_seqids -mask_data ${workD}/${sname}/${sname}_dust.asnb -out ${workD}/${sname}/${sname} -title "Caenorhabditis elegans, N2 pseudogenic_transcripts"

export sname="n2_transposon_transcript"
mkdir ${workD}/${sname}
dustmasker -in ${comD}/n2/c_elegans.PRJNA13758.WS274.transposon_transcripts.fa -infmt fasta -parse_seqids -outfmt maskinfo_asn1_bin -out ${workD}/${sname}/${sname}_dust.asnb
makeblastdb -in ${comD}/n2/c_elegans.PRJNA13758.WS274.transposon_transcripts.fa -input_type fasta -dbtype nucl -parse_seqids -mask_data ${workD}/${sname}/${sname}_dust.asnb -out ${workD}/${sname}/${sname} -title "Caenorhabditis elegans, N2 transposon_transcripts"
