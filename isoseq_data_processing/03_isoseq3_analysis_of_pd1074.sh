#!/bin/bash
# using pacbio ccs 4.2.0
# using pacbio lima 1.11.0 (commit v1.11.0)
# using pacbio isoseq3 3.3.0 (commit v3.3.0)

mkdir /mnt/3main/lotusbo/isoseq.analysis/pd1074/
mkdir /mnt/3main/lotusbo/isoseq.analysis/pd1074/isoseq3.one
export workDirectory="/mnt/3main/lotusbo/isoseq.analysis/pd1074/isoseq3.one"
export sname="pd1074"

## Isoseq3 analysis
# Reference: https://github.com/PacificBiosciences/IsoSeq_SA3nUP/wiki/Tutorial:-Installing-and-Running-Iso-Seq-3-using-Conda

# 01. Generate CCS reads
# When I made CCS read with min-passes 3, many reads were omitted. So, in this script, I changed the value of min-passes to 1.
mkdir ${workDirectory}/ccs
ccs /mnt/3main/junkim/isoseq/PD/m54229_191004_010418.subreads.bam ${workDirectory}/ccs/${sname}.isoseq.ccs.bam --min-rq 0.8 --min-passes 1 --num-threads 17
mv ccs_report.txt ${workDirectory}/ccs

# 02. Classify full-length read: library-prep-primer and barcode removal and demultiplexing
# Because my sample does not have barcode, I removed only library-prep-primer.
mkdir ${workDirectory}/demux
lima --isoseq --dump-clips --peek-guess -j 17 ${workDirectory}/ccs/${sname}.isoseq.ccs.bam  /mnt/3main/lotusbo/isoseq.analysis/resources/primer.libraryprep.fasta ${workDirectory}/demux/${sname}.isoseq.demux.bam

# 03. Classify full-length read: remove polyA tails and artificial concatamers
mkdir ${workDirectory}/flnc
isoseq3 refine --require-polya -j 17 ${workDirectory}/demux/${sname}.isoseq.demux.Clontech_5p--NEB_Clontech_3p.subreadset.xml /mnt/3main/lotusbo/isoseq.analysis/resources/primer.libraryprep.fasta ${workDirectory}/flnc/${sname}.flnc.bam

# 04. Cluster FLNC reads
mkdir ${workDirectory}/polish
isoseq3 cluster ${workDirectory}/flnc/${sname}.flnc.bam ${workDirectory}/polish/${sname}.polished.bam --verbose --use-qvs -j 17

# 05. Align reads to a genome using minimap2-2.17
mkdir ${workDirectory}/minimap2
minimap2 -t 17 -ax splice -uf --secondary=no -C5 -O6,24 -B4 /mnt/3main/lotusbo/isoseq.analysis/resources/pd1074/c_elegans.PRJEB28388.WS274.genomic.rename.fa ${workDirectory}/polish/${sname}.polished.hq.fasta.gz > ${workDirectory}/minimap2/${sname}.hq_isoforms.fasta.sam 2> ${workDirectory}/minimap2/${sname}.hq_isoforms.fasta.sam.log
sort -k 3,3 -k 4,4n ${workDirectory}/minimap2/${sname}.hq_isoforms.fasta.sam > ${workDirectory}/minimap2/${sname}.hq_isoforms.fasta.sorted.sam

# 06. Collapse redundant isoforms
mkdir ${workDirectory}/collapse
gunzip ${workDirectory}/polish/${sname}.polished.hq.fasta.gz
/mnt/main/tools/cDNA_cupcake/bin/collapse_isoforms_by_sam.py --input ${workDirectory}/polish/${sname}.polished.hq.fasta -s ${workDirectory}/minimap2/${sname}.hq_isoforms.fasta.sorted.sam --dun-merge-5-shorter -o ${workDirectory}/collapse/${sname} > ${workDirectory}/collapse/log.txt

# 07. Obtain associated count information
/mnt/main/tools/cDNA_cupcake/bin/get_abundance_post_collapse.py ${workDirectory}/collapse/${sname}.collapsed ${workDirectory}/polish/${sname}.polished.cluster_report.csv

# 08. Filter away 5' degraded isoforms
# Since my sample treated with Clontech SMARTer cDNA kit, I conducted this process.
/mnt/main/tools/cDNA_cupcake/bin/filter_away_subset.py ${workDirectory}/collapse/${sname}.collapsed
