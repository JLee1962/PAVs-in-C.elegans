#!/bin/bash

# Basecalling of nanopore sequencing (ONT) of the CB4856 genomic DNA
/mnt/main/tools/ont-guppy-cpu/bin/guppy_basecaller \
--input_path /mnt/3main/lotusbo/cb4856.genome.assembly/raw_data/cb4856/20190917_0636_MN28424_FAL17593_179f4c20/fast5/ \
--save_path /mnt/3main/lotusbo/cb4856.genome.assembly/raw_data/basecall.fastq/ --flowcell FLO-MIN106 --kit SQK-LSK109 \
--cpu_threads_per_caller 10 > /mnt/3main/lotusbo/cb4856.genome.assembly/raw_data/log.basecalling
