#!/bin/bash
#PBS -q normal
#PBS -l walltime=01:00:00
#PBS -l mem=64GB
#PBS -l ncpus=16
#PBS -j oe

#script copied from kango2/ausarg

module load minimap2 samtools pythonlib

# Check if the task has already been completed
checkpoint_file="${outputdir}/${sampleid}.ribocop"

if [ -f "$checkpoint_file.done" ];
then
    echo "Task already completed. Exiting."
    exit 0
else
    touch "${checkpoint_file}.running"
fi

set -ex
# Define the deuterostomia rDNA file as a variable
rdnalibfa="/g/data/xl04/hrp561/rdnalib/deuterostomia.18s28s.rDNA.fasta"

# Define the sample id as a variable
# sampleid="POGVITdef.h1.corrected_YAHS"

# Define the input fasta file as a variable
# inputfasta="POGVITdef.h1.corrected_YAHS.fasta"

# Define the output directory as a variable
# outputdir="/g/data/xl04/genomeprojects/Pogona_vitticeps/tmp/rdna/"

mkdir -p $outputdir
cd $outputdir

# Run minimap2 to generate primary PAF file
minimap2 -t ${PBS_NCPUS} --secondary=no -o ${sampleid}.primary.paf ${rdnalibfa} ${inputfasta}

# Extract the most frequent 18S and 28S sequences
eighteen=$(awk '$11>=$7*0.9' ${sampleid}.primary.paf | cut -f6 | sort | uniq -c | grep _18S | sort -k1,1n | tail -n 1 | sed -r 's/^\s+//g' | cut -f2 -d ' ')
twoeight=$(awk '$11>=$7*0.9' ${sampleid}.primary.paf | cut -f6 | sort | uniq -c | grep _28S | sort -k1,1n | tail -n 1 | sed -r 's/^\s+//g' | cut -f2 -d ' ')

# Extract sequences from the fasta file
samtools faidx $rdnalibfa $eighteen > ${sampleid}.ssulsurdna.fa
samtools faidx $rdnalibfa $twoeight >> ${sampleid}.ssulsurdna.fa

# Run minimap2 to generate refined PAF file
minimap2 -t ${PBS_NCPUS} --secondary=no -o ${sampleid}.refined.paf ${sampleid}.ssulsurdna.fa ${inputfasta}

# Run the rdnabuilder.py script with the refined PAF file
python3 /g/data/te53/hrp561/ausarg/scripts/python/rdnabuilder.py -p ${sampleid}.refined.paf -f ${inputfasta} -o ./ -s ${sampleid}

# Map assembly to the rDNA morphs to generate bed file for annotations
# This won't be perfect as the morphs may be incomplete at the edge of contigs
# Map assembly to rDNA morphs
# minimap2 was failing due to memory issues. may be due to too many identical sequences
# select one morph sequence and use that for alignments

# Select the morph sequence
samtools faidx ${sampleid}.rDNA.morphs.fasta
refrdnamorph=$(perl -lne '@a=split("\t",$_); push(@l,$a[1]); $s{$a[1]}=$a[0]; END { @l = sort {$a<=>$b} @l; print $s{@l[int($#l/2)]}}' ${sampleid}.rDNA.morphs.fasta.fai)
samtools faidx ${sampleid}.rDNA.morphs.fasta $refrdnamorph > ${sampleid}.rDNA.refmorph.fasta
# map against the selected morph
minimap2 -t ${PBS_NCPUS} --secondary=no -o ${sampleid}.asm2refmorph.paf ${sampleid}.rDNA.refmorph.fasta ${inputfasta}

# Create a unique temporary file
tempscript=$(mktemp tempfile.XXXXXXXX)

# Create R script for processing the PAF output
cat << 'EOF' > ${tempscript}.R
library(tidyverse)
library(GenomicRanges)
args <- commandArgs(trailingOnly = TRUE)
paf_file <- args[1]
output_file <- args[2]
alnratio <- 0.8
minalnresidue <- 500
edge <- 100
gapwidth <- 100
tmp <- read_delim(paf_file, delim = "\t", col_names = F, col_select = c(1:12))
minregionlength <- tmp %>% filter(X10>minalnresidue & X10/X11>alnratio) %>% pull(X10) %>% median() / 2
gr <- makeGRangesFromDataFrame(tmp %>% 
    filter(X10>minalnresidue & X10/X11>alnratio) %>% 
    mutate(X3=case_when(X3<edge~0,TRUE~X3), 
           X4=case_when(X2-X4<edge~X2,TRUE~X4)) %>% 
    select(chr=X1,start=X3,end=X4))

write_delim(
    data.frame(GenomicRanges::reduce(gr, min.gapwidth=gapwidth)) %>% 
    filter(width > minregionlength) %>%
    select(c(1:3)), 
    output_file, 
    col_names = F, 
    delim = "\t"
)
EOF

# Run R script
module load Rlib/4.3.1
Rscript ${tempscript}.R "${sampleid}.asm2refmorph.paf" "${sampleid}.rdnaregions.bed"
rm ${tempscript}.R
# Clean up
# rm getrDNAregions.R

# Mark the task as completed
touch "${checkpoint_file}.done"
rm -f "${checkpoint_file}.running"

exit 0

# TODO:
# - create a contained repo for this. perhaps include these commands in rdnabuilder.py script
# - pangenome graph for rDNA morphs
# - align assembly to morphs to ensure that all morphs are captured
# - annotate rDNA regions with repeatmasker
# - create reporting stats such as 
#   - length distribution, 
#   - completeness of morphs,
#   - head-to-tail organisation,
#   - number of scaffolds with rDNA seq only and their names
#   - number of rDNA copies in the assembly
#   - sequence divergence metrics between morphs
#   - multiple sequence alignments
#   - phylogenetic trees
