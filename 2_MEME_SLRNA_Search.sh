#!/usr/bin/env bash

# This script is meant top assist in the manual curation of SL RNA sequences.
# It is not a stand alone procedure.

# To run this script first it is necesary to preselect promising SLRNA sequences
# and store them in fasta files in a folder named Manual_Curation located in
# Analysis_Results

archivos=$(ls "Analysis_Results/Manual_Curation_MEME/"*".fasta" | sed "s/.fasta//" | sed "s/.*\///")
for file in $archivos
do
  rm -r "Analysis_Results/Manual_Curation_MEME/"$file"_meme_oops"
  meme "Analysis_Results/Manual_Curation_MEME/"$file".fasta" -rna -oc "Analysis_Results/Manual_Curation_MEME/"$file"_meme_oops" -time 30000 -nostatus -mod oops -nmotifs 8 -minw 5 -maxw 50
done
