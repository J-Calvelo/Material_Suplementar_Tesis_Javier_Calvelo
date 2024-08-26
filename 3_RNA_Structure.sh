#!/usr/bin/env bash

# This script is meant top assist in the manual curation of SL RNA sequences.
# It is not a stand alone procedure.

# To run it it requires to create a folder named Structure_RNA_Manual located
# in Analysis_Results. And then create the input files with the sequence and
# location of the SM site in 2 columns and at least 2 lines
# Example:
#Clonorchis_sinensis_1_sec       AAUGUUCGGUUUUCUGCCGUGUAUAUUAGUGCACGGUAAUAAUCGACUCCGACCUAUGGUCGGAUGAAUUCUUUGGCUAGCCCACC
#Clonorchis_sinensis_1_sm1       ..................................................................xxxxxxxxxx..........
# Multiples sm sites can be provided
# File names must end with "_Input.tab"


dir_corrida=Analysis_Results/Structure_RNA_Manual
rm -r "Analysis_Results/Structure_RNA_Manual/RNAfold_"*

input_files=$(ls Analysis_Results/Structure_RNA_Manual"/"*"_Input.tab" | sed "s/.*\///" | sed "s/_Input.tab//")

for file in $input_files
do
  echo $file

  mkdir "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file
  mkdir "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/SM_Sites"
  mkdir "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/PS_Files"
  mkdir "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/Temporales"

  echo "SL Len GC sm num_sm MFE_energy centroid_energy centroid_distance MEA_energy MEA_Value mfe_structure_frec ensemble_diversity print_hairpin" | tr " " "\t" >> "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/RNAfold_Summary.tab"
  sl_rnas=$(awk -F "\t" '{print $1}' "Analysis_Results/Structure_RNA_Manual/"$file"_Input.tab" | grep _sec$ | sed "s/_sec//")

  echo "Analysis_Results/Structure_RNA_Manual/"$file"_Input.tab"
  for sl in $sl_rnas
  do
    echo $sl
    grep -P $sl"_sec\t" "Analysis_Results/Structure_RNA_Manual/"$file"_Input.tab" | sed "s/^/>/" | tr "\t" "\n" >> "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/Temporales/"$sl".fasta"
    seq_len=$(seqkit fx2tab -l -n "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/Temporales/"$sl".fasta" | tail -n 1 | awk -F "\t" '{print $2}')
    seq_gc=$(seqkit fx2tab -g -n "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/Temporales/"$sl".fasta" | tail -n 1 | awk -F "\t" '{print $2}')

    sm_sites=$(grep $sl "Analysis_Results/Structure_RNA_Manual/"$file"_Input.tab" | awk -F "\t" '{print $1}' | grep -v _sec$ | sed "s/$sl//" | sed "s/_//")
    num_sm_sites=$(grep $sl "Analysis_Results/Structure_RNA_Manual/"$file"_Input.tab" | awk -F "\t" '{print $1}' | grep -v _sec$ | grep -c .)

    for sm in $sm_sites
    do
      echo "Procesando sm: "$sm
      seqkit seq -w 0 "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/Temporales/"$sl".fasta" >> "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/SM_Sites/"$sl"_"$sm".work"
      grep -P $sl"_"$sm "Analysis_Results/Structure_RNA_Manual/"$file"_Input.tab" | awk -F "\t" '{print $2}' >> "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/SM_Sites/"$sl"_"$sm".work"

      RNAfold -p --MEA -C < "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/SM_Sites/"$sl"_"$sm".work" > "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/SM_Sites/"$sl"_"$sm"_work_estructura.out"

      ps2pdf $sl"_sec_dp.ps" "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/PS_Files/"$sl"_"$sm"_sec_dp.pdf"
      ps2pdf $sl"_sec_ss.ps" "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/PS_Files/"$sl"_"$sm"_sec_ss.pdf"

      rm $sl"_"*".ps"

      cat "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/SM_Sites/"$sl"_"$sm".work" | sed "/>/s/_sec/_$sm/" >> "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/Registro_estructura.txt"
      sed -n 6p "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/SM_Sites/"$sl"_"$sm"_work_estructura.out" | awk -F " " '{print $1" MEA"}' >> "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/Registro_estructura.txt"
      sed -n 3p "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/SM_Sites/"$sl"_"$sm"_work_estructura.out" | awk -F " " '{print $1" MFE"}' >> "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/Registro_estructura.txt"
      sed -n 4p "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/SM_Sites/"$sl"_"$sm"_work_estructura.out" | awk -F " " '{print $1" Calidad"}' >> "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/Registro_estructura.txt"
      echo "" >> "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/Registro_estructura.txt"

      MFE_energy=$(sed -n 3p "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/SM_Sites/"$sl"_"$sm"_work_estructura.out" | sed "s/ -/-/"| awk -F " " '{print $2}' | sed "s/(//" | sed "s/)//")
      pairwise_free_energy=$(sed -n 4p "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/SM_Sites/"$sl"_"$sm"_work_estructura.out" | sed "s/ -/-/" | awk -F " " '{print $2}' | sed "s/\[//" | sed "s/\]//")
      centroid_energy=$(sed -n 5p "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/SM_Sites/"$sl"_"$sm"_work_estructura.out" | sed "s/{ /{/" | awk -F " " '{print $2}' | sed "s/{//")
      centroid_distance=$(sed -n 5p "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/SM_Sites/"$sl"_"$sm"_work_estructura.out" | sed "s/{ /{/" | awk -F " " '{print $3}' | sed "s/}//" | sed "s/d=//")
      MEA_energy=$(sed -n 6p "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/SM_Sites/"$sl"_"$sm"_work_estructura.out" | sed "s/{ /{/" | sed "s/  / /g" | sed "s/{//" | awk -F " " '{print $2}')
      MEA_Value=$(sed -n 6p "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/SM_Sites/"$sl"_"$sm"_work_estructura.out" | sed "s/.*MEA=//" | sed "s/}$//")

      mfe_structure_frec=$(sed -n 7p "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/SM_Sites/"$sl"_"$sm"_work_estructura.out" | awk -F ";" '{print $1}' | sed "s/.* //")
      ensemble_diversity=$(sed -n 7p "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/SM_Sites/"$sl"_"$sm"_work_estructura.out" | awk -F ";" '{print $2}' | sed "s/ //g" | sed "s/ensemblediversity//")

      sed -n 6p "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/SM_Sites/"$sl"_"$sm"_work_estructura.out" | awk -F " " '{print $1}' | grep -o . >> "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/Temporales/"$sl"_"$sm"_work.temp"

      check_abort_1=$(grep -c "|" "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/Temporales/"$sl"_"$sm"_work.temp")
      check_abort_2=$(grep -c "{" "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/Temporales/"$sl"_"$sm"_work.temp")
      check_abort_3=$(grep -c "}" "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/Temporales/"$sl"_"$sm"_work.temp")

      if [ $check_abort_1 -eq 0 ] && [ $check_abort_2 -eq 0 ] && [ $check_abort_3 -eq 0 ]
      then
        print_hairpin=Har
        hairpin_num=1
        start_run=RUN_START
        count=1
        recorrer=$(grep -c . "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/Temporales/"$sl"_"$sm"_work.temp")
        count_hairpin=0
        len_hairpin=0

        while [ $count -le $recorrer ]
        do
          base_info=$(sed -n $count"p" "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/Temporales/"$sl"_"$sm"_work.temp")
          if [ "$base_info" == "(" ]
          then
            start_run=CONTEO
            count_hairpin=$(($count_hairpin + 1))
            len_hairpin=$(($len_hairpin + 1))
          elif [ "$base_info" == ")" ]
          then
            count_hairpin=$(($count_hairpin - 1))
            len_hairpin=$(($len_hairpin + 1))
          elif [ "$start_run" == "CONTEO" ]
          then
            len_hairpin=$(($len_hairpin + 1))
          fi

          if [ $count_hairpin -eq 0 ] && [ $start_run == "CONTEO" ]
          then
            print_hairpin=$(echo $print_hairpin$hairpin_num"_"$len_hairpin";Har")
            hairpin_num=$(($hairpin_num + 1))
            len_hairpin=0
            start_run=NUEVO
          fi
          count=$(($count + 1))
        done
      else
        print_hairpin=Manual_Check
      fi
      print_hairpin=$(echo $print_hairpin | sed "s/;Har$//" )
      echo $sl" "$seq_len" "$seq_gc" "$sm" "$num_sm_sites" "$MFE_energy" "$centroid_energy" "$centroid_distance" "$MEA_energy" "$MEA_Value" "$mfe_structure_frec" "$ensemble_diversity" "$print_hairpin | tr " " "\t" >> "Analysis_Results/Structure_RNA_Manual/RNAfold_"$file"/RNAfold_Summary.tab"
    done
  done
done
