#!/usr/bin/env bash

unique_SLRNA_structure_input=Selected_Sequences_Work/Unique_SL_Structure_Input.tab
unique_donnor_sites=Selected_Sequences_Work/Unique_SL_Donnor.fasta
sm_sites=sm

rm -r "Analysis_Results/Unique_SL_Structure"
mkdir "Analysis_Results/Unique_SL_Structure"
mkdir "Analysis_Results/Unique_SL_Structure/PS_Files"
mkdir "Analysis_Results/Unique_SL_Structure/SM_Sites"
mkdir "Analysis_Results/Unique_SL_Structure/Temporales"

echo "SL Len GC sm num_sm MFE_energy centroid_energy centroid_distance MEA_energy MEA_Value mfe_structure_frec ensemble_diversity print_hairpin Donnor_Site SM-Like_Site" | tr " " "\t" >> "Analysis_Results/Unique_SL_Structure/RNAfold_Summary.tab"

sl_rnas=$(awk -F "\t" '{print $1}' $unique_SLRNA_structure_input | grep _sec$ | sed "s/_sec//")
num_sm_sites=1
for sl in $sl_rnas
do
  echo $sl
  grep -P $sl"_sec\t" $unique_SLRNA_structure_input | sed "s/^/>/" | tr "\t" "\n" >> "Analysis_Results/Unique_SL_Structure/Temporales/"$sl".fasta"
  seq_len=$(seqkit fx2tab -l -n "Analysis_Results/Unique_SL_Structure/Temporales/"$sl".fasta" | tail -n 1 | awk -F "\t" '{print $2}')
  seq_gc=$(seqkit fx2tab -g -n "Analysis_Results/Unique_SL_Structure/Temporales/"$sl".fasta" | tail -n 1 | awk -F "\t" '{print $2}')
  donnor_site=$(seqkit grep -p $sl $unique_donnor_sites | seqkit seq -s)

  for sm in $sm_sites
  do
    echo "Procesando sm: "$sm
    sm_start=$(grep $sl"_"$sm $unique_SLRNA_structure_input | awk -F "\t" '{print $2}' | sed "s/./&\n/g" | grep -n x | head -n 1 | awk -F ":" '{print $1}')
    sm_end=$(grep $sl"_"$sm $unique_SLRNA_structure_input  | awk -F "\t" '{print $2}' | sed "s/./&\n/g" | grep -n x | tail -n 1 | awk -F ":" '{print $1}')

    sm_sequence=$(grep $sl"_sec" $unique_SLRNA_structure_input | awk -F "\t" '{print ">Temp\n"$2}' | seqkit subseq -r $sm_start":"$sm_end | tail -n 1)

    seqkit seq -w 0 "Analysis_Results/Unique_SL_Structure/Temporales/"$sl".fasta" >> "Analysis_Results/Unique_SL_Structure/SM_Sites/"$sl"_"$sm".work"
    grep -P $sl"_"$sm $unique_SLRNA_structure_input | awk -F "\t" '{print $2}' >> "Analysis_Results/Unique_SL_Structure/SM_Sites/"$sl"_"$sm".work"

    RNAfold -p --MEA -C < "Analysis_Results/Unique_SL_Structure/SM_Sites/"$sl"_"$sm".work" > "Analysis_Results/Unique_SL_Structure/SM_Sites/"$sl"_"$sm"_work_estructura.out"
    mv $sl"_"*"ps" "Analysis_Results/Unique_SL_Structure/PS_Files"

    cat "Analysis_Results/Unique_SL_Structure/SM_Sites/"$sl"_"$sm".work" | sed "/>/s/_sec/_$sm/" >> "Analysis_Results/Unique_SL_Structure/Registro_estructura.txt"
    sed -n 6p "Analysis_Results/Unique_SL_Structure/SM_Sites/"$sl"_"$sm"_work_estructura.out" | awk -F " " '{print $1" MEA"}' >> "Analysis_Results/Unique_SL_Structure/Registro_estructura.txt"
    sed -n 3p "Analysis_Results/Unique_SL_Structure/SM_Sites/"$sl"_"$sm"_work_estructura.out" | awk -F " " '{print $1" MFE"}' >> "Analysis_Results/Unique_SL_Structure/Registro_estructura.txt"
    sed -n 4p "Analysis_Results/Unique_SL_Structure/SM_Sites/"$sl"_"$sm"_work_estructura.out" | awk -F " " '{print $1" Calidad"}' >> "Analysis_Results/Unique_SL_Structure/Registro_estructura.txt"
    echo "" >> "Analysis_Results/Unique_SL_Structure/Registro_estructura.txt"

    MFE_energy=$(sed -n 3p "Analysis_Results/Unique_SL_Structure/SM_Sites/"$sl"_"$sm"_work_estructura.out" | sed "s/ -/-/"| awk -F " " '{print $2}' | sed "s/(//" | sed "s/)//")
    pairwise_free_energy=$(sed -n 4p "Analysis_Results/Unique_SL_Structure/SM_Sites/"$sl"_"$sm"_work_estructura.out" | sed "s/ -/-/" | awk -F " " '{print $2}' | sed "s/\[//" | sed "s/\]//")
    centroid_energy=$(sed -n 5p "Analysis_Results/Unique_SL_Structure/SM_Sites/"$sl"_"$sm"_work_estructura.out" | sed "s/{ /{/" | awk -F " " '{print $2}' | sed "s/{//")
    centroid_distance=$(sed -n 5p "Analysis_Results/Unique_SL_Structure/SM_Sites/"$sl"_"$sm"_work_estructura.out" | sed "s/{ /{/" | awk -F " " '{print $3}' | sed "s/}//" | sed "s/d=//")
    MEA_energy=$(sed -n 6p "Analysis_Results/Unique_SL_Structure/SM_Sites/"$sl"_"$sm"_work_estructura.out" | sed "s/{ /{/" | sed "s/  / /g" | sed "s/{//" | awk -F " " '{print $2}')
    MEA_Value=$(sed -n 6p "Analysis_Results/Unique_SL_Structure/SM_Sites/"$sl"_"$sm"_work_estructura.out" | sed "s/.*MEA=//" | sed "s/}$//")

    mfe_structure_frec=$(sed -n 7p "Analysis_Results/Unique_SL_Structure/SM_Sites/"$sl"_"$sm"_work_estructura.out" | awk -F ";" '{print $1}' | sed "s/.* //")
    ensemble_diversity=$(sed -n 7p "Analysis_Results/Unique_SL_Structure/SM_Sites/"$sl"_"$sm"_work_estructura.out" | awk -F ";" '{print $2}' | sed "s/ //g" | sed "s/ensemblediversity//")

    sed -n 6p "Analysis_Results/Unique_SL_Structure/SM_Sites/"$sl"_"$sm"_work_estructura.out" | awk -F " " '{print $1}' | grep -o . >> "Analysis_Results/Unique_SL_Structure/Temporales/"$sl"_"$sm"_work.temp"

    check_abort_1=$(grep -c "|" "Analysis_Results/Unique_SL_Structure/Temporales/"$sl"_"$sm"_work.temp")
    check_abort_2=$(grep -c "{" "Analysis_Results/Unique_SL_Structure/Temporales/"$sl"_"$sm"_work.temp")
    check_abort_3=$(grep -c "}" "Analysis_Results/Unique_SL_Structure/Temporales/"$sl"_"$sm"_work.temp")

    if [ $check_abort_1 -eq 0 ] && [ $check_abort_2 -eq 0 ] && [ $check_abort_3 -eq 0 ]
    then
      print_hairpin=Har
      hairpin_num=1
      start_run=RUN_START
      count=1
      recorrer=$(grep -c . "Analysis_Results/Unique_SL_Structure/Temporales/"$sl"_"$sm"_work.temp")
      count_hairpin=0
      len_hairpin=0

      while [ $count -le $recorrer ]
      do
        base_info=$(sed -n $count"p" "Analysis_Results/Unique_SL_Structure/Temporales/"$sl"_"$sm"_work.temp")
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
    echo $sl" "$seq_len" "$seq_gc" "$sm" "$num_sm_sites" "$MFE_energy" "$centroid_energy" "$centroid_distance" "$MEA_energy" "$MEA_Value" "$mfe_structure_frec" "$ensemble_diversity" "$print_hairpin" "$donnor_site" "$sm_sequence | tr " " "\t" >> "Analysis_Results/Unique_SL_Structure/RNAfold_Summary.tab"
  done
done
