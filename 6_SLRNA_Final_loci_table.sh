#!/usr/bin/env bash

selected_loci_work=Selected_Sequences_Work/Selected_work_loci.txt
selected_reference_work=Selected_Sequences_Work/Selected_work_reference_loci.txt
selected_on_recovery_blast=Selected_Sequences_Work/Selected_recovery_blast.txt
unique_SLRNA=Selected_Sequences_Work/Unique_SL_Sequences.fasta
unique_structure=Selected_Sequences_Work/Unique_SL_RNA_structure_registry.txt

FINAL_LOCI_TABLE=FALSE
FINAL_SLRNA_ALN=TRUE

if [ "$FINAL_LOCI_TABLE" == "TRUE" ]
then
  count=2
  recorrer=$(grep -c . "Analysis_Results/Raw_Table_SLRNA_Loci/Raw_SLRNA_loci.txt")

  rm -r "Analysis_Results/Final_SLRNA_information"
  mkdir "Analysis_Results/Final_SLRNA_information"
  mkdir "Analysis_Results/Final_SLRNA_information/Temp"

  echo "Species Loci_ID Identification_Method Chromosome Strand Search_Start Search_End SL_TAG Unknown_Bases Repeats CD-HIT_Cluster TreeCluster Unique_SLRNA_ID Refined_Start Refined_End Hairpin_Structure" | tr " " "\t" >> "Analysis_Results/Final_SLRNA_information/Final_loci_Table.txt"

  while [ $count -le $recorrer ]
  do
    echo $count" // "$recorrer

    gru=$(sed -n $count"p" "Analysis_Results/Raw_Table_SLRNA_Loci/Raw_SLRNA_loci.txt" | awk -F "\t" '{print $1}')
    spe=$(sed -n $count"p" "Analysis_Results/Raw_Table_SLRNA_Loci/Raw_SLRNA_loci.txt" | awk -F "\t" '{print $2}')
    raw_cdhit_cluster=$(sed -n $count"p" "Analysis_Results/Raw_Table_SLRNA_Loci/Raw_SLRNA_loci.txt" | awk -F "\t" '{print $3}')
    chr=$(sed -n $count"p" "Analysis_Results/Raw_Table_SLRNA_Loci/Raw_SLRNA_loci.txt" | awk -F "\t" '{print $5}')
    strand=$(sed -n $count"p" "Analysis_Results/Raw_Table_SLRNA_Loci/Raw_SLRNA_loci.txt" | awk -F "\t" '{print $6}')
    start=$(sed -n $count"p" "Analysis_Results/Raw_Table_SLRNA_Loci/Raw_SLRNA_loci.txt" | awk -F "\t" '{print $7}')
    end=$(sed -n $count"p" "Analysis_Results/Raw_Table_SLRNA_Loci/Raw_SLRNA_loci.txt" | awk -F "\t" '{print $8}')
    total_n_bases=$(sed -n $count"p" "Analysis_Results/Raw_Table_SLRNA_Loci/Raw_SLRNA_loci.txt" | awk -F "\t" '{print $9}')
    repeats=$(sed -n $count"p" "Analysis_Results/Raw_Table_SLRNA_Loci/Raw_SLRNA_loci.txt" | awk -F "\t" '{print $10}')
    identification_method=$(sed -n $count"p" "Analysis_Results/Raw_Table_SLRNA_Loci/Raw_SLRNA_loci.txt" | awk -F "\t" '{print $11}')
    raw_filocluster=$(sed -n $count"p" "Analysis_Results/Raw_Table_SLRNA_Loci/Raw_SLRNA_loci.txt" | awk -F "\t" '{print $12}')
    raw_sequence_ID=$(sed -n $count"p" "Analysis_Results/Raw_Table_SLRNA_Loci/Raw_SLRNA_loci.txt" | awk -F "\t" '{print $13}')
    Definitive_loci=$(sed -n $count"p" "Analysis_Results/Raw_Table_SLRNA_Loci/Raw_SLRNA_loci.txt" | awk -F "\t" '{print $14}')

    if [ $gru == "Cestodes" ]
    then
      short_gru=Ces
    else
      short_gru=Tre
    fi

    replace=$(echo $short_gru"_cdhit_Clus-")
    cdhit_cluster=$(echo $raw_cdhit_cluster | sed "s/Cluster_/$replace/")

    filocluster_add=$(echo $short_gru"_Tree_clus")
    if [ ! "$raw_filocluster" == "NA" ]
    then
      if [ $raw_filocluster -ge 1 ]
      then
        filocluster=$(echo $filocluster_add"-"$raw_filocluster)
      else
        filocluster=NA
      fi
    else
      if [ ! "$cdhit_cluster" == "." ]
      then
        check_cluster=$(grep $gru "Analysis_Results/Raw_Table_SLRNA_Loci/Raw_SLRNA_loci.txt" | awk -F "\t" -v clus=$raw_cdhit_cluster '{if ($3==clus) print $12}' | grep -v -c NA)
        if [ $check_cluster -gt 0 ]
        then
          raw_filocluster=$(grep $gru "Analysis_Results/Raw_Table_SLRNA_Loci/Raw_SLRNA_loci.txt" | awk -F "\t" -v clus=$raw_cdhit_cluster '{if ($3==clus) print $12}' | grep -v NA)
          filocluster=$(echo $filocluster_add"-"$raw_filocluster)
        else
          filocluster=NA
        fi
      fi
    fi

    if [ ! "$identification_method" == "RBL" ]
    then
      cat "Analysis_Results/Alineamientos_pSL/"$gru"/Problematic_SLs.fasta" "Analysis_Results/Alineamientos_pSL/"$gru"/"$spe"_SLs_todos.fasta" "Analysis_Results/Alineamientos_pSL/"$gru"/"$spe"_SLs_Truncado.fasta" | seqkit grep -p $raw_sequence_ID  >> "Analysis_Results/Final_SLRNA_information/Temp/"$raw_sequence_ID"_raw.fasta"
      check_recovery_blast=$(grep -w -c $raw_sequence_ID $selected_on_recovery_blast)
      if [ $check_recovery_blast -gt 0 ]
      then
        identification_method=RBL
      fi
    else
      seqkit grep -p $raw_sequence_ID "Analysis_Results/Extra_SLs_Busqueda_BLAST/New_Sequences/All_New_Seq.fasta" >> "Analysis_Results/Final_SLRNA_information/Temp/"$raw_sequence_ID"_raw.fasta"
    fi

    check_SL_tag=$(seqkit seq --rna2dna "Analysis_Results/Final_SLRNA_information/Temp/"$raw_sequence_ID"_raw.fasta" | seqkit locate -f "Selected_Sequences_Work/"$gru"_SL_Tags.fasta" | grep -c .)
    if [ $check_SL_tag -gt 1 ]
    then
      SL_TAG=$(seqkit seq --rna2dna "Analysis_Results/Final_SLRNA_information/Temp/"$raw_sequence_ID"_raw.fasta" | seqkit locate -f "Selected_Sequences_Work/"$gru"_SL_Tags.fasta" | sed 1d | awk -F "\t" '{print $2}' | sed "s/ .*//" | tr "\n" ";" | sed "s/;$/\n/")
    else
      SL_TAG=NA
    fi

    check_selected_work=$(grep -F -w -c $raw_sequence_ID $selected_loci_work)
    if [ $check_selected_work -gt 0 ]
    then
      unique_ID=$(seqkit locate -f $unique_SLRNA "Analysis_Results/Final_SLRNA_information/Temp/"$raw_sequence_ID"_raw.fasta" | sed 1d | awk -F "\t" '{print $2}')
      unique_Structure=$(grep -w $unique_ID $unique_structure | awk -F "\t" '{print $2}')

      temp_start=$(seqkit locate -f $unique_SLRNA "Analysis_Results/Final_SLRNA_information/Temp/"$raw_sequence_ID"_raw.fasta" | sed 1d | awk -F "\t" '{print $5}')
      temp_end=$(seqkit locate -f $unique_SLRNA "Analysis_Results/Final_SLRNA_information/Temp/"$raw_sequence_ID"_raw.fasta" | sed 1d | awk -F "\t" '{print $6}')

      fixed_start=$(echo $start" "$end | awk -F " " -v correct=$temp_start '{print $1 + correct - 1}')
      fixed_end=$(echo $start" "$end | awk -F " " -v correct=$temp_end '{print $2 - (($2 - $1 + 1) - correct)}')
    else
      unique_ID=NA
      unique_Structure=NA

      fixed_start=NA
      fixed_end=NA
    fi

    echo $spe" "$Definitive_loci" "$identification_method" "$chr" "$strand" "$start" "$end" "$SL_TAG" "$total_n_bases" "$repeats" "$cdhit_cluster" "$filocluster" "$unique_ID" "$fixed_start" "$fixed_end" "$unique_Structure | tr " " "\t" >> "Analysis_Results/Final_SLRNA_information/Final_loci_Table.txt"
    count=$(($count + 1))
  done
fi

if [ "$FINAL_SLRNA_ALN" == "TRUE" ]
then
  rm "Analysis_Results/Final_SLRNA_information/Final_SLRNA_phy.fasta"
  rm "Analysis_Results/Final_SLRNA_information/Final_SLRNA_phy.aln"

  species=$(sed 1d "Analysis_Results/Final_SLRNA_information/Final_loci_Table.txt" | awk -F "\t" '{print $1}' | sort -u)

  for spe in $species
  do
    uniq_sl_rna=$(grep -w $spe "Analysis_Results/Final_SLRNA_information/Final_loci_Table.txt" |  awk -F "\t" '{print $13}' | sort -u)
    spe_count=1

    for uniq in $uniq_sl_rna
    do
      check_number=$(grep -w $spe "Analysis_Results/Final_SLRNA_information/Final_loci_Table.txt" | grep -w -c $uniq )

      if [ $check_number -gt 1 ]
      then
        new_id=$(echo ">"$spe"_SL_"$spe_count"__"$uniq"+")
      else
        new_id=$(echo ">"$spe"_SL_"$spe_count"__"$uniq)
      fi

      seqkit grep -p $uniq $unique_SLRNA | sed "s/>.*/$new_id/" >> "Analysis_Results/Final_SLRNA_information/Final_SLRNA_phy.fasta"

      spe_count=$(($spe_count + 1))
    done
  done

  species_referencia=$(awk -F "\t" '{print $3}' "Selected_Sequences_Work/Selected_work_reference_loci.txt" | sort -u)
  for spe in $species_referencia
  do
    uniq_sl_rna=$(grep -w $spe "Selected_Sequences_Work/Selected_work_reference_loci.txt" |  awk -F "\t" '{print $2}' | sort -u)
    for uniq in $uniq_sl_rna
    do
      check_number=$(grep -w $spe "Selected_Sequences_Work/Selected_work_reference_loci.txt" | grep -w -c $uniq )
      if [ $check_number -gt 1 ]
      then
        new_id=$(echo ">"$spe"_Reference_SL__"$uniq"+")
      else
        new_id=$(echo ">"$spe"_Reference_SL__"$uniq)
      fi
      seqkit grep -p $uniq $unique_SLRNA | sed "s/>.*/$new_id/" >> "Analysis_Results/Final_SLRNA_information/Final_SLRNA_phy.fasta"
    done
  done

  mafft --localpair --maxiterate 1000 --thread 20 "Analysis_Results/Final_SLRNA_information/Final_SLRNA_phy.fasta" | seqkit seq -u > "Analysis_Results/Final_SLRNA_information/Final_SLRNA_phy.aln"
fi
