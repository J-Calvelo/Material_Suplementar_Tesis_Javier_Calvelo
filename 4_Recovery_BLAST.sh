#!/usr/bin/env bash

# Once the initial candidates have been selected, this script conducts a new
# blast search, trying to identify more SLRNA sequences that managed to avoid
# detection so far by a full BLAST search with preselected and reference
# sequences.

# Other Variables
# Archivos SLs info
selected_sls_sec=Selected_Sequences_Work/Selected_SLRNA_preblast.fasta
selected_sls_ref_sec=Selected_Sequences_Work/Selected_SLRNA_referencia.fasta
SM_like_sites=Selected_Sequences_Work/SM_Patterns.fasta

# cutoffs
cov_cutoff=85 # aprox 13 missing bases
iden_cutoff=90 # aprox 9 mismatches
Rango_Region=500

## Variables corrida de modulos
RECOVERY_BLAST=FALSE
GENERATE_RAW_TABLE=TRUE

################################################################################

### Functions
extract_new_sls ()
{
  for nl in $New_Loci
  do
    rename=$(echo $especie"_"$nl)
    new_extract=$(grep -w $nl $store_blast"/"$especie"_loci.tab" | awk -F "\t" '{print $1" "$4-30"-"$5+30" "$6}')
    printf "%s %s %s %s\n%s %s %s\n" $new_extract | blastdbcmd -db $store_ref"/"$especie"_ref" -entry_batch - | sed "s/>.*/>$rename/" | tr -d " " | seqkit seq -w 0 --dna2rna  >> "Analysis_Results/Extra_SLs_Busqueda_BLAST/New_Sequences/"$especie"_"$grupo".fasta"
  done
}

extract_sucio_sls ()
{
  for nl in $New_Loci
  do
    rename=$(echo $especie"_"$nl)
    get_chr=$(grep -w $especie "Analysis_Results/pSL_Overlapping_Features/Matching_Info_Overlap.tab" | grep -w $nl | awk -F "\t" '{print $4}' | sort -u)
    get_coord=$(grep -w $especie "Analysis_Results/pSL_Overlapping_Features/Matching_Info_Overlap.tab" | grep -w $nl | awk -F "\t" '{print $5}' | sort -u | awk -F "_" '{print $1"-"$2}')
    strand=$(grep -w $nl $store_blast"/"$especie"_loci.tab" | awk -F "\t" '{print $6}')
    echo $nl": "$get_chr" "$get_coord" "$strand

    printf "%s %s %s %s\n%s %s %s\n" $get_chr" "$get_coord" "$strand | blastdbcmd -db $store_ref"/"$especie"_ref" -entry_batch - | sed "s/>.*/>$rename/" | tr -d " " | seqkit seq -w 0 --dna2rna  >> "Analysis_Results/Extra_SLs_Busqueda_BLAST/New_Sequences/"$especie"_"$grupo".fasta"
  done
}

################################################################################
################################################################################

# short cuts
store_ref=Analysis_Results/Extra_SLs_Busqueda_BLAST/Referencias_BLAST
store_blast=Analysis_Results/Extra_SLs_Busqueda_BLAST/BLAST_Files
store_temp=Analysis_Results/Extra_SLs_Busqueda_BLAST/Temporales


if [ "$RECOVERY_BLAST" == "TRUE" ]
then
  rm -r "Analysis_Results/Extra_SLs_Busqueda_BLAST"
  mkdir "Analysis_Results/Extra_SLs_Busqueda_BLAST"
  mkdir "Analysis_Results/Extra_SLs_Busqueda_BLAST/New_Sequences"
  mkdir $store_ref
  mkdir $store_blast
  mkdir $store_temp

  cat $selected_sls_sec $selected_sls_ref_sec >> "Analysis_Results/Extra_SLs_Busqueda_BLAST/SLs_For_search.fasta"

  especies_genome=$(grep -w -v Run_Status "Analysis_Results/Summary_SLFinder.tab" | awk -F "\t" '{print $1"\t"$15}' | grep -w -v NA | awk -F "\t" '{print $1}')

  echo "Especie Nuevos_Hits RepMasker_Base RepMasker_Custom" | tr " " "\t" >> "Analysis_Results/Extra_SLs_Busqueda_BLAST/Summary.tab"

  for especie in $especies_genome
  do
    new_loci_counter=1
    genome_file_blast=$(grep -w $especie "Analysis_Results/Summary_SLFinder.tab" | awk -F "\t" '{print $15}')
    genome_file=$(grep -w $especie "Analysis_Results/Summary_SLFinder.tab" | awk -F "\t" '{print $15}')
    grupo=$(grep -w $especie "Analysis_Results/Summary_SLFinder.tab" | awk -F "\t" '{print $2}')

    echo "Especie: "$especie" Archivo: "$genome_file
    makeblastdb -in "Analysis_Results/Reference_Genome/"$genome_file_blast -dbtype nucl -parse_seqids -input_type fasta -out $store_ref"/"$especie"_ref"
    echo "BLAST"

    blastn -query "Analysis_Results/Extra_SLs_Busqueda_BLAST/SLs_For_search.fasta" -db $store_ref"/"$especie"_ref" -outfmt "6 std nident qlen" -out $store_blast"/"$especie"_raw.blastn" -qcov_hsp_perc $cov_cutoff -perc_identity $iden_cutoff
    echo "Filtrando"
    echo "Cromosoma Ident Matches Start End Strand SL_Match Hit_num Loci RepMasker_Base RepMasker_Custom" | tr " " "\t" >> $store_blast"/"$especie"_loci.tab"

    #Chr Iden Qcov Start End Strand SL_ID
    awk -F "\t" '{if ($10>$9) print $2"\t"$3"\t"$13"\t"$9"\t"$10"\tplus\t"$1; else print $2"\t"$3"\t"$13"\t"$10"\t"$9"\tminus\t"$1}' $store_blast"/"$especie"_raw.blastn" | sort -u >> $store_temp"/"$especie"_temp.blastn"

    chromosomes=$(awk -F "\t" '{print $1}' $store_temp"/"$especie"_temp.blastn" | sort -u)
    loci_num=1

    for chr in $chromosomes
    do
      # Getting info locis coord
      grep -w $especie "Analysis_Results/pSL_Overlapping_Features/Matching_Info_Overlap.tab" | grep -w $chr | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$5}' | sort -u > $store_temp"/"$especie"_"$chr"_coord.temp"
      grep -w $especie "Analysis_Results/pSL_Overlapping_Features/Matching_Info_Region.tab" | grep -w $chr | grep RepeatMasker_Basic | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$9}' | grep -v No_Data |  grep -v Unknown | sort -u >> $store_temp"/"$especie"_"$chr"_transposon_basic.temp"
      grep -w $especie "Analysis_Results/pSL_Overlapping_Features/Matching_Info_Region.tab" | grep -w $chr | grep RepeatMasker_Custom | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$9}' | grep -v No_Data |  grep -v Unknown | sort -u >> $store_temp"/"$especie"_"$chr"_transposon_custom.temp"

      total_coords=$(grep -c . $store_temp"/"$especie"_"$chr"_coord.temp")

      repeat_masker_file1=$(echo $store_temp"/"$especie"_"$chr"_full_transposon_base.temp")
      repeat_masker_file2=$(echo $store_temp"/"$especie"_"$chr"_full_transposon_custom.temp")

      grep -w $chr "Analysis_Results/Info_Transposones/RepeatMasker_Runs/"$especie"/"$genome_file".work" | grep -v Simple_repeat >> $repeat_masker_file1
      grep -w $chr "Analysis_Results/Info_Transposones/RepeatMasker_Custom_Runs/"$especie"/"$genome_file".work" | grep -v Simple_repeat >> $repeat_masker_file2

      grep -w $chr $store_temp"/"$especie"_temp.blastn" | sort -n -r -k 4 >> $store_temp"/"$especie"_"$chr"_temp2.blastn"

      start_hit=$(sed -n 1p $store_temp"/"$especie"_"$chr"_temp2.blastn" | awk -F "\t" '{print $4}')
      end_hit=$(sed -n 1p $store_temp"/"$especie"_"$chr"_temp2.blastn" | awk -F "\t" '{print $5}')

      while_stop=$(grep -c . $store_temp"/"$especie"_"$chr"_temp2.blastn")
      count=2

      sed -n 1p $store_temp"/"$especie"_"$chr"_temp2.blastn" | awk -F "\t" -v loci=$loci_num '{print $0"\tWork_"loci}' >> $store_temp"/"$especie"_"$chr"_temp3.blastn"

      while [ $count -le $while_stop ]
      do
        test_start=$(sed -n $count"p" $store_temp"/"$especie"_"$chr"_temp2.blastn" | awk -F "\t" '{print $4}')
        test_end=$(sed -n $count"p" $store_temp"/"$especie"_"$chr"_temp2.blastn" | awk -F "\t" '{print $5}')

        if [ $test_start -le $end_hit -a $start_hit -le $test_end ]
        then
          sed -n $count"p" $store_temp"/"$especie"_"$chr"_temp2.blastn" | awk -F "\t" -v loci=$loci_num '{print $0"\tWork_"loci}' >> $store_temp"/"$especie"_"$chr"_temp3.blastn"
        else
          loci_num=$(($loci_num + 1))
          sed -n $count"p" $store_temp"/"$especie"_"$chr"_temp2.blastn" | awk -F "\t" -v loci=$loci_num '{print $0"\tWork_"loci}' >> $store_temp"/"$especie"_"$chr"_temp3.blastn"

          start_hit=$test_start
          end_hit=$test_end
        fi
        count=$(($count + 1))
      done
      loci_num=$(($loci_num + 1))

      # Determinando_Hits_previos:
      loci_ids=$(awk -F "\t" '{print $8}' $store_temp"/"$especie"_"$chr"_temp3.blastn" | sort -u)
      for locus in $loci_ids
      do
        best_coverage=$(grep $locus$ $store_temp"/"$especie"_"$chr"_temp3.blastn" | awk -F "\t" '{print $3}' | sort -n | tail -n 1)
        best_ident=$(grep $locus$ $store_temp"/"$especie"_"$chr"_temp3.blastn" | awk -F "\t" -v cov=$best_coverage '{if ($3==cov) print $2*1000}' | sort -n | tail -n 1 | awk '{print $1/1000}')

        start_best_loci=$(grep $locus$ $store_temp"/"$especie"_"$chr"_temp3.blastn" | awk -F "\t" -v cov=$best_coverage -v iden=$best_ident '{if ($3==cov && $2==iden) print}' | head -n 1 | awk -F "\t" '{print $4}')
        end_best_loci=$(grep $locus$ $store_temp"/"$especie"_"$chr"_temp3.blastn" | awk -F "\t" -v cov=$best_coverage -v iden=$best_ident '{if ($3==cov && $2==iden) print}' | head -n 1 | awk -F "\t" '{print $5}')

        print_info=$(grep $locus$ $store_temp"/"$especie"_"$chr"_temp3.blastn" | awk -F "\t" -v cov=$best_coverage -v iden=$best_ident '{if ($3==cov && $2==iden) print}' | head -n 1)

        count_loci=1

        Known_Hit=FALSE
        while [ $count_loci -le $total_coords ]
        do
          loci_start=$(sed -n $count_loci"p" $store_temp"/"$especie"_"$chr"_coord.temp" | awk -F "\t" '{print $4}' | awk -F "_" '{print $1}')
          loci_end=$(sed -n $count_loci"p" $store_temp"/"$especie"_"$chr"_coord.temp" | awk -F "\t" '{print $4}' | awk -F "_" '{print $2}')

          if [ $start_best_loci -le $loci_end -a $loci_start -le $end_best_loci ]
          then
            Known_Hit=TRUE
            loci_ID=$(sed -n $count_loci"p" $store_temp"/"$especie"_"$chr"_coord.temp" | awk -F "\t" '{print $3}')
            check_transposon_1=$(grep -c -w $loci_ID $store_temp"/"$especie"_"$chr"_transposon_basic.temp")
            if [ $check_transposon_1 -gt 0 ]
            then
              transposon_info1=$(grep -w $loci_ID $store_temp"/"$especie"_"$chr"_transposon_basic.temp" | awk -F ";" '{print $2}' | tr "\n" ";" | sed "s/;$//")
            else
              transposon_info1=NA
            fi

            check_transposon_2=$(grep -c -w $loci_ID $store_temp"/"$especie"_"$chr"_transposon_custom.temp")
            if [ $check_transposon_2 -gt 0 ]
            then
              transposon_info2=$(grep -w $loci_ID $store_temp"/"$especie"_"$chr"_transposon_custom.temp" | awk -F ";" '{print $2}' | tr "\n" ";" | sed "s/;$//")
            else
              transposon_info2=NA
            fi
            count=$(($total_coords + $total_coords))
          fi
          count_loci=$(($count_loci + 1))
        done

        if [ "$Known_Hit" == "FALSE" ]
        then
          loci_ID=$(echo "NEW_"$new_loci_counter)
          transposon_info1=$(awk -F "\t" -v start=$loci_start -v end=$loci_end '{if (start <= $7 && $6 <= end) print $10";"$11}' $repeat_masker_file1 | sort -u | tr "\n" "|" | sed "s/|/__/g")
          transposon_info2=$(awk -F "\t" -v start=$loci_start -v end=$loci_end '{if (start <= $7 && $6 <= end) print $10";"$11}' $repeat_masker_file2 | sort -u | tr "\n" "|" | sed "s/|/__/g")

          new_loci_counter=$(($new_loci_counter + 1))
        fi
        echo $print_info" "$loci_ID" "$transposon_info1" "$transposon_info2 | tr " " "\t" >> $store_blast"/"$especie"_loci.tab"
      done
      rm $store_temp"/"$especie"_"$chr"_temp2.blastn"
      rm $store_temp"/"$especie"_"$chr"_temp3.blastn"
    done

    # Get New SL candidates
    check_New=$(awk -F "\t" '{print $9}' $store_blast"/"$especie"_loci.tab" | grep -v RepMasker_Base | grep -c NEW_)
    if [ $check_New -gt 0 ]
    then
      New_Loci=$(awk -F "\t" '{print $9}' $store_blast"/"$especie"_loci.tab" | grep NEW_)
      grupo=Nuevos
      extract_new_sls
    fi

    check_dirty1=$(awk -F "\t" '{print $9"\t"$10}' $store_blast"/"$especie"_loci.tab" | grep -v RepMasker_Base | grep -v NEW | grep -c -w -v NA)
    if [ $check_dirty1 -gt 0 ]
    then
      New_Loci=$(awk -F "\t" '{print $9"\t"$10}' $store_blast"/"$especie"_loci.tab" | grep -v RepMasker_Base | grep -v NEW | grep -w -v NA | awk -F "\t" '{print $1}')
      grupo=Transposon_Base
      extract_sucio_sls
    fi

    check_dirty2=$(awk -F "\t" '{print $9"\t"$11}' $store_blast"/"$especie"_loci.tab" | grep -v RepMasker_ | grep -v NEW | grep -c -w -v NA)
    if [ $check_dirty2 -gt 0 ]
    then
      New_Loci=$(awk -F "\t" '{print $9"\t"$11}' $store_blast"/"$especie"_loci.tab" | grep -v RepMasker_ | grep -v NEW | grep -w -v NA | awk -F "\t" '{print $1}')
      grupo=Transposon_Custom
      extract_sucio_sls
    fi

    check_dirty3=$(awk -F "\t" '{print $9"\t"$10"_"$11}' $store_blast"/"$especie"_loci.tab" | grep -v NEW | grep -w -c NA_NA)
    if [ $check_dirty3 -gt 0 ]
    then
      New_Loci=$(awk -F "\t" '{print $9"\t"$10"_"$11}' $store_blast"/"$especie"_loci.tab" | grep -v NEW | grep -w NA_NA | awk -F "\t" '{print $1}')

      grupo=Asociated
      extract_sucio_sls
    fi

    echo $especie" "$check_New" "$check_dirty1" "$check_dirty2 | tr " " "\t" >> "Analysis_Results/Extra_SLs_Busqueda_BLAST/Summary.tab"
  done

  cat "Analysis_Results/Extra_SLs_Busqueda_BLAST/New_Sequences/"*"_Nuevos.fasta" >> "Analysis_Results/Extra_SLs_Busqueda_BLAST/New_Sequences/All_New_Seq.fasta"
  cat "Analysis_Results/Extra_SLs_Busqueda_BLAST/New_Sequences/"*"_Asociated.fasta" >> "Analysis_Results/Extra_SLs_Busqueda_BLAST/New_Sequences/All_Ass_Seq.fasta"

  seqkit seq -n "Analysis_Results/Extra_SLs_Busqueda_BLAST/New_Sequences/All_Ass_Seq.fasta" >> "Analysis_Results/Extra_SLs_Busqueda_BLAST/New_Sequences/All_Ass_Seq.ids"
  seqkit locate -m 2 -P -f $SM_like_sites "Analysis_Results/Extra_SLs_Busqueda_BLAST/New_Sequences/All_New_Seq.fasta" >> "Analysis_Results/Extra_SLs_Busqueda_BLAST/New_Sequences/All_New_Seq_sm.coord"
  seqkit locate -m 2 -P -f $SM_like_sites "Analysis_Results/Extra_SLs_Busqueda_BLAST/New_Sequences/All_Ass_Seq.fasta" >> "Analysis_Results/Extra_SLs_Busqueda_BLAST/New_Sequences/All_Ass_Seq_sm.coord"
fi

################################################################################

if [ "$GENERATE_RAW_TABLE" == "TRUE" ]
then
  grupos=$(echo "Cestodes Trematoda")

  rm -r "Analysis_Results/Raw_Table_SLRNA_Loci"
  mkdir "Analysis_Results/Raw_Table_SLRNA_Loci"
  mkdir "Analysis_Results/Raw_Table_SLRNA_Loci/Temp"

  echo "Group Species Cluster_CDHIT Raw_Loci Chromosome Strand Scan_start Scan_End Total_N_Bases Repeats Identification Filo_Cluster Raw_loci Definitive_loci" | tr " " "\t" > "Analysis_Results/Raw_Table_SLRNA_Loci/Raw_SLRNA_loci.txt"

  for gru in $grupos
  do
    species=$(ls -d $gru"/"*"/" | awk -F "/" '{print $2}')

    echo $gru": "$species
    for spe in $species
    do
      genome_file=$(grep -w $spe "Analysis_Results/Summary_SLFinder.tab" | awk -F "\t" '{print $15}')
      chromosomes=$(seqkit seq -n -i $gru"/"$spe"/"$genome_file )
      Definitive_loci_add=0

      for chr in $chromosomes
      do
        check_standard_search=$(grep -c -F -w $chr "Analysis_Results/Alineamientos_pSL/"$gru"/"$spe"_SLs.coords")
        if [ $check_standard_search -gt 0 ]
        then
          awk -F "\t" -v chr=$chr '{if ($4==chr) print}' "Analysis_Results/Alineamientos_pSL/"$gru"/"$spe"_SLs.coords" > "Analysis_Results/Raw_Table_SLRNA_Loci/Temp/standard_loci.tmp"

          count=1
          recorrer=$(grep -c . "Analysis_Results/Raw_Table_SLRNA_Loci/Temp/standard_loci.tmp")

          while [ $count -le $recorrer ]
          do
            raw_loci=$(sed -n $count"p" "Analysis_Results/Raw_Table_SLRNA_Loci/Temp/standard_loci.tmp" | awk -F "\t" '{print $2}')
            scan_start=$(sed -n $count"p" "Analysis_Results/Raw_Table_SLRNA_Loci/Temp/standard_loci.tmp" | awk -F "\t" '{print $5}' | awk -F "-" '{print $1}')
            scan_end=$(sed -n $count"p" "Analysis_Results/Raw_Table_SLRNA_Loci/Temp/standard_loci.tmp" | awk -F "\t" '{print $5}' | awk -F "-" '{print $2}')
            strand=$(sed -n $count"p" "Analysis_Results/Raw_Table_SLRNA_Loci/Temp/standard_loci.tmp" | awk -F "\t" '{print $6}')

            echo "Original Loci: "$count" // "$recorrer" || "$spe" "$raw_loci

            clusterID=$(grep -w $spe "Analysis_Results/pSL_Overlapping_Features/Matching_Info_Region.tab" | awk -F "\t" -v loci=$raw_loci '{if (loci==$3) print $6}' | sort -u)
            check_repeats=$(grep -w $spe "Analysis_Results/pSL_Overlapping_Features/Matching_Info_Region.tab" | grep "RepeatMasker_" | awk -F "\t" -v loci=$raw_loci '{if (loci==$3) print $9}' | grep -v "No_Data" | grep -v "Unknown" | grep -c .)
            if [ $check_repeats -gt 0 ]
            then
              repeats=$(grep -w $spe "Analysis_Results/pSL_Overlapping_Features/Matching_Info_Region.tab" | grep "RepeatMasker_" | awk -F "\t" -v loci=$raw_loci '{if (loci==$3) print $9}' | awk -F ";" '{print $2}' | sort -u | tr "\n" ";" | sed "s/;$/\n/" | sed "s/^;//")
            else
              repeats=NA
            fi

            check_detection_both=$(cat "Analysis_Results/Alineamientos_pSL/"$gru"/Problematic_SLs.fasta" "Analysis_Results/Alineamientos_pSL/"$gru"/"$spe"_SLs_todos.fasta" "Analysis_Results/Alineamientos_pSL/"$gru"/"$spe"_SLs_Truncado.fasta" | grep "_"$spe"_" | grep BLAST_SLF_ | grep -c $raw_loci$)
            check_detection_SLF=$(cat "Analysis_Results/Alineamientos_pSL/"$gru"/Problematic_SLs.fasta" "Analysis_Results/Alineamientos_pSL/"$gru"/"$spe"_SLs_todos.fasta" "Analysis_Results/Alineamientos_pSL/"$gru"/"$spe"_SLs_Truncado.fasta"  | grep "_"$spe"_"  | grep SLFinder_Loci_ | grep -c $raw_loci$)

            if [ $check_detection_both -gt 0 ]
            then
              identification_method=SLF+BL
              raw_sequence_ID=$(cat "Analysis_Results/Alineamientos_pSL/"$gru"/Problematic_SLs.fasta" "Analysis_Results/Alineamientos_pSL/"$gru"/"$spe"_SLs_todos.fasta" "Analysis_Results/Alineamientos_pSL/"$gru"/"$spe"_SLs_Truncado.fasta" | grep "_"$spe"_" | grep BLAST_SLF_ | grep $raw_loci$ | sed "s/>//")
            elif [ $check_detection_SLF -gt 0 ]
            then
              identification_method=SLF
              raw_sequence_ID=$(cat "Analysis_Results/Alineamientos_pSL/"$gru"/Problematic_SLs.fasta" "Analysis_Results/Alineamientos_pSL/"$gru"/"$spe"_SLs_todos.fasta" "Analysis_Results/Alineamientos_pSL/"$gru"/"$spe"_SLs_Truncado.fasta" | grep "_"$spe"_" | grep SLFinder_Loci_ | grep $raw_loci$ | sed "s/>//")
            else
              identification_method=BL
              raw_sequence_ID=$(cat "Analysis_Results/Alineamientos_pSL/"$gru"/Problematic_SLs.fasta" "Analysis_Results/Alineamientos_pSL/"$gru"/"$spe"_SLs_todos.fasta" "Analysis_Results/Alineamientos_pSL/"$gru"/"$spe"_SLs_Truncado.fasta" | grep "_"$spe"_" | grep BLAST_New_ | grep "_"$raw_loci"_" | sed "s/>//")
            fi

            total_N=$(cat "Analysis_Results/Alineamientos_pSL/"$gru"/Problematic_SLs.fasta" "Analysis_Results/Alineamientos_pSL/"$gru"/"$spe"_SLs_todos.fasta" "Analysis_Results/Alineamientos_pSL/"$gru"/"$spe"_SLs_Truncado.fasta" | seqkit grep -p $raw_sequence_ID | seqkit seq -s | grep -o N | grep -c .)

            # Check truncated loci
            check_truncated=$(cat "Analysis_Results/Alineamientos_pSL/"$gru"/Problematic_SLs.fasta" "Analysis_Results/Alineamientos_pSL/"$gru"/"$spe"_SLs_Truncado.fasta" | grep -w -c $raw_sequence_ID)
            if [ $check_truncated -gt 0 ]
            then
              identification_method=$(echo $identification_method"*")
            fi

            # Check the exact ID used for Treeclus
            ######################################################################
            check_filocluster=$(grep -c -F $raw_sequence_ID "Analysis_Results/MEME_Profile/Filogeny_Clustering/General/"$gru"_filo_sec_Region.fasta")
            if [ $check_filocluster -eq 0 ]
            then
              filocluster_ID=NA
            elif [ $check_filocluster -gt 1 ]
            then
              check_noise=$(grep -c -F $raw_sequence_ID"_" "Analysis_Results/MEME_Profile/Filogeny_Clustering/General/"$gru"_filo_sec_Region.fasta")
              if [ $check_noise -gt 0 ]
              then
                filocluster_ID=$(grep -F $raw_sequence_ID"_" "Analysis_Results/MEME_Profile/Filogeny_Clustering/General/"$gru"_filo_sec_Region.fasta" | sed "s/>//")
              else
                filocluster_ID=$(grep $raw_sequence_ID$ "Analysis_Results/MEME_Profile/Filogeny_Clustering/General/"$gru"_filo_sec_Region.fasta" | sed "s/>//")
              fi
            else
              filocluster_ID=$(grep -F $raw_sequence_ID "Analysis_Results/MEME_Profile/Filogeny_Clustering/General/"$gru"_filo_sec_Region.fasta" | sed "s/>//")
            fi

            if [ $check_filocluster -eq 0 ]
            then
              filo_cluster=NA
            else
              filo_cluster=$(grep -w -F $filocluster_ID "Analysis_Results/MEME_Profile/Filogeny_Clustering/General/"$gru"_filo_sec_Region.clus" | awk -F "\t" '{print $2}')
            fi
            echo $gru" "$spe" "$clusterID" "$raw_loci" "$chr" "$strand" "$scan_start" "$scan_end" "$total_N" "$repeats" "$identification_method" "$filo_cluster" "$raw_sequence_ID | tr " " "\t" >> "Analysis_Results/Raw_Table_SLRNA_Loci/Temp/current_loci_temp.txt"
            count=$(($count + 1))
          done
        fi

        check_last_blast_search=$(grep -F -w $chr "Analysis_Results/Extra_SLs_Busqueda_BLAST/BLAST_Files/"$spe"_loci.tab" | grep -c "NEW_")
        if [ $check_last_blast_search -gt 0 ]
        then
          awk -F "\t" -v chr=$chr '{if ($1==chr) print}' "Analysis_Results/Extra_SLs_Busqueda_BLAST/BLAST_Files/"$spe"_loci.tab" | grep "NEW_" > "Analysis_Results/Raw_Table_SLRNA_Loci/Temp/new_loci.tmp"
          count=1
          recorrer=$(grep -c . "Analysis_Results/Raw_Table_SLRNA_Loci/Temp/new_loci.tmp")

          while [ $count -le $recorrer ]
          do
            echo "New Loci: "$count" // "$recorrer
            raw_loci=$(sed -n $count"p" "Analysis_Results/Raw_Table_SLRNA_Loci/Temp/new_loci.tmp" | awk -F "\t" '{print $2}')
            scan_start=$(sed -n $count"p" "Analysis_Results/Raw_Table_SLRNA_Loci/Temp/new_loci.tmp" | awk -F "\t" '{print $5}' | awk -F "-" '{print $1}')
            scan_end=$(sed -n $count"p" "Analysis_Results/Raw_Table_SLRNA_Loci/Temp/new_loci.tmp" | awk -F "\t" '{print $5}' | awk -F "-" '{print $2}')
            strand=$(sed -n $count"p" "Analysis_Results/Raw_Table_SLRNA_Loci/Temp/new_loci.tmp" | awk -F "\t" '{print $6}')
            raw_sequence_ID=$(sed -n $count"p" "Analysis_Results/Raw_Table_SLRNA_Loci/Temp/new_loci.tmp" | awk -F "\t" -v spe=$spe '{print spe"_"$9}')

            total_N=$(seqkit grep -p $raw_sequence_ID "Analysis_Results/Extra_SLs_Busqueda_BLAST/New_Sequences/All_New_Seq.fasta" | seqkit seq -s | grep -o N | grep -c .)

            clusterID=NA
            identification_method=RBL
            filo_cluster=NA

            # Get repeats near new SLRNA
            repeat_range_start=$(($scan_start - $Rango_Region))
            repeat_range_end=$(($scan_end + $Rango_Region))

            repeat_basic=$(grep -w $chr "Analysis_Results/Info_Transposones/RepeatMasker_Runs/"$spe"/"$genome_file".work" | grep -v Simple_repeat | awk -F "\t" -v start=$repeat_range_start -v end=$repeat_range_end '{if (start <= $7 && $6 <= end) print $11}' | sort -u)
            repeat_custom=$(grep -w $chr "Analysis_Results/Info_Transposones/RepeatMasker_Custom_Runs/"$spe"/"$genome_file".work" | grep -v Simple_repeat | awk -F "\t" -v start=$repeat_range_start -v end=$repeat_range_end '{if (start <= $7 && $6 <= end) print $11}' | sort -u)

            check_repeats=$(echo $repeat_basic" "$repeat_custom | tr " " "\n" | grep -c .)
            if [ $check_repeats -eq 0 ]
            then
              repeats=NA
            else
              repeats=$(echo $repeat_basic" "$repeat_custom | tr " " "\n" | sort -u | tr "\n" ";" | sed "s/;$/\n/")
            fi
            echo $gru" "$spe" "$clusterID" "$raw_loci" "$chr" "$strand" "$scan_start" "$scan_end" "$total_N" "$repeats" "$identification_method" "$filo_cluster" "$raw_sequence_ID | tr " " "\t" >> "Analysis_Results/Raw_Table_SLRNA_Loci/Temp/current_loci_temp.txt"
            count=$(($count + 1))
          done
        fi
        if [ -f "Analysis_Results/Raw_Table_SLRNA_Loci/Temp/current_loci_temp.txt" ]
        then
          sort -n -k 7 "Analysis_Results/Raw_Table_SLRNA_Loci/Temp/current_loci_temp.txt" | awk -F "\t" -v add=$Definitive_loci_add '{print $0"\tDef_Loci-"NR+add}' >> "Analysis_Results/Raw_Table_SLRNA_Loci/Raw_SLRNA_loci.txt"
          Definitive_loci_add=$(awk -F "\t" -v add=$Definitive_loci_add '{print NR+add}' "Analysis_Results/Raw_Table_SLRNA_Loci/Temp/current_loci_temp.txt" | tail -n 1)
          rm "Analysis_Results/Raw_Table_SLRNA_Loci/Temp/current_loci_temp.txt"
        fi
      done
    done
  done
fi
