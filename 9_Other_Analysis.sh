#!/usr/bin/env bash

work_directory=Analysis_Results/SL_Acceptor_Sites
Result_Storage=Defenitive_Orthogroups
sitios_interes=/home/amanda/PROYECTS/SLFinder_Publication/Resultados_chi_Cuadrado
operones_interes=/home/amanda/PROYECTS/SLFinder_Publication/Interest_Operon_Pairs.in
interesting_extra_annotation=get_sl_counts_hog.txt

# Module variables
RUN_ACCEPTOR_INDIVIDUAL_SAMPLES_FREC=TRUE
RUN_ACCEPTOR_PHO_EXPRESION=TRUE
INTEREST_OPERONS=TRUE
INTEREST_OPERONS_CONTROL=TRUE
RUN_INTEREST_ANNOTATION=TRUE

# Functions
extract_info_interpro () {
  echo "HOG Species Gen InterprotID Description" | tr " " "\t" >> $work_directory"/"$Result_Storage"/Interprot_Results/Interest_PHO/"$group_id"_pho_annotation.txt"
  for pho in $selected_PHOs
  do
    grep -w -F $pho $PHO_File | tr "\t" "\n" | tr "," "\n" | grep . | tr -d " " | sed 1,3d > $work_directory"/"$Result_Storage"/Interprot_Results/Interest_PHO/Temp/genes_"$pho".tmp"
    for spe in $Especies
    do
      echo $spe": "$pho
      awk -F "\t" -v pho=$pho -v spe=$spe '{print pho";"spe";"$1";"$12";"$13}' $work_directory"/"$Result_Storage"/Interprot_Results/"$spe".tsv" | grep -w -F -f $work_directory"/"$Result_Storage"/Interprot_Results/Interest_PHO/Temp/genes_"$pho".tmp" | grep ";IPR" >> $work_directory"/"$Result_Storage"/Interprot_Results/Interest_PHO/Temp/genes_current_interpro.tmp"
    done
    sort -u $work_directory"/"$Result_Storage"/Interprot_Results/Interest_PHO/Temp/genes_current_interpro.tmp" | grep . | tr ";" "\t" >> $work_directory"/"$Result_Storage"/Interprot_Results/Interest_PHO/"$group_id"_pho_annotation.txt"
    rm $work_directory"/"$Result_Storage"/Interprot_Results/Interest_PHO/Temp/genes_current_interpro.tmp"
  done
}

##################################################################################################################################################
##################################################################################################################################################
###################################################################################################################################################

if [ "$RUN_ACCEPTOR_PHO_EXPRESION" == "TRUE" ]
then
  rm -r $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Sitios_Interes"
  mkdir $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Sitios_Interes"

  # Hacer los headers
  Especies=$(awk -F "\t" '{print $2}' $sitios_interes"/Sitios_significativos_todos.txt" | sort -u | grep -v Filtrado.Aceptor | tr -d '"')
  for spe in $Especies
  do
    echo $spe" header"
    count_data_file=$(ls $work_directory"/"$Result_Storage"/Mediciones_Expresion/Conteos/"*$spe"_All_htseq.counts")
    header_chisquare=$(sed -n 1p  $sitios_interes"/"$spe"_chisquare.txt" | tr -d "\r" | tr -d '"' | tr "\t" ";")
    header_expresion=$(sed -n 1p  $count_data_file | tr "\t" "\n" | sed 1d | tr "\n" ";" | sed "s/;$//")

    echo "Gen;Transcript;"$header_chisquare";Operon;Chimera;"$header_expresion | tr ";" "\t" >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Sitios_Interes/"$spe"_extra_info_chisquare.txt"
  done

  count=2
  recorrer=$(grep -c . $sitios_interes"/Sitios_significativos_todos.txt")

  while [ $count -le $recorrer ]
  do
    spe=$(sed -n $count"p" $sitios_interes"/Sitios_significativos_todos.txt" | awk -F "\t" '{print $2}' | tr -d '"')
    count_data_file=$(ls $work_directory"/"$Result_Storage"/Mediciones_Expresion/Conteos/"*$spe"_All_htseq.counts")
    site_info=$(sed -n $count"p" $sitios_interes"/Sitios_significativos_todos.txt" | awk -F "\t" '{print $3}' | tr -d '"')

    full_line=$(grep -w -F $site_info $sitios_interes"/"$spe"_chisquare.txt" | tr -d "\r" | sed 's/.*"\t"/"/' | tr "\t" ";")

    ace=$(echo $site_info | awk -F "__" '{print $1}')
    chr=$(echo $site_info | awk -F "__" '{print $2}')
    strand=$(echo $site_info | awk -F "__" '{print $3}')

    original_gen_id=$(awk -F "\t" -v ace=$ace -v chr=$chr -v strand=$strand '{if ($2==ace && $3==chr && $5=="Concordante" && $6!="Trans_End") print $1}' $work_directory"/SL_Read_Counts/"$spe"_tier3.counts" | sed "s/^gene-//")

    echo $spe" "$original_gen_id" "$site_info

    is_operon=$(grep -w -F $original_gen_id $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $10}' | grep -c Operon | awk '{if ($1>0) print "X"; else print "."}')
    is_chimera=$(grep -w -F $original_gen_id $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $7}' | grep -c Chimeric | awk '{if ($1>0) print "X"; else print "."}')

    check_trans=$(grep -w -F -c $original_gen_id $work_directory"/"$Result_Storage"/Mediciones_Expresion/Background_TPM/"$spe"_per_gene.txt")
    if [ $check_trans -gt 0 ]
    then
      transcripts=$(grep -w -F $original_gen_id  $work_directory"/"$Result_Storage"/Mediciones_Expresion/Background_TPM/"$spe"_per_gene.txt" | awk -F "\t" '{print $1}')
      for trans in $transcripts
      do
        count_data=$(grep -w -F $trans $count_data_file | tr "\t" "\n" | sed 1d | tr "\n" ";" | sed "s/;$//")
        echo $original_gen_id";"$trans";"$full_line";"$is_operon";"$is_chimera";"$count_data | tr ";" "\t" >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Sitios_Interes/"$spe"_extra_info_chisquare.txt"
      done
    else
      count_data=$(sed -n 1p  $count_data_file | tr "\t" "\n" | sed 1d | sed "s/.*/NA/")
      echo $original_gen_id";"$trans";"$full_line";"$is_operon";"$is_chimera";"$count_data | tr ";" "\t" >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Sitios_Interes/"$spe"_extra_info_chisquare.txt"
    fi

    count=$(($count + 1))
  done
fi

##################################################################################################################################################
##################################################################################################################################################
###################################################################################################################################################

if [ "$RUN_ACCEPTOR_INDIVIDUAL_SAMPLES_FREC" == "TRUE" ]
then
  rm -r $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Frecuencias_per_Sample"
  mkdir $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Frecuencias_per_Sample"

  Especies=$(ls $work_directory"/"$Result_Storage"/"*pep | sed "s/.*\///" | sed "s/.pep//")

  for spe in $Especies
  do
    if [ -f $sitios_interes"/"$spe"_chisquare.txt" ]
    then
      samples=$(head -n 1 $work_directory"/SL_Read_Counts/"$spe"_tier3.counts" | tr "\t" "\n" | sed 1,7d)
      sl_tags=$(head -n 1 $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Count_Table_Per_Tag/"$spe"_SL_Work_Counts.tab" | tr "\t" "\n" | sed 1d)

      site_info=$(cat $sitios_interes"/"$spe"_chisquare.txt" | tr -d "\r" | tr -d '"' | awk -F "\t" '{print $2}' | sed 1d )

      for sam in $samples
      do
        for tag in $sl_tags
        do
          check_data_tag=$(grep -c $sam $work_directory"/SL_Read_Counts/"$spe"_"$tag"_read.counts" )
          total_for_tag=0
          if [ $check_data_tag -gt 0 ]
          then
            sample_col=$(head -n 1 $work_directory"/SL_Read_Counts/"$spe"_"$tag"_read.counts"  | tr "\t" "\n" | grep -n $sam | awk -F ":" '{print $1}')
            for site in $site_info
            do
              ace=$(echo $site | awk -F "__" '{print $1}')
              chr=$(echo $site | awk -F "__" '{print $2}')
              strand=$(echo $site | awk -F "__" '{print $3}')
              read_counts=0

              gen_id=$(grep -w Concordante $work_directory"/SL_Read_Counts/"$spe"_tier3.counts" | grep -w -v Trans_End | awk -F "\t" -v chr=$chr -v ace=$ace -v strand=$strand '{if ($3==chr && $4==strand && $2==ace) print $1}' )

              check_read_counts=$(grep -w Concordante $work_directory"/SL_Read_Counts/"$spe"_"$tag"_read.counts" | grep -w -v Trans_End | awk -F "\t" -v gen=$gen_id -v chr=$chr -v ace=$ace -v strand=$strand '{if ($1==gen && $3==chr && $4==strand && $5==ace) print}' | grep -c .)
              if [ $check_read_counts -gt 0 ]
              then
                read_counts=$(grep -w Concordante $work_directory"/SL_Read_Counts/"$spe"_"$tag"_read.counts" | grep -w -v Trans_End | awk -F "\t" -v gen=$gen_id -v chr=$chr -v ace=$ace -v strand=$strand -v N=$sample_col '{if ($1==gen && $3==chr && $4==strand && $5==ace) print $N}' | sort -u)
              fi

              echo "Debug: "$spe" "$tag" "$ace": "$read_counts
              total_for_tag=$(($total_for_tag + $read_counts))
            done
          fi
          echo $sam" "$tag" "$total_for_tag | tr " " "\t" >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Frecuencias_per_Sample/"$spe"_reads_per_sample.txt"
        done
      done
    fi
  done
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$INTEREST_OPERONS" == "TRUE" ]
then
  rm -r $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_Arragements"
  mkdir $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_Arragements"

  results_orthofinder=$(ls -d $work_directory"/"$Result_Storage"/OrthoFinder/Results_"*)
  PHO_File=$(echo $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")

  interest_pairs=$(awk -F "\t" '{print $3}' $operones_interes | sort -u)
  species_order=$(echo "Mesocestoides_corti Hmicrostoma Hdiminuta Taenia_multiceps Taenia_saginata Taenia_asiatica Taenia_solium Egranulosus Emultilocularis Spirometra_erinaceieuropaei Schistocephalus_solidus Sparganum_proliferum Trichobilharzia_regenti Schistosoma_japonicum Schistosoma_mansoni Schistosoma_haematobium Schistosoma_bovis Paragonimus_heterotremus Paragonimus_westermani Clonorchis_sinensis Opisthorchis_felineus Fasciolopsis_buski Fasciola_gigantica Fasciola_hepatica")

  for pair in $interest_pairs
  do
    interest_ID=$(grep $pair $operones_interes | awk -F "\t" '{print $1"__"$2}')

    echo "Running: "$pair" --- "$interest_ID
    echo "Species 1st_Gene Union 2nd_Gene Warn1 Warn2 Order" | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_Arragements/"$interest_ID"_"$pair"_Registry.txt"

    first_pho=$(echo $pair | awk -F "__" '{print $1}')
    second_pho=$(echo $pair | awk -F "__" '{print $2}')
    for spe in $species_order
    do
      spe_column=$(head -n 1 $PHO_File | tr "\t" "\n" | grep -n $spe | awk -F ":" '{print $1}')
      abemus_data=FALSE
      check_fist_info=$(grep -c $first_pho $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" )
      if [ $check_fist_info -gt 0 ]
      then
        first_gene_IDs=$(grep $first_pho $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $1}')
      else
        first_gene_IDs=Missing
      fi

      check_second_info=$(grep -c $second_pho $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt")
      if [ $check_second_info -gt 0 ]
      then
        second_gene_IDs=$(grep $second_pho $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $1}')
      else
        second_gene_IDs=Missing
      fi

      check_missing_genes_first=$(grep $first_pho $PHO_File | awk -F "\t" -v N=$spe_column '{print $N}' | tr -d " " | tr "," "\n" |grep -c .)
      check_missing_genes_second=$(grep $second_pho $PHO_File | awk -F "\t" -v N=$spe_column '{print $N}' | tr -d " " | tr "," "\n" |grep -c .)

      if [ $check_fist_info -eq $check_missing_genes_first ]
      then
        warning_first=$(echo .)
      else
        warning_first=X
      fi

      if [ $check_second_info -eq $check_missing_genes_second ]
      then
        warning_second=$(echo .)
      else
        warning_second=X
      fi

      abemus_data=TRUE

      for first_gene in $first_gene_IDs
      do
        if [ $first_gene != "Missing" ]
        then
          first_gene_chr=$(grep -F -w $first_gene $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $2}')
          first_gene_start=$(grep -F -w $first_gene $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $3}')
          first_gene_end=$(grep -F -w $first_gene $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $4}')
          first_gene_strand=$(grep -F -w $first_gene $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $5}')
          first_gene_SL=$(grep -F -w $first_gene $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $6}')
          first_gene_chimera=$(grep -F -w $first_gene $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $7}')
          first_gene_operon=$(grep -F -w $first_gene $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $10}')
        else
          first_gene_chr=Skip
          first_gene_start=Skip
          first_gene_end=Skip
          first_gene_strand=Skip
          first_gene_SL=Skip
          first_gene_chimera=Skip
          first_gene_operon=Skip
        fi

        for second_gene in $second_gene_IDs
        do
          control_order=NA
          if [ $second_gene != "Missing" ]
          then
            second_gene_chr=$(grep -F -w $second_gene $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $2}')
            second_gene_start=$(grep -F -w $second_gene $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $3}')
            second_gene_end=$(grep -F -w $second_gene $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $4}')
            second_gene_strand=$(grep -F -w $second_gene $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $5}')
            second_gene_SL=$(grep -F -w $second_gene $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $6}')
            second_gene_chimera=$(grep -F -w $second_gene $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $7}')
            second_gene_operon=$(grep -F -w $second_gene $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $10}')
          else
            second_gene_chr=Skip
            second_gene_start=Skip
            second_gene_end=Skip
            second_gene_strand=Skip
            second_gene_SL=Skip
            second_gene_chimera=Skip
            second_gene_operon=Skip
          fi

          if [ "$first_gene" == "Missing" ] || [ "$second_gene" == "Missing" ]
          then
            union=Missing
          elif [ "$first_gene" == "$second_gene" ]
          then
            if [ $first_gene_chimera == "Chimeric" ]
            then
              union=Chimera
              control_order=Chimera
            else
              union=Skip
            fi
          elif [ "$first_gene_chr" != "$second_gene_chr" ]
          then
            union=Dif_Chr
          elif [ "$first_gene_strand" != "$second_gene_strand" ]
          then
            distance=$(echo $first_gene" "$first_gene_start" "$first_gene_end"_ooo_"$second_gene" "$second_gene_start" "$second_gene_end | sed "s/_ooo_/\n/" | tr " " "\t" | sort -n -k 2 | awk -F "\t" '{print $2" "$3}' | tr " " "\n" | sed 1d | sed 3d | tr "\n" " " | awk -F " " '{print $2-$1}')
            union=$(echo "Dif_Strand_("$distance")")
          else
            distance=$(echo $first_gene" "$first_gene_start" "$first_gene_end"_ooo_"$second_gene" "$second_gene_start" "$second_gene_end | sed "s/_ooo_/\n/" | tr " " "\t" | sort -n -k 2 | awk -F "\t" '{print $2" "$3}' | tr " " "\n" | sed 1d | sed 3d | tr "\n" " " | awk -F " " '{print $2-$1}')
            if [ "$first_gene_operon" == "$second_gene_operon" ] && [ "$first_gene_operon" != "NA" ]
            then
              union=$(echo "Operon_("$distance")")
            else
              union=$(echo "Colinear_("$distance")")
            fi
            if [ "$first_gene_start" == "Skip" ] || [ "$second_gene_start" == "Skip" ]
            then
              control_order=Skip_1
            else
              if [ "$first_gene_strand" == "$second_gene_strand" ]
              then
                if [ "$first_gene_strand" == "+" ]
                then
                  if [ $first_gene_end  -lt $second_gene_start ]
                  then
                    control_order=$(echo .)
                  else
                    control_order=$(echo X)
                  fi
                else
                  if [ $first_gene_start  -gt $second_gene_end ]
                  then
                    control_order=$(echo .)
                  else
                    control_order=$(echo X)
                  fi
                fi
              else
                control_order=Skip_2
              fi
            fi
          fi

          if [ $first_gene_SL == "SL" ]
          then
            print_first_gene=$(echo $first_gene"(SL)")
          else
            print_first_gene=$first_gene
          fi

          if [ $second_gene_SL == "SL" ]
          then
            print_second_gene=$(echo $second_gene"(SL)")
          else
            print_second_gene=$second_gene
          fi

          echo $spe" "$print_first_gene" "$union" "$print_second_gene" "$warning_first" "$warning_second" "$control_order | tr " " "\t" >>  $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_Arragements/"$interest_ID"_"$pair"_Registry.txt"
        done
      done
    done
  done
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ $INTEREST_OPERONS_CONTROL == TRUE ]
then
  rm -r $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_Arragements/Control_multi_chr"
  mkdir $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_Arragements/Control_multi_chr"
  mkdir $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_Arragements/Control_multi_chr/Hog_CHR"
  mkdir $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_Arragements/Control_multi_chr/Temp"

  results_orthofinder=$(ls -d $work_directory"/"$Result_Storage"/OrthoFinder/Results_"*)
  PHO_File=$(echo $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")

  interest_structures=$(awk -F "\t" '{print $1}' $operones_interes | sort -u)
  species_order=$(echo "Mesocestoides_corti Hmicrostoma Hdiminuta Taenia_multiceps Taenia_saginata Taenia_asiatica Taenia_solium Egranulosus Emultilocularis Spirometra_erinaceieuropaei Schistocephalus_solidus Sparganum_proliferum Trichobilharzia_regenti Schistosoma_japonicum Schistosoma_mansoni Schistosoma_haematobium Schistosoma_bovis Paragonimus_heterotremus Paragonimus_westermani Clonorchis_sinensis Opisthorchis_felineus Fasciolopsis_buski Fasciola_gigantica Fasciola_hepatica")

  for int in $interest_structures
  do
    hogs_int=$(grep -w $int $operones_interes | awk -F "\t" '{print $3}' | sed "s/__/\n/g" | sort -u)
    echo "Especie HOG Self-check "$hogs_int | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_Arragements/Control_multi_chr/"$int"_chr_check.txt"
    for spe in $species_order
    do
      echo $spe" "$int
      spe_column=$(head -n 1 $PHO_File | tr "\t" "\n" | grep -n $spe | awk -F ":" '{print $1}')

      for hog in $hogs_int
      do
        check_data=$(grep -w $hog $PHO_File | awk -F "\t" -v N=$spe_column '{print $N}' | grep -c .)
        if [ $check_data -gt 0 ]
        then
          grep -w $hog $PHO_File | awk -F "\t" -v N=$spe_column '{print $N}' | tr -d " " | tr "," "\n" | sort -u > $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_Arragements/Control_multi_chr/Temp/current_genes.temp"
          grep -w -f $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_Arragements/Control_multi_chr/Temp/current_genes.temp" $work_directory"/"$Result_Storage"/Modified_GTF_Files/"$spe"_with_Chimera_genes.gtf" | awk -F '\t' '{print $1}' >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_Arragements/Control_multi_chr/Hog_CHR/"$spe"_"$hog"_chr.temp"
        fi
      done

      for hog_first in $hogs_int
      do
        if [ -f $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_Arragements/Control_multi_chr/Hog_CHR/"$spe"_"$hog_first"_chr.temp" ]
        then
          check_self=$(sort $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_Arragements/Control_multi_chr/Hog_CHR/"$spe"_"$hog_first"_chr.temp" | uniq -d | grep -c . | awk '{if ($1==0) print "."; else print "X"}')
        else
          check_self=NA
        fi
        registro=$(echo $spe" "$hog_first" "$check_self)

        for hog_second in $hogs_int
        do
          if [ "$hog_first" != "$hog_second" ]
          then
            if [ -f $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_Arragements/Control_multi_chr/Hog_CHR/"$spe"_"$hog_first"_chr.temp" ] && [ -f $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_Arragements/Control_multi_chr/Hog_CHR/"$spe"_"$hog_second"_chr.temp" ]
            then
              shared_chromosomes=$(grep -w -F -c -f $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_Arragements/Control_multi_chr/Hog_CHR/"$spe"_"$hog_first"_chr.temp" $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_Arragements/Control_multi_chr/Hog_CHR/"$spe"_"$hog_second"_chr.temp" | awk '{if ($1==0) print "."; else print "X"}')
            else
              shared_chromosomes=NA
            fi
          else
            shared_chromosomes=Skip
          fi
          registro=$(echo $registro" "$shared_chromosomes)
        done

        echo $registro | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_Arragements/Control_multi_chr/"$int"_chr_check.txt"
      done
      rm $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_Arragements/Control_multi_chr/Hog_CHR/"*"_chr.temp"
    done
  done
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$RUN_INTEREST_ANNOTATION" == "TRUE" ]
then
  rm -r $work_directory"/"$Result_Storage"/Interprot_Results/Interest_PHO"
  mkdir $work_directory"/"$Result_Storage"/Interprot_Results/Interest_PHO"
  mkdir $work_directory"/"$Result_Storage"/Interprot_Results/Interest_PHO/Temp"

  results_orthofinder=$(ls -d $work_directory"/"$Result_Storage"/OrthoFinder/Results_"*)
  PHO_File=$(echo $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")
  Especies=$(ls $work_directory"/"$Result_Storage"/"*pep | sed "s/.*\///" | sed "s/.pep//")

  # 1) PHOs elegidos por consiste4ncia de SL Trans-Splicing
  group_id=Selected
  selected_PHOs=$(sed 1d $work_directory"/"$Result_Storage"/PHO_Selection_TPM/General_Selection_PHO.txt" | awk -F "\t" '{print $1}' | sort -u)
  extract_info_interpro

  group_id=Operon_Good_Pairs
  selected_PHOs=$(awk -F "\t" '{print $2}' $work_directory"/"$Result_Storage"/Operons/Questions/Work_gene_pairs.txt" | sed "s/__/\n/" | sort -u)
  extract_info_interpro

  group_id=Anexo_Construccion
  selected_PHOs=$(cat $interesting_extra_annotation)
  extract_info_interpro
fi
