#!/usr/bin/env bash

# Directory locations and file IDs
work_directory=Analysis_Results/SL_Acceptor_Sites
Result_Storage=Defenitive_Orthogroups
interprot_path=/home/javier/PROYECTOS_2/Interpro/interproscan-5.62-94.0
change_name_itol=/home/amanda/PROYECTS/SLFinder_Publication/Selected_Sequences_Work/change_names_itol.tmp
info_busco_gene_markers=Selected_Sequences_Work/info_mappings_all_busco_datasets_odb10.txt



# Chimera search
NCBI_GFF_format=$(echo "Fasciolopsis_buski Sparganum_proliferum Paragonimus_heterotremus")
Multiple_isoforms=$(echo "Egranulosus Emultilocularis Hdiminuta Hmicrostoma Mesocestoides_corti Opisthorchis_felineus Schistosoma_mansoni Schistosoma_japonicum")
qcov=10

# Mapping Variables
STAR_genomeChrBinNbits=12
STAR_genomeSAindexNbases=12

# Find Operons
operon_range=300
min_read_coverage=5
warning_too_much_reads=100

# Special variables
threads=30
max_sl_rna_loci_distance=5000

# Run sections
PREPARE_TROUBLE_MAKERS=FALSE
PREPARE_ORTHOFINDER=FALSE
RUN_ORTHOFINDER=FALSE
RUN_MAPEOS=FALSE # run "ulimit -n 2048" on the terminal before
PREPARE_CHIMERIC_GTF=FALSE
RUN_CONTEOS=FALSE
RUN_TPM=FALSE
RUN_BASELINE_TPM=FALSE
RUN_MAPEOS_SECOND_PASS=FALSE # run "ulimit -n 2048" on the terminal before
RUN_INTERPRO=FALSE
RUN_OPERON_SEARCH=FALSE
RUN_OPERON_CONSERVATION_QUESTIONS_1=FALSE
RUN_OPERON_CONSERVATION_QUESTIONS_2=FALSE
RUN_OPERON_CONSERVATION_QUESTIONS_3=FALSE
RUN_OPERON_CONSERVATION_QUESTIONS_4=FALSE
RUN_OPERON_CONSERVATION_QUESTIONS_5=FALSE
RUN_OPERON_CONSERVATION_QUESTIONS_6=FALSE
RUN_OPERON_CONSERVATION_QUESTIONS_7=FALSE
RUN_OPERON_CONSERVATION_QUESTIONS_8=FALSE
RUN_HOG_INFORMATION=FALSE
RUN_HOG_SELECTION=FALSE
RUN_HOG_SEQUENCES=FALSE
RUN_HOG_PHYLOGENY=FALSE
RUN_HOG_ITOL_DATASETS=FALSE
RUN_ACCEPTOR_TAG_COUNTS=FALSE
RUN_ACCEPTOR_COMPETENCE_INITIAL=FALSE
RUN_ACCEPTOR_COMPETENCE_SECOND=FALSE
RUN_ACCEPTOR_COMPETENCE_OPERON=FALSE
RUN_ACCEPTOR_COMPETENCE_FINAL=FALSE
RUN_ACCEPTOR_HOG_HEATMAP=FALSE
RUN_BUSCO=FALSE
RUN_BUSCO_ANOTATION=FALSE
RUN_INTEREST_ANNOTATION=FALSE
RUN_TABLAS_PAPER=FALSE
DEFINE_CLUSTER=FALSE
CONTEXT_CLUSTER=FALSE
UP_SET_EXPLORATION=FALSE
UP_SET_INPUT=FALSE
UP_SET_GROUP_QUESTIONS=FALSE
SITE_MAIN_SUMMARY=FALSE
SUP_TABLE_DATA_USED=FALSE
CHIMERA_SUMARY=FALSE
NEW_N0_TABLE=FALSE
ITOL_SPECIES_TREE_DATA=FALSE
SL_READS_FIGURE_TAB_MAKE=FALSE
RAREFACTION_CURVES=FALSE

####################
# Functions:

set_references () {
  if [ $filo == "Cestodes" ]
  then
    comparison_species=$(echo "Emultilocularis Hmicrostoma Sparganum_proliferum Taenia_multiceps")
  else
    comparison_species=$(echo "Schistosoma_mansoni Fasciola_hepatica Clonorchis_sinensis Opisthorchis_felineus")
  fi
}

get_chimera () {
  rename_sequence=$(echo ">"$seq"_Chimera_"$chimeric_interval"__"$id_start_interval"_"$id_end_interval)

  seqkit subseq -r $start_interval":"$end_interval $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/Genefiles/"$seq".fasta" | sed "s/>.*/$rename_sequence/" > $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/Genefiles/current_interval.fasta"
  cat $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/Genefiles/current_interval.fasta" >> $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/"$spe"_chimeric_intervals.fasta"
  getorf -minsize 10 -find 0 -reverse no -sequence $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/Genefiles/current_interval.fasta" -outseq $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/Genefiles/current_interval.peptide"
  best_peptide=$(seqkit fx2tab -n -i -l  $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/Genefiles/current_interval.peptide" | sort -n -k 2 | tail -n 1 | awk -F "\t" '{print $1}')

  seqkit grep -j $threads -p $best_peptide $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/Genefiles/current_interval.peptide" | seqkit seq -i  | sed "s/>.*/$rename_sequence/" >> $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/"$spe"_chimeric_intervals.aa"

  rm $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/Genefiles/current_interval.peptide" $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/Genefiles/current_interval.fasta"
  chimeric_interval=$(($chimeric_interval + 1))
}

competicion_summary () {
  echo "Especie N_Samples Ace Multiple_Genes Total_Competition Total_Pure Percentage_Pure_Total Operon_Competition Operon_Pure Percentage_Operon N_1_Donnor N_2_Donnor N_3_Donnor N_More_Donnor" | tr " " "\t" > $storage"/Summary_ace_per_species.tab"
  for spe in $Especies
  do
    if [ "$is_this_first" == "TRUE" ]
    then
      grupo=$(ls $sj_info"/"*"_"$spe"_"*"_SJ.out.tab" | sed "s/.*\///" | awk -F "_" '{print $1}' | sort -u)
      sj_file=$(echo $sj_info"/"$grupo"_"$spe)
    else
      sj_file=$(echo $sj_info"/"$spe)
    fi

    samples=$(ls $sj_file"_"*"_SJ.out.tab" | sed "s/.*\///" | sed "s/.*$spe//" | awk -F "_" '{print $2}')

    total_samples=$(echo $samples | tr " " "\n" | grep -c .)

    all_sl_data=$(grep -w Concordante $work_directory"/SL_Read_Counts/"$spe"_tier3.counts" | grep -w -v Trans_End | awk -F "\t" '{print $2";"$3";"$4}' | sort -u)
    echo "Aceptor Chromosome Strand Multiple_Genes GenIDs TranscriptIDs Is_Operon Site_Status Total_SJ_Reads Mean_SJ_Reads("$total_samples") Donnor_Sites "$samples | tr " " "\t" >> $storage"/"$spe"_Cis_SJ_on_SL_ace.tab"

    for sl_data in $all_sl_data
    do
      echo "Competicion "$spe": "$sl_data
      basic_info_ace_counts

      total_uniq_reads=0
      registro=sample

      donnor_sites=$(awk -F "\t" -v chr=$chr -v strand=$check_strand '{if ($1==chr && $4==strand && $7>0) print}' $sj_file"_"*"_SJ.out.tab" | awk -F "\t" -v ace=$ace '{if ($2==ace || $3==ace) print $1,$2,$3}' | sort -u | grep -c .)

      for sam in $samples
      do
        sample_uniq_reads=$(awk -F "\t" -v chr=$chr -v strand=$check_strand '{if ($1==chr && $4==strand && $7>0) print}' $sj_file"_"$sam"_SJ.out.tab" | awk -F "\t" -v ace=$ace '{if ($2==ace || $3==ace) print $7}' | awk -F "\t" '{ sum += $1 } END {print sum}')
        check_sj_data=$(echo $sample_uniq_reads | grep -c .)

        if [ $check_sj_data -eq 0 ]
        then
          sample_uniq_reads=0
        fi
        registro=$(echo $registro" "$sample_uniq_reads)
        total_uniq_reads=$(($sample_uniq_reads + $total_uniq_reads))
      done
      registro=$(echo $registro | tr " " "\n" | grep -v sample)

      if [ $total_uniq_reads -gt 0 ]
      then
        site_status=Competition
        mean_uniq_reads=$(echo $total_uniq_reads" "$total_samples | awk -F " " '{printf "%.2f\n", $1/$2}')
      else
        site_status=Pure
        mean_uniq_reads=NA
      fi

      print_gene_id=$(echo $gen_id | tr " " ";" | sed "s/;$//")

      echo $ace" "$chr" "$strand" "$check_double_gene_nonsense" "$print_gene_id" "$trans_ids" "$is_operon" "$site_status" "$total_uniq_reads" "$mean_uniq_reads" "$donnor_sites" "$registro | tr " " "\t" >> $storage"/"$spe"_Cis_SJ_on_SL_ace.tab"
    done

    total_ace=$(grep -c . $storage"/"$spe"_Cis_SJ_on_SL_ace.tab" | awk '{print $1-1}')
    multiple_genes=$(grep -w -c Multiple $storage"/"$spe"_Cis_SJ_on_SL_ace.tab")
    total_competition=$(grep -w -c Competition $storage"/"$spe"_Cis_SJ_on_SL_ace.tab")
    total_pure=$(grep -w -c Pure $storage"/"$spe"_Cis_SJ_on_SL_ace.tab")
    percentage_total=$(echo $total_pure" "$total_competition | awk -F " " '{printf "%.2f\n", $1/($1+$2)}')

    competition_operon=$(awk -F "\t" '{if ($7=="TRUE") print}' $storage"/"$spe"_Cis_SJ_on_SL_ace.tab" | grep -w -c Competition)
    pure_operon=$(awk -F "\t" '{if ($7=="TRUE") print}' $storage"/"$spe"_Cis_SJ_on_SL_ace.tab" | grep -w -c Pure)
    percentage_operon=$(echo $pure_operon" "$competition_operon | awk -F " " '{printf "%.2f\n", $1/($1+$2)}')

    ace_1_donnor=$(awk -F "\t" '{if ($11==1) print }' $storage"/"$spe"_Cis_SJ_on_SL_ace.tab" | grep -c .)
    ace_2_donnor=$(awk -F "\t" '{if ($11==2) print }' $storage"/"$spe"_Cis_SJ_on_SL_ace.tab" | grep -c .)
    ace_3_donnor=$(awk -F "\t" '{if ($11==3) print }' $storage"/"$spe"_Cis_SJ_on_SL_ace.tab" | grep -c .)
    ace_more_donnor=$(awk -F "\t" '{if ($11>3) print }' $storage"/"$spe"_Cis_SJ_on_SL_ace.tab" | grep -v Donnor_Sites| grep -c .)

    echo $spe" "$total_samples" "$total_ace" "$multiple_genes" "$total_competition" "$total_pure" "$percentage_total" "$competition_operon" "$pure_operon" "$percentage_operon" "$ace_1_donnor" "$ace_2_donnor" "$ace_3_donnor" "$ace_more_donnor | tr " " "\t" >> $storage"/Summary_ace_per_species.tab"
  done
}

upset_exploration () {
  cp $data_path"/"*"_"$data_id $work_directory"/"$Result_Storage"/Upset_graph/Temp"
  for ex in $excluded_species
  do
    rm $work_directory"/"$Result_Storage"/Upset_graph/Temp/"$ex"_"*$data_id
    rm $work_directory"/"$Result_Storage"/Upset_graph/Temp/"$ex"_"*$data_id
    rm $work_directory"/"$Result_Storage"/Upset_graph/Temp/"$ex"_"*$data_id
  done

  all_items=$(cat $work_directory"/"$Result_Storage"/Upset_graph/Temp/"*"_"$data_id | sort -u)
  N_all_tiems=$(echo $all_items | tr " " "\n" | grep -c .)

  echo "Hit Total GroupID" | tr " " "\t" >> $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_exploration.txt"
  for item in $all_items
  do
    hits=$(grep -w $item $work_directory"/"$Result_Storage"/Upset_graph/Temp/"*"_"$data_id | awk -F ":" '{print $1}' | sed "s/.*\///" | sed "s/_$data_id//" )
    total_hits=$(echo $hits | tr " " "\n" | grep -c .)
    group_id=$(echo $hits | tr " " ";")

    echo $analysis": "$item" --- "$total_hits

    echo $item" "$total_hits" "$group_id | tr " " "\t" >> $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_exploration.txt"
  done

  N_grupos=$(sed 1d  $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_exploration.txt"  | awk -F "\t" '{print $3}' | sort -u | grep -c .)
  N_grupos_over5=$(sed 1d  $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_exploration.txt"  | awk -F "\t" '{if ($2>=5) print $3}' | sort -u | grep -c .)

  echo "Grupo analisis: "$analysis" ("$N_all_tiems")" >> $work_directory"/"$Result_Storage"/Upset_graph/Resumen_"$analysis"_exploration.txt"
  echo "Total Grupos: "$N_grupos "("$N_grupos_over5")"  >> $work_directory"/"$Result_Storage"/Upset_graph/Resumen_"$analysis"_exploration.txt"
  echo "" >> $work_directory"/"$Result_Storage"/Upset_graph/Resumen_"$analysis"_exploration.txt"
  echo "Top 50 Items with 5 hits:" >> $work_directory"/"$Result_Storage"/Upset_graph/Resumen_"$analysis"_exploration.txt"
  sed "1d" $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_exploration.txt"  | awk -F "\t" '{if ($2>=5) print $1"\t"$2}' | sort -n -r -k 2 | head -n 50 >> $work_directory"/"$Result_Storage"/Upset_graph/Resumen_"$analysis"_exploration.txt"
  echo "" >> $work_directory"/"$Result_Storage"/Upset_graph/Resumen_"$analysis"_exploration.txt"
  echo "Top 50 groups:" >> $work_directory"/"$Result_Storage"/Upset_graph/Resumen_"$analysis"_exploration.txt"
  sed "1d" $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_exploration.txt" | awk -F "\t" '{if ($2>=5) print $3}' | sort | uniq -c | sort -n -r -k 1 | awk '{print $2" "$1}' | head -n 50 >> $work_directory"/"$Result_Storage"/Upset_graph/Resumen_"$analysis"_exploration.txt"
}

confirm_missing_genes () {
  check_missing=$(grep -w $confirm_pho $PHO_File | awk -F "\t" -v N=$spe_column '{print $N}' | grep -c .)
  if [ $check_missing -eq 0 ]
  then
    found_genes=Missing
  else
    found_genes=$(grep -w $confirm_pho $PHO_File | awk -F "\t" -v N=$spe_column '{print $N}'  | tr -d "," | tr " " "\n" | awk '{print $1"[M]"}')
  fi
}

skip_gene_overlap_shenanigans () {
  echo $spe" "$next_feature_ID" "$skip_because >> $work_directory"/"$Result_Storage"/Operons/Skipped_Genes.txt"
  extract=$(($extract + 1))
  count=$(($count + 1))
  if [ $extract -gt $recorrer ]
  then
    use_this_distance=TRUE
    distance=End
  fi
}

get_chimeric_ace () {
  total_SL_Reads=0
  for ace in $work_aces
  do
    temp_count=$(awk -F "\t" -v gen=$extract_gene_id -v N=$read_count_column -v ace=$ace '{if ($1==gen && $2==ace) print $N}' $work_directory"/SL_Read_Counts/"$spe"_total.counts" )

    echo $pho" "$spe" "$transcript" "$ace" "$temp_count" "$strand | tr " " "\t" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Control_Calidad.txt"
    total_SL_Reads=$(($total_SL_Reads + $temp_count))
  done
}

binary_data_itol () {
  echo "DATASET_BINARY" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
  echo "SEPARATOR TAB" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
  echo "" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
  echo "DATASET_LABEL "$dataset_ID | tr " " "\t" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
  echo "COLOR	#000000" | tr " " "\t" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
  echo "MARGIN "$margin | tr " " "\t" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
  echo "HEIGHT_FACTOR "$height_factor | tr " " "\t" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
  echo "FIELD_LABELS "$field_labels | tr " " "\t" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
  echo "FIELD_SHAPES "$field_itol_shape | tr " " "\t" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
  echo "FIELD_COLORS "$field_itol_color | tr " " "\t" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
  echo "DATA" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
}

operon_heatmap () {
  echo $header_operon | tr " " "\t" > $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Heatmap_percentage_"$individual_hog_file
  echo $header_pair | tr " " "\t" > $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Heatmap_percentage_"$pair_hog_file

  echo $header_operon | tr " " "\t" > $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Heatmap_raw_"$individual_hog_file
  echo $header_pair | tr " " "\t" > $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Heatmap_raw_"$pair_hog_file

  for spe1 in $species_order
  do
    N_base_pho_Operon=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe1"_"$individual_hog_file)
    N_base_pho_Operon_pair=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe1"_"$pair_hog_file)

    registro_pho_operon_raw=$(echo $spe1"("$N_base_pho_Operon")")
    registro_pho_operon_pair_raw=$(echo $spe1"("$N_base_pho_Operon_pair")")

    registro_pho_operon_percentage=$(echo $spe1"("$N_base_pho_Operon")")
    registro_pho_operon_pair_percentage=$(echo $spe1"("$N_base_pho_Operon_pair")")

    for spe2 in $species_order
    do
      echo $spe1" vs "$spe2

      if [ "$spe1" == "$spe2" ]
      then
        registro_pho_operon_raw=$(echo $registro_pho_operon_raw" NA")
        registro_pho_operon_pair_raw=$(echo $registro_pho_operon_pair_raw" NA")

        registro_pho_operon_percentage=$(echo $registro_pho_operon_percentage" NA")
        registro_pho_operon_pair_percentage=$(echo $registro_pho_operon_pair_percentage" NA")
      else
        N_pho_operon=$(grep -w -c -f $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe1"_"$individual_hog_file $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe2"_"$individual_hog_file)
        N_pair_operon=$(grep -w -c -f $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe1"_"$pair_hog_file $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe2"_"$pair_hog_file)

        all_pho_operon=$(cat $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe1"_"$individual_hog_file $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe2"_"$individual_hog_file | sort -u | grep -c .)
        all_pair_operon=$(cat $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe1"_"$pair_hog_file $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe2"_"$pair_hog_file | sort -u | grep -c .)

        echo "Total: "$N_pho_operon >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Pairwise/"$spe1"_vs"$spe2"_"$individual_hog_file
        grep -w -f $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe1"_"$individual_hog_file $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe2"_"$individual_hog_file  >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Pairwise/"$spe1"_vs"$spe2"_"$individual_hog_file

        echo "Total: "$N_pair_operon >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Pairwise/"$spe1"_vs"$spe2"_"$pair_hog_file
        grep -w -f $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe1"_"$pair_hog_file $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe2"_"$pair_hog_file >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Pairwise/"$spe1"_vs"$spe2"_"$pair_hog_file

        percentage_pho_operon=$(echo $N_pho_operon" "$all_pho_operon | awk -F " " '{printf "%.10f\n", $1/($1+$2)}')
        percentage_pair_operon=$(echo $N_pair_operon" "$all_pair_operon | awk -F " " '{printf "%.10f\n", $1/($1+$2)}')

        registro_pho_operon_raw=$(echo $registro_pho_operon_raw" "$N_pho_operon)
        registro_pho_operon_pair_raw=$(echo $registro_pho_operon_pair_raw" "$N_pair_operon)

        registro_pho_operon_percentage=$(echo $registro_pho_operon_percentage" "$percentage_pho_operon)
        registro_pho_operon_pair_percentage=$(echo $registro_pho_operon_pair_percentage" "$percentage_pair_operon)
      fi
    done

    echo $registro_pho_operon_raw | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Heatmap_raw_"$individual_hog_file
    echo $registro_pho_operon_pair_raw | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Heatmap_raw_"$pair_hog_file

    echo $registro_pho_operon_percentage | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Heatmap_percentage_"$individual_hog_file
    echo $registro_pho_operon_pair_percentage | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Heatmap_percentage_"$pair_hog_file
  done
}

operon_heatmap_asimetrico () {
  echo $header_operon | tr " " "\t" > $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Heatmaps_Asimetricos/Heatmap_asimetrico_percentage_"$individual_hog_file
  echo $header_pair | tr " " "\t" > $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Heatmaps_Asimetricos/Heatmap_asimetrico_percentage_"$pair_hog_file

  echo $header_operon | tr " " "\t" > $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Heatmaps_Asimetricos/Heatmap_asimetrico_raw_"$individual_hog_file
  echo $header_pair | tr " " "\t" > $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Heatmaps_Asimetricos/Heatmap_asimetrico_raw_"$pair_hog_file

  for spe1 in $species_order
  do
    N_base_pho_Operon=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe1"_"$individual_hog_file)
    N_base_pho_Operon_pair=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe1"_"$pair_hog_file)

    registro_pho_operon_raw=$(echo $spe1"("$N_base_pho_Operon")")
    registro_pho_operon_pair_raw=$(echo $spe1"("$N_base_pho_Operon_pair")")

    registro_pho_operon_percentage=$(echo $spe1"("$N_base_pho_Operon")")
    registro_pho_operon_pair_percentage=$(echo $spe1"("$N_base_pho_Operon_pair")")

    for spe2 in $species_order
    do
      echo $spe1" vs "$spe2

      if [ "$spe1" == "$spe2" ]
      then
        registro_pho_operon_raw=$(echo $registro_pho_operon_raw" NA")
        registro_pho_operon_pair_raw=$(echo $registro_pho_operon_pair_raw" NA")

        registro_pho_operon_percentage=$(echo $registro_pho_operon_percentage" NA")
        registro_pho_operon_pair_percentage=$(echo $registro_pho_operon_pair_percentage" NA")
      else
        N_pho_operon=$(grep -w -c -f $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe1"_"$individual_hog_file $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe2"_"$individual_hog_file)
        N_pair_operon=$(grep -w -c -f $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe1"_"$pair_hog_file $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe2"_"$pair_hog_file)

        all_pho_operon=$N_base_pho_Operon
        all_pair_operon=$N_base_pho_Operon_pair

        percentage_pho_operon=$(echo $N_pho_operon" "$all_pho_operon | awk -F " " '{printf "%.10f\n", $1/($1+$2)}')
        percentage_pair_operon=$(echo $N_pair_operon" "$all_pair_operon | awk -F " " '{printf "%.10f\n", $1/($1+$2)}')

        registro_pho_operon_raw=$(echo $registro_pho_operon_raw" "$N_pho_operon)
        registro_pho_operon_pair_raw=$(echo $registro_pho_operon_pair_raw" "$N_pair_operon)

        registro_pho_operon_percentage=$(echo $registro_pho_operon_percentage" "$percentage_pho_operon)
        registro_pho_operon_pair_percentage=$(echo $registro_pho_operon_pair_percentage" "$percentage_pair_operon)
      fi
    done

    echo $registro_pho_operon_raw | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Heatmaps_Asimetricos/Heatmap_asimetrico_raw_"$individual_hog_file
    echo $registro_pho_operon_pair_raw | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Heatmaps_Asimetricos/Heatmap_asimetrico_raw_"$pair_hog_file

    echo $registro_pho_operon_percentage | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Heatmaps_Asimetricos/Heatmap_asimetrico_percentage_"$individual_hog_file
    echo $registro_pho_operon_pair_percentage | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Heatmaps_Asimetricos/Heatmap_asimetrico_percentage_"$pair_hog_file
  done
}

interest_pho_sequences () {
  for pho in $selected_PHOs
  do
    echo "Sequence Average_TPM N_SL_Reads Species_Median_TPM ALL_Species_SL_Reads" | tr " " "\t" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/"$pho"_expresion_data.txt"
    echo "Species Gen Transcript Ace Strand N_SL_Reads" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/"$pho"_SL_ace_information.txt"

    genes_pho=$(grep -F -w $pho $work_directory"/"$Result_Storage"/PHO_Selection_TPM/PHO_Gene_Registry.txt" | awk -F "\t" '{print $3}')
    for transcript in $genes_pho
    do
      spe=$(grep -w -F $transcript $work_directory"/"$Result_Storage"/PHO_Selection_TPM/PHO_Gene_Registry.txt" | awk -F "\t" '{print $2}')
      read_count_column=$(head -n 1 $work_directory"/SL_Read_Counts/"$spe"_total.counts" | tr "\t" "\n" | grep -n -w Total  )
      transcriptID=$(echo $transcript | sed "s/_Chimera.*//")

      species_SL_median=$(grep -w $spe $work_directory"/"$Result_Storage"/Mediciones_Expresion/Background_TPM/SL_Reference_median.txt" | awk -F "\t" '{print $2}')
      Species_SL_Reads=$(awk -F "\t" -v N=$read_count_column '{if ($5=="Concordante" && $6!="Trans_End") print $N}' $work_directory"/SL_Read_Counts/"$spe"_total.counts" | awk -F "\t" '{ sum += $1 } END {print sum}')

      # Retrieve ITOL_Data
      average_TPM=$(grep -w -F $transcript $work_directory"/"$Result_Storage"/PHO_Selection_TPM/PHO_Gene_Registry.txt" | awk -F "\t" '{print $4}')

      # Incluir los genes quiméricos a complejizado todos los conteos de splicing
      # Este bloque revisa que sitios aceptores son incluidos dentro de cada intervalo y los suma para el registro
      # Una consecuencia de las decisiones tomadas es que algunos sitios aceptores secundarios pueden, teoricamente, no ser incluidos.
      #   --->>> Aceptar su perdida hacia el 3', extender selección hacia el 5'

      if [ $spe == "Fasciolopsis_buski" ] || [ $spe == "Paragonimus_heterotremus" ]
      then
        extract_trans_id=$(echo "gene-"$transcriptID)
      else
        extract_trans_id=$transcriptID
      fi

      echo $pho" "$transcript" "$spe
      check_SL=$(grep -w -F $extract_trans_id $work_directory"/SL_Read_Counts/"$spe"_"*"_read.counts" | awk -F "\t" '{if ($6=="Concordante" && $7!="Trans_End") print}' | grep -c .)
      total_SL_Reads=0

      if [ $check_SL -gt 0 ]
      then
        extract_gene_id=$(grep -w -F $transcriptID $work_directory"/SL_Read_Counts/"$spe"_"*"_read.counts" | awk -F "\t" '{print $1}' | awk -F ":" '{print $2}' | sort -u)

        check_chimera=$(echo $transcript | grep -c "Chimera")
        check_chimera_First=$(echo $transcript | grep -c "_Start_")
        check_chimera_End=$(echo $transcript | grep -c "_End$")

        if [ $check_chimera -eq 0 ]
        then
          work_aces=$(grep -w -F $transcript $work_directory"/SL_Read_Counts/"$spe"_"*"_read.counts" | awk -F "\t" '{if ($6=="Concordante" && $7!="Trans_End") print $5}' | sort -u)
          for ace in $work_aces
          do
            awk -F "\t" -v species=$spe -v gen=$extract_gene_id -v trans=$transcript -v ace=$ace -v N=$read_count_column '{if ($1==gen && $2==ace) print species"\t"$1"\t"trans"\t"$2"\t"$4"\t"$N}' $work_directory"/SL_Read_Counts/"$spe"_total.counts" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/"$pho"_SL_ace_information.txt"
            count_aces=$(awk -F "\t" -v gen=$extract_gene_id -v ace=$ace -v N=$read_count_column '{if ($1==gen && $2==ace) print $N}' $work_directory"/SL_Read_Counts/"$spe"_total.counts" )
            total_SL_Reads=$(($total_SL_Reads + $count_aces))
          done
        else
          strand=$(awk -F "\t" -v gen=$extract_gene_id '{if ($1==gen) print $4}' $work_directory"/SL_Read_Counts/"$spe"_total.counts" | sort -u)

          awk -F "\t" -v species=$spe -v gen=$extract_gene_id -v trans=$transcriptID -v ace=$ace -v N=$read_count_column '{if ($1==gen) print species"\t"$1"\t"trans"\t"$2"\t"$4"\t"$N}' $work_directory"/SL_Read_Counts/"$spe"_total.counts" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/"$pho"_SL_ace_information.txt"

          if [ $check_chimera_First -gt 0 ]
          then
            # Contar reads SL en Chimera 1
            border_ace=$(echo $transcript | sed "s/.*_Start_//" )
            if [ "$strand" == "+" ]
            then
              work_aces=$(grep -w -F $transcriptID $work_directory"/SL_Read_Counts/"$spe"_"*"_read.counts" | awk -F "\t" -v border=$border_ace '{if ($5<border && $6=="Concordante" && $7!="Trans_End") print $5}' | sort -u)
              get_chimeric_ace
            else
              work_aces=$(grep -w -F $transcriptID $work_directory"/SL_Read_Counts/"$spe"_"*"_read.counts" | awk -F "\t" -v border=$border_ace '{if ($5>border && $6=="Concordante" && $7!="Trans_End") print $5}' | sort -u)
              get_chimeric_ace
            fi
          elif [ $check_chimera_End -gt 0 ]
          then
            extra_aces=$(echo $transcript | sed "s/.*__//" | sed "s/_End//" )
            Control_1=$extra_aces
            Control_2=End
            expand_border
            if [ "$strand" == "+" ]
            then
              border_ace=$(echo $extra_aces | tr " " "\n" | sort -n -u | head -n 1)
              work_aces=$(grep -w -F $transcriptID $work_directory"/SL_Read_Counts/"$spe"_"*"_read.counts" | awk -F "\t" -v border=$border_ace '{if ($5>=border && $6=="Concordante" && $7!="Trans_End") print $5}' | sort -u)
              get_chimeric_ace
            else
              border_ace=$(echo $extra_aces | tr " " "\n" | sort -n -u | tail -n 1)
              work_aces=$(grep -w -F $transcriptID $work_directory"/SL_Read_Counts/"$spe"_"*"_read.counts" | awk -F "\t" -v border=$border_ace '{if ($5<=border && $6=="Concordante" && $7!="Trans_End") print $5}' | sort -u)
              get_chimeric_ace
            fi
          else
            border_first_ace=$(echo $transcript | sed "s/.*__//" | awk -F "_" '{print $1}')
            border_last_ace=$(echo $transcript | sed "s/.*__//" | awk -F "_" '{print $2}')

            Control_1=$border_first_ace
            Control_2=$border_last_ace

            if [ "$strand" == "+" ]
            then
              extra_aces=$(echo $Control_1" "$Control_2 | tr " " "\n" | sort -n | head -n 1)
              expand_border
              border_first_ace=$(echo $extra_aces | tr " " "\n" | sort -n -u | head -n 1)
              work_aces=$(grep -w -F $transcriptID $work_directory"/SL_Read_Counts/"$spe"_"*"_read.counts" | awk -F "\t" -v first_border=$border_first_ace -v second_border=$border_last_ace '{if ($5>=first_border && $5<second_border && $6=="Concordante" && $7!="Trans_End") print $5}' | sort -u)
              get_chimeric_ace
            else
              extra_aces=$(echo $Control_1" "$Control_2 | tr " " "\n" | sort -n | tail -n 1)
              expand_border
              border_first_ace=$(echo $extra_aces | tr " " "\n" | sort -n -u | tail -n 1)
              work_aces=$(grep -w -F $transcriptID $work_directory"/SL_Read_Counts/"$spe"_"*"_read.counts" | awk -F "\t" -v first_border=$border_first_ace -v second_border=$border_last_ace '{if ($5<=first_border && $5>second_border && $6=="Concordante" && $7!="Trans_End") print $5}' | sort -u)
              get_chimeric_ace
            fi
          fi
        fi
      fi
      echo $spe"_"$transcript" "$average_TPM" "$total_SL_Reads" "$species_SL_median" "$Species_SL_Reads | tr " " "\t" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/"$pho"_expresion_data.txt"

      # Retrieve Sequences
      rename=$(echo ">"$spe"_")
      seqkit grep -p $transcript $work_directory"/"$Result_Storage"/"$spe".pep" | sed "s/>/$rename/" >> $work_directory"/"$Result_Storage"/"$storage_place"/Sequences/"$pho"_sequences.aa"
    done
  done
}

interest_pho_phylogeny () {
  for pho in $selected_PHOs
  do
    echo $pho" MAFFT"
    mafft --quiet --maxiterate 1000 --localpair --thread $threads $work_directory"/"$Result_Storage"/"$storage_place"/Sequences/"$pho"_sequences.aa" | seqkit seq -w 60 -u >> $work_directory"/"$Result_Storage"/"$storage_place"/Phylogeny/Alignment/"$pho".aln"
    echo $pho" iqtree"
    iqtree --seqtype "AA" --quiet -T $threads -s $work_directory"/"$Result_Storage"/"$storage_place"/Phylogeny/Alignment/"$pho".aln" --prefix $work_directory"/"$Result_Storage"/"$storage_place"/Phylogeny/"$pho -m MFP -B 1000
  done
}

interest_pho_itol () {
  echo "Egranulosus Emultilocularis Hdiminuta Hmicrostoma Mesocestoides_corti Sparganum_proliferum Spirometra_erinaceieuropaei Taenia_asiatica Taenia_multiceps Taenia_saginata Taenia_solium Schistocephalus_solidus Clonorchis_sinensis Fasciola_gigantica Fasciola_hepatica Fasciolopsis_buski Opisthorchis_felineus Paragonimus_heterotremus Paragonimus_westermani Schistosoma_bovis Schistosoma_haematobium Schistosoma_japonicum Schistosoma_mansoni Trichobilharzia_regenti" | tr " " "\n" > $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/find_species.txt"

  # 1) Información sobre el TPM Relativo a la mediana de detección
  dataset_ID=SL_Trans-Splicing
  width=100
  margin=0
  upload_order=1
  height_factor=$(echo "0.5")
  field_itol_shape=3
  field_labels=$(echo "f1")
  field_itol_color=$(echo "#000000")

  for pho in $selected_PHOs
  do
    binary_data_itol
    sequences=$(seqkit seq -n -i $work_directory"/"$Result_Storage"/"$storage_place"/Sequences/"$pho"_sequences.aa")
    for seq in $sequences
    do
      echo $seq" "$dataset_ID
      temp_value=$(grep -w -F $seq $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/"$pho"_expresion_data.txt" | awk -F "\t" '{print $3}')
      if [ $temp_value -eq 0 ]
      then
        value=$(echo "-1")
      elif [ $temp_value -lt 4 ]
      then
        value=$(echo "0")
      else
        value=$(echo "1")
      fi
      echo $seq" "$value | tr " " "\t" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
    done
  done

  # 2) Información sobre el TPM Relativo a la mediana de detección
  dataset_ID=Relative_TPM
  width=100
  margin=0
  itol_color=$(echo "#00ff16")
  barshift=$(echo "1.5")
  upload_order=2
  value_line=1
  ref_line=$(echo $value_line"-1x Relative TPM-"$itol_color"-1-1-1")

  for pho in $selected_PHOs
  do
    simplebar_data_itol

    sequences=$(seqkit seq -n -i $work_directory"/"$Result_Storage"/"$storage_place"/Sequences/"$pho"_sequences.aa")
    for seq in $sequences
    do
      echo $seq" "$dataset_ID
      check_saturated_value=$(grep -w -F $seq $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/"$pho"_expresion_data.txt" | awk -F "\t" '{if ($2/$4 > 10) print "TRUE"; else print "FALSE"}')
      if [ "$check_saturated_value" == "TRUE" ]
      then
        value=10
      else
        value=$(grep -w -F $seq $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/"$pho"_expresion_data.txt" | awk -F "\t" '{printf "%.2f\n", $2/$4}')
      fi

      echo $seq" "$value | tr " " "\t" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
    done
  done

  # 3) Información sobre el SL Reads a la mediana de detección
  dataset_ID=Log10_SL_Reads
  itol_color=$(echo "#ff0000")
  margin=$(echo "-"$width)
  barshift=$(echo "-1.5")
  upload_order=3
  value_line=$(echo 0.69897)
  ref_line=$(echo $value_line"-4 SL Reads-"$itol_color"-1-1-1")

  for pho in $selected_PHOs
  do
    simplebar_data_itol

    sequences=$(seqkit seq -n -i $work_directory"/"$Result_Storage"/"$storage_place"/Sequences/"$pho"_sequences.aa")
    for seq in $sequences
    do
      echo $seq" "$dataset_ID
      temp_value=$(grep -w -F $seq $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/"$pho"_expresion_data.txt" | awk -F "\t" '{print $3+1}')

      if [ $temp_value -gt 100 ]
      then
        value=$(echo 100 | awk '{print log($1)/log(10)}')
      else
        value=$(echo $temp_value | awk '{print log($1)/log(10)}')
      fi

      echo $seq" "$value | tr " " "\t" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
    done
  done

  # 4) Control sobre la saturación
  dataset_ID=Saturation
  width=100
  margin=10
  field_itol_color=$(echo "#00ff16 #ff0000")
  field_itol_shape=$(echo "2 2")
  field_labels=$(echo "f1 f2")
  height_factor=$(echo "0.5")
  upload_order=4

  for pho in $selected_PHOs
  do
    binary_data_itol
    sequences=$(seqkit seq -n -i $work_directory"/"$Result_Storage"/"$storage_place"/Sequences/"$pho"_sequences.aa")
    for seq in $sequences
    do
      echo $pho": "$seq" "$dataset_ID
      check_saturation_TPM=$(grep -w -F $seq $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/"$pho"_expresion_data.txt" | awk -F "\t" '{if ($2/$4 > 10) print "TRUE"; else print "FALSE"}')
      if [ $check_saturation_TPM == "TRUE" ]
      then
        TMP_value=$(echo "1")
      else
        TMP_value=$(echo "-1")
      fi

      check_saturation_SL=$(grep -w -F $seq $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/"$pho"_expresion_data.txt" |  awk -F "\t" '{if ($3+1 > 100) print "TRUE"; else print "FALSE"}')
      if [ $check_saturation_SL == "TRUE" ]
      then
        SL_value=$(echo "1")
      else
        SL_value=$(echo "-1")
      fi

      echo $seq" "$TMP_value" "$SL_value | tr " " "\t" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
    done
  done

  # 5) Parte de un Operon
  dataset_ID=Operon_Members
  width=100
  margin=0
  upload_order=5
  height_factor=$(echo "0.5")
  field_itol_shape=6
  field_labels=$(echo "f1")
  field_itol_color=$(echo "#0000ff")

  for pho in $selected_PHOs
  do
    binary_data_itol
    sequences=$(seqkit seq -n -i $work_directory"/"$Result_Storage"/"$storage_place"/Sequences/"$pho"_sequences.aa")
    for seq in $sequences
    do
      current_species=$(echo $seq | grep -o -f $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/find_species.txt")
      grupo=$(ls $work_directory"/"$Result_Storage"/Mediciones_Expresion/Mapeos/"*"_"$current_species"_"*"_Aligned.sortedByCoord.out.bam" | sed "s/.*\///" | sed "s/$current_species.*//" | sed "s/_$//" | sort -u)
      gff_file=$(ls $grupo"/"$current_species"/"*"annotations.gtf")

      search_seq=$(echo $seq | sed "s/_Chimera.*//" | sed "s/$current_species//" | sed "s/^_//")

      check_trouble=$(echo "Sparganum_proliferum Paragonimus_heterotremus Fasciolopsis_buski" | grep -c $current_species)
      check_trouble2=$(echo "Taenia_multiceps Schistosoma_mansoni" | grep -c $current_species)

      if [ $check_trouble -gt 0 ]
      then
        gene_ID=$(grep $search_seq'";' $gff_file | awk -F "\t" '{print $9}' | tr ";" "\n" | grep gene_id | sed 's/gene_id "//' | tr -d '"' | tr -d ' ' | sort -u | sed "s/^gene-//")
      elif [ $check_trouble2 -gt 0 ]
      then
        gene_ID=$(grep $search_seq'";' $gff_file | awk -F "\t" '{print $9}' |  tr ";" "\n" | grep gene_id | sed 's/gene_id "gene://' | sed 's/gene_id "transcript://' | tr -d '"' | tr -d ' ' | sort -u)
      else
        gene_ID=$(grep $search_seq'";' $gff_file | awk -F "\t" '{print $9}' |  tr ";" "\n" | grep gene_id | sed 's/gene_id "gene://' | tr -d '"' | tr -d ' ' | sort -u)
      fi

      echo $pho": "$seq" "$dataset_ID" || "$gene_ID
      check_is_operon=$(grep -w -F $gene_ID $work_directory"/"$Result_Storage"/Operons/"$current_species"_gene.txt" | grep -c Operon)

      if [ $check_is_operon -gt 0 ]
      then
        echo $seq" 1" | tr " " "\t" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
      else
        echo $seq" -1" | tr " " "\t" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
      fi
    done
  done

  # 5) Parte de un Operon
  dataset_ID=Rename_branches
  upload_order=6

  echo ""
  echo "##################"
  echo ""

  for pho in $selected_PHOs
  do
    echo "LABELS" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
    echo "SEPARATOR TAB" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
    echo "" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
    echo "DATA" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"

    sequences=$(seqkit seq -n -i $work_directory"/"$Result_Storage"/"$storage_place"/Sequences/"$pho"_sequences.aa")
    for seq in $sequences
    do
      current_species=$(echo $seq | grep -o -f $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/find_species.txt")
      change_seq=$(echo $current_species"_")
      replace_name=$(grep -w $current_species $change_name_itol | awk -F "\t" '{print $2}')
      new_name=$(echo $seq | sed "s/$change_seq/$replace_name /")

      echo "seq:"$seq
      echo $change_seq" || "$replace_name
      echo "new_name:"$new_name

      echo $seq";"$new_name | tr ";" "\t" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
    done
  done
}

simplebar_data_itol () {
  echo "DATASET_SIMPLEBAR" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
  echo "SEPARATOR TAB" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
  echo "" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
  echo "DATASET_LABEL "$dataset_ID | tr " " "\t" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
  echo "COLOR	"$itol_color | tr " " "\t" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
  echo "WIDTH "$width | tr " " "\t" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
  echo "MARGIN "$margin | tr " " "\t" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
  echo "HEIGHT_FACTOR 0.1" | tr " " "\t" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
  echo "BAR_SHIFT "$barshift | tr " " "\t" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
  echo "DATASET_SCALE;"$ref_line | tr ";" "\t" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
  echo "" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
  echo "DATA" >> $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets/"$pho"_"$upload_order"_"$dataset_ID"_Dataset.txt"
}

expand_border () {
  # En la mayoiría de los casos encontrar el sitio ACE apropiado para contar es trivial. Sin embargo cuando se combina:
    # a) genes quimericos
    # b) Multiples sitios aceptores internos conde 2 o más intervalos fueron combinados en uno solo
  # En estos casos los puntos de control pueden aparecer en más de una linea en el archivo "All_cut_points.txt" y lo que quiero es la menor existente para las comparaciones

  first_find_hit_search=$(grep -w -F $transcriptID $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/All_cut_points.txt" | grep -n . | awk -F "\t" -v control=$Control_1 '{if ($5==control || $6==control) print}' | awk -F ":" '{print $1}')
  second_find_hit_search=$(grep -w -F $transcriptID $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/All_cut_points.txt" | grep -n . | awk -F "\t" -v control=$Control_2 '{if ($5==control || $6==control) print}' | awk -F ":" '{print $1}')
  current_line=$(echo $first_find_hit_search" "$second_find_hit_search | tr " " "\n" | sort -n -u | head -n 1)

  while [ $current_line -gt 1 ]
  do
    check_annotation_presence=$(grep -w -F $transcriptID $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/All_cut_points.txt" | sed -n $current_line"p" | awk -F "\t" '{print $4}' | grep -c .)
    check_annotation_abort=$(grep -w -F $transcriptID $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/All_cut_points.txt" | sed -n $current_line"p" | awk -F "\t" '{print $4}')
    if [ $check_annotation_presence -eq 0 ] || [ "$check_annotation_abort" == "ABORT" ]
    then
      new_aces=$(grep -w -F $transcriptID $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/All_cut_points.txt" | sed -n $current_line"p" | awk -F "\t" '{print $5,$6}' | tr " " "\n" | grep -v Start | grep -v End)
      extra_aces=$(echo $extra_aces" "$new_aces)
      current_line=$(($current_line - 1))
    else
      current_line=0
    fi
  done
}

check_orthogroups () {
  if [ $artifact != "Chimeric" ]
  then
    retrive_transcript_id=$(echo "Schistosoma_mansoni Schistosoma_haematobium Hdiminuta Hmicrostoma Emultilocularis Egranulosus Taenia_saginata Taenia_asiatica Opisthorchis_felineus" | grep -c $spe)

    if [ $retrive_transcript_id -gt 0 ]
    then
      cp $work_directory"/"$Result_Storage"/Operons/Temporal/"$grupo"_"$spe"_transcript.temp" $work_directory"/"$Result_Storage"/Operons/Temporal/transcript.id"
    elif [ $check_trouble3 -gt 0 ]
    then
      echo $gene_id | sed "s/^gene-//" > $work_directory"/"$Result_Storage"/Operons/Temporal/transcript.id"
    else
      echo $gene_id > $work_directory"/"$Result_Storage"/Operons/Temporal/transcript.id"
    fi

    check_orthofinder_data=$(grep -c -F -w -f $work_directory"/"$Result_Storage"/Operons/Temporal/transcript.id" $PHO_File)
    if [ $check_orthofinder_data -eq 1 ]
    then
      pho=$(grep -F -w -f $work_directory"/"$Result_Storage"/Operons/Temporal/transcript.id" $PHO_File | awk -F "\t" '{print $1}')
    elif [ $check_orthofinder_data -eq 0 ]
    then
      pho=No_Hits
    else
      echo "Please Reconsider"
      echo $spe" "$gene_id
      exit
    fi
  else
    pho=XX_
    for chim in $chimeric_ids
    do
      check_pho_hit=$(grep -c -F $chim"__" $PHO_File)
      if [ $check_pho_hit -gt 0 ]
      then
        temp_pho=$(grep -F $chim"__" $PHO_File | awk -F "\t" '{print $1}')
      else
        temp_pho=No_Hits
      fi
      pho=$(echo $pho"__"$temp_pho)
    done
    pho=$(echo $pho | sed "s/^XX___//")
  fi
}

check_read_depth () {
  depth_star=$(($end_feature + 1))
  depth_end=$(($next_feature_start - 1))

  should_have_data=$(($depth_end - $depth_star + 1))

  for read in $read_files
  do
    samtools depth -r $chr":"$depth_star"-"$depth_end $work_directory"/"$Result_Storage"/Mediciones_Expresion/Mapeos/"$grupo"_"$spe"_"$read"_Aligned.sortedByCoord.out.bam" > $work_directory"/"$Result_Storage"/Operons/Temporal/"$read"_coverage.tmp"
  done

  if [ $N_read_files -gt 1 ]
  then
    check_data=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Temporal/"*"_coverage.tmp" | awk -F ":" '{print $2}' | sort -n -u | tail -n 1)
  else
    check_data=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Temporal/"$read"_coverage.tmp")
  fi

  if [ $check_data -gt 0 ]
  then
    if [ $check_data -lt $should_have_data ]
    then
      read_coverage=interrupted_missing
    else
      internal_positions=$(awk -F "\t" '{print $2}' $work_directory"/"$Result_Storage"/Operons/Temporal/"*"_coverage.tmp" | sort -u)

      for int in $internal_positions
      do
        suma=$(grep -w $int $work_directory"/"$Result_Storage"/Operons/Temporal/"*"_coverage.tmp" | awk -F '\t' '{ sum += $3 } END { if (NR > 0) print sum/NR}')
        echo $int" "$suma | tr " " "\t" >>  $work_directory"/"$Result_Storage"/Operons/Temporal/All_work_coverage.tmp"
      done

      low_cov_bases=$(awk -F "\t" -v too_low=$min_read_coverage '{if ($2<=too_low) print}' $work_directory"/"$Result_Storage"/Operons/Temporal/All_work_coverage.tmp" | grep -c . )
      if [ $low_cov_bases -gt 0 ]
      then
        read_coverage=interrupted_low_cov
      else
        read_coverage=$(awk '{ sum += $2 } END { if (NR > 0) print sum / NR }' $work_directory"/"$Result_Storage"/Operons/Temporal/All_work_coverage.tmp" | awk -F "\t" -v high=$warning_too_much_reads '{if ($1>=high) print "high"; else print "good" }')
      fi
    fi
  else
    read_coverage=critical_missing
  fi
  rm $work_directory"/"$Result_Storage"/Operons/Temporal/"*"_coverage.tmp"
}


calculate_counts () {
  for spe in $Especies
  do
    read_files=$(ls $grupo"/"$spe"/"*"_work_read_1P.gz" | sed "s/.*\///" | sed "s/_work_read_1P.gz//")
    gff_file=$work_directory"/"$Result_Storage"/Modified_GTF_Files/"$spe"_with_Chimera_genes.gtf"

    for read in $read_files
    do
      echo $spe" "$read
      htseq-count -r pos -s no --idattr=ID --format=bam --type=transcript -q $work_directory"/"$Result_Storage"/Mediciones_Expresion/Mapeos/"$grupo"_"$spe"_"$read"_Aligned.sortedByCoord.out.bam" $gff_file > $storage_count"/"$grupo"_"$spe"_"$read"_htseq.counts"
    done

    count_da_genes=$(awk -F "\t" '{print $1}' $storage_count"/"$grupo"_"$spe"_"$read"_htseq.counts" | grep -v "^__" )

    echo "GenID "$read_files| tr " " "\t" >> $storage_count"/"$grupo"_"$spe"_All_htseq.counts"
    for gen in $count_da_genes
    do
      print_gene=$(echo $gen | sed "s/^gene://")

      echo "Counting "$grupo" "$spe": "$gen
      record=$print_gene

      for read in $read_files
      do
        counting=$(grep -w -F $gen $storage_count"/"$grupo"_"$spe"_"$read"_htseq.counts" | awk -F "\t" '{print $2}')
        record=$(echo $record" "$counting)
      done
      echo $record | tr " " "\t" >> $storage_count"/"$grupo"_"$spe"_All_htseq.counts"
    done
  done
}

calculate_TPM () {
  for spe in $Especies
  do
    # Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
    # Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
    # Divide the RPK values by the “per million” scaling factor. This gives you TPM.

    echo "RPK Making "$grupo" "$spe

    read_files=$(ls $grupo"/"$spe"/"*"_work_read_1P.gz" | sed "s/.*\///" | sed "s/_work_read_1P.gz//")
    extract_cols=$(ls $grupo"/"$spe"/"*"_work_read_1P.gz" | grep -c . | awk '{print $1+1}')

    seqkit fx2tab -n -i -l $work_directory"/"$Result_Storage"/"$spe".pep" | awk -F "\t" '{print $1"\t"$2*3}' >> $storage_tpm"/Temporales/"$spe"_temp.len"
    head -n 1 $storage_count"/"$grupo"_"$spe"_All_htseq.counts" >> $storage_tpm"/RPK/"$grupo"_"$spe"_RPK.counts"

    count_da_genes=$(awk -F "\t" '{print $1}' $storage_count"/"$grupo"_"$spe"_All_htseq.counts" | grep -v "GenID" )

    for gen in $count_da_genes
    do
      gen_len=$(grep -w -F $gen $storage_tpm"/Temporales/"$spe"_temp.len" | awk -F "\t" '{printf "%.10f\n", $2/1000}')

      echo "RPK :"$grupo"_"$spe": "$gen"("$gen_len")"
      grep -w -F $gen $storage_count"/"$grupo"_"$spe"_All_htseq.counts" | sed "s/\t/\n/" | tail -n 1
      RPK=$(grep -w -F $gen $storage_count"/"$grupo"_"$spe"_All_htseq.counts" | sed "s/\t/\n/" | tail -n 1 | tr "\t" "\n" | awk -F "\t" -v len=$gen_len '{printf "%.10f\n", $1/len}' | tr "\n" " ")
      echo "Done RPK"

      echo $gen" "$RPK | tr " " "\t" >> $storage_tpm"/RPK/"$grupo"_"$spe"_RPK.counts"
    done

    echo "TPM Making "$grupo" "$spe
    extr=2

    awk -F "\t" '{print $1}' $storage_tpm"/RPK/"$grupo"_"$spe"_RPK.counts" > $storage_tpm"/Temporales/making_tpm"

    while [ $extr -le $extract_cols ]
    do
      read_id=$(head -n 1 $storage_tpm"/RPK/"$grupo"_"$spe"_RPK.counts" | awk -F "\t" -v col=$extr '{print $col}')
      factor=$(grep -v "GenID" $storage_tpm"/RPK/"$grupo"_"$spe"_RPK.counts" | awk -F "\t" -v col=$extr '{print $col}' | awk -F "\t" '{ sum += $1 } END {print sum/1000000}')

      echo $read_id" "$factor | tr " " "\t" >> $storage_tpm"/"$grupo"_"$spe"_factors.txt"

      echo "TPM "$spe": "$read_id" ("$factor")"
      grep -v "GenID" $storage_tpm"/RPK/"$grupo"_"$spe"_RPK.counts" | awk -F "\t" -v col=$extr '{print $col}' | awk -v factor=$factor '{printf "%.3f\n", $1/factor}' | tr "\n" " " | awk -v read=$read_id '{print read" "$0}' | tr " " "\n" | grep . > $storage_tpm"/Temporales/count_tpm"
      echo "Done TPM"
      paste $storage_tpm"/Temporales/making_tpm" $storage_tpm"/Temporales/count_tpm" >> $storage_tpm"/Temporales/add_tpm"
      mv $storage_tpm"/Temporales/add_tpm" $storage_tpm"/Temporales/making_tpm"

      extr=$(($extr + 1))
    done
    cp $storage_tpm"/Temporales/making_tpm" $storage_tpm"/"$grupo"_"$spe"_TPM.counts"
  done
}

mapeo () {
  for spe in $Especies
  do
    mkdir $work_directory"/"$Result_Storage"/Mediciones_Expresion/Reference"
    read_files=$(ls $grupo"/"$spe"/"*"_work_read_1P.gz" | sed "s/.*\///" | sed "s/_work_read_1P.gz//")

    genome_file=$(ls $grupo"/"$spe"/"*"genomic.fa")
    gff_file=$(ls $grupo"/"$spe"/"*"annotations.gtf")

    echo "STAR Reference"
    STAR --runMode genomeGenerate --runThreadN $threads --genomeDir $work_directory"/"$Result_Storage"/Mediciones_Expresion/Reference" --genomeFastaFiles $genome_file --genomeChrBinNbits $STAR_genomeChrBinNbits--genomeSAindexNbases $STAR_genomeSAindexNbases --sjdbGTFfile $gff_file

    for read in $read_files
    do
      echo $read
      echo "Unzip 1P"
      seqkit seq $grupo"/"$spe"/"$read"_work_read_1P.gz" >> $work_directory"/"$Result_Storage"/Mediciones_Expresion/Read_temps/"$read"_work_read_1P"
      echo "Unzip 2P"
      seqkit seq $grupo"/"$spe"/"$read"_work_read_2P.gz" >> $work_directory"/"$Result_Storage"/Mediciones_Expresion/Read_temps/"$read"_work_read_2P"

      echo "Mapping"
      STAR --runMode alignReads --sjdbGTFfile $gff_file --runThreadN $threads --genomeSAindexNbases $STAR_genomeSAindexNbases --outFileNamePrefix $work_directory"/"$Result_Storage"/Mediciones_Expresion/Mapeos/"$grupo"_"$spe"_"$read"_" --outSAMtype BAM SortedByCoordinate --genomeDir $work_directory"/"$Result_Storage"/Mediciones_Expresion/Reference" --readFilesIn $work_directory"/"$Result_Storage"/Mediciones_Expresion/Read_temps/"$read"_work_read_1P" $work_directory"/"$Result_Storage"/Mediciones_Expresion/Read_temps/"$read"_work_read_2P"

      echo "Conteo"
      samtools index $work_directory"/"$Result_Storage"/Mediciones_Expresion/Mapeos/"$grupo"_"$spe"_"$read"_Aligned.sortedByCoord.out.bam"

      rm $work_directory"/"$Result_Storage"/Mediciones_Expresion/Read_temps/"$read"_work_read_1P"
      rm $work_directory"/"$Result_Storage"/Mediciones_Expresion/Read_temps/"$read"_work_read_2P"
    done
    rm -r $work_directory"/"$Result_Storage"/Mediciones_Expresion/Reference"
  done
}

check_hog_order () {
  gen_num=1
  for gen in $interval_genes
  do
    check_pho=$(grep -w -F -c $gen $PHO_File)
    if [ $check_pho -eq 0 ]
    then
      HOG=NA
    else
      HOG=$(grep -w -F $gen $PHO_File | awk -F "\t" '{print $1}')
    fi

    echo $location" "$HOG" "$gen" "$gen_num | tr " " "\t" >> $work_directory"/"$Result_Storage"/Context_Clusters_HOG/"$cluster_id"_hog.txt"
    gen_num=$(($gen_num + 1))
  done
}

check_gene_sl_function () {
  if [ "$wanted_gene" != "Missing" ] && [ $check_manual_review -eq 0 ]
  then
    check_gene_SL=$(grep -F -w $wanted_gene $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $6}')

    if [ $check_gene_SL == "SL" ]
    then
      wanted_gene=$(echo $wanted_gene"(SL)")
    fi
  fi
}

cluster_heatmap () {
  for cluster_id in $all_cluster
  do
    grep $search_term $work_directory"/"$Result_Storage"/Context_Clusters_HOG/"$cluster_id"_hog.txt" | awk -F "\t" '{print $2}' | grep -v -w NA | sort -u >> $work_directory"/"$Result_Storage"/Context_Clusters_HOG/Temp/"$cluster_id"_"$search_term".temp"
    total_terms=$(grep -c . $work_directory"/"$Result_Storage"/Context_Clusters_HOG/Temp/"$cluster_id"_"$search_term".temp")
    header=$(echo $header" "$cluster_id"("$total_terms")")
  done

  echo $header | tr " " "\t" >> $work_directory"/"$Result_Storage"/Context_Clusters_HOG/Cluster_heatmap_"$file_name".txt"

  for cluster1 in $all_cluster
  do
    total_terms=$(grep -c . $work_directory"/"$Result_Storage"/Context_Clusters_HOG/Temp/"$cluster1"_"$search_term".temp")
    registro=$(echo $cluster1"("$total_terms")")
    for cluster2 in $all_cluster
    do
      if [ $cluster1 == $cluster2 ]
      then
        registro=$(echo $registro" NA")
      else
        count=$(cat $work_directory"/"$Result_Storage"/Context_Clusters_HOG/Temp/"$cluster1"_"$search_term".temp" $work_directory"/"$Result_Storage"/Context_Clusters_HOG/Temp/"$cluster2"_"$search_term".temp" | sort | uniq -d | grep -c .)
        registro=$(echo $registro" "$count)
      fi
    done
    echo $registro | tr " " "\t" >> $work_directory"/"$Result_Storage"/Context_Clusters_HOG/Cluster_heatmap_"$file_name".txt"
  done
}

mapeo_second_pass () {
  mkdir $work_directory"/"$Result_Storage"/Mapeos_Second_Pass"
  mkdir $work_directory"/"$Result_Storage"/Mapeos_Second_Pass/Read_temps"

  for spe in $Especies
  do
    sj_pass=$(ls $work_directory"/"$Result_Storage"/Mediciones_Expresion/Mapeos/"$grupo"_"$spe"_"*"_SJ.out.tab" | sort)

    mkdir $work_directory"/"$Result_Storage"/Mapeos_Second_Pass/Reference"
    read_files=$(ls $grupo"/"$spe"/"*"_work_read_1P.gz" | sed "s/.*\///" | sed "s/_work_read_1P.gz//")

    genome_file=$(ls $grupo"/"$spe"/"*"genomic.fa")
    gff_file=$(ls $grupo"/"$spe"/"*"annotations.gtf")

    echo "STAR Reference"
    STAR --runMode genomeGenerate --runThreadN $threads --genomeDir $work_directory"/"$Result_Storage"/Mapeos_Second_Pass/Reference" --genomeFastaFiles $genome_file --genomeChrBinNbits $STAR_genomeChrBinNbits --genomeSAindexNbases $STAR_genomeSAindexNbases --sjdbGTFfile $gff_file

    for read in $read_files
    do
      echo $read
      echo "Unzip 1P"
      seqkit seq $grupo"/"$spe"/"$read"_work_read_1P.gz" >> $work_directory"/"$Result_Storage"/Mapeos_Second_Pass/Read_temps/"$read"_work_read_1P"
      echo "Unzip 2P"
      seqkit seq $grupo"/"$spe"/"$read"_work_read_2P.gz" >> $work_directory"/"$Result_Storage"/Mapeos_Second_Pass/Read_temps/"$read"_work_read_2P"

      echo "Mapping"
      STAR --runMode alignReads --sjdbFileChrStartEnd $sj_pass --sjdbGTFfile $gff_file --runThreadN $threads --genomeSAindexNbases $STAR_genomeSAindexNbases --outFileNamePrefix $work_directory"/"$Result_Storage"/Mapeos_Second_Pass/"$spe"_"$read"_" --outSAMtype BAM SortedByCoordinate --genomeDir $work_directory"/"$Result_Storage"/Mapeos_Second_Pass/Reference" --readFilesIn $work_directory"/"$Result_Storage"/Mapeos_Second_Pass/Read_temps/"$read"_work_read_1P" $work_directory"/"$Result_Storage"/Mapeos_Second_Pass/Read_temps/"$read"_work_read_2P"

      rm $work_directory"/"$Result_Storage"/Mapeos_Second_Pass/Read_temps/"$read"_work_read_1P"
      rm $work_directory"/"$Result_Storage"/Mapeos_Second_Pass/Read_temps/"$read"_work_read_2P"
    done
    rm -r $work_directory"/"$Result_Storage"/Mapeos_Second_Pass/Reference"
  done
}

new_GTF () {
  for spe in $Especies
  do
    gff_file=$(ls $grupo"/"$spe"/"*"annotations.gtf")

    check_bastard=$(echo $NCBI_GFF_format | grep -c $spe )

    echo $spe": "$gff_file
    echo "Bastard: "$check_bastard

    if [ $check_bastard -eq 0 ]
    then
      grep -v "#" $gff_file | awk -F "\t" '{if ($3=="transcript") print}' > $work_directory"/"$Result_Storage"/Modified_GTF_Files/Temporales/current.gtf"

      count=1
      recorrer=$(grep -c . $work_directory"/"$Result_Storage"/Modified_GTF_Files/Temporales/current.gtf")
      echo $spe": "$count" --- "$recorrer
      while [ $count -le $recorrer ]
      do
        transcript=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Modified_GTF_Files/Temporales/current.gtf" | awk -F "\t" '{print $9}' | tr ";" "\n" | grep transcript_id | tr -d '"' | sed "s/transcript_id transcript://")
        internal_new_GTF
        count=$(($count + 1))
      done
    else
      grep -v "#" $gff_file | awk -F "\t" '{if ($3=="transcript") print}' > $work_directory"/"$Result_Storage"/Modified_GTF_Files/Temporales/current.gtf"
      count=1
      recorrer=$(grep -c . $work_directory"/"$Result_Storage"/Modified_GTF_Files/Temporales/current.gtf")

      while [ $count -le $recorrer ]
      do
        transcript=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Modified_GTF_Files/Temporales/current.gtf" | awk -F "\t" '{print $9}' | tr ";" "\n" | grep gene_id | tr -d '"' | sed "s/gene_id //" | tr -d ' ' | sed "s/^gene-//")
        internal_new_GTF
        count=$(($count + 1))
      done
    fi
  done
}

internal_new_GTF () {
  check_standard_transcript=$(grep -c -F -w ">"$transcript $work_directory"/"$Result_Storage"/"$spe".pep")
  check_chimera_transcript=$(grep -c -F ">"$transcript"_Chimera" $work_directory"/"$Result_Storage"/"$spe".pep")
  echo "#############################################"
  echo $transcript "|| "$check_standard_transcript" || "$check_chimera_transcript
  echo
  echo "#############################################"
  if [ $check_standard_transcript -gt 0 ]
  then
    information=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Modified_GTF_Files/Temporales/current.gtf" | awk -F "\t" '{print $1,$2,$3,$4,$5,$6,$7,$8}')
    echo $information" ID="$transcript | tr " " "\t" >> $work_directory"/"$Result_Storage"/Modified_GTF_Files/"$spe"_with_Chimera_genes.gtf"
  elif [ $check_chimera_transcript -gt 0 ]
  then
    information=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Modified_GTF_Files/Temporales/current.gtf" | awk -F "\t" '{print $1,$2,$3}')
    gene_start=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Modified_GTF_Files/Temporales/current.gtf" | awk -F "\t" '{print $4}')
    gene_end=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Modified_GTF_Files/Temporales/current.gtf" | awk -F "\t" '{print $5}')
    gene_strand=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Modified_GTF_Files/Temporales/current.gtf" | awk -F "\t" '{print $7}')

    grep -F $transcript"_Chimera" $work_directory"/"$Result_Storage"/"$spe".pep" | sed "s/>//" > $work_directory"/"$Result_Storage"/Modified_GTF_Files/Temporales/current_chimera.temp"

    if [ $gene_strand == "+" ]
    then
      repleace_start=$gene_start
      repleace_end=$gene_end
      move_the_start=$(echo "1")
    else
      repleace_start=$gene_end
      repleace_end=$gene_start
      move_the_start=$(echo "-1")
    fi

    count_chimera=1
    recorrer_chimera=$(grep -c . $work_directory"/"$Result_Storage"/Modified_GTF_Files/Temporales/current_chimera.temp")

    while [ $count_chimera -le $recorrer_chimera ]
    do
      chimera_id=$(sed -n $count_chimera"p" $work_directory"/"$Result_Storage"/Modified_GTF_Files/Temporales/current_chimera.temp")

      check_coord_start=$(echo $chimera_id | grep -c "Start")

      if [ $check_coord_start -gt 0 ]
      then
        coordinates=$(echo $chimera_id | sed "s/__/\n/" | tail -n 1 | sed "s/Start/$repleace_start/" | tr "_" "\n" | sort -n | tr "\n" " " | sed "s/ $/\n/")
      else
        coordinates=$(echo $chimera_id | sed "s/__/\n/" | tail -n 1 | awk -F "_" -v move=$move_the_start '{print $1+move"_"$2 }' | sed "s/End/$repleace_end/"  | tr "_" "\n" | sort -n | tr "\n" " " | sed "s/ $/\n/")
      fi
      echo
      echo $information" "$coordinates" . "$gene_strand" . ID="$chimera_id
      echo $information" "$coordinates" . "$gene_strand" . ID="$chimera_id | tr " " "\t" >> $work_directory"/"$Result_Storage"/Modified_GTF_Files/"$spe"_with_Chimera_genes.gtf"
      count_chimera=$(($count_chimera + 1))
    done
  fi
}

check_first_chimera () {
  # Verificando si el Chimera1 cuenta
  if [ $strand == "+" ]
  then
    include_first_gene=$(grep -w -F $gen $work_directory"/SL_Read_Counts/"$spe"_tier3.counts" | awk -F "\t" -v ace=$ace '{if ($2<ace) print}' | grep -c .)
  else
    include_first_gene=$(grep -w -F $gen $work_directory"/SL_Read_Counts/"$spe"_tier3.counts" | awk -F "\t" -v ace=$ace '{if ($2>ace) print}' | grep -c .)
  fi
}

check_operon_info () {
  if [ $check_artifact -gt 0 ]
  then
    print_gen=$(echo $print_gen"[C]")
  fi

  check_SL=$(sed -n $recorrer"p" $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | grep -w -c "SL")

  if [ $check_SL -gt 0 ]
  then
    print_gen=$(echo $print_gen"*")
  fi

  if [ $op_strand == "+" ]
  then
    if [ $recorrer -lt $last_op_line  ]
    then
      registro_genes=$(echo $registro_genes")_"$print_gen"_("$print_dist)
    else
      registro_genes=$(echo $registro_genes")_"$print_gen)
    fi
  else
    if [ $recorrer -eq $last_op_line ]
    then
      registro_genes=$(echo $print_gen"_(")
    elif [ $recorrer -gt $first_op_line  ]
    then
      registro_genes=$(echo $registro_genes$print_dist")_"$print_gen"_(")
    else
      registro_genes=$(echo $registro_genes$print_dist")_"$print_gen)
    fi
  fi

  registro_pho=$(echo $registro_pho"__"$print_pho)
}

find_operons () {
  for spe in $Especies
  do
    if [ ! -f $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" ]
    then
      read_files=$(ls $grupo"/"$spe"/"*"_work_read_1P.gz" | sed "s/.*\///" | sed "s/_work_read_1P.gz//")
      N_read_files=$(echo $read_files | tr " " "\n" | grep -c .)

      gff_file=$(ls $grupo"/"$spe"/"*"annotations.gtf")
      chromosomes=$(grep -v "#" $gff_file | awk -F "\t" '{if ($3=="transcript") print $1}' | grep -v mitochondrion | grep -v _MITO$ | sort -u)

      echo "GenID Chromosome Start End Strand SL_Insertion Chimeric Dist_Next_Gene Read_Coverage OperonID PHO" | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt"
      check_trouble=$(echo "Sparganum_proliferum Paragonimus_heterotremus Fasciolopsis_buski" | grep -c $spe)
      check_trouble2=$(echo "Taenia_multiceps Schistosoma_mansoni" | grep -c $spe)
      check_trouble3=$(echo "Paragonimus_heterotremus Fasciolopsis_buski" | grep -c $spe)
      operon_ID=1

      for chr in $chromosomes
      do
        echo "Iniciando: "$spe" "$chr
        grep -F -w $chr $gff_file | awk -F "\t" '{if ($3=="transcript") print}' > $work_directory"/"$Result_Storage"/Operons/Temporal/temp_info.tab"
        grep -F -w $chr $gff_file | awk -F "\t" '{if ($3=="CDS") print}' > $work_directory"/"$Result_Storage"/Operons/Temporal/check_cds_gene.temp"

        check_gene_in_chr=$(awk -F "\t" '{print $3}' $work_directory"/"$Result_Storage"/Operons/Temporal/check_cds_gene.temp" | grep -w -c CDS)

        if [ $check_gene_in_chr -eq 0 ]
        then
          echo $grupo" "$spe" "$chr": does not have genes" >>  $work_directory"/"$Result_Storage"/Operons/Warnings.txt"
        else
          if [ $check_trouble -gt 0 ]
          then
            genes=$(awk -F "\t" '{print $9}' $work_directory"/"$Result_Storage"/Operons/Temporal/temp_info.tab" | tr ";" "\n" | grep gene_id | sed 's/gene_id "//' | sed "s/^gene-//" | tr -d '"' | tr -d ' ' | sort -u)
          elif [ $check_trouble2 -gt 0 ]
          then
            genes=$(awk -F "\t" '{print $9}' $work_directory"/"$Result_Storage"/Operons/Temporal/temp_info.tab" | tr ";" "\n" | grep gene_id | sed 's/gene_id "gene://' | sed 's/gene_id "transcript://' | tr -d '"' | tr -d ' ' | sort -u)
          else
            genes=$(awk -F "\t" '{print $9}' $work_directory"/"$Result_Storage"/Operons/Temporal/temp_info.tab" | tr ";" "\n" | grep gene_id | sed 's/gene_id "gene://' | tr -d '"' | tr -d ' ' | sort -u)
          fi

          for gen in $genes
          do
            check_protein_coding_gene=$(grep -F $gen'";' $work_directory"/"$Result_Storage"/Operons/Temporal/check_cds_gene.temp" | awk -F "\t" '{print $3}' | grep -c CDS )
            if [ $check_protein_coding_gene -gt 0 ]
            then
              gen_start=$(grep -F $gen'";' $work_directory"/"$Result_Storage"/Operons/Temporal/temp_info.tab" | awk -F "\t" '{print $4}' | sort -n -u | head -n 1)
              gen_end=$(grep -F $gen'";' $work_directory"/"$Result_Storage"/Operons/Temporal/temp_info.tab" | awk -F "\t" '{print $5}' | sort -n -u | tail -n 1)
              gen_strand=$(grep -F $gen'";' $work_directory"/"$Result_Storage"/Operons/Temporal/temp_info.tab" | awk -F "\t" '{print $7}' | sort -u)
              echo $gen" "$chr" "$gen_start" "$gen_end" "$gen_strand | tr " " "\t"  >> $work_directory"/"$Result_Storage"/Operons/Temporal/"$grupo"_"$spe"_gene.temp"
            fi
          done

          sort -k3,3n -k4,4r -s $work_directory"/"$Result_Storage"/Operons/Temporal/"$grupo"_"$spe"_gene.temp" > $work_directory"/"$Result_Storage"/Operons/Temporal/"$grupo"_"$spe"_gene.temp2"
          rm $work_directory"/"$Result_Storage"/Operons/Temporal/"$grupo"_"$spe"_gene.temp"

          count=1
          recorrer=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Temporal/"$grupo"_"$spe"_gene.temp2")
          continue_operon=FALSE
          member=NA

          while [ $count -le $recorrer ]
          do
            chimera_found=FALSE
            pho=NA
            read_coverage=NA

            if [ $check_trouble3 -gt 0 ]
            then
              info_line=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Operons/Temporal/"$grupo"_"$spe"_gene.temp2" | sed "s/^gene-//")
            else
              info_line=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Operons/Temporal/"$grupo"_"$spe"_gene.temp2")
            fi
            gene_id=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Operons/Temporal/"$grupo"_"$spe"_gene.temp2" | awk -F "\t" '{print $1}')
            transcript_id=$(grep -w $gene_id'$' $work_directory"/"$Result_Storage"/Mediciones_Expresion/Background_TPM/"$spe"_per_gene.txt" | grep -F $gene_id | awk -F "\t" '{print $1}')

            check_SL_Gene=$(grep -w Concordante $work_directory"/SL_Read_Counts/"$spe"_tier3.counts" | grep -w -v Trans_End | grep -c -w -F $gene_id)
            if [ $check_SL_Gene -gt 0 ]
            then
              SL_Gene=SL
            else
              SL_Gene=X
            fi

            if [ $check_trouble -gt 0 ]
            then
              if [ "$spe" == "Sparganum_proliferum" ]
              then
                grep $gene_id'";' $work_directory"/"$Result_Storage"/Operons/Temporal/temp_info.tab" | awk -F "\t" '{print $9}' | tr ";" "\n" | grep gene_id | sed 's/gene_id "//' | tr -d '"' | tr -d ' ' | sort -u  > $work_directory"/"$Result_Storage"/Operons/Temporal/"$grupo"_"$spe"_transcript.temp"
              else
                grep $gene_id'";' $work_directory"/"$Result_Storage"/Operons/Temporal/temp_info.tab" | awk -F "\t" '{print $9}' | tr ";" "\n" | grep transcript_id | sed 's/transcript_id "//' | tr -d '"' | tr -d ' ' | sort -u | awk -F "|" '{print $3}' > $work_directory"/"$Result_Storage"/Operons/Temporal/"$grupo"_"$spe"_transcript.temp"
              fi
            else
              grep $gene_id'";' $work_directory"/"$Result_Storage"/Operons/Temporal/temp_info.tab" | awk -F "\t" '{print $9}' | tr ";" "\n" | grep transcript_id | sed 's/transcript_id "transcript://' | tr -d '"' | tr -d ' ' | sort -u  > $work_directory"/"$Result_Storage"/Operons/Temporal/"$grupo"_"$spe"_transcript.temp"
            fi

            awk '{print $1"_Chimera"}' $work_directory"/"$Result_Storage"/Operons/Temporal/"$grupo"_"$spe"_transcript.temp" > $work_directory"/"$Result_Storage"/Operons/Temporal/"$grupo"_"$spe"_Chimera.temp"

            artifact=$(grep -c -F -f $work_directory"/"$Result_Storage"/Operons/Temporal/"$grupo"_"$spe"_Chimera.temp" $work_directory"/"$Result_Storage"/"$spe".pep" | awk '{if ($1>=1) print "Chimeric"; else print "NA"}' )
            if [ $artifact == "Chimeric" ]
            then
              chimeric_ids=$(grep -F -f $work_directory"/"$Result_Storage"/Operons/Temporal/"$grupo"_"$spe"_Chimera.temp" $work_directory"/"$Result_Storage"/"$spe".pep" | sed "s/>//"  | sed "s/__.*//" | sort)
            fi

            start_feature=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Operons/Temporal/"$grupo"_"$spe"_gene.temp2" | awk -F "\t" '{print $3}')
            end_feature=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Operons/Temporal/"$grupo"_"$spe"_gene.temp2" | awk -F "\t" '{print $4}')
            strand_feature=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Operons/Temporal/"$grupo"_"$spe"_gene.temp2" | awk -F "\t" '{print $5}')

            extract=$(($count + 1))
            if [ $extract -le $recorrer ]
            then
              use_this_distance=FALSE

              while [ "$use_this_distance" == "FALSE" ]
              do
                next_feature_ID=$(sed -n $extract"p" $work_directory"/"$Result_Storage"/Operons/Temporal/"$grupo"_"$spe"_gene.temp2" | awk -F "\t" '{print $1}')
                next_feature_start=$(sed -n $extract"p" $work_directory"/"$Result_Storage"/Operons/Temporal/"$grupo"_"$spe"_gene.temp2" | awk -F "\t" '{print $3}')
                next_feature_end=$(sed -n $extract"p" $work_directory"/"$Result_Storage"/Operons/Temporal/"$grupo"_"$spe"_gene.temp2" | awk -F "\t" '{print $4}')

                next_feature_strand=$(sed -n $extract"p" $work_directory"/"$Result_Storage"/Operons/Temporal/"$grupo"_"$spe"_gene.temp2" | awk -F "\t" '{print $5}')

                distance=$(($next_feature_start - $end_feature))

                if [ $distance -gt 0 ]
                then
                  use_this_distance=TRUE
                else
                  if [ $next_feature_end -gt $end_feature ]
                  then
                    first_gene_lenght=$(($end_feature - $start_feature + 1))
                    second_gene_lenght=$(($next_feature_end - $next_feature_start + 1))
                    overlap_lenght=$(($end_feature - $next_feature_start + 1))

                    test_overlap_first=$(echo $first_gene_lenght $overlap_lenght | awk -F " " '{if ($2/$1>=0.9) print "SKIP"; else print "WORK"}' )
                    test_overlap_second=$(echo $second_gene_lenght $overlap_lenght | awk -F " " '{if ($2/$1>=0.9) print "SKIP"; else print "WORK"}' )

                    if [ "$test_overlap_first" == "WORK" ] || [ "$test_overlap_second" == "WORK" ]
                    then
                      use_this_distance=TRUE
                    else
                      skip_because=Large_Overlap
                      skip_gene_overlap_shenanigans
                    fi
                  else
                    skip_because=Complete_Overlap
                    skip_gene_overlap_shenanigans
                  fi
                fi
              done
            else
              distance=End
            fi

            if [ "$distance" == "End" ]
            then
              if [ "$continue_operon" == "TRUE" ]
              then
                member=$(echo "Operon-"$operon_ID)
                operon_ID=$(($operon_ID + 1))
                continue_operon=FALSE
              elif [ $artifact == "Chimeric" ]
              then
                member=$(echo "Operon-"$operon_ID)
                read_coverage=Chimeric
                continue_operon=FALSE
                operon_ID=$(($operon_ID + 1))
              fi
            elif [ "$strand_feature" != "$next_feature_strand" ]
            then
              if [ $artifact == "Chimeric" ]
              then
                member=$(echo "Operon-"$operon_ID)
                read_coverage=Chimeric
                continue_operon=TRUE
              fi

              if [ "$continue_operon" == "TRUE" ]
              then
                operon_ID=$(($operon_ID + 1))
                continue_operon=FALSE
              fi
            elif [ $artifact == "Chimeric" ] && [ $continue_operon == "FALSE" ]
            then
              member=$(echo "Operon-"$operon_ID)
              read_coverage=Chimeric
              if [ $distance -le $operon_range ]
              then
                continue_operon=TRUE
              else
                operon_ID=$(($operon_ID + 1))
              fi
            elif [ $distance -le $operon_range ] && [ $continue_operon == "FALSE" ]
            then
              member=$(echo "Operon-"$operon_ID)
              continue_operon=TRUE
              if [ $distance -gt 1 ]
              then
                check_read_depth
              elif [ $distance -eq 1 ]
              then
                read_coverage=Adjacent
              else
                read_coverage=Overlap
              fi
            elif [ $distance -gt $operon_range ] && [ $continue_operon == "TRUE" ]
            then
              continue_operon=FALSE
              operon_ID=$(($operon_ID + 1))
            fi

            check_orthogroups
            echo $info_line" "$SL_Gene" "$artifact" "$distance" "$read_coverage" "$member" "$pho | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt"
            if [ $continue_operon == "FALSE" ]
            then
              member=NA
            fi
            count=$(($count + 1))
          done
        fi
      done

      # Make Operon Final Table
      echo "OpID Chromosome Strand Start End Genes PHOs Gene_Overlaps" | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/"$spe"_operon.txt"
      all_operons=$(grep -v -w "GenID" $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $10}' | grep Operon | awk -F "-" '{print $2}' | sort -u | sort -n)
      for ope in $all_operons
      do
        affected_by_Overlaps=$(grep -w -n "Operon-"$ope $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | grep -w -c Overlap | awk '{if ($1>0) print "X"; else print "."}' )

        echo "Operones: "$spe" "$ope
        first_op_line=$(grep -w -n "Operon-"$ope $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F ":" '{print $1}' | head -n 1)
        last_op_line=$(grep -w -n "Operon-"$ope $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F ":" '{print $1}' | tail -n 1)

        op_strand=$(sed -n $first_op_line"p" $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $5}')
        op_chromosome=$(sed -n $first_op_line"p" $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $2}')
        start_op=$(sed -n $first_op_line"p" $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $3}')
        end_op=$(sed -n $last_op_line"p" $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $4}')

        if [ $op_strand == "+" ]
        then
          recorrer=$first_op_line
          registro_genes=
          registro_pho=

          while [ $recorrer -le $last_op_line ]
          do
            print_gen=$(sed -n $recorrer"p" $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $1}')
            print_dist=$(sed -n $recorrer"p" $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $8}')
            check_artifact=$(sed -n $recorrer"p" $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $7}' | grep -c "Chimeric")
            print_pho=$(sed -n $recorrer"p" $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $11}')

            check_operon_info
            recorrer=$(($recorrer + 1))
          done
        else
          recorrer=$last_op_line
          registro_genes=
          registro_pho=

          while [ $recorrer -ge $first_op_line ]
          do
            print_gen=$(sed -n $recorrer"p" $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $1}')
            print_dist=$(sed -n $recorrer"p" $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $8}')
            check_artifact=$(sed -n $recorrer"p" $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $7}' | grep -c "Chimeric")
            print_pho=$(sed -n $recorrer"p" $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $11}')

            check_operon_info
            recorrer=$(($recorrer - 1))
          done
        fi
        registro_genes=$(echo $registro_genes | sed "s/^)_//" | sed 's/_($//')
        registro_pho=$(echo $registro_pho | sed "s/^__//")

        echo "Operon-"$ope" "$op_chromosome" "$op_strand" "$start_op" "$end_op" "$registro_genes" "$registro_pho" "$affected_by_Overlaps | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/"$spe"_operon.txt"
      done
    fi
  done
}

sl_background () {
  for spe in $Especies
  do
    check_points=$(echo $annoying_points | grep -c $spe)

    if [ $check_points -gt 0 ]
    then
      awk -F "\t" '{print $1}' $storage_tpm"/"$grupo"_"$spe"_TPM.counts" | grep -v GenID | sed "s/_Chimera.*//" | sort -u | awk -F "." '{print $1"\t"$0}' > $work_directory"/"$Result_Storage"/Mediciones_Expresion/Background_TPM/Temporales/current_spe.genes"
    elif [ "$spe" == "Mesocestoides_corti" ] || [ "$spe" == "Spirometra_erinaceieuropaei" ] || [ "$spe" == "Schistocephalus_solidus" ] || [ "$spe" == "Trichobilharzia_regenti" ]
    then
      awk -F "\t" '{print $1}' $storage_tpm"/"$grupo"_"$spe"_TPM.counts" |  grep -v GenID | sed "s/_Chimera.*//" | sort -u | awk -F "-" '{print $1"\t"$0}' > $work_directory"/"$Result_Storage"/Mediciones_Expresion/Background_TPM/Temporales/current_spe.genes"
    elif [ "$spe" == "Taenia_saginata" ] || [ "$spe" == "Taenia_asiatica" ]
    then
      awk -F "\t" '{print $1}' $storage_tpm"/"$grupo"_"$spe"_TPM.counts" |  grep -v GenID | sed "s/_Chimera.*//" | sort -u | awk -F "m" '{print $1"\t"$0}' > $work_directory"/"$Result_Storage"/Mediciones_Expresion/Background_TPM/Temporales/current_spe.genes"
    else
      awk -F "\t" '{print $1}' $storage_tpm"/"$grupo"_"$spe"_TPM.counts" | grep -v GenID | sed "s/_Chimera.*//" | sort -u | awk -F "\t" '{print $1"\t"$1}' | sed "s/-mRNA-1\t/\t/" > $work_directory"/"$Result_Storage"/Mediciones_Expresion/Background_TPM/Temporales/current_spe.genes"
    fi

    count=1
    recorrer=$(grep -c . $work_directory"/"$Result_Storage"/Mediciones_Expresion/Background_TPM/Temporales/current_spe.genes")

    while [ $count -le $recorrer ]
    do
      gen=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Mediciones_Expresion/Background_TPM/Temporales/current_spe.genes" | awk -F "\t" '{print $1}')
      transcript=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Mediciones_Expresion/Background_TPM/Temporales/current_spe.genes" | awk -F "\t" '{print $2}')

      check_sl_gene=$(grep -w Concordante $work_directory"/SL_Read_Counts/"$spe"_tier3.counts" | grep -w -v Trans_End | grep -w -F -c $gen )
      check_chimera=$(grep -c $transcript"_Chimera" $storage_tpm"/"$grupo"_"$spe"_TPM.counts")

      echo $spe" || "$gen" || "$transcript" || "$check_sl_gene" || "$check_chimera

      if [ $check_sl_gene -gt 0 ]
      then
        if [ $check_chimera -eq 0 ]
        then
          average=$(grep -w -F $transcript $storage_tpm"/"$grupo"_"$spe"_TPM.counts" | tr "\t" "\n" | sed 1d | awk '{ sum += $1 } END { if (NR > 0) print sum / NR}')
          echo $transcript" "$average" "$gen | tr " " "\t" >> $work_directory"/"$Result_Storage"/Mediciones_Expresion/Background_TPM/"$spe"_per_gene.txt"
        else
          ace=$(grep -F $transcript"_Chimera_1__Start_"  $storage_tpm"/"$grupo"_"$spe"_TPM.counts" | awk -F "\t" '{print $1}' | sed "s/.*_Chimera_1__Start_//")
          strand=$(grep -w -F $gen $work_directory"/SL_Read_Counts/"$spe"_tier3.counts" | awk -F "\t" '{print $4}' | sort -u)

          check_first_chimera

          if [ $include_first_gene -gt 0 ]
          then
            chimeric_intervals=$(grep -F $transcript"_Chimera" $storage_tpm"/"$grupo"_"$spe"_TPM.counts" | awk -F "\t" '{print $1}')
          else
            chimeric_intervals=$(grep -F $transcript"_Chimera" $storage_tpm"/"$grupo"_"$spe"_TPM.counts" | awk -F "\t" '{print $1}' | grep -v "_Chimera_1")
          fi

          for chimera in $chimeric_intervals
          do
            average=$(grep -w -F $chimera $storage_tpm"/"$grupo"_"$spe"_TPM.counts" | tr "\t" "\n" | sed 1d | awk '{ sum += $1 } END { if (NR > 0) print sum / NR}')
            echo $chimera" "$average" "$gen | tr " " "\t" >> $work_directory"/"$Result_Storage"/Mediciones_Expresion/Background_TPM/"$spe"_per_gene.txt"
          done
        fi
      fi
      count=$(($count + 1))
    done

    species_SL_median=$(awk -F "\t" '{print $2}' $work_directory"/"$Result_Storage"/Mediciones_Expresion/Background_TPM/"$spe"_per_gene.txt" | tr "." "," | datamash median 1 | tr "," ".")
    echo $spe" "$species_SL_median | tr " " "\t" >> $work_directory"/"$Result_Storage"/Mediciones_Expresion/Background_TPM/SL_Reference_median.txt"
  done
}

adjust_gene_ids () {
  print_new_gene_ID=$new_gene_ID
  if [ "$new_gene_chimera" == "Chimeric" ]
  then
    print_new_gene_ID=$(echo $print_new_gene_ID"[C]")
  fi

  if [ "$new_gene_sl" == "SL" ]
  then
    print_new_gene_ID=$(echo $print_new_gene_ID"*")
  fi
}

presence_absence (){
  echo "Pair Total Pair_Confirm "$Especies" N_Grupos_Cestoda N_Grupos_Cestoda_Derivados N_Grupos_Cestoda_Basal N_Grupos_Taenida N_Grupos_Hymenolepis N_Grupos_Trematoda N_Grupos_Schistosmidae N_Grupos_Trematoda_No_schistosoma N_Grupos_Fasciolidae N_Grupos_Trematoda_NoSchistoFasciola N_Grupos_Paragonimus N_Grupos_Schistosoma N_Grupos_Special" | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/"$presence_absence_results
  for pair in $pho_unit_of_interest
  do
    echo "checking: "$pair
    N_Grupos_Cestoda=0
    N_Grupos_Cestoda_Derivados=0
    N_Grupos_Cestoda_Basal=0
    N_Grupos_Taenida=0
    N_Grupos_Hymenolepis=0
    N_Grupos_Trematoda=0
    N_Grupos_Schistosmidae=0
    N_Grupos_Trematoda_No_schistosoma=0
    N_Grupos_Fasciolidae=0
    N_Grupos_Trematoda_NoSchistoFasciola=0
    N_Grupos_Paragonimus=0
    N_Grupos_Schistosoma=0
    N_Grupos_Special=0
    Total=0
    registro_conteos=$pair

    for spe in $Especies
    do
      check_pair_on_spe=$(grep -F -c $pair $work_directory"/"$Result_Storage"/Operons/"$spe"_operon.txt" | awk '{if ($1>0) print "X"; else print "."}')
      registro_conteos=$(echo $registro_conteos" "$check_pair_on_spe)

      if [ "$check_pair_on_spe" == "X" ]
      then
        Total=$(($Total + 1))

        # A que grupo pertenece este hit en particular?
        is_cestoda=$(echo $Grupos_Cestoda | grep -c $spe)
        if [ $is_cestoda -gt 0 ]
        then
          N_Grupos_Cestoda=$(($N_Grupos_Cestoda + 1))
        fi

        is_cestoda_derivados=$(echo $Grupos_Cestoda_Derivados | grep -c $spe)
        if [ $is_cestoda_derivados -gt 0 ]
        then
          N_Grupos_Cestoda_Derivados=$(($N_Grupos_Cestoda_Derivados + 1))
        fi

        is_cestoda_basal=$(echo $Grupos_Cestoda_Basal | grep -c $spe)
        if [ $is_cestoda_basal -gt 0 ]
        then
          N_Grupos_Cestoda_Basal=$(($N_Grupos_Cestoda_Basal + 1))
        fi

        is_taenida=$(echo $Grupos_Taenida | grep -c $spe)
        if [ $is_taenida -gt 0 ]
        then
          N_Grupos_Taenida=$(($N_Grupos_Taenida + 1))
        fi

        is_hymenolepis=$(echo $Grupos_Hymenolepis | grep -c $spe)
        if [ $is_taenida -gt 0 ]
        then
          N_Grupos_Hymenolepis=$(($N_Grupos_Hymenolepis + 1))
        fi

        is_trematoda=$(echo $Grupos_Trematoda | grep -c $spe)
        if [ $is_trematoda -gt 0 ]
        then
          N_Grupos_Trematoda=$(($N_Grupos_Trematoda + 1))
        fi

        is_schistosmidae=$(echo $Grupos_Schistosmidae | grep -c $spe)
        if [ $is_schistosmidae -gt 0 ]
        then
          N_Grupos_Schistosmidae=$(($N_Grupos_Schistosmidae + 1))
        fi

        is_not_schistosmidae=$(echo $Grupos_Trematoda_No_schistosoma | grep -c $spe)
        if [ $is_not_schistosmidae -gt 0 ]
        then
          N_Grupos_Trematoda_No_schistosoma=$(($N_Grupos_Trematoda_No_schistosoma + 1))
        fi

        is_fasciolidae=$(echo $Grupos_Fasciolidae | grep -c $spe)
        if [ $is_fasciolidae -gt 0 ]
        then
          N_Grupos_Fasciolidae=$(($N_Grupos_Fasciolidae + 1))
        fi

        is_other_trematoda=$(echo $Grupos_Trematoda_NoSchistoFasciola | grep -c $spe)
        if [ $is_other_trematoda -gt 0 ]
        then
          N_Grupos_Trematoda_NoSchistoFasciola=$(($N_Grupos_Trematoda_NoSchistoFasciola + 1))
        fi

        is_paragonimus=$(echo $Grupos_Paragonimus | grep -c $spe)
        if [ $is_paragonimus -gt 0 ]
        then
          N_Grupos_Paragonimus=$(($N_Grupos_Paragonimus + 1))
        fi

        is_schistosoma=$(echo $Grupos_Schistosoma | grep -c $spe)
        if [ $is_schistosoma -gt 0 ]
        then
          N_Grupos_Schistosoma=$(($N_Grupos_Schistosoma + 1))
        fi

        is_special=$(echo $Grupos_Special | grep -c $spe)
        if [ $is_schistosoma -gt 0 ]
        then
          N_Grupos_Special=$(($N_Grupos_Special + 1))
        fi

      fi
    done
    echo $pair" "$Total" "$registro_conteos" "$N_Grupos_Cestoda" "$N_Grupos_Cestoda_Derivados" "$N_Grupos_Cestoda_Basal" "$N_Grupos_Taenida" "$N_Grupos_Hymenolepis" "$N_Grupos_Trematoda" "$N_Grupos_Schistosmidae" "$N_Grupos_Trematoda_No_schistosoma" "$N_Grupos_Fasciolidae" "$N_Grupos_Trematoda_NoSchistoFasciola" "$N_Grupos_Paragonimus" "$N_Grupos_Schistosoma" "$N_Grupos_Special  | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/"$presence_absence_results
  done
}

gene_location_info (){
  if [ -f $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$pair"_"$intermediary_storage"_"$species_name"_gene_coordinates.txt" ]
  then
    rm $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$pair"_"$intermediary_storage"_"$species_name"_gene_coordinates.txt"
  fi

  for gen in $interest_genes
  do
    grep -w -F "ID="$gen $work_directory"/"$Result_Storage"/Modified_GTF_Files/"$species_name"_with_Chimera_genes.gtf" | awk -F "\t" '{print $9"\t"$1"\t"$4"\t"$5"\t"$7}' | sed "s/^ID=//" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$pair"_"$intermediary_storage"_"$species_name"_gene_coordinates.txt"
  done

  if [ ! -f $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$pair"_"$intermediary_storage"_"$species_name"_gene_coordinates.txt" ]
  then
    echo "" > $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$pair"_"$intermediary_storage"_"$species_name"_gene_coordinates.txt"
  fi
}

confirm_pair_synteny () {
  mkdir $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$result_storage"_details"
  echo "PHO_Pair Species Status Missing_PHOs Alternative_Operons_First Alternative_Operons_Second" | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/Summary_"$result_storage".txt"

  for pair in $interest_pairs
  do
    count_species=4 # Primer columna con especies

    echo $result_storage" Pair:"$pair
    first_pho=$(echo $pair | awk -F "__" '{print $1}')
    second_pho=$(echo $pair | awk -F "__" '{print $2}')

    # buscar otros pares de interes, en la especie
    echo "AVOID" > $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/first_pho_alternatives.temp"
    echo "AVOID" > $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/second_pho_alternatives.temp"

    awk -F "\t" '{print $2}' $work_directory"/"$Result_Storage"/Operons/Questions/Work_gene_pairs.txt" | grep $first_pho"_" | grep -v $pair >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/first_pho_alternatives.temp"
    awk -F "\t" '{print $2}' $work_directory"/"$Result_Storage"/Operons/Questions/Work_gene_pairs.txt" | grep "_"$second_pho | grep -v $pair >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/second_pho_alternatives.temp"

    while [ $count_species -le $recorrer_species ]
    do
      # Report Variables
      Missing_PHO=NA

      check_absence_presence=$(grep -w $pair $work_directory"/"$Result_Storage"/Operons/Questions/Presence_or_Absence_By_Pairs.txt" | awk -F "\t" -v N=$count_species '{print $N}'| grep -w -c "X")
      if [ $check_absence_presence -eq 0 ]
      then
        species_name=$(head -n 1 $work_directory"/"$Result_Storage"/Operons/Questions/Presence_or_Absence_By_Pairs.txt" | awk -F "\t" -v N=$count_species '{print $N}')
        spe_num=$(head -n 1 $PHO_File | tr "\t" "\n" | grep -n $species_name | awk -F ":" '{print $1}')
        check_first_alternative_pairs=$(grep -c -f $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/first_pho_alternatives.temp" $work_directory"/"$Result_Storage"/Operons/"$species_name"_operon.txt")

        if [ $check_first_alternative_pairs -gt 0 ]
        then
          first_alternative_pairs=$(grep -o -f $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/first_pho_alternatives.temp" $work_directory"/"$Result_Storage"/Operons/"$species_name"_operon.txt" | sort -u | tr "\n" ";" | sed "s/;$//")
        else
          first_alternative_pairs=NA
        fi

        check_second_alternative_pairs=$(grep -c -f $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/second_pho_alternatives.temp" $work_directory"/"$Result_Storage"/Operons/"$species_name"_operon.txt")
        if [ $check_second_alternative_pairs -gt 0 ]
        then
          second_alternative_pairs=$(grep -o -f $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/second_pho_alternatives.temp" $work_directory"/"$Result_Storage"/Operons/"$species_name"_operon.txt" | sort -u | tr "\n" ";" | sed "s/;$//")
        else
          second_alternative_pairs=NA
        fi

        More_than_one_Gene=FALSE
        echo $pair" --- "$species_name"("$spe_num")"

        interest_genes=$(grep -w $first_pho $PHO_File | awk -F "\t" -v N=$spe_num '{print $N}' | tr -d " " | tr "," "\n")
        intermediary_storage=First_PHO
        gene_location_info

        interest_genes=$(grep -w $second_pho $PHO_File | awk -F "\t" -v N=$spe_num '{print $N}' | tr -d " " | tr "," "\n")
        intermediary_storage=Second_PHO
        gene_location_info

        check_first_pho_data=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$pair"_First_PHO_"$species_name"_gene_coordinates.txt")
        check_second_pho_data=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$pair"_Second_PHO_"$species_name"_gene_coordinates.txt")

        # Control de sintenia: 1- El PHO esta presente en la especie
        if [ $check_first_pho_data -eq 0 ] || [ $check_second_pho_data -eq 0 ]
        then
          Status=Missing_Genes
          if [ $check_first_pho_data -eq 0 ]
          then
            Missing_PHO=$(echo $Missing_PHO"_"$first_pho)
          fi

          if [ $check_second_pho_data -eq 0 ]
          then
            Missing_PHO=$(echo $Missing_PHO"_"$second_pho)
          fi
          Missing_PHO=$(echo $Missing_PHO | sed "s/NA_//" | tr "_" "\n" | sort -u | tr "\n" "_" | sed "s/_$//")
        else
          # Control de sintenia: 2- Hay representantes del PHO en el mismo cromosoma
          awk -F "\t" '{print $2}' $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$pair"_First_PHO_"$species_name"_gene_coordinates.txt" | sort -u > $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/temp_first_chr"
          check_same_chromosome=$(grep -c -w -F -f $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/temp_first_chr" $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$pair"_Second_PHO_"$species_name"_gene_coordinates.txt")

          if [ $check_same_chromosome -eq 0 ]
          then
            Status=Dif_Chr
          else
            # Registro genes par a par
            matching_chr=$(grep -w -F -f $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/temp_first_chr" $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$pair"_Second_PHO_"$species_name"_gene_coordinates.txt" | awk -F "\t" '{print $2}' | sort -u)
            echo "Pair Chr Gen1 Gen2 Internal_Status N_interval_genes Genes_Internos" | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$result_storage"_details/"$pair"_"$species_name"_"$result_storage"_Gene_Details.txt"
            for chr in $matching_chr
            do
              first_pho_genes=$(grep -w -F $chr $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$pair"_First_PHO_"$species_name"_gene_coordinates.txt" | awk -F "\t" '{print $1}')
              second_pho_genes=$(grep -w -F $chr $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$pair"_Second_PHO_"$species_name"_gene_coordinates.txt" | awk -F "\t" '{print $1}')

              for gen1 in $first_pho_genes
              do
                gen1_start=$(grep -w -F $gen1 $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$pair"_First_PHO_"$species_name"_gene_coordinates.txt" | awk -F "\t" '{print $3}')
                gen1_end=$(grep -w -F $gen1 $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$pair"_First_PHO_"$species_name"_gene_coordinates.txt" | awk -F "\t" '{print $4}')
                gen1_strand=$(grep -w -F $gen1 $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$pair"_First_PHO_"$species_name"_gene_coordinates.txt" | awk -F "\t" '{print $5}')

                for gen2 in $second_pho_genes
                do
                  if [ $gen1 != $gen2 ]
                  then
                    More_than_one_Gene=TRUE
                    registro_genes_internos=NA

                    gen2_start=$(grep -w -F $gen2 $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$pair"_Second_PHO_"$species_name"_gene_coordinates.txt" | awk -F "\t" '{print $3}')
                    gen2_end=$(grep -w -F $gen2 $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$pair"_Second_PHO_"$species_name"_gene_coordinates.txt" | awk -F "\t" '{print $4}')
                    gen2_strand=$(grep -w -F $gen2 $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$pair"_Second_PHO_"$species_name"_gene_coordinates.txt" | awk -F "\t" '{print $5}')

                    interval_start=$(echo $gen1_start" "$gen1_end" "$gen2_start" "$gen2_end | tr " " "\n" | sort -n | head -n 1)
                    interval_end=$(echo $gen1_start" "$gen1_end" "$gen2_start" "$gen2_end | tr " " "\n" | sort -n | tail -n 1)

                    interval_genes=$(gffread -R -r $chr":"$interval_start".."$interval_end $work_directory"/"$Result_Storage"/Modified_GTF_Files/"$species_name"_with_Chimera_genes.gtf" | grep -v "#" | awk -F "\t" '{if ($3=="transcript") print $9}' | sed "s/^ID=//" | grep -v -w -F $gen1 | grep -v -w -F $gen2)
                    N_interval_genes=$(echo $interval_genes | tr " " "\n" | grep -c .)

                    if [ "$gen1_strand" != "$gen2_strand" ]
                    then
                      Internal_Status=Strand_Change
                    elif [ $N_interval_genes -eq 0 ]
                    then
                      Internal_Status=Viable_Operon
                    else
                      Internal_Status=Interruping_Genes
                    fi

                    if [ $N_interval_genes -gt 0 ] && [ $N_interval_genes -le 5 ]
                    then
                      for int_gen in $interval_genes
                      do
                        int_pho=$(grep -F -w $int_gen $PHO_File | awk -F "\t" '{print $1}')
                        registro_genes_internos=$(echo $registro_genes_internos";"$int_gen"("$int_pho")")
                      done
                      registro_genes_internos=$(echo $registro_genes_internos | sed "s/NA;//")
                    elif [ $N_interval_genes -gt 5 ]
                    then
                      registro_genes_internos=Too_Many
                    fi
                    echo $pair" "$chr" "$gen1"("$gen1_strand") "$gen2"("$gen2_strand") "$Internal_Status" "$N_interval_genes" "$registro_genes_internos | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$result_storage"_details/"$pair"_"$species_name"_"$result_storage"_Gene_Details.txt"
                  else
                    echo $pair" "$chr" "$gen1"("$gen1_strand") "$gen2"("$gen2_strand") Same_Gene! NA NA" | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$result_storage"_details/"$pair"_"$species_name"_"$result_storage"_Gene_Details.txt"
                  fi
                done
              done
            done

            if [ "$More_than_one_Gene" == "FALSE" ]
            then
              Status=Single
            else
              Status=$(grep -w -c "Viable_Operon" $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$result_storage"_details/"$pair"_"$species_name"_"$result_storage"_Gene_Details.txt" | awk -F "\t" '{if ($1>0) print "Viable"; else print "Compromised_Synteny"}')
            fi

          fi
          rm $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/temp_first_chr"
        fi

        echo $pair" "$species_name" "$Status" "$Missing_PHO" "$first_alternative_pairs" "$second_alternative_pairs | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/Summary_"$result_storage".txt"
      fi
      count_species=$(($count_species + 1))
    done
  done
}

basic_info_ace_counts (){
  ace=$(echo $sl_data | awk -F ";" '{print $1}')
  chr=$(echo $sl_data | awk -F ";" '{print $2}')
  strand=$(echo $sl_data | awk -F ";" '{print $3}')
  check_strand=$(echo $strand | awk '{if ($1=="+") print 1; else print 2}')

  if [ $spe == "Paragonimus_heterotremus" ] || [ $spe == "Fasciolopsis_buski" ]
  then
    gen_id=$(grep -w Concordante $work_directory"/SL_Read_Counts/"$spe"_tier3.counts" | grep -w -v Trans_End | awk -F "\t" -v chr=$chr -v ace=$ace '{if ($2==ace && $3==chr) print $1}' | sed "s/^gene-//")
  else
    gen_id=$(grep -w Concordante $work_directory"/SL_Read_Counts/"$spe"_tier3.counts" | grep -w -v Trans_End | awk -F "\t" -v chr=$chr -v ace=$ace '{if ($2==ace && $3==chr) print $1}')
  fi

  check_double_gene_nonsense=$(echo $gen_id | tr " " "\n" | grep -c . | awk '{if ($1>1) print "Multiple"; else print "Single"}')
  trans_ids=

  is_operon=$(echo .)

  for gen in $gen_id
  do
    trans=$(grep -w $gen'$' $work_directory"/"$Result_Storage"/Mediciones_Expresion/Background_TPM/"$spe"_per_gene.txt" | grep -F $gen | awk -F "\t" '{print $1}' | tr "\n" ";" | sed "s/;$/\n/")
    trans_ids=$(echo $trans_ids";"$trans)

    check_operon=$(grep -w -F $gen $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | grep -c Operon )
    if [ $check_operon -gt 0 ]
    then
      is_operon=TRUE
    fi

  done
  trans_ids=$(echo $trans_ids | sed "s/^;//")
}

get_pho_info_phos () {
  mkdir $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/"$group_name
  N_species=$(echo $wanted_species | tr " " "\n" | grep -c . | awk '{print $1-int($1/3)}')
  for spe in $wanted_species
  do
    echo "get info: "$group_name": "$spe
    awk -F "\t" '{if ($6=="SL") print $11}' $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | sed "s/__/\n/g" | sort -u | grep -v No_Hits >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/"$group_name"/pho_sl.temp"
    grep Operon $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | grep -v Read_Coverage |  awk -F "\t" '{print $11}' | sed "s/__/\n/g" | sort -u | grep -v No_Hits >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/"$group_name"/pho_operon.temp"
  done

  sort $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/"$group_name"/pho_sl.temp" | uniq -c | awk '{print $1";"$2}' | tr -d " " | awk -F ";" -v cutoff=$N_species '{if ($1>=cutoff) print $2}' >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/"$group_name"/PHO_SL.txt"
  sort $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/"$group_name"/pho_operon.temp" | uniq -c | awk '{print $1";"$2}' | tr -d " " | awk -F ";" -v cutoff=$N_species '{if ($1>=cutoff) print $2}' >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/"$group_name"/PHO_Operon.txt"
  rm $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/"$group_name"/"*".temp"
}

comparacion_phos_entre_grupos (){
  header_sl=Grupo
  header_operon=Grupo
  for int in $grupos_interes2
  do
    N_base_pho_SL=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/"$int"/PHO_SL.txt")
    N_base_pho_Operon=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/"$int"/PHO_Operon.txt")

    header_sl=$(echo $header_sl" "$int"("$N_base_pho_SL")")
    header_operon=$(echo $header_operon" "$int"("$N_base_pho_Operon")")
  done

  echo $header_sl | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/Heatmap_shared_SL_PHOs_"$storage"_raw.txt"
  echo $header_operon | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/Heatmap_shared_Operon_PHOs_"$storage"_raw.txt"

  echo $header_sl | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/Heatmap_shared_SL_PHOs_"$storage"_percentage.txt"
  echo $header_operon | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/Heatmap_shared_Operon_PHOs_"$storage"_percentage.txt"

  for interes1 in $grupos_interes1
  do
    N_base_pho_SL=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/"$interes1"/PHO_SL.txt")
    N_base_pho_Operon=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/"$interes1"/PHO_Operon.txt")

    registro_pho_SL_raw=$(echo $interes1"("$N_base_pho_SL")")
    registro_pho_operon_raw=$(echo $interes1"("$N_base_pho_Operon")")

    registro_pho_SL_percentage=$(echo $interes1"("$N_base_pho_SL")")
    registro_pho_operon_percentage=$(echo $interes1"("$N_base_pho_Operon")")

    for interes2 in $grupos_interes2
    do
      echo $interes1" vs "$interes2

      if [ "$interes1" == "$interes2" ]
      then
        registro_pho_SL_raw=$(echo $registro_pho_SL_raw" NA")
        registro_pho_operon_raw=$(echo $registro_pho_operon_raw" NA")

        registro_pho_SL_percentage=$(echo $registro_pho_SL_percentage" NA")
        registro_pho_operon_percentage=$(echo $registro_pho_operon_percentage" NA")
      else
        N_pho_SL=$(grep -w -c -f $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/"$interes1"/PHO_SL.txt" $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/"$interes2"/PHO_SL.txt")
        N_pho_operon=$(grep -w -c -f 2 $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/"$interes2"/PHO_Operon.txt")

        echo "Total:"$N_pho_SL >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/Pairwise_marches/"$interes1"_vs_"$interes2"_PHO_SL_"$storage".txt"
        grep -w -f $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/"$interes1"/PHO_SL.txt" $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/"$interes2"/PHO_SL.txt" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/Pairwise_marches/"$interes1"_vs_"$interes2"_PHO_SL_"$storage".txt"

        echo "Total:"$N_pho_operon >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/Pairwise_marches/"$interes1"_vs_"$interes2"_PHO_Operon_"$storage".txt"
        grep -w -f $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/"$interes1"/PHO_Operon.txt" $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/"$interes2"/PHO_Operon.txt" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/Pairwise_marches/"$interes1"_vs_"$interes2"_PHO_Operon_"$storage".txt"

        registro_pho_SL_raw=$(echo $registro_pho_SL_raw" "$N_pho_SL)
        registro_pho_operon_raw=$(echo $registro_pho_operon_raw" "$N_pho_operon)

        base_sl=$(cat $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/"$interes1"/PHO_SL.txt" $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/"$interes2"/PHO_SL.txt" | sort -u | grep -c .)
        base_operon=$(cat $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/"$interes1"/PHO_Operon.txt" $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/"$interes2"/PHO_Operon.txt"  | sort -u | grep -c .)

        percentage_sl=$(echo $N_pho_SL" "$base_sl | awk -F " " '{printf "%.10f\n", $1/($1+$2)}')
        percentage_operon=$(echo $N_pho_operon" "$base_sl | awk -F " " '{printf "%.10f\n", $1/($1+$2)}')

        registro_pho_SL_percentage=$(echo $registro_pho_SL_percentage" "$percentage_sl)
        registro_pho_operon_percentage=$(echo $registro_pho_operon_percentage" "$percentage_operon)
      fi
    done

    echo $registro_pho_SL_raw | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/Heatmap_shared_SL_PHOs_"$storage"_raw.txt"
    echo $registro_pho_operon_raw | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/Heatmap_shared_Operon_PHOs_"$storage"_raw.txt"

    echo $registro_pho_SL_percentage | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/Heatmap_shared_SL_PHOs_"$storage"_percentage.txt"
    echo $registro_pho_operon_percentage | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/Heatmap_shared_Operon_PHOs_"$storage"_percentage.txt"
  done
}

sj_ace_information () {
  count=1
  recorrer=$(grep -c . $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/Resgistros_por_sitio/"$spe"_"$ace_set"_ace.txt")

  N_sitios_SJ_over_4=0
  while [ $count -le $recorrer ]
  do
    ace=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/Resgistros_por_sitio/"$spe"_"$ace_set"_ace.txt" | awk -F "__" '{print $1}')
    chr=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/Resgistros_por_sitio/"$spe"_"$ace_set"_ace.txt" | awk -F "__" '{print $2}')
    sj_info_ace_info=$(awk -F "\t" -v ace=$ace -v chr=$chr '{if ($1==ace && $2==chr) print}' $storage"/"$spe"_Cis_SJ_on_SL_ace.tab" | awk -F "\t" '{print $9" "$11}')
    control_single=$(awk -F "\t" -v ace=$ace -v chr=$chr '{if ($1==ace && $2==chr) print $4}' $storage"/"$spe"_Cis_SJ_on_SL_ace.tab")

    if [ $control_single == "Single" ]
    then
      check_sj_reads=$(echo $sj_info_ace_info | awk -F " " '{print $1}')
      if [ $check_sj_reads -ge 4 ]
      then
        N_sitios_SJ_over_4=$(($N_sitios_SJ_over_4 + 1))
      fi

      registro_reads_sl=X
      todos_reads_SL=0

      for tag in $sl_tags
      do
        confirm_gen=$(grep -w Concordante $work_directory"/SL_Read_Counts/"$spe"_tier3.counts" | grep -w -v Trans_End | awk -F "\t" -v ace=$ace -v chr=$chr '{if ($2==ace && $3==chr) print $1}')
        echo "Debug: "$confirm_gen" -- "$ace_set

        check_tag_counts=$(grep -w Concordante $work_directory"/SL_Read_Counts/"$spe"_"$tag"_read.counts" | awk -F "\t" -v ace=$ace -v chr=$chr -v gen=$confirm_gen '{if ($1==gen && $5==ace && $3==chr) print}' | grep -w -v Trans_End | grep -c .)
        if [ $check_tag_counts -gt 0 ]
        then
          tag_counts=$(grep -w Concordante $work_directory"/SL_Read_Counts/"$spe"_"$tag"_read.counts" | awk -F "\t" -v ace=$ace -v chr=$chr -v gen=$confirm_gen '{if ($1==gen && $5==ace && $3==chr) print}' | grep -w -v Trans_End | tr "\t" "\n" | tail -n 1)
        else
          tag_counts=0
        fi
        todos_reads_SL=$(($todos_reads_SL + $tag_counts))
        registro_reads_sl=$(echo $registro_reads_sl" "$tag_counts)
      done
      registro_reads_sl=$(echo $registro_reads_sl | sed "s/X //")

      # echo $ace_set" "$ace" "$chr" "$sj_info_ace_info" || "$registro_reads_sl
      echo $ace_set" "$ace" "$chr" "$sj_info_ace_info" "$registro_reads_sl" "$todos_reads_SL | tr " " "\t" >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt"
    fi
    count=$(($count + 1))
  done
  echo $spe" "$ace_set" "$N_sitios_SJ_over_4" "$recorrer >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/Resumen_Buenos_sitios_competencia.txt"
}

final_site_table () {
  for spe in $Especies
  do
    gff_file=$(ls $grupo"/"$spe"/"*"annotations.gtf")
    SL_TAGS=$(ls $work_directory"/SL_Read_Counts/"$spe"_"*"_read.counts" | sed "s/.*\///" | sed "s/$spe//" | awk -F "_" '{print $2"_"$3}' )

    header_spe_tags=$(echo $SL_TAGS | tr " " ";")
    col_all_count=$(head -n 1 $work_directory"/SL_Read_Counts/"$spe"_total.counts" | tr "\t" "\n" | grep -n Total | awk -F ":" '{print $1}')
    sj_info=$work_directory"/"$Result_Storage"/Mapeos_Second_Pass/"$spe

    # Tengo que eliminar un puñado de pseudogenes que habñian roto las tablas pero nunca lo corregí
    site_interes=$(grep -v FGIG_13568 $work_directory"/SL_Read_Counts/"$spe"_total.counts" | grep -v Pseudo_ | awk -F "\t" '{print $2"__"$3}' | grep -v Aceptor__Chromosome | sort -u)

    total_sites=$(echo $site_interes | tr " " "\n" | grep -c .)
    current_site=1

    echo "Grupo Species Acceptor Chromosome Ace_Strand All_SL_Trans-Splicing_Reads Cis-Splcing_Reads Cis-Splicing_Donnor_Sites GenID Gene_Numer Gen_SL_Trans-Splicing_Reads "$header_spe_tags" Discarded_From_Operon Chimeric_Portion Matching_Strand Relative_Location Operon_Position Role_Site Role_Site_Analysis" | tr " " "\t" >> $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/"$spe"_SL_Site_Full_Summary.txt"

    for site in $site_interes
    do
      echo "Corriendo "$spe":"$current_site"/"$total_sites
      current_site=$(($current_site + 1))

      ace=$(echo $site | awk -F "__" '{print $1}')
      chr=$(echo $site | awk -F "__" '{print $2}')

      total_SL_reads_all=$(grep -w $ace $work_directory"/SL_Read_Counts/"$spe"_total.counts" | grep -F -w $chr  | awk -F "\t" -v N=$col_all_count '{print $N}' | awk -F "\t" '{ sum += $1 } END {print sum}')

      # 1) Grupo General Sitio
      check_grupo_general=$(grep -w $site $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/Resgistros_por_sitio/"$spe"_"*"_ace.txt" | grep -c .)
      if [ $check_grupo_general -gt 0 ]
      then
        grupo_general=$(grep -w $site $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/Resgistros_por_sitio/"$spe"_"*"_ace.txt" | sed "s/.*\///" | awk -F ":" '{print $1}' | sed "s/_ace.txt//" |  sed "s/$spe//" | sed "s/^_//" | tr -d "\n" | sed "s/$/\n/" )
      else
        grupo_general=Not_Included
      fi

      # 2) Información SL Tags y genes aceptores
      acceptor_genes=$(grep -w $ace $work_directory"/SL_Read_Counts/"$spe"_total.counts" | grep -F -w $chr | awk -F "\t" '{print $1}' | sort -u | sed "s/^gene-//")

      N_Genes=$(echo $acceptor_genes | tr " " "\n"  | grep -c .)
      gen_number=1
      for gen in $acceptor_genes
      do
        check_concordancia=$(grep -w $ace $work_directory"/SL_Read_Counts/"$spe"_total.counts" | grep -P $gen"\t" | grep -F $gen | awk -F "\t" '{print $5}' | grep -c "Concordante" | awk '{if ($1==0) print "."; else print "X"}')
        gen_strand=$(grep -w $ace $work_directory"/SL_Read_Counts/"$spe"_total.counts" | grep -P $gen"\t" | grep -F $gen  | awk -F "\t" '{print $4}')

        # 2.1) general
        gene_start=$(grep $gen'";' $gff_file | awk -F "\t" '{if ($3=="CDS") print $4}' | sort -u -n | head -n 1 )
        gene_end=$(grep $gen'";' $gff_file | awk -F "\t" '{if ($3=="CDS") print $5}' | sort -u -n | tail -n 1 )

        if [ $check_concordancia == "X" ]
        then
          ace_strand=$(echo $gen_strand | sed "s/+/Plus/" | sed "s/-/Minus/" )
        else
          if [ $gen_strand == "+" ]
          then
            ace_strand=Minus
          else
            ace_strand=Plus
          fi
        fi

        # 2.2 Posicion respecto a la maxima extencioon del gen (CDS)
        if [ "$gen_strand" == "+" ]
        then
          if [ $ace -le $gene_start ]
          then
            location_vs_gene_model=Start
          elif [ $ace -ge $gene_end ]
          then
            location_vs_gene_model=End
          else
            location_vs_gene_model=Internal
          fi
        else
          if [ $ace -ge $gene_end ]
          then
            location_vs_gene_model=Start
          elif [ $ace -le $gene_start ]
          then
            location_vs_gene_model=End
          else
            location_vs_gene_model=Internal
          fi
        fi

        # 3.3 Conteos SL
        total_SL_reads_all_gene=$(grep -w $ace $work_directory"/SL_Read_Counts/"$spe"_total.counts" | grep -P $gen"\t" | grep -F $gen | awk -F "\t" -v N=$col_all_count '{print $N}')
        registro_SL_Tags=X
        for tag in $SL_TAGS
        do
          check_data=$(grep -w $ace $work_directory"/SL_Read_Counts/"$spe"_"$tag"_read.counts" | grep -P $gen"\t" | grep -F $gen | grep -c .)
          if [ $check_data -gt 0 ]
          then
            column=$(head -n 1 $work_directory"/SL_Read_Counts/"$spe"_"$tag"_read.counts"  | tr "\t" "\n" | grep -n Total | awk -F ":" '{print $1}')
            count_tag=$(grep -w $ace $work_directory"/SL_Read_Counts/"$spe"_"$tag"_read.counts" | grep -P $gen"\t" | grep -F $gen | awk -F "\t" -v N=$column '{print $N}' | sort -u )
          else
            count_tag=NA
          fi
          registro_SL_Tags=$(echo $registro_SL_Tags";"$count_tag)
        done
        registro_SL_Tags=$(echo $registro_SL_Tags | sed "s/X;//")

        # 3.4 El gen forma parte de un operon y/o una quimera. Y cual porcion
        is_properly_scanned=$(grep -w $spe $work_directory"/"$Result_Storage"/Operons/Skipped_Genes.txt" | grep -w -F -c $gen)
        if [ $is_properly_scanned -gt 0 ]
        then
          is_chimera=NA
          is_operon=NA
          operon_position=NA
          chimera_part=NA
          Operon_Resolution=NA
          Skipped=$(grep -w $spe $work_directory"/"$Result_Storage"/Operons/Skipped_Genes.txt" | grep -w -F $gen | awk -F " " '{print $3}')
        else
          Skipped=$(echo ".")
          is_chimera=$(grep -w -F $gen $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{if ($7=="Chimeric") print "X"; else print "."}')
          is_operon=$(grep -w -F $gen $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $10}')

          check_non_operon=$(echo $is_operon | grep -w -c "NA")
          if [ $check_non_operon -gt 0 ]
          then
            operon_position=NA
          else
            genes_in_operon=$(grep -w $is_operon $work_directory"/"$Result_Storage"/Operons/"$spe"_operon.txt"  | awk -F "\t" '{print $6}' | sed "s/_(/\n/g" | sed "s/.*)_//" | tr -d "*")
            registro_genes_operon=NA

            for genOP in $genes_in_operon
            do
              check_chimera=$(echo $genOP | grep -c -F "[C]")
              if [ $check_chimera -gt 0 ]
              then
                grep -w $ace $work_directory"/SL_Read_Counts/"$spe"_"*"_read.counts" | grep -P $gen"\t" | grep -F $gen | awk -F "\t" '{print $2"_Chimera"}' | sed "s/^gene-//" > $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/Temporales/chimera_search.tmp"
                chimera_sections=$(grep -F -f $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/Temporales/chimera_search.tmp" $work_directory"/"$Result_Storage"/"$spe".pep" | sed "s/>//")
                registro_genes_operon=$(echo $registro_genes_operon" "$chimera_sections)
              else
                registro_genes_operon=$(echo $registro_genes_operon" "$genOP)
              fi
              registro_genes_operon=$(echo $registro_genes_operon | sed "s/NA //")
            done
          fi

          if [ $is_chimera == "." ]
          then
            chimera_part=NA
            if [ $check_non_operon -eq 0 ]
            then
              operon_position=$(echo $registro_genes_operon | tr " " "\n" | grep -F -w -n $gen | awk -F ":" '{print $1}')
            else
              operon_position=NA
            fi
          else
            grep -w $ace $work_directory"/SL_Read_Counts/"$spe"_"*"_read.counts" | grep -P $gen"\t" | grep -F $gen | awk -F "\t" '{print $2"_Chimera"}' | sed "s/^gene-//" > $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/Temporales/chimera_search.tmp"
            chimera_sections=$(grep -F -f $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/Temporales/chimera_search.tmp" $work_directory"/"$Result_Storage"/"$spe".pep" | sed "s/>//")
            recorrer_chimeras=$(echo $chimera_sections | tr " " "\n" | grep -c .)
            count_chimera=1

            chimera_part=Wololo

            while [ $count_chimera -le $recorrer_chimeras ]
            do
              chimeraID=$(echo $chimera_sections | tr " " "\n" | sed -n $count_chimera"p")
              check_primero=$(echo $chimeraID | grep -c _Start_)
              check_ultimo=$(echo $chimeraID | grep -c _End)

              if [ $check_primero -gt 0 ]
              then
                threshold=$(echo $chimeraID | sed "s/.*_Start_//")
                chimera_part=$(echo $ace" "$threshold | awk -F " " -v strand=$gen_strand '{if (strand=="+" && $1 < $2) print "Hit"; else if (strand=="-" && $1 > $2) print "Hit"; else print "NA"}' )
              elif [ $check_ultimo -gt 0 ]
              then
                threshold=$(echo $chimeraID | sed "s/.*__//" | sed "s/_End//")
                chimera_part=$(echo $ace" "$threshold | awk -F " " -v strand=$gen_strand '{if (strand=="+" && $1 >= $2) print "Hit"; else if (strand=="-" && $1 <= $2) print "Hit"; else print "NA"}' )
              else
                threshold_1=$(echo $chimeraID | sed "s/.*__//" | awk -F "_" '{print $1}')
                threshold_2=$(echo $chimeraID | sed "s/.*__//" | awk -F "_" '{print $2}')
                chimera_part=$(echo $ace" "$threshold_1" "$threshold_2 | awk -F " " -v strand=$gen_strand '{if (strand=="+" && $1 >= $2 && $1 < $3) print "Hit plus"; if (strand=="-" && $1 > $2 && $1 <= $3) print "Hit minus"; else if ($1 < $2 || $1 > $3) print "NA"}')
              fi

              if [ "$chimera_part" == "NA" ]
              then
                count_chimera=$(($count_chimera + 1))
              else
                chimera_part=$chimeraID
                count_chimera=9000
                operon_position=$(echo $registro_genes_operon | tr " " "\n" | grep -F -w -n $chimera_part | awk -F ":" '{print $1}')
              fi
            done
          fi
        fi

        # 3.5 Importancia en la resolución de operones
        if [ $Skipped != "." ]
        then
          Operon_Resolution=NA
        elif [ $check_concordancia != "X" ]
        then
          Operon_Resolution=Missmatching_Strand
        elif  [ $location_vs_gene_model == "End" ]
        then
          Operon_Resolution=Out_of_Gene_Model
        elif [ $total_SL_reads_all -le 3 ]
        then
          Operon_Resolution=Low_SL_Read_Counts
        else
          if [ "$operon_position" == "NA" ] || [ "$operon_position" == "1" ]
          then
            Operon_Resolution=No
          else
            Operon_Resolution=Yes
          fi
        fi

        # 4) informacion cis-splicing (repetitivo en este loop pero preciso la información de la hebra de los sitios)
        check_cis_splicing_info=$(grep -w $ace $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/"$spe"_Cis_SJ_on_SL_ace.tab" | grep -F -w -c $chr)
        if [ $check_cis_splicing_info -gt 0 ]
        then
          cis_splicing_reads=$(grep -w $ace $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/"$spe"_Cis_SJ_on_SL_ace.tab" | grep -F -w $chr  | awk -F "\t" '{print $9}')
          cis_splicing_donnor_site=$(grep -w $ace $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/"$spe"_Cis_SJ_on_SL_ace.tab" | grep -F -w $chr  | awk -F "\t" '{print $11}')
        else
          check_strand=$(echo $ace_strand | awk '{if ($1=="Plus") print 1; else print 2}')
          cis_splicing_donnor_site=$(awk -F "\t" -v chr=$chr -v strand=$check_strand '{if ($1==chr && $4==strand && $7>0) print}' $sj_info"_"*"_SJ.out.tab" | awk -F "\t" -v ace=$ace '{if ($2==ace || $3==ace) print $1,$2,$3}' | sort -u | grep -c .)
          if [ $cis_splicing_donnor_site -gt 0 ]
          then
            cis_splicing_reads=$(awk -F "\t" -v chr=$chr -v strand=$check_strand '{if ($1==chr && $4==strand && $7>0) print}' $sj_info"_"*"_SJ.out.tab" | awk -F "\t" -v ace=$ace '{if ($2==ace || $3==ace) print $7}' | awk -F "\t" '{ sum += $1 } END {print sum}')
          else
            cis_splicing_reads=0
          fi
        fi

        echo ""
        echo "grupo: "$grupo
        echo "spe "$spe
        echo "ace: "$ace
        echo "chr: "$chr
        echo "ace_strand: "$ace_strand
        echo "total_SL_reads_all: "$total_SL_reads_all
        echo "cis_splicing_reads: "$cis_splicing_reads
        echo "cis_splicing_donnor_site: "$cis_splicing_donnor_site
        echo "gen: "$gen
        echo "gen_number: "$gen_number"/"$N_Genes
        echo "Skipped: "$Skipped
        echo "chimera_part: "$chimera_part
        echo "check_concordancia: "$check_concordancia
        echo "location_vs_gene_model: "$location_vs_gene_model
        echo "operon_position: "$operon_position
        echo "Operon_Resolution: "$Operon_Resolution"("$grupo_general")"
        echo "total_SL_reads_all_gene: "$total_SL_reads_all_gene
        echo "registro_SL_Tags: "$registro_SL_Tags
        echo "############################################################"
        echo ""

        echo $grupo" "$spe" "$ace" "$chr" "$ace_strand" "$total_SL_reads_all" "$cis_splicing_reads" "$cis_splicing_donnor_site" "$gen" "$gen_number"_of_"$N_Genes" "$total_SL_reads_all_gene" "$registro_SL_Tags" "$Skipped" "$chimera_part" "$check_concordancia" "$location_vs_gene_model" "$operon_position" "$Operon_Resolution" "$grupo_general | tr " " "\t" >> $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/"$spe"_SL_Site_Full_Summary.txt"
        gen_number=$(($gen_number + 1))
      done
    done
  done
}

upset_input () {
  all_items=$(sed 1d $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_exploration.txt" | awk -F "\t" '{print $1}')
  rm -r $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput"
  mkdir $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput"

  echo "Item "$Especies | tr " " "\t" > $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_utest_matrix.txt"
  for item in $all_items
  do
    registro=$item
    for spe in $Especies
    do
      echo $spe" "$item

      check_Grupos_Cestoda_Derivados=$(echo $Grupos_Cestoda_Derivados | grep -c $spe)
      check_Grupos_Cestoda_Basal=$(echo $Grupos_Cestoda_Basal | grep -c $spe)
      check_Grupos_Schistosmidae=$(echo $Grupos_Schistosmidae | grep -c $spe)
      check_Grupos_Trematoda_No_schistosoma=$(echo $Grupos_Trematoda_No_schistosoma | grep -c $spe)

      check_data=$(grep -w -c $item $data_path"/"$spe"_"$data_id)
      registro=$(echo $registro" "$check_data)

      if [ $check_data -gt 0 ]
      then
        if [ $check_Grupos_Cestoda_Derivados -gt 0 ]
        then
          echo $item >> $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/Cyclophyllidea_graph.tmp"
        elif [ $check_Grupos_Cestoda_Basal -gt 0 ]
        then
          echo $item >> $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/Basal_Cestoda_graph.tmp"
        elif [ $check_Grupos_Schistosmidae -gt 0 ]
        then
          echo $item >> $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/Schistosomatids_graph.tmp"
        elif [ $check_Grupos_Trematoda_No_schistosoma -gt 0 ]
        then
          echo $item >> $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/Other_Trematoda_graph.tmp"
        fi
      fi
    done
    echo $registro  | tr " " "\t" >> $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_utest_matrix.txt"
  done

  sort -u $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/Cyclophyllidea_graph.tmp" >> $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/work_data.tmp"
  sort -u $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/Basal_Cestoda_graph.tmp" >> $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/work_data.tmp"
  sort -u $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/Schistosomatids_graph.tmp" >> $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/work_data.tmp"
  sort -u $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/Other_Trematoda_graph.tmp" >> $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/work_data.tmp"

  sort $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/work_data.tmp" | uniq -d >> $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/work_data.txt"

  grep -f $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/work_data.txt" $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/Cyclophyllidea_graph.tmp" | sort -u >> $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/Cyclophyllidea_graph.in"
  grep -f $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/work_data.txt" $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/Basal_Cestoda_graph.tmp" | sort -u >> $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/Basal_Cestoda_graph.in"
  grep -f $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/work_data.txt" $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/Schistosomatids_graph.tmp" | sort -u >> $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/Schistosomatids_graph.in"
  grep -f $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/work_data.txt" $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/Other_Trematoda_graph.tmp" | sort -u >> $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/Other_Trematoda_graph.in"

  ######################################################################################################################################################################################################################################################################################################

  grep -v -f $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/work_data.txt" $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/Cyclophyllidea_graph.tmp" | sort -u >> $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/Cyclophyllidea_uniq.id"
  grep -v -f $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/work_data.txt" $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/Basal_Cestoda_graph.tmp" | sort -u >> $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/Basal_Cestoda_uniq.id"
  grep -v -f $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/work_data.txt" $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/Schistosomatids_graph.tmp" | sort -u >> $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/Schistosomatids_uniq.id"
  grep -v -f $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/work_data.txt" $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/Other_Trematoda_graph.tmp" | sort -u >> $work_directory"/"$Result_Storage"/Upset_graph/"$analysis"_Intput/Other_Trematoda_uniq.id"
}

################################################################################################################################################
################################################################################################################################################

if [ "$PREPARE_TROUBLE_MAKERS" == "TRUE" ]
then
  rm -r $work_directory"/"$Result_Storage
  mkdir $work_directory"/"$Result_Storage
  mkdir $work_directory"/"$Result_Storage"/Storage_sequences"

  Filos=$(echo "Cestodes Trematoda")
  for filo in $Filos
  do
    species=$(grep -v N_Assemblies $work_directory"/"$filo"_Summary_hits.tab" | awk -F "\t" '{print $1}')
    for spe in $species
    do
      check_bastard=$(echo $NCBI_GFF_format | grep -c $spe )

      echo $spe": "$check_bastard
      rna_file=$(ls $filo"/"$spe"/"*".CDS_transcripts.fa")

      if [ $check_bastard -eq 0 ]
      then
        cp $rna_file $work_directory"/"$Result_Storage"/Storage_sequences/"$spe"_mRNA.fasta"
        seqkit grep -j $threads -v -s -p "X" "Ref_Proteins/"$spe".protein.fa" >> $work_directory"/"$Result_Storage"/Storage_sequences/"$spe".aa"
      else
        sed "s/.*locus_tag//" $rna_file | sed "s/].*//" | sed "s/=/>/" > $work_directory"/"$Result_Storage"/Storage_sequences/"$spe"_mRNA.fasta"
        rename_proteins=$(seqkit seq -n -i "Ref_Proteins/"$spe".protein.fa")
        for id in $rename_proteins
        do
          rename=$(grep -w -F $id "Ref_Proteins/ID_Cross_over/"$spe"_Anot" |  tr ";" "\n" | grep locus_tag= | sed "s/locus_tag=/>/" | sort -u )

          echo "renaming: "$id" to "$rename
          seqkit grep -j $threads -p $id "Ref_Proteins/"$spe".protein.fa" | sed "s/>.*/$rename/" | seqkit grep -j $threads -v -s -p "X" >> $work_directory"/"$Result_Storage"/Storage_sequences/"$spe".aa"
        done
      fi
    done
  done
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$PREPARE_ORTHOFINDER" == "TRUE" ]
then
  rm $work_directory"/"$Result_Storage"/"*"pep"
  rm -r $work_directory"/"$Result_Storage"/Chimera_Details"
  mkdir $work_directory"/"$Result_Storage"/Chimera_Details"
  mkdir $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details"

  Filos=$(echo "Cestodes Trematoda")

  for filo in $Filos
  do
    species=$(grep -v N_Assemblies $work_directory"/"$filo"_Summary_hits.tab" | awk -F "\t" '{print $1}')
    set_references

    for spe in $species
    do
      rm $work_directory"/"$Result_Storage"/"$spe".pep"
      mkdir $work_directory"/"$Result_Storage"/Chimera_Details/"$spe
      mkdir $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/Genefiles"

      check_bastard=$(echo $NCBI_GFF_format | grep -c $spe )

      prot_file=$work_directory"/"$Result_Storage"/Storage_sequences/"$spe".aa"
      rna_file=$work_directory"/"$Result_Storage"/Storage_sequences/"$spe"_mRNA.fasta"

      sequences=$(seqkit seq -n -i $prot_file | sort -u)

      echo $filo" "$spe" || "$check_bastard
      echo $rna_file
      echo $prot_file
      echo "------------------------------------------"

      for seq in $sequences
      do
        check_artifact=$(grep $seq"_" $work_directory"/Inserciones_Internas/Anotacion/"$spe"_registro_tier3.tab" | grep -c "Artifact")
        if [ $check_artifact -gt 0 ]
        then
          if [ "$spe" ==  "Paragonimus_heterotremus" ] || [ "$spe" ==  "Fasciolopsis_buski" ]
          then
            remove=$(echo "gene-"$seq"_")
            acceptor_sites=$(grep "gene-"$seq"_" $work_directory"/Inserciones_Internas/Anotacion/BLAST_Searches/"$spe"_"*"_registro_tier3_hits.tab" | grep Artifact | awk -F ":" '{print $2}' | awk -F "\t" '{print $1}' | sed "s/$remove//" | awk -F "_" '{print $1}' | sort -u )
          else
            remove=$(echo $seq"_")
            acceptor_sites=$(grep $seq"_" $work_directory"/Inserciones_Internas/Anotacion/BLAST_Searches/"$spe"_"*"_registro_tier3_hits.tab" | grep Artifact | awk -F ":" '{print $2}' | awk -F "\t" '{print $1}' | sed "s/$remove//" | awk -F "_" '{print $1}' | sort -u )
          fi

          seqkit grep -j $threads -p $seq $rna_file >> $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/Genefiles/"$seq".fasta"

          for comp in $comparison_species
          do
            if [ "$comp" != "$spe" ]
            then
              artifact_gene_ids=$(grep $seq"_" $work_directory"/Inserciones_Internas/Anotacion/BLAST_Searches/"$spe"_"$comp"_registro_tier3_hits.tab" | grep "Artifact" | awk -F "\t" '{print $5}' | grep -v -w NA)
              for arti in $artifact_gene_ids
              do
                seqkit grep -j $threads -p $arti "Ref_Proteins/"$comp".protein.fa" >> $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/Genefiles/"$seq"_chimera_hits.temp"
              done
            fi
          done

          seqkit rmdup --quiet -s $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/Genefiles/"$seq"_chimera_hits.temp" >> $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/Genefiles/"$seq"_chimera_hits.fasta"
          rm $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/Genefiles/"$seq"_chimera_hits.temp"

          makeblastdb -dbtype nucl -in  $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/Genefiles/"$seq".fasta" -out $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/Genefiles/"$seq"_ref"

          echo "EMERGENCY CONTROL:"$spe" || "$seq"_chimera_hits.fasta"
          tblastn -query  $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/Genefiles/"$seq"_chimera_hits.fasta" -db $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/Genefiles/"$seq"_ref" -outfmt "6 std qcovhsp" -out $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/Genefiles/"$seq".tblastn" -max_hsps 1
          echo "DONE tblastn"

          awk -F "\t" -v qcov=$qcov '{if ($9 < $10 && $13 <=qcov) print}' $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/Genefiles/"$seq".tblastn" >> $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/Genefiles/"$seq"_filtered.tblastn"

          rm $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/Genefiles/"$seq"_ref"*

          echo $seq" Start 0" | tr " " "\t" > $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/current_aceptor_site_locations.txt"
          for ace in $acceptor_sites
          do
            interest_header=$(grep $seq"_"$ace"_" $work_directory"/Inserciones_Internas/Anotacion/"$spe"_inserciones_internas_anot.fasta" | grep Mitad-A | sed "s/>//" )
            search_seq=$(seqkit grep -j $threads -p $interest_header $work_directory"/Inserciones_Internas/Anotacion/"$spe"_inserciones_internas_anot.fasta" | seqkit seq -w 0 -s)
            seqkit locate -j $threads -p $search_seq $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/Genefiles/"$seq".fasta" | grep -v seqID | awk -F "\t" -v ace=$ace '{print $1"\t"ace"\t"$6}' >> $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/current_aceptor_site_locations.txt"
          done
          gene_length=$(seqkit fx2tab -i -n -l $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/Genefiles/"$seq".fasta" | awk -F "\t" '{print $2}')
          echo $seq" End "$gene_length | tr " " "\t" >> $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/current_aceptor_site_locations.txt"

          sort -n -k 3 $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/current_aceptor_site_locations.txt" > $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/"$spe"_aceptor_site_locations.txt"
          cat $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/"$spe"_aceptor_site_locations.txt" >> $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/Aceptor_site_locations.txt"

          count=1
          check_cutting_places=$(grep -c . $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/"$spe"_aceptor_site_locations.txt")

          echo "Setting intervals"
          while [ $count -lt $check_cutting_places ]
          do
            next=$(($count + 1))
            interval_start=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/"$spe"_aceptor_site_locations.txt" | awk -F "\t" '{print $3+1}')
            interval_end=$(sed -n $next"p" $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/"$spe"_aceptor_site_locations.txt" | awk -F "\t" '{print $3}')

            id_start=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/"$spe"_aceptor_site_locations.txt" | awk -F "\t" '{print $2}')
            id_end=$(sed -n $next"p" $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/"$spe"_aceptor_site_locations.txt" | awk -F "\t" '{print $2}')

            if [ "$id_start" == "Start" ]
            then
              control_location=$interval_end
              if [ $interval_start -lt $interval_end ]
              then
                BLAST_Hits=$(awk -F "\t" -v start=$interval_start -v end=$interval_end -v control=$control_location -v qcov=$qcov '{if ($9 < $10 && $9<=end && start <=$10 && control-$9 >= 30 && $13 >=qcov) print $1}' $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/Genefiles/"$seq".tblastn" | sort -u | tr "\n" ";" | sed "s/;$/\n/")
              else
                BLAST_Hits=ABORT
              fi
            else
              control_location=$interval_start
              if [ $interval_start -lt $interval_end ]
              then
                BLAST_Hits=$(awk -F "\t" -v start=$interval_start -v end=$interval_end -v control=$control_location  -v qcov=$qcov '{if ($9 < $10 && $9<=end && start <=$10 && $10-control >= 30 && $13 >=qcov) print $1}' $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/Genefiles/"$seq".tblastn" | sort -u | tr "\n" ";" | sed "s/;$/\n/")
              else
                BLAST_Hits=ABORT
              fi
            fi
            echo $seq" "$interval_start" "$interval_end" "$BLAST_Hits" "$id_start" "$id_end" "$control_location | tr " " "\t" >> $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/current_intervals.temp"
            count=$(($count + 1))
          done
          cat $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/current_intervals.temp" >> $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/All_cut_points.txt"

          echo "Done setting intervals"
          count=1
          check_cutting_places=$(grep -c . $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/current_intervals.temp")
          chimeric_interval=1

          echo "Preparing sequences"
          while [ $count -le $check_cutting_places ]
          do
            next=$(($count + 1))
            start_interval=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/current_intervals.temp" | awk -F "\t" '{print $2}')
            end_interval=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/current_intervals.temp" | awk -F "\t" '{print $3}')
            blast_hits=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/current_intervals.temp" | awk -F "\t" '{print $4}')

            id_start_interval=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/current_intervals.temp" | awk -F "\t" '{print $5}')
            id_end_interval=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/current_intervals.temp" | awk -F "\t" '{print $6}')

            abemus_data=$(echo $blast_hits | grep -c .)
            if [ $abemus_data -gt 0 ] &&  [ "$blast_hits" != "ABORT" ]
            then
              stop_the_madness=FALSE

              while [ "$stop_the_madness" == "FALSE" ]
              do
                if [ $next -gt $check_cutting_places ]
                then
                  get_chimera
                  stop_the_madness=TRUE
                else
                  Next_intervals=$(sed -n $next"p" $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/current_intervals.temp" | awk -F "\t" '{print $4}')

                  echo $blast_hits | tr ";" "\n" > $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/compare_1.temp"
                  echo $Next_intervals | tr ";" "\n" > $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/compare_2.temp"
                  check_shared_hits=$(grep -c -w -F -f $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/compare_1.temp" $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/compare_2.temp")

                  if [ $check_shared_hits -eq 0 ]
                  then
                    get_chimera
                    stop_the_madness=TRUE
                  else
                    end_interval=$(sed -n $next"p" $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/current_intervals.temp" | awk -F "\t" '{print $3}')
                    id_end_interval=$(sed -n $next"p" $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/current_intervals.temp" | awk -F "\t" '{print $6}')

                    next=$(($count + 1))
                    count=$(($count + 1))
                  fi
                fi
              done
            fi
            count=$(($count + 1))
          done

          if [ -f $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/"$spe"_chimeric_intervals.fasta" ]
          then
            check_number_halfs=$(grep -c ">" $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/"$spe"_chimeric_intervals.fasta")
          else
            check_number_halfs=0
          fi

          echo "DEBUG: "$seq": "$check_number_halfs
          if [ $check_number_halfs -gt 1 ]
          then
            echo $seq >> $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/"$spe"_artifact"
            cat $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/"$spe"_chimeric_intervals.aa" >> $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/"$spe"_chimeric_intervals.aa"
            cat $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/"$spe"_chimeric_intervals.fasta" >>   $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/"$spe"_chimeric_intervals.fasta"
          else
            echo $seq >> $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/Trouble_makers.txt"
          fi

          if [ -f $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/"$spe"_chimeric_intervals.aa" ]
          then
            rm $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/"$spe"_chimeric_intervals.aa"
            rm $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/"$spe"_chimeric_intervals.fasta"
          fi
          rm $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/current_intervals.temp"
        fi
      done

      echo "Getting non Artifacts"
      seqkit grep -j $threads -v -f $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/"$spe"_artifact" $prot_file >> $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/"$spe".pep"

      echo "Done!!"

      check_duplicated=$(echo $Multiple_isoforms | grep -c $spe)
      if [ $check_duplicated -eq 0 ]
      then
        seqkit seq -i $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/"$spe".pep" >> $work_directory"/"$Result_Storage"/"$spe".pep"
        cat $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/"$spe"_chimeric_intervals.aa" >> $work_directory"/"$Result_Storage"/"$spe".pep"
      else
        echo "Measuring things"
        seqkit fx2tab -j $threads -n -l $prot_file >> $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/"$spe".lenght"
        echo "Done!!"

        final_extract=$(awk -F "\t" '{print $1}' $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/"$spe".lenght" | awk -F " " '{print $3}' | sort -u | sed "s/gene=//")
        for dup in $final_extract
        do
          grep -P $dup"\t" $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/"$spe".lenght" | awk -F " " '{print $1"_Chimera"}' > $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/last_minute_artifact.gene"
          check_artifact=$(grep -c -F -f $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/last_minute_artifact.gene" $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/"$spe"_chimeric_intervals.aa")

          if [ $check_artifact -eq 0 ]
          then
            long_isoform=$(grep -w -F $dup $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/"$spe".lenght"  | tr " " "\t" | sort -k 4 -n | tail -n 1 | awk -F "\t" '{print $1}')
            seqkit grep -j $threads -p $long_isoform $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/"$spe".pep" | seqkit seq -i >> $work_directory"/"$Result_Storage"/"$spe".pep"
          else
            grep -F -f $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/last_minute_artifact.gene" $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/"$spe"_chimeric_intervals.aa" | sed "s/>//" | sed "s/_Chimera.*//" | sort -u > $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/last_minute_artifact.ids"
            long_isoform=$(grep -w -F -f $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/last_minute_artifact.ids" $work_directory"/"$Result_Storage"/Chimera_Details/Temp_details/"$spe".lenght"  | tr " " "\t" | sort -k 4 -n | tail -n 1 | awk -F "\t" '{print $1}')
            seqkit grep -j $threads -r -p $long_isoform"_Ch" $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/"$spe"_chimeric_intervals.aa" | seqkit seq -i >> $work_directory"/"$Result_Storage"/"$spe".pep"
          fi
          echo "retrieving sequences:"$dup" ("$long_isoform"): "$check_artifact
        done
      fi
    done
  done
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$RUN_ORTHOFINDER" == "TRUE" ]
then
  orthofinder -f $work_directory"/"$Result_Storage -t $threads -a $threads
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$RUN_MAPEOS" == "TRUE" ]
then
  if [ ! -d $work_directory"/"$Result_Storage"/Mediciones_Expresion" ]
  then
    mkdir $work_directory"/"$Result_Storage"/Mediciones_Expresion"
    mkdir $work_directory"/"$Result_Storage"/Mediciones_Expresion/Mapeos"
    mkdir $work_directory"/"$Result_Storage"/Mediciones_Expresion/Read_temps"
  fi

  grupo=Cestodes
  Especies=$(echo "Egranulosus Emultilocularis Hdiminuta Hmicrostoma Mesocestoides_corti Sparganum_proliferum Spirometra_erinaceieuropaei Taenia_asiatica Taenia_multiceps Taenia_saginata Taenia_solium Schistocephalus_solidus")
  mapeo

  grupo=Trematoda
  Especies=$(echo "Clonorchis_sinensis Fasciola_gigantica Fasciola_hepatica Fasciolopsis_buski Opisthorchis_felineus Paragonimus_heterotremus Paragonimus_westermani Schistosoma_bovis Schistosoma_haematobium Schistosoma_japonicum Schistosoma_mansoni Trichobilharzia_regenti")
  mapeo
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$RUN_MAPEOS_SECOND_PASS" == "TRUE" ]
then
  if [ ! -d $work_directory"/"$Result_Storage"/Mapeos_Second_Pass" ]
  then
    mkdir $work_directory"/"$Result_Storage"/Mapeos_Second_Pass"
    mkdir $work_directory"/"$Result_Storage"/Mapeos_Second_Pass/Read_temps"
  fi

  grupo=Cestodes
  Especies=$(echo "Egranulosus Emultilocularis Hdiminuta Hmicrostoma Mesocestoides_corti Sparganum_proliferum Spirometra_erinaceieuropaei Taenia_asiatica Taenia_multiceps Taenia_saginata Taenia_solium Schistocephalus_solidus")
  mapeo_second_pass

  grupo=Trematoda
  Especies=$(echo "Clonorchis_sinensis Fasciola_gigantica Fasciola_hepatica Fasciolopsis_buski Opisthorchis_felineus Paragonimus_heterotremus Paragonimus_westermani Schistosoma_bovis Schistosoma_haematobium Schistosoma_japonicum Schistosoma_mansoni Trichobilharzia_regenti")

  mapeo_second_pass
fi



##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$PREPARE_CHIMERIC_GTF" == "TRUE" ]
then
  rm -r $work_directory"/"$Result_Storage"/Modified_GTF_Files"
  mkdir $work_directory"/"$Result_Storage"/Modified_GTF_Files"
  mkdir $work_directory"/"$Result_Storage"/Modified_GTF_Files/Temporales"

  grupo=Cestodes
  Especies=$(echo "Egranulosus Emultilocularis Hdiminuta Hmicrostoma Mesocestoides_corti Sparganum_proliferum Spirometra_erinaceieuropaei Taenia_asiatica Taenia_multiceps Taenia_saginata Taenia_solium Schistocephalus_solidus")
  new_GTF

  grupo=Trematoda
  Especies=$(echo "Clonorchis_sinensis Fasciola_gigantica Fasciola_hepatica Fasciolopsis_buski Opisthorchis_felineus Paragonimus_heterotremus Paragonimus_westermani Schistosoma_bovis Schistosoma_haematobium Schistosoma_japonicum Schistosoma_mansoni Trichobilharzia_regenti")
  new_GTF
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$RUN_CONTEOS" == "TRUE" ]
then
  storage_count=$work_directory"/"$Result_Storage"/Mediciones_Expresion/Conteos"

  rm -r $storage_count
  mkdir $storage_count

  grupo=Cestodes
  Especies=$(echo "Egranulosus Emultilocularis Hdiminuta Hmicrostoma Mesocestoides_corti Sparganum_proliferum Spirometra_erinaceieuropaei Taenia_asiatica Taenia_multiceps Taenia_saginata Taenia_solium Schistocephalus_solidus")
  calculate_counts

  grupo=Trematoda
  Especies=$(echo "Clonorchis_sinensis Fasciola_gigantica Fasciola_hepatica Fasciolopsis_buski Opisthorchis_felineus Paragonimus_heterotremus Paragonimus_westermani Schistosoma_bovis Schistosoma_haematobium Schistosoma_japonicum Schistosoma_mansoni Trichobilharzia_regenti")
  calculate_counts
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$RUN_TPM" == "TRUE" ]
then
  storage_tpm=$work_directory"/"$Result_Storage"/Mediciones_Expresion/Conteos/TPM"
  rm -r $storage_tpm
  mkdir $storage_tpm
  mkdir $storage_tpm"/RPK"
  mkdir $storage_tpm"/Temporales"

  storage_count=$work_directory"/"$Result_Storage"/Mediciones_Expresion/Conteos"

  grupo=Cestodes
  Especies=$(echo "Egranulosus Emultilocularis Hdiminuta Hmicrostoma Mesocestoides_corti Sparganum_proliferum Spirometra_erinaceieuropaei Taenia_asiatica Taenia_multiceps Taenia_saginata Taenia_solium Schistocephalus_solidus")
  calculate_TPM

  grupo=Trematoda
  Especies=$(echo "Clonorchis_sinensis Fasciola_gigantica Fasciola_hepatica Fasciolopsis_buski Opisthorchis_felineus Paragonimus_heterotremus Paragonimus_westermani Schistosoma_bovis Schistosoma_haematobium Schistosoma_japonicum Schistosoma_mansoni Trichobilharzia_regenti")
  calculate_TPM
fi

if [ $RUN_BASELINE_TPM == "TRUE" ]
then
  storage_tpm=$work_directory"/"$Result_Storage"/Mediciones_Expresion/Conteos/TPM"

  rm -r $work_directory"/"$Result_Storage"/Mediciones_Expresion/Background_TPM"
  mkdir $work_directory"/"$Result_Storage"/Mediciones_Expresion/Background_TPM"
  mkdir $work_directory"/"$Result_Storage"/Mediciones_Expresion/Background_TPM/Temporales"

  annoying_points=$(echo "Egranulosus Emultilocularis Hdiminuta Hmicrostoma Opisthorchis_felineus Schistosoma_haematobium Schistosoma_japonicum Schistosoma_mansoni") # Especies en los que hay informacion de isoformas que molesta para esto

  grupo=Cestodes
  Especies=$(echo "Egranulosus Emultilocularis Hdiminuta Hmicrostoma Mesocestoides_corti Sparganum_proliferum Spirometra_erinaceieuropaei Taenia_asiatica Taenia_multiceps Taenia_saginata Taenia_solium Schistocephalus_solidus")
  sl_background

  grupo=Trematoda
  Especies=$(echo "Clonorchis_sinensis Fasciola_gigantica Fasciola_hepatica Fasciolopsis_buski Opisthorchis_felineus Paragonimus_heterotremus Paragonimus_westermani Schistosoma_bovis Schistosoma_haematobium Schistosoma_japonicum Schistosoma_mansoni Trichobilharzia_regenti")
  sl_background
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$RUN_HOG_INFORMATION" == "TRUE" ]
then
  storage_place=PHO_Selection_TPM

  rm -r $work_directory"/"$Result_Storage"/"$storage_place""
  mkdir $work_directory"/"$Result_Storage"/"$storage_place""
  mkdir $work_directory"/"$Result_Storage"/"$storage_place"/Temporales"

  results_orthofinder=$(ls -d $work_directory"/"$Result_Storage"/OrthoFinder/Results_"*)
  PHO_File=$(echo $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")

  all_phos=$(awk -F "\t" '{print $1}' $PHO_File | grep -v -w "HOG" | sort -u)

  no_points=$(echo "Clonorchis_sinensis Fasciola_hepatica")
  Especies=$(echo "Egranulosus Emultilocularis Hdiminuta Hmicrostoma Mesocestoides_corti Sparganum_proliferum Spirometra_erinaceieuropaei Taenia_asiatica Taenia_multiceps Taenia_saginata Taenia_solium Schistocephalus_solidus Clonorchis_sinensis Fasciola_gigantica Fasciola_hepatica Fasciolopsis_buski Opisthorchis_felineus Paragonimus_heterotremus Paragonimus_westermani Schistosoma_bovis Schistosoma_haematobium Schistosoma_japonicum Schistosoma_mansoni Trichobilharzia_regenti")

  echo "PHO Species Gen Average_TPM SL Above_SL_Median" | tr " " "\t" >> $work_directory"/"$Result_Storage"/"$storage_place"/PHO_Gene_Registry.txt"

  for pho in $all_phos
  do
    for spe in $Especies
    do
      spe_num=$(head -n 1 $PHO_File | tr "\t" "\n" | grep -n $spe | awk -F ":" '{print $1}')
      species_SL_median=$(grep -w $spe $work_directory"/"$Result_Storage"/Mediciones_Expresion/Background_TPM/SL_Reference_median.txt" | awk -F "\t" '{print $2}')
      grupo=$(ls $work_directory"/"$Result_Storage"/Mediciones_Expresion/Mapeos/"*"_"$spe"_"*"_Aligned.sortedByCoord.out.bam" | sed "s/.*\///" | sed "s/$spe.*//" | sed "s/_$//" | sort -u)
      check_points=$(echo $no_points | grep -c $spe)
      transcript_pho=$(grep $pho $PHO_File | awk -v N=$spe_num -F "\t" '{print $N}' | tr -d " " | tr "," "\n")

      echo $pho" "$grupo" "$spe
      for transcript in $transcript_pho
      do
        check_SL=$(grep -w -F -c $transcript $work_directory"/"$Result_Storage"/Mediciones_Expresion/Background_TPM/"$spe"_per_gene.txt" | awk -F "\t" '{if ($1>=1) print "SL"; else print "X"}')
        results_gene=$(grep -w -F $transcript $work_directory"/"$Result_Storage"/Mediciones_Expresion/Conteos/TPM/"$grupo"_"$spe"_TPM.counts" | tr "\t" "\n" | sed 1d | awk '{ sum += $1 } END { if (NR > 0) print sum / NR}' | awk '{printf "%.2f\n", $1}' | awk -v spe_sl=$species_SL_median -v SL=$check_SL '{if ($1>=spe_sl) print $1" "SL" Yes"; else print $1" "SL" No"}' )

        echo $pho" "$spe" "$transcript" "$results_gene
        echo $pho" "$spe" "$transcript" "$results_gene | tr " " "\t" >> $work_directory"/"$Result_Storage"/"$storage_place"/PHO_Gene_Registry.txt"
      done
    done
  done
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$RUN_HOG_SELECTION" == "TRUE" ]
then
  storage_place=PHO_Selection_TPM
  rm $work_directory"/"$Result_Storage"/"$storage_place"/Selection_PHO_Table.txt"

  echo "PHO N_Species N_Wanted_Cestodes N_Wanted_Trematodes N_SL_Species N_SL_Species_High N_All_BAD N_Bad_No_SL N_Wanted_Cestodes_BAD N_Wanted_Trematodes_BAD" | tr " " "\t" >> $work_directory"/"$Result_Storage"/"$storage_place"/Selection_PHO_Table.txt"

  all_phos=$(awk '{print $1}' $work_directory"/"$Result_Storage"/"$storage_place"/PHO_Gene_Registry.txt" | grep -v "PHO" | sort -u)
  wanted_cestodes=$(echo "Emultilocularis Hmicrostoma Sparganum_proliferum Taenia_multiceps")
  wanted_trematodes=$(echo "Clonorchis_sinensis Fasciola_hepatica Opisthorchis_felineus Schistosoma_mansoni")

  echo $wanted_cestodes | tr " " "\n" > $work_directory"/"$Result_Storage"/"$storage_place"/Temporales/wanted_cestodes.tmp"
  echo $wanted_trematodes | tr " " "\n" > $work_directory"/"$Result_Storage"/"$storage_place"/Temporales/wanted_trematodes.tmp"

  for pho in $all_phos
  do
    grep -w $pho $work_directory"/"$Result_Storage"/"$storage_place"/PHO_Gene_Registry.txt" > $work_directory"/"$Result_Storage"/"$storage_place"/Temporales/current_pho.tmp"
    species_in_pho=$(awk -F "\t" '{print $2}' $work_directory"/"$Result_Storage"/"$storage_place"/Temporales/current_pho.tmp" | sort -u)
    N_species_in_pho=$(echo $species_in_pho | tr " " "\n" | grep -c .)
    N_wanted_cestodes=$(echo $species_in_pho | tr " " "\n" | grep -c -f $work_directory"/"$Result_Storage"/"$storage_place"/Temporales/wanted_cestodes.tmp")
    N_wanted_trematodes=$(echo $species_in_pho | tr " " "\n" | grep -c -f $work_directory"/"$Result_Storage"/"$storage_place"/Temporales/wanted_trematodes.tmp")
    N_SL_Species_All=$(awk -F "\t" '{if ($5=="SL") print $2}' $work_directory"/"$Result_Storage"/"$storage_place"/Temporales/current_pho.tmp" | sort -u | grep -c .)
    N_SL_Species_All_and_high=$(awk -F "\t" '{if ($5=="SL" && $6=="Yes") print $2}' $work_directory"/"$Result_Storage"/"$storage_place"/Temporales/current_pho.tmp" | sort -u | grep -c .)

    All_total_bad=0
    No_SL_And_BAD=0
    Wanted_cestodes_bad=0
    Wanted_trematodes_bad=0

    for spe in $species_in_pho
    do
      check_species=$(grep $spe $work_directory"/"$Result_Storage"/"$storage_place"/Temporales/current_pho.tmp" | awk -F "\t" '{if ($6=="Yes") print}' | grep -c .)
      check_wanted_cestodes=$(grep -c $spe $work_directory"/"$Result_Storage"/"$storage_place"/Temporales/wanted_cestodes.tmp")
      check_wanted_trematodes=$(grep -c $spe $work_directory"/"$Result_Storage"/"$storage_place"/Temporales/wanted_trematodes.tmp")

      check_species_no_sl_and_bad=$(grep $spe $work_directory"/"$Result_Storage"/"$storage_place"/Temporales/current_pho.tmp" | awk -F "\t" '{if ($5=="SL" && $6=="No") print}' | grep -c .)

      if [ $check_species -eq 0 ]
      then
        All_total_bad=$(($All_total_bad + 1))

        if [ $check_species_no_sl_and_bad -eq 0 ]
        then
          No_SL_And_BAD=$(($No_SL_And_BAD + 1))
          if [ $check_wanted_cestodes -gt 0 ]
          then
            Wanted_cestodes_bad=$(($Wanted_cestodes_bad + 1))
          elif [ $check_wanted_trematodes -gt 0 ]
          then
            Wanted_trematodes_bad=$(($Wanted_trematodes_bad + 1))
          fi
        fi
      fi
    done

    echo "pho:"$pho" || "$species_in_pho
    echo $pho" "$N_species_in_pho" "$N_wanted_cestodes" "$N_wanted_trematodes" "$N_SL_Species_All" "$N_SL_Species_All_and_high" "$All_total_bad" "$No_SL_And_BAD" "$Wanted_cestodes_bad" "$Wanted_trematodes_bad | tr " " "\t" >> $work_directory"/"$Result_Storage"/"$storage_place"/Selection_PHO_Table.txt"
  done
  # Criterio selección: 2 Cestodos; 2 Trematodos; 5 especies con SLs; Todas las especies con al menos un representante viable.
  head -n 1 $work_directory"/"$Result_Storage"/"$storage_place"/Selection_PHO_Table.txt" >> $work_directory"/"$Result_Storage"/"$storage_place"/General_Selection_PHO.txt"
  awk -F "\t" '{if ($3>=2 && $4>=2 && $5>=5 && $8==0) print}' $work_directory"/"$Result_Storage"/"$storage_place"/Selection_PHO_Table.txt" >> $work_directory"/"$Result_Storage"/"$storage_place"/General_Selection_PHO.txt"
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$RUN_HOG_SEQUENCES" == "TRUE" ]
then
  storage_place=PHO_Selection_TPM
  results_orthofinder=$(ls -d $work_directory"/"$Result_Storage"/OrthoFinder/Results_"*)
  PHO_File=$(echo $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")

  rm -r $work_directory"/"$Result_Storage"/"$storage_place"/Sequences"
  rm -r $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data"
  mkdir $work_directory"/"$Result_Storage"/"$storage_place"/Sequences"
  mkdir $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data"

  selected_PHOs=$(sed 1d $work_directory"/"$Result_Storage"/"$storage_place"/General_Selection_PHO.txt" | awk -F "\t" '{print $1}' | sort -u)
  interest_pho_sequences
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$RUN_HOG_PHYLOGENY" == "TRUE" ]
then
  storage_place=PHO_Selection_TPM

  rm -r $work_directory"/"$Result_Storage"/"$storage_place"/Phylogeny"
  mkdir $work_directory"/"$Result_Storage"/"$storage_place"/Phylogeny"
  mkdir $work_directory"/"$Result_Storage"/"$storage_place"/Phylogeny/Alignment"

  selected_PHOs=$(sed 1d $work_directory"/"$Result_Storage"/"$storage_place"/General_Selection_PHO.txt" | awk -F "\t" '{print $1}' | sort -u)
  interest_pho_phylogeny
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$RUN_INTERPRO" == "TRUE" ]
then
  # conda activate Interprot
  rm -r $work_directory"/"$Result_Storage"/Interprot_Results"
  mkdir $work_directory"/"$Result_Storage"/Interprot_Results"
  mkdir $work_directory"/"$Result_Storage"/Interprot_Results/Sequences"

  Especies=$(ls $work_directory"/"$Result_Storage"/"*pep | sed "s/.*\///" | sed "s/.pep//")
  for spe in $Especies
  do
    if [ ! -f $work_directory"/"$Result_Storage"/Interprot_Results/"$spe".tsv" ]
    then
      echo "Running: "$spe
      seqkit grep -v -s -p "*" $work_directory"/"$Result_Storage"/"$spe".pep" >> $work_directory"/"$Result_Storage"/Interprot_Results/Sequences/"$spe"_interprot.aa"

      $interprot_path"/interproscan.sh" -f TSV -pa -b $work_directory"/"$Result_Storage"/Interprot_Results/"$spe -dra -i $work_directory"/"$Result_Storage"/Interprot_Results/Sequences/"$spe"_interprot.aa"
      echo ""
      echo "################################################################################"
      echo ""
    fi
  done
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$RUN_OPERON_SEARCH" == "TRUE" ]
then
  results_orthofinder=$(ls -d $work_directory"/"$Result_Storage"/OrthoFinder/Results_"*)
  PHO_File=$(echo $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")

  rm -r $work_directory"/"$Result_Storage"/Operons"
  mkdir $work_directory"/"$Result_Storage"/Operons"
  mkdir $work_directory"/"$Result_Storage"/Operons/Temporal"

  grupo=Cestodes
  Especies=$(echo "Egranulosus Emultilocularis Hdiminuta Hmicrostoma Mesocestoides_corti Sparganum_proliferum Spirometra_erinaceieuropaei Taenia_asiatica Taenia_multiceps Taenia_saginata Taenia_solium Schistocephalus_solidus")
  find_operons

  grupo=Trematoda
  Especies=$(echo "Clonorchis_sinensis Fasciola_gigantica Fasciola_hepatica Fasciolopsis_buski Opisthorchis_felineus Paragonimus_heterotremus Paragonimus_westermani Schistosoma_bovis Schistosoma_haematobium Schistosoma_japonicum Schistosoma_mansoni Trichobilharzia_regenti")
  find_operons
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

# Requiere arreglos importantes en varias preguntas

if [ "$RUN_OPERON_CONSERVATION_QUESTIONS_1" == "TRUE" ]
then
  rm -r $work_directory"/"$Result_Storage"/Operons/Questions"
  mkdir $work_directory"/"$Result_Storage"/Operons/Questions"
  mkdir $work_directory"/"$Result_Storage"/Operons/Questions/Temp"
  storage_place=PHO_Selection_TPM

  Especies=$(echo "Clonorchis_sinensis Egranulosus Emultilocularis Fasciola_gigantica Fasciola_hepatica Fasciolopsis_buski Hdiminuta Hmicrostoma Mesocestoides_corti Opisthorchis_felineus Paragonimus_heterotremus Paragonimus_westermani Schistocephalus_solidus Schistosoma_bovis Schistosoma_haematobium Schistosoma_japonicum Schistosoma_mansoni Sparganum_proliferum Spirometra_erinaceieuropaei Taenia_asiatica Taenia_multiceps Taenia_saginata Taenia_solium Trichobilharzia_regenti")
  Grupos_Cestoda=$(echo "Egranulosus Emultilocularis Hdiminuta Hmicrostoma Mesocestoides_corti Schistocephalus_solidus Sparganum_proliferum Spirometra_erinaceieuropaei Taenia_asiatica Taenia_multiceps Taenia_saginata Taenia_solium")
  Grupos_Cestoda_Derivados=$(echo "Egranulosus Emultilocularis Hdiminuta Hmicrostoma Mesocestoides_corti Taenia_asiatica Taenia_multiceps Taenia_saginata Taenia_solium")
  Grupos_Cestoda_Basal=$(echo "Schistocephalus_solidus Sparganum_proliferum Spirometra_erinaceieuropaei")
  Grupos_Taenidae=$(echo "Egranulosus Emultilocularis Taenia_asiatica Taenia_multiceps Taenia_saginata Taenia_solium")
  Grupos_Hymenolepis=$(echo "Hdiminuta Hmicrostoma")

  ### Trematoda
  Grupos_Trematoda=$(echo "Clonorchis_sinensis Fasciola_gigantica Fasciola_hepatica Fasciolopsis_buski Opisthorchis_felineus Paragonimus_heterotremus Paragonimus_westermani Schistosoma_bovis Schistosoma_haematobium Schistosoma_japonicum Schistosoma_mansoni Trichobilharzia_regenti")
  Grupos_Schistosmidae=$(echo "Schistosoma_bovis Schistosoma_haematobium Schistosoma_japonicum Schistosoma_mansoni Trichobilharzia_regenti")
  Grupos_Trematoda_No_schistosoma=$(echo "Clonorchis_sinensis Fasciola_gigantica Fasciola_hepatica Fasciolopsis_buski Opisthorchis_felineus Paragonimus_heterotremus Paragonimus_westermani")
  Grupos_Fasciolidae=$(echo "Fasciola_gigantica Fasciola_hepatica Fasciolopsis_buski")
  Grupos_Trematoda_NoSchistoFasciola=$(echo "Clonorchis_sinensis Opisthorchis_felineus Paragonimus_heterotremus Paragonimus_westermani")
  Grupos_Paragonimus=$(echo "Paragonimus_heterotremus Paragonimus_westermani")
  Grupos_Schistosoma=$(echo "Schistosoma_bovis Schistosoma_haematobium Schistosoma_japonicum Schistosoma_mansoni")
  Grupos_Special=$(echo "Clonorchis_sinensis Fasciola_gigantica Fasciola_hepatica Fasciolopsis_buski Opisthorchis_felineus Paragonimus_heterotremus Paragonimus_westermani Schistosoma_bovis Schistosoma_haematobium Schistosoma_japonicum Schistosoma_mansoni Trichobilharzia_regenti Schistocephalus_solidus Sparganum_proliferum Spirometra_erinaceieuropaei")

  # 0) Generando set de palabras basico

  awk -F "\t" '{print $7}' $work_directory"/"$Result_Storage"/Operons/"*"_operon.txt" | grep -v PHOs | sort -u >> $work_directory"/"$Result_Storage"/Operons/Questions/Temp/raw_gene_pairs.temp"

  count=1
  recorrer=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Questions/Temp/raw_gene_pairs.temp")

  while [ $count -le $recorrer ]
  do
    first_of_pair=1
    last_of_pair=2
    all_members=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Operons/Questions/Temp/raw_gene_pairs.temp" | sed "s/__/\n/g" | grep -c .)

    while [ $last_of_pair -le $all_members ]
    do
      first_one=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Operons/Questions/Temp/raw_gene_pairs.temp" | awk -F "__" -v N=$first_of_pair '{print $N}')
      last_one=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Operons/Questions/Temp/raw_gene_pairs.temp" | awk -F "__" -v N=$last_of_pair '{print $N}')

      if [ "$first_one" != "No_Hits" ] && [ "$last_one" != "No_Hits" ]
      then
        echo $first_one"__"$last_one >> $work_directory"/"$Result_Storage"/Operons/Questions/Temp/refined_gene_pairs.temp"
      fi
      first_of_pair=$(($first_of_pair + 1))
      last_of_pair=$(($last_of_pair + 1))
    done
    count=$(($count + 1))
  done

  sort -u $work_directory"/"$Result_Storage"/Operons/Questions/Temp/refined_gene_pairs.temp" | awk '{print "Op_PHO_Pair-"NR"\t"$0}' >> $work_directory"/"$Result_Storage"/Operons/Questions/Work_gene_pairs.txt"

  #####################################################################################################################################
  #####################################################################################################################################
  #####################################################################################################################################

  sed "s/\t/__/" $work_directory"/"$Result_Storage"/Operons/Questions/Work_gene_pairs.txt" | awk -F "__" '{if ($2==$3) print}' | sed "s/__/\t/" > $work_directory"/"$Result_Storage"/Operons/Questions/Tandem_repeat_pairs.txt"
  # 1) Pairs with tandem_repeats:
  tandem_pairs=$(awk -F "\t" '{print $2}' $work_directory"/"$Result_Storage"/Operons/Questions/Tandem_repeat_pairs.txt")
  echo "Pairs in Tandem"
  echo "Tandem_Pair "$Especies | tr " " "\t">> $work_directory"/"$Result_Storage"/Operons/Questions/Tandem_repeat_pairs_Operons.txt"
  for tandem_pair in $tandem_pairs
  do
    registro_repetidos=$tandem_pair
    for spe in $Especies
    do
      N_operons=$(grep -c $tandem_pair $work_directory"/"$Result_Storage"/Operons/"$spe"_operon.txt")
      registro_repetidos=$(echo $registro_repetidos" "$N_operons)
    done
    echo $registro_repetidos | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Tandem_repeat_pairs_Operons.txt"
  done

  awk '{print $2}' $work_directory"/"$Result_Storage"/Operons/Questions/Tandem_repeat_pairs.txt" >> $work_directory"/"$Result_Storage"/Operons/Questions/Tandem_repeat_pairs.temp"
  echo "Species N_Tandem_operons" | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Tandem_repeat_counts_Operons.txt"
  for spe in $Especies
  do
    tandem_counts=$(grep -c -f $work_directory"/"$Result_Storage"/Operons/Questions/Tandem_repeat_pairs.temp" $work_directory"/"$Result_Storage"/Operons/"$spe"_operon.txt")
    echo $spe" tandem_counts:"$tandem_counts
    echo $spe" "$tandem_counts | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Tandem_repeat_counts_Operons.txt"
  done
  rm $work_directory"/"$Result_Storage"/Operons/Questions/Tandem_repeat_pairs.temp"

  #####################################################################################################################################
  #####################################################################################################################################
  #####################################################################################################################################

  # 2) Absence/Presence of words per species
  pho_unit_of_interest=$(awk -F "\t" '{print $2}' $work_directory"/"$Result_Storage"/Operons/Questions/Work_gene_pairs.txt" | sort -u)
  presence_absence_results=Presence_or_Absence_By_Pairs.txt
  presence_absence

  pho_unit_of_interest=$(awk -F "\t" '{print $2}' $work_directory"/"$Result_Storage"/Operons/Questions/Work_gene_pairs.txt" | sed  "s/__/\n/" | sort -u)
  presence_absence_results=Presence_or_Absence_By_PHO.txt
  presence_absence

  #####################################################################################################################################
  #####################################################################################################################################
  #####################################################################################################################################

  # 3) Posiciones particulares de PHOs en los pares (la primer posición es muy wibbly wobbly pero debería tener una señal más clara en la segunda)
  all_phos=$(awk -F "\t" '{print $2}' $work_directory"/"$Result_Storage"/Operons/Questions/Work_gene_pairs.txt" | sed "s/__/\n/g" | sort -u)
  echo "PHO Total Primero Segundo Double Percentage_First Percentage_Second"  | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/PHO_in_the_pairs.txt"
  for pho in $all_phos
  do
    total_pairs=$(grep -c $pho $work_directory"/"$Result_Storage"/Operons/Questions/Work_gene_pairs.txt")
    first_only=$(awk -F "\t" '{print $2}' $work_directory"/"$Result_Storage"/Operons/Questions/Work_gene_pairs.txt" | awk -F "__" '{if ($1!=$2) print $1}' | grep -c $pho)
    second_only=$(awk -F "\t" '{print $2}' $work_directory"/"$Result_Storage"/Operons/Questions/Work_gene_pairs.txt" | awk -F "__" '{if ($1!=$2) print $2}' | grep -c $pho)
    double_tap=$(awk -F "\t" '{print $2}' $work_directory"/"$Result_Storage"/Operons/Questions/Work_gene_pairs.txt" | awk -F "__" '{if ($1==$2) print}' | grep -c $pho)

    if [ $total_pairs -eq $double_tap ]
    then
      percentage_first=NA
      percentage_second=NA
    else
      percentage_first=$(echo $first_only" "$total_pairs" "$double_tap | awk -F " " '{printf "%.2f\n", ($1/($2-$3))*100}')
      percentage_second=$(echo $second_only" "$total_pairs" "$double_tap | awk -F " " '{printf "%.2f\n", ($1/($2-$3))*100}')
    fi
    echo $pho" "$total_pairs" "$first_only" "$second_only" "$double_tap" "$percentage_first" "$percentage_second
    echo $pho" "$total_pairs" "$first_only" "$second_only" "$double_tap" "$percentage_first" "$percentage_second | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/PHO_in_the_pairs.txt"
  done

  # 4) Selected PHOs in Operons
   selected_PHOs=$(sed 1d $work_directory"/"$Result_Storage"/"$storage_place"/General_Selection_PHO.txt" | awk -F "\t" '{print $1}' | sort -u)
   echo "Selected_PHO "$Especies" All" | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Selected_PHOs_In_operons.txt"
   for pho in $selected_PHOs
   do
     echo $pho" Selected in Operon"
     registro_conteos=$pho
     all_spe=0
     for spe in $Especies
     do
       count_pho=$(grep -c $pho $work_directory"/"$Result_Storage"/Operons/"$spe"_operon.txt" | awk '{if ($1>0) print 1; else print 0}')
       registro_conteos=$(echo $registro_conteos" "$count_pho)
       all_spe=$(($all_spe + $count_pho))
     done
     echo $registro_conteos" "$all_spe | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Selected_PHOs_In_operons.txt"
   done
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$RUN_OPERON_CONSERVATION_QUESTIONS_2" == "TRUE" ]
then
  results_orthofinder=$(ls -d $work_directory"/"$Result_Storage"/OrthoFinder/Results_"*)
  PHO_File=$(echo $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")

  rm -r $work_directory"/"$Result_Storage"/Operons/Questions/Specific"
  mkdir $work_directory"/"$Result_Storage"/Operons/Questions/Specific"
  mkdir $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check"

  # 1) Uniq PHOs pairs
  count_species=4 # Primer columna con especies
  recorrer_species=27 # ultima columna con especies

  echo "Species Uniq_PHO_Pairs" | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Total_uniq_PHO_pairs.txt"
  while [ $count_species -le $recorrer_species ]
  do
    species_name=$(head -n 1 $work_directory"/"$Result_Storage"/Operons/Questions/Presence_or_Absence_By_Pairs.txt" | awk -F "\t" -v N=$count_species '{print $N}' )
    uniq_pairs=$(awk -F "\t" -v N=$count_species '{if ($2==1) print $N}'  $work_directory"/"$Result_Storage"/Operons/Questions/Presence_or_Absence_By_Pairs.txt" | grep -w -c "X")

    echo $species_name" "$uniq_pairs | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Total_uniq_PHO_pairs.txt"

    count_species=$(($count_species + 1))
  done

  ###########################################################################################################################################

  # 2) Group specific PHOs RAW # It is a basic count of how many species have the pair, regardless of the PHO
  echo "PHO_Pair Total_Species Group" | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_pairs_raw.txt"
  awk -F "\t" '{if ($28>=5 && $33==0) print $1"\t"$2"\tCestode_Exclusive"}' $work_directory"/"$Result_Storage"/Operons/Questions/Presence_or_Absence_By_Pairs.txt" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_pairs_raw.txt"
  awk -F "\t" '{if ($28==0 && $33>=5) print $1"\t"$2"\tTrematode_Exclusive"}' $work_directory"/"$Result_Storage"/Operons/Questions/Presence_or_Absence_By_Pairs.txt" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_pairs_raw.txt"
  awk -F "\t" '{if ($28>=2 && $33>=2 && $28+$33>=10) print $1"\t"$2"\tBalanced"}' $work_directory"/"$Result_Storage"/Operons/Questions/Presence_or_Absence_By_Pairs.txt" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_pairs_raw.txt"

  ###########################################################################################################################################

  # 3) Group specific PHOs Corrected by Phylogeny # It is a basic count of how many species have the pair, regardless of the PHO
  echo "PHO_Pair Total_Cestodes Total_Trematodes Echinococus Hymenolepis Mesocestoides Schistocephalus Sparganum Spirometra Taenia Clonorchis Fasciola Fasciolopis Opisthorchis Paragonimus Schistosoma Trichobilharzia Selected_Reason" | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_pairs_Only_Genus.txt"
  interest_pair=$(grep -v "PHO_Pair" $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_pairs_raw.txt" | awk -F "\t" '{print $1}')
  for pair in $interest_pair
  do
    selected_reason=$(grep $pair $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_pairs_raw.txt" | awk -F "\t" '{print $3}')

    Presence_in_Echinococus=$(grep -w $pair $work_directory"/"$Result_Storage"/Operons/Questions/Presence_or_Absence_By_Pairs.txt" | awk -F "\t" '{print $5" "$6}' | grep -c "X")
    Presence_in_Hymenolepis=$(grep -w $pair $work_directory"/"$Result_Storage"/Operons/Questions/Presence_or_Absence_By_Pairs.txt" | awk -F "\t" '{print $10" "$11}' | grep -c "X")
    Presence_in_Mesocestoides_corti=$(grep -w $pair $work_directory"/"$Result_Storage"/Operons/Questions/Presence_or_Absence_By_Pairs.txt" | awk -F "\t" '{print $12}' | grep -c "X")
    Presence_in_Schistocephalus_solidus=$(grep -w $pair $work_directory"/"$Result_Storage"/Operons/Questions/Presence_or_Absence_By_Pairs.txt" | awk -F "\t" '{print $16}' | grep -c "X")
    Presence_in_Sparganum_proliferum=$(grep -w $pair $work_directory"/"$Result_Storage"/Operons/Questions/Presence_or_Absence_By_Pairs.txt" | awk -F "\t" '{print $21}' | grep -c "X")
    Presence_in_Spirometra_erinaceieuropaei=$(grep -w $pair $work_directory"/"$Result_Storage"/Operons/Questions/Presence_or_Absence_By_Pairs.txt" | awk -F "\t" '{print $22}' | grep -c "X")
    Presence_in_Taenia=$(grep -w $pair $work_directory"/"$Result_Storage"/Operons/Questions/Presence_or_Absence_By_Pairs.txt" | awk -F "\t" '{print $23" "$24" "$25" "$26}' | grep -c "X")
    Total_Cestodos=$(($Presence_in_Echinococus + $Presence_in_Hymenolepis + $Presence_in_Mesocestoides_corti + $Presence_in_Schistocephalus_solidus + $Presence_in_Sparganum_proliferum + $Presence_in_Spirometra_erinaceieuropaei + $Presence_in_Taenia))

    Presence_in_Clonorchis_sinensis=$(grep -w $pair $work_directory"/"$Result_Storage"/Operons/Questions/Presence_or_Absence_By_Pairs.txt" | awk -F "\t" '{print $4}' | grep -c "X")
    Presence_in_Fasciola=$(grep -w $pair $work_directory"/"$Result_Storage"/Operons/Questions/Presence_or_Absence_By_Pairs.txt" | awk -F "\t" '{print $7" "$8}' | grep -c "X")
    Presence_in_Fasciolopis=$(grep -w $pair $work_directory"/"$Result_Storage"/Operons/Questions/Presence_or_Absence_By_Pairs.txt" | awk -F "\t" '{print $9}' | grep -c "X")
    Presence_in_Opisthorchis_felineus=$(grep -w $pair $work_directory"/"$Result_Storage"/Operons/Questions/Presence_or_Absence_By_Pairs.txt" | awk -F "\t" '{print $13}' | grep -c "X")
    Presence_in_Paragonimus=$(grep -w $pair $work_directory"/"$Result_Storage"/Operons/Questions/Presence_or_Absence_By_Pairs.txt" | awk -F "\t" '{print $14" "$15}' | grep -c "X")
    Presence_in_Schistosoma=$(grep -w $pair $work_directory"/"$Result_Storage"/Operons/Questions/Presence_or_Absence_By_Pairs.txt" | awk -F "\t" '{print $17" "$18" "$19" "$20}' | grep -c "X")
    Presence_in_Trichobilharzia_regenti=$(grep -w $pair $work_directory"/"$Result_Storage"/Operons/Questions/Presence_or_Absence_By_Pairs.txt" | awk -F "\t" '{print $27}' | grep -c "X")
    Total_Trematodos=$(($Presence_in_Clonorchis_sinensis + $Presence_in_Fasciola + $Presence_in_Fasciolopis + $Presence_in_Opisthorchis_felineus + $Presence_in_Paragonimus + $Presence_in_Schistosoma + $Presence_in_Trichobilharzia_regenti))

    echo $pair" "$Total_Cestodos" "$Total_Trematodos" "$Presence_in_Echinococus" "$Presence_in_Hymenolepis" "$Presence_in_Mesocestoides_corti" "$Presence_in_Schistocephalus_solidus" "$Presence_in_Sparganum_proliferum" "$Presence_in_Spirometra_erinaceieuropaei" "$Presence_in_Taenia" "$Presence_in_Clonorchis_sinensis" "$Presence_in_Fasciola" "$Presence_in_Fasciolopis" "$Presence_in_Opisthorchis_felineus" "$Presence_in_Paragonimus" "$Presence_in_Schistosoma" "$Presence_in_Trichobilharzia_regenti" "$selected_reason | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_pairs_Only_Genus.txt"
  done

  ###########################################################################################################################################

  # 4) Evaluar cambios sintenicos en de los operones de interes cuando no son encontrados
  recorrer_species=27 # ultima columna con especies

  result_storage=Cestode_Exclusive
  interest_pairs=$(grep -w "Cestode_Exclusive" $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_pairs_Only_Genus.txt" | awk -F "\t" '{if ($2+$3>=4) print $1}')
  confirm_pair_synteny

  result_storage=Trematode_Exclusive
  interest_pairs=$(grep -w "Trematode_Exclusive" $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_pairs_Only_Genus.txt" | awk -F "\t" '{if ($2+$3>=4) print $1}')
  confirm_pair_synteny

  result_storage=Balanced
  interest_pairs=$(grep -w "Balanced" $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Interest_Operon_pairs_Only_Genus.txt" | awk -F "\t" '{if ($2+$3>=4) print $1}')
  confirm_pair_synteny

  result_storage=Ganancia_Hymenolepis
  interest_pairs=$(echo "N0.HOG0004578__N0.HOG0008654 N0.HOG0004578__N0.HOG0004577 N0.HOG0004578__N0.HOG0008554 N0.HOG0012528__N0.HOG0004577")
  confirm_pair_synteny

  result_storage=Operones_Paper_Hymenolepis
  interest_pairs=$(echo "N0.HOG0003395__N0.HOG0003396 N0.HOG0003757__N0.HOG0003758 N0.HOG0008711__N0.HOG0008107 N0.HOG0011155__N0.HOG0014589 N0.HOG0011749__N0.HOG0009965 N0.HOG0012340__N0.HOG0009961 N0.HOG0012359__N0.HOG0006427 N0.HOG0013560__N0.HOG0009782")
  confirm_pair_synteny

  result_storage=PHOs_Seleccionados
  sed 1d $work_directory"/"$Result_Storage"/PHO_Selection_TPM/General_Selection_PHO.txt" | awk -F "\t" '{print $1}' | sort -u >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Selected_PHO.temp"
  interest_pairs=$(grep -f $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Selected_PHO.temp" $work_directory"/"$Result_Storage"/Operons/Questions/Work_gene_pairs.txt" | awk -F "\t" '{print $2}')
  rm $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Selected_PHO.temp"
  confirm_pair_synteny

  result_storage=Selected_By_Association
  interest_pairs=$(echo "N0.HOG0011496__N0.HOG0010248 N0.HOG0010460__N0.HOG0015076 N0.HOG0002658__N0.HOG0008527 N0.HOG0025669__N0.HOG0004267 N0.HOG0006436__N0.HOG0019938 N0.HOG0014388__N0.HOG0011302 N0.HOG0000286__N0.HOG0012130 N0.HOG0010740__N0.HOG0015684 N0.HOG0011332__N0.HOG0003821 N0.HOG0019220__N0.HOG0006207 N0.HOG0003428__N0.HOG0000049 N0.HOG0003427__N0.HOG0008150 N0.HOG0003427__N0.HOG0003428 N0.HOG0005019__N0.HOG0012341 N0.HOG0000286__N0.HOG0009961 N0.HOG0014511__N0.HOG0009668 N0.HOG0019884__N0.HOG0010962 N0.HOG0015070__N0.HOG0006584 N0.HOG0015076__N0.HOG0010460 N0.HOG0003427__N0.HOG0008149 N0.HOG0003395__N0.HOG0019656 N0.HOG0004577__N0.HOG0009987 N0.HOG0009961__N0.HOG0000286 N0.HOG0005523__N0.HOG0005058 N0.HOG0005523__N0.HOG0000286")
  confirm_pair_synteny
fi

#######################################################################################################################################################################################################
#######################################################################################################################################################################################################
#######################################################################################################################################################################################################

if [ "$RUN_OPERON_CONSERVATION_QUESTIONS_3" == "TRUE" ]
then
  results_orthofinder=$(ls -d $work_directory"/"$Result_Storage"/OrthoFinder/Results_"*)
  PHO_File=$(echo $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")

  rm -r $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Final_Sumary"
  mkdir $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Final_Sumary"

  # Hacer una tabla con cada par revisado y...
  # A) Registrar si hay SLs en alguno de los genes identificados como viables, y solapados con genes (y registrar quienes)
  # B) Confirmar que el orden de PHOs revisado sea valido
  # C) Registrar si el par fue utilizado en otro grupo de pares de interes

  explore_sumaries=$(ls $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/Summary_"*".txt" | sed "s/.*\///" | sed "s/Summary_//" | sed "s/.txt//")
  Especies=$(ls $work_directory"/"$Result_Storage"/"*pep | sed "s/.*\///" | sed "s/.pep//")

  echo "Grupo Par N_Especies_SL_Original Especies_SL_Original Genes_SL_Original Chimeric_Original N_Especies_SL_Extendido Especies_SL_Extendido Genes_SL_Extendido Chimeric_Extendido Cambios_Orden Usado_Tambien_En" | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Final_Sumary/Operon_Pairs_Final.txt"
  for summary in $explore_sumaries
  do
    pho_pairs=$(awk -F "\t" '{print $1}' $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/Summary_"$summary".txt" | grep -v -w "PHO_Pair" | sort -u)

    for pair in $pho_pairs
    do
      echo $summary" "$pair
      first_pho=$(echo $pair | awk -F "__" '{print $1}')
      second_pho=$(echo $pair | awk -F "__" '{print $2}')
      # Variables de Registro
      registro_otros_sumarios=NA

      ajustar_SL_original=FALSE
      ajustar_chimera_original=FALSE
      registro_especies_SL_original=NA
      registro_genes_SL_original=NA
      chimera_original=NA
      total_SL_original=0

      ajustar_SL_extendido=FALSE
      ajustar_chimera_extendido=FALSE
      registro_especies_SL_extendido=NA
      registro_genes_SL_extendido=NA
      chimera_extendido=NA
      total_SL_extendido=0

      ajustar_cambio_orden=FALSE
      confirmar_cambio_orden=NA

      check_other_summay=$(grep -w -c $pair $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/Summary_"*".txt" | grep -v "Summary_"$summary".txt" | awk -F ":" '{if ($2>0) print}' | grep -c .)
      if [ $check_other_summay -gt 0 ]
      then
        registro_otros_sumarios=$(grep -w -c $pair $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/Summary_"*".txt" | grep -v "Summary_"$summary".txt" | awk -F ":" '{if ($2>0) print $1}' | sed "s/.*\///" | sed "s/Summary_//" | sed "s/.txt//" | tr "\n" ";" | sed "s/;$//")
      fi

      for spe in $Especies
      do
        # 1) Fijarse si es uno de los pares originales
        check_original_pair=$(grep -c $pair $work_directory"/"$Result_Storage"/Operons/"$spe"_operon.txt")

        if [ $check_original_pair -gt 0 ]
        then
          operon_ids=$(grep $pair $work_directory"/"$Result_Storage"/Operons/"$spe"_operon.txt" | awk -F "\t" '{print $1}' | sort -u)
          for op in $operon_ids
          do
            First_current_genes_SL=$(grep -w -F $op $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | grep $first_pho | awk -F "\t" '{if ($6=="SL") print $1}')
            First_current_chimera_SL=$(grep -w -F $op $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | grep $first_pho | awk -F "\t" '{if ($6=="SL") print $1}')

            Second_current_genes_SL=$(grep -w -F $op $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | grep $second_pho | awk -F "\t" '{if ($6=="SL") print $1}')
            registro_genes_SL_original=$(echo $registro_genes_SL_original $First_current_genes_SL $Second_current_genes_SL)

            check_chimera=$(grep -w -F $op $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | grep -c -w Chimeric )
            if [ $check_chimera -gt 0 ]
            then
              ajustar_chimera_original=TRUE
              chimera_original=$(echo $chimera_original" "$spe )
            fi
          done
          check_genes_SL=$(echo $registro_genes_SL_original | tr " " "\n" | grep -w -v NA | grep -c .)
          if [ $check_genes_SL -gt 0 ]
          then
            ajustar_SL_original=TRUE
            registro_especies_SL_original=$(echo $registro_especies_SL_original" "$spe)
          fi
        elif [ -f $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$summary"_details/"$pair"_"$spe"_"$summary"_Gene_Details.txt" ]
        then
          count=2
          recorrer=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$summary"_details/"$pair"_"$spe"_"$summary"_Gene_Details.txt")
          while [ $count -le $recorrer ]
          do
            check_gene_pairing=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$summary"_details/"$pair"_"$spe"_"$summary"_Gene_Details.txt" | awk -F "\t" '{print $5}')
            if [ $check_gene_pairing == "Interruping_Genes" ] || [ $check_gene_pairing == "Viable_Operon" ]
            then
              temp_first_gene=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$summary"_details/"$pair"_"$spe"_"$summary"_Gene_Details.txt" | awk -F "\t" '{print $3}' | sed "s/(.*//")
              temp_second_gene=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$summary"_details/"$pair"_"$spe"_"$summary"_Gene_Details.txt" | awk -F "\t" '{print $4}' | sed "s/(.*//")

              check_chimera=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$summary"_details/"$pair"_"$spe"_"$summary"_Gene_Details.txt" | grep -c _Chimera_)

              if [ $check_chimera -gt 0 ]
              then
                ajustar_chimera_extendido=TRUE
                chimera_extendido=$(echo $chimera_extendido" "$spe)
              fi

              # 2) Confirmar SLs
              check_gene_sl=$(grep -w -F $temp_first_gene $work_directory"/"$Result_Storage"/PHO_Selection_TPM/PHO_Gene_Registry.txt" | grep -w -c SL)
              if [ $check_gene_sl -gt 0 ]
              then
                ajustar_SL_extendido=TRUE
                registro_genes_SL_extendido=$(echo $registro_genes_SL_extendido" "$temp_first_gene)
                registro_especies_SL_extendido=$(echo $registro_especies_SL_extendido" "$spe)
              fi

              check_gene_sl=$(grep -w -F $temp_second_gene $work_directory"/"$Result_Storage"/PHO_Selection_TPM/PHO_Gene_Registry.txt" | grep -w -c SL)

              if [ $check_gene_sl -gt 0 ]
              then
                ajustar_SL_extendido=TRUE
                registro_genes_SL_extendido=$(echo $registro_genes_SL_extendido" "$temp_second_gene)
                registro_especies_SL_extendido=$(echo $registro_especies_SL_extendido" "$spe)
              fi

              # 3) Confirmar orden en el par
              grep -w -F $temp_first_gene $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$pair"_"*"_"$spe"_gene_coordinates.txt" | awk -F ":" '{print $2}' >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Final_Sumary/temp_gene_coords"
              grep -w -F $temp_second_gene $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Synteny_Check/"$pair"_"*"_"$spe"_gene_coordinates.txt" | awk -F ":" '{print $2}' >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Final_Sumary/temp_gene_coords"

              strand=$(awk -F "\t" '{print $5}' $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Final_Sumary/temp_gene_coords" | sort -u)

              if [ $strand == "+" ]
              then
                check_pho_first_gene=$(sort -n -k 3 $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Final_Sumary/temp_gene_coords" | head -n 1 | awk -F "\t" '{print $1}')
              else
                check_pho_first_gene=$(sort -n -r -k 4 $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Final_Sumary/temp_gene_coords" | head -n 1 | awk -F "\t" '{print $1}')
              fi

              check_pho_order=$(grep -w -F $check_pho_first_gene $PHO_File | awk -F "\t" '{print $1}')
              if [ $check_pho_order != $first_pho ]
              then
                ajustar_cambio_orden=TRUE
                confirmar_cambio_orden=$(echo $confirmar_cambio_orden" "$spe )
              fi
              rm $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Final_Sumary/temp_gene_coords"
            fi
            count=$(($count + 1))
          done
        fi
      done

      if [ "$ajustar_SL_original" == "TRUE" ]
      then
        registro_genes_SL_original=$(echo $registro_genes_SL_original | tr " " "\n" | grep -w -v NA | grep . | sort -u | tr "\n" ";" | sed "s/;$//")
        registro_especies_SL_original=$(echo $registro_especies_SL_original | tr " " "\n" | grep -w -v NA | grep . | sort -u | tr "\n" ";" | sed "s/;$//")
        total_SL_original=$(echo $registro_especies_SL_original | tr ";" "\n" | grep -c . )
      fi

      if [ "$ajustar_chimera_original" == "TRUE" ]
      then
        chimera_original=$(echo $chimera_original | tr " " "\n" | grep -w -v NA | grep . | sort -u | tr "\n" ";" | sed "s/;$//")
      fi

      if [ "$ajustar_SL_extendido" == "TRUE" ]
      then
        registro_genes_SL_extendido=$(echo $registro_genes_SL_extendido | tr " " "\n" | grep -w -v NA | grep . | sort -u | tr "\n" ";" | sed "s/;$//")
        registro_especies_SL_extendido=$(echo $registro_especies_SL_extendido | tr " " "\n" | grep -w -v NA | grep . | sort -u | tr "\n" ";" | sed "s/;$//")
        total_SL_extendido=$(echo $registro_especies_SL_extendido | tr ";" "\n" | grep -c . )
      fi
      if [ "$ajustar_chimera_extendido" == "TRUE" ]
      then
        chimera_extendido=$(echo $chimera_extendido  | tr " " "\n" | grep -w -v NA | grep . | sort -u | tr "\n" ";" | sed "s/;$//")
      fi

      if [ "$ajustar_cambio_orden" == "TRUE" ]
      then
        confirmar_cambio_orden=$(echo $confirmar_cambio_orden | tr " " "\n" | grep -w -v NA | grep . | sort -u | tr "\n" ";" | sed "s/;$//" )
      fi
      echo $summary" "$pair" "$total_SL_original" "$registro_especies_SL_original" "$registro_genes_SL_original" "$chimera_original" "$total_SL_extendido" "$registro_especies_SL_extendido" "$registro_genes_SL_extendido" "$chimera_extendido" "$confirmar_cambio_orden" "$registro_otros_sumarios | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Final_Sumary/Operon_Pairs_Final.txt"
    done
  done
fi

#######################################################################################################################################################################################################
#######################################################################################################################################################################################################
#######################################################################################################################################################################################################

if [ $RUN_OPERON_CONSERVATION_QUESTIONS_4 == "TRUE" ]
then
  # Intentar confirmar el rango de PHOs que se comparten entre los distintos grupos
  rm -r $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO"
  mkdir $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO"
  mkdir $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/Pairwise_marches"

  Grupos_Cestoda=$(echo "Egranulosus Emultilocularis Hdiminuta Hmicrostoma Mesocestoides_corti Schistocephalus_solidus Sparganum_proliferum Spirometra_erinaceieuropaei Taenia_asiatica Taenia_multiceps Taenia_saginata Taenia_solium")
  Grupos_Cestoda_Derivados=$(echo "Egranulosus Emultilocularis Hdiminuta Hmicrostoma Mesocestoides_corti Taenia_asiatica Taenia_multiceps Taenia_saginata Taenia_solium")
  Grupos_Cestoda_Basal=$(echo "Schistocephalus_solidus Sparganum_proliferum Spirometra_erinaceieuropaei")
  Grupos_Taenidae=$(echo "Egranulosus Emultilocularis Taenia_asiatica Taenia_multiceps Taenia_saginata Taenia_solium")
  Grupos_Hymenolepis=$(echo "Hdiminuta Hmicrostoma")
  Grupos_Echinococus=$(echo "Egranulosus Emultilocularis")
  Grupos_Taenia=$(echo "Taenia_asiatica Taenia_multiceps Taenia_saginata Taenia_solium")

  Grupos_Trematoda=$(echo "Clonorchis_sinensis Fasciola_gigantica Fasciola_hepatica Fasciolopsis_buski Opisthorchis_felineus Paragonimus_heterotremus Paragonimus_westermani Schistosoma_bovis Schistosoma_haematobium Schistosoma_japonicum Schistosoma_mansoni Trichobilharzia_regenti")
  Grupos_Schistosmidae=$(echo "Schistosoma_bovis Schistosoma_haematobium Schistosoma_japonicum Schistosoma_mansoni Trichobilharzia_regenti")
  Grupos_Trematoda_No_schistosoma=$(echo "Clonorchis_sinensis Fasciola_gigantica Fasciola_hepatica Fasciolopsis_buski Opisthorchis_felineus Paragonimus_heterotremus Paragonimus_westermani")
  Grupos_Fasciolidae=$(echo "Fasciola_gigantica Fasciola_hepatica Fasciolopsis_buski")
  Grupos_Trematoda_NoSchistoFasciola=$(echo "Clonorchis_sinensis Opisthorchis_felineus Paragonimus_heterotremus Paragonimus_westermani")
  Grupos_Trematoda_Clonorchis_Opisthorchis=$(echo "Clonorchis_sinensis Opisthorchis_felineus")
  Grupos_Paragonimus=$(echo "Paragonimus_heterotremus Paragonimus_westermani")
  Grupos_Schistosoma=$(echo "Schistosoma_bovis Schistosoma_haematobium Schistosoma_japonicum Schistosoma_mansoni")
  Grupos_Fasciola=$(echo "Fasciola_gigantica Fasciola_hepatica")

  wanted_species=$Grupos_Cestoda
  group_name=Grupos_Cestoda
  get_pho_info_phos

  wanted_species=$Grupos_Cestoda_Derivados
  group_name=Grupos_Cestoda_Derivados
  get_pho_info_phos

  wanted_species=$Grupos_Cestoda_Basal
  group_name=Grupos_Cestoda_Basal
  get_pho_info_phos

  wanted_species=$Grupos_Taenidae
  group_name=Grupos_Taenidae
  get_pho_info_phos

  wanted_species=$Grupos_Hymenolepis
  group_name=Grupos_Hymenolepis
  get_pho_info_phos

  wanted_species=$Grupos_Echinococus
  group_name=Grupos_Echinococus
  get_pho_info_phos

  wanted_species=$Grupos_Taenia
  group_name=Grupos_Taenia
  get_pho_info_phos

  ########################################################################

  wanted_species=$Grupos_Trematoda
  group_name=Grupos_Trematoda
  get_pho_info_phos

  wanted_species=$Grupos_Schistosmidae
  group_name=Grupos_Schistosmidae
  get_pho_info_phos

  wanted_species=$Grupos_Trematoda_No_schistosoma
  group_name=Grupos_Trematoda_No_schistosoma
  get_pho_info_phos

  wanted_species=$Grupos_Fasciolidae
  group_name=Grupos_Fasciolidae
  get_pho_info_phos

  wanted_species=$Grupos_Trematoda_NoSchistoFasciola
  group_name=Grupos_Trematoda_NoSchistoFasciola
  get_pho_info_phos

  wanted_species=$Grupos_Trematoda_Clonorchis_Opisthorchis
  group_name=Grupos_Trematoda_Clonorchis_Opisthorchis
  get_pho_info_phos

  wanted_species=$Grupos_Paragonimus
  group_name=Grupos_Paragonimus
  get_pho_info_phos

  wanted_species=$Grupos_Schistosoma
  group_name=Grupos_Schistosoma
  get_pho_info_phos

  wanted_species=$Grupos_Fasciola
  group_name=Grupos_Fasciola
  get_pho_info_phos

  Especies=$(ls $work_directory"/"$Result_Storage"/"*pep | sed "s/.*\///" | sed "s/.pep//")
  for spe in $Especies
  do
    wanted_species=$spe
    group_name=$(echo "Individual_"$spe)
    get_pho_info_phos
  done

  ########################################################################

  grupos_interes1=$(echo "Grupos_Cestoda Grupos_Cestoda_Basal Grupos_Cestoda_Derivados Grupos_Taenidae Grupos_Taenia Grupos_Echinococus Grupos_Hymenolepis Grupos_Trematoda Grupos_Trematoda_No_schistosoma Grupos_Trematoda_NoSchistoFasciola Grupos_Trematoda_Clonorchis_Opisthorchis Grupos_Fasciolidae Grupos_Fasciola Grupos_Paragonimus Grupos_Schistosmidae Grupos_Schistosoma")
  grupos_interes2=$(echo "Grupos_Cestoda Grupos_Cestoda_Basal Grupos_Cestoda_Derivados Grupos_Taenidae Grupos_Taenia Grupos_Echinococus Grupos_Hymenolepis Grupos_Trematoda Grupos_Trematoda_No_schistosoma Grupos_Trematoda_NoSchistoFasciola Grupos_Trematoda_Clonorchis_Opisthorchis Grupos_Fasciolidae Grupos_Fasciola Grupos_Paragonimus Grupos_Schistosmidae Grupos_Schistosoma")
  storage=per_group
  comparacion_phos_entre_grupos

  grupos_interes1=$(ls -d $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Group_Specific_PHO/"*"/" | sed "s/\/$//" | sed "s/.*\///" | grep "Individual_" | grep -v "Pairwise_marches")
  grupos_interes2=$(echo "Grupos_Cestoda Grupos_Cestoda_Basal Grupos_Cestoda_Derivados Grupos_Taenidae Grupos_Taenia Grupos_Echinococus Grupos_Hymenolepis Grupos_Trematoda Grupos_Trematoda_No_schistosoma Grupos_Trematoda_NoSchistoFasciola Grupos_Trematoda_Clonorchis_Opisthorchis Grupos_Fasciolidae Grupos_Fasciola Grupos_Paragonimus Grupos_Schistosmidae Grupos_Schistosoma")
  storage=individual_species
  comparacion_phos_entre_grupos

  ##################################################################################################################################################################

  rm -r $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps"
  mkdir $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps"
  mkdir $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo"
  mkdir $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Pairwise"

  species_order=$(echo "Mesocestoides_corti Hmicrostoma Hdiminuta Taenia_multiceps Taenia_saginata Taenia_asiatica Taenia_solium Egranulosus Emultilocularis Spirometra_erinaceieuropaei Schistocephalus_solidus Sparganum_proliferum Trichobilharzia_regenti Schistosoma_japonicum Schistosoma_mansoni Schistosoma_haematobium Schistosoma_bovis Paragonimus_heterotremus Paragonimus_westermani Clonorchis_sinensis Opisthorchis_felineus Fasciolopsis_buski Fasciola_gigantica Fasciola_hepatica")
  pho_unit_of_interest=$(awk -F "\t" '{print $2}' $work_directory"/"$Result_Storage"/Operons/Questions/Work_gene_pairs.txt" | sort -u)

  ###############################################
  ########### Preparar heatmap basico ###########
  ###############################################

  header_operon=Especie
  header_pair=Especie

  individual_hog_file=operon_pho.txt
  pair_hog_file=operon_pair.txt

  for spe in $species_order
  do
    echo $spe" checking_1"
    grep Operon $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | grep -v Read_Coverage |  awk -F "\t" '{print $11}' | sed "s/__/\n/g" | sort -u | grep -v No_Hits >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe"_"$individual_hog_file

    for pair in $pho_unit_of_interest
    do
      check_pair=$(grep -c $pair $work_directory"/"$Result_Storage"/Operons/"$spe"_operon.txt")
      if [ $check_pair -gt 0 ]
      then
        echo $pair >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe"_"$pair_hog_file
      fi
    done
    N_base_pho_Operon=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe"_"$individual_hog_file)
    N_base_pho_Operon_pair=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe"_"$pair_hog_file)

    header_operon=$(echo $header_operon" "$spe"("$N_base_pho_Operon")")
    header_pair=$(echo $header_pair" "$spe"("$N_base_pho_Operon_pair")")
  done

  operon_heatmap

  #################################################
  ########### Preparar heatmap solo SLs ###########
  #################################################

  header_operon=Especie
  header_pair=Especie

  individual_hog_file=operon_pho_sl_only.txt
  pair_hog_file=operon_pair_sl_only.txt

  for spe in $species_order
  do
    echo $spe" checking_2"
    grep -F "*" $work_directory"/"$Result_Storage"/Operons/"$spe"_operon.txt" | awk -F "\t" '{print $7}' | sed "s/__/\n/g" | sort -u | grep -v No_Hits >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe"_"$individual_hog_file
    for pair in $pho_unit_of_interest
    do
      check_pair=$(grep -F "*" $work_directory"/"$Result_Storage"/Operons/"$spe"_operon.txt" | grep -c $pair)
      if [ $check_pair -gt 0 ]
      then
        echo $pair >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe"_"$pair_hog_file
      fi
    done
    N_base_pho_Operon=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe"_"$individual_hog_file)
    N_base_pho_Operon_pair=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe"_"$pair_hog_file)

    header_operon=$(echo $header_operon" "$spe"("$N_base_pho_Operon")")
    header_pair=$(echo $header_pair" "$spe"("$N_base_pho_Operon_pair")")
  done
  operon_heatmap
fi

if [ "$RUN_OPERON_CONSERVATION_QUESTIONS_6" == "TRUE" ]
then
  rm -r $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Heatmaps_Asimetricos"
  mkdir $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Heatmaps_Asimetricos"

  species_order=$(echo "Mesocestoides_corti Hmicrostoma Hdiminuta Taenia_multiceps Taenia_saginata Taenia_asiatica Taenia_solium Egranulosus Emultilocularis Spirometra_erinaceieuropaei Schistocephalus_solidus Sparganum_proliferum Trichobilharzia_regenti Schistosoma_japonicum Schistosoma_mansoni Schistosoma_haematobium Schistosoma_bovis Paragonimus_heterotremus Paragonimus_westermani Clonorchis_sinensis Opisthorchis_felineus Fasciolopsis_buski Fasciola_gigantica Fasciola_hepatica")

  header_operon=Especie
  header_pair=Especie

  individual_hog_file=operon_pho.txt
  pair_hog_file=operon_pair.txt

  for spe in $species_order
  do
    N_base_pho_Operon=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe"_"$individual_hog_file)
    N_base_pho_Operon_pair=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe"_"$pair_hog_file)

    header_operon=$(echo $header_operon" "$spe"("$N_base_pho_Operon")")
    header_pair=$(echo $header_pair" "$spe"("$N_base_pho_Operon_pair")")
  done

  operon_heatmap_asimetrico

  ###############################################################################################################

  header_operon=Especie
  header_pair=Especie

  individual_hog_file=operon_pho_sl_only.txt
  pair_hog_file=operon_pair_sl_only.txt

  for spe in $species_order
  do
    N_base_pho_Operon=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe"_"$individual_hog_file)
    N_base_pho_Operon_pair=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe"_"$pair_hog_file)

    header_operon=$(echo $header_operon" "$spe"("$N_base_pho_Operon")")
    header_pair=$(echo $header_pair" "$spe"("$N_base_pho_Operon_pair")")
  done

  operon_heatmap_asimetrico
fi

#######################################################################################################################################################################################################
#######################################################################################################################################################################################################
#######################################################################################################################################################################################################

if [ "$RUN_OPERON_CONSERVATION_QUESTIONS_5" == "TRUE" ]
then
  rm -r $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Operon_overall_Sumary"
  mkdir $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Operon_overall_Sumary"

  Especies=$(ls $work_directory"/"$Result_Storage"/"*pep | sed "s/.*\///" | sed "s/.pep//")

  individual_pairs=$(awk -F "\t" '{print $2}' $work_directory"/"$Result_Storage"/Operons/Questions/Work_gene_pairs.txt")

  echo "Species N_Candidates N_Candidates_with_SL N_Uniq_HOGs N_Uniq_HOGs_OpSL N_Chimeric N_HOG_Pair N_HOG_Pair_SL" | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Operon_overall_Sumary/Resumen_Operones.txt"
  for spe in $Especies
  do
    N_Operon_Candidates=$(awk -F "\t" '{print $1}' $work_directory"/"$Result_Storage"/Operons/"$spe"_operon.txt" | grep -v -c OpID)
    N_Operon_Candidates_SL=$(grep -F "*" $work_directory"/"$Result_Storage"/Operons/"$spe"_operon.txt" | awk -F "\t" '{print $1}' | grep -c .)

    N_Uniq_HOG=$(awk -F "\t" '{print $7}' $work_directory"/"$Result_Storage"/Operons/"$spe"_operon.txt" | grep -v -w PHOs | sed "s/__/\n/g" | grep -v No_Hits | sort -u | grep -c .)
    N_Uniq_HOG_in_op_SL=$(grep -F "*" $work_directory"/"$Result_Storage"/Operons/"$spe"_operon.txt" | awk -F "\t" '{print $7}' | grep -v -w PHOs | sed "s/__/\n/g" | grep -v No_Hits | sort -u | grep -c .)

    N_Chimeric=$(grep -F "[C]" $work_directory"/"$Result_Storage"/Operons/"$spe"_operon.txt" | awk -F "\t" '{print $1}' | sort -u | grep -c .)

    total_pairs_all=0
    total_pairs_sl=0

    for pair in $individual_pairs
    do
      echo "Run: "$spe" "$pair
      check_pair=$(grep -c $pair $work_directory"/"$Result_Storage"/Operons/"$spe"_operon.txt")
      if [ $check_pair -gt 0 ]
      then
        total_pairs_all=$(($total_pairs_all + 1))
      fi

      check_pair_sl=$(grep -F "*" $work_directory"/"$Result_Storage"/Operons/"$spe"_operon.txt" | grep -c $pair)
      if [ $check_pair_sl -gt 0 ]
      then
        total_pairs_sl=$(($total_pairs_sl + 1))
      fi
    done
    echo $spe" "$N_Operon_Candidates" "$N_Operon_Candidates_SL" "$N_Uniq_HOG" "$N_Uniq_HOG_in_op_SL" "$N_Chimeric" "$total_pairs_all" "$total_pairs_sl | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/Operon_overall_Sumary/Resumen_Operones.txt"
  done
fi
#######################################################################################################################################################################################################
#######################################################################################################################################################################################################
#######################################################################################################################################################################################################

if [ "$RUN_OPERON_CONSERVATION_QUESTIONS_7" == "TRUE" ]
then
  results_orthofinder=$(ls -d $work_directory"/"$Result_Storage"/OrthoFinder/Results_"*)
  PHO_File=$(echo $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")

  rm -r $work_directory"/"$Result_Storage"/Operons/Questions/Specific/HOG_vs_OG"
  mkdir $work_directory"/"$Result_Storage"/Operons/Questions/Specific/HOG_vs_OG"

  HOG_interest=$(awk -F "\t" '{print $2}' $work_directory"/"$Result_Storage"/Operons/Questions/Work_gene_pairs.txt" | sed "s/__/\n/" | sort -u)
  echo "HOG N_Original_Pairs N_Related_HOGs Original_OG Related_HOGs N_Related_HOG_Pairs" | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/HOG_vs_OG/Original_OG_related_HOG_pairs.txt"

  for hog in $HOG_interest
  do
    N_original_pairs=$(grep $hog $work_directory"/"$Result_Storage"/Operons/Questions/Work_gene_pairs.txt" | awk -F "\t" '{print $2}' | grep -c .)

    original_og=$(grep -w $hog $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv" | awk -F "\t" '{print $2}' | sort -u)
    check_original_og=$(awk -F "\t" -v ori=$original_og '{if ($2==ori) print $1}' $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv" | grep -c . | awk '{print $1-1}')

    echo $hog" "$original_og" "$check_original_og
    if [ $check_original_og -eq 0 ]
    then
      status=$(echo ".")
      print_related_hogs=NA
    else
      status=$(echo "X")

      related_hogs=$(awk -F "\t" -v ori=$original_og '{if ($2==ori) print $1}' $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv" | grep -v $hog)
      total_related_hog_pairs=0
      registro_related_hog_pairs=NA
      for related in $related_hogs
      do
        check_related_hog_pairs=$(grep -c $related $work_directory"/"$Result_Storage"/Operons/Questions/Work_gene_pairs.txt")
        if [ $check_related_hog_pairs -gt 0 ]
        then
          total_related_hog_pairs=$(($total_related_hog_pairs + $check_related_hog_pairs))
        fi
      done
      N_related_hog_pairs=$(echo $registro_related_hog_pairs | tr ";" "\n" | grep -c -v NA)
      print_related_hogs=$(echo $related_hogs | tr " " ";" | sed "s/;$//")
    fi

    echo $hog" "$N_original_pairs" "$check_original_og" "$original_og" "$print_related_hogs | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/HOG_vs_OG/Original_OG_related_HOG_pairs.txt"
  done

  Interest_HOG_Pairs=$(awk -F "\t" '{print $2}' $work_directory"/"$Result_Storage"/Operons/Questions/Work_gene_pairs.txt" | sort -u)
  for pair in $Interest_HOG_Pairs
  do
    echo "check: "$pair
    first_HOG=$(echo $pair | awk -F "__" '{print $1}')
    second_HOG=$(echo $pair | awk -F "__" '{print $2}')

    presence_species=$(grep -w $pair $work_directory"/"$Result_Storage"/Operons/Questions/Presence_or_Absence_By_Pairs.txt" | awk -F "\t" '{print $2}')

    check_first=$(awk -F "\t" -v hog=$first_HOG '{if ($1==hog) print $3}' $work_directory"/"$Result_Storage"/Operons/Questions/Specific/HOG_vs_OG/Original_OG_related_HOG_pairs.txt" | awk '{if ($1==0) print "X"; else print "."}')
    check_second=$(awk -F "\t" -v hog=$second_HOG '{if ($1==hog) print $3}' $work_directory"/"$Result_Storage"/Operons/Questions/Specific/HOG_vs_OG/Original_OG_related_HOG_pairs.txt" | awk '{if ($1==0) print "X"; else print "."}')

    if [ $check_first == "X" ] && [ $check_second == "X" ]
    then
      echo $pair" "$presence_species | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/HOG_vs_OG/Specific_pairs.txt"
    else
      echo $pair" "$presence_species | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/HOG_vs_OG/Compromised_pairs.txt"
    fi
  done

  # Making Summary of the problem
  # Specific_HOGs
  Specific_HOGs=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Questions/Specific/HOG_vs_OG/Specific_pairs.txt")
  All_HOGs=$(grep -v N_Original_Pairs $work_directory"/"$Result_Storage"/Operons/Questions/Specific/HOG_vs_OG/Original_OG_related_HOG_pairs.txt" | grep -c .)
  specific_hogs_percentage=$(echo $Specific_HOGs" "$All_HOGs | awk -F " " '{printf "%.2f\n", $1/$2*100}' | awk '{print $1"%"}')

  # HOGs with compromised pairs
  Compromised_HOGs=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Questions/Specific/HOG_vs_OG/Compromised_pairs.txt")
  compromised_hogs_percentage=$(echo $Compromised_HOGs" "$All_HOGs | awk -F " " '{printf "%.2f\n", $1/$2*100}' | awk '{print $1"%"}')

  tier_5_specific_pairs=$(awk -F "\t" '{if ($2>=5) print}' $work_directory"/"$Result_Storage"/Operons/Questions/Specific/HOG_vs_OG/Specific_pairs.txt" | grep -c .)
  mean_specific_pairs_spe=$(awk -F "\t" '{print $2}' $work_directory"/"$Result_Storage"/Operons/Questions/Specific/HOG_vs_OG/Specific_pairs.txt" | datamash mean 1 )

  tier_5_compromised_pairs=$(awk -F "\t" '{if ($2>=5) print}' $work_directory"/"$Result_Storage"/Operons/Questions/Specific/HOG_vs_OG/Compromised_pairs.txt" | grep -c .)
  mean_compromised_pairs_spe=$(awk -F "\t" '{print $2}' $work_directory"/"$Result_Storage"/Operons/Questions/Specific/HOG_vs_OG/Compromised_pairs.txt"  | datamash mean 1 )

  echo "Summary of the HOG vs OG dispute" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/HOG_vs_OG/Summary.txt"
  echo "Specific HOGs: "$Specific_HOGs"/"$All_HOGs" ("$specific_hogs_percentage")" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/HOG_vs_OG/Summary.txt"
  echo "Compromised HOGs with Pairs: "$Compromised_HOGs"/"$All_HOGs" ("$compromised_hogs_percentage")" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/HOG_vs_OG/Summary.txt"
  echo "............................................." >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/HOG_vs_OG/Summary.txt"
  echo "" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/HOG_vs_OG/Summary.txt"
  echo "Total specific pairs: "$total_specific_pairs >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/HOG_vs_OG/Summary.txt"
  echo "Mean species on specific pairs: "$mean_specific_pairs_spe >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/HOG_vs_OG/Summary.txt"
  echo "More than 5: "$tier_5_specific_pairs >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/HOG_vs_OG/Summary.txt"
  echo "" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/HOG_vs_OG/Summary.txt"
  echo "Total compromised pairs: "$total_compromised_pairs >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/HOG_vs_OG/Summary.txt"
  echo "Mean species on compromised pairs: "$mean_compromised_pairs_spe >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/HOG_vs_OG/Summary.txt"
  echo "More than 5: "$tier_5_compromised_pairs >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/HOG_vs_OG/Summary.txt"
fi



if [ "$RUN_OPERON_CONSERVATION_QUESTIONS_8" == "TRUE" ]
then
  results_orthofinder=$(ls -d $work_directory"/"$Result_Storage"/OrthoFinder/Results_"*)
  PHO_File=$(echo $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")

  species_order=$(echo "Mesocestoides_corti Hmicrostoma Hdiminuta Taenia_multiceps Taenia_saginata Taenia_asiatica Taenia_solium Egranulosus Emultilocularis Spirometra_erinaceieuropaei Schistocephalus_solidus Sparganum_proliferum Trichobilharzia_regenti Schistosoma_japonicum Schistosoma_mansoni Schistosoma_haematobium Schistosoma_bovis Paragonimus_heterotremus Paragonimus_westermani Clonorchis_sinensis Opisthorchis_felineus Fasciolopsis_buski Fasciola_gigantica Fasciola_hepatica")
  pho_unit_of_interest=$(awk -F "\t" '{print $2}' $work_directory"/"$Result_Storage"/Operons/Questions/Work_gene_pairs.txt" | sort -u)

  header_operon=Especie
  header_pair=Especie

  individual_hog_file=operon_pair_sl_second_in_pair_only_hogs.txt
  pair_hog_file=operon_pair_sl_second_in_pair_only_pairs.txt

  rm -r $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Temp_Files_better_heatmap"
  mkdir $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Temp_Files_better_heatmap"

  for spe in $species_order
  do
    rm $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe"_"$individual_hog_file
    rm $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe"_"$pair_hog_file

    operon_ids=$(grep -F "*" $work_directory"/"$Result_Storage"/Operons/"$spe"_operon.txt" | awk -F "\t" '{print $1}')
    for opid in $operon_ids
    do
      echo $spe": "$opid
      grep -w $opid $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" > $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Temp_Files_better_heatmap/temp_gene.tmp"
      count=1
      recorrer=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Temp_Files_better_heatmap/temp_gene.tmp" )

      while [ $count -le $recorrer ]
      do
        check_chimera=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Temp_Files_better_heatmap/temp_gene.tmp"  | awk -F "\t" '{print $7}')
        if [ $check_chimera == "Chimeric" ]
        then
          operon_id=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Temp_Files_better_heatmap/temp_gene.tmp"  | awk -F "\t" '{print $10}')
          gen_id=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Temp_Files_better_heatmap/temp_gene.tmp"  | awk -F "\t" '{print $1}')
          grep -w -F $gen_id $work_directory"/SL_Read_Counts/"$spe"_"*"_read.counts" | awk -F "\t" '{print $2"_Chimera_"}' | sort -u > $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Temp_Files_better_heatmap/temp_transcript_chimeric.tmp"

          chimeric_parts=$(grep -f $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Temp_Files_better_heatmap/temp_transcript_chimeric.tmp" $work_directory"/"$Result_Storage"/"$spe".pep" | sed "s/>//")
          for part in $chimeric_parts
          do
            check_chim_hog=$(grep -c -w -F $part $PHO_File)
            if [ $check_chim_hog -gt 0 ]
            then
              hog=$(grep -w -F $part $PHO_File | awk -F "\t" '{print $1}')
            else
              hog=No_Hits
            fi

            check_sl=$(grep -w -F -c $part $work_directory"/"$Result_Storage"/Mediciones_Expresion/Background_TPM/"$spe"_per_gene.txt")
            if [ $check_sl -gt 0 ]
            then
              sl=SL
            else
              sl=X
            fi
            echo $operon_id" "$part" "$sl" "$hog | tr " " "\t" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Temp_Files_better_heatmap/"$spe"_genes_in_operon.txt"
          done
        else
           sed -n $count"p" $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Temp_Files_better_heatmap/temp_gene.tmp" | awk -F "\t" '{print $10"\t"$1"\t"$6"\t"$11}' >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Temp_Files_better_heatmap/"$spe"_genes_in_operon.txt"
        fi
        count=$(($count + 1))
      done
    done

    for op_pair in $operon_ids
    do
      grep -w -F $op_pair $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Temp_Files_better_heatmap/"$spe"_genes_in_operon.txt" > $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Temp_Files_better_heatmap/temp.tmp"
      count=2
      recorrer=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Temp_Files_better_heatmap/temp.tmp")

      while [ $count -le $recorrer ]
      do
        previous=$(($count - 1))
        skip_no_hit=$(grep -n . $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Temp_Files_better_heatmap/temp.tmp" | sed -n $previous","$count"p" | grep -c No_Hits)
        check_sl=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Temp_Files_better_heatmap/temp.tmp" | awk -F "\t" '{print $3}')

        if [ $skip_no_hit -eq 0 ]
        then
          if [ $check_sl == "SL" ]
          then
             sed -n $previous","$count"p" $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Temp_Files_better_heatmap/temp.tmp" | awk -F "\t" '{print $4}' | tr "\n" "_" | sed "s/_/__/" | sed "s/_$/\n/" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Temp_Files_better_heatmap/"$spe"_work_pairs.tmp"
          fi
        fi
        count=$(($count + 1))
      done
    done
    sort -u $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Temp_Files_better_heatmap/"$spe"_work_pairs.tmp" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe"_"$pair_hog_file
    sort -u $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Temp_Files_better_heatmap/"$spe"_work_pairs.tmp" >> $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe"_"$individual_hog_file

    N_base_pho_Operon=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe"_"$individual_hog_file)
    N_base_pho_Operon_pair=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe"_"$pair_hog_file)

    header_operon=$(echo $header_operon" "$spe"("$N_base_pho_Operon")")
    header_pair=$(echo $header_pair" "$spe"("$N_base_pho_Operon_pair")")
  done

  operon_heatmap

  header_operon=Especie
  header_pair=Especie

  for spe in $species_order
  do
    N_base_pho_Operon=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe"_"$individual_hog_file)
    N_base_pho_Operon_pair=$(grep -c . $work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo/"$spe"_"$pair_hog_file)

    header_operon=$(echo $header_operon" "$spe"("$N_base_pho_Operon")")
    header_pair=$(echo $header_pair" "$spe"("$N_base_pho_Operon_pair")")
  done

  operon_heatmap_asimetrico
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$RUN_HOG_ITOL_DATASETS" == "TRUE" ]
then
  storage_place=PHO_Selection_TPM

  rm -r $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets"
  mkdir $work_directory"/"$Result_Storage"/"$storage_place"/Other_Data/Itol_Datasets"

  selected_PHOs=$(sed 1d $work_directory"/"$Result_Storage"/"$storage_place"/General_Selection_PHO.txt" | awk -F "\t" '{print $1}' | sort -u)
  interest_pho_itol
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$RUN_ACCEPTOR_TAG_COUNTS" == "TRUE" ]
then
  rm -r $work_directory"/"$Result_Storage"/Accepor_Site_Analisis"
  mkdir $work_directory"/"$Result_Storage"/Accepor_Site_Analisis"
  mkdir $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Count_Table_Per_Tag"

  Especies=$(ls $work_directory"/"$Result_Storage"/"*pep | sed "s/.*\///" | sed "s/.pep//")

  for spe in $Especies
  do
    tags_sl=$(ls $work_directory"/SL_Read_Counts/"$spe"_"*"_read.counts" | sed "s/_read.counts//" | sed "s/.*\/$spe//" | sed "s/^_//")
    grupo=$(ls $work_directory"/"$Result_Storage"/Mediciones_Expresion/Mapeos/"*"_"$spe"_"*"_Aligned.sortedByCoord.out.bam" | sed "s/.*\///" | awk -F "_" '{print $1}' | sort -u)
    samples=$(ls $work_directory"/"$Result_Storage"/Mediciones_Expresion/Mapeos/"*"_"$spe"_"*"_Aligned.sortedByCoord.out.bam" | sed "s/.*\///" | sed "s/.*$spe//" | awk -F "_" '{print $2}')
    total_samples=$(echo $samples | tr " " "\n" | grep -c .)

    all_sl_data=$(grep -w Concordante $work_directory"/SL_Read_Counts/"$spe"_tier3.counts" | grep -w -v Trans_End | awk -F "\t" '{print $2";"$3";"$4}' | sort -u)
    is_operon=$(echo .)

    echo "Ace Chromosome Strand Gen_IDs Trans_IDs Is_Operon N_Samples Total_SL_Reads "$tags_sl | tr " " "\t" >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Count_Table_Per_Tag/"$spe"_SL_Tags_Counts.tab"

    for sl_data in $all_sl_data
    do
      basic_info_ace_counts

      if [ "$check_double_gene_nonsense" == "Single" ]
      then
        total_sl_reads=0
        registro=sample

        for tag in $tags_sl
        do
          echo "TAGS: "$spe" || "$sl_data" || "$tag
          check_tag_counts=$(grep -w Concordante $work_directory"/SL_Read_Counts/"$spe"_"$tag"_read.counts" | awk -F "\t" -v ace=$ace -v chr=$chr '{if ($5==ace && $3==chr) print}' | grep -w -v -c Trans_End)

          if [ $check_tag_counts -gt 0 ]
          then
            tag_counts=$(grep -w Concordante $work_directory"/SL_Read_Counts/"$spe"_"$tag"_read.counts" | awk -F "\t" -v ace=$ace -v chr=$chr '{if ($5==ace && $3==chr) print}' | grep -w -v Trans_End | tr "\t" "\n" | tail -n 1)
          else
            tag_counts=0
          fi
          total_sl_reads=$(($total_sl_reads + $tag_counts))
          registro=$(echo $registro" "$tag_counts)
        done
        registro=$(echo $registro | tr " " "\n" | grep -v sample)

        print_genID=$(echo $gen_id | tr "\n" ";" | sed "s/;$//")

        echo $ace" "$chr" "$strand" "$print_genID" "$trans_ids" "$is_operon" "$total_samples" "$total_sl_reads" "$registro | tr " " "\t" >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Count_Table_Per_Tag/"$spe"_SL_Tags_Counts.tab"
      fi
    done

    total_reads=$(grep -v -w "N_Samples" $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Count_Table_Per_Tag/"$spe"_SL_Tags_Counts.tab" | awk -F "\t" '{print $8}' | awk -F "\t" '{ sum += $1 } END {print sum*0.05}')
    awk -F "\t" '{print $1"__"$2"__"$3}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Count_Table_Per_Tag/"$spe"_SL_Tags_Counts.tab" >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Count_Table_Per_Tag/"$spe"_table_make.temp"

    for tag in $tags_sl
    do
      column_num=$(head -n 1 $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Count_Table_Per_Tag/"$spe"_SL_Tags_Counts.tab" | tr "\t" "\n" | grep -n $tag)
      check=$(grep -v -w "N_Samples" $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Count_Table_Per_Tag/"$spe"_SL_Tags_Counts.tab" | awk -F "\t" -v N=$column_num '{print $N}' | awk '{ sum += $1 } END {print sum}' | awk -v cutooff=$total_reads '{if ($1>=cutooff) print "TRUE"; else print "FALSE"}')

      if [ $check == "TRUE" ]
      then
        awk -F "\t" -v N=$column_num '{print $N}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Count_Table_Per_Tag/"$spe"_SL_Tags_Counts.tab" >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Count_Table_Per_Tag/"$spe"_"$tag"_tag.temp"
      fi
    done

    paste $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Count_Table_Per_Tag/"$spe"_table_make.temp" $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Count_Table_Per_Tag/"$spe"_"*"_tag.temp" >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Count_Table_Per_Tag/"$spe"_SL_Work_Counts.tab"
    rm $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Count_Table_Per_Tag/"$spe"_"*"_tag.temp"
    rm $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Count_Table_Per_Tag/"$spe"_table_make.temp"
  done
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$RUN_ACCEPTOR_COMPETENCE_INITIAL" == "TRUE" ]
then
  rm -r $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition"
  mkdir $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition"

  Especies=$(ls $work_directory"/"$Result_Storage"/"*pep | sed "s/.*\///" | sed "s/.pep//")

  storage=$work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition"
  sj_info=$work_directory"/"$Result_Storage"/Mediciones_Expresion/Mapeos"
  is_this_first=TRUE
  competicion_summary
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$RUN_ACCEPTOR_COMPETENCE_SECOND" == "TRUE" ]
then
  rm -r $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass"
  mkdir $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass"

  Especies=$(ls $work_directory"/"$Result_Storage"/"*pep | sed "s/.*\///" | sed "s/.pep//")
  storage=$work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass"
  sj_info=$work_directory"/"$Result_Storage"/Mapeos_Second_Pass"
  is_this_first=FALSE
  competicion_summary
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$RUN_ACCEPTOR_COMPETENCE_OPERON" == "TRUE" ]
then
  rm -r $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron"
  mkdir $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron"
  mkdir $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/Resgistros_por_sitio"

  storage=$work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition"
  sj_info=$work_directory"/"$Result_Storage"/Mapeos_Second_Pass"

  Especies=$(ls $work_directory"/"$Result_Storage"/"*pep | sed "s/.*\///" | sed "s/.pep//")
  echo "Especie Grupo Sitios_SJ_Competencia Todos" >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/Resumen_Buenos_sitios_competencia.txt"

  Cestoda=$(echo "Egranulosus Emultilocularis Hdiminuta Hmicrostoma Mesocestoides_corti Sparganum_proliferum Spirometra_erinaceieuropaei Taenia_asiatica Taenia_multiceps Taenia_saginata Taenia_solium Schistocephalus_solidus")

  for spe in $Especies
  do
    check_cestodes=$(echo $Cestoda | grep -c $spe)
    if [ $check_cestodes -gt 0 ]
    then
      grupo=Cestodes
    else
      grupo=Trematoda
    fi

    sl_tags=$(head -n 1 $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Count_Table_Per_Tag/"$spe"_SL_Work_Counts.tab" | tr "\t" "\n" | sed 1d)
    gff_file=$(ls $grupo"/"$spe"/"*"annotations.gtf")

    echo "Ace_Set Ace Chr N_Reads N_Donadores "$sl_tags" Todos" | tr " " "\t" >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt"

    count=2
    recorrer=$(grep -c . $work_directory"/"$Result_Storage"/Operons/"$spe"_operon.txt")
    while [ $count -le $recorrer ]
    do
      genes_operon=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Operons/"$spe"_operon.txt" | awk -F "\t" '{print $6}' | sed "s/_(/\n/g" | sed "s/.*)_//")
      pos_operon=1

      for gen in $genes_operon
      do
        echo "Recuperando aces: "$gen" --- "$pos_operon

        check_chimera=$(echo $gen | grep -c -F "[C]")
        if [ $check_chimera -gt 0 ]
        then
          if [ $spe == "Fasciolopsis_buski" ] || [ $spe == "Paragonimus_heterotremus" ]
          then
            check_gene=$(echo $gen | awk -F "[" '{print $1}' | awk '{print "gene-"$0}')
          else
            check_gene=$(echo $gen | awk -F "[" '{print $1}')
          fi

          grep -w -F $check_gene $work_directory"/SL_Read_Counts/"$spe"_tier3.counts" | grep -w Concordante | grep -w -v Trans_End | awk -F "\t" '{print $2"__"$3}' | sort -u >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/Resgistros_por_sitio/"$spe"_chimera_ace.txt"
          pos_operon=$(($pos_operon + 1))
        else
          if [ ! $pos_operon -eq 1 ]
          then
            if [ $spe == "Fasciolopsis_buski" ] || [ $spe == "Paragonimus_heterotremus" ]
            then
              check_gene=$(echo $gen | sed "s/\*//" | awk '{print "gene-"$0}')
            else
              check_gene=$(echo $gen | sed "s/\*//")
            fi
            strand=$(grep $gen'";' $gff_file | awk -F "\t" '{print $7}' | sort -u)
            if [ $strand == "+" ]
            then
              start_codon=$(grep $gen'";' $gff_file| awk -F "\t" '{if ($3=="CDS") print $4}' | sort -n  | head -n 1)
              grep -w -F $check_gene $work_directory"/SL_Read_Counts/"$spe"_tier3.counts" | grep -w Concordante | grep -w -v Trans_End | awk -F "\t" -v CDS=$start_codon '{if ($2 < CDS) print $2"__"$3}' | sort -u >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/Resgistros_por_sitio/"$spe"_operon_ace.txt"
            else
              start_codon=$(grep $gen'";' $gff_file| awk -F "\t" '{if ($3=="CDS") print $5}' | sort -n  | tail -n 1)
              grep -w -F $check_gene $work_directory"/SL_Read_Counts/"$spe"_tier3.counts" | grep -w Concordante | grep -w -v Trans_End | awk -F "\t" -v CDS=$start_codon '{if ($2 > CDS) print $2"__"$3}' | sort -u >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/Resgistros_por_sitio/"$spe"_operon_ace.txt"
            fi
          fi
        fi

        pos_operon=$(($pos_operon + 1))
      done
      count=$(($count + 1))
    done

    echo "Done... Time for all others"
    grep -w Concordante $work_directory"/SL_Read_Counts/"$spe"_tier3.counts" | grep -w -v Trans_End | awk -F "\t" '{print $2"__"$3}' | sort -u > $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/Resgistros_por_sitio/Temp.txt"
    grep -v -F -w -f $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/Resgistros_por_sitio/"$spe"_operon_ace.txt" $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/Resgistros_por_sitio/Temp.txt" | grep -v -F -w -f $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/Resgistros_por_sitio/"$spe"_chimera_ace.txt" | sort -u > $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/Resgistros_por_sitio/"$spe"_other_ace.txt"

    # Agregar sitios al inicio del gen, estricto
    count_estricto=1
    recorrer_estricto=$(grep -c . $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/Resgistros_por_sitio/"$spe"_other_ace.txt")

    while [ $count_estricto -le $recorrer_estricto ]
    do
      ace=$(sed -n $count_estricto"p" $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/Resgistros_por_sitio/"$spe"_other_ace.txt" | awk -F "__" '{print $1}')
      chr=$(sed -n $count_estricto"p" $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/Resgistros_por_sitio/"$spe"_other_ace.txt" | awk -F "__" '{print $2}')

      check_single=$(grep -w $ace $storage"/"$spe"_Cis_SJ_on_SL_ace.tab" | grep -w -F $chr | grep -w -c Single)

      if [ $check_single -eq 1 ]
      then
        echo "check better: "$ace" -- "$chr
        gen=$(grep -w $ace $storage"/"$spe"_Cis_SJ_on_SL_ace.tab" | grep -w -F $chr | awk -F "\t" '{print $5}' )

        gen_start=$(grep -w -F $gen $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $3}')
        gen_end=$(grep -w -F $gen $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $4}' )
        gen_strand=$(grep -w -F $gen $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $5}')
        grep -w -F $gen $work_directory"/SL_Read_Counts/"$spe"_tier3.counts" | grep -w $ace | grep -w Concordante | grep -w -v Trans_End | awk -F "\t" -v strand=$gen_strand -v start=$gen_start -v end=$gen_end '{if ($4=="+" && start-$2>0) print $2"__"$3; else if ($4=="-" && $2-end >0) print $2"__"$3}' | sort -u >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/Resgistros_por_sitio/"$spe"_strict_start_ace.txt"
      fi
      count_estricto=$(($count_estricto + 1))
    done

    grep -v -F -w -f $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/Resgistros_por_sitio/"$spe"_strict_start_ace.txt" $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/Resgistros_por_sitio/"$spe"_other_ace.txt" | sort -u > $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/Resgistros_por_sitio/"$spe"_other_strict_ace.txt"

    ace_set=operon
    sj_ace_information

    ace_set=other_strict
    sj_ace_information

    ace_set=strict_start
    sj_ace_information

    ace_set=chimera
    sj_ace_information

    grep -v -w chimera $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads_work.txt"
  done
fi

if [ "$RUN_ACCEPTOR_COMPETENCE_FINAL" == "TRUE" ]
then
  rm $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Final_Summary.txt"

  Especies=$(ls $work_directory"/"$Result_Storage"/"*pep | sed "s/.*\///" | sed "s/.pep//")
  echo "Species Group_SL-Aces Total SL_Only Competition Competition_Good N1_Donnor_Site N2_Donnor_Site N3_Donnor_Site More_Donnor_Site" | tr " " "\t" > $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Final_Summary.txt"

  for spe in $Especies
  do
    # Conteos Other
    All_other=$(awk -F "\t" '{if ($1=="other_strict") print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)
    other_SL_only=$(awk -F "\t" '{if ($1=="other_strict" && $4==0) print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)
    other_SL_competition=$(awk -F "\t" '{if ($1=="other_strict" && $4>0) print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)
    other_SL_competition_good=$(awk -F "\t" '{if ($1=="other_strict" && $4>3) print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)
    other_donnor_site_1=$(awk -F "\t" '{if ($1=="other_strict" && $5==1) print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)
    other_donnor_site_2=$(awk -F "\t" '{if ($1=="other_strict" && $5==2) print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)
    other_donnor_site_3=$(awk -F "\t" '{if ($1=="other_strict" && $5==3) print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)
    other_donnor_site_more=$(awk -F "\t" '{if ($1=="other_strict" && $5>3) print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)

    All_operon=$(awk -F "\t" '{if ($1=="operon") print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)
    operon_SL_only=$(awk -F "\t" '{if ($1=="operon" && $4==0) print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)
    operon_SL_competition=$(awk -F "\t" '{if ($1=="operon" && $4>0) print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)
    operon_SL_competition_good=$(awk -F "\t" '{if ($1=="operon" && $4>3) print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)
    operon_donnor_site_1=$(awk -F "\t" '{if ($1=="operon" && $5==1) print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)
    operon_donnor_site_2=$(awk -F "\t" '{if ($1=="operon" && $5==2) print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)
    operon_donnor_site_3=$(awk -F "\t" '{if ($1=="operon" && $5==3) print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)
    operon_donnor_site_more=$(awk -F "\t" '{if ($1=="operon" && $5>3) print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)

    All_start=$(awk -F "\t" '{if ($1=="strict_start") print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)
    start_SL_only=$(awk -F "\t" '{if ($1=="strict_start" && $4==0) print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)
    start_SL_competition=$(awk -F "\t" '{if ($1=="strict_start" && $4>0) print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)
    start_SL_competition_good=$(awk -F "\t" '{if ($1=="strict_start" && $4>3) print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)
    start_donnor_site_1=$(awk -F "\t" '{if ($1=="strict_start" && $5==1) print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)
    start_donnor_site_2=$(awk -F "\t" '{if ($1=="strict_start" && $5==2) print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)
    start_donnor_site_3=$(awk -F "\t" '{if ($1=="strict_start" && $5==3) print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)
    start_donnor_site_more=$(awk -F "\t" '{if ($1=="strict_start" && $5>3) print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)

    All_chimera=$(awk -F "\t" '{if ($1=="chimera") print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)
    chimera_SL_only=$(awk -F "\t" '{if ($1=="chimera" && $4==0) print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)
    chimera_SL_competition=$(awk -F "\t" '{if ($1=="chimera" && $4>0) print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)
    chimera_SL_competition_good=$(awk -F "\t" '{if ($1=="chimera" && $4>3) print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)
    chimera_donnor_site_1=$(awk -F "\t" '{if ($1=="chimera" && $5==1) print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)
    chimera_donnor_site_2=$(awk -F "\t" '{if ($1=="chimera" && $5==2) print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)
    chimera_donnor_site_3=$(awk -F "\t" '{if ($1=="chimera" && $5==3) print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)
    chimera_donnor_site_more=$(awk -F "\t" '{if ($1=="chimera" && $5>3) print}' $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Operon_vs_monocistron/"$spe"_Registro_Splice_reads.txt" | grep -c .)

    echo $spe" Start "$All_start" "$start_SL_only" "$start_SL_competition" "$start_SL_competition_good" "$start_donnor_site_1" "$start_donnor_site_2" "$start_donnor_site_3" "$start_donnor_site_more | tr " " "\t" >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Final_Summary.txt"
    echo $spe" Operon_Resolution "$All_operon" "$operon_SL_only" "$operon_SL_competition" "$operon_SL_competition_good" "$operon_donnor_site_1" "$operon_donnor_site_2" "$operon_donnor_site_3" "$operon_donnor_site_more | tr " " "\t" >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Final_Summary.txt"
    echo $spe" Other "$All_other" "$other_SL_only" "$other_SL_competition" "$other_SL_competition_good" "$other_donnor_site_1" "$other_donnor_site_2" "$other_donnor_site_3" "$other_donnor_site_more | tr " " "\t" >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Final_Summary.txt"
    echo $spe" Chimeric "$All_chimera" "$chimera_SL_only" "$chimera_SL_competition" "$chimera_SL_competition_good" "$chimera_donnor_site_1" "$chimera_donnor_site_2" "$chimera_donnor_site_3" "$chimera_donnor_site_more | tr " " "\t" >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Competition_Second_Pass/Final_Summary.txt"
  done
fi



##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$RUN_ACCEPTOR_HOG_HEATMAP" == "TRUE" ]
then
  rm -r $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Heatmap_SL_HOG"
  mkdir $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Heatmap_SL_HOG"
  mkdir $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Heatmap_SL_HOG/Species_List"

  results_orthofinder=$(ls -d $work_directory"/"$Result_Storage"/OrthoFinder/Results_"*)
  PHO_File=$(echo $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")

  species_order=$(echo "Mesocestoides_corti Hmicrostoma Hdiminuta Taenia_multiceps Taenia_saginata Taenia_asiatica Taenia_solium Egranulosus Emultilocularis Spirometra_erinaceieuropaei Schistocephalus_solidus Sparganum_proliferum Trichobilharzia_regenti Schistosoma_japonicum Schistosoma_mansoni Schistosoma_haematobium Schistosoma_bovis Paragonimus_heterotremus Paragonimus_westermani Clonorchis_sinensis Opisthorchis_felineus Fasciolopsis_buski Fasciola_gigantica Fasciola_hepatica")

  header=Especie
  for spe in $species_order
  do
    header=$(echo $header" "$spe"(replace_"$spe")")
    awk -F "\t" '{if ($6=="SL" && $7=="NA" ) print $11}' $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | sort -u | grep -v No_Hits >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Heatmap_SL_HOG/Species_List/"$spe"_sl_pho.temp"

    chimeric_genes=$(awk -F "\t" '{if ($6=="SL" && $7=="Chimeric" ) print $1}' $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | sort -u)
    for chim_gen in $chimeric_genes
    do
      check_data=$(grep -w -F -c $chim_gen $work_directory"/"$Result_Storage"/Mediciones_Expresion/Background_TPM/"$spe"_per_gene.txt")
      if [ $check_data -eq 0 ]
      then
        echo "ERROR ERROR: "$chim_gen
      else
        get_chim=$(grep -w -F $chim_gen $work_directory"/"$Result_Storage"/Mediciones_Expresion/Background_TPM/"$spe"_per_gene.txt" | awk -F "\t" '{print $1}')
        for part in $get_chim
        do
          hog_chim=$(grep -w -F $part $PHO_File | awk -F "\t" '{print $1}' | sort -u | grep .)
          echo $part" --- "$hog_chim
          echo $hog_chim >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Heatmap_SL_HOG/Species_List/"$spe"_sl_pho.temp"
        done
      fi
    done
    sort -u $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Heatmap_SL_HOG/Species_List/"$spe"_sl_pho.temp" >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Heatmap_SL_HOG/Species_List/"$spe"_sl_pho.ids"
  done

  echo $header | tr " " "\t" >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Heatmap_SL_HOG/Heatmap_SL_HOG_raw.txt"
  echo $header | tr " " "\t" >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Heatmap_SL_HOG/Heatmap_SL_HOG_percentage.txt"
  echo $header | tr " " "\t" >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Heatmap_SL_HOG/Heatmap_SL_HOG_asimetric.txt"

  for spe1 in $species_order
  do
    change=$(echo "replace_"$spe1)
    total_spe1=$(grep -c . $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Heatmap_SL_HOG/Species_List/"$spe1"_sl_pho.ids")

    sed -i "s/$change/$total_spe1/" $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Heatmap_SL_HOG/Heatmap_SL_HOG_raw.txt"
    sed -i "s/$change/$total_spe1/" $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Heatmap_SL_HOG/Heatmap_SL_HOG_percentage.txt"
    sed -i "s/$change/$total_spe1/" $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Heatmap_SL_HOG/Heatmap_SL_HOG_asimetric.txt"

    registro_raw=$(echo $spe1"("$total_spe1")")
    registro_percentage=$(echo $spe1"("$total_spe1")")
    registro_asimetric=$(echo $spe1"("$total_spe1")")

    for spe2 in $species_order
    do
      if [ "$spe1" == "$spe2" ]
      then
        registro_raw=$(echo $registro_raw" NA")
        registro_percentage=$(echo $registro_percentage" NA")
        registro_asimetric=$(echo $registro_asimetric" NA")
      else
        # Raw heatmap
        N_pho_compartidos=$(grep -w -F -c -f $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Heatmap_SL_HOG/Species_List/"$spe1"_sl_pho.ids" $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Heatmap_SL_HOG/Species_List/"$spe2"_sl_pho.ids")
        registro_raw=$(echo $registro_raw" "$N_pho_compartidos)

        # Percetage heatmap symetric
        all_phos=$(cat $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Heatmap_SL_HOG/Species_List/"$spe1"_sl_pho.ids" $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Heatmap_SL_HOG/Species_List/"$spe2"_sl_pho.ids" | sort -u | grep -c .)
        shared_percetage=$(echo $N_pho_compartidos" "$all_phos | awk -F " " '{printf "%.10f\n", $1/$2}')
        registro_percentage=$(echo $registro_percentage" "$shared_percetage)

        # Percetage heatmap asymetric
        spe1_total=$(grep -c . $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Heatmap_SL_HOG/Species_List/"$spe1"_sl_pho.ids")
        asimetric_percentage=$(echo $N_pho_compartidos" "$spe1_total | awk -F " " '{printf "%.10f\n", $1/$2}')
        registro_asimetric=$(echo $registro_asimetric" "$asimetric_percentage)
      fi
    done
    echo $registro_raw | tr " " "\t" >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Heatmap_SL_HOG/Heatmap_SL_HOG_raw.txt"
    echo $registro_percentage | tr " " "\t" >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Heatmap_SL_HOG/Heatmap_SL_HOG_percentage.txt"
    echo $registro_asimetric | tr " " "\t" >> $work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Heatmap_SL_HOG/Heatmap_SL_HOG_asimetric.txt"
  done
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$RUN_BUSCO" == "TRUE" ]
then
  rm -r $work_directory"/"$Result_Storage"/BUSCO_RUNS"
  mkdir $work_directory"/"$Result_Storage"/BUSCO_RUNS"

  Especies=$(ls $work_directory"/"$Result_Storage"/"*pep | sed "s/.*\///" | sed "s/.pep//")
  for spe in $Especies
  do
    busco -i $work_directory"/"$Result_Storage"/"$spe".pep" -o $spe"_Busco" -m prot -l metazoa
    mv $spe"_Busco" $work_directory"/"$Result_Storage"/BUSCO_RUNS/"$spe
  done
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$RUN_BUSCO_ANOTATION" == "TRUE" ]
then
  rm -r $work_directory"/"$Result_Storage"/BUSCO_Genes_SL"
  mkdir $work_directory"/"$Result_Storage"/BUSCO_Genes_SL"
  mkdir $work_directory"/"$Result_Storage"/BUSCO_Genes_SL/Temp"

  results_orthofinder=$(ls -d $work_directory"/"$Result_Storage"/OrthoFinder/Results_"*)
  PHO_File=$(echo $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")

  Especies=$(ls $work_directory"/"$Result_Storage"/"*pep | sed "s/.*\///" | sed "s/.pep//")
  echo "Especie N_Genes_SL N_Genes_Busco N_Genes_Busco_SL N_Genes_Operon N_Genes_Operon_Busco" | tr " " "\t" >> $work_directory"/"$Result_Storage"/BUSCO_Genes_SL/Summary_SL_genes_busco.txt"

  for spe in $Especies
  do
    genes=$(awk -F "\t" -v N=$read_count_column '{if ($5=="Concordante" && $6!="Trans_End") print $1}' $work_directory"/SL_Read_Counts/"$spe"_tier3.counts" | sed "s/^gene-//" | sort -u)
    echo "GenID TransID Chimeric PHO_ID Is_Operon OperonID Busco_OG_ID Status N_Genes_Busco" | tr " " "\t" >> $work_directory"/"$Result_Storage"/BUSCO_Genes_SL/"$spe"_SL_genes_busco.txt"
    for gen in $genes
    do
      echo "Debug "$spe": "$gen
      check_data=$(grep -c -w -F $gen $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt")

      if [ $check_data -gt 0 ]
      then
        chimera=$(grep -w -F $gen $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $7}')
        is_operon=$(grep -w -F $gen $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{if ($10 !="NA") print "is_operon"; else print "NA"}' )
        ID_operon=$(grep -w -F $gen $work_directory"/"$Result_Storage"/Operons/"$spe"_gene.txt" | awk -F "\t" '{print $10}')

        if [ $chimera == "Chimeric" ]
        then
          grep -P $gen"\t" $work_directory"/SL_Read_Counts/"$spe"_"*"_read.counts" | awk -F "\t" '{print $2}' | sort -u | sed "s/^gene-//" | awk -F "\t" '{print $1"_Chimera_"}' > $work_directory"/"$Result_Storage"/BUSCO_Genes_SL/Temp/trans.temp"
          transcripts=$(grep -F -f $work_directory"/"$Result_Storage"/BUSCO_Genes_SL/Temp/trans.temp" $work_directory"/"$Result_Storage"/"$spe".pep" | sed "s/>//")
        else
          grep -P $gen"\t" $work_directory"/SL_Read_Counts/"$spe"_"*"_read.counts" | awk -F "\t" '{print $2}' | sort -u | sed "s/^gene-//" > $work_directory"/"$Result_Storage"/BUSCO_Genes_SL/Temp/trans.temp"
          transcripts=$(grep -w -F -f $work_directory"/"$Result_Storage"/BUSCO_Genes_SL/Temp/trans.temp" $work_directory"/"$Result_Storage"/"$spe".pep" | sed "s/>//")
        fi
      else
        grep -P $gen"\t" $work_directory"/SL_Read_Counts/"$spe"_"*"_read.counts" | awk -F "\t" '{print $2}' | sort -u | sed "s/^gene-//" > $work_directory"/"$Result_Storage"/BUSCO_Genes_SL/Temp/trans.temp"
        transcripts=$(grep -F -f $work_directory"/"$Result_Storage"/BUSCO_Genes_SL/Temp/trans.temp" $work_directory"/"$Result_Storage"/"$spe".pep" | sed "s/>//")
        chimera=Skipped_Overlap
        is_operon=Skipped_Overlap
        ID_operon=Skipped_Overlap
      fi

      check_ortho_data=$(echo $transcripts | grep -c .)
      if [ $check_ortho_data -gt 0 ]
      then
        for trans in $transcripts
        do
          echo "Debug: "$gen" -- "$trans
          check_SL=$(grep -w -c $trans $work_directory"/"$Result_Storage"/Mediciones_Expresion/Background_TPM/"$spe"_per_gene.txt")
          if [ $check_SL -gt 0 ]
          then
            pho_id=$(grep -w -F $trans $PHO_File | awk -F "\t" '{print $1}')
            check_busco_targets=$(grep -c -w -F $trans $work_directory"/"$Result_Storage"/BUSCO_RUNS/"$spe"/run_metazoa_odb10/full_table.tsv")

            if [ $check_busco_targets -gt 0 ]
            then
              busco_id=$(grep -w -F $trans $work_directory"/"$Result_Storage"/BUSCO_RUNS/"$spe"/run_metazoa_odb10/full_table.tsv" | awk -F "\t" '{print $1}')
              busco_status=$(grep -w -F $trans $work_directory"/"$Result_Storage"/BUSCO_RUNS/"$spe"/run_metazoa_odb10/full_table.tsv" | awk -F "\t" '{print $2}')
              N_genes_busco=$(grep -c -w -F $busco_id $info_busco_gene_markers)
            else
              busco_id=NA
              busco_status=NA
              N_genes_busco=0
            fi
            echo $gen" "$trans" "$chimera" "$pho_id" "$is_operon" "$ID_operon" "$busco_id" "$busco_status" "$N_genes_busco | tr " " "\t" >> $work_directory"/"$Result_Storage"/BUSCO_Genes_SL/"$spe"_SL_genes_busco.txt"
          fi
        done
      else
        pho_id=Skipped
        chimera=NA
        busco_id=Skipped
        busco_status=Skipped
        N_genes_busco=Skipped

        echo $gen" "$trans" "$chimera" "$pho_id" "$is_operon" "$ID_operon" "$busco_id" "$busco_status" "$N_genes_busco | tr " " "\t" >> $work_directory"/"$Result_Storage"/BUSCO_Genes_SL/"$spe"_SL_genes_busco.txt"
      fi
    done

    N_genes=$(grep -v "N_Genes_Busco" $work_directory"/"$Result_Storage"/BUSCO_Genes_SL/"$spe"_SL_genes_busco.txt" | awk -F "\t" '{print $2}' | grep -c .)
    N_Busco_markers=$(grep -v Missing $work_directory"/"$Result_Storage"/BUSCO_RUNS/"$spe"/run_metazoa_odb10/full_table.tsv" | grep -v -F "#" | grep -c .)
    N_Busco_markers_SL=$(grep -v "N_Genes_Busco" $work_directory"/"$Result_Storage"/BUSCO_Genes_SL/"$spe"_SL_genes_busco.txt" | awk -F "\t" '{if ($7!="NA") print}' | grep -c .)

    N_genes_in_Operon=$(grep -w -c "is_operon" $work_directory"/"$Result_Storage"/BUSCO_Genes_SL/"$spe"_SL_genes_busco.txt")
    N_Busco_markers_in_Operon=$(grep -w "is_operon" $work_directory"/"$Result_Storage"/BUSCO_Genes_SL/"$spe"_SL_genes_busco.txt" | awk -F "\t" '{if ($7!="NA") print}' | grep -c .)

    echo $spe" "$N_genes" "$N_Busco_markers" "$N_Busco_markers_SL" "$N_genes_in_Operon" "$N_Busco_markers_in_Operon | tr " " "\t" >> $work_directory"/"$Result_Storage"/BUSCO_Genes_SL/Summary_SL_genes_busco.txt"
  done
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$RUN_TABLAS_PAPER" == "TRUE" ]
then
  rm -r $work_directory"/"$Result_Storage"/Tablas_Sup_Paper"
  mkdir $work_directory"/"$Result_Storage"/Tablas_Sup_Paper"

  # 1) Tabla general con los resultados
  groups=$(echo "Cestodes Trematoda")

  echo "Lineage Species SRA Trinity_Transcripts N50 Trinity_Genes" | tr " " "\t" >> $work_directory"/"$Result_Storage"/Tablas_Sup_Paper/Table_Data.txt"
  for lineage in $groups
  do
    species=$(ls -d $lineage"/"*"/" | awk -F "/" '{print $2}')
    for spe in $species
    do
      seqkit stats -T -a $work_directory"/"$Result_Storage"/Tablas_Sup_Paper/"* > $work_directory"/"$Result_Storage"/Tablas_Sup_Paper/seqkit.temp"

      rnaseq=$(ls $lineage"/"$spe"/"*"_Trinity.Trinity.fasta.gz" | sed "s/.*\///" | sed "s/_Trinity.Trinity.fasta.gz//")
      for rna in $rnaseq
      do
        Total_Sequences=$(grep $rna"_Trinity.Trinity.fasta.gz" $work_directory"/"$Result_Storage"/Tablas_Sup_Paper/seqkit.temp" | awk -F "\t" '{print $4}')
        N50=$(grep $rna"_Trinity.Trinity.fasta.gz" $work_directory"/"$Result_Storage"/Tablas_Sup_Paper/seqkit.temp" | awk -F "\t" '{print $13}')
        Trinity_Genes=$(seqkit seq -i -n $lineage"/"$spe"/"$rna"_Trinity.Trinity.fasta.gz" | sed "s/_i.*//" | sort -u | grep -c .)

        echo $lineage" "$spe" "$rna" "$Total_Sequences" "$N50" "$Trinity_Genes | tr " " "\t" >> $work_directory"/"$Result_Storage"/Tablas_Sup_Paper/Table_Data.txt"
      done
    done
  done
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$DEFINE_CLUSTER" == "TRUE" ]
then
  rm -r $work_directory"/"$Result_Storage"/Cluster_Final_Definition"
  mkdir $work_directory"/"$Result_Storage"/Cluster_Final_Definition"
  mkdir $work_directory"/"$Result_Storage"/Cluster_Final_Definition/temp"

  Cestoda=$(echo "Egranulosus Emultilocularis Hdiminuta Hmicrostoma Mesocestoides_corti Sparganum_proliferum Spirometra_erinaceieuropaei Taenia_asiatica Taenia_multiceps Taenia_saginata Taenia_solium Schistocephalus_solidus")

  echo "ClusterID Group Species Chr Strand Start End N_loci Average_Distance LocisID SL_Tags Unique_SLRNA Hairpin_Structure" | tr " " "\t" >> $work_directory"/"$Result_Storage"/Cluster_Final_Definition/SL_Cluster_coordinates.txt"
  awk -F "\t" '{if ($13!="NA") print}' "Analysis_Results/Final_SLRNA_information/Final_loci_Table.txt" >> $work_directory"/"$Result_Storage"/Cluster_Final_Definition/Work_SL_RNA_Loci.txt"

  N_current_cluster=1
  check_grupo=$(sed -n 2p $work_directory"/"$Result_Storage"/Cluster_Final_Definition/Work_SL_RNA_Loci.txt" | awk -F "\t" '{print $1}')
  check_species=$(sed -n 2p $work_directory"/"$Result_Storage"/Cluster_Final_Definition/Work_SL_RNA_Loci.txt" | awk -F "\t" '{print $2}')
  check_loci=$(sed -n 2p $work_directory"/"$Result_Storage"/Cluster_Final_Definition/Work_SL_RNA_Loci.txt" | awk -F "\t" '{print $3}')
  check_chr=$(sed -n 2p $work_directory"/"$Result_Storage"/Cluster_Final_Definition/Work_SL_RNA_Loci.txt" | awk -F "\t" '{print $6}')
  check_strand=$(sed -n 2p $work_directory"/"$Result_Storage"/Cluster_Final_Definition/Work_SL_RNA_Loci.txt" | awk -F "\t" '{print $7}')
  check_start=$(sed -n 2p $work_directory"/"$Result_Storage"/Cluster_Final_Definition/Work_SL_RNA_Loci.txt" | awk -F "\t" '{print $12}')
  check_end=$(sed -n 2p $work_directory"/"$Result_Storage"/Cluster_Final_Definition/Work_SL_RNA_Loci.txt" | awk -F "\t" '{print $13}')
  is_cluster=FALSE

  count=3
  recorrer=$(grep -c . $work_directory"/"$Result_Storage"/Cluster_Final_Definition/Work_SL_RNA_Loci.txt")

  while [ $count -le $recorrer ]
  do
    echo $count" // "$recorrer

    new_species=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Cluster_Final_Definition/Work_SL_RNA_Loci.txt" | awk -F "\t" '{print $1}')
    new_loci=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Cluster_Final_Definition/Work_SL_RNA_Loci.txt" | awk -F "\t" '{print $2}')
    new_chr=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Cluster_Final_Definition/Work_SL_RNA_Loci.txt" | awk -F "\t" '{print $4}')
    new_strand=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Cluster_Final_Definition/Work_SL_RNA_Loci.txt" | awk -F "\t" '{print $5}')
    new_start=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Cluster_Final_Definition/Work_SL_RNA_Loci.txt" | awk -F "\t" '{print $14}')
    new_end=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Cluster_Final_Definition/Work_SL_RNA_Loci.txt" | awk -F "\t" '{print $15}')

    check_cestodes=$(echo $Cestoda | grep -c $spe)
    if [ $check_cestodes -gt 0 ]
    then
      new_grupo=Cestodes
    else
      new_grupo=Trematoda
    fi

    distance_previous=$(($new_start - $check_end))

    if [ "$new_species" == "$check_species" ] && [ "$new_chr" == "$check_chr" ] && [ "$new_strand" == "$check_strand" ] && [ "$distance_previous" -le "$max_sl_rna_loci_distance" ] && [ "$distance_previous" -ge 0 ]
    then
      if [ $is_cluster == "FALSE" ]
      then
        grupo_cluster=$check_grupo
        spe_cluster=$check_species
        chr_cluster=$check_chr
        strand_cluster=$check_strand
        start_cluster=$check_start
        cluster_id=$(echo "SL_RNA_Cluster-"$N_current_cluster)
        registro_loci=$(echo $check_loci";"$new_loci)
        N_SL_RNA_loci=2
        registro_separador=$distance_previous

        is_cluster=TRUE
      else
        N_SL_RNA_loci=$(($N_SL_RNA_loci + 1))
        registro_loci=$(echo $registro_loci";"$new_loci)
        registro_separador=$(echo $registro_separador" "$distance_previous)
      fi
    else
      if [ $is_cluster == "TRUE" ]
      then
        is_cluster=FALSE
        end_cluster=$check_end
        promedio_distancia=$(echo $registro_separador | tr " " "\n" | awk '{ sum += $1 } END { if (NR > 0) print sum/NR}' )

        N_current_cluster=$(($N_current_cluster + 1))

        echo $registro_loci | tr ";" "\n" > $work_directory"/"$Result_Storage"/Cluster_Final_Definition/temp/loci_ids.tmp"
        SL_tags=$(grep -w $spe_cluster $work_directory"/"$Result_Storage"/Cluster_Final_Definition/Work_SL_RNA_Loci.txt" | grep -w -f $work_directory"/"$Result_Storage"/Cluster_Final_Definition/temp/loci_ids.tmp" | awk -F "\t" '{print $8}'| sort | uniq -c | awk '{print $2"("$1")"}' | tr "\n" ";" | sed "s/;$/\n/")
        Unique_SLRNA=$(grep -w $spe_cluster $work_directory"/"$Result_Storage"/Cluster_Final_Definition/Work_SL_RNA_Loci.txt" | grep -w -f $work_directory"/"$Result_Storage"/Cluster_Final_Definition/temp/loci_ids.tmp" | awk -F "\t" '{print $13}'| sort | uniq -c | awk '{print $2"("$1")"}' | tr "\n" ";" | sed "s/;$/\n/")
        Hairpin_Structure=$(grep -w $spe_cluster $work_directory"/"$Result_Storage"/Cluster_Final_Definition/Work_SL_RNA_Loci.txt" | grep -w -f $work_directory"/"$Result_Storage"/Cluster_Final_Definition/temp/loci_ids.tmp" | awk -F "\t" '{print $16}'| sort | uniq -c | awk '{print $2"("$1")"}' | tr "\n" ";" | sed "s/;$/\n/")


        echo $cluster_id" "$grupo_cluster" "$spe_cluster" "$chr_cluster" "$strand_cluster" "$start_cluster" "$end_cluster" "$N_SL_RNA_loci" "$promedio_distancia" "$registro_loci" "$SL_tags" "$Unique_SLRNA" "$Hairpin_Structure | tr " " "\t" >> $work_directory"/"$Result_Storage"/Cluster_Final_Definition/SL_Cluster_coordinates.txt"
      fi
    fi
    check_grupo=$new_grupo
    check_species=$new_species
    check_loci=$new_loci
    check_chr=$new_chr
    check_strand=$new_strand
    check_start=$new_start
    check_end=$new_end

    count=$(($count + 1))
  done
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$CONTEXT_CLUSTER" == "TRUE" ]
then
  rm -r $work_directory"/"$Result_Storage"/Context_Clusters_HOG"
  mkdir $work_directory"/"$Result_Storage"/Context_Clusters_HOG"
  mkdir $work_directory"/"$Result_Storage"/Context_Clusters_HOG/Temp"

  results_orthofinder=$(ls -d $work_directory"/"$Result_Storage"/OrthoFinder/Results_"*)
  PHO_File=$(echo $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")

  interval=100000

  count=2
  recorrer=$(grep -c . $work_directory"/"$Result_Storage"/Cluster_Final_Definition/SL_Cluster_coordinates.txt")

  while [ $count -le $recorrer ]
  do
    echo "Discovery fase: "$count" // "$recorrer

    start_compromised=FALSE
    end_compromised=FALSE

    cluster_id=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Cluster_Final_Definition/SL_Cluster_coordinates.txt" | awk -F "\t" '{print $1}')
    group=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Cluster_Final_Definition/SL_Cluster_coordinates.txt" | awk -F "\t" '{print $2}')
    spe=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Cluster_Final_Definition/SL_Cluster_coordinates.txt" | awk -F "\t" '{print $3}')
    chr=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Cluster_Final_Definition/SL_Cluster_coordinates.txt" | awk -F "\t" '{print $4}')
    strand=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Cluster_Final_Definition/SL_Cluster_coordinates.txt" | awk -F "\t" '{print $5}')
    start=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Cluster_Final_Definition/SL_Cluster_coordinates.txt" | awk -F "\t" '{print $6}')
    end=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Cluster_Final_Definition/SL_Cluster_coordinates.txt" | awk -F "\t" '{print $7}')
    N_SL_Loci=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Cluster_Final_Definition/SL_Cluster_coordinates.txt" | awk -F "\t" '{print $8}')
    N_average_distance=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Cluster_Final_Definition/SL_Cluster_coordinates.txt" | awk -F "\t" '{print $9}')
    SL_tags=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Cluster_Final_Definition/SL_Cluster_coordinates.txt" | awk -F "\t" '{print $11}')

    genome_file=$(ls $group"/"$spe"/"*"genomic.fa")
    chr_end=$(seqkit grep -p $chr $genome_file | seqkit fx2tab -n -i -l | awk -F "\t" '{print $2}')

    work_start=$(($start - $interval))
    work_end=$(($end + $interval))

    if [ $work_start -le 0 ]
    then
      work_start=1
      start_compromised=TRUE
    fi

    if [ $work_end -gt $chr_end ]
    then
      work_end=$chr_end
      end_compromised=TRUE
    fi

    echo "## Especie: "$spe >> $work_directory"/"$Result_Storage"/Context_Clusters_HOG/"$cluster_id"_hog.txt"
    echo "## Location: "$chr":"$start".."$end"("$strand")" >> $work_directory"/"$Result_Storage"/Context_Clusters_HOG/"$cluster_id"_hog.txt"
    echo "## Nº Loci: "$N_SL_Loci"("$N_average_distance ")" >> $work_directory"/"$Result_Storage"/Context_Clusters_HOG/"$cluster_id"_hog.txt"
    echo "## SL Tags: "$SL_tags >> $work_directory"/"$Result_Storage"/Context_Clusters_HOG/"$cluster_id"_hog.txt"
    echo "" >> $work_directory"/"$Result_Storage"/Context_Clusters_HOG/"$cluster_id"_hog.txt"
    echo "## Review: "$chr":"$work_start".."$work_end"("$strand")" >> $work_directory"/"$Result_Storage"/Context_Clusters_HOG/"$cluster_id"_hog.txt"
    echo "## Start Compromised: "$start_compromised >> $work_directory"/"$Result_Storage"/Context_Clusters_HOG/"$cluster_id"_hog.txt"
    echo "## End Compromised: "$end_compromised >> $work_directory"/"$Result_Storage"/Context_Clusters_HOG/"$cluster_id"_hog.txt"
    echo "" >> $work_directory"/"$Result_Storage"/Context_Clusters_HOG/"$cluster_id"_hog.txt"
    echo "Location HOG Gen Position" | tr " " "\t" >> $work_directory"/"$Result_Storage"/Context_Clusters_HOG/"$cluster_id"_hog.txt"

    #1) upstream genes
    location=upstream
    if [ $strand == "plus" ]
    then
      interval_start=$work_start
      interval_end=$start
      echo $cluster_id" ("$location") "$chr":"$interval_start".."$interval_end" || "$strand
      interval_genes=$(gffread -r $chr":"$interval_start".."$interval_end $work_directory"/"$Result_Storage"/Modified_GTF_Files/"$spe"_with_Chimera_genes.gtf" | grep -v "#" | awk -F "\t" '{if ($3=="transcript") print $4"\t"$9}' | sort -n -k 1 | awk -F "ID=" '{print $2}')
    else
      interval_start=$end
      interval_end=$work_end
      echo $cluster_id" ("$location") "$chr":"$interval_start".."$interval_end" || "$strand
      interval_genes=$(gffread -r $chr":"$interval_start".."$interval_end $work_directory"/"$Result_Storage"/Modified_GTF_Files/"$spe"_with_Chimera_genes.gtf" | grep -v "#" | awk -F "\t" '{if ($3=="transcript") print $4"\t"$9}' | sort -n -r -k 1 | awk -F "ID=" '{print $2}')
    fi

    check_hog_order

    #3) Internal genes
    location=internal
    interval_start=$start
    interval_end=$end
    echo $cluster_id" ("$location") "$chr":"$interval_start".."$interval_end" || "$strand
    interval_genes=$(gffread -R -r $chr":"$interval_start".."$interval_end $work_directory"/"$Result_Storage"/Modified_GTF_Files/"$spe"_with_Chimera_genes.gtf" | grep -v "#" | awk -F "\t" '{if ($3=="transcript") print $4"\t"$9}' | sort -n -k 1 | awk -F "ID=" '{print $2}')

    check_hog_order

    #3) downstream genes
    location=downstream

    if [ $strand == "plus" ]
    then
      interval_start=$end
      interval_end=$work_end
      echo $cluster_id" ("$location") "$chr":"$interval_start".."$interval_end" || "$strand
      interval_genes=$(gffread -r $chr":"$interval_start".."$interval_end $work_directory"/"$Result_Storage"/Modified_GTF_Files/"$spe"_with_Chimera_genes.gtf" | grep -v "#" | awk -F "\t" '{if ($3=="transcript") print $4"\t"$9}' | sort -n -k 1 | awk -F "ID=" '{print $2}')
    else
      interval_start=$work_start
      interval_end=$start
      echo $cluster_id" ("$location") "$chr":"$interval_start".."$interval_end" || "$strand
      interval_genes=$(gffread -r $chr":"$interval_start".."$interval_end $work_directory"/"$Result_Storage"/Modified_GTF_Files/"$spe"_with_Chimera_genes.gtf" | grep -v "#" | awk -F "\t" '{if ($3=="transcript") print $4"\t"$9}' | sort -n -r -k 1 | awk -F "ID=" '{print $2}')
    fi

    check_hog_order

    count=$(($count + 1))
  done

  # Heatmaps:
  all_cluster=$(awk -F "\t" '{print $1}' $work_directory"/"$Result_Storage"/Cluster_Final_Definition/SL_Cluster_coordinates.txt" | sed 1d)

  #1) Upstreap
  search_term=upstream
  header=Cluster
  file_name=Upstream_HOGs
  cluster_heatmap

  #2) Downstream
  search_term=downstream
  header=Cluster
  file_name=Downstream_HOGs
  cluster_heatmap

  #3) All
  search_term=stream
  header=Cluster
  file_name=All_HOGs
  cluster_heatmap
fi

################################################################################
################################################################################
################################################################################

if [ "$SUP_TABLE_DATA_USED" == "TRUE" ]
then
  rm -r $work_directory"/"$Result_Storage"/Used_Data_Table"
  mkdir $work_directory"/"$Result_Storage"/Used_Data_Table"

  grupos=$(echo "Cestodes Trematoda")
  echo "Grupo Species GenomeID RNA_ID Work_Reads"  | tr " " "\t" >> $work_directory"/"$Result_Storage"/Used_Data_Table/Sup_Table_Work_Data.txt"
  for gru in $grupos
  do
    species=$(ls -d $gru"/"*"/" | awk -F "/" '{print $2}')
    for spe in $species
    do
      check_genome_file=$(ls $gru"/"$spe"/"* | grep -c "genomic.fa")
      if [ $check_genome_file -gt 0 ]
      then
        check_wormbase=$(ls $gru"/"$spe"/"* | grep -c "WBPS15.genomic.fa")
        if [ $check_wormbase -gt 0 ]
        then
          echo "Roto 1":
          genomeID=$(ls $gru"/"$spe"/"*"genomic.fa" | sed "s/.*\///" | awk -F "." '{print $2}')
          echo "Done"
        else
          echo "Roto 2":
          genomeID=$(ls $gru"/"$spe"/"*"genomic.fa" | sed "s/.*\///" | awk -F "_" '{print $1"_"$2}')
          echo "Done"
        fi
        echo "samples"
        samples=$(ls $gru"/"$spe"/"*"_work_read_1P.gz" | sed "s/.*\///" | sed "s/_work_read_1P.gz//")

        echo "done"
        for sam in $samples
        do
          work_reads=$(seqkit stats -T $gru"/"$spe"/"$sam"_work_read_1P.gz" | tail -n 1 | awk -F "\t" '{print $4}')
          echo $gru" "$spe" "$genomeID" "$sam" "$work_reads
          echo $gru" "$spe" "$genomeID" "$sam" "$work_reads | tr " " "\t" >> $work_directory"/"$Result_Storage"/Used_Data_Table/Sup_Table_Work_Data.txt"
        done
      fi
    done
  done
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$SITE_MAIN_SUMMARY" == "TRUE" ]
then
  rm -r $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY"
  mkdir $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY"
  mkdir $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/Temporales"

  grupo=Cestodes
  Especies=$(echo "Egranulosus Emultilocularis Hdiminuta Hmicrostoma Mesocestoides_corti Sparganum_proliferum Spirometra_erinaceieuropaei Taenia_asiatica Taenia_multiceps Taenia_saginata Taenia_solium Schistocephalus_solidus")
  final_site_table

  grupo=Trematoda
  Especies=$(echo "Clonorchis_sinensis Fasciola_gigantica Fasciola_hepatica Fasciolopsis_buski Opisthorchis_felineus Paragonimus_heterotremus Paragonimus_westermani Schistosoma_bovis Schistosoma_haematobium Schistosoma_japonicum Schistosoma_mansoni Trichobilharzia_regenti")
  final_site_table
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$CHIMERA_SUMARY" == "TRUE" ]
then
  rm -r $work_directory"/"$Result_Storage"/Chimeric_summary"
  mkdir $work_directory"/"$Result_Storage"/Chimeric_summary"

  Especies=$(ls $work_directory"/"$Result_Storage"/"*pep | sed "s/.*\///" | sed "s/.pep//")
  echo "Species Gene_Model Transcript_Start Transcript_End Marching_Hits Start End Length Complete Chimeric_Portion" | tr " " "\t" >> $work_directory"/"$Result_Storage"/Chimeric_summary/Summary_Chimera.txt"
  for spe in $Especies
  do
    chimeric_genes=$(grep "Chimera_" $work_directory"/"$Result_Storage"/"$spe".pep" | sed "s/_Chimera.*//" | tr -d ">" | sort -u)
    for gen in $chimeric_genes
    do
      echo $spe" "$gen
      chimeric_portions=$(grep $gen"_Chimera_" $work_directory"/"$Result_Storage"/"$spe".pep" | tr -d ">")
      grep -w $gen $work_directory"/"$Result_Storage"/Chimera_Details/"$spe"/All_cut_points.txt" > $work_directory"/"$Result_Storage"/Chimeric_summary/cut_points_tmp.txt"

      count=1
      recorrer=$(grep -c . $work_directory"/"$Result_Storage"/Chimeric_summary/cut_points_tmp.txt")

      continue_hit=PLACE_HOLDER
      while [ $count -le $recorrer ]
      do
        Complete_hit=$(echo NA)
        information1=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Chimeric_summary/cut_points_tmp.txt" | awk -F "\t" '{print $1" "$2" "$3}')
        information2=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Chimeric_summary/cut_points_tmp.txt" | awk -F "\t" '{print $5" "$6" "$7}')

        check_matching_genes=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Chimeric_summary/cut_points_tmp.txt" | awk -F "\t" '{print $4}' | grep -c .)
        if [ $check_matching_genes -gt 0 ]
        then
          matching_genes=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Chimeric_summary/cut_points_tmp.txt" | awk -F "\t" '{print $4}')
        else
          matching_genes=NA
        fi

        start=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Chimeric_summary/cut_points_tmp.txt" | awk -F "\t" '{print $5}')
        end=$(sed -n $count"p" $work_directory"/"$Result_Storage"/Chimeric_summary/cut_points_tmp.txt" | awk -F "\t" '{print $6}')

        #1) Complete_hit
        check_complete=$(echo $chimeric_portions | grep -c "__"$start"_"$end)
        if [ $check_complete -gt 0 ]
        then
          chimeric_part=$(echo $chimeric_portions | tr " " "\n" | grep "__"$start"_"$end)
          Complete_hit=$(echo X)
          chimeric_portions=$(echo $chimeric_portions | tr " " "\n" | grep -v -w -F $chimeric_part)
        else
          # Confirmar hit parcial
          check_start=$(echo $chimeric_portions | grep -c "__"$start"_")
          if [ $check_start -gt 0 ]
          then
            chimeric_part=$(echo $chimeric_portions | grep "__"$start"_")
            Complete_hit=$(echo .)
            continue_hit=TRUE
          fi

          check_end=$(echo $chimeric_portions | grep -c "_"$end$)
          if [ $check_end -gt 0 ] && [ $continue_hit == "TRUE" ]
          then
            chimeric_part=$(echo $chimeric_portions | grep "_"$end$)
            Complete_hit=$(echo .)
            continue_hit=FALSE
          fi
        fi

        if [ $Complete_hit == "NA" ]
        then
          chimeric_part=NA
        fi

        echo "information1: "$information1
        echo "matching_genes: "$matching_genes
        echo "information2: "$information2
        echo "Complete_hit: "$Complete_hit
        echo "chimeric_part: "$chimeric_part
        echo ""
        echo "###############################################"
        echo ""

        echo $spe" "$information1" "$matching_genes" "$information2" "$Complete_hit" "$chimeric_part | sed "s/ABORT/Skip/" | tr " " "\t" >> $work_directory"/"$Result_Storage"/Chimeric_summary/Summary_Chimera.txt"

        count=$(($count + 1))
      done
    done
  done
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ $NEW_N0_TABLE == "TRUE" ]
then
  rm -r $work_directory"/"$Result_Storage"/Orthofinder_Sup_Table"
  mkdir $work_directory"/"$Result_Storage"/Orthofinder_Sup_Table"
  mkdir $work_directory"/"$Result_Storage"/Orthofinder_Sup_Table/Temp"

  results_orthofinder=$(ls -d $work_directory"/"$Result_Storage"/OrthoFinder/Results_"*)
  PHO_File=$(echo $results_orthofinder"/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")

  Especies=$(ls $work_directory"/"$Result_Storage"/"*pep | sed "s/.*\///" | sed "s/.pep//")

  header=$(echo "HOG OG Gene_Tree_Parent_Clade Interpro_IDs")
  for spe in $Especies
  do
    add=$(echo $spe"_Total_Genes "$spe"_Genes_SL "$spe"_Gene_IDs" )
    header=$(echo $header" "$add)
  done

  echo $header | tr " " "\t" > $work_directory"/"$Result_Storage"/Orthofinder_Sup_Table/New_N0_Sup_Table.txt"

  count=2
  recorrer=$(grep -c . $PHO_File)

  while [ $count -le $recorrer ]
  do
    hog=$(sed -n $count"p" $PHO_File | awk -F "\t" '{print $1}')
    HOG_info=$(sed -n $count"p" $PHO_File | awk -F "\t" '{print $1" "$2" "$3}')
    registro_HOG=XXXX

    echo $hog": "$count"/"$recorrer

    for spe in $Especies
    do
      spe_col=$(head -n 1 $PHO_File | tr "\t" "\n" | grep -n $spe | awk -F ":" '{print $1}')
      N_spe_genes=$(grep -w $hog $PHO_File | awk -F "\t" -v N=$spe_col '{print $N}' | tr -d " " | tr "," "\n" | grep -c .)

      echo "NA" > $work_directory"/"$Result_Storage"/Orthofinder_Sup_Table/Temp/current_interprot_list.tmp"
      if [ $N_spe_genes -gt 0 ]
      then
        grep -w $hog $PHO_File | awk -F "\t" -v N=$spe_col '{print $N}' | tr -d " " | tr "," "\n" > $work_directory"/"$Result_Storage"/Orthofinder_Sup_Table/Temp/current_gene_list.tmp"
        grep -w -F -f $work_directory"/"$Result_Storage"/Orthofinder_Sup_Table/Temp/current_gene_list.tmp" $work_directory"/"$Result_Storage"/Interprot_Results/"$spe".tsv" | awk -F "\t" '{print $12}' | grep IPR | sort -u >> $work_directory"/"$Result_Storage"/Orthofinder_Sup_Table/Temp/current_interprot_list.tmp"

        genes_ids=$(grep -w $hog $PHO_File | awk -F "\t" -v N=$spe_col '{print $N}' | tr -d " ")
        N_spe_SL_genes=$(grep -w -F -c -f $work_directory"/"$Result_Storage"/Orthofinder_Sup_Table/Temp/current_gene_list.tmp" $work_directory"/"$Result_Storage"/Mediciones_Expresion/Background_TPM/"$spe"_per_gene.txt")
      else
        genes_ids=NA
        N_spe_SL_genes=0
      fi
      registro_HOG=$(echo $registro_HOG" "$N_spe_genes" "$N_spe_SL_genes" "$genes_ids)
    done

    registro_HOG=$(echo $registro_HOG | sed "s/XXXX //")

    check_interpro=$(grep -v -w -c NA $work_directory"/"$Result_Storage"/Orthofinder_Sup_Table/Temp/current_interprot_list.tmp")
    if [ $check_interpro -gt 0 ]
    then
      interpro_annotation=$(sort -u $work_directory"/"$Result_Storage"/Orthofinder_Sup_Table/Temp/current_interprot_list.tmp" | grep -v NA | tr "\n" "," | sed "s/,$//")
    else
      interpro_annotation=NA
    fi

    echo $HOG_info" "$interpro_annotation" "$registro_HOG | tr " " "\t" >> $work_directory"/"$Result_Storage"/Orthofinder_Sup_Table/New_N0_Sup_Table.txt"
    count=$(($count +1))
  done
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ $ITOL_SPECIES_TREE_DATA == "TRUE" ]
then
  rm -r $work_directory"/"$Result_Storage"/ITOL_Species_Tree_Data"
  mkdir $work_directory"/"$Result_Storage"/ITOL_Species_Tree_Data"

  Especies=$(ls $work_directory"/"$Result_Storage"/"*pep | sed "s/.*\///" | sed "s/.pep//")
  Cestoda=$(echo "Egranulosus Emultilocularis Hdiminuta Hmicrostoma Mesocestoides_corti Sparganum_proliferum Spirometra_erinaceieuropaei Taenia_asiatica Taenia_multiceps Taenia_saginata Taenia_solium Schistocephalus_solidus")
  echo "Especie N_Genes N_Genes_SL N_Operones N_Operones_Chimera N_Gene_Model" | tr " " "\t" >> $work_directory"/"$Result_Storage"/ITOL_Species_Tree_Data/Itol_data.txt"

  for spe in $Especies
  do
    echo $spe
    check_cestodes=$(echo $Cestoda | grep -c $spe)
    if [ $check_cestodes -gt 0 ]
    then
      grupo=Cestodes
    else
      grupo=Trematoda
    fi

    gff_file=$(ls $grupo"/"$spe"/"*"annotations.gtf")

    all_genes=$(awk -F "\t" '{if ($3=="CDS") print $9}' $gff_file | tr ";" "\n" | grep -w gene_id | sort -u  | grep -c .)
    genes_SL=$(grep -w Concordante $work_directory"/SL_Read_Counts/"$spe"_tier3.counts" | grep -w -v Trans_End | awk -F "\t" '{print $1}' | sort -u | grep -c .)
    all_operones=$(grep -c -v -w "OpID" $work_directory"/"$Result_Storage"/Operons/"$spe"_operon.txt")
    chimera_operones=$(grep -c -F "[C]" $work_directory"/"$Result_Storage"/Operons/"$spe"_operon.txt")
    operon_gene_models=$(grep -v -w "OpID" $work_directory"/"$Result_Storage"/Operons/"$spe"_operon.txt" | awk -F "\t" '{print $6}' | sed "s/_(/\n/g" | grep -c .)

    echo $spe" "$all_genes" "$genes_SL" "$all_operones" "$chimera_operones" "$operon_gene_models
    echo $spe" "$all_genes" "$genes_SL" "$all_operones" "$chimera_operones" "$operon_gene_models | tr " " "\t" >> $work_directory"/"$Result_Storage"/ITOL_Species_Tree_Data/Itol_data.txt"
  done
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$SL_READS_FIGURE_TAB_MAKE" == TRUE ]
then
  rm -r $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/Info_Figura"
  mkdir $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/Info_Figura"

  Especies=$(ls $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/"*"_SL_Site_Full_Summary_work.txt" | sed "s/.*\///" | sed "s/_SL_Site_Full_Summary_work.txt//")
  echo "Especie Grupo N_Sitios" | tr " " "\t" >> $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/Info_Figura/Grafico_stacked_cols_Sitios.txt"
  echo "Especie Grupo N_SL_Reads" | tr " " "\t" >> $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/Info_Figura/Grafico_stacked_cols_SLReads.txt"
  echo "Especie Grupo N_SL_splicing" | tr " " "\t" >> $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/Info_Figura/Grafico_stacked_cols_splicing.txt"

  for spe in $Especies
  do
    # Sitios:
    echo $spe
    N_Start_sites=$(awk -F "\t" '{if ($11=="Start" && $8=="1_of_1" && $10 =="X") print $1"__"$2}' $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/"$spe"_SL_Site_Full_Summary_work.txt" | sort -u | grep -c .)
    N_Internal_sites=$(awk -F "\t" '{if ($11=="Internal" && $8=="1_of_1" && $10 =="X") print $1"__"$2}' $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/"$spe"_SL_Site_Full_Summary_work.txt" | sort -u | grep -c .)
    N_End_sites=$(awk -F "\t" '{if ($11=="End" && $8=="1_of_1" && $10 =="X") print $1"__"$2}' $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/"$spe"_SL_Site_Full_Summary_work.txt" | sort -u | grep -c .)

    echo $spe" Upstream "$N_Start_sites | tr " " "\t" >> $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/Info_Figura/Grafico_stacked_cols_Sitios.txt"
    echo $spe" Internal "$N_Internal_sites | tr " " "\t" >> $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/Info_Figura/Grafico_stacked_cols_Sitios.txt"
    echo $spe" Downstream "$N_End_sites | tr " " "\t" >> $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/Info_Figura/Grafico_stacked_cols_Sitios.txt"

    ##############################

    N_Start_SL_reads=$(awk -F "\t" '{if ($11=="Start" && $8=="1_of_1" && $10 =="X") print $4 }' $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/"$spe"_SL_Site_Full_Summary_work.txt" | awk -F "\t" '{ sum += $1 } END {print sum}')
    N_Internal_SL_reads=$(awk -F "\t" '{if ($11=="Internal" && $8=="1_of_1" && $10 =="X") print $4 }' $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/"$spe"_SL_Site_Full_Summary_work.txt" | awk -F "\t" '{ sum += $1 } END {print sum}')

    if [ $N_End_sites -gt 0 ]
    then
      N_End_SL_reads=$(awk -F "\t" '{if ($11=="End" && $8=="1_of_1" && $10 =="X") print $4 }' $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/"$spe"_SL_Site_Full_Summary_work.txt" | awk -F "\t" '{ sum += $1 } END {print sum}')
    else
      N_End_SL_reads=0
    fi

    echo $spe" Upstream "$N_Start_SL_reads | tr " " "\t" >> $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/Info_Figura/Grafico_stacked_cols_SLReads.txt"
    echo $spe" Internal "$N_Internal_SL_reads | tr " " "\t" >> $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/Info_Figura/Grafico_stacked_cols_SLReads.txt"
    echo $spe" Downstream "$N_End_SL_reads | tr " " "\t" >> $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/Info_Figura/Grafico_stacked_cols_SLReads.txt"

    ##############################

    N_Start_SL_splicing=$(awk -F "\t" '{if ($11=="Start" && $8=="1_of_1" && $10 =="X") print $5 }' $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/"$spe"_SL_Site_Full_Summary_work.txt" | awk -F "\t" '{ sum += $1 } END {print sum}')
    N_Internal_SL_splicing=$(awk -F "\t" '{if ($11=="Internal" && $8=="1_of_1" && $10 =="X") print $5 }' $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/"$spe"_SL_Site_Full_Summary_work.txt" | awk -F "\t" '{ sum += $1 } END {print sum}')

    if [ $N_End_sites -gt 0 ]
    then
      N_End_SL_splicing=$(awk -F "\t" '{if ($11=="End" && $8=="1_of_1" && $10 =="X") print $5 }' $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/"$spe"_SL_Site_Full_Summary_work.txt" | awk -F "\t" '{ sum += $1 } END {print sum}')
    else
      N_End_SL_splicing=0
    fi

    echo $spe" Upstream "$N_Start_SL_splicing | tr " " "\t" >> $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/Info_Figura/Grafico_stacked_cols_splicing.txt"
    echo $spe" Internal "$N_Internal_SL_splicing | tr " " "\t" >> $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/Info_Figura/Grafico_stacked_cols_splicing.txt"
    echo $spe" Downstream "$N_End_SL_splicing | tr " " "\t" >> $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/Info_Figura/Grafico_stacked_cols_splicing.txt"
  done

  outputs=$(ls $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/Info_Figura/Grafico_"*)

  for out in $outputs
  do
    sed -i "s/Egranulosus/Echinococcus granulosus/" $out
    sed -i "s/Emultilocularis/Echinococcus multilocularis/" $out
    sed -i "s/Emultilocularis/Echinococcus multilocularis/" $out
    sed -i "s/Hdiminuta/Hymenolepis diminuta/" $out
    sed -i "s/Hmicrostoma/Hymenolepis microstoma/" $out
    sed -i "s/_/ /g" $out
  done
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$UP_SET_EXPLORATION" == "TRUE" ]
then
  rm -r $work_directory"/"$Result_Storage"/Upset_graph"
  mkdir $work_directory"/"$Result_Storage"/Upset_graph"
  mkdir $work_directory"/"$Result_Storage"/Upset_graph/Temp"

  analysis=SL_HOG
  data_path=$work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Heatmap_SL_HOG/Species_List"
  data_id=sl_pho.ids
  excluded_species=$(echo "Mesocestoides_corti Spirometra_erinaceieuropaei")
  upset_exploration

  analysis=Op_Pair_second_with_SL
  data_path=$work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo"
  data_id=operon_pair_sl_second_in_pair_only_pairs.txt
  excluded_species=$(echo "Mesocestoides_corti Spirometra_erinaceieuropaei Trichobilharzia_regenti Paragonimus_heterotremus Fasciolopsis_buski Opisthorchis_felineus")
  upset_exploration
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$UP_SET_INPUT" == "TRUE" ]
then
  Grupos_Cestoda_Derivados=$(echo "Egranulosus Emultilocularis Hdiminuta Hmicrostoma Mesocestoides_corti Taenia_asiatica Taenia_multiceps Taenia_saginata Taenia_solium")
  Grupos_Cestoda_Basal=$(echo "Schistocephalus_solidus Sparganum_proliferum Spirometra_erinaceieuropaei")

  ### Trematoda
  Grupos_Schistosmidae=$(echo "Schistosoma_bovis Schistosoma_haematobium Schistosoma_japonicum Schistosoma_mansoni Trichobilharzia_regenti")
  Grupos_Trematoda_No_schistosoma=$(echo "Clonorchis_sinensis Fasciola_gigantica Fasciola_hepatica Fasciolopsis_buski Opisthorchis_felineus Paragonimus_heterotremus Paragonimus_westermani")

  analysis=SL_HOG
  Especies=$(ls $work_directory"/"$Result_Storage"/"*pep | sed "s/.*\///" | sed "s/.pep//" | grep -v "Spirometra_erinaceieuropaei" | grep -v "Mesocestoides_corti")
  data_path=$work_directory"/"$Result_Storage"/Accepor_Site_Analisis/Heatmap_SL_HOG/Species_List"
  data_id=sl_pho.ids
  upset_input

  Especies=$(ls $work_directory"/"$Result_Storage"/"*pep | sed "s/.*\///" | sed "s/.pep//" | grep -v "Spirometra_erinaceieuropaei" | grep -v "Mesocestoides_corti" | grep -v "Trichobilharzia_regenti" | grep -v "Paragonimus_heterotremus" | grep -v "Fasciolopsis_buski" | grep -v "Opisthorchis_felineus")
  analysis=Op_Pair_second_with_SL
  data_path=$work_directory"/"$Result_Storage"/Operons/Questions/Specific/More_Operon_Heatmaps/Registro_conteo"
  data_id=operon_pair_sl_second_in_pair_only_pairs.txt
  upset_input
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$UP_SET_GROUP_QUESTIONS" == "TRUE" ]
then
  rm -r $work_directory"/"$Result_Storage"/Upset_graph/Operon_question"
  mkdir $work_directory"/"$Result_Storage"/Upset_graph/Operon_question"

  groups=$(ls $work_directory"/"$Result_Storage"/Upset_graph/SL_HOG_Intput/"*"_graph.in" | sed "s/.*\///" | sed "s/_graph.in//" )
  for gru in $groups
  do
    echo $gru
    total_uniq_sl_hogs=$(grep -c . $work_directory"/"$Result_Storage"/Upset_graph/SL_HOG_Intput/"$gru"_uniq.id")
    total_shared_sl_hogs=$(grep -c . $work_directory"/"$Result_Storage"/Upset_graph/SL_HOG_Intput/"$gru"_graph.in" )
    all_sl_hogs=$(cat $work_directory"/"$Result_Storage"/Upset_graph/SL_HOG_Intput/"$gru"_uniq.id" $work_directory"/"$Result_Storage"/Upset_graph/SL_HOG_Intput/"$gru"_graph.in" | sort -u | grep -c .)
    per_uniq_sl_hogs=$(echo $total_uniq_sl_hogs" "$all_sl_hogs | awk -F " " '{printf "%.2f\n", ($1/$2)*100}')

    total_operon_uniq=$(grep -c . $work_directory"/"$Result_Storage"/Upset_graph/Op_Pair_second_with_SL_Intput/"$gru"_uniq.id")
    total_operon_shared=$(grep -c . $work_directory"/"$Result_Storage"/Upset_graph/Op_Pair_second_with_SL_Intput/"$gru"_graph.in")
    all_operon=$(cat $work_directory"/"$Result_Storage"/Upset_graph/Op_Pair_second_with_SL_Intput/"$gru"_uniq.id" $work_directory"/"$Result_Storage"/Upset_graph/Op_Pair_second_with_SL_Intput/"$gru"_graph.in" | sort -u | grep -c .)
    per_operon_uniq=$(echo $total_operon_uniq" "$all_operon | awk -F " " '{printf "%.2f\n", ($1/$2)*100}')

    hogs_in_uniq_operons=$(grep -o -f  $work_directory"/"$Result_Storage"/Upset_graph/SL_HOG_Intput/"$gru"_uniq.id" $work_directory"/"$Result_Storage"/Upset_graph/Op_Pair_second_with_SL_Intput/"$gru"_uniq.id" | sort -u | grep -c .)
    per_hogs_in_uniq_operons=$(echo $hogs_in_uniq_operons" "$total_uniq_sl_hogs | awk -F " " '{printf "%.2f\n", ($1/$2)*100}')

    hogs_in_shared_operons=$(grep -o -f  $work_directory"/"$Result_Storage"/Upset_graph/SL_HOG_Intput/"$gru"_uniq.id" $work_directory"/"$Result_Storage"/Upset_graph/Op_Pair_second_with_SL_Intput/"$gru"_graph.in" | sort -u | grep -c .)
    per_hogs_in_shared_operons=$(echo $hogs_in_shared_operons" "$total_uniq_sl_hogs | awk -F " " '{printf "%.2f\n", ($1/$2)*100}')

    echo $gru >> $work_directory"/"$Result_Storage"/Upset_graph/Operon_question/Summary_Hogs_operones.txt"
    echo "All HOGs with SL: "$all_sl_hogs >> $work_directory"/"$Result_Storage"/Upset_graph/Operon_question/Summary_Hogs_operones.txt"
    echo "Unique: "$total_uniq_sl_hogs" ("$per_uniq_sl_hogs"%)" >> $work_directory"/"$Result_Storage"/Upset_graph/Operon_question/Summary_Hogs_operones.txt"
    echo "Shared: "$total_shared_sl_hogs >> $work_directory"/"$Result_Storage"/Upset_graph/Operon_question/Summary_Hogs_operones.txt"
    echo "" >> $work_directory"/"$Result_Storage"/Upset_graph/Operon_question/Summary_Hogs_operones.txt"
    echo "All HOG pairs in Operons: "$all_operon >> $work_directory"/"$Result_Storage"/Upset_graph/Operon_question/Summary_Hogs_operones.txt"
    echo "Unique: "$total_operon_uniq "("$per_operon_uniq"%)" >> $work_directory"/"$Result_Storage"/Upset_graph/Operon_question/Summary_Hogs_operones.txt"
    echo "Shared: "$total_operon_shared >> $work_directory"/"$Result_Storage"/Upset_graph/Operon_question/Summary_Hogs_operones.txt"
    echo "" >> $work_directory"/"$Result_Storage"/Upset_graph/Operon_question/Summary_Hogs_operones.txt"
    echo "HOG Pairs in unique operons: "$hogs_in_uniq_operons" ("$per_hogs_in_uniq_operons"%)" >> $work_directory"/"$Result_Storage"/Upset_graph/Operon_question/Summary_Hogs_operones.txt"
    echo "HOG Pairs in shared operons: "$hogs_in_shared_operons" ("$per_hogs_in_shared_operons"%)" >> $work_directory"/"$Result_Storage"/Upset_graph/Operon_question/Summary_Hogs_operones.txt"
    echo "" >> $work_directory"/"$Result_Storage"/Upset_graph/Operon_question/Summary_Hogs_operones.txt"
    echo "#########################" >> $work_directory"/"$Result_Storage"/Upset_graph/Operon_question/Summary_Hogs_operones.txt"
    echo "" >> $work_directory"/"$Result_Storage"/Upset_graph/Operon_question/Summary_Hogs_operones.txt"
  done
fi

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if [ "$RAREFACTION_CURVES" == TRUE ]
then
  rm -r $work_directory"/"$Result_Storage"/Rarefaction_curves"
  mkdir $work_directory"/"$Result_Storage"/Rarefaction_curves"

  Especies=$(ls $work_directory"/"$Result_Storage"/"*pep | sed "s/.*\///" | sed "s/.pep//")
  for spe in $Especies
  do
    genes=$(awk -F "\t" '{if ($10 == "X" && $8 == "1_of_1" && $11 != "End") print $9}' $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/"$spe"_SL_Site_Full_Summary_work.txt" | sort -u)
    echo "Especie "$genes | tr " " "\t" >>  $work_directory"/"$Result_Storage"/Rarefaction_curves/"$spe"_rarefaction.in"
    registro=$spe
    for gen in $genes
    do
      echo $spe" "$gen
      SL_read_count=$(grep -w $gen  $work_directory"/"$Result_Storage"/SL_ACE_SUMMARY/"$spe"_SL_Site_Full_Summary_work.txt" | awk -F "\t" '{if ($10 == "X" && $11 != "End") print $4}' | awk -F "\t" '{ sum += $1 } END {print sum}')
      registro=$(echo $registro" "$SL_read_count)
    done
    echo $registro | tr " " "\t" >>  $work_directory"/"$Result_Storage"/Rarefaction_curves/"$spe"_rarefaction.in"
  done
fi
