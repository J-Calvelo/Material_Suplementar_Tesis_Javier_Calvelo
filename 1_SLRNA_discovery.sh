#!/usr/bin/env bash
home_direcory=$(pwd)
grupos=$(echo "Cestodes Trematoda")
uniprot_reference=/home/amanda/PROYECTS/SLFinder_Publication/UNIPROT_29-Junio-2020/uniprot.fasta
use_TruSeq2=$(echo "Taenia_multiceps Taenia_solium Trichobilharzia_regenti Schistosoma_mansoni Schistosoma_japonicum")
trimmomatic_installation=/home/amanda/PROYECTS/SLFinder_Publication/Trimmomatic-0.36
path_to_ninja=/home/amanda/PROYECTS/SLFinder_Publication/ninja_1.2.2
threads=24
SLs_Conocidos_Fasta=SLs_Conocidos.fasta
Adaptadores_Ilumina_Fasta=Adaptadores_Ilumina.fasta
path_to_ninja=/home/amanda/PROYECTS/SLFinder_Publication/ninja_1.2.2
IDEN_HOOK=80
Min_SL_Hit=10
Min_Adapt_Hit=10
NMotifs=20
MEME_Min_Cluster=4
SM_like_sites=Selected_Sequences_Work/SM_Patterns.fasta
Rango_Region=500
TreeCluster_dist=0.7
TreeCluster_Mode=med_clade
Full_SLs=SLs_Conocidos_Full.fasta

## Variables corrida de modulos
PREPARE_READS=FALSE
ASSEMBLY=FALSE
RUN_SLFINDER1=FALSE
RUN_SLFINDER2=FALSE
SLRNA_RETRIVAL=FALSE
RUN_REPEATMASKER=FALSE
SLRNA_INITIAL_FILTERING=TRUE

################################################################################
# Functions for SL-RNA identification and initial evaluation:

get_hook_variants () {
  Total_hooks=$(grep -c ">" $gru"/"$spec"/SL-analysis/Results/Best-Hook.fasta")
  Work_Hook_Variants=$(grep -c ">" $gru"/"$spec"/SL-analysis/Results/Hook_variants/Hook_all_variants.fa")

  Hooks_IDs=$(seqkit seq -n -i $gru"/"$spec"/SL-analysis/Results/Best-Hook.fasta")

  for hook in $Hooks_IDs
  do
    echo "Procesando "$hook
    seqkit grep -r -p $hook $gru"/"$spec"/SL-analysis/Results/Best-Hook.fasta" > "Analysis_Results/Alineamientos_Hooks/"$gru"/Secuencias.temp"

    check_hook_forward=$(ls $gru"/"$spec"/SL-analysis/Results/Hook_variants/Hook_Individual/" | grep forward | grep -c $hook)
    check_hook_reverse=$(ls $gru"/"$spec"/SL-analysis/Results/Hook_variants/Hook_Individual/" | grep reverse | grep -c $hook)

    if [ $check_hook_forward -gt 0 ]
    then
      cat $gru"/"$spec"/SL-analysis/Results/Hook_variants/Hook_Individual/"$hook"-"*"forward-SL.fa" | seqkit seq -u >> "Analysis_Results/Alineamientos_Hooks/"$gru"/Secuencias.temp"
    fi

    if [ $check_hook_reverse -gt 0 ]
    then
      cat $gru"/"$spec"/SL-analysis/Results/Hook_variants/Hook_Individual/"$hook"-"*"reverse-SL.fa" | seqkit seq -u -p -r -v -t DNA >> "Analysis_Results/Alineamientos_Hooks/"$gru"/Secuencias.temp"
    fi

    if [ $check_hook_forward -gt 0 ] || [ $check_hook_reverse -gt 0 ]
    then
      mafft --globalpair --quiet --maxiterate 1000 "Analysis_Results/Alineamientos_Hooks/"$gru"/Secuencias.temp" >> "Analysis_Results/Alineamientos_Hooks/"$gru"/"$spec"_"$hook".aln"
    fi
    rm "Analysis_Results/Alineamientos_Hooks/"$gru"/Secuencias.temp"
  done

  # BLAST chequeo
  blastn -task blastn-short -query $gru"/"$spec"/SL-analysis/Results/Best-Hook.fasta" -db "Analysis_Results/BLAST_Identificacion/Referencias/SLs_Conocidos_Fasta_ref" -outfmt "6 std" -out "Analysis_Results/BLAST_Identificacion/"$gru"/"$spec"_hook_sl_conocidos.blastn" -perc_identity $IDEN_HOOK
  hit_sl_hook=$(awk -F "\t" '{print $1}' "Analysis_Results/BLAST_Identificacion/"$gru"/"$spec"_hook_sl_conocidos.blastn" | sort -u | grep -c .)

  blastn -task blastn-short -query $gru"/"$spec"/SL-analysis/Results/Best-Hook.fasta" -db "Analysis_Results/BLAST_Identificacion/Referencias/Adaptadores_Ilumina_Fasta_ref" -outfmt "6 std" -out "Analysis_Results/BLAST_Identificacion/"$gru"/"$spec"_hook_adaptadores.blastn" -perc_identity $IDEN_HOOK
  hit_adaptadores_hook=$(awk -F "\t" '{print $1}' "Analysis_Results/BLAST_Identificacion/"$gru"/"$spec"_hook_adaptadores.blastn" | sort -u | grep -c .)

  for hook in $Hooks_IDs
  do
    check_SL_Conocido=$(grep $hook "Analysis_Results/BLAST_Identificacion/"$gru"/"$spec"_hook_sl_conocidos.blastn" | awk -F "\t" -v min=$Min_SL_Hit '{if ($4>=min) print}' | grep -c .)
    if [ $check_SL_Conocido -eq 0 ]
    then
      SL_Conocido_Hook=$(echo ".")
      SL_Conocido_Hook_len=$(echo ".")
      SL_Conocido_Hook_iden=$(echo ".")
    else
      SL_Conocido_Hook=$(grep $hook "Analysis_Results/BLAST_Identificacion/"$gru"/"$spec"_hook_sl_conocidos.blastn" | awk -F "\t" -v min=$Min_SL_Hit '{if ($4>=min) print $2"\t"$12*10"\t"$4"\t"$3}' | sort -n -r -k 2 | head -n 1 | awk -F "\t" '{print $1}')
      SL_Conocido_Hook_len=$(grep $hook "Analysis_Results/BLAST_Identificacion/"$gru"/"$spec"_hook_sl_conocidos.blastn" | awk -F "\t" -v min=$Min_SL_Hit '{if ($4>=min) print $2"\t"$12*10"\t"$4"\t"$3}' | sort -n -r -k 2 | head -n 1 | awk -F "\t" '{print $3}')
      SL_Conocido_Hook_iden=$(grep $hook "Analysis_Results/BLAST_Identificacion/"$gru"/"$spec"_hook_sl_conocidos.blastn" | awk -F "\t" -v min=$Min_SL_Hit '{if ($4>=min) print $2"\t"$12*10"\t"$4"\t"$3}' | sort -n -r -k 2 | head -n 1 | awk -F "\t" '{print $4}')
    fi

    check_adaptador=$(grep $hook "Analysis_Results/BLAST_Identificacion/"$gru"/"$spec"_hook_adaptadores.blastn" | awk -F "\t" -v min=$Min_Adapt_Hit '{if ($4>=min) print}' | grep -c .)
    if [ $check_adaptador -eq 0 ]
    then
      Adaptador_Hook=$(echo ".")
      Adaptador_Hook_len=$(echo ".")
      Adaptador_Hook_iden=$(echo ".")
    else
      Adaptador_Hook=$(grep $hook "Analysis_Results/BLAST_Identificacion/"$gru"/"$spec"_hook_adaptadores.blastn" | awk -F "\t" -v min=$Min_Adapt_Hit '{if ($4>=min) print $2"\t"$12*10"\t"$4"\t"$3}' | sort -n -r -k 2 | head -n 1 | awk -F "\t" '{print $1}')
      Adaptador_Hook_len=$(grep $hook "Analysis_Results/BLAST_Identificacion/"$gru"/"$spec"_hook_adaptadores.blastn" | awk -F "\t" -v min=$Min_Adapt_Hit '{if ($4>=min) print $2"\t"$12*10"\t"$4"\t"$3}' | sort -n -r -k 2 | head -n 1 | awk -F "\t" '{print $3}')
      Adaptador_Hook_iden=$(grep $hook "Analysis_Results/BLAST_Identificacion/"$gru"/"$spec"_hook_adaptadores.blastn" | awk -F "\t" -v min=$Min_Adapt_Hit '{if ($4>=min) print $2"\t"$12*10"\t"$4"\t"$3}' | sort -n -r -k 2 | head -n 1 | awk -F "\t" '{print $4}')
    fi
    echo $spec" "$hook" "$SL_Conocido_Hook" "$SL_Conocido_Hook_len" "$SL_Conocido_Hook_iden" "$Adaptador_Hook" "$Adaptador_Hook_iden" "$Adaptador_Hook_len | tr " " "\t" >> "Analysis_Results/Summary_Hook_Anot.tab"
  done
}

coord_check () {
  if [ $Cor_start -lt 1 ]
  then
    print_sequence=FALSE
    Cor_start=1
  fi

  chr_end=$(grep -w $chr "Analysis_Results/Reference_Genome/"$reference_genome"_chr" | awk -F "\t" '{print $2}')

  if [ $Cor_end -gt $chr_end ]
  then
    print_sequence=FALSE
    Cor_end=$chr_end
  fi
  coord=$(echo $Cor_start"-"$Cor_end)
}

### Functions for the initial SL-RNA selection
gff_file_scan ()
{
  if [ -f $gff_file ]
  then

    grep ^$chr $gff_file | awk -F "\t" -v start=$Cor_start -v end=$Cor_end '{if (start <= $5  && $4 <= end) print $4"_"$5";;"$3";"$9}' | sed "s/ /_/g" | sort -u >> "Analysis_Results/pSL_Overlapping_Features/Temporales/GFF_Info.temp"
    scan_info=$(grep -c . "Analysis_Results/pSL_Overlapping_Features/Temporales/GFF_Info.temp")

    if [ $scan_info -gt 0 ]
    then
      conteo=1
      while [ $conteo -le $scan_info ]
      do
        gff_info=$(sed -n $conteo"p" "Analysis_Results/pSL_Overlapping_Features/Temporales/GFF_Info.temp" | sed "s/;;/ /")
        echo $filo" "$spec" "$loci_name" "$chr" "$Cor_start"_"$Cor_end" "$locus_cdhit" GFF "$gff_info | tr " " "\t" >> "Analysis_Results/pSL_Overlapping_Features/Matching_Info_"$Extencion".tab"
        conteo=$(($conteo + 1))
      done
    else
      echo $filo" "$spec" "$loci_name" "$chr" "$Cor_start"_"$Cor_end" "$locus_cdhit" GFF NA No_Data" | tr " " "\t" >> "Analysis_Results/pSL_Overlapping_Features/Matching_Info_"$Extencion".tab"
    fi
    rm "Analysis_Results/pSL_Overlapping_Features/Temporales/GFF_Info.temp"
  fi
}

RepeatMasker_scan ()
{
  if [ -f $repeat_masker_file ]
  then
    grep -w $chr $repeat_masker_file | grep -v Simple_repeat | awk -F "\t" -v start=$Cor_start -v end=$Cor_end '{if (start <= $7 && $6 <= end) print $6"_"$7" "$10";"$11}' >> "Analysis_Results/pSL_Overlapping_Features/Temporales/RepeatMasker_Info.temp"
    scan_info=$(grep -c . "Analysis_Results/pSL_Overlapping_Features/Temporales/RepeatMasker_Info.temp")
    if [ $scan_info -gt 0 ]
    then
      conteo=1
      while [ $conteo -le $scan_info ]
      do
        repeatmasker_basic_info=$(sed -n $conteo"p" "Analysis_Results/pSL_Overlapping_Features/Temporales/RepeatMasker_Info.temp")
        echo $filo" "$spec" "$loci_name" "$chr" "$Cor_start"_"$Cor_end" "$locus_cdhit" "$repeat_masker_type" "$repeatmasker_basic_info | tr " " "\t" >> "Analysis_Results/pSL_Overlapping_Features/Matching_Info_"$Extencion".tab"
        conteo=$(($conteo + 1))
      done
    else
      echo $filo" "$spec" "$loci_name" "$chr" "$Cor_start"_"$Cor_end" "$locus_cdhit" "$repeat_masker_type" NA No_Data" | tr " " "\t" >> "Analysis_Results/pSL_Overlapping_Features/Matching_Info_"$Extencion".tab"
    fi
    rm "Analysis_Results/pSL_Overlapping_Features/Temporales/RepeatMasker_Info.temp"
  fi
}

get_sec_filo_runs ()
{
  echo $filo" Recuperando secuencias"

  for filo_seq in $sequence_filo
  do
    check_origin_loci=$(echo $filo_seq | grep -c ^BLAST_New_)
    if [ $check_origin_loci -eq 0 ]
    then
      search_locus=$(echo $filo_seq | tr "_" "\n" | tail -n 1)
      search_spec=$(echo $filo_seq | sed "s/.*_Loci_//" | sed "s/_Locus.*//")
    else
      search_locus=$(echo $filo_seq | sed "s/_Loci_.*//" | sed "s/BLAST_New_//")
      search_spec=$(echo $filo_seq | sed "s/.*_Loci_//" | sed "s/_SLF_NA//")
    fi

    loci_type=clear
    rename_seq=$filo_seq

    check_cluster=$(grep $search_spec "Analysis_Results/pSL_Overlapping_Features/Matching_Info_"$Extencion".tab" | grep -w $search_locus | awk -F "\t" '{print $6}' | sort -u)
    if [ "$check_cluster" != "." ]
    then
      echo "|"$check_cluster"| :"$rename_seq
      rename_seq=$(echo $check_cluster"_"$rename_seq)
    fi

    check_transposon=$(grep $search_spec "Analysis_Results/pSL_Overlapping_Features/Matching_Info_"$Extencion".tab" | grep -w $search_locus | grep RepeatMasker_ | grep -v Unknown | grep -v No_Data | grep -c .)
    if [ $check_transposon -gt 0 ]
    then
      loci_type=transposon
      transposones=$(grep $search_spec "Analysis_Results/pSL_Overlapping_Features/Matching_Info_"$Extencion".tab" | grep -w $search_locus | grep RepeatMasker_ | grep -v Unknown | grep -v No_Data | awk -F "\t" '{print $9}' | awk -F ";" '{print $2}' | sort -u | tr "\n" "_" | sed "s/\//__/g" | sed "s/_$//")
      rename_seq=$(echo $rename_seq"_"$transposones)
    fi

    seqkit grep -w 0 -p $filo_seq Analysis_Results/MEME_Profile/$filo"_All_pSLs.work" | sed "s/>.*/>$rename_seq/" >> "Analysis_Results/MEME_Profile/Filogeny_Clustering/General/"$filo"_filo_sec_"$Extencion".fasta"

  done

  Filo_grupo=General
  seqkit grep -r -p $filo $Full_SLs | seqkit seq -w 0 --dna2rna >> "Analysis_Results/MEME_Profile/Filogeny_Clustering/"$Filo_grupo"/"$filo"_filo_sec_"$Extencion".fasta"

  check_seq=$(grep -c ">" "Analysis_Results/MEME_Profile/Filogeny_Clustering/"$Filo_grupo"/"$filo"_filo_sec_"$Extencion".fasta")

  if [ $check_seq -gt 2 ]
  then
    mafft --globalpair --quiet --maxiterate 1000 "Analysis_Results/MEME_Profile/Filogeny_Clustering/"$Filo_grupo"/"$filo"_filo_sec_"$Extencion".fasta" >> "Analysis_Results/MEME_Profile/Filogeny_Clustering/"$Filo_grupo"/"$filo"_filo_sec_"$Extencion".aln"
    $path_to_ninja/ninja "Analysis_Results/MEME_Profile/Filogeny_Clustering/"$Filo_grupo"/"$filo"_filo_sec_"$Extencion".aln" | grep . > "Analysis_Results/MEME_Profile/Filogeny_Clustering/"$Filo_grupo"/"$filo"_filo_sec_"$Extencion".nwk"
    TreeCluster.py -i "Analysis_Results/MEME_Profile/Filogeny_Clustering/"$Filo_grupo"/"$filo"_filo_sec_"$Extencion".nwk" -o "Analysis_Results/MEME_Profile/Filogeny_Clustering/"$Filo_grupo"/"$filo"_filo_sec_"$Extencion".clus" -t $TreeCluster_dist -m $TreeCluster_Mode

    meme_clusters=$(awk -F "\t" '{print $2}' "Analysis_Results/MEME_Profile/Filogeny_Clustering/"$Filo_grupo"/"$filo"_filo_sec_"$Extencion".clus" | grep -v ClusterNumber | sort | uniq -c | awk -v min=$MEME_Min_Cluster '{if ($1>=min) print $2}' | awk '{if ($1>0) print}')
    for mclus in $meme_clusters
    do
      meme_clus_seq=$(awk -F "\t" -v clus=$mclus '{if ($2==clus) print $1}' "Analysis_Results/MEME_Profile/Filogeny_Clustering/"$Filo_grupo"/"$filo"_filo_sec_"$Extencion".clus" )
      for mseq in $meme_clus_seq
      do
        echo $mclus": "$mseq
        seqkit grep -p $mseq "Analysis_Results/MEME_Profile/Filogeny_Clustering/"$Filo_grupo"/"$filo"_filo_sec_"$Extencion".fasta" >> "Analysis_Results/MEME_Profile/Filogeny_Clustering/"$Filo_grupo"/"$filo"_cluster_"$mclus"_"$Extencion".fasta"
      done
      meme "Analysis_Results/MEME_Profile/Filogeny_Clustering/"$Filo_grupo"/"$filo"_cluster_"$mclus"_"$Extencion".fasta" -rna -oc "Analysis_Results/MEME_Profile/Filogeny_Clustering/"$Filo_grupo"/"$filo"_cluster_"$mclus"_oops" -time 30000 -nostatus -mod oops -nmotifs $NMotifs -minw 5 -maxw 50
    done
  fi
}

cluster_check ()
{
  check_origin_clus=$(echo $loci_name | grep -c "KLoc-")
  if [ $check_origin_clus -eq 0 ]
  then
    search_loci_clus=$(echo $spec"_"$loci_name";")
  else
    search_loci_clus=$(echo "_New_"$loci_name"_Loci_")
  fi

  check_cluster=$(grep $spec "Analysis_Results/pSL_Overlapping_Features/Temporales/"$filo"_Cluster.temp" | grep $search_loci_clus | tr ";" "\n" | grep -c .)
  uggly_scan=FALSE
  if [ $check_cluster -gt 2 ]
  then
    uggly_scan=TRUE
    locus_cdhit=$(grep $spec "Analysis_Results/pSL_Overlapping_Features/Temporales/"$filo"_Cluster.temp" | grep $search_loci_clus | awk -F ";" '{print $1}')
  else
    locus_cdhit=$(echo ".")
  fi
}

################################################################################

if [ "$PREPARE_READS" == "TRUE" ]
then
  for gru in $grupos
  do
    especies=$(ls -d $gru"/"*"/" | sed "s/\/$//" | sed "s/.*\///")
    for spe in $especies
    do
      echo $gru": "$spe

      rm -r $gru"/"$spe"/raw_reads"
      rm -r $gru"/"$spe"/trimmed_reads"

      mkdir $gru"/"$spe"/raw_reads"
      mkdir $gru"/"$spe"/trimmed_reads"

      echo "#!/usr/bin/env bash" > $gru"/"$spe"/Run_trimmomatic.sh"

      check_truseq2=$(echo $use_TruSeq2 | grep -c $spe)
      if [ $check_truseq2 -gt 0 ]
      then
        adapter_file=$(echo $trimmomatic_installation"/adapters/TruSeq2-PE.fa:2")
      else
        adapter_file=$(echo $trimmomatic_installation"/adapters/TruSeq3-PE-2.fa")
      fi

      sra_files=$(ls $gru"/"$spe"/"*".sra" | sed "s/.*\///" | sed "s/.sra//")
      for read in $sra_files
      do
        echo $read
        fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files $gru"/"$spe"/"$read".sra" -O $gru"/"$spe"/raw_reads"
        echo "Comprimiendo"
        gzip $gru"/"$spe"/raw_reads/"*".fastq"
        echo "Prepare Trimming"
        echo "trimmomatic PE "$gru"/"$spe"/raw_reads/"$read'_1.fastq.gz' $gru"/"$spe"/raw_reads/"$read'_2.fastq.gz' -baseout $gru"/"$spe"/trimmed_reads/"$read'_trimmed_read.gz' -threads $threads "ILLUMINACLIP:"$adapter_file":2:30:10 SLIDINGWINDOW:5:20 MINLEN:25" >> $gru"/"$spe"/Run_trimmomatic.sh"
      done
      chmod +x $gru"/"$spe"/Run_trimmomatic.sh"
      $gru"/"$spe"/Run_trimmomatic.sh"

      fastqc -t $threads $gru"/"$spe"/raw_reads/"*
      fastqc -t $threads $gru"/"$spe"/trimmed_reads/"*
    done
  done
fi

################################################################################
################################################################################

if [ "$ASSEMBLY" == "TRUE" ]
then
  for gru in $grupos
  do
    especies=$(ls -d $gru"/"*"/" | sed "s/\/$//" | sed "s/.*\///")
    for spe in $especies
    do
      read_files=$(ls $gru"/"$spe"/trimmed_reads/"*"_trimmed_read_1P.gz"| sed 's/_trimmed_read_1P.gz//' | sed "s/.*\///")
      for read in $read_files
      do
        echo "Running: "$read
        Trinity --seqType fq --CPU $threads --left $gru"/"$spe"/trimmed_reads/"$read"_trimmed_read_1P.gz" --right $gru"/"$spe"/trimmed_reads/"$read"_trimmed_read_2P.gz" --no_normalize_reads --full_cleanup --output $read"_Trinity" --max_memory 30G
      done
    done
  done
fi

################################################################################
################################################################################

if [ "$RUN_SLFINDER1" == "TRUE" ]
then
  for gru in $grupos
  do
    especies=$(ls -d $gru"/"*"/" | sed "s/\/$//" | sed "s/.*\///")
    for spe in $especies
    do
      cd $gru"/"$spe

      rm -r Filtered_Assemblies
      rm -r SL-analysis

      ./SLFinder-step0 -a _Trinity.Trinity.fasta
      ./SLFinder-step1 -t $threads -me F -a _Trinity.Trinity.fasta
      ./SLFinder-step2 -t $threads

      cd $home_direcory
    done
  done
fi

################################################################################
################################################################################

if [ "$RUN_SLFINDER2" == "TRUE" ]
then
  for gru in $grupos
  do
    especies=$(ls -d $gru"/"*"/" | sed "s/\/$//" | sed "s/.*\///")
    for spe in $especies
    do
      cd $gru"/"$spe
      genome_file=$(ls *"genomic.fa")

      ./SLFinder-step3 -t $threads -g $genome_file -br $uniprot_reference

      cd $home_direcory
    done
  done
fi

################################################################################
################################################################################

if [ "$SLRNA_RETRIVAL" == "TRUE" ]
then
  # SLFinder_Classification.sh
  rm -r "Analysis_Results"

  mkdir "Analysis_Results"
  mkdir "Analysis_Results/Alineamientos_Hooks"
  mkdir "Analysis_Results/Alineamientos_pSL"
  mkdir "Analysis_Results/SLFinder_Runs"
  mkdir "Analysis_Results/SLFinder_Runs/Raw_Runs"
  mkdir "Analysis_Results/BLAST_Identificacion"
  mkdir "Analysis_Results/BLAST_Identificacion/Referencias"
  mkdir "Analysis_Results/Structure"
  mkdir "Analysis_Results/MEME_Profile"
  mkdir "Analysis_Results/Known_SLs"
  mkdir "Analysis_Results/Reference_Genome"

  makeblastdb -in $SLs_Conocidos_Fasta -dbtype nucl -parse_seqids -input_type fasta -out "Analysis_Results/BLAST_Identificacion/Referencias/SLs_Conocidos_Fasta_ref"
  makeblastdb -in $Adaptadores_Ilumina_Fasta -dbtype nucl -parse_seqids -input_type fasta -out "Analysis_Results/BLAST_Identificacion/Referencias/Adaptadores_Ilumina_Fasta_ref"

  echo "Species Group Assemblies Run_Status Total_Hooks Hook_Variants SL_Classes Unclear_SL_Classes Total_Loci Clear Unclear Other Hit_SL_conocidos Hit_Adaptadores Genome_File" | tr " " "\t" >> "Analysis_Results/Summary_SLFinder.tab"
  echo "Especie Hook SL_Reportado SL_Hit_Len SL_Iden Adaptador Adaptador_Iden Adaptador_Hit_Len" | tr " " "\t" >> "Analysis_Results/Summary_Hook_Anot.tab"
  echo "Filo Especie Loci_BLAST SL_Ref Ident Start End Bitscore Sense Chr Similar_SLs SLFinder_Locus" | tr " " "\t"  >> "Analysis_Results/Known_SLs/Registro_Locis_SLs_conocidos.tab"

  echo "Filos seleccionados: "$grupos

  for gru in $grupos
  do
    run_meme=FALSE
    echo "Processing: "$gru
    mkdir "Analysis_Results/Alineamientos_Hooks/"$gru
    mkdir "Analysis_Results/Alineamientos_pSL/"$gru
    mkdir "Analysis_Results/BLAST_Identificacion/"$gru
    mkdir "Analysis_Results/Structure/"$gru
    mkdir "Analysis_Results/Known_SLs/"$gru
    mkdir "Analysis_Results/Known_SLs/"$gru"/BLAST_Search"
    mkdir "Analysis_Results/Known_SLs/"$gru"/BLAST_Search/References"
    mkdir "Analysis_Results/Known_SLs/"$gru"/Temporales"
    mkdir "Analysis_Results/Known_SLs/"$gru"/Potential_Locus"

    echo "Species Chr Start-End Strand Loci_ID" | tr " " "\t" >> "Analysis_Results/Alineamientos_pSL/"$gru"/Problematic_SLs.coords"
    echo "" > "Analysis_Results/Alineamientos_pSL/"$gru"/Problematic_SLs.fasta"

    species=$(ls -d $gru"/"*"/" | sed "s/\/$//" | sed "s/.*\///")

    for spec in $species
    do
      echo "Species Chr Start-End Strand Loci_ID" | tr " " "\t" >> "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs_Truncado.coord"
      echo "" >> "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs_Truncado.fasta"

      cp -r $gru"/"$spec"/SL-analysis" "Analysis_Results/SLFinder_Runs/"$spec"_Results"
      rm -r "Analysis_Results/SLFinder_Runs/"$spec"_Results/temp_files"

      Summary_Files=$(ls $gru"/"$spec"/SL-analysis/"*"_Summary" | grep -c .)
      assemblies=$(ls $gru"/"$spec"/"*"_Trinity.Trinity.fasta" | grep -c .)
      missing_hook=NA
      reference_genome=NA
      Total_hooks=NA
      total_loci=NA
      Work_Hook_Variants=NA
      sl_classes=NA
      unclear_sl_classes=NA
      clear=NA
      unclear=NA
      other=NA

      echo $gru" "$spec" "$Summary_Files" "$assemblies

      if [ $Summary_Files -le 1 ]
      then
        SL_Status=No_Candidates
      elif [ $Summary_Files -eq 2 ]
      then
        SL_Status=Unverified
        get_hook_variants
      elif [ $Summary_Files -eq 3 ]
      then
        run_meme=TRUE
        reference_genome=$(grep "Genome Reference file: "  $gru"/"$spec"/SL-analysis/Step3_Summary" | sed "s/Genome Reference file: //")
        cp $gru"/"$spec"/"$reference_genome "Analysis_Results/Reference_Genome/"$reference_genome
        seqkit fx2tab -i -l -n  $gru"/"$spec"/"$reference_genome >> "Analysis_Results/Reference_Genome/"$reference_genome"_chr"

        SL_Status=Candidates
        get_hook_variants

        #########################
        #### Check pSL Locis ####
        #########################

        echo "Recuperando SL Locis"
        recorrer_loci=$(grep -c . $gru"/"$spec"/SL-analysis/Results/Step3_Loci.tab")
        count=2
        clear=0
        unclear=0
        other=0
        total_loci=$(($recorrer_loci - 1))

        while [ $count -le $recorrer_loci ]
        do
          temp_loci_name=$(sed -n $count"p" $gru"/"$spec"/SL-analysis/Results/Step3_Loci.tab" | awk -F "\t" '{print $1}')
          loci_name=$(echo "SLFinder_Loci_"$spec"_"$temp_loci_name)

          sense=$(sed -n $count"p" $gru"/"$spec"/SL-analysis/Results/Step3_Loci.tab" | awk -F "\t" '{print $4}')
          type=$(sed -n $count"p" $gru"/"$spec"/SL-analysis/Results/Step3_Loci.tab" | awk -F "\t" '{print $6}' | grep -c "prima")
          chr=$(sed -n $count"p" $gru"/"$spec"/SL-analysis/Results/Step3_Loci.tab" | awk -F "\t" '{print $2}')

          # Conteos info
          hook_hit=$(sed -n $count"p" $gru"/"$spec"/SL-analysis/Results/Step3_Loci.tab" | awk -F "\t" '{print $7}' | awk -F "-" '{print $1"-"$2}')
          check_adaptador=$(grep $spec "Analysis_Results/Summary_Hook_Anot.tab" | grep -w $hook_hit | awk -F "\t" '{print $6}' | grep -c -v \\.)
          adaptador_print=$(grep -w $hook_hit "Analysis_Results/BLAST_Identificacion/"$gru"/"$spec"_hook_adaptadores.blastn" | awk -F "\t" -v min=$Min_Adapt_Hit '{if ($4>=min) print $2}' | sort -u | tr "\n" ";" | sed "s/;$//")

          if [ $check_adaptador -gt 0 ]
          then
            echo "Adaptador encontrado"
            echo $spec" "$temp_loci_name" "$grupo" "$chr" "$adaptador_print >> "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLp_Loci_Descatados.tab"
          else
            if [ $type -gt 0 ]
            then
              grupo=Donor_Site
              clear=$(($clear + 1))
            else
              type=$(sed -n $count"p" $gru"/"$spec"/SL-analysis/Results/Step3_Loci.tab" | awk -F "\t" '{print $6}' | grep -c "Unclear")
              if [ $type -gt 0 ]
              then
                grupo=Unclear
                unclear=$(($unclear + 1))
              else
                grupo=Other
                other=$(($other + 1))
              fi
            fi

            # Conteos SL Classes
            if [ -f $gru"/"$spec"/SL-analysis/Results/SL_Clases.fa" ]
            then
              sl_classes=$(grep -c ">" $gru"/"$spec"/SL-analysis/Results/SL_Clases.fa")
            else
              sl_classes=0
            fi
            if [ -f $gru"/"$spec"/SL-analysis/Results/SL_Clases-Unclear.fa" ]
            then
              unclear_sl_classes=$(grep -c ">" $gru"/"$spec"/SL-analysis/Results/SL_Clases-Unclear.fa")
            else
              unclear_sl_classes=0
            fi

            print_sequence=TRUE
            # Extraer secuencias SLs
            if [ "$sense" == "Forward" ]
            then
              Cor_start=$(sed -n $count"p" $gru"/"$spec"/SL-analysis/Results/Step3_Loci.tab" | awk -F "\t" '{print $8-20}')
              Cor_end=$(sed -n $count"p" $gru"/"$spec"/SL-analysis/Results/Step3_Loci.tab" | awk -F "\t" '{print $9+100}')
              coord_check

              if [ "$print_sequence" == "TRUE" ]
              then
                printf "%s %s %s %s\n%s %s %s\n" $chr" "$coord" plus" | blastdbcmd -db $gru"/"$spec"/SL-analysis/temp_files/Reference/Reference_genome-ref" -entry_batch - | sed "s/>.*/>$loci_name/" | seqkit seq -w 0 --dna2rna >> "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs_"$grupo".fasta"
                printf "%s %s %s %s\n%s %s %s\n" $chr" "$coord" plus" | blastdbcmd -db $gru"/"$spec"/SL-analysis/temp_files/Reference/Reference_genome-ref" -entry_batch - | sed "s/>.*/>$loci_name/" | seqkit seq -w 0 --dna2rna >> "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs_todos.fasta"
              else
                printf "%s %s %s %s\n%s %s %s\n" $chr" "$coord" plus" | blastdbcmd -db $gru"/"$spec"/SL-analysis/temp_files/Reference/Reference_genome-ref" -entry_batch - | sed "s/>.*/>$loci_name/" | seqkit seq -w 0 --dna2rna >> "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs_Truncado.fasta"
                echo $spec" "$chr" "$Cor_start"-"$Cor_end" plus "$loci_name >> "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs_Truncado.coord"
              fi

              print_coord=$(echo $Cor_start"-"$Cor_end)
              echo $spec" "$temp_loci_name" "$grupo" "$chr" "$print_coord" plus" | tr " " "\t" >> "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs.coords"
            elif [ "$sense" == "Reverse" ]
            then
              Cor_start=$(sed -n $count"p" $gru"/"$spec"/SL-analysis/Results/Step3_Loci.tab" | awk -F "\t" '{print $8-100}')
              Cor_end=$(sed -n $count"p" $gru"/"$spec"/SL-analysis/Results/Step3_Loci.tab" | awk -F "\t" '{print $9+20}')
              coord_check

              if [ "$print_sequence" == "TRUE" ]
              then
                printf "%s %s %s %s\n%s %s %s\n" $chr" "$coord" minus" | blastdbcmd -db $gru"/"$spec"/SL-analysis/temp_files/Reference/Reference_genome-ref" -entry_batch - | sed "s/>.*/>$loci_name/" | seqkit seq -w 0 --dna2rna >> "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs_"$grupo".fasta"
                printf "%s %s %s %s\n%s %s %s\n" $chr" "$coord" minus" | blastdbcmd -db $gru"/"$spec"/SL-analysis/temp_files/Reference/Reference_genome-ref" -entry_batch - | sed "s/>.*/>$loci_name/" | seqkit seq -w 0 --dna2rna >> "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs_todos.fasta"
              else
                printf "%s %s %s %s\n%s %s %s\n" $chr" "$coord" minus" | blastdbcmd -db $gru"/"$spec"/SL-analysis/temp_files/Reference/Reference_genome-ref" -entry_batch - | sed "s/>.*/>$loci_name/" | seqkit seq -w 0 --dna2rna >> "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs_Truncado.fasta"
                echo $spec" "$chr" "$Cor_start"-"$Cor_end" minus "$loci_name | tr " " "\t" >> "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs_Truncado.coord"
              fi

              print_coord=$(echo $Cor_start"-"$Cor_end)
              echo $spec" "$temp_loci_name" "$grupo" "$chr" "$print_coord" minus" | tr " " "\t" >> "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs.coords"
            else
              Cor_start=$(sed -n $count"p" $gru"/"$spec"/SL-analysis/Results/Step3_Loci.tab" | awk -F "\t" '{print $8-100}')
              Cor_end=$(sed -n $count"p" $gru"/"$spec"/SL-analysis/Results/Step3_Loci.tab" | awk -F "\t" '{print $9+100}')

              coord_check
              printf "%s %s %s %s\n%s %s %s\n" $chr" "$coord" plus" | blastdbcmd -db $gru"/"$spec"/SL-analysis/temp_files/Reference/Reference_genome-ref" -entry_batch - | sed "s/>.*/>$loci_name/" | seqkit seq -w 0 --dna2rna >> "Analysis_Results/Alineamientos_pSL/"$gru"/Problematic_SLs.fasta"
              echo $spec" "$chr" "$Cor_start" "$Cor_end" "$loci_name | tr " " "\t" >> "Analysis_Results/Alineamientos_pSL/"$gru"/Problematic_SLs.coords"
              echo $spec" "$chr" "$Cor_start"-"$Cor_end" plus "$loci_name | tr " " "\t" >> "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs_Truncado.coord"

              print_coord=$(echo $Cor_start"-"$Cor_end)
              echo $spec" "$temp_loci_name" "$grupo" "$chr" "$print_coord" plus" | tr " " "\t" >> "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs.coords"
            fi
          fi
          count=$(($count + 1))
        done

        #########################
        #### Known SL Search ####
        #########################

        makeblastdb -in "Analysis_Results/Reference_Genome/"$reference_genome -dbtype nucl -parse_seqids -input_type fasta -out "Analysis_Results/Known_SLs/"$gru"/BLAST_Search/References/"$spec"_blast-ref"
        blastn -task blastn-short -query $SLs_Conocidos_Fasta -db "Analysis_Results/Known_SLs/"$gru"/BLAST_Search/References/"$spec"_blast-ref" -outfmt "6 std sstrand" -out "Analysis_Results/Known_SLs/"$gru"/BLAST_Search/"$spec"_raw.blastn" -num_threads $threads -perc_identity 95 -ungapped -qcov_hsp_perc 75

        check_known_blast_results=$(grep -c . "Analysis_Results/Known_SLs/"$gru"/BLAST_Search/"$spec"_raw.blastn")
        if [ $check_known_blast_results -eq 0 ]
        then
          N_known_SL_Hits=0
        else
          known_sls_chr=$(awk -F "\t" '{print $2}' "Analysis_Results/Known_SLs/"$gru"/BLAST_Search/"$spec"_raw.blastn" | sort -u)
          klocus=0
          for kchr in $known_sls_chr
          do
            echo "Comprobando SLs Conocidos: "$kchr
            grep -w $kchr "Analysis_Results/Known_SLs/"$gru"/BLAST_Search/"$spec"_raw.blastn" | sort -n -k 9 | awk -F "\t" '{print $1"\t"$3"\t"$9"\t"$10"\t"$12*10"\t"$13}' >> "Analysis_Results/Known_SLs/"$gru"/Temporales/"$kchr"_Potential"
            klocus=$(($klocus + 1))
            kcount=2
            khits=$(grep -c -w $kchr "Analysis_Results/Known_SLs/"$gru"/BLAST_Search/"$spec"_raw.blastn")

            main_sl=$(sed -n 1p "Analysis_Results/Known_SLs/"$gru"/Temporales/"$kchr"_Potential" | awk -F "\t" '{print $1}')
            iden_main_sl=$(sed -n 1p "Analysis_Results/Known_SLs/"$gru"/Temporales/"$kchr"_Potential" | awk -F "\t" '{print $2}')
            kstart=$(sed -n 1p "Analysis_Results/Known_SLs/"$gru"/Temporales/"$kchr"_Potential" | awk -F "\t" '{print $3"\n"$4}' | sort -n | head -n 1)
            kend=$(sed -n 1p "Analysis_Results/Known_SLs/"$gru"/Temporales/"$kchr"_Potential" | awk -F "\t" '{print $3"\n"$4}' | sort -n | tail -n 1)
            main_sl_bitscore=$(sed -n 1p "Analysis_Results/Known_SLs/"$gru"/Temporales/"$kchr"_Potential" | awk -F "\t" '{print $5}')
            main_sl_strand=$(sed -n 1p "Analysis_Results/Known_SLs/"$gru"/Temporales/"$kchr"_Potential" | awk -F "\t" '{if ($6=="minus") print "Reverse"; else print "Forward"}' )
            echo "KLoc-"$klocus" "$main_sl" "$iden_main_sl" "$kstart" "$kend" "$main_sl_bitscore" "$main_sl_strand" "$kchr | tr " " "\t" >> "Analysis_Results/Known_SLs/"$gru"/Potential_Locus/"$spec"_pot_locus.hits"

            while [ $kcount -le $khits ]
            do
              main_sl=$(sed -n $kcount"p" "Analysis_Results/Known_SLs/"$gru"/Temporales/"$kchr"_Potential" | awk -F "\t" '{print $1}')
              iden_main_sl=$(sed -n $kcount"p" "Analysis_Results/Known_SLs/"$gru"/Temporales/"$kchr"_Potential" | awk -F "\t" '{print $2}')
              main_sl_bitscore=$(sed -n $kcount"p" "Analysis_Results/Known_SLs/"$gru"/Temporales/"$kchr"_Potential" | awk -F "\t" '{print $5}')
              main_sl_strand=$(sed -n $kcount"p" "Analysis_Results/Known_SLs/"$gru"/Temporales/"$kchr"_Potential" | awk -F "\t" '{if ($6=="minus") print "Reverse"; else print "Forward"}' )

              test_kstart=$(sed -n $kcount"p" "Analysis_Results/Known_SLs/"$gru"/Temporales/"$kchr"_Potential" | awk -F "\t" '{print $3"\n"$4}' | sort -n | head -n 1)
              test_kend=$(sed -n $kcount"p" "Analysis_Results/Known_SLs/"$gru"/Temporales/"$kchr"_Potential" | awk -F "\t" '{print $3"\n"$4}' | sort -n | tail -n 1)

              if ! [ $kstart -le $test_kend -a $test_kstart -le $kend ]
              then
                kstart=$test_kstart
                kend=$test_kend
                klocus=$(($klocus + 1))
              fi

              echo "KLoc-"$klocus" "$main_sl" "$iden_main_sl" "$test_kstart" "$test_kend" "$main_sl_bitscore" "$main_sl_strand" "$kchr | tr " " "\t" >> "Analysis_Results/Known_SLs/"$gru"/Potential_Locus/"$spec"_pot_locus.hits"
              kcount=$(($kcount + 1))
            done
          done

          klocu_num=$(awk -F "\t" '{print $1}' "Analysis_Results/Known_SLs/"$gru"/Potential_Locus/"$spec"_pot_locus.hits" | sort -u)
          for kl in $klocu_num
          do
            best_bitscore=$(grep -w $kl "Analysis_Results/Known_SLs/"$gru"/Potential_Locus/"$spec"_pot_locus.hits" | awk -F "\t" '{print $6}' | sort -n | tail -n 1)
            ksl_hits=$(grep -w $kl "Analysis_Results/Known_SLs/"$gru"/Potential_Locus/"$spec"_pot_locus.hits" | awk -F "\t" -v bit=$best_bitscore '{if ($6=bit) print $2}' | tr "\n" "_" | sed "s/_$//")

            grep -w $kl "Analysis_Results/Known_SLs/"$gru"/Potential_Locus/"$spec"_pot_locus.hits" | awk -F "\t" -v bit=$best_bitscore -v sl_hits=$ksl_hits '{if ($6=bit) print $1" "$2" "$3" "$4" "$5" "$6/10" "$7" "$8" "sl_hits}' | head -n 1 | tr " " "\t" >> "Analysis_Results/Known_SLs/"$gru"/Potential_Locus/"$spec"_pot_locus.tab"
            echo $kl" -- "$best_bitscore" -- "$ksl_hits
          done

          # Determinar que Loci recuperé de novo y como nombrarlos

          awk -F "\t" '{print $2" "$4" "$5}' "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs.coords" | tr "-" " " | tr " " "\t"  | sed "s/Locus\t/Locus-/" | sed "s/F_gigantica\t/F_gigantica-/" >> "Analysis_Results/Known_SLs/"$gru"/Temporales/"$spec"_SLFinder_Loci_coords.temp"
          for kl in $klocu_num
          do
            main_sl_kl=$(grep -w $kl "Analysis_Results/Known_SLs/"$gru"/Potential_Locus/"$spec"_pot_locus.tab" | awk -F "\t" '{print $1"_"$2}')
            start_kl=$(grep -w $kl "Analysis_Results/Known_SLs/"$gru"/Potential_Locus/"$spec"_pot_locus.tab" | awk -F "\t" '{print $4}')
            sentido_kl=$(grep -w $kl "Analysis_Results/Known_SLs/"$gru"/Potential_Locus/"$spec"_pot_locus.tab" | awk -F "\t" '{print $7}')
            end_kl=$(grep -w $kl "Analysis_Results/Known_SLs/"$gru"/Potential_Locus/"$spec"_pot_locus.tab" | awk -F "\t" '{print $5}')
            chr_kl=$(grep -w $kl "Analysis_Results/Known_SLs/"$gru"/Potential_Locus/"$spec"_pot_locus.tab" | awk -F "\t" '{print $8}')

            check_SLFinder_loci=$(awk -F "\t" -v chr=$chr_kl '{ if ($2==chr) print}' "Analysis_Results/Known_SLs/"$gru"/Temporales/"$spec"_SLFinder_Loci_coords.temp" | awk -F "\t" -v start=$start_kl -v end=$end_kl '{if (start <= $4 && $3 <= end) print $1}' | grep -c .)

            if [ $check_SLFinder_loci -eq 0 ]
            then
              grep -w $kl "Analysis_Results/Known_SLs/"$gru"/Potential_Locus/"$spec"_pot_locus.tab" | awk -F "\t" -v filo=$gru -v spec=$spec '{print filo"\t"spec"\t"$0"\tNew"}' >> "Analysis_Results/Known_SLs/Registro_Locis_SLs_conocidos.tab"

              print_sequence=TRUE
              modificar=$(echo 'BLAST_New_'$main_sl_kl'_Loci_'$spec'_SLF_NA')
              chr=$chr_kl

              if [ "$sentido_kl" == "Forward" ]
              then
                echo "BLAST:"$main_sl_kl" --- "$sentido_kl
                Cor_start=$(echo $start_kl | awk '{print $1-20}')
                Cor_end=$(echo $end_kl | awk '{print $1+100}')
                coord_check

                if [ "$print_sequence" == "TRUE" ]
                then
                  printf "%s %s %s %s\n%s %s %s\n" $chr" "$coord" plus" | blastdbcmd -db $gru"/"$spec"/SL-analysis/temp_files/Reference/Reference_genome-ref" -entry_batch - | sed "s/>.*/>$modificar/" | seqkit seq -w 0 --dna2rna >> "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs_todos.fasta"
                else
                  printf "%s %s %s %s\n%s %s %s\n" $chr" "$coord" plus" | blastdbcmd -db $gru"/"$spec"/SL-analysis/temp_files/Reference/Reference_genome-ref" -entry_batch - | sed "s/>.*/>$modificar/" | seqkit seq -w 0 --dna2rna >> "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs_Truncado.fasta"
                fi
                print_coord=$(echo $Cor_start"-"$Cor_end)
                echo $spec" "$main_sl_kl" Nuevos_BLAST "$chr" "$print_coord" plus" | tr " " "\t" >> "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs.coords"
              else
                echo "BLAST:"$main_sl_kl" --- "$sentido_kl
                Cor_start=$(echo $start_kl | awk '{print $1-100}')
                Cor_end=$(echo $end_kl | awk '{print $1+20}')
                coord_check

                if [ "$print_sequence" == "TRUE" ]
                then
                  printf "%s %s %s %s\n%s %s %s\n" $chr" "$coord" minus" | blastdbcmd -db $gru"/"$spec"/SL-analysis/temp_files/Reference/Reference_genome-ref" -entry_batch - | sed "s/>.*/>$modificar/" | seqkit seq -w 0 --dna2rna >> "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs_todos.fasta"
                else
                  printf "%s %s %s %s\n%s %s %s\n" $chr" "$coord" minus" | blastdbcmd -db $gru"/"$spec"/SL-analysis/temp_files/Reference/Reference_genome-ref" -entry_batch - | sed "s/>.*/>$modificar/" | seqkit seq -w 0 --dna2rna >> "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs_Truncado.fasta"
                fi
                print_coord=$(echo $Cor_start"-"$Cor_end)
                echo $spec" "$main_sl_kl" Nuevos_BLAST "$chr" "$print_coord" minus" | tr " " "\t" >> "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs.coords"
              fi
            else
              SLFinder_loci=$(awk -F "\t" -v chr=$chr_kl '{ if ($2==chr) print}' "Analysis_Results/Known_SLs/"$gru"/Temporales/"$spec"_SLFinder_Loci_coords.temp" | awk -F "\t" -v start=$start_kl -v end=$end_kl '{if (start <= $4 && $3 <= end) print $1}')
              for slf_loci in $SLFinder_loci
              do
                cambiar=$(echo 'SLFinder_Loci_'$spec'_'$slf_loci)
                modificar=$(echo 'BLAST_SLF_'$main_sl_kl'_Loci_'$spec'_'$slf_loci)

                sed -i "/$cambiar$/s/>$cambiar/>$modificar/"  "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs_todos.fasta"
                if [ -f "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs_Truncado.fasta" ]
                then
                  sed -i "/$cambiar$/s/>$cambiar/>$modificar/" "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs_Truncado.fasta"
                fi
                sed -i "/$slf_loci\t/s/$/\tConocido/" "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs.coords"
                grep -w $kl "Analysis_Results/Known_SLs/"$gru"/Potential_Locus/"$spec"_pot_locus.tab" | awk -F "\t" -v filo=$gru -v spec=$spec -v slfinder=$slf_loci '{print filo"\t"spec"\t"$0"\t"slfinder}' >> "Analysis_Results/Known_SLs/Registro_Locis_SLs_conocidos.tab"
              done
              if [ $check_SLFinder_loci -gt 1 ]
              then
                echo $spec" - "$main_sl_kl" ($chr_kl:$start_kl/$end_kl) SLFinder: "$SLFinder_loci >> "Analysis_Results/Known_SLs/Problemas.tab"
              fi
            fi
          done
          grep -v Conocido "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs.coords" | grep -v Nuevos_BLAST >> "Analysis_Results/Known_SLs/Registro_Locis_SLs_Exclusivos_SLFinder.tab"
        fi

        ###########################
        #### SM site locations ####
        ###########################
        seqkit locate -m 2 -P -f $SM_like_sites "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs_todos.fasta" >> "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs_sm.coord"
        sm_loci=$(seqkit seq -n "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs_todos.fasta")

        for loc in $sm_loci
        do
          seqkit grep -p $loc "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs_todos.fasta" | seqkit fx2tab >> "Analysis_Results/Structure/"$gru"/"$spec"_estructura.tab"
          check_sm_hit=$(grep -c -w $loc "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs_sm.coord")
          seq_len=$(seqkit grep -p $loc "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs_todos.fasta" | seqkit fx2tab -n -l | awk -F "\t" '{print $2}')
          sm_count=1

          echo $loc" --- "$check_sm_hit
          if [ $check_sm_hit -eq 0 ]
          then
            seq_estructura=
            print_base=1
            while [ $print_base -le $seq_len ]
            do
              seq_estructura=$(echo $seq_estructura".")
              print_base=$(($print_base + 1))
            done
            echo "No_SM "$seq_estructura | tr " " "\t" >> "Analysis_Results/Structure/"$gru"/"$spec"_estructura.tab"
          else
            grep -w $loc "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs_sm.coord" | awk -F "\t" '{print $5"\t"$6"\t"$7}' | sort -u | awk -F "\t" '{print "SM_Site-"NR"\t"$1"\t"$2"\t"$3}' >> "Analysis_Results/Structure/SM_site.temp"
            check_sm_hit=$(grep -c . "Analysis_Results/Structure/SM_site.temp")

            while [ $sm_count -le $check_sm_hit ]
            do
              seq_estructura=
              print_base=1
              sm_site=$(sed -n $sm_count"p" "Analysis_Results/Structure/SM_site.temp" | awk -F "\t" '{print $1}')
              sm_start=$(sed -n $sm_count"p" "Analysis_Results/Structure/SM_site.temp" | awk -F "\t" '{print $2}')
              sm_end=$(sed -n $sm_count"p" "Analysis_Results/Structure/SM_site.temp" | awk -F "\t" '{print $3}')

              while [ $print_base -le $seq_len ]
              do
                if [ $print_base -ge $sm_start ] && [ $print_base -le $sm_end ]
                then
                  seq_estructura=$(echo $seq_estructura"x")
                else
                  seq_estructura=$(echo $seq_estructura".")
                fi
                print_base=$(($print_base + 1))
              done
              echo $sm_site"_"$sm_start"-"$sm_end " "$seq_estructura | tr " " "\t" >> "Analysis_Results/Structure/"$gru"/"$spec"_estructura.tab"
              sm_count=$(($sm_count + 1))
            done
            rm "Analysis_Results/Structure/SM_site.temp"
          fi
        done
        cat "Analysis_Results/Alineamientos_pSL/"$gru"/"$spec"_SLs_todos.fasta" >> "Analysis_Results/MEME_Profile/"$gru"_All_pSLs.fasta"
      fi
      echo $spec" "$gru" "$assemblies" "$SL_Status" "$Total_hooks" "$Work_Hook_Variants" "$sl_classes" "$unclear_sl_classes" "$total_loci" "$clear" "$unclear" "$other" "$hit_sl_hook" "$hit_adaptadores_hook" "$reference_genome | tr " " "\t" >> "Analysis_Results/Summary_SLFinder.tab"
    done

    if [ $run_meme == TRUE ]
    then
      echo "comprobando MEME"
      echo "Corrida General: "$gru
      seqkit grep -v -s -r -p NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN "Analysis_Results/MEME_Profile/"$gru"_All_pSLs.fasta" >> "Analysis_Results/MEME_Profile/"$gru"_All_pSLs.no_missing.fasta"
      cd-hit-est -i "Analysis_Results/MEME_Profile/"$gru"_All_pSLs.no_missing.fasta" -d 0 -o "Analysis_Results/MEME_Profile/"$gru"_All_pSLs.work"
      meme "Analysis_Results/MEME_Profile/"$gru"_All_pSLs.work" -rna -oc "Analysis_Results/MEME_Profile/"$gru"_All_Results_oops" -time 30000 -nostatus -mod oops -nmotifs $NMotifs -minw 5 -maxw 50
    fi

  done
fi

################################################################################
################################################################################

if [ "$RUN_REPEATMASKER" == "TRUE" ]
then
  rm -r "Analysis_Results/Info_Transposones"
  mkdir "Analysis_Results/Info_Transposones"
  mkdir "Analysis_Results/Info_Transposones/RepeatMasker_Runs"
  mkdir "Analysis_Results/Info_Transposones/RepeatMasker_Custom_Runs"
  mkdir "Analysis_Results/Info_Transposones/RepeatModeler_Runs"
  mkdir "Analysis_Results/Info_Transposones/References"
  mkdir "Analysis_Results/Info_Transposones/chr_info"

  for gru in $grupos
  do
    especies=$(ls -d $gru"/"*"/" | sed "s/\/$//" | sed "s/.*\///")
    for spe in $especies
    do
      echo "First masker -- "$spe
      genome_file=$(ls $gru"/"$spe"/"*"genomic.fa"  | sed "s/.*\///")
      RepeatMasker -no_is $gru"/"$spe"/"$genome_file -dir "Analysis_Results/Info_Transposones/RepeatMasker_Runs/"$spe
      cat "Analysis_Results/Info_Transposones/RepeatMasker_Runs/"$spe"/"$genome_file".out" | tr -s " " | sed "s/^ //" | tr " " "\t" >> "Analysis_Results/Info_Transposones/RepeatMasker_Runs/"$spe"/"$genome_file".work"

      echo "Modeler -- "$spe
      BuildDatabase $gru"/"$spe"/"$genome_file -name "Analysis_Results/Info_Transposones/References/"$spe"_Ref"
      cd Analysis_Results/Info_Transposones/References
      RepeatModeler -database $spe"_Ref"
      rep_mod_run=$(ls -d RM_*)
      mv $rep_mod_run "../RepeatModeler_Runs/"$spe"_RepMod"
      cd ../../../

      echo "Second RepeatMasker -- "$spe
      RepeatMasker -lib "Analysis_Results/Info_Transposones/RepeatModeler_Runs/"$spe"_RepMod/consensi.fa.classified" -no_is $gru"/"$spe"/"$genome_file -dir "Analysis_Results/Info_Transposones/RepeatMasker_Custom_Runs/"$spe
      cat "Analysis_Results/Info_Transposones/RepeatMasker_Custom_Runs/"$spe"/"$genome_file".out" | tr -s " " | sed "s/^ //" | tr " " "\t" >> "Analysis_Results/Info_Transposones/RepeatMasker_Custom_Runs/"$spe"/"$genome_file".work"
    done
  done
fi

################################################################################
################################################################################

if [ "$SLRNA_INITIAL_FILTERING" == "TRUE" ]
then
  rm -r "Analysis_Results/pSL_Overlapping_Features"
  rm -r "Analysis_Results/MEME_Profile/Filogeny_Clustering"

  mkdir "Analysis_Results/pSL_Overlapping_Features"
  mkdir "Analysis_Results/pSL_Overlapping_Features/Temporales"
  mkdir "Analysis_Results/MEME_Profile/Filogeny_Clustering"
  mkdir "Analysis_Results/MEME_Profile/Filogeny_Clustering/General"

  echo "Filo Especie Loci Chr Coord_pSL Cluster Fuente Coord_Ref Info" | tr " " "\t" >> "Analysis_Results/pSL_Overlapping_Features/Matching_Info_Overlap.tab"
  echo "Filo Especie Loci Chr Coord_pSL Cluster Fuente Coord_Ref Info" | tr " " "\t" >> "Analysis_Results/pSL_Overlapping_Features/Matching_Info_Region.tab"

  for filo in $grupos
  do
    if [ -f "Analysis_Results/MEME_Profile/"$filo"_All_pSLs.work.clstr" ]
    then
      sed "s/\.\.\..*$//" "Analysis_Results/MEME_Profile/"$filo"_All_pSLs.work.clstr" | sed "s/.*nt, >//" | sed "s/$/&;/" | tr -d "\n" | tr ">" "\n" | tr " " "_" | grep . >> "Analysis_Results/pSL_Overlapping_Features/Temporales/"$filo"_Cluster.temp"
    else
      echo "" >> "Analysis_Results/pSL_Overlapping_Features/Temporales/"$filo"_Cluster.temp"
    fi

    especies=$(grep $filo "Analysis_Results/Summary_SLFinder.tab" | awk -F "\t" '{print $1"\t"$15}' | grep -v -w "NA" | awk -F "\t" '{print $1}')

    for spec in $especies
    do
      mkdir "Analysis_Results/pSL_Overlapping_Features/"$spec
      genome=$(grep $spec "Analysis_Results/Summary_SLFinder.tab" | awk -F "\t" '{print $15}')
      gff_file=$(ls $filo"/"$spec"/"*"annotations.gtf")
      seqkit fx2tab -i -l -n "Analysis_Results/Reference_Genome/"$genome > "Analysis_Results/Reference_Genome/"$genome"_chr"
      echo $filo" "$spec" "$genome

      psl_data=$(grep -c . "Analysis_Results/Alineamientos_pSL/"$filo"/"$spec"_SLs.coords")
      count=1

      while [ $count -le $psl_data ]
      do
        Extencion=Overlap
        loci_name=$(sed -n $count"p" "Analysis_Results/Alineamientos_pSL/"$filo"/"$spec"_SLs.coords" | awk -F "\t" '{print $2}')
        Cor_start=$(sed -n $count"p" "Analysis_Results/Alineamientos_pSL/"$filo"/"$spec"_SLs.coords" | awk -F "\t" '{print $5}' | awk -F "-" '{print $1}')
        Cor_end=$(sed -n $count"p" "Analysis_Results/Alineamientos_pSL/"$filo"/"$spec"_SLs.coords" | awk -F "\t" '{print $5}' | awk -F "-" '{print $2}')
        chr=$(sed -n $count"p" "Analysis_Results/Alineamientos_pSL/"$filo"/"$spec"_SLs.coords" | awk -F "\t" '{print $4}')

        echo $loci_name" -- "$chr": "$Cor_start" - "$Cor_end
        check_origin_clus=$(echo $loci_name | grep -c "KLoc-")

        if [ $check_origin_clus -eq 0 ]
        then
          search_loci_clus=$(echo $spec"_"$loci_name";")
        else
          search_loci_clus=$(echo "_New_"$loci_name"_Loci_")
        fi

        check_cluster=$(grep $spec "Analysis_Results/pSL_Overlapping_Features/Temporales/"$filo"_Cluster.temp" | grep $search_loci_clus | tr ";" "\n" | grep -c .)
        uggly_scan=FALSE
        if [ $check_cluster -gt 2 ]
        then
          uggly_scan=TRUE
          locus_cdhit=$(grep $spec "Analysis_Results/pSL_Overlapping_Features/Temporales/"$filo"_Cluster.temp" | grep $search_loci_clus | awk -F ";" '{print $1}')
        else
          locus_cdhit=$(echo ".")
        fi

        cluster_check # traer
        ######## Adjust cluser names
        check_origin_clus=$(echo $loci_name | grep -c "KLoc-")
        if [ $check_origin_clus -eq 0 ]
        then
          search_loci_clus=$(echo $spec"_"$loci_name";")
        else
          search_loci_clus=$(echo "_New_"$loci_name"_Loci_")
        fi
        check_cluster=$(grep $spec "Analysis_Results/pSL_Overlapping_Features/Temporales/"$filo"_Cluster.temp" | grep $search_loci_clus | tr ";" "\n" | grep -c .)
        uggly_scan=FALSE
        if [ $check_cluster -gt 2 ]
        then
          uggly_scan=TRUE
          locus_cdhit=$(grep $spec "Analysis_Results/pSL_Overlapping_Features/Temporales/"$filo"_Cluster.temp" | grep $search_loci_clus | awk -F ";" '{print $1}')
        else
          locus_cdhit=$(echo ".")
        fi

        gff_file_scan
        repeat_masker_file="Analysis_Results/Info_Transposones/RepeatMasker_Runs/"$spec"/"$genome".work"
        repeat_masker_type=RepeatMasker_Basic
        RepeatMasker_scan
        repeat_masker_file="Analysis_Results/Info_Transposones/RepeatMasker_Custom_Runs/"$spec"/"$genome".work"
        repeat_masker_type=RepeatMasker_Custom
        RepeatMasker_scan

        Cor_start=$(($Cor_start - $Rango_Region))
        Cor_end=$(($Cor_end + $Rango_Region))

        ####################################################################################
        # Verify if coordinates remain within chr
        if [ $Cor_start -lt 1 ]
        then
          Cor_start=1
        fi
        chr_end=$(grep -w $chr "Analysis_Results/Reference_Genome/"$genome"_chr" | awk -F "\t" '{print $2}')
        if [ $Cor_end -gt $chr_end ]
        then
          Cor_end=$chr_end
        fi
        ####################################################################################

        Extencion=Region
        gff_file_scan
        repeat_masker_file="Analysis_Results/Info_Transposones/RepeatMasker_Runs/"$spec"/"$genome".work"
        repeat_masker_type=RepeatMasker_Basic
        RepeatMasker_scan
        repeat_masker_file="Analysis_Results/Info_Transposones/RepeatMasker_Custom_Runs/"$spec"/"$genome".work"
        repeat_masker_type=RepeatMasker_Custom
        RepeatMasker_scan
        count=$(($count + 1))
      done

      # Extract repeat family
      Extencion=Overlap
      repeat_masker_type=RepeatMasker_Custom
      rmd_family=$(grep $repeat_masker_type "Analysis_Results/pSL_Overlapping_Features/Matching_Info_"$Extencion".tab" | grep $spec | awk -F "\t" '{print $9}' | grep -v No_Data | grep -v Unknown | awk -F ";" '{print $1}' | sort -u)
      for rmd in $rmd_family
      do
        seqkit grep -r -p $rmd"#" "Analysis_Results/Info_Transposones/RepeatModeler_Runs/"$spec"_RepMod/consensi.fa.classified" >> "Analysis_Results/pSL_Overlapping_Features/"$spec"/"$rmd"_matching_loci.fasta"
        rmd_loci=$(grep $repeat_masker_type "Analysis_Results/pSL_Overlapping_Features/Matching_Info_"$Extencion".tab" | grep $spec | grep $rmd | awk -F "\t" '{print $3}' | sort -u)
        for rloc in $rmd_loci
        do
          seqkit grep -r -p $rloc "Analysis_Results/MEME_Profile/"$filo"_All_pSLs.fasta" "Analysis_Results/Info_Transposones/RepeatModeler_Runs/"$spec"_RepMod/consensi.fa.classified" | seqkit seq --rna2dna >> "Analysis_Results/pSL_Overlapping_Features/"$spec"/"$rmd"_matching_loci.fasta"
        done
      done
    done

    sequence_filo=$(seqkit seq -n "Analysis_Results/MEME_Profile/"$filo"_All_pSLs.work")
    Extencion=Overlap
    get_sec_filo_runs

    Extencion=Region
    get_sec_filo_runs
  done
fi

################################################################################
################################################################################
