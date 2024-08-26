#!/usr/bin/env bash

# Once identified the SL leaders, it is necesary to define SL tags. Short
# portions of the leader sequence that are going to be used to identify the
# SL bearing reads. Once done this scripts produce some initial exploratory
# counts and classification, and starts the identification of putative chimeras.

### Variables
grupos=$(echo "Cestodes Trematoda")
threads=24
Skip_UTR_pre_classification=$(echo "Sparganum_proliferum") # Due to differences in format, the UTR localization for this species was ignored in the initial serch.
cestoda_sl_tags=/home/amanda/PROYECTS/SLFinder_Publication/Selected_Sequences_Work/Cestodes_SL_Tags.fasta
trematopda_sl_tags=/home/amanda/PROYECTS/SLFinder_Publication/Selected_Sequences_Work/Trematoda_SL_Tags.fasta
qcov=40 # Query coverage for the annotation of potential chimeric gene models

############################################################################################

## Controle Modules to run
RUN_SLFINDER_GENES=FALSE
RUN_ACE_PRE_CLASSIFICATION=FALSE
RUN_ACE_INTERNAL_SEC=TRUE
RUN_BLAST_INTERNAL_STIES=FALSE
RUN_CLASSIFY_INTERNAL_STIES=FALSE

############################################################################################

### Functions
gene_basic_info () {
  if [ $check_wormbase -gt 0 ]
  then
    chr=$(grep -F 'gene:'$gen'";' $gtf_file | grep -F -w transcript | awk -F "\t" '{print $1}' | sort -u)
    strand=$(grep -F 'gene:'$gen'";' $gtf_file | grep -F -w CDS | awk -F "\t" '{print $7}' | sort -u)
  else
    chr=$(grep -F $gen'";' $gtf_file | grep -F -w transcript | awk -F "\t" '{print $1}' | sort -u)
    strand=$(grep  -F $gen'";' $gtf_file | grep -F -w CDS | awk -F "\t" '{print $7}' | sort -u)
  fi
}

reads_counts () {
  echo "Conteo reads "$spe
  echo "Gen_ID Trans Chromosome Strand Aceptor Concordance Status Exon_ID "$archivos" Total" | tr " " "\t" >> "../Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts/"$spe"_"$sl"_read.counts"
  for gen_ace in $genes_aces
  do
    gen=$(echo $gen_ace | sed "s/|||.*//")
    ace=$(echo $gen_ace | sed "s/.*|||//")
    gene_basic_info

    if [ $check_wormbase -gt 0 ]
    then
      gene_trans=$(grep -F 'gene_id "gene:'$gen'";' $gtf_file | sed "s/.*transcript://" | sed 's/".*//' | sort -u)
    else
      gene_trans=$gen
    fi
    for trans in $gene_trans
    do
      ace_clasificator
      print_conteo=$(echo $gen" "$trans" "$chr" "$strand" "$ace" "$ace_concordance" "$status" "$exon)
      total=0
      for file in $archivos
      do
        conteo=$(grep -F -w $gen $spe"/SL_Genes/Results/Potential_aceptor_sites/"$file"_"$sl"_aceptor_sites.tab" | awk -F "\t" -v ace=$ace '{if ($2 == "Y" && $3 == ace) print}' | grep -c .)
        total=$(($total + $conteo))

        echo $sl" "$conteo" "$file" "$total | tr " " "\t" >> "../Analysis_Results/SL_Acceptor_Sites/Temporal/"$spe"_SL_read_counts"
        print_conteo=$(echo $print_conteo" "$conteo)
      done
      echo $print_conteo" "$total | tr " " "\t" >> "../Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts/"$spe"_"$sl"_read.counts"
    done
  done
  echo "Completado Conteo "$spe
}

check_non_coding () {
  check_coding_gene=$(grep -c . "../Analysis_Results/SL_Acceptor_Sites/Temporal/Trans_Exon_temp_"$trans)
  if [ $check_coding_gene -eq 0 ]
  then
    if [ "$strand" == "+" ]
    then
      if [ $check_wormbase -gt 0 ]
      then
        grep -F 'transcript_id "transcript:'$trans'";' $gtf_file  | awk -F "\t" '{if ($3=="exon") print $4"\t"$5}' | grep . | sort -n -k 1 | awk -F "\t" '{print $0"\tEx-"NR}' >> "../Analysis_Results/SL_Acceptor_Sites/Temporal/Trans_Exon_temp_"$trans
      else
        grep $trans'";' $gtf_file | awk -F "\t" '{if ($3=="exon") print $4"\t"$5}' | grep . | sort -n -k 1 | awk -F "\t" '{print $0"\tEx-"NR}' >> "../Analysis_Results/SL_Acceptor_Sites/Temporal/Trans_Exon_temp_"$trans
      fi
    else
      if [ $check_wormbase -gt 0 ]
      then
        grep -F 'transcript_id "transcript:'$trans'";' $gtf_file | awk -F "\t" '{if ($3=="exon") print $5"\t"$4}' | grep . | sort -r -n -k 1 | awk -F "\t" '{print $0"\tEx-"NR}' >> "../Analysis_Results/SL_Acceptor_Sites/Temporal/Trans_Exon_temp_"$trans
      else
        grep $trans'";' $gtf_file  | awk -F "\t" '{if ($3=="exon") print $5"\t"$4}' | grep . | sort -r -n -k 1 | awk -F "\t" '{print $0"\tEx-"NR}' >> "../Analysis_Results/SL_Acceptor_Sites/Temporal/Trans_Exon_temp_"$trans
      fi
    fi
  fi
}

complete_sl_sum () {
  echo "Suma Total "$spe

  check_files_to_count=$(ls "../Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts/"$spe"_"*"_read.counts" | grep -c .)
  if [ $check_files_to_count -eq 1 ]
  then
    genes_ace_suma=$(grep -F -v Gen_ID "../Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts/"$spe"_"*"_read.counts" | awk -F "\t" '{print $1"|||"$5}' | sort -u)
  else
    genes_ace_suma=$(grep -F -v Gen_ID "../Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts/"$spe"_"*"_read.counts" | awk -F ":" '{print $2}' | awk -F "\t" '{print $1"|||"$5}' | sort -u)
  fi

  archivos_2=$(ls $spe"/SL_Genes/Results/Potential_aceptor_sites/"*"_aceptor_sites.tab" | sed "s/.*\///g" | awk -F "_" '{print $1}' | sort -u)
  echo "Gen_ID Aceptor Chromosome Strand Concordance Status Exon_ID "$archivos_2 | tr " " "\t" >> "../Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts/"$spe"_tier3.counts"
  echo "Gen_ID Aceptor Chromosome Strand Concordance Status Exon_ID "$archivos_2 | tr " " "\t" >> "../Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts/"$spe"_tier10.counts"
  echo "Gen_ID Aceptor Chromosome Strand Concordance Status Exon_ID "$archivos_2" Total Tier_3_reads Tier_10_reads" | tr " " "\t" >> "../Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts/"$spe"_total.counts"

  echo "Conteo final: "
  for gen_ace2 in $genes_ace_suma
  do
    gen=$(echo $gen_ace2 | sed "s/|||.*//")
    ace=$(echo $gen_ace2 | sed "s/.*|||//")

    gene_basic_info
    site_concord=$(grep -F -w $gen "../Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts/"$spe"_"*"_read.counts" | grep -F -w $ace | awk -F "\t" '{print $6}' | sort -u | tr "\n" "_"| sed "s/_$//")
    site_status=$(grep -F -w $gen "../Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts/"$spe"_"*"_read.counts" | grep -F -w $ace | awk -F "\t" '{print $7}' | sort -u | tr "\n" "_"| sed "s/_$//")
    site_exon=$(grep -F -w $gen "../Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts/"$spe"_"*"_read.counts" | grep -F -w $ace | awk -F "\t" '{print $8}' | sort -u | tr "\n" "_"| sed "s/_$//")

    print_conteo=$(echo $gen" "$ace" "$chr" "$strand" "$site_concord" "$site_status" "$site_exon)

    registro_conteo=0
    for file2 in $archivos_2
    do
      conteo=$(grep -F -w $gen $spe"/SL_Genes/Results/Potential_aceptor_sites/"$file2"_"*"_aceptor_sites.tab" | awk -F "\t" -v ace=$ace '{if ($2 == "Y" && $3 == ace) print}' | grep -c .)

      print_conteo=$(echo $print_conteo" "$conteo)

      registro_conteo=$(($registro_conteo + $conteo))
    done

    ### AD OC FIX!!!!!
    if [ "$spe" == "Fasciolopsis_buski" ]
    then
      registro_conteo=$(($registro_conteo / 2))
    fi

    print_conteo=$(echo $print_conteo" "$registro_conteo)

    if [ $registro_conteo -gt 3 ]
    then
      tier_read=$(echo "X")
      echo $print_conteo | tr " " "\t" >> "../Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts/"$spe"_tier3.counts"
    else
      tier_read=$(echo ".")
    fi

    if [ $registro_conteo -gt 10 ]
    then
      tier_read=$(echo $tier_read" X")
      echo $print_conteo | tr " " "\t" >> "../Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts/"$spe"_tier10.counts"
    else
      tier_read=$(echo $tier_read" .")
    fi

    echo $print_conteo" "$tier_read | tr " " "\t" >> "../Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts/"$spe"_total.counts"
  done
  echo "Completado "$spe
}

summary_types () {
  sumary_genes_aces=$(grep -F -v Gen_ID $file_count | awk -F "\t" '{print $1"|||"$2}' | sort -u)
  num_sumary_genes_aces=$(grep -F -v Gen_ID $file_count | awk -F "\t" '{print $1"|||"$2}' | sort -u | grep -c .)

  total_isoform=0
  total_discordante=0
  total_end_of_line=0
  total_start=0
  total_utr=0
  total_gen_end=0
  total_ex=0
  total_intron_can=0
  total_intron_alt=0

  for sum_ace in $sumary_genes_aces
  do
    gen=$(echo $sum_ace | sed "s/|||.*//")
    ace=$(echo $sum_ace | sed "s/.*|||//")

    check_isoform_shenanigans=$(grep -F -w $gen "../Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts/"$spe"_"*"_read.counts" | grep -F -w $ace | awk -F "\t" '{print $6"_"$7}' | sort -u | grep -c .)
    if [ $check_isoform_shenanigans -gt 1 ]
    then
      total_isoform=$(($total_isoform + 1))
    else
      check_end_of_line=$(grep -F -w $gen "../Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts/"$spe"_"*"_read.counts" | grep -F -w $ace | awk -F "\t" '{print $6}' | grep -c "coting_limit")
      if [ $check_end_of_line -gt 0 ]
      then
        total_end_of_line=$(($total_end_of_line + 1))
      else
        check_discord=$(grep -F -w $gen "../Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts/"$spe"_"*"_read.counts" | grep -F -w $ace |  awk -F "\t" '{print $6}' | grep -c "Inverso")
        if [ $check_discord -gt 0 ]
        then
          total_discordante=$(($total_discordante + 1))
        else
          count_start=$(grep -F -w $gen "../Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts/"$spe"_"*"_read.counts" | grep -w "Concordante" | grep -F -w $ace | awk -F "\t" '{print $7}' | sort -u | grep -c "Trans_Start")
          count_utr=$(grep -F -w $gen "../Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts/"$spe"_"*"_read.counts" | grep -w "Concordante" | grep -F -w $ace | awk -F "\t" '{print $7}' | sort -u | grep -c "Trans_UTR")
          count_end=$(grep -F -w $gen "../Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts/"$spe"_"*"_read.counts" | grep -w "Concordante" | grep -F -w $ace | awk -F "\t" '{print $7}' | sort -u | grep -c "Trans_End")
          count_ex=$(grep -F -w $gen "../Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts/"$spe"_"*"_read.counts" | grep -w "Concordante" | grep -F -w $ace | awk -F "\t" '{print $7}' | sort -u | grep -c "Exon")
          count_inron_can=$(grep -F -w $gen "../Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts/"$spe"_"*"_read.counts" | grep -w "Concordante" | grep -F -w $ace | awk -F "\t" '{print $7}' | sort -u | grep -c "Intron_Canonical")
          count_intron_alt=$(grep -F -w $gen "../Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts/"$spe"_"*"_read.counts" | grep -w "Concordante" | grep -F -w $ace | awk -F "\t" '{print $7}' | sort -u | grep -c "Intron_Alt")

          total_start=$(($total_start + $count_start))
          total_utr=$(($total_utr + $count_utr))
          total_gen_end=$(($total_gen_end + $count_end))
          total_ex=$(($total_ex + $count_ex))
          total_intron_can=$(($total_intron_can + $count_inron_can))
          total_intron_alt=$(($total_intron_alt + $count_intron_alt))
        fi
      fi
    fi
  done
  echo $spe" "$num_sumary_genes_aces" "$total_isoform" "$total_end_of_line" "$total_discordante" "$total_start" "$total_utr" "$total_gen_end" "$total_ex" "$total_intron_can" "$total_intron_alt | tr " " "\t" >> "../Analysis_Results/SL_Acceptor_Sites/"$filo"_Summary_"$sumary_tag"_tipos.tab"
}

summary_sl_tags_info ()
{
  echo "Resumen Tags "$spe
  print_info=$(echo $spe" "$ensamblados_todos)
  for sl_filo in $sl_tags_filo
  do
    check_valid_hits=$(ls $spe"/SL_Genes/Results/Potential_aceptor_sites/"*"_aceptor_sites.tab" | grep -F -c $sl_filo)

    if [ $check_valid_hits -gt 0 ]
    then
      SL_tag_gen_N=$(awk -F "\t" '{if ($2=="Y") print $1}' $spe"/SL_Genes/Results/Potential_aceptor_sites/"*"_"$sl_filo"_aceptor_sites.tab" | sort -u | grep -c .)
      if [ $SL_tag_gen_N -gt 0 ]
      then
        SL_tag_read_N=$(grep -F -w $sl_filo "../Analysis_Results/SL_Acceptor_Sites/Temporal/"$spe"_SL_read_counts" | awk '{ sum += $2; n++ } END { if (n > 0) print sum; }')
        SL_tag_info=$(echo $SL_tag_gen_N"|"$SL_tag_read_N)
      else
        SL_tag_info=0
      fi
    else
      SL_tag_info=NA
    fi
    print_info=$(echo $print_info" "$SL_tag_info)
  done

  # Sitios no canonicos

  no_canon_genes=$(awk -F "\t" '{if ($2=="N") print $1}' $spe"/SL_Genes/Results/Potential_aceptor_sites/"*"_aceptor_sites.tab" | sort -u | grep -c .)
  no_canon_reads=$(awk -F "\t" '{if ($2=="N") print $1}' $spe"/SL_Genes/Results/Potential_aceptor_sites/"*"_aceptor_sites.tab" | grep -c .)
  no_canon_info=$(echo $no_canon_genes"|"$no_canon_reads)

  couting_limit_genes=$(awk -F "\t" '{if ($2=="Coting_Boundary") print $1}' $spe"/SL_Genes/Results/Potential_aceptor_sites/"*"_aceptor_sites.tab" | sort -u | grep -c .)
  couting_limit_reads=$(awk -F "\t" '{if ($2=="Coting_Boundary") print $1}' $spe"/SL_Genes/Results/Potential_aceptor_sites/"*"_aceptor_sites.tab" | grep -c .)
  coting_limit_info=$(echo $couting_limit_genes"|"$couting_limit_reads)

  echo $print_info" "$no_canon_info" "$coting_limit_info | tr " " "\t" >> "../Analysis_Results/SL_Acceptor_Sites/"$filo"_Summary_hits.tab"
  echo "Done Hit Sumary!!!"

  echo "Read Summary"

  total_reads=$(seqkit stats -T $spe"/"*"1P.gz" | grep -F -v -w num_seqs | awk '{ sum += $4; n++ } END { if (n > 0) print sum; }' )
  raw_sl_reads=$(seqkit stats -T $spe"/SL_Genes/Results/read_SL/"*"_reads_filtered_1.fq" | grep -F -v -w num_seqs | awk '{ sum += $4; n++ } END { if (n > 0) print sum; }' )
  mapped_reads=$(awk -F "\t" '{if ($2=="Y") print $1}' $spe"/SL_Genes/Results/Potential_aceptor_sites/"*"_aceptor_sites.tab" | grep -c .)

  raw_percentage_total=$(echo $raw_sl_reads" "$total_reads | awk -F " " '{printf "%.4f\n", 100*$1/$2}')
  map_percentage_total=$(echo $mapped_reads" "$total_reads | awk -F " " '{printf "%.4f\n", 100*$1/$2}')

  echo $spe" "$ensamblados_todos" "$total_reads" "$raw_sl_reads"("$raw_percentage_total"%) "$mapped_reads"("$map_percentage_total"%)"

  print_info=$(echo $spe" "$ensamblados_todos" "$total_reads" "$raw_sl_reads"("$raw_percentage_total"%) "$mapped_reads"("$map_percentage_total"%)")

  for sl_filo in $sl_tags_filo
  do
    check_valid_hits2=$(ls $spe"/SL_Genes/Results/Potential_aceptor_sites/"*"_aceptor_sites.tab" | grep -c $sl_filo)
    if [ $check_valid_hits2 -gt 0 ]
    then
      sl_reads=$(awk -F "\t" '{if ($2=="Y") print $1}' $spe"/SL_Genes/Results/Potential_aceptor_sites/"*"_"$sl_filo"_aceptor_sites.tab" | grep -c .)
      if [ $sl_reads -eq 0 ]
      then
        print_info=$(echo $print_info" 0(0%)")
      else
        per_sl_reads=$(echo $sl_reads" "$mapped_reads | awk -F " " '{printf "%.4f\n", 100*$1/$2}')
        print_info=$(echo $print_info" "$sl_reads"("$per_sl_reads"%)")
      fi
    else
      print_info=$(echo $print_info" 0(NA)")
    fi
  done
  echo $print_info | tr " " "\t" >> "../Analysis_Results/SL_Acceptor_Sites/"$filo"_Summary_reads.tab"
  echo "Read Summary Done"

  sumary_tag=All
  file_count="../Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts/"$spe"_total.counts"
  summary_types

  sumary_tag=tier_3
  file_count="../Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts/"$spe"_tier3.counts"
  summary_types

  sumary_tag=tier_10
  file_count="../Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts/"$spe"_tier10.counts"
  summary_types
  echo "Tipo de ace Sumary"
}

coord_check () {
  if [ $coord -lt 1 ]
  then
    print=false
  else
    chr_end=$(seqkit grep -p $chr $genome_fasta | seqkit fx2tab -n -i -l | awk -F "\t" '{print $2}')

    if [ $coord -gt $chr_end ]
    then
      print=false
    fi
  fi
}

start_classifier () {
  check_banned=$(echo $Skip_UTR_pre_classification | grep -c $spe)
  if [ $check_banned -gt 0 ]
  then
    status=Trans_Start
    exon=Start
  else
    if [ $check_wormbase -gt 0 ]
    then
      UTR=$(grep -F 'transcript_id "transcript:'$trans'";' $gtf_file | awk -F "\t" '{if ($3=="exon") print $4"\t"$5}' | grep . | sort -n -k 1 | awk -F "\t" -v ace=$ace '{if ($1<=ace && $2>=ace) print}' | grep -c .)
    else
      UTR=$(grep -F $trans'";' $gtf_file | awk -F "\t" '{if ($3=="exon") print $4"\t"$5}' | grep . | sort -n -k 1 | awk -F "\t" -v ace=$ace '{if ($1<=ace && $2>=ace) print}' | grep -c .)
    fi

    if [ $UTR -eq 0 ]
    then
      status=Trans_Start
      exon=Start
    else
      status=Trans_UTR
      exon=UTR
    fi
  fi
}

ace_clasificator () {
  if [ "$strand" == "+" ]
  then
    if [ $check_wormbase -gt 0 ]
    then
      grep -F 'transcript_id "transcript:'$trans'";' $gtf_file  | awk -F "\t" '{if ($3=="CDS") print $4"\t"$5}' | grep . | sort -n -k 1 | awk -F "\t" '{print $0"\tEx-"NR}' >> "../Analysis_Results/SL_Acceptor_Sites/Temporal/Trans_Exon_temp_"$trans
    else
      grep -F $trans'";' $gtf_file | awk -F "\t" '{if ($3=="CDS") print $4"\t"$5}' | grep . | sort -n -k 1 | awk -F "\t" '{print $0"\tEx-"NR}' >> "../Analysis_Results/SL_Acceptor_Sites/Temporal/Trans_Exon_temp_"$trans
    fi

    print=true
    coord=$(echo $ace | awk '{print $1-1}')
    coord_check
    if [ "$print" == "true" ]
    then
      ace_concordance=$(seqkit subseq -r $coord":"$ace --chr $chr --gtf $gtf_file $genome_fasta | grep -c "AG" | awk '{if ($1==1) print "Concordante"; else print "Inverso"}' )
    else
      ace_concordance=coting_limit
    fi

    gene_start=$(head -n 1 "../Analysis_Results/SL_Acceptor_Sites/Temporal/Trans_Exon_temp_"$trans | awk -F "\t" '{print $1}')
    gene_end=$(tail -n 1 "../Analysis_Results/SL_Acceptor_Sites/Temporal/Trans_Exon_temp_"$trans | awk -F "\t" '{print $2}')

    if [ $ace -gt $gene_end ] # 1) aceptor por fuera del gen
    then
      status=Trans_End
      exon=NA
    elif [ $ace -lt $gene_start ] # 2) aceptor delante del gen
    then
      start_classifier
    else # 3) Aceptor en medio del gen
      recorrer=$(grep -c . "../Analysis_Results/SL_Acceptor_Sites/Temporal/Trans_Exon_temp_"$trans)
      count=1
      while [ $count -le $recorrer ]
      do
        start_ex=$(sed -n $count"p" "../Analysis_Results/SL_Acceptor_Sites/Temporal/Trans_Exon_temp_"$trans | awk -F "\t" '{print $1}')
        end_ex=$(sed -n $count"p" "../Analysis_Results/SL_Acceptor_Sites/Temporal/Trans_Exon_temp_"$trans | awk -F "\t" '{print $2}')
        id_ex=$(sed -n $count"p" "../Analysis_Results/SL_Acceptor_Sites/Temporal/Trans_Exon_temp_"$trans | awk -F "\t" '{print $3}')

        if [ $ace -gt $end_ex ] # El SL está asociado a otro exon
        then
          count=$(($count + 1))
        elif [ $ace -ge $start_ex ] && [ $ace -le $end_ex ] # El SL está en este exon
        then
          status=Exon
          exon=$id_ex
          count=$(($recorrer + $recorrer))
        elif [ $ace -lt $start_ex ] # El SL está en un intron
        then
          distance=$(($start_ex - $ace))
          prev=$(($count - 1))
          prev_ex=$(sed -n $prev"p" "../Analysis_Results/SL_Acceptor_Sites/Temporal/Trans_Exon_temp_"$trans | awk -F "\t" '{print $3}')
          count=$(($recorrer + $recorrer))

          exon=$(echo $prev_ex"_"$id_ex)
          if [ $distance -eq 1 ]
          then
            status=Intron_Canonical
          else
            status=$(echo "Intron_Alt-"$distance)
          fi
        else
          echo "BANANA"
        fi
      done
    fi
  else
    if [ $check_wormbase -gt 0 ]
    then
      grep -F 'transcript_id "transcript:'$trans'";' $gtf_file | awk -F "\t" '{if ($3=="CDS") print $5"\t"$4}' | grep . | sort -r -n -k 1 | awk -F "\t" '{print $0"\tEx-"NR}' >> "../Analysis_Results/SL_Acceptor_Sites/Temporal/Trans_Exon_temp_"$trans
    else
      grep -F $trans'";' $gtf_file  | awk -F "\t" '{if ($3=="CDS") print $5"\t"$4}' | grep . | sort -r -n -k 1 | awk -F "\t" '{print $0"\tEx-"NR}' >> "../Analysis_Results/SL_Acceptor_Sites/Temporal/Trans_Exon_temp_"$trans
    fi
    check_non_coding

    print=true
    coord=$(echo $ace | awk '{print $1+1}')
    coord_check
    if [ "$print" == "true" ]
    then
      ace_concordance=$(seqkit subseq -r $ace":"$coord --chr $chr --gtf $gtf_file $genome_fasta | grep -c "CT" | awk '{if ($1==1) print "Concordante"; else print "Inverso"}' )
    else
      ace_concordance=coting_limit
    fi

    gene_start=$(head -n 1 "../Analysis_Results/SL_Acceptor_Sites/Temporal/Trans_Exon_temp_"$trans | awk -F "\t" '{print $1}')
    gene_end=$(tail -n 1 "../Analysis_Results/SL_Acceptor_Sites/Temporal/Trans_Exon_temp_"$trans | awk -F "\t" '{print $2}')
    echo $trans" "$ace": "$gene_start"_"$gene_end

    if [ $ace -lt $gene_end ] # 1) aceptor por fuera del gen
    then
      status=Trans_End
      exon=NA
    elif [ $ace -gt $gene_start ] # 2) aceptor delante del gen
    then
      start_classifier
    else # 3) Aceptor en medio del gen
      recorrer=$(grep -c . "../Analysis_Results/SL_Acceptor_Sites/Temporal/Trans_Exon_temp_"$trans)
      count=1
      while [ $count -le $recorrer ]
      do
        start_ex=$(sed -n $count"p" "../Analysis_Results/SL_Acceptor_Sites/Temporal/Trans_Exon_temp_"$trans | awk -F "\t" '{print $1}')
        end_ex=$(sed -n $count"p" "../Analysis_Results/SL_Acceptor_Sites/Temporal/Trans_Exon_temp_"$trans | awk -F "\t" '{print $2}')
        id_ex=$(sed -n $count"p" "../Analysis_Results/SL_Acceptor_Sites/Temporal/Trans_Exon_temp_"$trans | awk -F "\t" '{print $3}')

        if [ $ace -lt $end_ex ] # El SL está asociado a otro exon
        then
          count=$(($count + 1))
        elif [ $ace -le $start_ex ] && [ $ace -ge $end_ex ] # El SL está en este exon
        then
          status=Exon
          exon=$id_ex
          count=$(($recorrer + $recorrer))
        elif [ $ace -gt $start_ex ] # El SL está en este exon
        then
          distance=$(($ace - $start_ex))
          prev=$(($count - 1))
          prev_ex=$(sed -n $prev"p" "../Analysis_Results/SL_Acceptor_Sites/Temporal/Trans_Exon_temp_"$trans | awk -F "\t" '{print $3}')
          count=$(($recorrer + $recorrer))

          exon=$(echo $prev_ex"_"$id_ex)

          if [ $distance -eq 1 ]
          then
            status=Intron_Canonical
          else
            status=$(echo "Intron_Alt-"$distance)
          fi
        else
          echo "BANANA"
        fi
      done
    fi
  fi
  rm "../Analysis_Results/SL_Acceptor_Sites/Temporal/Trans_Exon_temp_"$trans
}

extract_info () {
  transcripts=$(grep -w $gen $work_directory"/SL_Read_Counts/"$spe"_"*"_read.counts" | grep $get_info | sed "s/.*://" | awk -F "\t" '{print $2}' | sort -u)
  chr=$(grep -w $gen $work_directory"/SL_Read_Counts/"$spe"_total.counts" | awk -F "\t" '{print $3}' | sort -u)
  strand=$(grep -w $gen $work_directory"/SL_Read_Counts/"$spe"_total.counts" | awk -F "\t" '{print $4}' | sort -u)
  echo $gen" "$chr" "$strand
}

get_tier_data () {
  count=2
  recorrer=$(grep -c . $work_directory"/SL_Read_Counts/"$spe"_"$tier".counts")

  while [ $count -le $recorrer ]
  do
    gen=$(sed -n $count"p" $work_directory"/SL_Read_Counts/"$spe"_"$tier".counts" | awk -F "\t" '{print $1}')
    ace=$(sed -n $count"p" $work_directory"/SL_Read_Counts/"$spe"_"$tier".counts" | awk -F "\t" '{print $2}')

    transcripts=$(grep -w -F $gen $work_directory"/SL_Read_Counts/"$spe"_"*"_read.counts" | awk -F "\t" '{print $2}' | sort -u)
    for trans in $transcripts
    do
      grep -F $trans"_"$ace"_" $work_directory"/Inserciones_Internas/"$spe"_inserciones_internas.tab" >> $work_directory"/Inserciones_Internas/"$spe"_inserciones_internas_"$tier".tab"
      grep -F $trans"_" $work_directory"/Inserciones_Internas/"$spe"_intrones.tab" >> $work_directory"/Inserciones_Internas/"$spe"_intrones_"$tier".tab"
    done
    count=$(($count + 1))
  done
}

separator () {
  recorrer_exones=$(grep -c . $file_exon_data)
  count=1
  mitad=A

  while [ $count -le $recorrer_exones ]
  do
    exon_turno=$(sed -n $count"p" $file_exon_data | awk -F "\t" '{print $3}')
    coord=$(sed -n $count"p" $file_exon_data | awk -F "\t" '{print $1"-"$2" "$4}')

    if [ "$exon_turno" == "$local_ace_exon" ]
    then
      mitad=D
    fi
    exon_name=$(echo $trans"_"$ace"_"$exon_turno"_Mitad-"$mitad)
    printf "%s %s %s %s\n%s %s %s\n" $chr" "$coord | blastdbcmd -db $reference_blast"/"$spe"_ref" -entry_batch - | sed "s/>.*/>$exon_name/" | seqkit seq -w 0 >> $work_directory"/Inserciones_Internas/Registro_Exones/"$spe"/"$trans"_"$ace"_"$tipo"_"$mitad".fasta"

    count=$(($count + 1))
  done
  seq_name=$(echo $trans"_"$ace"_"$local_ace_exon"_"$tipo"_Mitad-A___")
  seqkit seq -s $work_directory"/Inserciones_Internas/Registro_Exones/"$spe"/"$trans"_"$ace"_"$tipo"_A.fasta" | tr -d "\n" | sed "s/^/>$seq_name/" | sed "s/$/\n/" | sed "s/___/\n/" >> $work_directory"/Inserciones_Internas/"$spe"_inserciones_internas.fasta"

  seq_name=$(echo $trans"_"$ace"_"$local_ace_exon"_"$tipo"_Mitad-D___")
  seqkit seq -s $work_directory"/Inserciones_Internas/Registro_Exones/"$spe"/"$trans"_"$ace"_"$tipo"_D.fasta" | tr -d "\n" | sed "s/^/>$seq_name/" | sed "s/$/\n/" | sed "s/___/\n/" >> $work_directory"/Inserciones_Internas/"$spe"_inserciones_internas.fasta"
}

intron_secs () {
  recorrer_exones=$(grep -c . $file_exon_data)
  count=1

  while [ $count -lt $recorrer_exones ]
  do
    sense=$(sed -n $count"p" $file_exon_data | awk -F "\t" '{print $4}')
    if [ "$strand" == "+" ]
    then
      exon1=$(sed -n $count"p" $file_exon_data | awk -F "\t" '{print $3}')
      start=$(sed -n $count"p" $file_exon_data | awk -F "\t" '{print $2+1}')
      count=$(($count + 1))
      end=$(sed -n $count"p" $file_exon_data | awk -F "\t" '{print $1-1}')
      exon2=$(sed -n $count"p" $file_exon_data | awk -F "\t" '{print $3}')
    else
      exon1=$(sed -n $count"p" $file_exon_data | awk -F "\t" '{print $3}')
      end=$(sed -n $count"p" $file_exon_data | awk -F "\t" '{print $1-1}')
      count=$(($count + 1))
      start=$(sed -n $count"p" $file_exon_data | awk -F "\t" '{print $2+1}')
      exon2=$(sed -n $count"p" $file_exon_data | awk -F "\t" '{print $3}')
    fi

    ace_num=1
    while [ $ace_num -le $N_ace_exones ]
    do
      ace_exon=$(sed -n $ace_num"p" $work_directory"/Inserciones_Internas/Temporal/ACE_Exones")
      if [ "$exon2" == "$ace_exon" ]
      then
        exon2=$(echo $exon2"_ACE")
        ace_num=$(($N_ace_exones + $N_ace_exones))
      else
        ace_num=$(($ace_num + 1))
      fi
    done

    intron_name=$(echo $trans"_intron_"$exon1"-"$exon2)
    printf "%s %s %s %s\n%s %s %s\n" $chr" "$start"-"$end" "$sense | blastdbcmd -db $reference_blast"/"$spe"_ref" -entry_batch - | sed "s/>.*/>$intron_name/" | seqkit seq -w 0 >> $work_directory"/Inserciones_Internas/Intrones/"$spe"_intrones.fasta"
  done
}

sec_to_extract () {
  if [ $check_nonsense -gt 0 ]
  then
    exon_code=CDS
  else
    exon_code=exon
  fi
}

print_data_blast () {
  check_hit=$(grep -c $hit"_"$mitad $work_directory"/Anotacion/BLAST_Searches/"$spe"_"$spe2".blastx")

  if [ $check_hit -gt 0 ]
  then
    potential_hits=$(grep -c $hit"_"$mitad $work_directory"/Anotacion/BLAST_Searches/"$spe"_"$spe2".blastx")
    best_hit=$(grep $hit"_"$mitad $work_directory"/Anotacion/BLAST_Searches/"$spe"_"$spe2".blastx" | awk -F "\t" '{print $2"\t"$12*10}' | sort -n -k 2 | tail -n 1 | awk -F "\t" '{print $1}')
    best_iden=$(grep $hit"_"$mitad $work_directory"/Anotacion/BLAST_Searches/"$spe"_"$spe2".blastx" | awk -F "\t" '{print $3"\t"$12*10}' | sort -n -k 2 | tail -n 1 | awk -F "\t" '{print $1}')
    best_len=$(grep $hit"_"$mitad $work_directory"/Anotacion/BLAST_Searches/"$spe"_"$spe2".blastx" | awk -F "\t" '{print $4"\t"$12*10}' | sort -n -k 2 | tail -n 1 | awk -F "\t" '{print $1}')
    best_qcov=$(grep $hit"_"$mitad $work_directory"/Anotacion/BLAST_Searches/"$spe"_"$spe2".blastx" | awk -F "\t" '{print $13"\t"$12*10}' | sort -n -k 2 | tail -n 1 | awk -F "\t" '{print $1}')
    best_hit_bitscore=$(grep $hit"_"$mitad $work_directory"/Anotacion/BLAST_Searches/"$spe"_"$spe2".blastx" | awk -F "\t" '{print $2"\t"$12*10}' | sort -n -k 2 | tail -n 1 | awk -F "\t" '{print $2/10}')
  else
    potential_hits=0
    best_hit=$(echo "X")
    best_iden=$(echo "X")
    best_len=$(echo "X")
    best_qcov=$(echo "X")
    best_hit_bitscore=$(echo "X")
  fi
  length=$(seqkit grep -p $hit"_"$mitad $work_directory"/"$spe"_inserciones_internas.fasta" | seqkit fx2tab -n -l | awk -F "\t" '{print $2}')

  echo $hit" "$mitad" "$length" "$potential_hits" "$best_hit" "$best_iden" "$best_len" "$best_qcov" "$best_hit_bitscore" XSTATUSX" | tr " " "\t" >> $work_directory"/Anotacion/BLAST_Searches/"$spe"_"$spe2"_registro_hits.tab"
}

Sumary_Maker () {
  for spe2 in $references
  do
    if [ ! "$spe" == "$spe2" ]
    then
      count_All=$(awk -F "\t" '{print $1"\t"$10}' $work_directory"/Anotacion/BLAST_Searches/"$spe"_"$spe2"_"$data_file | sort -u | grep -c .)
      count_Internal=$(awk -F "\t" '{print $1"\t"$10}' $work_directory"/Anotacion/BLAST_Searches/"$spe"_"$spe2"_"$data_file | sort -u | grep -c Internal)
      count_Corto=$(awk -F "\t" '{print $1"\t"$10}' $work_directory"/Anotacion/BLAST_Searches/"$spe"_"$spe2"_"$data_file | sort -u | grep -c Corto)
      count_Missing=$(awk -F "\t" '{print $1"\t"$10}' $work_directory"/Anotacion/BLAST_Searches/"$spe"_"$spe2"_"$data_file | sort -u | grep -c Missing_Hit)
      count_Artifact=$(awk -F "\t" '{print $1"\t"$10}' $work_directory"/Anotacion/BLAST_Searches/"$spe"_"$spe2"_"$data_file | sort -u | grep -w -c Artifact)
      count_Unreliable=$(awk -F "\t" '{print $1"\t"$10}' $work_directory"/Anotacion/BLAST_Searches/"$spe"_"$spe2"_"$data_file | sort -u | grep -w -c U_Internal)
      count_Unreliable_Artifact=$(awk -F "\t" '{print $1"\t"$10}' $work_directory"/Anotacion/BLAST_Searches/"$spe"_"$spe2"_"$data_file | sort -u | grep -c U_Artifact)
    else
      count_All=o
      count_Internal=o
      count_Corto=o
      count_Missing=o
      count_Artifact=o
      count_Unreliable=o
      count_Unreliable_Artifact=o
    fi
    print_all=$(echo $print_all" "$count_All)
    print_Internal=$(echo $print_Internal" "$count_Internal)
    print_Corto=$(echo $print_Corto" "$count_Corto)
    print_Missing=$(echo $print_Missing" "$count_Missing)
    print_Artifact=$(echo $print_Artifact" "$count_Artifact)
    print_Unreliable=$(echo $print_Unreliable" "$count_Unreliable)
    print_Unreliable_Artifact=$(echo $print_Unreliable_Artifact" "$count_Unreliable_Artifact)
  done
  echo $print_all | tr " " "\t" >> $work_directory"/Anotacion/"$spe"_"$tier"_summary.tab"
  echo $print_Internal | tr " " "\t" >> $work_directory"/Anotacion/"$spe"_"$tier"_summary.tab"
  echo $print_Corto | tr " " "\t" >> $work_directory"/Anotacion/"$spe"_"$tier"_summary.tab"
  echo $print_Missing | tr " " "\t" >> $work_directory"/Anotacion/"$spe"_"$tier"_summary.tab"
  echo $print_Artifact | tr " " "\t" >> $work_directory"/Anotacion/"$spe"_"$tier"_summary.tab"
  echo $print_Unreliable | tr " " "\t" >> $work_directory"/Anotacion/"$spe"_"$tier"_summary.tab"
  echo $print_Unreliable_Artifact | tr " " "\t" >> $work_directory"/Anotacion/"$spe"_"$tier"_summary.tab"
}

get_tier_data2 () {
  count=2
  recorrer=$(grep -c . "Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts/"$spe"_"$tier".counts")

  while [ $count -le $recorrer ]
  do
    gen=$(sed -n $count"p" "Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts/"$spe"_"$tier".counts" | awk -F "\t" '{print $1}')
    ace=$(sed -n $count"p" "Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts/"$spe"_"$tier".counts" | awk -F "\t" '{print $2}')

    echo "Tier Nonsense!!!"
    echo $gen" "$ace

    transcripts=$(grep -w -F $gen "Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts/"$spe"_"*"_read.counts" | awk -F "\t" '{print $2}' | sort -u)
    for trans in $transcripts
    do
      for spe2 in $references
      do
        if [ ! "$spe" == "$spe2" ]
        then
          grep -F $trans"_"$ace"_" $work_directory"/Anotacion/BLAST_Searches/"$spe"_"$spe2"_registro_hits.tab" >> $work_directory"/Anotacion/BLAST_Searches/"$spe"_"$spe2"_registro_"$tier"_hits.tab"
        fi
      done
      grep -F $trans"_"$ace"_" $work_directory"/Anotacion/"$spe"_registro_global.tab" >> $work_directory"/Anotacion/"$spe"_registro_"$tier".tab"
      grep -F $trans"_"$ace"_" $work_directory"/"$spe"_inserciones_internas.tab" >> $work_directory"/"$spe"_inserciones_internas_"$tier".tab"
      grep -F $trans"_" $work_directory"/"$spe"_intrones.tab" >> $work_directory"/"$spe"_intrones_"$tier".tab"
    done
    count=$(($count + 1))
  done
}

################################################################################
################################################################################
################################################################################

if [ "$RUN_SLFINDER_GENES" == "TRUE" ]
then
  for gru in $grupos
  do
    especies=$(ls -d $gru"/"*"/" | sed "s/\/$//" | sed "s/.*\///")
    for spe in $especies
    do
      echo $spe
      cd $gru"/"$spe
      rm -r SL_Genes
      rm *seqkit.fai
      if [ $gru == "Cestodes" ]
      then
        cp $cestoda_sl_tags SL_tags.fa
      else
        cp $trematopda_sl_tags SL_tags.fa
      fi

      read_files=$(ls "trimmed_reads/"*"_trimmed_read_1P.gz" | sed "s/.*\///" | sed "s/_trimmed_read_1P.gz//")
      for read in $read_files
      do
        echo "Reformating: "$read
        seqkit seq "trimmed_reads/"$read"_trimmed_read_1P.gz" >> $read"_work_read_1P"
        gzip $read"_work_read_1P"

        seqkit seq "trimmed_reads/"$read"_trimmed_read_2P.gz" >> $read"_work_read_2P"
        gzip $read"_work_read_2P"
      done

      genome_file=$(ls *"genomic.fa")
      gff_file=$(ls *"annotations.gtf")

      echo "Debug SLFinder-Genes -t $threads -g $genome_file -gt $gff_file -r1 _work_read_1P.gz -r2 _work_read_2P.gz"
      ./SLFinder-Genes -t $threads -g $genome_file -gt $gff_file -r1 "_work_read_1P.gz" -r2 "_work_read_2P.gz"

      rm *"_work_read_1P.gz"
      rm *"_work_read_2P.gz"

      cd ../../
    done
  done
fi

################################################################################
################################################################################

if [ $RUN_ACE_PRE_CLASSIFICATION == "TRUE" ]
then
  rm -r "Analysis_Results/SL_Acceptor_Sites"

  mkdir "Analysis_Results/SL_Acceptor_Sites"
  mkdir "Analysis_Results/SL_Acceptor_Sites/SL_Read_Counts"
  mkdir "Analysis_Results/SL_Acceptor_Sites/Temporal"

  for filo in $grupos
  do
    echo $filo
    cd $filo
    species=$(ls -d *"/" | sed "s/\/$//" | sed "s/.*\///")

    if [ $gru == "Cestodes" ]
    then
      sl_tag_file=$cestoda_sl_tags
    else
      sl_tag_file=$trematopda_sl_tags
    fi

    sl_tags_filo=$(seqkit seq -i -n $sl_tag_file)

    echo "Especie N_Assemblies "$sl_tags_filo" Non_Canon Coting_limit" | tr " " "\t" >> ../Analysis_Results/SL_Acceptor_Sites/$filo"_Summary_hits.tab"
    echo "Especie N_Assemblies Total_reads Raw_SL_reads Mapped_SL_reads "$sl_tags_filo  | tr " " "\t" >> ../Analysis_Results/SL_Acceptor_Sites/$filo"_Summary_reads.tab"

    echo "Especie Total_sitios Iso_Complicada End_Of_Seq Discordantes Comienzo UTR Final Exon Intron_Can Intron_Alt" | tr " " "\t" >> ../Analysis_Results/SL_Acceptor_Sites/$filo"_Summary_All_tipos.tab"
    echo "Especie Total_sitios Iso_Complicada End_Of_Seq Discordantes Comienzo UTR Final Exon Intron_Can Intron_Alt" | tr " " "\t" >> ../Analysis_Results/SL_Acceptor_Sites/$filo"_Summary_tier_3_tipos.tab"
    echo "Especie Total_sitios Iso_Complicada End_Of_Seq Discordantes Comienzo UTR Final Exon Intron_Can Intron_Alt" | tr " " "\t" >> ../Analysis_Results/SL_Acceptor_Sites/$filo"_Summary_tier_10_tipos.tab"

    for spe in $species
    do
      if [ -d $spe"/SL_Genes" ]
      then
        gtf_file=$(ls $spe"/"*".annotations.gtf")
        genome_fasta=$(ls $spe"/"*"genomic.fa")

        check_wormbase=$(grep -c 'transcript_id "transcript:' $gtf_file)
        echo $spe" "$gtf_file
        echo "wormbase: "$check_wormbase

        seqkit stats -T $spe"/"*"1P.gz" | awk -F "\t" '{print $1"\t"$4}' >> ../Analysis_Results/SL_Acceptor_Sites/Temporal/$spe"_read_counts.temp"
        seqkit stats -T $spe"/SL_Genes/Results/read_SL/"*"_reads_filtered_1.fq" | awk -F "\t" '{print $1"\t"$4}' >> ../Analysis_Results/SL_Acceptor_Sites/Temporal/$spe"_sl_read_counts.temp"

        ensamblados_todos=$(ls $spe"/"*_Trinity.Trinity.fasta | grep -c .)

        SL_Tags_Hits=$(ls $spe"/SL_Genes/Results/Potential_aceptor_sites/"*"_aceptor_sites.tab" | sed "s/.*\///g" | sed "s/_/_x_/" | sed "s/.*_x_//" | sed "s/_aceptor_sites.tab//" | sort -u)
        for sl in $SL_Tags_Hits
        do
          no_known_genes=$(awk -F "\t" '{if ($2=="Y") print $1}' $spe"/SL_Genes/Results/Potential_aceptor_sites/"*"_"$sl"_aceptor_sites.tab" | grep -c .)

          if [ ! $no_known_genes -eq 0 ]
          then
            echo "Processing: "$sl
            archivos=$(ls $spe"/SL_Genes/Results/Potential_aceptor_sites/"*"_"$sl"_aceptor_sites.tab" | sed "s/.*\///g" | awk -F "_" '{print $1}' | sort -u)
            genes_aces=$(awk -F "\t" '{if ($2=="Y") print $1"|||"$3}' $spe"/SL_Genes/Results/Potential_aceptor_sites/"*"_"$sl"_aceptor_sites.tab" | sort -u)
            reads_counts
          else
            echo "No good alignments: "$sl
          fi
        done
        complete_sl_sum
        summary_sl_tags_info
      fi
    done
    echo "Completado "$filo
    cd ../
  done
fi

################################################################################
################################################################################

if [ $RUN_ACE_INTERNAL_SEC == "TRUE" ]
then
  work_directory="Analysis_Results/SL_Acceptor_Sites"
  reference_blast="Analysis_Results/Extra_SLs_Busqueda_BLAST/Referencias_BLAST"

  rm -r $work_directory"/Inserciones_Internas"
  rm $work_directory"/Inserciones_Internas/extract_secs.id"
  mkdir $work_directory"/Inserciones_Internas"
  mkdir $work_directory"/Inserciones_Internas/Registro_Exones"
  mkdir $work_directory"/Inserciones_Internas/Intrones"
  mkdir $work_directory"/Inserciones_Internas/Anotacion"
  mkdir $work_directory"/Inserciones_Internas/Temporal"

  for filo in $grupos
  do
    species=$(grep -F -v Iso_Complicada $work_directory"/"$filo"_Summary_All_tipos.tab" | awk -F "\t" '{print $1}')

    echo $species
    for spe in $species
    do
      echo $filo" "$spe
      mkdir $work_directory"/Inserciones_Internas/Registro_Exones/"$spe

      get_info=Intron_
      genes=$(grep -F -w -v Gen_ID $work_directory"/SL_Read_Counts/"$spe"_"*"_read.counts" | grep -F -w Concordante | grep -F $get_info | sed "s/.*://" | awk -F "\t" '{print $1}' | sort -u)

      gtf_file=$(ls $filo"/"$spe"/"*".gtf")
      check_wormbase=$(grep -c -F 'transcript_id "transcript:' $gtf_file)

      for gen in $genes
      do
        extract_info
        for trans in $transcripts
        do
          file_exon_data=$work_directory"/Inserciones_Internas/Temporal/Trans_Exon_temp_"$trans
          aceptor_sites=$(grep -F -w $trans $work_directory"/SL_Read_Counts/"$spe"_"*"_read.counts" | grep -F $trans | grep -F -w Concordante | grep Intron_ | awk -F "\t" '{print $5}' | sort -u)

          if [ $check_wormbase -gt 0 ]
          then
            if [ "$strand" == "+" ]
            then
              check_nonsense=$(grep -F 'transcript_id "transcript:'$trans'";' $gtf_file  | awk -F "\t" -v exon=CDS '{if ($3==exon) print}' | grep -c .)
              sec_to_extract
              grep -F 'transcript_id "transcript:'$trans'";' $gtf_file  | awk -F "\t" -v exon=$exon_code '{if ($3==exon) print $4"\t"$5}' | grep . | sort -n -k 1 | awk -F "\t" '{print $0"\tEx-"NR"\tplus"}' >> $file_exon_data
            else
              check_nonsense=$(grep -F 'transcript_id "transcript:'$trans'";' $gtf_file | awk -F "\t" -v exon=CDS '{if ($3==exon) print}' | grep -c .)
              sec_to_extract
              grep -F 'transcript_id "transcript:'$trans'";' $gtf_file  | awk -F "\t" -v exon=$exon_code '{if ($3==exon) print $5"\t"$4}' | grep . | sort -r -n -k 1 | awk -F "\t" '{print $2"\t"$1"\tEx-"NR"\tminus"}' >> $file_exon_data
            fi
          else
            if [ "$strand" == "+" ]
            then
              check_nonsense=$(grep -F $trans'";' $gtf_file | awk -F "\t" -v exon=CDS '{if ($3==exon) print}' | grep -c .)
              sec_to_extract
              grep -F $trans'";' $gtf_file  | awk -F "\t" -v exon=$exon_code '{if ($3==exon) print $4"\t"$5}' | grep . | sort -n -k 1 | awk -F "\t" '{print $0"\tEx-"NR"\tplus"}' >> $file_exon_data
            else
              check_nonsense=$(grep -F $trans'";' $gtf_file | awk -F "\t" -v exon=CDS '{if ($3==exon) print}' | grep -c .)
              sec_to_extract
              grep -F $trans'";' $gtf_file  | awk -F "\t" -v exon=$exon_code '{if ($3==exon) print $5"\t"$4}' | grep . | sort -r -n -k 1 | awk -F "\t" '{print $2"\t"$1"\tEx-"NR"\tminus"}' >> $file_exon_data
            fi
          fi

          for ace in $aceptor_sites
          do
            local_ace_exon=$(grep -F -w $trans $work_directory"/SL_Read_Counts/"$spe"_"*"_read.counts" | grep -F -w $ace | awk -F "\t" '{print $8}' | awk -F "_" '{print $2}' | sort -u)
            grep -F -w $trans $work_directory"/SL_Read_Counts/"$spe"_"*"_read.counts" | grep -F -w $ace | awk -F "\t" '{print $8}' | awk -F "_" '{print $2}' | sort -u >> $work_directory"/Inserciones_Internas/Temporal/ACE_Exones"
            tipo=Intron
            separator
          done
          N_ace_exones=$(grep -c . $work_directory"/Inserciones_Internas/Temporal/ACE_Exones")

          intron_secs
          rm $file_exon_data
          rm $work_directory"/Inserciones_Internas/Temporal/ACE_Exones"
        done
      done

      for gen in $genes
      do
        extract_info
        for trans in $transcripts
        do
          if [ $check_wormbase -gt 0 ]
          then
            if [ "$strand" == "+" ]
            then
              grep -F 'transcript_id "transcript:'$trans'";' $gtf_file  | awk -F "\t"  -v exon=$exon_code '{if ($3==exon) print $4"\t"$5}' | grep . | sort -n -k 1 | awk -F "\t" '{print $0"\tEx-"NR"\tplus"}' >> $work_directory"/Inserciones_Internas/Temporal/Trans_Exon_temp1_"$trans
            else
              grep -F 'transcript_id "transcript:'$trans'";' $gtf_file  | awk -F "\t" -v exon=$exon_code '{if ($3==exon) print $5"\t"$4}' | grep . | sort -r -n -k 1 | awk -F "\t" '{print $2"\t"$1"\tEx-"NR"\tminus"}' >> $work_directory"/Inserciones_Internas/Temporal/Trans_Exon_temp1_"$trans
            fi
          else
            if [ "$strand" == "+" ]
            then
              grep -F $trans'";' $gtf_file | awk -F "\t"  -v exon=$exon_code '{if ($3==exon) print $4"\t"$5}' | grep . | sort -n -k 1 | awk -F "\t" '{print $0"\tEx-"NR"\tplus"}' >> $work_directory"/Inserciones_Internas/Temporal/Trans_Exon_temp1_"$trans
            else
              grep -F $trans'";' $gtf_file | awk -F "\t" -v exon=$exon_code '{if ($3==exon) print $5"\t"$4}' | grep . | sort -r -n -k 1 | awk -F "\t" '{print $2"\t"$1"\tEx-"NR"\tminus"}' >> $work_directory"/Inserciones_Internas/Temporal/Trans_Exon_temp1_"$trans
            fi
          fi

          # Comprobar inserciones Anomaly

          aceptor_sites=$(grep -F -w $trans $work_directory"/SL_Read_Counts/"$spe"_"*"_read.counts" | grep -F $trans | grep -w Concordante | grep Exon | awk -F "\t" '{print $5}' | sort -u)
          for ace in $aceptor_sites
          do
            check_anomaly=$(grep -w -c $ace $work_directory"/Inserciones_Internas/Temporal/Trans_Exon_temp1_"$trans)
            if [ $check_anomaly -gt 0 ]
            then
              info_anomaly=$(grep -w $ace $work_directory"/Inserciones_Internas/Temporal/Trans_Exon_temp1_"$trans)
              echo $filo" "$spe" "$trans" "$info_anomaly | tr " " "\t" >> $work_directory"/Inserciones_Internas/Strange_exons.tab"
            else
              file_exon_data=$work_directory"/Inserciones_Internas/Temporal/Trans_Exon_temp_"$trans"_"$ace
              local_ace_exon=$(grep -F -w $trans $work_directory"/SL_Read_Counts/"$spe"_"*"_read.counts" | grep -F -w $ace | awk -F "\t" '{print $8}' | sort -u)
              echo "Exon_info "$gen" "$strand": "$trans" --- "$ace" || "$local_ace_exon

              recorrer_exones=$(grep -c . $work_directory"/Inserciones_Internas/Temporal/Trans_Exon_temp1_"$trans)
              count=1

              while [ $count -le $recorrer_exones ]
              do
                exon_turno=$(sed -n $count"p" $work_directory"/Inserciones_Internas/Temporal/Trans_Exon_temp1_"$trans | awk -F "\t" '{print $3}')

                if [ "$exon_turno" == "$local_ace_exon" ]
                then
                  start_ex=$(sed -n $count"p" $work_directory"/Inserciones_Internas/Temporal/Trans_Exon_temp1_"$trans | awk -F "\t" '{print $1}')
                  end_ex=$(sed -n $count"p" $work_directory"/Inserciones_Internas/Temporal/Trans_Exon_temp1_"$trans | awk -F "\t" '{print $2}')
                  sentido=$(sed -n $count"p" $work_directory"/Inserciones_Internas/Temporal/Trans_Exon_temp1_"$trans | awk -F "\t" '{print $4}')

                  if [ "$strand" == "+" ]
                  then
                    echo $start_ex" "$ace" "$exon_turno" "$sentido | tr " " "\t" >> $file_exon_data
                    cut_off=$(($ace + 1))
                    echo $cut_off" "$end_ex" Cut_"$exon_turno" "$sentido | tr " " "\t" >> $file_exon_data
                  else
                    echo $ace" "$end_ex" "$exon_turno" "$sentido | tr " " "\t" >> $file_exon_data
                    cut_off=$(($ace - 1))
                    echo $start_ex" "$cut_off" Cut_"$exon_turno" "$sentido | tr " " "\t" >> $file_exon_data
                  fi
                else
                  sed -n $count"p" $work_directory"/Inserciones_Internas/Temporal/Trans_Exon_temp1_"$trans >> $file_exon_data
                fi
                count=$(($count + 1))
              done
              local_ace_exon=$(echo "Cut_"$local_ace_exon)
              tipo=Exon
              separator
              rm $file_exon_data
            fi
          done
          rm $work_directory"/Inserciones_Internas/Temporal/Trans_Exon_temp1_"$trans
        done
      done

      seqkit fx2tab -n -i -l -g $work_directory"/Inserciones_Internas/"$spe"_inserciones_internas.fasta" >> $work_directory"/Inserciones_Internas/"$spe"_inserciones_internas.tab"
      seqkit fx2tab -n -i -l -g $work_directory"/Inserciones_Internas/Intrones/"$spe"_intrones.fasta" >> $work_directory"/Inserciones_Internas/"$spe"_intrones.tab"

      tier=tier3
      get_tier_data

      tier=tier10
      get_tier_data
    done
  done
fi

################################################################################
################################################################################


if [ "$RUN_BLAST_INTERNAL_STIES" == "TRUE" ]
then
  work_directory=Analysis_Results/SL_Acceptor_Sites/Inserciones_Internas

  rm -r $work_directory"/Anotacion"
  rm $work_directory"/"*".tab"
  mkdir $work_directory"/Anotacion"
  mkdir $work_directory"/Anotacion/References_Prot"
  mkdir $work_directory"/Anotacion/BLAST_Searches"

  for filo in $grupos
  do
    especies=$(awk -F "\t" '{print $1}' "Analysis_Results/SL_Acceptor_Sites/"$filo"_Summary_hits.tab" | grep -v -w "Especie")
    set_references

    for spe in $references
    do
      echo "Referencia Blast: "$spe
      makeblastdb -dbtype prot -in "Ref_Proteins/"$spe".protein.fa" -out $work_directory"/Anotacion/References_Prot/"$spe"_ref"
    done

    for spe in $especies
    do
      seqkit fx2tab -n -i -l -g $work_directory"/"$spe"_inserciones_internas.fasta" >> $work_directory"/"$spe"_inserciones_internas.tab"
      seqkit fx2tab -n -i -l -g $work_directory"/Intrones/"$spe"_intrones.fasta" >> $work_directory"/"$spe"_intrones.tab"

      awk -F "\t" '{if ($2>=60) print $1}' $work_directory"/"$spe"_inserciones_internas.tab" | sed "s/_Ex-.*//" | sort | uniq -d >> $work_directory"/extract_secs.id"
      seqkit grep -r -f $work_directory"/extract_secs.id" $work_directory"/"$spe"_inserciones_internas.fasta" >> $work_directory"/Anotacion/"$spe"_inserciones_internas_anot.fasta"
      rm $work_directory"/extract_secs.id"
    done

    for spe in $especies
    do
      for spe2 in $references
      do
        if [ ! "$spe" == "$spe2" ]
        then
          echo "Realizando BLAST"
          echo $spe" vs "$spe2

          blastx -query $work_directory"/Anotacion/"$spe"_inserciones_internas_anot.fasta" -db $work_directory"/Anotacion/References_Prot/"$spe2"_ref" -outfmt "6 std qcovs" -out $work_directory"/Anotacion/BLAST_Searches/"$spe"_"$spe2".blastx" -num_threads 24 -qcov_hsp_perc $qcov
        fi
      done
    done
  done
fi

if [ "$RUN_CLASSIFY_INTERNAL_STIES" == "TRUE" ]
then
  work_directory=Analysis_Results/SL_Acceptor_Sites/Inserciones_Internas

  for filo in $grupos
  do
    especies=$(awk -F "\t" '{print $1}' "Analysis_Results/SL_Acceptor_Sites/"$filo"_Summary_hits.tab" | grep -v -w "Especie")
    echo $especies

    set_references

    for spe in $especies
    do
      print_all=All
      print_Internal=Internal
      print_Corto=Corto
      print_Missing=Missing
      print_Artifact=Artifact
      print_Unreliable=U_Internal
      print_Unreliable_Artifact=U_Artifact

      hits_analizar=$(seqkit seq -n $work_directory"/"$spe"_inserciones_internas.fasta" | grep "_Mitad-A" | sed "s/_Mitad-A//" | sort -u)

      for spe2 in $references
      do
        if [ ! "$spe" == "$spe2" ]
        then
          echo $spe" vs "$spe2
          echo "Gen Half Len N_Hits Best_Hit Ident Matches Qcov Bitscore Status" | tr " " "\t" >> $work_directory"/Anotacion/BLAST_Searches/"$spe"_"$spe2"_registro_hits.tab"

          for hit in $hits_analizar
          do
            check_hit_half=$(grep -c $hit $work_directory"/Anotacion/"$spe"_inserciones_internas_anot.fasta")
            if [ $check_hit_half -eq 0 ]
            then
              mitad=Mitad-A
              length=$(seqkit grep -p $hit"_"$mitad $work_directory"/"$spe"_inserciones_internas.fasta" | seqkit fx2tab -n -l | awk -F "\t" '{print $2}')
              echo $hit" "$mitad" "$length" NA NA NA NA NA NA Corto" | tr " " "\t" >> $work_directory"/Anotacion/BLAST_Searches/"$spe"_"$spe2"_registro_hits.tab"

              mitad=Mitad-D
              length=$(seqkit grep -p $hit"_"$mitad $work_directory"/"$spe"_inserciones_internas.fasta" | seqkit fx2tab -n -l | awk -F "\t" '{print $2}')
              echo $hit" "$mitad" "$length" NA NA NA NA NA NA Corto" | tr " " "\t" >> $work_directory"/Anotacion/BLAST_Searches/"$spe"_"$spe2"_registro_hits.tab"
            else
              mitad=Mitad-A
              print_data_blast
              best_hit_1=$best_hit

              mitad=Mitad-D
              print_data_blast
              best_hit_2=$best_hit

              if [ "$best_hit_1" == "X" ] || [ "$best_hit_2" == "X" ]
              then
                sed -i "/$hit/s/XSTATUSX$/Missing_Hit/" $work_directory"/Anotacion/BLAST_Searches/"$spe"_"$spe2"_registro_hits.tab"
              elif [ "$best_hit_1" == "$best_hit_2" ]
              then
                sed -i "/$hit/s/XSTATUSX$/Internal/" $work_directory"/Anotacion/BLAST_Searches/"$spe"_"$spe2"_registro_hits.tab"
              else
                if [ $potential_hits -gt 1 ]
                then
                  hits_A=$(grep $hit"_Mitad-A" $work_directory"/Anotacion/BLAST_Searches/"$spe"_"$spe2".blastx" | awk -F "\t" '{print $2}' | sort -u)
                  hits_D=$(grep $hit"_Mitad-D" $work_directory"/Anotacion/BLAST_Searches/"$spe"_"$spe2".blastx" | awk -F "\t" '{print $2}' | sort -u)

                  check_unreliable=$(echo $hits_A" "$hits_D | tr " " "\n" | sort | uniq -d | grep -c .)
                  if [ $check_unreliable -gt 0 ]
                  then
                    sed -i "/$hit/s/XSTATUSX$/U_Internal/" $work_directory"/Anotacion/BLAST_Searches/"$spe"_"$spe2"_registro_hits.tab"
                  else
                    sed -i "/$hit/s/XSTATUSX$/U_Artifact/" $work_directory"/Anotacion/BLAST_Searches/"$spe"_"$spe2"_registro_hits.tab"
                  fi
                else
                  sed -i "/$hit/s/XSTATUSX$/Artifact/" $work_directory"/Anotacion/BLAST_Searches/"$spe"_"$spe2"_registro_hits.tab"
                fi
              fi
            fi
          done
        fi
      done

      data_file=registro_hits.tab
      tier=All
      echo "Tipo "$references | tr " " "\t" >> $work_directory"/Anotacion/"$spe"_"$tier"_summary.tab"
      echo "Hit "$references | tr " " "\t" >> $work_directory"/Anotacion/"$spe"_registro_global.tab"
      Sumary_Maker

      for hit in $hits_analizar
      do
        print_result=$hit
        for spe2 in $references
        do
          if [ "$spe" == "$spe2" ]
          then
            resultado=$(echo "o")
          else
            check_status=$(grep -w $hit $work_directory"/Anotacion/BLAST_Searches/"$spe"_"$spe2"_registro_hits.tab" | awk -F "\t" '{print $10}' | sort -u)

            if [ "$check_status" == "Internal" ]
            then
              resultado=$(echo ".")
            elif [ "$check_status" == "Corto" ]
            then
              resultado=$(echo "na")
            elif [ "$check_status" == "Missing_Hit" ]
            then
              resultado=$(echo "m")
            elif [ "$check_status" == "U_Internal" ]
            then
              resultado=$(echo "u")
            elif [ "$check_status" == "U_Artifact" ]
            then
              resultado=$(echo "U_Artifact")
            elif [ "$check_status" == "Artifact" ]
            then
              resultado=$(echo "Artifact")
            fi
          fi
          print_result=$(echo $print_result" "$resultado)
        done
        echo $print_result | tr " " "\t" >> $work_directory"/Anotacion/"$spe"_registro_global.tab"
      done

      tier=tier3
      data_file="registro_"$tier"_hits.tab"
      echo "Tipo "$references | tr " " "\t" > $work_directory"/Anotacion/"$spe"_"$tier"_summary.tab"
      echo "Hit "$references | tr " " "\t" > $work_directory"/Anotacion/"$spe"_registro_"$tier".tab"

      print_all=All
      print_Internal=Internal
      print_Corto=Corto
      print_Missing=Missing
      print_Artifact=Artifact
      print_Unreliable=U_Internal
      print_Unreliable_Artifact=U_Artifact

      get_tier_data2
      Sumary_Maker

      tier=tier10
      data_file="registro_"$tier"_hits.tab"


      echo "Tipo "$references | tr " " "\t" > $work_directory"/Anotacion/"$spe"_"$tier"_summary.tab"
      echo "Hit "$references | tr " " "\t" > $work_directory"/Anotacion/"$spe"_registro_"$tier".tab"

      print_all=All
      print_Internal=Internal
      print_Corto=Corto
      print_Missing=Missing
      print_Artifact=Artifact
      print_Unreliable=U_Internal
      print_Unreliable_Artifact=U_Artifact

      get_tier_data2
      Sumary_Maker
    done
  done
fi
