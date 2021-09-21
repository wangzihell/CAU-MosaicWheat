#!/bin/bash
set -euxo pipefail

WD=/data2/rawdata2/tetraintro/

# aggregate 3B PAH from various dataset
gawk -F"\t" -vOFS="\t" '$5=="CL"{print $1 > "CL.txt"} $5=="NCL"{print $1 > "NCL.txt"}' /data2/rawdata2/sample_metadata/tetra_intro/metadata_210112_spCVLR.txt
tail -n +297 /data2/rawdata2/sample_metadata/tetra_intro/metadata_NG20.txt|gawk -F"\t" -vOFS="\t" '{if($3=="Landrace"){if($2=="China"){print $1>>"CL.txt"}else{print $1>>"NCL.txt"}}}'
gawk -F"\t" -vOFS="\t" '$4=="CL"{print $1>>"CL.txt"} $5=="NCL"{print $1>>"NCL.txt"}' /data2/rawdata2/sample_metadata/tetra_intro/metadata_MP20.txt
gawk -F"\t" -vOFS="\t" '{if($5=="LANDRACE"){if($6=="CHINA"){print $1>>"CL.txt"}else{print $1>>"NCL.txt"}}}' /data2/rawdata2/sample_metadata/tetra_intro/PassportData_160809.txt
gawk -F"\t" -vOFS="\t" 'BEGIN{print "Accession\thap1\thap2\thap3"} ARGIND==1{a[$1]} ARGIND==2{if($1 in a){printf $1;for(i=1;i<=3;i++){if($2==i){printf "\t"1}else{printf "\t"0}};printf "\n"}}' <(cat CL.txt NCL.txt) <(cut -f 1,11 whole.cluster) > clade_map_input.txt

# infer the spreading route of 3 main 3B PAH
for H in 1 2 3;do
  # select samples for each haplotype
  gawk -vH=$H 'BEGIN{print "Accession"} ARGIND==1{a[$1]} ARGIND==2{b[$1]} ARGIND==3{if(($1 in a)&&($1 in b)&&($11==H)){print $1}}' sample_location_Eurasia.txt LR_samplelist.txt allsample.centromo.cluster > hap${H}_samplist.txt
  
  # merge samples according to their locations, and filter sites with too less samples (<3?)
  Rscript merge_adjacent_location.R -g sample_location_tetraintro_WEC_201221.txt -i hap${H}_samplist.txt -r 1 -m 3 -o geog_hap${H}.txt
  gawk -F"\t" -vOFS="\t" '{print "chr3B",280000000,415000000,$1 > $2".bed"}' geog_hap${H}.txt
  
  ## 
  tail -n +2 geog_hap${H}.txt |cut -f 2|sort|uniq > geog_list_hap${H}.txt

  parallel -j 10 --joblog parallel.log python seg_diversity_across.py -d /data2/rawdata2/variant_density/6.cross_sample_5M_200923/raw/ -a {1}.bed -b {2}.bed -s name.txt -o {1}_{2}.dist :::: geog_list_hap${H}.txt :::: geog_list_hap${H}.txt
  
  while read i;do
  while read j;do
    gawk '{sum+=$1} END{printf sum/NR"\t"}' ${i}_${j}.dist
  done < geog_list_hap${H}.txt
  printf "\n"
  done < geog_list_hap${H}.txt > hap${H}.dist
  
  Rscript minimum_spanning_tree.R -i hap${H}.dist -o hap${H}.mst -g geog_list_hap${H}.txt
  gawk -F"\t" -vOFS="\t" 'BEGIN{print "x1\ty1\tx3\ty3"}ARGIND==1{lon[$2]=$3;lat[$2]=$4} ARGIND==2{if(FNR==1){for(i=1;i<=NF;i++){a[i]=$i}}else{for(i=2;i<=NF;i++){if($i==1){print lon[a[FNR-1]],lat[a[FNR-1]],lon[a[i-1]],lat[a[i-1]]}}}}' geog_hap${H}.txt hap${H}.mst > hap${H}.linedata.txt
done

# convert to plot format
gawk -F"\t" -vOFS="\t" 'NR==1{$5="color";$6="width";$7="direct"} NR>1{$5="darkblue";$6=2;$7="down"} {print}' hap1.linedata.txt > hap123.linedata.txt
gawk -F"\t" -vOFS="\t" 'NR>1{$5="orange";$6=2;$7="up";print}' hap2.linedata.txt >> hap123.linedata.txt
gawk -F"\t" -vOFS="\t" 'NR>1{$5="darkred";$6=2;$7="up";print}' hap3.linedata.txt >> hap123.linedata.txt

# to plot
Rscript hapmap_mst.R -g sample_location_tetraintro_WEC_201221.txt --GeoRange="30,120,20,50" -i clade_map_input.txt -c /data2/rawdata2/hapmap/tetraintro/201111_china_local/colorF.txt -o chr3B_china_mst_all.pdf -m 3 -s 0.25 -l hap123.linedata.txt
