#!/bin/bash
set -euxo pipefail

WD=/data2/rawdata2/tetraintro/
for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B;do
  (if true;then
  # number of type switchs between adjacent bin
  gawk -F"\t" -vOFS="\t" 'function abs(v) {return v<0?-v:v} {if(NR==2){for(i=4;i<=NF;i++){a[i]=abs($i)}};if(NR>2){for(i=4;i<=NF;i++){if(abs($i)!=a[i]){sum+=1};a[i]=abs($i)};print $1,$2,$3,sum;sum=0}}' ../../data/${CHR}.grp > ${CHR}.simibins
  # smoothing using median of window size of 8
  Rscript ${WD}/bin/color_switch_smoothing.R -i ${CHR}.simibins

  if [[ $CHR =~ "A" ]];then thr=30; else thr=50;fi
  detect START and END postion, by continuous beyond thres
  gawk -F"\t" -vOFS="\t" -vthr=30 'BEGIN{sig=0} ARGIND==1{if($4<thr){posi+=1;neg=0}else{posi=0;neg+=1};if(posi>=5&&sig==0){s=NR-4;sig=1};if(neg>=5&&sig==1){e=NR-5;nextfile}} ARGIND==2{if((FNR>=s+2&&FNR<=e+2)||FNR==1){print}}' ${CHR}.simibins ../../data/${CHR}.grp > ${CHR}.centromo
  # mannual
  gawk -F"\t" -vOFS="\t" -vs=$s -ve=$e '{if(($2>=s&&$3<=e)||NR==1){print}}' ../../data/${CHR}.grp > ${CHR}.centromo
    
  # clustering based on centromo
  Rscript ${WD}/bin/centromere_clustering.R -i ${CHR}.centromo -o ${CHR}.centromo.cluster.bk
  # centromo contribution according to their ordering in .fi file
  Rscript ${WD}/bin/centromere_contribution.R -i ${CHR}.centromo.cluster.bk -s ../../data/${CHR}.fi -o ${CHR}.contri
  # change group identifier number to its order 
  gawk 'ARGIND==1{if($1>1){a[$2]=NR}else{a[$2]="unique"}} ARGIND==2{print $1"\t"a[$2]}' <(cut -f 2 ${CHR}.centromo.cluster.bk|sort |uniq -c |sort -k1,1nr) ${CHR}.centromo.cluster.bk | sponge ${CHR}.centromo.cluster
  fi
  plotH=$(gawk 'NR==1{print (NF-3)/5}' ../../data/${CHR}.grp)
  Rscript ${WD}/bin/draw_haplo_block_v2.R -H ${plotH} -m /data2/rawdata2/sample_metadata/tetra_intro/metadata_201221.txt -i ../../data/${CHR}.contri.grp -o ${CHR}_centromere.pdf -s ../../data/${CHR}.contri.fi -n 20 -c ${WD}/bin/20_distinct_colors2.txt -g ${CHR}.centromo.cluster -C "$(gawk 'NR==2{print $2}' ${CHR}.centromo),$(gawk '{a=$3} END{print a}' ${CHR}.centromo),$(gawk -vCHR=$CHR '$1==CHR{print $3}' /data2/rawdata2/readDepth//DP_bySampel/CS_karyotype.txt)") &
done
