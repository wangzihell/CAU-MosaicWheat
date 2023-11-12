#!/bin/bash
set -euxo pipefail

for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B;do

  # number of type switchs between adjacent bin
  gawk -F"\t" -vOFS="\t" 'function abs(v) {return v<0?-v:v} {if(NR==2){for(i=4;i<=NF;i++){a[i]=abs($i)}};if(NR>2){for(i=4;i<=NF;i++){if(abs($i)!=a[i]){sum+=1};a[i]=abs($i)};print $1,$2,$3,sum;sum=0}}' ${CHR}.grp > ${CHR}.simibins

  # smoothing using median of window size of 8
  Rscript color_switch_smoothing.R -i ${CHR}.simibins

  # detect START and END postion, by continuous beyond thres
  gawk -F"\t" -vOFS="\t" -vthr=30 'BEGIN{sig=0} ARGIND==1{if($4<thr){posi+=1;neg=0}else{posi=0;neg+=1};if(posi>=5&&sig==0){s=NR-4;sig=1};if(neg>=5&&sig==1){e=NR-5;nextfile}} ARGIND==2{if((FNR>=s+2&&FNR<=e+2)||FNR==1){print}}' ${CHR}.simibins ${CHR}.grp > ${CHR}.centromo
    
  # clustering based on centromo
  Rscript ${WD}/bin/centromere_clustering.R -i ${CHR}.centromo -o ${CHR}.centromo.cluster

done
