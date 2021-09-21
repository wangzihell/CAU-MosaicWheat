#!/bin/bash
set -euxo pipefail

if true;then

# for each of four taxonomic groups
for GROUP in CV LR DO WE;do
n=$(wc -l < ${GROUP}.order)
t=100  # rounds

# shuf ${GROUP} sample order
for j in $(seq $t);do
(shuf ${GROUP}.order > ${GROUP}.order.shuf${j}
for i in $(seq $n);do
for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B;do
  gawk -F"\t" -vOFS="\t" 'ARGIND==1{a[$i+3]} ARGIND==2{if(FNR>1){for(i=4;i<=NF;i++){if(i in a){b[FNR][$i]}}}} END{for(i in b){sum+=length(b[i])};print sum}' <(head -n $i ${GROUP}.order.shuf${j} ) ../../../data/${CHR}.grp
done > shuf_${j}_${i}.txt
done) &
sleep 5
done

wait

> saturation_100_${GROUP}.txt
for j in $(seq $t);do
for i in $(seq $n);do
gawk -vi=$i '{sum+=$1} END{print i,sum}' shuf_${j}_${i}.txt
done
done >> saturation_100_${GROUP}.txt

gawk -vOFS="\t" '{$2=$2/2030;print}' saturation_100_${GROUP}.txt| sort -k1,1n | datamash -g 1 perc:5 2 mean 2 perc:95 2 > saturation_100_${GROUP}.percentile.txt
done
fi
