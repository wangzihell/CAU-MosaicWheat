#!/bin/bash
set -euxo pipefail

# A&B subgenomes
for CHR in chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B;do
  n=$(wc -l < /data2/rawdata2/variant_density/CNV_masker/201224/${CHR}.CNVfilter.5M.txt)
  > ${CHR}_AABB_HH.dist; > ${CHR}_AABB_DT.dist; > ${CHR}_AABB_all.dist
  (for i in $(seq 1 $n);do
     # between hexaploid samples
     gawk -f dist_pair_sample.awk -vBin=$i /data2/rawdata2/variant_density/CNV_masker/201224/${CHR}.CNVfilter.5M.txt AABB_HH_linen_rand1000.txt ../raw/${CHR}.${i}.5M.dist.txt >> ${CHR}_AABB_HH.dist

     # between domesticated tetraploid samples
     gawk -f dist_pair_sample.awk -vBin=$i /data2/rawdata2/variant_density/CNV_masker/201224/${CHR}.CNVfilter.5M.txt AABB_DT_linen_rand1000.txt ../raw/${CHR}.${i}.5M.dist.txt >> ${CHR}_AABB_DT.dist

     # between all samples
     gawk -f dist_pair_sample.awk -vBin=$i /data2/rawdata2/variant_density/CNV_masker/201224/${CHR}.CNVfilter.5M.txt AABB_all_linen_rand1000.txt ../raw/${CHR}.${i}.5M.dist.txt >> ${CHR}_AABB_all.dist
done ) &
done

cat chr??_AABB_HH.dist |shuf -n 100000 | sort -n > AABB_HH_sample_by_count_flat.txt
cat chr??_AABB_DT.dist |shuf -n 100000 | sort -n > AABB_DT_sample_by_count_flat.txt

# D subgenome
for CHR in chr1D chr2D chr3D chr4D chr5D chr6D chr7D;do
  (n=$(wc -l < /data2/rawdata2/variant_density/CNV_masker/201224/${CHR}.CNVfilter.5M.txt)
  > ${CHR}_DD_HH.dist;
  for i in $(seq 1 $n);do
     # between hexaploid samples
     gawk -f dist_pair_sample_DD.awk -vBin=$i /data2/rawdata2/variant_density/CNV_masker/201224/${CHR}.CNVfilter.5M.txt DD_HH_linen_rand1000.txt ../raw/${CHR}.${i}.5M.dist.txt >> ${CHR}_DD_HH.dist
done ) &
done

cat chr?D_DD_HH.dist |shuf -n 100000 | sort -n > DD_HH_sample_by_count_flat.txt

# modeling
Rscript distribution_ABD.R
