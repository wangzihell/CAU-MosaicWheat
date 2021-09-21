#!/bin/bash

if false;then
echo > dist.grp
for i in ../../data/chr??.grp;do
  tail -n +2 $i|cut -f 1-3 --complement >> dist.grp
done

Rscript AHG_based_dist.R -i dist.grp -o dist_mat_normedCNV.txt

Rscript /data2/rawdata2/tree/script/NJtree.R -d dist_mat_normedCNV.txt -s samplelist.txt -a /data2/rawdata2/sample_metadata/tetra_intro/metadata_210112.txt -p samplelist.txt -c colorF_9cols.txt -P T -T tiplab.txt -o dist_mat_normedCNV.pdf

Rscript /data2/rawdata2/tree/script/NJtree_iTOL.R -d dist_mat_normedCNV.txt -s samplelist.txt -a /data2/rawdata2/sample_metadata/tetra_intro/metadata_210112_itol.txt -p samplelist.txt -c ../colorF.txt -o dist_mat_normedCNV

cp /data2/rawdata2/tree/script/dataset_symbols_template.txt dataset_symbols.txt
# gawk -vOFS="," 'ARGIND==1{a[$1]=$2} ARGIND==2{b[$1]=$3;c[$1]=$6} ARGIND==3{print b[$1],2,1,a[c[$1]],1,1}' geog_color.txt /data2/rawdata2/sample_metadata/tetra_intro/metadata_210112_itol.txt samplelist.txt >> dataset_symbols.txt
gawk -vOFS="," 'ARGIND==1{a[$1]=$2} ARGIND==2{b[$1]=$3;c[$1]=$6} ARGIND==3{print b[$1],2,1,a[$2],1,1}' clade_color.txt /data2/rawdata2/sample_metadata/tetra_intro/metadata_210112_itol.txt WE_founder_clade.txt >> dataset_symbols.txt

# correspondence between pericentromeric and whole genome tree.
Rscript /data2/rawdata2/tree/script/NJtree.R -d dist_mat_normedCNV.txt -s samplelist.txt -a /data2/rawdata2/sample_metadata/tetra_intro/metadata_210112.txt -p samplelist.txt -c ../colorF.txt -T tiplab.txt -o WGS_wholegenome.pdf -l rectangular -W 7 -H 60

gawk 'ARGIND==1{a[$1]=$5} ARGIND==2{if((FNR==1)||(a[$1]!="CV")){print}}' /data2/rawdata2/sample_metadata/tetra_intro/metadata_210112.txt samplelist.txt > samplelist_WEDOLR.txt
Rscript /data2/rawdata2/tree/script/NJtree.R -d dist_mat_normedCNV.txt -s samplelist.txt -a /data2/rawdata2/sample_metadata/tetra_intro/metadata_210112.txt -p samplelist_WEDOLR.txt -c ../colorF.txt -o WGS_WEDOLR_whole.pdf -l rectangular -W 5 -H 30
fi

# large block detection between Western and Eastern LR using heatmap sorting.
# gawk -F"\t" -vOFS="\t" 'ARGIND==1{a[$1]=$5;b[$1]=$6;c[$1]=$4} ARGIND==2{if(a[$1]=="LR"&&c[$1]!="SP"){if(b[$1]=="Europe"||b[$1]=="WestAsia"){print FNR,$1,"Western"}else if(b[$1]=="EastAsia"||b[$1]=="CentralAsia"){print FNR,$1,"Eastern"}}}' /data2/rawdata2/sample_metadata/tetra_intro/metadata_210112.txt samplelist.txt > LR_order_WestEast.txt
gawk -F"\t" -vOFS="\t" 'ARGIND==1{a[$1]=$5;b[$1]=$6;c[$1]=$4} ARGIND==2{if(a[$1]=="LR"){if(b[$1]=="Europe"||b[$1]=="WestAsia"){print FNR,$1,"Western"}else if(b[$1]=="EastAsia"||b[$1]=="CentralAsia"){print FNR,$1,"Eastern"}}}' /data2/rawdata2/sample_metadata/tetra_intro/metadata_210112.txt samplelist.txt > LR_order_WestEast.txt
gawk 'ARGIND==1{a[$1]} ARGIND==2{if(FNR in a){for(i=1;i<=NF;i++){if(i in a){printf $i"\t"}};printf "\n"}}' LR_order_WestEast.txt dist_mat_normedCNV.txt|sed 's/\t$//'  > LR_WestEast.dist.txt

# add CV accessions
gawk -F"\t" -vOFS="\t" 'ARGIND==1{a[$1]=$5;b[$1]=$6;c[$1]=$4} ARGIND==2{if(a[$1]=="LR"&&c[$1]!="SP"){if(b[$1]=="Europe"||b[$1]=="WestAsia"){print FNR,$1,"Western","LR"}else if(b[$1]=="EastAsia"||b[$1]=="CentralAsia"){print FNR,$1,"Eastern","LR"}};if(a[$1]=="CV"){if(b[$1]=="Europe"||b[$1]=="WestAsia"){print FNR,$1,"Western","CV"}else if(b[$1]=="EastAsia"||b[$1]=="CentralAsia"){print FNR,$1,"Eastern","CV"}}}' /data2/rawdata2/sample_metadata/tetra_intro/metadata_210112.txt samplelist.txt > LRCV_order_WestEast.txt
gawk 'ARGIND==1{a[$1]} ARGIND==2{if(FNR in a){for(i=1;i<=NF;i++){if(i in a){printf $i"\t"}};printf "\n"}}' LR_order_WestEast.txt dist_mat_normedCNV.txt|sed 's/\t$//'  > LRCV_WestEast.dist.txt

# all accessions
Rscript heatmap.R -d dist_mat_normedCNV.txt -s samplelist.txt -c colorfile.txt -W 7 -H 4 -o whole_genome_heatmap.pdf

# all but CV accessions
gawk -F"\t" -vOFS="\t" 'ARGIND==1{a[$1]=$5} ARGIND==2{if(a[$1]!="CV"){print FNR,$1}}' /data2/rawdata2/sample_metadata/tetra_intro/metadata_210112.txt samplelist.txt > WEDOLR_ordername.txt
gawk 'ARGIND==1{a[$1]} ARGIND==2{if(FNR in a){for(i=1;i<=NF;i++){if(i in a){printf $i"\t"}};printf "\n"}}' WEDOLR_ordername.txt dist_mat_normedCNV.txt|sed 's/\t$//'  > WEDOLR.dist.txt
cut -f 2 WEDOLR_ordername.txt > WEDOLR_samplelist.txt
