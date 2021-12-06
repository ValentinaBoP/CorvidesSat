/proj/sllstore2017073/private/scripts/RM2BED.sh -i lycPyr_PB.fasta.out -o lycPyr_PB.fasta.out.bed

grep 'Satellite' lycPyr_PB.fasta.out.bed > lycPyr_PB_sats.bed
ml bioinfo-tools BEDTools R_packages/3.6.0
bedtools sort -i lycPyr_PB_sats.bed > temp

#awk '$5 <= 10 {print $0}' temp > lycPyr_PB_sats_div10.bed
#bedtools merge -d 100 -c 4 -o collapse,count -i lycPyr_PB_sats_div10_rename.bed > lycPyr_PB_sats_div10_rename_long.bed


awk '{sub("_.*", "", $4); print $0}' OFS="\t" lycPyr_PB_sats.bed | bedtools sort -i stdin > lycPyr_PB_sats_rename.bed

bedtools merge -d 100 -c 4 -o collapse,count -i lycPyr_PB_sats_rename.bed | awk '{print $0, $3-$2+1}' OFS="\t" > lycPyr_PB_sats_rename_merged

awk '$5 > 1 {print $0}' lycPyr_PB_sats_rename_merged > lycPyr_PB_sats_rename_merged_long

Rscript --vanilla $CODE/array_rmsk_assemblies.R lycPyr_PB_sats_rename_merged_long

awk '$7 > 10000 {print $0}' total_arrays_composition_orientation.txt | cut -f 5 | sort | uniq -c | tr -s ' ' | sed 's/^ //' | awk '{print $2, $1}' OFS="\t" | sort -r -n -k2,2 > array_types.txt




SAMPLE=ptiInt

/proj/sllstore2017073/private/scripts/RM2BED.sh -i ptiInt_PB.fasta.out -o ptiInt_PB.fasta.out.bed

grep 'Satellite' ${SAMPLE}_PB.fasta.out.bed > ${SAMPLE}_PB_sats.bed
ml bioinfo-tools BEDTools R_packages/3.6.0
bedtools sort -i ${SAMPLE}_PB_sats.bed > temp

#awk '$5 <= 10 {print $0}' temp > ${SAMPLE}_PB_sats_div10.bed
#bedtools merge -d 100 -c 4 -o collapse,count -i ${SAMPLE}_PB_sats_div10_rename.bed > ${SAMPLE}_PB_sats_div10_rename_long.bed


awk '{sub("_.*", "", $4); print $0}' OFS="\t" ${SAMPLE}_PB_sats.bed | bedtools sort -i stdin > ${SAMPLE}_PB_sats_rename.bed

bedtools merge -d 100 -c 4 -o collapse,count -i ${SAMPLE}_PB_sats_rename.bed | awk '{print $0, $3-$2+1}' OFS="\t" > ${SAMPLE}_PB_sats_rename_merged

awk '$5 > 1 {print $0}' ${SAMPLE}_PB_sats_rename_merged > ${SAMPLE}_PB_sats_rename_merged_long

Rscript --vanilla $CODE/array_rmsk_assemblies.R ${SAMPLE}_PB_sats_rename_merged_long

awk '$7 > 10000 {print $0}' total_arrays_composition_orientation.txt | cut -f 5 | sort | uniq -c | tr -s ' ' | sed 's/^ //' | awk '{print $2, $1}' OFS="\t" | sort -r -n -k2,2 > array_types.txt







for FILE in $( ls corCor* )
do
 OUT=`echo $FILE | sed 's/corCor/corCon/'`
 mv $FILE $OUT
done


