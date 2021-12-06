#!/bin/bash -l
#SBATCH -J RMSK_LANDSCAPE_corvide
#SBATCH -o RMSK_LANDSCAPE_corvide.output
#SBATCH -e RMSK_LANDSCAPE_corvide.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 06:00:00
#SBATCH -A snic2020-5-680
#SBATCH -p core
#SBATCH -n 2

# produce the repeat landscape of Lycocorax assemblies using only the RepBase library

# load the module
module load bioinfo-tools RepeatMasker/4.0.7

# go to folder
WD=/proj/sllstore2017073/private/Valentina/2021SatCorvides/Intermediate/RE2/RMSK
cd $WD

ls -d */ | sed 's%/%%g' > list_dirs

for DIR in $( cat list_dirs )
do
 
 # go to directory 
 cd $DIR
 
 # Convert the .align file into a .tbl file with Kimura 2-parameter distances (excluding CpG sites)
 perl /proj/sllstore2017073/private/scripts/RM_divergence/calcDivergenceFromAlign.pl -noCpG ${DIR}.fasta.align > ${DIR}.k2p.noCpG

 # Add sizes of the individual TE copies
 perl /proj/sllstore2017073/private/scripts/RM_divergence/RM.addsize.div.pl ${DIR}.k2p.noCpG > ${DIR}.k2p.noCpG.size

 # Remove the term "kimura=" from the last column
 sed -i 's/kimura=//g' ${DIR}.k2p.noCpG.size

 cd $WD

done


# go to folder
WD=/proj/sllstore2017073/private/Valentina/2021SatCorvides/Intermediate/Abundance
cd $WD

for ALIGN in $( ls *.align )
do
 
 # Convert the .align file into a .tbl file with Kimura 2-parameter distances (excluding CpG sites)
 perl /proj/sllstore2017073/private/scripts/RM_divergence/calcDivergenceFromAlign.pl -noCpG $ALIGN > $ALIGN.k2p.noCpG

 cat $ALIGN.k2p.noCpG | tr -s ' ' | sed 's/^ //' | tr -s ' ' '\t' | sed 's/kimura=//' | awk '{print $0, $7-$6+1}' OFS="\t" | grep -v 'sat6_b_parRub' >$ALIGN.k2p.noCpG.size

 Rscript --vanilla $CODE/abundance_sat.R $ALIGN.k2p.noCpG.size
 cd $WD

done


ml R_packages/3.6.0
CODE=/home/vpeona/2021SatCorvides/Code
WD=/proj/sllstore2017073/private/Valentina/2021SatCorvides/Intermediate/RE2/RMSK
cd $WD
for DIR in $( cat list_dirs_all )
do
 
 # go to directory 
 cd $DIR
 
 # Convert the .align file into a .tbl file with Kimura 2-parameter distances (excluding CpG sites)
 perl /proj/sllstore2017073/private/scripts/RM_divergence/calcDivergenceFromAlign.pl -noCpG ${DIR}.fasta.align > ${DIR}.k2p.noCpG
 cat ${DIR}.k2p.noCpG | tr -s ' ' | sed 's/^ //' | tr -s ' ' '\t' | sed 's/kimura=//' | awk '{print $0, $7-$6+1}' OFS="\t"  > ${DIR}.k2p.noCpG.size

 # Add sizes of the individual TE copies
 Rscript --vanilla $CODE/abundance_sat.R ${DIR}.k2p.noCpG.size

 cd $WD

done



RepeatMasker -pa 16 -a -xsmall -gccalc -excln -lib $LIB ${SAMPLE}.fasta
cp -t $DIR/ ${SAMPLE}.fasta.align ${SAMPLE}.fasta.out ${SAMPLE}.fasta.tbl

cd $DIR
ml R_packages/3.6.0

perl /proj/sllstore2017073/private/scripts/RM_divergence/calcDivergenceFromAlign.pl -noCpG ${SAMPLE}.fasta.align > ${SAMPLE}.fasta.align.k2p.noCpG
cat ${SAMPLE}.fasta.align.k2p.noCpG | tr -s ' ' | sed 's/^ //' | tr -s ' ' '\t' | sed 's/kimura=//' | awk '{print $0, $7-$6+1}' OFS="\t" | grep -v 'sat6_b_parRub' > ${SAMPLE}.fasta.align.k2p.noCpG.size
Rscript --vanilla $CODE/abundance_sat.R ${SAMPLE}.fasta.align.k2p.noCpG.size







# for ALIGN in $( ls *.fasta.align )
# do
#   # Convert the .align file into a .tbl file with Kimura 2-parameter distances (excluding CpG sites)
#  perl /proj/sllstore2017073/private/scripts/RM_divergence/calcDivergenceFromAlign.pl -noCpG $ALIGN > $ALIGN.k2p.noCpG

#  # Add sizes of the individual TE copies
#  perl /proj/sllstore2017073/private/scripts/RM_divergence/RM.addsize.div.pl $ALIGN.k2p.noCpG > $ALIGN.k2p.noCpG.size

#  # Remove the term "kimura=" from the last column
#  sed -i 's/kimura=//g' $ALIGN.k2p.noCpG.size
# done

ml R_packages/3.6.0
CODE=/home/vpeona/2021SatCorvides/Code
WD=/home/vpeona/2021SatCorvides/Intermediate/lycPyr_PB/RMSK
cd $WD
for DIR in $( ls lycPyr_PB_reads_*.fasta.align )
do
 
 # Convert the .align file into a .tbl file with Kimura 2-parameter distances (excluding CpG sites)
 perl /proj/sllstore2017073/private/scripts/RM_divergence/calcDivergenceFromAlign.pl -noCpG $DIR > ${DIR}.k2p.noCpG
 cat ${DIR}.k2p.noCpG | tr -s ' ' | sed 's/^ //' | tr -s ' ' '\t' | sed 's/kimura=//' | awk '{print $0, $7-$6+1}' OFS="\t"  > ${DIR}.k2p.noCpG.size

 # Add sizes of the individual TE copies
 Rscript --vanilla $CODE/abundance_sat.R ${DIR}.k2p.noCpG.size
done



ml R_packages/3.6.0
CODE=/home/vpeona/2021SatCorvides/Code
WD=/home/vpeona/2021SatCorvides/Intermediate/ptiInt_PB/RMSK
cd $WD
for DIR in $( ls ptiIntPB_reads_*.fasta.align )
do
 
 # Convert the .align file into a .tbl file with Kimura 2-parameter distances (excluding CpG sites)
 perl /proj/sllstore2017073/private/scripts/RM_divergence/calcDivergenceFromAlign.pl -noCpG $DIR > ${DIR}.k2p.noCpG
 cat ${DIR}.k2p.noCpG | tr -s ' ' | sed 's/^ //' | tr -s ' ' '\t' | sed 's/kimura=//' | awk '{print $0, $7-$6+1}' OFS="\t"  > ${DIR}.k2p.noCpG.size

 # Add sizes of the individual TE copies
 Rscript --vanilla $CODE/abundance_sat.R ${DIR}.k2p.noCpG.size
done
