# RepeatMasker on reads
# example with Astrapia rothschildi 

#!/bin/bash
#SBATCH -J RMSK_sat_astRot
#SBATCH -o RMSK_sat_astRot.output
#SBATCH -e RMSK_sat_astRot.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 20:00:00
#SBATCH -A snic2019-8-371
#SBATCH -p core
#SBATCH -n 16

# load modules
ml bioinfo-tools RepeatMasker/4.0.8

# species ID
SAMPLE=astRot

# set output directory
DIR=/proj/sllstore2017073/private/Valentina/2021SatCorvides/Intermediate/RE2/RMSK/$SAMPLE
mkdir -p $DIR

# copy files into SNIC_TMP
cp /proj/sllstore2017073/private/Valentina/2021SatCorvides/Data/Illumina/BOPs/Filtered/${SAMPLE}_RE2_ready.fasta.gz $SNIC_TMP/${SAMPLE}.fasta.gz
cp /proj/sllstore2017073/private/Valentina/2021SatCorvides/Results/corvides_sat1.2_dimers.fasta $SNIC_TMP/corvides_sat1.2_dimers.lib

# move to SNIC_TMP
cd $SNIC_TMP

LIB=`ls *.lib`
gunzip ${SAMPLE}.fasta.gz

RepeatMasker -pa 16 -a -xsmall -gccalc -excln -lib $LIB ${SAMPLE}.fasta
cp -t $DIR/ ${SAMPLE}.fasta.align ${SAMPLE}.fasta.out ${SAMPLE}.fasta.tbl


cd /home/vpeona/sllstore2017073/private/Valentina/2021ChrW/Code/RMSK_sat_jobs

SAMPLE=astRot

for SUB in cinMag cinReg dreAlb epiMey astSte lycPyr manCha manKer parBre parLaw parRag parRub ptiInt ptiMag
do
 sed "s/${SAMPLE}/${SUB}/g" RMSK_sat_${SAMPLE}.sh > RMSK_sat_${SUB}.sh
done



for JOB in $(ls RMSK*.sh )
do
 chmod +x $JOB
 sbatch $JOB
done
