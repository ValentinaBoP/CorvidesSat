#!/bin/bash
#SBATCH -J RMSK_lycPyr_PB_reads_01
#SBATCH -o RMSK_lycPyr_PB_reads_01.output
#SBATCH -e RMSK_lycPyr_PB_reads_01.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 24:00:00
#SBATCH -A snic2019-8-371
#SBATCH -p core
#SBATCH -n 20

# load modules
ml bioinfo-tools RepeatMasker/4.1.0 seqtk

# set output directory
DIR=/proj/sllstore2017073/private/Valentina/2021SatCorvides/Intermediate/lycPyr_PB/RMSK
mkdir -p $DIR

# copy files to temporary directory
cp /proj/sllstore2017073/private/Valentina/2021SatCorvides/Results/corvides_sat1.3_dimers.fasta $SNIC_TMP/corvides_sat1.3_dimers.lib

cd $SNIC_TMP
LIB=`ls *.lib`

# link subreads
READS=lycPyr_filtered_subreads_01.fastq
ln -s /proj/sllstore2017073/private/BirdsOfParadise/lycPyr_PacBio/subreads/$READS .

seqtk seq -A $READS > temp

NUM=`echo $READS | sed 's/lycPyr_filtered_subreads_//' | sed 's/.fastq//'`
FASTA=`echo $READS | sed 's/filtered_sub/PB_/' | sed 's/fastq/fasta/'`
awk -v num=$NUM '{if($1 ~ />/){print ">batch" num "_" NR} else {print $0}}' temp > $FASTA

rm temp

RepeatMasker -pa 20 -a -xsmall -gccalc -excln -lib $LIB $FASTA
cp -t $DIR/ $FASTA.align $FASTA.out $FASTA.tbl $FASTA





for NUM in _02 _03 _04 _05 _06 _07 _08 _09 _10 _11 _12 _13 _14
do
 sed "s/_01/${NUM}/" RMSK_lycPyr_PB_reads_01.sh > RMSK_lycPyr_PB_reads${NUM}.sh
done
