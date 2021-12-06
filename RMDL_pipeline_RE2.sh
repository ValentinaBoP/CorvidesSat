## RepeatModeler pipeline to align genomic hits to the consensus sequences of the satDNAs produced by RepeatExplorer2 (on Galaxy)
## job example for Astrapia rothschildi but it was done for all the species for which there was an assembly available

#!/bin/bash
#SBATCH -J RMDL_pipeline_RE2_astRot
#SBATCH -o RMDL_pipeline_RE2_astRot.output
#SBATCH -e RMDL_pipeline_RE2_astRot.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 02:00:00
#SBATCH -A snic2020-5-680
#SBATCH -p core
#SBATCH -n 2

SAMPLE=astRot

DIR=/proj/sllstore2017073/private/Valentina/2021SatCorvides/Intermediate/RE2/RMDL/$SAMPLE
mkdir -p $DIR
cd $DIR

REF=/home/vpeona/sllstore2017073/private/Valentina/2021SatCorvides/Data/Genomes/${SAMPLE}.fasta
CONS=/home/vpeona/sllstore2017073/private/Valentina/2021SatCorvides/Intermediate/RE2/Consensi/${SAMPLE}_RE2_satellite.fasta

ml purge
ml bioinfo-tools blast/2.9.0+ MAFFT/7.310

unset $MAFFT_BINARIES

makeblastdb -in $REF -out $REF -dbtype nucl -parse_seqids

awk '{print$1;}' $CONS > temp

CONS=`basename $CONS`
CONS=${CONS%.*}.clean.fasta
mv temp $CONS

sed -i 's/\//_/g' $CONS

perl /proj/sllstore2017073/private/scripts/repeatModelerPipeline4.pl $REF $REF $CONS

module purge

module load bioinfo-tools T-Coffee/11.00.8cbe486

export CACHE_4_TCOFFEE=/home/vpeona/bin/

rm *emp.out

mkdir final

cd aligned; for i in $(ls *.fa); do name=`ls $i | cut -f1 -d "."`; cat $i | perl -ne 'chomp;s/>\s+/>/;if(/>(\S+)/){$id{$1}++;$id2=$1;}if($id{$id2}==1){print "$_\n"}' >../final/$name.fa; done; cd ../

cd final; for i in $(ls *.fa); do name=`ls $i | cut -f1 -d "."`; t_coffee -other_pg seq_reformat -in $i -action +rm_gap 95 >$name.gaps95.fa; done; cd ../



cd /home/vpeona/sllstore2017073/private/Valentina/2021SatCorvides/Code/RMDL_pipeline_jobs

SAMPLE=astRot

for SUB in dreAlb manKer parBre parLaw parRub ptiInt ptiMag
do
 sed "s/${SAMPLE}/${SUB}/g" RMDL_pipeline_RE2_${SAMPLE}.sh > RMDL_pipeline_RE2_${SUB}.sh
done

for JOB in $(ls RMSK*.sh )
do
 chmod +x $JOB
 sbatch $JOB
done
