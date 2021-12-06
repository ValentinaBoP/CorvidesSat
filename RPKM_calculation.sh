# RPKM calculation for Corvus cornix, Corvus corone and Lycocorax pyrrhopterus

#!/bin/bash -l 
#SBATCH -J RPKM_RNAseq_sat_corCon
#SBATCH -o RPKM_RNAseq_sat_corCon.output
#SBATCH -e RPKM_RNAseq_sat_corCon.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL                           
#SBATCH -t 00:15:00
#SBATCH -A snic2020-5-680
#SBATCH -p core
#SBATCH -n 4

# Calculate RPKM for each satellite families in lycPyr
# FOLLOW FROM BWA_RNAseq_lycPyrOnly_sat.sh

#modules
ml bioinfo-tools samtools/1.12 R_packages/3.6.0

WD=/proj/sllstore2017073/private/Valentina/2021SatCorvides/Intermediate/Transcription/corCon_lib
cd $WD

CODE=/proj/sllstore2017073/private/Valentina/2021SatCorvides/Code

# STEP1: count the total reads in the library and divide by 1,000,000 – this is our “per million” scaling factor
R1=/home/vpeona/sllstore2017073/private/Valentina/2021SatCorvides/Data/Illumina/Corvus/corCon_RNA_R1.fastq.gz
R2=/home/vpeona/sllstore2017073/private/Valentina/2021SatCorvides/Data/Illumina/Corvus/corCon_RNA_R2.fastq.gz

# read contained in half the library
zcat $R1 | grep '>' | wc -l | awk '{print $1/2}' > tot_reads
# not going to count the number of reads in R2 (assume that R1 and R2 are the same size -- they actually are)

FACTOR=`awk '{print $1/1000000}' tot_reads`

# STEP2: divide the read counts by the scaling factor
# STEP2.1: count the reads for each satellite family from the bam file with samtools view
BAM=/proj/sllstore2017073/private/Valentina/2021SatCorvides/Intermediate/Transcription/corCon_lib/corCon.bam
cp /proj/sllstore2017073/private/Valentina/2021SatCorvides/Results/corvides_sat1.4_dimers.fasta .

LIB=corvides_sat1.4_dimers.fasta

grep 'corCon\|crowSat' $LIB | cut -c2- > list_seq
#printf "bopSat1_dimer#Satellite\nbopSat2_dimer#Satellite\nbopSat3_dimer#Satellite\n" >> list_seq
perl /proj/sllstore2017073/private/scripts/extractFromFasta.pl $LIB list list_seq > temp 
mv temp $LIB

samtools faidx $LIB
#cat $LIB.fai | awk '{print $1, "0", $2}' OFS="\t" > regions.bed
cat $LIB.fai | awk '{print $1 ":0-" $2}' > regions

samtools view -h -b -o temp -@ 4 -q 30 $BAM
BAM=${BAM%.*}_filterQ30.bam
mv temp $BAM
samtools index -b $BAM

for REGION in $( cat regions )
do
 samtools view -c -o temp -@ 4 $BAM $REGION
 NUM=$(cat temp)
 printf "$REGION\t$NUM\n" >> summary_temp
done

cat summary_temp | tr -s ':' '\t' | tr -s '-' '\t' > corvides_sat1.4_dimers.counts

rm summary_temp
rm temp

# STEP2.2: divide the read count for each satellite family by the scaling factor
# reads per million (RPM)
awk -v factor=$FACTOR '{print $0, $3/1000, $4/factor}' OFS="\t" corvides_sat1.4_dimers.counts > corvides_sat1.4_dimers.counts.rpm

# STEP3: Divide the RPM values by the length of the gene, in kilobases. This gives you RPKM.
awk '{print $0, $6/$5}' OFS="\t" corvides_sat1.4_dimers.counts.rpm > corvides_sat1.4_dimers.counts.rpkm

printf "Satellite\tStart\tEnd\tReadCount\tSizeKb\trpm\trpkm\n" > header

cat header corvides_sat1.4_dimers.counts.rpkm > temp
mv temp corvides_sat1.4_dimers.counts.rpkm

# STEP4: aggregate by family in R
Rscript --vanilla $CODE/filterFamiliesRPKM.R corvides_sat1.4_dimers.counts.rpkm

# STEP5: barplot of the rpkm values for each satellite family
Rscript --vanilla $CODE/plotFamiliesRPKM.R corvides_sat1.4_dimers.counts.rpkm

### LYCPYR
#!/bin/bash -l 
#SBATCH -J RPKM_RNAseq_sat
#SBATCH -o RPKM_RNAseq_sat.output
#SBATCH -e RPKM_RNAseq_sat.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL                           
#SBATCH -t 00:15:00
#SBATCH -A snic2020-5-680
#SBATCH -p core
#SBATCH -n 4

# Calculate RPKM for each satellite families in lycPyr
# FOLLOW FROM BWA_RNAseq_lycPyrOnly_sat.sh

#modules
ml bioinfo-tools samtools/1.12 R_packages/3.6.0

WD=/proj/sllstore2017073/private/Valentina/2021SatCorvides/Intermediate/Transcription/lycPyr_lib
cd $WD

CODE=/proj/sllstore2017073/private/Valentina/2021SatCorvides/Code

# STEP1: count the total reads in the library and divide by 1,000,000 – this is our “per million” scaling factor
R1=/proj/sllstore2017073/private/01_raw_data/SH-2274/191218_A00181_0135_AHHKFMDRXX/Sample_SH-2274-THY-410-25-11-13/SH-2274-THY-410-25-11-13_S9_L001_R1_001.fastq.gz
R2=/proj/sllstore2017073/private/01_raw_data/SH-2274/191218_A00181_0135_AHHKFMDRXX/Sample_SH-2274-THY-410-25-11-13/SH-2274-THY-410-25-11-13_S9_L001_R2_001.fastq.gz

# read contained in half the library
zcat $R1 | grep '>' | wc -l | awk '{print $1/2}' > tot_reads
# not going to count the number of reads in R2 (assume that R1 and R2 are the same size -- they actually are)

FACTOR=`awk '{print $1/1000000}' tot_reads`

# STEP2: divide the read counts by the scaling factor
# STEP2.1: count the reads for each satellite family from the bam file with samtools view
BAM=/proj/sllstore2017073/private/Valentina/2021SatCorvides/Intermediate/Transcription/lycPyr.bam
cp /proj/sllstore2017073/private/Valentina/2021SatCorvides/Results/corvides_sat1.3_dimers.fasta .

LIB=corvides_sat1.3_dimers.fasta
grep 'lycPyr' $LIB | cut -c2- > list_seq
printf "bopSat1_dimer#Satellite\nbopSat2_dimer#Satellite\nbopSat3_dimer#Satellite\n" >> list_seq
perl /proj/sllstore2017073/private/scripts/extractFromFasta.pl $LIB list list_seq > temp 
mv temp $LIB

samtools faidx $LIB
#cat $LIB.fai | awk '{print $1, "0", $2}' OFS="\t" > regions.bed
cat $LIB.fai | awk '{print $1 ":0-" $2}' > regions

for REGION in $( cat regions )
do
 samtools view -c -o temp -@ 4 $BAM $REGION
 NUM=$(cat temp)
 printf "$REGION\t$NUM\n" >> summary_temp
done

cat summary_temp | tr -s ':' '\t' | tr -s '-' '\t' > corvides_sat1.3_dimers.counts

rm summary_temp
rm temp

# STEP2.2: divide the read count for each satellite family by the scaling factor
# reads per million (RPM)
awk -v factor=$FACTOR '{print $0, $3/1000, $4/factor}' OFS="\t" corvides_sat1.3_dimers.counts > corvides_sat1.3_dimers.counts.rpm

# STEP3: Divide the RPM values by the length of the gene, in kilobases. This gives you RPKM.
awk '{print $0, $6/$5}' OFS="\t" corvides_sat1.3_dimers.counts.rpm > corvides_sat1.3_dimers.counts.rpkm

printf "Satellite\tStart\tEnd\tReadCount\tSizeKb\trpm\trpkm\n" > header

cat header corvides_sat1.3_dimers.counts.rpkm > temp
mv temp corvides_sat1.3_dimers.counts.rpkm

# STEP4: aggregate by family in R
Rscript --vanilla $CODE/filterFamiliesRPKM.R corvides_sat1.3_dimers.counts.rpkm

# STEP5: barplot of the rpkm values for each satellite family
Rscript --vanilla $CODE/plotFamiliesRPKM.R corvides_sat1.3_dimers.counts.rpkm

## CORCOR
#!/bin/bash -l 
#SBATCH -J RPKM_RNAseq_sat_corCor
#SBATCH -o RPKM_RNAseq_sat_corCor.output
#SBATCH -e RPKM_RNAseq_sat_corCor.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL                           
#SBATCH -t 00:15:00
#SBATCH -A snic2020-5-680
#SBATCH -p core
#SBATCH -n 4

# Calculate RPKM for each satellite families in lycPyr
# FOLLOW FROM BWA_RNAseq_lycPyrOnly_sat.sh

#modules
ml bioinfo-tools samtools/1.12 R_packages/3.6.0

WD=/proj/sllstore2017073/private/Valentina/2021SatCorvides/Intermediate/Transcription/corCor_lib
cd $WD

CODE=/proj/sllstore2017073/private/Valentina/2021SatCorvides/Code

# STEP1: count the total reads in the library and divide by 1,000,000 – this is our “per million” scaling factor
R1=/home/vpeona/sllstore2017073/private/Valentina/2021SatCorvides/Data/Illumina/Corvus/corCor_RNA_R1.fastq.gz
R2=/home/vpeona/sllstore2017073/private/Valentina/2021SatCorvides/Data/Illumina/Corvus/corCor_RNA_R2.fastq.gz

# read contained in half the library
zcat $R1 | grep '>' | wc -l | awk '{print $1/2}' > tot_reads
# not going to count the number of reads in R2 (assume that R1 and R2 are the same size -- they actually are)

FACTOR=`awk '{print $1/1000000}' tot_reads`

# STEP2: divide the read counts by the scaling factor
# STEP2.1: count the reads for each satellite family from the bam file with samtools view
BAM=/proj/sllstore2017073/private/Valentina/2021SatCorvides/Intermediate/Transcription/corCor_lib/corCor.bam
cp /proj/sllstore2017073/private/Valentina/2021SatCorvides/Results/corvides_sat1.4_dimers.fasta .

LIB=corvides_sat1.4_dimers.fasta
grep 'corCor\|crowSat' $LIB | cut -c2- > list_seq
#printf "bopSat1_dimer#Satellite\nbopSat2_dimer#Satellite\nbopSat3_dimer#Satellite\n" >> list_seq
perl /proj/sllstore2017073/private/scripts/extractFromFasta.pl $LIB list list_seq > temp 
mv temp $LIB

samtools faidx $LIB
#cat $LIB.fai | awk '{print $1, "0", $2}' OFS="\t" > regions.bed
cat $LIB.fai | awk '{print $1 ":0-" $2}' > regions

samtools view -h -b -o temp -@ 4 -q 30 $BAM
BAM=${BAM%.*}_filterQ30.bam
mv temp $BAM
samtools index -b $BAM

for REGION in $( cat regions )
do
 samtools view -c -o temp -@ 4 $BAM $REGION
 NUM=$(cat temp)
 printf "$REGION\t$NUM\n" >> summary_temp
done

cat summary_temp | tr -s ':' '\t' | tr -s '-' '\t' > corvides_sat1.4_dimers.counts

rm summary_temp
rm temp

# STEP2.2: divide the read count for each satellite family by the scaling factor
# reads per million (RPM)
awk -v factor=$FACTOR '{print $0, $3/1000, $4/factor}' OFS="\t" corvides_sat1.4_dimers.counts > corvides_sat1.4_dimers.counts.rpm

# STEP3: Divide the RPM values by the length of the gene, in kilobases. This gives you RPKM.
awk '{print $0, $6/$5}' OFS="\t" corvides_sat1.4_dimers.counts.rpm > corvides_sat1.4_dimers.counts.rpkm

printf "Satellite\tStart\tEnd\tReadCount\tSizeKb\trpm\trpkm\n" > header

cat header corvides_sat1.4_dimers.counts.rpkm > temp
mv temp corvides_sat1.4_dimers.counts.rpkm

# STEP4: aggregate by family in R
Rscript --vanilla $CODE/filterFamiliesRPKM.R corvides_sat1.4_dimers.counts.rpkm

# STEP5: barplot of the rpkm values for each satellite family
Rscript --vanilla $CODE/plotFamiliesRPKM.R corvides_sat1.4_dimers.counts.rpkm
