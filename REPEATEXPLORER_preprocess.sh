#!/bin/bash -l 
#SBATCH -J REPEATEXPLORER2_BOP_preprocess
#SBATCH -o REPEATEXPLORER2_BOP_preprocess.output
#SBATCH -e REPEATEXPLORER2_BOP_preprocess.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL                           
#SBATCH -t 12:00:00
#SBATCH -A snic2020-5-680
#SBATCH -p core
#SBATCH -n 10

# modules
ml R_packages/3.6.0
ml bioinfo-tools
ml seqtk/1.2-r101
ml cutadapt/2.3

# working directory
DIR=/home/vpeona/sllstore2017073/private/Valentina/2021SatCorvides/Data/Illumina/BOPs/
mkdir -p $DIR && cd $DIR
mkdir Filtered

# manCha female
ln -s /proj/sllstore2017073/private/01_raw_data/A.Suh_17_02/fastq/HGYLKALXX/outs/fastq_path/A_Suh_17_02/P7359_101/P7359_101_S1_L006_R1_001.fastq.gz ./manCha_R1.fastq.gz
ln -s /proj/sllstore2017073/private/01_raw_data/A.Suh_17_02/fastq/HGYLKALXX/outs/fastq_path/A_Suh_17_02/P7359_101/P7359_101_S1_L006_R2_001.fastq.gz ./manCha_R2.fastq.gz

# parHel female
ln -s /proj/sllstore2017073/private/01_raw_data/A.Suh_18_01/P9957/P9957_102/02-FASTQ/180302_ST-E00214_0206_AHJ5WKCCXY/P9957_102_S2_L002_R1_001.fastq.gz ./parHel_R1.fastq.gz
ln -s /proj/sllstore2017073/private/01_raw_data/A.Suh_18_01/P9957/P9957_102/02-FASTQ/180302_ST-E00214_0206_AHJ5WKCCXY/P9957_102_S2_L002_R2_001.fastq.gz ./parHel_R2.fastq.gz

# cinMag female
ln -s /proj/sllstore2017073/private/01_raw_data/A.Suh_19_10/P13554/P13554_109/02-FASTQ/191004_A00187_0202_BHNF7GDSXX/P13554_109_S2_L003_R1_001.fastq.gz ./cinMag_R1.fastq.gz
ln -s /proj/sllstore2017073/private/01_raw_data/A.Suh_19_10/P13554/P13554_109/02-FASTQ/191004_A00187_0202_BHNF7GDSXX/P13554_109_S2_L003_R2_001.fastq.gz ./cinMag_R2.fastq.gz

# cinReg female
ln -s /proj/sllstore2017073/private/01_raw_data/A.Suh_19_10/P13554/P13554_114/02-FASTQ/190904_A00621_0120_BHMFMYDSXX/P13554_114_S12_L003_R1_001.fastq.gz ./cinReg_R1.fastq.gz
ln -s /proj/sllstore2017073/private/01_raw_data/A.Suh_19_10/P13554/P13554_114/02-FASTQ/190904_A00621_0120_BHMFMYDSXX/P13554_114_S12_L003_R2_001.fastq.gz ./cinReg_R2.fastq.gz

# epiMey female
ln -s /proj/sllstore2017073/private/01_raw_data/A.Suh_17_02/fastq/HGYLKALXX/outs/fastq_path/A_Suh_17_02/P7359_102/P7359_102_S2_L006_R1_001.fastq.gz ./epiMey_R1.fastq.gz
ln -s /proj/sllstore2017073/private/01_raw_data/A.Suh_17_02/fastq/HGYLKALXX/outs/fastq_path/A_Suh_17_02/P7359_102/P7359_102_S2_L006_R2_001.fastq.gz ./epiMey_R2.fastq.gz

# lycPyr102 female reference
ln -s /proj/sllstore2017073/private/BirdsOfParadise/lycPyr_10XG/H5JJHALXX/outs/fastq_path/A_Suh_16_01/Sample_P6156_102/P6156_102_S8_L008_R1_001.fastq.gz ./lycPyr_R1.fastq.gz
ln -s /proj/sllstore2017073/private/BirdsOfParadise/lycPyr_10XG/H5JJHALXX/outs/fastq_path/A_Suh_16_01/Sample_P6156_102/P6156_102_S8_L008_R2_001.fastq.gz ./lycPyr_R2.fastq.gz

# astRot male 10X
ln -s /crex/proj/sllstore2017073/private/01_raw_data/A.Suh_19_10/P13554/P13554_105/02-FASTQ/190904_A00621_0120_BHMFMYDSXX/P13554_105_S5_L002_R1_001.fastq.gz ./astRot_R1.fastq.gz
ln -s /crex/proj/sllstore2017073/private/01_raw_data/A.Suh_19_10/P13554/P13554_105/02-FASTQ/190904_A00621_0120_BHMFMYDSXX/P13554_105_S5_L002_R2_001.fastq.gz ./astRot_R2.fastq.gz

# dreAlb male
ln -s /crex/proj/sllstore2017073/private/01_raw_data/A.Suh_19_10/P13554/P13554_107/02-FASTQ/190904_A00621_0120_BHMFMYDSXX/P13554_107_S7_L002_R1_001.fastq.gz ./dreAlb_R1.fastq.gz
ln -s /crex/proj/sllstore2017073/private/01_raw_data/A.Suh_19_10/P13554/P13554_107/02-FASTQ/190904_A00621_0120_BHMFMYDSXX/P13554_107_S7_L002_R2_001.fastq.gz ./dreAlb_R2.fastq.gz

# parBre male
ln -s /crex/proj/sllstore2017073/private/01_raw_data/A.Suh_19_10/P13554/P13554_106/02-FASTQ/190904_A00621_0120_BHMFMYDSXX/P13554_106_S6_L002_R1_001.fastq.gz ./parBre_R1.fastq.gz
ln -s /crex/proj/sllstore2017073/private/01_raw_data/A.Suh_19_10/P13554/P13554_106/02-FASTQ/190904_A00621_0120_BHMFMYDSXX/P13554_106_S6_L002_R2_001.fastq.gz ./parBre_R2.fastq.gz

# parRub female
ln -s /crex/proj/sllstore2017073/private/01_raw_data/A.Suh_19_09/P13361/P13361_134/02-FASTQ/190920_A00621_0128_BHMH2JDSXX/P13361_134_S9_L002_R1_001.fastq.gz ./parRub_R1.fastq.gz
ln -s /crex/proj/sllstore2017073/private/01_raw_data/A.Suh_19_09/P13361/P13361_134/02-FASTQ/190920_A00621_0128_BHMH2JDSXX/P13361_134_S9_L002_R2_001.fastq.gz ./parRub_R2.fastq.gz

# ptiInt female
ln -s /crex/proj/sllstore2017073/private/01_raw_data/A.Suh_18_01/P9957/P9957_101/02-FASTQ/180302_ST-E00214_0206_AHJ5WKCCXY/P9957_101_S1_L002_R1_001.fastq.gz ./ptiInt_R1.fastq.gz
ln -s /crex/proj/sllstore2017073/private/01_raw_data/A.Suh_18_01/P9957/P9957_101/02-FASTQ/180302_ST-E00214_0206_AHJ5WKCCXY/P9957_101_S1_L002_R2_001.fastq.gz ./ptiInt_R2.fastq.gz

# ptiMag female
ln -s /crex/proj/sllstore2017073/private/01_raw_data/A.Suh_19_10/P13554/P13554_108/02-FASTQ/191004_A00187_0202_BHNF7GDSXX/P13554_108_S1_L003_R1_001.fastq.gz ./ptiMag_R1.fastq.gz
ln -s /crex/proj/sllstore2017073/private/01_raw_data/A.Suh_19_10/P13554/P13554_108/02-FASTQ/191004_A00187_0202_BHNF7GDSXX/P13554_108_S1_L003_R2_001.fastq.gz ./ptiMag_R2.fastq.gz

# parLaw female
ln -s /crex/proj/sllstore2017073/private/01_raw_data/A.Suh_19_10/P13554/P13554_116/02-FASTQ/190904_A00621_0120_BHMFMYDSXX/P13554_116_S14_L003_R1_001.fastq.gz ./parLaw_R1.fastq.gz
ln -s /crex/proj/sllstore2017073/private/01_raw_data/A.Suh_19_10/P13554/P13554_116/02-FASTQ/190904_A00621_0120_BHMFMYDSXX/P13554_116_S14_L003_R2_001.fastq.gz ./parLaw_R2.fastq.gz

# manKer female
ln -s /crex/proj/sllstore2017073/private/01_raw_data/A.Suh_19_10/P13554/P13554_115/02-FASTQ/190904_A00621_0120_BHMFMYDSXX/P13554_115_S13_L003_R1_001.fastq.gz ./manKer_R1.fastq.gz
ln -s /crex/proj/sllstore2017073/private/01_raw_data/A.Suh_19_10/P13554/P13554_115/02-FASTQ/190904_A00621_0120_BHMFMYDSXX/P13554_115_S13_L003_R2_001.fastq.gz ./manKer_R2.fastq.gz

# astSte female
ln -s /proj/sllstore2017055/private/181015_BOP/1.clean_reads/5.quality_trimmed/AsteF677/AsteF677_L001/AsteF677_L001_R1.fastq.gz ./astSte_R1.fastq.gz
ln -s /proj/sllstore2017055/private/181015_BOP/1.clean_reads/5.quality_trimmed/AsteF677/AsteF677_L001/AsteF677_L001_R2.fastq.gz ./astSte_R2.fastq.gz

# parRag male
ln -s /proj/sllstore2017055/private/181015_BOP/1.clean_reads/5.quality_trimmed/PragM588/PragM588_L005/PragM588_L005_R1.fastq.gz ./parRag_R1.fastq.gz
ln -s /proj/sllstore2017055/private/181015_BOP/1.clean_reads/5.quality_trimmed/PragM588/PragM588_L005/PragM588_L005_R2.fastq.gz ./parRag_R2.fastq.gz

# go to temporary directory
cp *_R1.fastq.gz $SNIC_TMP
cp *_R2.fastq.gz $SNIC_TMP
cd $SNIC_TMP

ls *_R1.fastq.gz > list_fastq

NUMREADS=1000000

for R1 in $( cat list_fastq )
do
 SAMPLE=`echo ${R1%.*.*} | sed 's/_R1//'`
 #ID=`echo $SAMPLE | cut -c 1,4,7,8`
 R2=${SAMPLE}_R2.fastq.gz

 OUT1=${SAMPLE}_R1_filtered.fasta
 OUT2=${SAMPLE}_R2_filtered.fasta
 /home/vpeona/re_utilities/paired_fastq_filtering.R -a $R1 -b $R2 -x $OUT1 -y $OUT2 -R -f fasta
 #  -s 123 -n $NUMREADS
 
 seqtk sample -s 123 $OUT1 $NUMREADS > ${SAMPLE}_R1_filtered_sampled.fasta
 seqtk sample -s 123 $OUT2 $NUMREADS > ${SAMPLE}_R2_filtered_sampled.fasta

 # Interlace reads
 /home/vpeona/re_utilities/fasta_interlacer.py -a ${SAMPLE}_R1_filtered_sampled.fasta -b ${SAMPLE}_R2_filtered_sampled.fasta -p ${SAMPLE}_INTERLACED.fasta -x ${SAMPLE}_NOPAIR.fasta

 # Affix identifier to the reads
 /home/vpeona/re_utilities/fasta_affixer.py -f ${SAMPLE}_INTERLACED.fasta -o ${SAMPLE}_RE2_ready.fasta -p $SAMPLE

 cp ${SAMPLE}_RE2_ready.fasta $DIR/Filtered
done

cat *_RE2_ready.fasta > BOPS_RE2_ready.fasta

cp BOPS_RE2_ready.fasta $DIR/Filtered
