#!/bin/bash -l 
#SBATCH -J SAMPLES_download
#SBATCH -o SAMPLES_download.output
#SBATCH -e SAMPLES_download.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL                           
#SBATCH -t 12:00:00
#SBATCH -A snic2020-5-680
#SBATCH -p core
#SBATCH -n 1


# Download all the Corvus Illumina libraries needed

DIR=/proj/sllstore2017073/private/Valentina/2021SatCorvides/Data/Illumina
mkdir -p $DIR && cd $DIR

# Corvus hawaiensis
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR345/005/ERR3458785/ERR3458785_1.fastq.gz -O corHaw_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR345/005/ERR3458785/ERR3458785_2.fastq.gz -O corHaw_R2.fastq.gz

# Corvus brachyrhynchos
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR944/SRR944630/SRR944630_1.fastq.gz -O corBra_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR944/SRR944630/SRR944630_2.fastq.gz -O corBra_R2.fastq.gz

# Corvus monedula / Coloeus monedula
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR136/006/SRR1361906/SRR1361906_1.fastq.gz -O colMon_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR136/006/SRR1361906/SRR1361906_2.fastq.gz -O colMon_R2.fastq.gz

# Corvus kubaryi
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR345/006/ERR3458786/ERR3458786_1.fastq.gz -O corKub_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR345/006/ERR3458786/ERR3458786_2.fastq.gz -O corKub_R2.fastq.gz

# Corvus corax
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR345/003/ERR3458783/ERR3458783_1.fastq.gz -O corCax_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR345/003/ERR3458783/ERR3458783_2.fastq.gz -O corCax_R2.fastq.gz

# Corvus dauuricus
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR850/ERR850077/ERR850077_1.fastq.gz -O corDau_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR850/ERR850077/ERR850077_2.fastq.gz -O corDau_R2.fastq.gz

# Corvus frugileus
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR136/009/SRR1368309/SRR1368309_1.fastq.gz -O corFru_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR136/009/SRR1368309/SRR1368309_2.fastq.gz -O corFru_R2.fastq.gz

# Corvus moneduloides
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR345/004/ERR3457574/ERR3457574_1.fastq.gz -O corMon_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR345/004/ERR3457574/ERR3457574_2.fastq.gz -O corMon_R2.fastq.gz

# Corvus pectoralis
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR290/003/ERR2900313/ERR2900313_1.fastq.gz -O corPec_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR290/003/ERR2900313/ERR2900313_2.fastq.gz -O corPec_R2.fastq.gz

# Corvus woodfordi
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR345/000/ERR3457550/ERR3457550_1.fastq.gz -O corWoo_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR345/000/ERR3457550/ERR3457550_2.fastq.gz -O corWoo_R2.fastq.gz

# Corvus splendens
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR345/008/ERR3457558/ERR3457558_1.fastq.gz -O corSpl_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR345/008/ERR3457558/ERR3457558_2.fastq.gz -O corSpl_R2.fastq.gz

# Corvus cryptoleucus
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR633/000/SRR6337630/SRR6337630.fastq.gz -O corCry_R.fastq.gz

# Corvus tasmanicus
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR345/004/ERR3458784/ERR3458784_1.fastq.gz -O corTas_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR345/004/ERR3458784/ERR3458784_2.fastq.gz -O corTas_R2.fastq.gz

# Corvus macrorhynchos
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR250/DRR250114/DRR250114_1.fastq.gz -O corMac_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR250/DRR250114/DRR250114_2.fastq.gz -O corMac_R2.fastq.gz


# corvus corone and corvus cornix are to be downloaded yet; i?m confused about the samples
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/001/SRR1284441/SRR1284441_1.fastq.gz -O corCor_1_R1.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/009/SRR1284299/SRR1284299_1.fastq.gz -O corCor_2_R1.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/000/SRR1284320/SRR1284320_1.fastq.gz -O corCor_3_R1.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR210/008/SRR2107328/SRR2107328_1.fastq.gz -O corCor_RNA_R1.fastq.gz
		
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/000/SRR1265090/SRR1265090_1.fastq.gz -O corCon_1_R1.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/000/SRR1265100/SRR1265100_1.fastq.gz -O corCon_2_R1.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/005/SRR1265105/SRR1265105_1.fastq.gz -O corCon_3_R1.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR210/003/SRR2107373/SRR2107373_1.fastq.gz -O corCon_RNA_R1.fastq.gz
		
		
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR290/005/ERR2900305/ERR2900305_1.fastq.gz -O hyb_M1_R1.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR290/007/ERR2900307/ERR2900307_1.fastq.gz -O hyb_M2_R1.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR290/003/ERR2900303/ERR2900303_1.fastq.gz -O hyb_F1_R1.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR290/006/ERR2900306/ERR2900306_1.fastq.gz -O hyb_F2_R1.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR290/002/ERR2900302/ERR2900302_1.fastq.gz -O hyb_M3_R1.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR290/004/ERR2900304/ERR2900304_1.fastq.gz -O hyb_F3_R1.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/001/SRR1284441/SRR1284441_2.fastq.gz -O corCor_1_R2.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/009/SRR1284299/SRR1284299_2.fastq.gz -O corCor_2_R2.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/000/SRR1284320/SRR1284320_2.fastq.gz -O corCor_3_R2.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR210/008/SRR2107328/SRR2107328_2.fastq.gz -O corCor_RNA_R2.fastq.gz
		
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/000/SRR1265090/SRR1265090_2.fastq.gz -O corCon_1_R2.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/000/SRR1265100/SRR1265100_2.fastq.gz -O corCon_2_R2.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/005/SRR1265105/SRR1265105_2.fastq.gz -O corCon_3_R2.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR210/003/SRR2107373/SRR2107373_2.fastq.gz -O corCon_RNA_R2.fastq.gz
		
		
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR290/005/ERR2900305/ERR2900305_2.fastq.gz -O hyb_M1_R2.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR290/007/ERR2900307/ERR2900307_2.fastq.gz -O hyb_M2_R2.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR290/003/ERR2900303/ERR2900303_2.fastq.gz -O hyb_F1_R2.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR290/006/ERR2900306/ERR2900306_2.fastq.gz -O hyb_F2_R2.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR290/002/ERR2900302/ERR2900302_2.fastq.gz -O hyb_M3_R2.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR290/004/ERR2900304/ERR2900304_2.fastq.gz -O hyb_F3_R2.fastq.gz
