# Working directory
mkdir -p /proj/sllstore2017073/private/Valentina/2021SatCorvides/Intermediate/Simulation/Input
cd /proj/sllstore2017073/private/Valentina/2021SatCorvides/Intermediate/Simulation/Input

# Collect genome and RepeatMasker annotation and sat library
ln -s /proj/sllstore2017073/private/Bird_references_2.0/assemblies/lycPyr7.4.fasta .
ln -s /proj/sllstore2017073/private/Valentina/2021SatCorvides/Intermediate/Abundance/lycPyr_PB_chr.annotation.light ./lycPyr7.4.fasta.out.bed
ln -s /proj/sllstore2017073/private/Valentina/2021SatCorvides/Results/corvides_sat1.4.fasta .


# Keep only chromosomes in the fasta
grep 'chr' lycPyr7.4.fasta > list_chr
perl /proj/sllstore2017073/private/scripts/extractFromFasta.pl lycPyr7.4.fasta list list_chr > lycPyr7.4_onlyChr.fasta

# hard mask with bedtools
ml bioinfo-tools BEDTools

bedtools maskfasta -fi lycPyr7.4_onlyChr.fasta -bed lycPyr7.4.fasta.out.bed -fo lycPyr7.4_onlyChr.fasta.masked

# remove all Ns from the genome
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < lycPyr7.4_onlyChr.fasta.masked | sed 's/N//g' | tr "\t" "\n" > lycPyr_clean.fasta

# generate fai for the new genome
ml samtools

samtools faidx lycPyr_clean.fasta



## THE SAT LIB MUST BE SUBSET FOR THE LYCPYR SATS
