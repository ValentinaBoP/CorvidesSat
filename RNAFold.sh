#!/bin/bash -l 
#SBATCH -J RNAFold_corvides
#SBATCH -o RNAFold_corvides.output
#SBATCH -e RNAFold_corvides.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL                           
#SBATCH -t 24:00:00
#SBATCH -A snic2020-5-680
#SBATCH -p core
#SBATCH -n 2

ml bioinfo-tools RNAfold/2.4.17

WD=/proj/sllstore2017073/private/Valentina/2021SatCorvides/Intermediate/RNAfold
mkdir -p $WD && cd $WD

mkdir Output
mkdir Plots

cd Output
ln -s /proj/sllstore2017073/private/Valentina/2021SatCorvides/Results/corvides_sat1.3.fasta .
RNAfold --noGU --noconv -P dna_mathews2004.par -p -g corvides_sat1.3.fasta

for DP in $( ls *_dp.ps )
do
 SS=`echo $DP | sed 's/_dp.ps/_ss.ps/'`
 OUT=`echo $DP | sed 's/_dp.ps/_rss.ps/'`
 /sw/bioinfo/RNAfold/2.4.17/rackham/share/ViennaRNA/bin/relplot.pl $SS $DP > ../$OUT
done
