# make landscape plots for Corvus genre
ml R_packages/3.6.0

CODE=/home/vpeona/2021SatCorvides/Code
WD=/home/vpeona/2021SatCorvides/Intermediate/RE2/Landscape
mkdir -p $WD && cd $WD

# collect out RMSK files
cd ../RMSK/

for DIR in $( cat list_dirs )
do
 cp $DIR/$DIR.k2p.noCpG.size $WD/
done

cd $WD

cp ../PresAbs/list_sats .
ln -s ../../LibSizes/GenomeSize_summary.txt .

ls *k2p.noCpG.size > list_files

for FILE in $( cat list_files )
do
 SPECIES=${FILE%.*.*.*}
 Rscript --vanilla $CODE/landscape_plots.R $FILE GenomeSize_summary.txt list_sats ${SPECIES}_landscape1.1
done

mkdir -p Plots
mv *.pdf Plots



# diversity of the landscapes
mkdir Summary
cd Summary

mv ../*_divergence_abundance.txt .

ls cor*_divergence_abundance.txt > list_corvus

ls *_divergence_abundance.txt | grep -v 'cor' > list_bops

for FILE in $( cat list_corvus )
do   
 SPECIES=${FILE%_*_*}
 awk 'NR > 1 && $4 <= 2 {print $0}' $FILE | cut -f1 | sed 's/\$.*//' | sort -u | wc -l >> num
 echo $SPECIES >> species
done

paste species num > summary_corvus

rm species
rm num

for FILE in $( cat list_bops )
do   
 SPECIES=${FILE%_*_*}
 awk 'NR > 1 && $4 <= 2 {print $0}' $FILE | cut -f1 | sed 's/\$.*//' | sort -u | wc -l >> num
 echo $SPECIES >> species
done

paste species num > summary_bops



ls cor*_divergence_abundance_summary.txt > list_corvus

ls *_divergence_abundance_summary.txt | grep -v 'cor' > list_bops
rm species
rm num
for FILE in $( cat list_corvus )
do   
 SPECIES=${FILE%_*_*_*}
 awk '$2 >= 1 {print $1}' $FILE | sort -u | wc -l >> num
 echo $SPECIES >> species
done

paste species num > summary2_corvus

rm species
rm num

for FILE in $( cat list_bops )
do   
 SPECIES=${FILE%_*_*}
 awk '$2 >= 1 {print $1}' $FILE | sort -u | wc -l >> num
 echo $SPECIES >> species
done

paste species num > summary2_bops

