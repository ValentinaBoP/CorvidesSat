ml R_packages/3.6.0

CODE=/home/vpeona/2021SatCorvides/Code
WD=/home/vpeona/2021SatCorvides/Intermediate/RE2/RMSK
cd $WD

ls -d */ | sed 's%/%%g' > list_dirs

for DIR in $( cat list_dirs )
do
 
 cd $DIR

 # reshape the out file
 cat ${DIR}.k2p.noCpG | tr -s ' ' | sed 's/^ //' | tr -s ' ' '\t' | sed 's/kimura=//' | awk '{print $0, $7-$6+1}' OFS="\t" | grep -v 'sat6_b_parRub' > ${DIR}.k2p.noCpG.size
 Rscript --vanilla $CODE/abundance_sat.R ${DIR}.k2p.noCpG.size

 cd $WD
 
done

# abundance for lycPyr individuals
CODE=/home/vpeona/2021SatCorvides/Code
WD=/home/vpeona/2021SatCorvides/Intermediate/RE2/RMSK
cd $WD

ls -d */ | sed 's%/%%g' | grep 'lycPyr1' > list_dirs_lc

for DIR in $( cat list_dirs_lc )
do
 
 cd $DIR

 # reshape the out file
 cat ${DIR}.k2p.noCpG | tr -s ' ' | sed 's/^ //' | tr -s ' ' '\t' | sed 's/kimura=//' | awk '{print $0, $7-$6+1}' OFS="\t" > ${DIR}.k2p.noCpG.size
 Rscript --vanilla $CODE/abundance_sat_lc.R ${DIR}.k2p.noCpG.size

 cd $WD
 
done


# ABUNDANCE SHORT READ ASSEMBLIES
WD=/home/vpeona/2021SatCorvides/Intermediate/Abundance
mkdir -p $WD && cd $WD
CODE=/home/vpeona/2021SatCorvides/Code
for FILE in $( ls *k2p.noCpG )
do
 
# reshape the out file
 cat $FILE | tr -s ' ' | sed 's/^ //' | tr -s ' ' '\t' | sed 's/kimura=//' | awk '{print $0, $7-$6+1}' OFS="\t" | grep -v 'sat6_b_parRub' > $FILE.size
 Rscript --vanilla $CODE/abundance_sat.R $FILE.size
 
done

Rscript --vanilla $CODE/abundance_proportions.R list_sats 
