ml R_packages/3.6.0

CODE=/home/vpeona/2021SatCorvides/Code
WD=/home/vpeona/2021SatCorvides/Intermediate/RE2/PresAbs
mkdir -p $WD && cd $WD

# collect the abundance cluster files from ABUNDANCE_corvides.sh output
cd ../RMSK/
for DIR in $( cat list_dirs )
do
 cp ${DIR}/${DIR}.abundance.clusters.filter.txt $WD/
done

cd $WD

cat *.abundance.clusters.filter.txt | cut -f 1 | sort -u | grep -v 'Cluster' > list_sats

ln -s /home/vpeona/2021SatCorvides/Intermediate/LibSizes/GenomeSize_summary.txt .

Rscript --vanilla $CODE/presence_sat_threshold.R list_sats GenomeSize_summary.txt Figure_presAbs3.3



WD=/home/vpeona/2021SatCorvides/Intermediate/RE2/PresAbs/lycPyr_individuals
mkdir -p $WD && cd $WD

cd ../../RMSK/
for DIR in $( cat list_dirs_lc )
do
 
 cp ${DIR}/${DIR}.abundance.clusters.filter.txt $WD/
done

cd $WD

cat *.abundance.clusters.filter.txt | grep -v 'Cluster' | cut -f1 | sort -u > list_sats

ln -s /home/vpeona/2021SatCorvides/Intermediate/LibSizes/GenomeSize_summary_lc.txt ./GenomeSize_summary_lc.txt

Rscript --vanilla $CODE/presence_sat_threshold_lc.R list_sats GenomeSize_summary_lc.txt Figure_presAbs_lc2.2


CODE=/home/vpeona/2021SatCorvides/Code
WD=/home/vpeona/2021SatCorvides/Intermediate/RE2/PresAbs/Hybrids
mkdir -p $WD && cd $WD

cd ../../RMSK/

ls | grep 'hyb\|corCon_\|corCor_' > list_corvus

for DIR in $( cat list_corvus )
do
 cp ${DIR}/${DIR}.abundance.clusters.filter.txt $WD/
done

cd $WD

cat *.abundance.clusters.filter.txt | grep -v 'Cluster' | cut -f1 | sort -u > list_sats

ln -s /home/vpeona/2021SatCorvides/Intermediate/LibSizes/GenomeSize_summary_corvus.txt ./GenomeSize_summary_corvus.txt

Rscript --vanilla $CODE/presence_sat_threshold_corvus.R list_sats GenomeSize_summary_corvus.txt Figure_presAbs_hybs_1.1
