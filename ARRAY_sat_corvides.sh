ml bioinfo-tools blast/2.11.0+ R_packages/3.6.0 BEDTools/2.29.2

CODE=/home/vpeona/2021SatCorvides/Code
DATA=/home/vpeona/2021SatCorvides/Data/Genomes
LIB=/home/vpeona/2021SatCorvides/Results/corvides_sat1.3_dimers.fasta

WD=/home/vpeona/2021SatCorvides/Intermediate/RE2/Arrays
mkdir -p $WD && cd $WD

# run blast for all the corvus assemblies
for FASTA in $( ls $DATA/cor*.rename.fasta )
do
 OUT=`basename $FASTA`
 OUT=${OUT%.*.*}.blast.out
 blastn -db $FASTA -query $LIB -outfmt '6 qseqid qstart qend qlen sseqid sstart send slen sstrand evalue length pident' -out $WD/$OUT
done

# run blast for all the BOP assemblies
for FASTA in $( ls $DATA/*.rename.fasta | grep -v 'cor' )
do
 OUT=`basename $FASTA`
 OUT=${OUT%.*.*}.blast.out
 #makeblastdb -in $FASTA -out $FASTA -dbtype nucl -parse_seqids
 blastn -db $FASTA -query $LIB -outfmt '6 qseqid qstart qend qlen sseqid sstart send slen sstrand evalue length pident' -out $WD/$OUT
done

# run blast for the PB assemblies of ptiInt and lycPyr
for FASTA in $( ls $DATA/*_PB.fasta )
do
 OUT=`basename $FASTA`
 OUT=${OUT%.*}.blast.out
 makeblastdb -in $FASTA -out $FASTA -dbtype nucl -parse_seqids
 blastn -db $FASTA -query $LIB -outfmt '6 qseqid qstart qend qlen sseqid sstart send slen sstrand evalue length pident' -out $WD/$OUT
done

# identify the arrays
for BLAST in $( ls *.blast.out )
do
 # filter for perc identity and make BED file
 OUT=${BLAST%.*}.bed
 awk '$12 > 0.8 {print $5, $6, $7, $1, $2, $3}' OFS="\t" $BLAST | awk '{if($2 > $3) {print $1, $3, $2, $4, $5, $6} else {print $0}}' OFS="\t" | bedtools sort -i stdin > $OUT
 MERGE=${BLAST%.*}.merge.bed
 bedtools merge -i $OUT -c 4 -o distinct,count -d 1000 | awk '{print $0, $3-$2+1}' OFS="\t" > $MERGE
done

# select longer arrays with the most number of repeats


####### TEST
perl /proj/sllstore2017073/private/scripts/extractFromFasta.pl $DATA/ptiMag.rename.fasta list ptiMag.list > ptiMag.list.fasta
