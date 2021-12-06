# find palindromes with EMBOSS palindrome

module load bioinfo-tools emboss/6.6.0

WD=/proj/sllstore2017073/private/Valentina/2021SatCorvides/Intermediate/Palindrome
mkdir -p $WD && cd $WD

ln -s /proj/sllstore2017073/private/Valentina/2021SatCorvides/Results/corvides_sat1.4.fasta .

palindrome -sequence corvides_sat1.3.fasta -minpallen 10 -maxpallen 100 -gaplimit 20 -nummismatches 1 -outfile corvides_sat1.3.fasta.pal -overlap Y
