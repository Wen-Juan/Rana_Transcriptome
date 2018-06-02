module add UHTS/Analysis/kallisto/0.43.0
for f in A10*_pairedR1.fq.gz
do
kallisto quant -i TvAmKj_trans -o ./output/${f%%_pairedR1.fq.gz} -b 1000 ${f%%_pairedR1.fq.gz}_pairedR1.fq.gz ${f%%_pairedR1.fq.gz}_pairedR2.fq.gz
done
