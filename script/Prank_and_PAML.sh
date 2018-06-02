# run prank with codon function #using bash loop scripts
module add SequenceAnalysis/MultipleSequenceAlignment/prank/140603

for i in *.fa;
do
        prank -d=$i -o=$i -f=phylips -F -codon

done

# run PAML with codeml function
