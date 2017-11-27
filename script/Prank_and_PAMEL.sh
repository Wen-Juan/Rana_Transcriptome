#Prank
#for codon alignment, it requires the sequences to be length which can multiply by 3.

module add UHTS/Assembler/TransDecoder/2.0.1
TransDecoder.LongOrfs -t file.fasta > file_output.cds
TransDecoder.Predict -t target_transcripts.fasta #the final cds or .pep file are not uniqe, need to remove duplicates with shorter sequence lengths

Note: TransDecoder version above v4.0 will generate the longest ORFs per transcript, but not the earlier version, e.g. in cluster the version 2.0.1 generate multiple ORFs.

#modify the header of the transdecoder.cds file.
awk '{print $1}' test.cds
#check the duplicate names, e.g. nameip1.p1 needs to be modified before using the below code, as it will remove the ip1 in the main naming as well. my solution.
#manually replace .p1 (in nameip1.p1) as _ttt, then manually modify the name back, this will work for a few duplicate names.

sed 's/.p1//g' test.cds > test1.cds

# note TransDecoder normally generates more than one open reading frames per transcript.
#select the longest ORF per transcript
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  input.fasta  # unwrap fasta sequences
awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' #calculate size
sort -t ' ' -k2,2 -k1,1nr  #sort on name, inverse length
sort -k1,1 -u -s) #sort on name, unique, stable sort
sed 's/    /./'  #restore name
cut -f 1,2 #cut name, sequence
tr "\t" "\n" < linearized.fasta #go back to fasta
tr "\t" "\n" < linearized.fasta | fold -w 60 #pretty fasta

# run prank with codon function #using bash loop scripts
module add SequenceAnalysis/MultipleSequenceAlignment/prank/140603

for i in *.fa;
do
        prank -d=$i -o=$i -f=phylips -F -codon

done

# run PAML with codeml function
