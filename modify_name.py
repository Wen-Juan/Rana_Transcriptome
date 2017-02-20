import os
import sys

args = sys.argv ### gets a list of command args
test2 = args[1] ## take the 2nd one

file1 = open(test2) ## open file names 

read1_all = ""
read2_all =  ""

for line in file1:
        line = line.rstrip("\n") 
        line = line.split("R2_pair.fastq")[0]
        read1 = line + "R1_pair.fastq"
        read2 = line + "R2_pair.fastq"
 
        read1_all = read1_all + "," + read1
        read2_all = read2_all + "," + read2

read1_all = read1_all.lstrip(",")
read2_all = read2_all.lstrip(",")

print(read1_all)
print(read2_all)
