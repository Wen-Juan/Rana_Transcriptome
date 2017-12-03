infile1 = open("xtro_id40_hsp50.txt","r")
infile2 = open("Xtropicalisv9.0.Named.primaryTrs.fa","r")
outfile = open("xtro_id40_hsp50.fasta","w")

table = {}
keep = "no"
nothing = 0

for line1 in infile1:
	line1 = line1.strip("\n")
	line1 = line1.split()
	seq_I_want = line1[0]
	table[seq_I_want] = seq_I_want
        print(line1)

for line2 in infile2:
	if line2.startswith(">"):
		line2 = line2.strip("\n")
		fullline2 = line2
		line2 = line2.split(">")
		seq_name = line2[1]
		if table.has_key(seq_name):
			outfile.write(fullline2 + "\n")
			keep = "yes"
		else:
			keep = "no"
	else:
		if keep == "yes":
			outfile.write(line2)
		else:
			nothing += 1
