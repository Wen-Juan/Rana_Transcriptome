module add UHTS/Assembler/trinityrnaseq/2.4.0
module add UHTS/Analysis/kallisto/0.43.0
module add R/3.3.2

/software/UHTS/Assembler/trinityrnaseq/2.1.1/util/abundance_estimates_to_matrix.pl \
--est_method kallisto  --out_prefix Amm_all_ \
--name_sample_by_basedir \
Am2_231_L1/abundance.tsv \
Am2_233_L1/abundance.tsv \
Am4_231_L3/abundance.tsv \
Am4_232_L4/abundance.tsv \
Am5_231/abundance.tsv \
Am5_232_L6/abundance.tsv \
Am2_274_L1/abundance.tsv \
Am2_275_L1/abundance.tsv \
Am4_271_L2/abundance.tsv \
Am4_272_L3/abundance.tsv \
Am5_273_L4/abundance.tsv \
Am5_274_L5/abundance.tsv \
Am2_314_all/abundance.tsv \
Am2_315_L6/abundance.tsv \
Am2_316_L3/abundance.tsv \
Am2_318_L2/abundance.tsv \
Am4_314_L7/abundance.tsv \
Am4_317_L7/abundance.tsv \
Am5_311_L6/abundance.tsv \
Am5_312_L8/abundance.tsv \
Am5_313_L7/abundance.tsv \
Am1_434_L5/abundance.tsv \
Am2_433_all/abundance.tsv \
Am2_434_L4/abundance.tsv \
Am4_434_L3/abundance.tsv \
Am4_435_L4/abundance.tsv \
Am5_433_L8/abundance.tsv \
Am2_463_L3/abundance.tsv \
Am2_464_L2/abundance.tsv \
Am4_461_L7/abundance.tsv \
Am5_461_L7/abundance.tsv \
Am6_462_L6/abundance.tsv \
Am6_464_L8/abundance.tsv \
A10FO1/abundance.tsv \
A15MT1/abundance.tsv \
A17FO1/abundance.tsv \
A17MT1/abundance.tsv \
A2FO2/abundance.tsv \
A6FO1/abundance.tsv \
A6MT1/abundance.tsv \
A8FO2/abundance.tsv \
A8MT1/abundance.tsv \
A12MT1/abundance.tsv \
A10FB/abundance.tsv \
A10MB/abundance.tsv \
A12FB/abundance.tsv \
A15FB/abundance.tsv \
A15MB/abundance.tsv \
A16FB/abundance.tsv \
A16MB/abundance.tsv \
A17FB/abundance.tsv \
A17MB/abundance.tsv \
A8MB/abundance.tsv \
A12FL1/abundance.tsv \
A12ML1/abundance.tsv \
A15ML1/abundance.tsv \
A17FL1/abundance.tsv \
A17ML1/abundance.tsv \
A2FL1/abundance.tsv \
A7FL1/abundance.tsv \
A8FL1/abundance.tsv \
A8ML1/abundance.tsv

