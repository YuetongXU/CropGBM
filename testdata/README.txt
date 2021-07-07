These files contain all the data required by the run_test.py

configfile.params is the configuration file of the cropgbm


genofile.map, genofile.ped.zip are the genotype files before preprocessing (ped format).

genofile.bed, genofile.bim, genofile.fam are the genotype files before preprocessing (bed format).


ksampleid_file.txt content is sample ID to be extracted.

ksnpid_file.txt content is SNP ID to be extracted.

rsampleid_file.txt content is sample ID to be removed.

rsnpid_file.txt content is SNP ID to be removed.


phefile.txt is the phenotype file.
The file has 4 columns in total, the first column is the sample ID, the second column is the paternal ID, the third column is the strain, and the fourth column is the phenotype value.


train.geno is the genotype data used for training.

train.phe is a continuous phenotype data corresponding to train.geno.
the phenotype is extracted from the DTT information in the phefile.txt
The file has 2 columns in total, the first column is the sample ID, the second column is the phenotype value.

train_class.phe is a discrete phenotype data corresponding to train.geno.
the phenotype is extracted from the pop_class information in the phefile.txt

valid.geno is the genotype data used for verification.

valid.phe is a continuous phenotype data corresponding to valid.geno.

valid_class.phe is a discrete phenotype data corresponding to valid.geno.

test.geno is the genotype data used for testing.

test.phe is a continuous phenotype data corresponding to valid.geno.

test_class.phe is a discrete phenotype data corresponding to valid.geno.

train.lgb_model is the model file after training, and the tree structure is recorded in the file.


