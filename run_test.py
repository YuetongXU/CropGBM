import os

print('Test 1')
os.system('cropgbm -pg all '
          '--fileprefix ./testdata/genofile')

print('Test 2')
os.system('cropgbm -o ./ -pg all '
          '--fileprefix ./testdata/genofile '
          '--fileformat bed '
          '--snpmaxmiss 0.1 '
          '--samplemaxmiss 0.1 '
          '--maf-max 0.05 '
          '--r2-cutoff 0.7')
#
print('Test 3')
os.system('cropgbm -o ./ -pg all '
          '-c ./testdata/configfile.params')

print('Test 4')
os.system('cropgbm -o ./ -pg filter '
          '--fileprefix ./testdata/genofile '
          '--remove-sampleid-path ./testdata/ksampleid_file.txt '
          '--exclude-snpid-path ./testdata/ksnpid_file.txt')

print('Test 5')
os.system("cropgbm -pp --phe-plot "
          "--phefile-path ./testdata/phefile.txt "
          "--phe-name DTT")

print('Test 6')
os.system("cropgbm -pp "
          "--phefile-path ./testdata/phefile.txt "
          "--phe-name DTT "
          "--phefile-sep ',' "
          "--ppexsampleid-path ./testdata/ksampleid_file.txt "
          "--phe-norm")

print('Test 7')
os.system("cropgbm -pp --phe-plot "
          "--phefile-path ./testdata/phefile.txt "
          "--phe-name DTT "
          "--ppgroupfile-path ./testdata/phefile.txt "
          "--ppgroupfile-sep ',' "
          "--ppgroupid-name paternal_line")

print('Test 8')
os.system("cropgbm -pp "
          "--phefile-path ./testdata/phefile.txt "
          "--phe-name DTT "
          "--phe-norm "
          "--ppgroupfile-path ./testdata/phefile.txt "
          "--ppgroupfile-sep ',' "
          "--ppgroupid-name paternal_line")

print('Test 9')
os.system("cropgbm -pp "
          "--phefile-path ./testdata/phefile.txt "
          "--phe-name pop_class "
          "--phe-recode word2num")

print('Test 10')
os.system("cropgbm -pp "
          "--phefile-path ./preprocessed/phefile.numphe "
          "--phe-name phe "
          "--phe-recode num2word "
          "--num2wordfile-path ./preprocessed/phefile.word2num")

print('Test 11')
os.system("cropgbm -s "
          "--genofile-path ./preprocessed/genofile_filter.geno "
          "--structure-plot "
          "--redim-mode pca "
          "--pca-explained-var 0.98 "
          "--cluster-mode kmeans "
          "--n-clusters 30")

print('Test 12')
os.system("cropgbm -s "
          "--genofile-path ./preprocessed/genofile_filter.geno "
          "--structure-plot "
          "--redim-mode tsne "
          "--window-size 5 "
          "--cluster-mode optics "
          "--optics-min-sample 0.025 "
          "--optics-xi 0.01 "
          "--optics-min-cluster-size 0.03")

print('Test 13')
os.system("cropgbm -s "
          "--genofile-path ./preprocessed/genofile_filter.geno "
          "--structure-plot "
          "--redim-mode tsne "
          "--window-size 5 "
          "--cluster-mode optics "
          "--optics-min-sample 0.025 "
          "--optics-xi 0.01 "
          "--optics-min-cluster-size 0.03 "
          "--sgroupfile-path ./testdata/phefile.txt "
          "--sgroupfile-sep ',' "
          "--sgroupid-name paternal_line")

print('Test 14')
os.system("cropgbm -e -cv "
          "--traingeno ./testdata/train.geno "
          "--trainphe ./testdata/train.phe "
          "--cv-nfold 5 "
          "--min-detal 0.5")

print('Test 15')
os.system("cropgbm -e -cv "
          "--traingeno ./testdata/train.geno "
          "--trainphe ./testdata/train_class.phe "
          "--cv-nfold 5 "
          "--min-detal 0.5 "
          "--objective multiclass "
          "--num-class 6 "
          "--num-boost-round 20 "
          "--num-leaves 5")

print('Test 16')
os.system("cropgbm -e -t "
          "--traingeno ./testdata/train.geno "
          "--trainphe ./testdata/train.phe "
          "--validgeno ./testdata/valid.geno "
          "--validphe ./testdata/valid.phe")

print('Test 17')
os.system("cropgbm -e -t -sf "
          "--bygain-boxplot "
          "--traingeno ./testdata/train.geno "
          "--trainphe ./testdata/train.phe "
          "--min-gain 0.05 "
          "--max-colorbar 0.6 "
          "--cv-times 3")

print('Test 18')
os.system("cropgbm -e -p "
          "--testgeno ./testdata/test.geno "
          "--modelfile-path ./engine/train.lgb_model")

print('Test 19')
os.system("cropgbm -e -t -sf "
          "--bygain-boxplot "
          "--traingeno ./testdata/train.geno "
          "--trainphe ./testdata/train_class.phe "
          "--min-gain 0.05 "
          "--max-colorbar 0.6 "
          "--cv-times 3 "
          "--objective multiclass "
          "--num-class 6 "
          "--num-boost-round 20 "
          "--num-leaves 5")

print('Test 20')
os.system("cropgbm -e -p "
          "--testgeno ./testdata/test.geno "
          "--modelfile-path ./engine/train.lgb_model "
          "--objective multiclass")

print('Test is complete!')
print('Delete temporary files used for test')
os.system('rm -rf preprocessed/ structure/ engine/')


