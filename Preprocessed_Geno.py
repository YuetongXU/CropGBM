# -*- coding: utf-8 -*-
import os
import subprocess
import pandas as pd
from pandas import read_csv
from cropgbm import Visualize


# one-hot code
genolist = ['00', '01', '10', '11']
onehotlist = ['0', '1', '1', '2']
codedict = dict(zip(genolist, onehotlist))


def recode012(fileprefix):
    ped_path = fileprefix + '.ped'
    map_path = fileprefix + '.map'
    geno_path = fileprefix + '.geno'
    geno_file = open(geno_path, 'w+')

    with open(map_path) as map_file:
        snpid_list = []
        for i, row in enumerate(map_file):
            row = row.strip().split('\t')
            snpid_list.append(row[1])
        header = ','.join(snpid_list)
        geno_file.write(header + '\n')

    with open(ped_path) as ped_file:
        for i, row in enumerate(ped_file):
            row = row.strip().split()
            geno = row[6:]
            geno_list = [row[1]]
            for j in range(0, len(geno)):
                diploid = codedict[geno[j]]
                geno_list.append(diploid)
            genos = ','.join(geno_list)
            geno_file.write(genos + '\n')
    geno_file.close()


def exid(user_params, fpf, spf, savedir):

    extract_snpid = user_params['extract_snpid_path']
    exclude_snpid = user_params['exclude_snpid_path']
    keep_sampleid = user_params['keep_sampleid_path']
    remove_sampleid = user_params['remove_sampleid_path']

    plink_path = user_params['plink_path']
    plink_result = os.popen('which ' + plink_path).read().strip()
    if not os.path.exists(plink_result):
        raise IOError("Can't find plink by --plink-path")

    if extract_snpid:
        if keep_sampleid:
            process = subprocess.Popen(plink_path + ' --bfile ' + fpf + ' --out ' + spf + ' --extract ' + extract_snpid +
                                       ' --keep ' + keep_sampleid + ' --recode compound-genotypes 01 '
                                                                    '--output-missing-genotype 3 >> '
                                       + savedir + '_preprocessed.log', shell=True)

        elif remove_sampleid:
            process = subprocess.Popen(plink_path + ' --bfile ' + fpf + ' --out ' + spf + ' --extract ' + extract_snpid +
                                       ' --remove ' + remove_sampleid + ' --recode compound-genotypes 01 '
                                                                        '--output-missing-genotype 3 >> '
                                       + savedir + '_preprocessed.log', shell=True)
        else:
            process = subprocess.Popen(plink_path + ' --bfile ' + fpf + ' --out ' + spf + ' --extract ' + extract_snpid +
                                       ' --recode compound-genotypes 01 --output-missing-genotype 3 >> '
                                       + savedir + '_preprocessed.log', shell=True)
    elif exclude_snpid:
        if keep_sampleid:
            process = subprocess.Popen(plink_path + ' --bfile ' + fpf + ' --out ' + spf + ' --exclude ' + exclude_snpid +
                                       ' --keep ' + keep_sampleid + ' --recode compound-genotypes 01 '
                                                                    '--output-missing-genotype 3 >> '
                                       + savedir + '_preprocessed.log', shell=True)
        elif remove_sampleid:
            process = subprocess.Popen(plink_path + ' --bfile ' + fpf + ' --out ' + spf + ' --exclude ' + exclude_snpid +
                                       ' --remove ' + remove_sampleid + ' --recode compound-genotypes 01 '
                                                                        '--output-missing-genotype 3 >> '
                                       + savedir + '_preprocessed.log', shell=True)
        else:
            process = subprocess.Popen(plink_path + ' --bfile ' + fpf + ' --out ' + spf + ' --exclude ' + exclude_snpid +
                                       ' --recode compound-genotypes 01 --output-missing-genotype 3 >> '
                                       + savedir + '_preprocessed.log', shell=True)
    else:
        if keep_sampleid:
            process = subprocess.Popen(plink_path + ' --bfile ' + fpf + ' --out ' + spf + ' --keep ' +
                                       keep_sampleid + ' --recode compound-genotypes 01 --output-missing-genotype 3 >> '
                                       + savedir + '_preprocessed.log', shell=True)
        elif remove_sampleid:
            process = subprocess.Popen(plink_path + ' --bfile ' + fpf + ' --out ' + spf + ' --remove ' + remove_sampleid
                                       + ' --recode compound-genotypes 01 --output-missing-genotype 3 >> '
                                       + savedir + '_preprocessed.log', shell=True)
        else:
            process = subprocess.Popen(plink_path + ' --bfile ' + fpf + ' --out ' + spf +
                                       ' --recode compound-genotypes 01 --output-missing-genotype 3 >> '
                                       + savedir + '_preprocessed.log', shell=True)
    process.wait()


def analyze_genotype(user_params, fpf, spf):

    fileformat = user_params['fileformat']
    snp_miss = user_params['snpmaxmiss']
    sample_miss = user_params['samplemaxmiss']
    maf = user_params['maf_max']
    r2 = user_params['r2_cutoff']

    plink_path = user_params['plink_path']
    plink_result = os.popen('which ' + plink_path).read().strip()
    if not os.path.exists(plink_result):
        raise IOError("Can't find plink by --plink-path")

    process = subprocess.Popen(plink_path + fileformat + fpf + ' --out ' + spf +
                               ' --make-bed --freqx --missing --geno ' + str(snp_miss) +
                               ' --mind ' + str(sample_miss) + ' --maf ' + str(maf) + ' >> '
                               + spf + '_preprocessed.log', shell=True)
    process.wait()

    freqx_data = read_csv(spf + '.frqx', sep='\t')
    sample_num = read_csv(spf + '.fam', sep=' ').shape[0]

    maf_data = pd.DataFrame({'maf': (freqx_data['C(HOM A1)'] * 2 + freqx_data['C(HET)']) / (sample_num * 2)})
    het_rate_data = pd.DataFrame(freqx_data['C(HET)'] / sample_num)

    Visualize.plot_hist(maf_data, spf + '_maf.pdf', title='MAF')
    Visualize.plot_hist(het_rate_data, spf + '_het.pdf', title='Het Rate')

    with open(spf + '.imiss') as imiss_file:
        header = imiss_file.readline()
        imiss_rate = []
        for i, row in enumerate(imiss_file):
            imiss_rate.append(float(row.strip()[-1]))
    imiss_data = pd.DataFrame(imiss_rate)

    with open(spf + '.lmiss') as lmiss_file:
        header = lmiss_file.readline()
        lmiss_rate = []
        for i, row in enumerate(lmiss_file):
            lmiss_rate.append(float(row.strip()[-1]))
    lmiss_data = pd.DataFrame(lmiss_rate)

    Visualize.plot_hist(imiss_data, spf + '_imiss.pdf', title='Sample Missing Rate')
    Visualize.plot_hist(lmiss_data, spf + '_lmiss.pdf', title='Snp Missing Rate')

    # fill the missing snp
    spf2 = spf + '_f'
    process = subprocess.Popen(plink_path + ' --bfile ' + spf + ' --out ' + spf2 +
                               ' --make-bed --fill-missing-a2 >> '
                               + spf + '_preprocessed.log', shell=True)
    process.wait()

    # indep and recode
    spf3 = spf + '_r'
    process = subprocess.Popen(plink_path + ' --bfile ' + spf2 + ' --out ' + spf +
                               ' --indep-pairwise 50 10 ' + str(r2) +
                               ' >> ' + spf + '_preprocessed.log', shell=True)
    process.wait()

    process = subprocess.Popen(plink_path + ' --bfile ' + spf2 + ' --out ' + spf3 + ' --extract ' + spf +
                               '.prune.in --make-bed >> ' + spf + '_preprocessed.log', shell=True)
    process.wait()
    return spf2, spf3



