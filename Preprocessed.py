# -*- coding: utf-8 -*-
import subprocess
import numpy as np
import pandas as pd
from pandas import read_csv
from sklearn import preprocessing


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


def exid(extract_snpid, exclude_snpid, keep_sampleid, remove_sampleid, fpf, spf, savedir, plink_path):
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
            process = subprocess.Popen(plink_path + ' --bfile ' + fpf + ' --out ' + spf + ' --remove ' +
                                       remove_sampleid + ' --recode compound-genotypes 01 --output-missing-genotype 3 >> '
                                       + savedir + '_preprocessed.log', shell=True)
        else:
            process = subprocess.Popen(plink_path + ' --bfile ' + fpf + ' --out ' + spf +
                                       ' --make-bed --recode compound-genotypes 01 --output-missing-genotype 3 >> '
                                       + savedir + '_preprocessed.log', shell=True)
    process.wait()


def normphe(phe_data, savepath_prefix, user_params):
    """Normalize and standardize phenotypic data.

    Parameters:
        phe_data: pandas Dataframe, phenotype data
        savepath_prefix: str
            Indicate the processed phenotype file storage path.
        user_params: dict
            params_name as key, params_value as value

    Return: array
        transphe_array
    """
    if user_params['norm_mode']:
        norm_mode = user_params['norm_mode']
    else:
        norm_mode = 'z-score'

    if norm_mode in ['0-1', 'max-min']:
        transphe_array = preprocessing.MinMaxScaler().fit_transform(phe_data)
    elif norm_mode == 'z-score':
        transphe_array = preprocessing.scale(phe_data)
    else:
        raise KeyError("The parameter of norm_mode is error. "
                       "Alternate parameters are ['0-1', 'max-min', 'z-score']")
    transphe_data = pd.DataFrame({'sampleid': phe_data.index.values, 'phe_norm': transphe_array.flatten()})
    transphe_data.to_csv(savepath_prefix + '.phenorm', header=True, index=None)
    return transphe_data


def recodephe(phe_data, savepath_prefix, user_params):
    """Recode phenotypic data.

    Parameters:
        phe_data: pandas Dataframe, phenotype data
        savepath_prefix: str
            Indicate the processed phenotype file storage path.
        user_params: dict
            params_name as key, params_value as value

    Return: None
    """
    phe_recode = user_params['phe_recode']
    phe_array = phe_data.values
    if phe_recode == 'word2num':
        phe_word = list(set(phe_array))
        phe_num = np.array(range(len(phe_word)))
        word2num = pd.DataFrame({'num': phe_num, 'word': phe_word})
        word2num.to_csv(savepath_prefix + '.word2num', header=True, index=None)
        word_num_dict = dict(zip(phe_word, phe_num))
        num_array = [word_num_dict[i] for i in phe_array]
        num_data = pd.DataFrame({'sampleid': phe_data.index.values, 'phe': num_array})
        num_data.to_csv(savepath_prefix + '.numphe', header=True, index=None)
    elif phe_recode == 'num2word':
        num2wordfile_path = user_params['num2wordfile_path']
        num2wordfile_data = read_csv(num2wordfile_path, header=0, index_col=None)
        phe_word = num2wordfile_data['word'].values
        phe_num = num2wordfile_data['num'].values
        num_word_dict = dict(zip(phe_num, phe_word))
        word_array = [num_word_dict[i] for i in phe_array]
        word_data = pd.DataFrame({'sampleid': phe_data.index.values, 'phe': word_array})
        word_data.to_csv(savepath_prefix + '.wordphe', header=True, index=None)
    else:
        raise KeyError("The parameter of phe_recode is error. "
                       "Alternate parameters are ['word2num', 'num2word']")
