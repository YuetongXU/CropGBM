# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from pandas import read_csv
from sklearn import preprocessing
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from cropgbm import Parameters as Params


def ex_sample(phe_data, user_params):

    exsampleid_path = user_params['ppexsampleid_path']

    id_array = phe_data.index.values
    if exsampleid_path:
        exsampleid_array = read_csv(exsampleid_path, header=None, index_col=None, sep='\t').iloc[:, 0].values
        ar1ar2, ar1_index, ar2_index = np.intersect1d(ar1=id_array, ar2=exsampleid_array, return_indices=True)
        phe_data = pd.DataFrame({'sampleid': id_array[ar1_index], 'phe': phe_data.values[ar1_index]})
        phe_data.index = id_array[ar1_index]
    else:
        phe_data = pd.DataFrame({'sampleid': id_array, 'phe': phe_data.values})
        phe_data.index = id_array

    return phe_data


def ex_gruopid(phe_data, groupfile_path, user_params):

    groupid_name = Params.check_params(user_params, 'ppgroupid_name')
    groupfile_sep = user_params['ppgroupfile_sep']

    groups_data = read_csv(groupfile_path, header=0, index_col=0, sep=groupfile_sep)
    groups_data = groups_data[groupid_name].to_frame()
    groups_data.columns = ['group']

    gphe_data = phe_data.merge(groups_data, left_index=True, right_index=True)
    groupname_list = sorted(list(set(gphe_data['group'].values)))

    return gphe_data, groupname_list


def normphe(phe_data, savepath_prefix):
    """Normalize and standardize phenotypic data.

    Parameters:
        phe_data: pandas Dataframe, phenotype data
        savepath_prefix: str
            Indicate the processed phenotype file storage path.

    Return: dataframe
    """

    transphe_array = preprocessing.scale(phe_data)
    transphe_data = pd.DataFrame({'sampleid': phe_data.index.values, 'phe_norm': transphe_array.flatten()})
    transphe_data.to_csv(savepath_prefix + '.phenorm', header=True, index=False)
    return transphe_data


def plot_phenodist_scatter(gphe_data, groupname_list, spf, user_params):

    phe_norm = user_params['phe_norm']

    pdf = PdfPages(spf + '_scatter.pdf')
    plt.figure(figsize=(15, 15))

    for index, group_name in enumerate(groupname_list):
        group_data = gphe_data[gphe_data['group'] == group_name].loc[:, ['sampleid', 'phe']]
        noise = np.random.uniform(-0.05, 0.05, group_data.shape[0])

        if phe_norm:
            normphe_data = normphe(group_data['phe'], spf + '_group' + str(group_name), user_params)
            plt.scatter(x=noise + index, y=normphe_data['phe_norm'].values)
            plt.boxplot(normphe_data['phe_norm'].values, positions=[index], widths=0.5,
                        medianprops={'color': 'black'})
        else:
            plt.scatter(x=noise + index, y=group_data['phe'].values)
            plt.boxplot(group_data['phe'].values, positions=[index], widths=0.5,
                        medianprops={'color': 'black'})
            group_data.to_csv(spf + '_group' + str(group_name) + '.phe', header=True, index=None)

    plt.xticks(np.array(range(0, len(groupname_list))), groupname_list, rotation=90)
    plt.tick_params(labelsize=18)
    plt.tight_layout()
    pdf.savefig()
    plt.close()
    pdf.close()


def plot_phenodist_hist(phe_data, spf, user_params):

    phe_norm = user_params['phe_norm']

    pdf = PdfPages(spf + '_plotdd.pdf')
    plt.figure(figsize=(15, 15))

    if phe_norm:
        normphe_data = normphe(phe_data['phe'], spf, user_params)
        normphe_data['phe_norm'].plot(kind='hist', bins=30)
    else:
        phe_data['phe'].plot(kind='hist', bins=30)
        phe_data.to_csv(spf + '.phe', header=True, index=False)

    plt.tick_params(labelsize=20)
    plt.ylabel('')
    plt.tight_layout()
    pdf.savefig()
    plt.close()
    pdf.close()


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
        word2num.to_csv(savepath_prefix + '.word2num', header=True, index=False)
        word_num_dict = dict(zip(phe_word, phe_num))
        num_array = [word_num_dict[i] for i in phe_array]
        num_data = pd.DataFrame({'sampleid': phe_data.index.values, 'phe': num_array})
        num_data.to_csv(savepath_prefix + '.numphe', header=True, index=False)

    elif phe_recode == 'num2word':
        num2wordfile_path = Params.check_params(user_params, 'num2wordfile_path')
        num2wordfile_data = read_csv(num2wordfile_path, header=0, index_col=False)
        phe_word = num2wordfile_data['word'].values
        phe_num = num2wordfile_data['num'].values
        num_word_dict = dict(zip(phe_num, phe_word))
        word_array = [num_word_dict[i] for i in phe_array]
        word_data = pd.DataFrame({'sampleid': phe_data.index.values, 'phe': word_array})
        word_data.to_csv(savepath_prefix + '.wordphe', header=True, index=False)
    else:
        raise KeyError("The parameter of phe_recode is error. "
                       "Alternate parameters are ['word2num', 'num2word']")
