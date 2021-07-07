# -*- coding: utf-8 -*-
import pandas as pd
from pandas import read_csv
from cropgbm import Engine
from cropgbm import Visualize
from cropgbm import Parameters as Params


def extree_info(model_path, num_boost_round, objective, num_class):
    """Extract feature utilization in each tree by creating n dictionaries, which n = num_boost_round

    Parameter:
        model_path: string, Path for storing the lightgbm model.
        num_boost_round: int, Number of tree in lightgbm model.

    Return: dict
        tree_info
    """
    tree_info_dict = {}
    if objective == 'regression':
        for tree_index in range(0, num_boost_round):
            tree_info_dict['tree_' + str(tree_index)] = {}
    elif objective == 'multiclass':
        for tree_index in range(0, num_boost_round*num_class):
            tree_info_dict['tree_' + str(tree_index)] = {}
    for i, row in enumerate(open(model_path)):
        if row.find('feature_names') != -1:
            features_name_list = row.strip().split('=')[-1].split(' ')
            continue
        if row.find('Tree=') != -1:
            tree_index = row.strip().split('=')[-1]
            continue
        if row.find('split_feature') != -1:
            features_index_list = row.strip().split('=')[-1].split(' ')
            features_index_list = [int(i) for i in features_index_list]
            continue
        if row.find('split_gain') != -1:
            features_gain_list = row.strip().split('=')[-1].split(' ')
            features_gain_list = [float(i) for i in features_gain_list]
            seq_index = 0
            for index in features_index_list:
                feature_name = features_name_list[index]
                feature_gain = features_gain_list[seq_index]
                try:
                    tree_info_dict['tree_' + tree_index][feature_name].append(feature_gain)
                except KeyError:
                    tree_info_dict['tree_' + tree_index][feature_name] = [feature_gain]
                seq_index += 1

    return tree_info_dict


def exfeature_by_classification(tree_info_dict, num_boost_round, num_class, save_path=None):
    """Extract information about each tree in the lightgbm classification model.
    And summarize the contribution of each feature to the model.
    The more important the feature, the larger the Gain value corresponding to its feature.

    For multi-category tasks, lightgbm's approach is to adopt a one-to-many strategy,
    which refers to one category as a positive class and the remaining categories as a negative class.
    There are K categories and then K classifiers will be generated.
    Suppose there are K categories, and K categories begin to fit the second tree after fitting the first tree.
    It is not allowed to learn the M trees of a certain category first, and then learn another category.
    M refers to the number of iterations of each classifier set in the multi-category task.
    After the training, a total of M*K trees will be generated.

    Parameter:
        tree_info_dict: dict, Stores information about each tree model in the lightgbm model
        num_boost_round: int, int, Number of tree in lightgbm model.
        num_class: int, Number of categories in a multi-category question
        save_path: string, default=None
            The storage path of feature information extracted from the lightgbm classification model.
            A total of n files will be generated, which n = cat_num.

    Return:
        None
    """
    for icat in range(0, num_class):
        icat_alltree = pd.DataFrame()
        for icat_ktree in range(icat, num_boost_round*num_class, num_class):
            feature_id, feature_gain = [], []
            ktree_info_dict = tree_info_dict['tree_' + str(icat_ktree)]
            for feature_name in ktree_info_dict:
                feature_id.append(feature_name)
                feature_gain.append(sum(ktree_info_dict[feature_name]))
            icat_ktree_df = pd.DataFrame({'featureid': feature_id, 'tree_' + str(icat_ktree): feature_gain})
            try:
                icat_alltree = icat_alltree.merge(icat_ktree_df, how='outer', on='featureid')
            except KeyError:
                icat_alltree = icat_ktree_df

        icat_alltree.fillna(0, inplace=True)
        feature_gain_sum = icat_alltree.sum(axis=1)
        icat_alltree.insert(1, 'featureGain_sum', feature_gain_sum)
        icat_alltree = icat_alltree.sort_values(by='featureGain_sum', axis=0, ascending=False)
        icat_alltree.to_csv(save_path + '.cat_' + str(icat), index=None)


def exfeature_by_regression(tree_info_dict, num_boost_round, save_path=None):
    """Extract information about each tree in the lightgbm regression model.
    And summarize the contribution of each feature to the model.
    The more important the feature, the larger the Gain value corresponding to its feature.

    Parameter:
        tree_info_dict:
        num_boost_round:
        save_path: string, default=None
            The storage path of feature information extracted from the lightgbm classification model.

    Return:
        None
    """
    alltree = pd.DataFrame()
    for itree in range(0, num_boost_round):
        feature_id, feature_gain = [], []
        itree_info_dict = tree_info_dict['tree_' + str(itree)]
        for feature_name in itree_info_dict:
            feature_id.append(feature_name)
            feature_gain.append(sum(itree_info_dict[feature_name]))
        itree_df = pd.DataFrame({'featureid': feature_id, 'tree_' + str(itree): feature_gain})
        try:
            alltree = alltree.merge(itree_df, how='outer', on='featureid')
        except KeyError:
            alltree = itree_df

    alltree.fillna(0, inplace=True)
    feature_gain_sum = alltree.sum(axis=1)
    alltree.insert(1, 'featureGain_sum', feature_gain_sum)
    alltree = alltree.sort_values(by='featureGain_sum', axis=0, ascending=False)
    alltree.to_csv(save_path, index=None)


def exfeature(traingeno_data, trainphe_data, savedir, params_dict, user_params):
    """Perform the select feature with given parameters.

    Parameters:
        traingeno_data: pandas Dataframe,
            The genotype data for the training set sample.
        trainphe_data: pandas Dataframe,
            The phenotype data for the training set sample.
        savedir: str,
        params_dict: dict,
            Parameters for Booster.
        user_params: dict,
            params_name as key, params_value as value

    Return: None
    """

    traingeno = Params.check_params(user_params, 'traingeno')

    bygain_boxplot = user_params['bygain_boxplot']

    cv_times = user_params['cv_times']
    num_boost_round = user_params['num_boost_round']
    objective = user_params['objective']
    num_class = user_params['num_class']
    gainmin = user_params['min_gain']
    colorbar_max = user_params['max_colorbar']

    trainfile_name = traingeno.strip().split('/')[-1].split('.')[:-1]
    trainfile_name = '.'.join(trainfile_name)
    savepath_prefix = savedir + trainfile_name

    if objective == 'regression':
        tree_info_dict = extree_info(savepath_prefix + '.lgb_model', num_boost_round, objective, num_class)
        exfeature_by_regression(tree_info_dict, num_boost_round, savepath_prefix + '.feature')
    elif objective == 'multiclass':
        tree_info_dict = extree_info(savepath_prefix + '.lgb_model', num_boost_round, objective, num_class)
        exfeature_by_classification(tree_info_dict, num_boost_round, num_class, savepath_prefix + '.feature')
    else:
        raise KeyError("The parameter of fileformat is error. Alternate parameters are ['regression', 'multiclass']")

    print('feature extraction is OK')

    feature_data = read_csv(savepath_prefix + '.feature', header=0, index_col=0)
    gainmax = feature_data.iloc[0, 0]
    feature_data = feature_data[feature_data['featureGain_sum'] >= (gainmax * gainmin)]

    if bygain_boxplot:
        bygain_feature_array = feature_data.index.values
        Engine.lgb_iter_feature(bygain_feature_array, traingeno_data, trainphe_data, params_dict,
                                cv_times, num_boost_round, savepath_prefix)

    feature_data.drop('featureGain_sum', axis=1, inplace=True)
    vmax = feature_data.iloc[0, 0] * colorbar_max
    Visualize.plot_heatmap(feature_data, savepath_prefix + '_heatmap.pdf', vmax=vmax)
