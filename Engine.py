# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import lightgbm as lgb
from pandas import read_csv
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from cropgbm import Parameters as Params


def get_params(user_params):
    """Get the values of the parameters necessary for lightgbm training
    from a configuration file or command line and return them as a dictionary.

    Parameters:
        user_params: dict
            params_name as key, params_value as value

    Return: dict
        lgb_params_dict
    """

    lgb_params_dict = {
        'learning_rate': user_params['learning_rate'],
        'num_leaves': user_params['num_leaves'],
        'num_threads': user_params['num_threads'],
        'min_data_in_leaf': user_params['min_data_in_leaf'],
        'objective': user_params['objective'],
        'device_type': user_params['device_type'],
        'max_depth': user_params['max_depth'],
        'feature_fraction': user_params['feature_fraction'],
        'verbosity': user_params['verbosity'],
        'num_class': user_params['num_class'],
    }

    return lgb_params_dict


def lgb_train(params_dict, savedir, user_params):
    """Perform the training with given parameters.

    Return: pandas Dataframe
        trainphe_data, train_boost
    """
    # load params
    traingeno = Params.check_params(user_params, 'traingeno')
    trainphe = Params.check_params(user_params, 'trainphe')

    validgeno = user_params['validgeno']
    init_model_path = user_params['init_model_path']

    num_boost_round = user_params['num_boost_round']
    verbose_eval = user_params['verbose_eval']
    early_stopping_rounds = user_params['early_stopping_rounds']

    trainfile_name = traingeno.strip().split('/')[-1].split('.')[:-1]
    trainfile_name = '.'.join(trainfile_name)
    model_savepath = savedir + trainfile_name + '.lgb_model'

    # train set
    traingeno_data = read_csv(traingeno, header=0, index_col=0)
    trainphe_data = read_csv(trainphe, header=0, index_col=0).dropna(axis=0)
    traingeno_data = traingeno_data.loc[trainphe_data.index.values, :]
    trainsnpid = traingeno_data.columns.values
    train_set = lgb.Dataset(traingeno_data, label=trainphe_data)

    # valid set
    if validgeno:
        validphe = Params.check_params(user_params, 'validphe')

        validgeno_data = read_csv(validgeno, header=0, index_col=0)
        validphe_data = read_csv(validphe, header=0, index_col=0).dropna(axis=0)

        validgeno_data = validgeno_data.loc[validphe_data.index.values, :]
        validsnpid = validgeno_data.columns.values
        ar1ar2, ar1_index, ar2_index = np.intersect1d(ar1=validsnpid, ar2=trainsnpid, return_indices=True)

        if len(validsnpid[ar1_index]) != len(trainsnpid):
            raise KeyError('Part of the snpid in the validgeno missing traingeno.')
        else:
            if (validsnpid == trainsnpid).all():
                pass
            else:
                print('Warning: The snpid order in the validgeno does not match the snpid of the traingeno.')
                validgeno_data = validgeno_data.loc[:, trainsnpid]

        valid_set = lgb.Dataset(validgeno_data, label=validphe_data)
        train_boost = lgb.train(params_dict, train_set, num_boost_round, valid_set, init_model=init_model_path,
                                early_stopping_rounds=early_stopping_rounds, verbose_eval=verbose_eval)

    else:
        train_boost = lgb.train(params_dict, train_set, num_boost_round, init_model=init_model_path)

    train_boost.save_model(model_savepath)
    return traingeno_data, trainphe_data


def lgb_cv(params_dict, user_params):
    """Perform the cross-validation with given paramaters.

    Parameters:
        params_dict: dict, Parameters for Booster.
        user_params: dict
            params_name as key, params_value as value

    Return: eval_hist – Evaluation history.
        The dictionary has the following format: {‘metric1-mean’: [values], ‘metric1-stdv’: [values],
        ‘metric2-mean’: [values], ‘metric2-stdv’: [values], …}.
    """
    traingeno = Params.check_params(user_params, 'traingeno')
    trainphe = Params.check_params(user_params, 'trainphe')

    cv_nfold = user_params['cv_nfold']
    verbose_eval = user_params['verbose_eval']
    num_boost_round = user_params['num_boost_round']
    early_stopping_rounds = user_params['early_stopping_rounds']
    min_detal = user_params['min_detal']
    objective = user_params['objective']
    init_model_path = user_params['init_model_path']

    traingeno_data = read_csv(traingeno, header=0, index_col=0)
    trainphe_data = read_csv(trainphe, header=0, index_col=0).dropna(axis=0)
    traingeno_data = traingeno_data.loc[trainphe_data.index.values, :]
    train_set = lgb.Dataset(traingeno_data, label=trainphe_data)
    eval_dict = lgb.cv(params_dict, train_set, num_boost_round, nfold=cv_nfold, verbose_eval=verbose_eval,
                       init_model=init_model_path, early_stopping_rounds=early_stopping_rounds,
                       stratified=False, shuffle=True)

    if objective == 'regression':
        eval_array1 = np.array(eval_dict['l2-mean'])
        eval_array2 = np.append(np.delete(eval_array1, 0, axis=0), [eval_array1[-1]], axis=0)
        detal = eval_array1 - eval_array2
        drate = np.true_divide(detal, eval_array1) * 100
        exdrate = drate[drate >= min_detal]
        best_iter = len(exdrate)
        best_l2 = round(eval_array1[best_iter-1], 3)
        print('When iterating to the %dth round, the improvement of the accuracy of the model '
              'is less than %.3f for each iteration.' % (best_iter, best_l2))


def lgb_predict(user_params, savedir):
    """Perform the predict with given parameters.

    Return: None
    """
    testgeno = Params.check_params(user_params, 'testgeno')

    objective = user_params['objective']
    modelfile_path = user_params['modelfile_path']

    testfile_name = testgeno.strip().split('/')[-1].split('.')[:-1]
    testfile_name = '.'.join(testfile_name)
    predict_savepath = savedir + testfile_name + '.predict'

    test_data = read_csv(testgeno, header=0, index_col=0)
    train_boost = lgb.Booster(model_file=modelfile_path)
    predict = train_boost.predict(test_data)

    if objective == 'regression':
        pass
    elif objective == 'multiclass':
        print('The prediction result of the category is output in the form of [0,N) integer encoding, '
              'and the result can be converted into the original encoding method by using '
              'the --phe-recode parameter of the Preprocessed module.')
        predict = np.argmax(predict, axis=1)

    predict = pd.DataFrame({'sampleid': test_data.index.values, 'predict_phe': predict})
    predict.to_csv(predict_savepath, header=True, index=False)


def lgb_iter_feature(bygain_feature_array, trainfile_data, trainphe_data, params_dict,
                     cv_times, num_boost_round, savepath_prefix):
    """This function adds features to the training model one by one according to
    the sequence of bygain_feature_array, until all features in the sequence are added to the training model.
    And observe the influence of the features continuously added to the model on the accuracy of the model.

    The prediction accuracy of lightgbm model increases as the number of available features
    in the model increases (the feature richness increases).
    When features are continuously added to the model through iteration,
    the improvement of model prediction accuracy is related to the increase of feature richness
    and the ability of features to distinguish between groups.

    The function will return three matrices, which are the cv results of bygain, random.
    Bygain means that the feature sequence is formed in descending order of the gain value in the training model.
    Random means that the features and sequences in the feature sequence are randomly sampled in the total features.
    By comparing bygain and random, observe the difference between the features selected by lightgbm
    and the randomly selected features.
    By comparing bygain, observe the influence of the nth feature in the bygain sequence on the accuracy.

    Parameters:
        bygain_feature_array: array,
        trainfile_data: string, numpy array, pandas DataFrame
            Data source of Dataset to lgb_cv. If string, it represents the path to txt file.
        trainphe_data: list, numpy 1-D array, pandas Series/one-column DataFrame
            Label of the data to lgb_cv.
        params_dict: dict, Parameters for lgb_cv Booster.
        cv_times: int, repeat times of Cross-validation
        num_boost_round: int,
            Number of boosting iterations. default = 100
        savepath_prefix: str,

    """
    print('lgb_iter_feature_cv_all is in progress. This process will take a long time, please wait patiently')

    exfeature_num = bygain_feature_array.shape[0]
    params_dict['verbosity'] = -1

    # bygain
    bygain_mse = []
    for i in range(1, exfeature_num + 1):
        features = bygain_feature_array[:i]
        features_gene = trainfile_data.loc[:, features]
        train_set = lgb.Dataset(features_gene, label=trainphe_data)
        cv_result = lgb.cv(params_dict, train_set, num_boost_round=num_boost_round, stratified=False)
        if params_dict['objective'] == 'regression':
            imse = cv_result['l2-mean'][-1]
        else:
            imse = cv_result['multi_logloss-mean'][-1]
        bygain_mse.append(imse)

    with PdfPages(savepath_prefix + '_bygain.pdf') as pdf:
        plt.figure(figsize=(9, 9))
        plt.scatter(range(1, exfeature_num + 1), bygain_mse)
        plt.xticks(range(1, len(bygain_feature_array) + 1), bygain_feature_array, rotation=90, size=6)
        plt.title('ByGain')
        plt.ylabel('MSE')
        plt.tight_layout()
        pdf.savefig()
        plt.close()

    # random
    random_imse = []
    for i in range(cv_times):
        snpid_array = trainfile_data.columns.values
        random_array = np.random.choice(range(trainfile_data.shape[1]), exfeature_num, replace=False)
        random_feature_array = snpid_array[random_array]
        random_jmse = []
        for j in range(1, exfeature_num + 1):
            features = random_feature_array[:j]
            features_gene = trainfile_data.loc[:, features]
            train_set = lgb.Dataset(features_gene, label=trainphe_data)
            cv_result = lgb.cv(params_dict, train_set, num_boost_round=num_boost_round, stratified=False)
            if params_dict['objective'] == 'regression':
                jmse = cv_result['l2-mean'][-1]
            else:
                jmse = cv_result['multi_logloss-mean'][-1]
            random_jmse.append(jmse)
        random_imse.append(random_jmse)
    random_mse_data = pd.DataFrame(random_imse)

    with PdfPages(savepath_prefix + '_random.pdf') as pdf:
        fig, axes = plt.subplots(figsize=(9, 9))
        random_mse_data.boxplot(ax=axes, rot=90, fontsize=6, grid=False)
        plt.xticks(range(1, random_mse_data.shape[1] + 1), range(1, random_mse_data.shape[1] + 1))
        plt.title('Random')
        plt.xlabel('Snp Number')
        plt.ylabel('MSE')
        plt.tight_layout()
        pdf.savefig()
        plt.close()
