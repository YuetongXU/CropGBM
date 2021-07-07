#  -*- coding: utf-8 -*-
import codecs


def import_config_params(config, user_params: dict):

    if config is not None:
        all_sections = config.sections()
        for section in all_sections:
            for param in config[section]:
                if not user_params[param]:
                    param_value = config[section][param]
                    if param_value not in ['', 'None', 'false', 'False', 'FALSE']:
                        user_params[param] = config[section][param]

    return user_params


def fill_bool_params(user_params: dict, bool_key: str):

    if user_params[bool_key]:
        pass
    elif not user_params[bool_key]:
        pass
    elif user_params[bool_key] in ['True', 'TRUE', 'true']:
        user_params[bool_key] = True
    else:
        raise ValueError("The optional value of the '" + bool_key + "' should be one of ['True', 'False', None]")


def fill_params_by_default(user_params: dict):
    if not user_params['output_folder']:
        user_params['output_folder'] = './'

    # Preprocessed params: genotype data
    if not user_params['fileformat']:
        user_params['fileformat'] = ' --bfile '
    elif user_params['fileformat'] == 'ped':
        user_params['fileformat'] = ' --file '
    elif user_params['fileformat'] == 'bed':
        user_params['fileformat'] = ' --bfile '
    else:
        raise ValueError("The optional value of the 'fileformat' should be one of ['ped', 'bed', None]")

    if not user_params['plink_path']:
        user_params['plink_path'] = 'plink'

    if not user_params['snpmaxmiss']:
        user_params['snpmaxmiss'] = 0.05
    else:
        user_params['snpmaxmiss'] = float(user_params['snpmaxmiss'])

    if not user_params['samplemaxmiss']:
        user_params['samplemaxmiss'] = 0.05
    else:
        user_params['samplemaxmiss'] = float(user_params['samplemaxmiss'])

    if not user_params['maf_max']:
        user_params['maf_max'] = 0.01
    else:
        user_params['maf_max'] = float(user_params['maf_max'])

    if not user_params['r2_cutoff']:
        user_params['r2_cutoff'] = 0.8
    else:
        user_params['r2_cutoff'] = float(user_params['r2_cutoff'])

    # Preprocessed params: phenotype data
    fill_bool_params(user_params, 'phe_norm')
    fill_bool_params(user_params, 'phe_plot')

    if not user_params['phefile_sep']:
        user_params['phefile_sep'] = ','
    else:
        user_params['phefile_sep'] = user_params['phefile_sep'].strip('\'')
        user_params['phefile_sep'] = codecs.decode(user_params['phefile_sep'], 'unicode_escape')

    if not user_params['ppgroupfile_sep']:
        user_params['ppgroupfile_sep'] = ','
    else:
        user_params['ppgroupfile_sep'] = user_params['ppgroupfile_sep'].strip('\'')
        user_params['ppgroupfile_sep'] = codecs.decode(user_params['ppgroupfile_sep'], 'unicode_escape')

    # Structure params
    fill_bool_params(user_params, 'structure_plot')

    if not user_params['sgroupfile_sep']:
        user_params['sgroupfile_sep'] = ','
    else:
        user_params['sgroupfile_sep'] = user_params['sgroupfile_sep'].strip('\'')
        user_params['sgroupfile_sep'] = codecs.decode(user_params['sgroupfile_sep'], 'unicode_escape')

    if not user_params['redim_mode']:
        user_params['redim_mode'] = 'pca'
    elif user_params['redim_mode'] not in ['tsne', 'pca']:
        raise ValueError("The optional value of the 'redim_mode' should be one of ['tsne', 'pca']")

    if not user_params['pca_explained_var']:
        user_params['pca_explained_var'] = 0.95
    else:
        if float(user_params['pca_explained_var']) < 1:
            user_params['pca_explained_var'] = float(user_params['pca_explained_var'])
        else:
            user_params['pca_explained_var'] = int(user_params['pca_explained_var'])

    if not user_params['window_size']:
        user_params['window_size'] = 20
    else:
        user_params['window_size'] = int(user_params['window_size'])

    if not user_params['cluster_mode']:
        user_params['cluster_mode'] = 'kmeans'
    elif user_params['cluster_mode'] not in ['kmeans', 'optics']:
        raise ValueError("The optional value of the 'cluster_mode' should be one of ['kmeans', 'optics']")

    if not user_params['optics_min_samples']:
        user_params['optics_min_samples'] = 0.025
    else:
        if float(user_params['optics_min_samples']) < 1:
            user_params['optics_min_samples'] = float(user_params['optics_min_samples'])
        else:
            user_params['optics_min_samples'] = int(user_params['optics_min_samples'])

    if not user_params['optics_xi']:
        user_params['optics_xi'] = 0.05
    else:
        if float(user_params['optics_xi']) < 1:
            user_params['optics_xi'] = float(user_params['optics_xi'])
        else:
            user_params['optics_xi'] = int(user_params['optics_xi'])

    if not user_params['optics_min_cluster_size']:
        user_params['optics_min_cluster_size'] = 0.03
    else:
        if float(user_params['optics_min_cluster_size']) < 1:
            user_params['optics_min_cluster_size'] = float(user_params['optics_min_cluster_size'])
        else:
            user_params['optics_min_cluster_size'] = int(user_params['optics_min_cluster_size'])

    # Engine params
    fill_bool_params(user_params, 'bygain_boxplot')

    if not user_params['cv_times']:
        user_params['cv_times'] = 5
    else:
        user_params['cv_times'] = int(user_params['cv_times'])

    if not user_params['cv_nfold']:
        user_params['cv_nfold'] = 5
    else:
        user_params['cv_nfold'] = int(user_params['cv_nfold'])

    if not user_params['min_detal']:
        user_params['min_detal'] = 0.05
    else:
        user_params['min_detal'] = float(user_params['min_detal'])

    if not user_params['min_gain']:
        user_params['min_gain'] = 0.05
    else:
        user_params['min_gain'] = float(user_params['min_gain'])

    if not user_params['max_colorbar']:
        user_params['max_colorbar'] = 0.6
    else:
        user_params['max_colorbar'] = float(user_params['max_colorbar'])

    # Engine params: LGB Hyperparameter
    if not user_params['learning_rate']:
        user_params['learning_rate'] = 0.1
    else:
        user_params['learning_rate'] = float(user_params['learning_rate'])

    if not user_params['num_leaves']:
        user_params['num_leaves'] = 10
    else:
        user_params['num_leaves'] = int(user_params['num_leaves'])

    if not user_params['num_threads']:
        user_params['num_threads'] = 0
    else:
        user_params['num_threads'] = int(user_params['num_threads'])

    if not user_params['min_data_in_leaf']:
        user_params['min_data_in_leaf'] = 1
    else:
        user_params['min_data_in_leaf'] = int(user_params['min_data_in_leaf'])

    if not user_params['objective']:
        user_params['objective'] = 'regression'
    elif user_params['objective'] not in ['regression', 'multiclass']:
        raise ValueError("The optional value of the 'objective' should be one of ['regression', 'multiclass']")

    if not user_params['device_type']:
        user_params['device_type'] = 'cpu'
    elif user_params['device_type'] not in ['cpu', 'gpu']:
        raise ValueError("The optional value of the 'device_type' should be one of ['cpu', 'gpu']")

    if not user_params['max_depth']:
        user_params['max_depth'] = -1
    else:
        user_params['max_depth'] = int(user_params['max_depth'])

    if not user_params['feature_fraction']:
        user_params['feature_fraction'] = 1
    else:
        user_params['feature_fraction'] = float(user_params['feature_fraction'])

    if not user_params['verbosity']:
        user_params['verbosity'] = 0
    else:
        user_params['verbosity'] = int(user_params['verbosity'])

    if not user_params['num_class']:
        user_params['num_class'] = 1
    else:
        user_params['num_class'] = int(user_params['num_class'])

    if not user_params['num_boost_round']:
        user_params['num_boost_round'] = 100
    else:
        user_params['num_boost_round'] = int(user_params['num_boost_round'])

    if not user_params['early_stopping_rounds']:
        user_params['early_stopping_rounds'] = 20
    else:
        user_params['early_stopping_rounds'] = int(user_params['early_stopping_rounds'])

    if not user_params['verbose_eval']:
        user_params['verbose_eval'] = 10
    else:
        user_params['verbose_eval'] = int(user_params['verbose_eval'])

    return user_params


def check_params(user_params: dict, param: str):

    if not user_params[param]:
        raise ValueError("Parameter '" + param + "' is not assigned ")

    return user_params[param]








