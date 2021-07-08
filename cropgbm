#!/usr/bin/env python
#  -*- coding: utf-8 -*-
import os
import argparse
import configparser
from pandas import read_csv
from cropgbm import Engine
from cropgbm import Feature
from cropgbm import Structure
from cropgbm import Visualize
from cropgbm import Parameters as Params
from cropgbm import Preprocessed_Geno as Pg
from cropgbm import Preprocessed_Pheno as Pp


parser = argparse.ArgumentParser()
parser.description = 'Crop Genomic Breeding machine (CropGBM) is a program that integrates pre-processing, ' \
                     'population structure analysis, snp screening, phenotyping prediction, visualization. ' \
                     'The currently accepted geno file format for the program is: ped, bed'

parser.add_argument('-c', '--config_path', type=str,
                    help="The path to the config file, which indicates the value of each parameter "
                         "that the program needs to call. When the parameters in the file are repeated with "
                         "the command line parameters, the program will use the command line value as the standard. "
                         "It is easy to assign values to a large number of parameters through parameter files. "
                         "Using the command line to call specific parameters makes it easy to assign a small number of "
                         "parameters based on the parameter file, and observe result.")
parser.add_argument('-o', '--output-folder', type=str, help="Folder path to store output files. default = './'")


# Preprocessed params: genotype data
parser.add_argument('-pg', '--preprocessed-geno', type=str, choices=['all', 'filter'],
                    help="Start the Preprocessed Module to process genotype data.")
pg = parser.add_argument_group('preprocessed_geno',
                               'The preprocessed_geno module is used to process genotype data. '
                               'For different format data, extract genotype data according to sampleid or snp, '
                               'convert the coding mode of genotype data, etc. to provide the required data for '
                               'downstream. The following parameters are only called in the preprocessed_geno module')
pg.add_argument('--fileprefix', type=str, help="The prefix of the genofile.")
pg.add_argument('--fileformat', type=str, choices=['ped', 'bed'],
                help="The format of the genofile. default = 'bed'")
pg.add_argument('--keep-sampleid-path', type=str,
                help="The sampleid file path, a space/tab-delimited text file with family IDs in the first column "
                     "and within-family IDs in the second column, and removes all unlisted samples from the "
                     "current analysis")
pg.add_argument('--remove-sampleid-path', type=str,
                help="The sampleid file path, a space/tab-delimited text file with family IDs in the first column "
                     "and within-family IDs in the second column, and removes all listed samples from the "
                     "current analysis")
pg.add_argument('--extract-snpid-path', type=str,
                help="The snpid file path, a text file with a list of snp IDs (one per line), "
                     "based on the snpid contained in the file, extracts the genetic "
                     "information of the corresponding sample in the genotype file as an output.")
pg.add_argument('--exclude-snpid-path', type=str,
                help="The snpid file path, a text file with a list of snp IDs (one per line), "
                     "based on the snpid contained in the file, exclude the genetic "
                     "information of the corresponding sample in the genotype file as an output.")
pg.add_argument('--plink-path', type=str, help="The sampleid file path. default = plink")
pg.add_argument('--snpmaxmiss', type=float,
                help="Filter out all snp with missing rates exceeding the value to be removed. default = 0.05")
pg.add_argument('--samplemaxmiss', type=float,
                help="Filter out all sample with missing rates exceeding the value to be removed. default = 0.05")
pg.add_argument('--maf-max', type=float,
                help="Filters out all snp with minor allele frequency below the threshold. default = 0.01")
pg.add_argument('--r2-cutoff', type=float,
                help="Equivalent to the --indep-pairwise parameter in plink. default = 0.8")


# Preprocessed params: phenotype data
parser.add_argument('-pp', '--preprocessed-phe', action="store_true",
                    help="Start the Preprocessed Module to process phenotype data.")
pp = parser.add_argument_group(
    'preprocessed_phe', 'The preprocessed_phe module is used to process phenotype data. '
                        'For different format data, extract phenotype data according to sampleid, '
                        'and normalized phenotypic data to provide the required data for downstream. '
                        'Simultaneously, this module also has the function of visualizing the distribution '
                        'of phenotype data in the form of box plot or density map. '
                        'The following parameters are only called in the preprocessed_phe module')
pp.add_argument('--phe-norm', action="store_true",
                help="Start the Preprocessed Module to normalize phenotypic data.")
pp.add_argument('--phe-recode', type=str, choices=['word2num', 'num2word'],
                help="Recode phenotypic data. Lightgbm only accepts consecutive integers with sample labels "
                     "[0,N) when performing classification tasks. If the training concentration sample comes "
                     "from 5 groups, [0, 1, 2, 3, 4] needs to be used as the reference of 5 groups, "
                     "but this generally does not match the actual group naming. With this parameter, the program can "
                     "achieve reversible conversion between sample tags and [0,N) consecutive integers. "
                     "Provides compliant phenotypic data for downstream classification tasks.")
pp.add_argument('--phe-plot', action="store_true",
                help="Plot the density distribution graph and boxplot for phenotypic data.")
pp.add_argument('--phefile-path', type=str, help="The path to the phenotype file.")
pp.add_argument('--phefile-sep', type=str, help="The separator of phefile. default = ','")
pp.add_argument('--phe-name', type=str, help="Specify the column name of phe to be processed")
pg.add_argument('--ppexsampleid-path', type=str,
                help="The sampleid file path, based on the sampleid contained in the file, extracts the genetic "
                     "information of the corresponding sample in the genotype file as an output.")
pp.add_argument('--ppgroupfile-path', type=str, help="The path to the file containing the groupid information.")
pp.add_argument('--ppgroupfile-sep', type=str, help="The separator of groupfile. default = ','")
pp.add_argument('--ppgroupid-name', type=str, help="Specify the column name of groupid")
pp.add_argument('--num2wordfile-path', type=str,
                help="This parameter is only called when the --phe-recode=num2word. "
                     "The file contains the correspondence between two coding methods of the numerical values "
                     "and category of the phenotype data.")


# Structure params
parser.add_argument('-s', '--structure', action="store_true",
                    help="Start the Structure Module to process genotype data.")
s = parser.add_argument_group(
    'structure', 'The structure module is used to analyze the population structure based on the sample genotype data, '
                 'use the tsne or pca method to reduce the genotype, and then use the optics method to cluster '
                 'and predict the population structure. Simultaneously, the program also supports visualizing '
                 'the predicted population structure and the real population structure in the form of scatter plots. '
                 'The following parameters are only called in the structure module')
s.add_argument('--structure-plot', action="store_true",
               help="Output 3 scatter plots showing the dimensionality reduction and clustering results. "
                    "'_plotredim.pdf' shows the 2-dimensional distribution of the reduced-dimensional data; "
                    "'_plotpredict.pdf' shows the predicted category information of the clustered samples "
                    "based on the reduced-dimensional data. The same category sample points have the same color "
                    "(if The user provides the real category information of the sample, and the actual category "
                    "information of the sample will be displayed in the same way on '_plotredim.pdf'); "
                    "'_plotreach.pdf' shows the reachable distance of the clustering result, "
                    "the same category of points have the same color")
s.add_argument('--genofile-path', type=str, help="The path to the genofile.")
s.add_argument('--sgroupfile-path', type=str, help="The path to the file containing the groupid information.")
s.add_argument('--sgroupid-name', type=str, help="Specify the column name of groupid")
s.add_argument('--sgroupfile-sep', type=str, help="The separator of groupfile. default = ','")
s.add_argument('--redim-mode', type=str, choices=['tsne', 'pca'], help="Dimension reduction method. default = 'pca'")
s.add_argument('--pca-explained-var', type=float,
               help="This parameter is only called when --redim-mode = pca. "
                    "int > 1 or float 0-1. The amount of variance that needs to be explained. default = 0.95")
s.add_argument('--window-size', type=int,
               help="This parameter is only called when --redim-mode = tsne. Window size. default = 20")
s.add_argument('--cluster-mode', type=str, choices=['kmeans', 'optics'], help="Cluster method. default = 'kmeans'")
s.add_argument('--n-clusters', type=int,
               help="This parameter is only called when --cluster-mode = kmeans. "
                    "The number of clusters to form as well as the number of centroids to generate.")
s.add_argument('--optics-min-samples', type=float,
               help="This parameter is only called when --cluster-mode = optics. "
                    "The parameter value is an integer greater than 1 or a decimal between 0-1. "
                    "When the parameter value is an integer, it indicates the minimum number of samples required to "
                    "form the core point; when the parameter value is a decimal, it indicates the ratio of "
                    "the minimum number of samples required to form the core point to the total number of samples. "
                    "default = 0.025")
s.add_argument('--optics-xi', type=float,
               help="This parameter is only called when --cluster-mode = optics. "
                    "The minimum of the reachable distance gradient required for the clustering class. default = 0.05")
s.add_argument('--optics-min-cluster-size', type=float,
               help="This parameter is only called when --cluster-mode = optics. "
                    "The parameter value is an integer greater than 1 or a decimal between 0-1. "
                    "When the parameter value is an integer, it indicates the minimum number of samples required for "
                    "the aggregation class; when the parameter value is a decimal, it indicates the ratio of "
                    "the minimum number of samples required for the aggregation class to the total number of samples. "
                    "default = 0.03")

# Engine params
parser.add_argument('-e', '--engine', action="store_true", help='Start the Engine Module to process genotype data.')
e = parser.add_argument_group(
    'engine', 'The engine module is used to train the inputdata set using the lightGBM model, predict '
              'the phenotype of unknown phenotypic samples, and screen out the snp sites that are decisive for '
              'the phenotype. At the same time, the program supports visualizing the importance of the selected snp '
              'in model in the form of a boxplot and a heatmap. '
              'The following parameters are only called in the engine module')
e.add_argument('-t', '--train', action="store_true", help='Perform train')
e.add_argument('-cv', action="store_true", help='Perform cross validation')
e.add_argument('-p', '--predict', action="store_true", help='Perform predict')
e.add_argument('-sf', '--select-feature', action="store_true",
               help='Perform select features. The result is output as a csv file, boxplot, and heatmap.')
e.add_argument('--bygain-boxplot', action="store_true",
               help='Draw a boxplot with the accuracy of the model added with snp. '
                    'This parameter is only used when -sf.')
e.add_argument('--traingeno', type=str,
               help="The genotype data file path for the training set sample. The file separator is ',', "
                    "the first row is snpid information, and the first column is sampleid information."
                    "This parameter is only used when -t or -cv.")
e.add_argument('--trainphe', type=str,
               help="The phenotype data file path for the training set sample. "
                    "The phenotype file format needs to have a column name and a header. The separator is ','. "
                    "It is recommended to input the output of the original phenotype file "
                    "through the program preprocessing module as input. This parameter is only used when -t or -cv.")
e.add_argument('--validgeno', type=str,
               help="The genotype data file path for the valid set sample. The file separator is ',', "
                    "the first row is snpid information, and the first column is sampleid information."
                    "This parameter is only used when -t.")
e.add_argument('--validphe', type=str,
               help="The phenotype data file path for the valid set sample. "
                    "The phenotype file format needs to have a column name and a header. The separator is ','. "
                    "It is recommended to input the output of the original phenotype file "
                    "through the program preprocessing module as input. This parameter is only used when -t.")
e.add_argument('--testgeno', type=str,
               help="The genotype data file path for the testing set sample. The file separator is ',', "
                    "the first row is snpid information, and the first column is sampleid information."
                    "This parameter is only used when -p.")
e.add_argument('--cv-times', type=int,
               help="Number of cross-validations. This parameter is only used when -sf. default = 5")
e.add_argument('--cv-nfold', type=int,
               help="Number of folds in cross-validation. This parameter is only used when -cv. default = 5")
e.add_argument('--min-detal', type=float,
               help="Minimum percentage value of model accuracy improvement per iteration. "
                    "This parameter is only used when -cv. default = 0.05")
e.add_argument('--min-gain', type=float,
               help="Minimum percentage of gain. The program calculates the total gain of each snp in the model, "
                    "taking the product of the maximum gain and Gainmin as the threshold. This snp will not be "
                    "selected when the total gain of snp in the model is less than the threshold."
                    "This parameter is only used when -sf. default = 0.05")
e.add_argument('--max-colorbar', type=float,
               help="The product of the maximum gain value of each snp in each tree and --max-colorbar "
                    "will be the maximum value of the colorbar in the heatmap. "
                    "This parameter is only used when -sf. default = 0.6")
e.add_argument('--learning-rate', type=float, help="Learning rates for LightGBM. default = 0.1")
e.add_argument('--num-leaves', type=int, help="Max number of leaves in one tree. default = 10")
e.add_argument('--num-threads', type=int,
               help="Number of threads for LightGBM. 0 means default number of threads in OpenMP. default = 0")
e.add_argument('--min-data-in-leaf', type=int,
               help="Minimal number of data in one leaf. Can be used to deal with over-fitting. default = 1")
e.add_argument('--objective', type=str, choices=['regression', 'multiclass'],
               help="The objective of LightGBM. default = 'regression'")
e.add_argument('--device-type', type=str,
               help="Device for the tree learning, you can use GPU to achieve the faster learning. default = 'cpu'")
e.add_argument('--max-depth', type=int,
               help="Limit the max depth for tree model. This is used to deal with over-fitting"
                    " when data is small. Tree still grows leaf-wise. <= 0 means no limit. default = -1")
e.add_argument('--feature-fraction', type=float,
               help="LightGBM will randomly select part of features on each iteration (tree) "
                    "if feature_fraction smaller than 1.0. For example, if you set it to 0.8, "
                    "LightGBM will select 80%% of features before training each tree. default = 1")
e.add_argument('--verbosity', type=int,
               help="Controls the level of LightGBM’s verbosity.  default = 0"
                    "< 0: Fatal, = 0: Error (Warning), = 1: Info, > 1: Debug")
e.add_argument('--num-class', type=int,
               help="This parameter is only called in multi-class "
                    "classification application. The number of categories in classification prediction. default = 1")
e.add_argument('--num-boost-round', type=int, help="Number of boosting iterations. default = 100")
e.add_argument('--early-stopping-rounds', type=int,
               help="This parameter is only called when --validgeno have a value. "
                    "Activates early stopping. The model will train until the validation score stops improving. "
                    "Validation score needs to improve at least every. default = 20")
e.add_argument('--verbose-eval', type=int,
               help="This parameter is only called when --validgeno have a value. "
                    "The eval metric on the valid set is printed at every verbose_eval boosting stage. "
                    "The last boosting stage is also printed. default = 10")
e.add_argument('--init-model-path', type=str,
               help="Filename of LightGBM model used for continue training. Used in the train.")
e.add_argument('--modelfile-path', type=str,
               help="The path of LightGBM model file for feature selection and prediction. Used in the predict.")


def main():
    args = parser.parse_args()
    user_params = vars(args)

    if args.config_path:
        config = configparser.ConfigParser()
        config.read(args.config_path, encoding="utf-8")
    else:
        config = None

    user_params = Params.import_config_params(config, user_params)
    user_params = Params.fill_params_by_default(user_params)

    if args.preprocessed_geno:
        # Preparation parameters
        fpf = Params.check_params(user_params, 'fileprefix')

        savepath = user_params['output_folder']

        try:
            os.mkdir(savepath + 'preprocessed')
        except FileExistsError:
            pass
        savedir = savepath + 'preprocessed/'
        filename = fpf.split('/')[-1]
        spf = savedir + filename

        if args.preprocessed_geno == 'all':

            # analyze
            spf2, spf3 = Pg.analyze_genotype(user_params, fpf, spf)

            # filter
            spf4 = spf + '_filter'
            Pg.exid(user_params, spf3, spf4, spf)

            # INFO summery
            log_list = [x + '.log' for x in [spf, spf2, spf3, spf4]]
            os.system('cat ' + ' '.join(log_list) + ' > ' + spf + '_plink.log')

            try:
                for pf in [spf, spf2, spf3, spf4]:
                    if pf != spf4:
                        for sf in ['.bed', '.bim', '.fam', '.nosex', '.log']:
                            os.remove(pf + sf)
                    else:
                        for sf in ['.nosex', '.log']:
                            os.remove(pf + sf)
            except FileNotFoundError:
                raise FileNotFoundError('plink call failed, please check the log file for the reason ')

        elif args.preprocessed_geno == 'filter':
            spf4 = spf + '_filter'
            Pg.exid(user_params, fpf, spf4, spf)
            os.system('cat ' + spf4 + '.log > ' + spf4 + '_plink.log')

            try:
                os.remove(spf4 + '.nosex')
                os.remove(spf4 + '.log')
            except FileNotFoundError:
                raise FileNotFoundError('plink call failed, please check the log file for the reason ')

        else:
            raise ValueError("The optional value of the 'preprocessed_geno' should be one of ['all', 'filter']")

        Pg.recode012(spf + '_filter')

    if args.preprocessed_phe:

        # Preparation parameters
        phefile_path = Params.check_params(user_params, 'phefile_path')
        phe_name = Params.check_params(user_params, 'phe_name')

        phe_norm = user_params['phe_norm']
        phe_plot = user_params['phe_plot']
        phe_recode = user_params['phe_recode']
        groupfile_path = user_params['ppgroupfile_path']
        phefile_sep = user_params['phefile_sep']
        savepath = user_params['output_folder']

        try:
            os.mkdir(savepath + 'preprocessed')
        except FileExistsError:
            pass
        savedir = savepath + 'preprocessed/'
        filename = phefile_path.strip().split('/')[-1].split('.')[:-1]
        filename = '.'.join(filename)
        spf = savedir + filename

        phe_data = read_csv(phefile_path, header=0, index_col=0, sep=phefile_sep)
        phe_data = phe_data.loc[:, phe_name].dropna(axis=0)

        # extract phe by sampleid
        phe_data = Pp.ex_sample(phe_data, user_params)

        # recode phenotype
        if phe_recode:
            Pp.recodephe(phe_data['phe'], spf, user_params)

        else:
            # extract phe by groupid
            if groupfile_path:
                gphe_data, groupname_list = Pp.ex_gruopid(phe_data, groupfile_path, user_params)

                if phe_plot:
                    Pp.plot_phenodist_scatter(gphe_data, groupname_list, spf, user_params)

                else:
                    for group_name in groupname_list:
                        group_data = gphe_data[gphe_data['group'] == group_name].loc[:, ['sampleid', 'phe']]
                        if phe_norm:
                            Pp.normphe(group_data['phe'], spf + '_group' + str(group_name))
                        else:
                            group_data.to_csv(spf + '_group' + str(group_name) + '.phe', header=True, index=None)

            else:
                if phe_norm:
                    phe_data = Pp.normphe(phe_data['phe'], spf)
                else:
                    phe_data.to_csv(spf + '.phe', header=True, index=False)

                if phe_plot:
                    Pp.plot_phenodist_hist(phe_data, spf, user_params)

    if args.structure:
        filepath = Params.check_params(user_params, 'genofile_path')

        structure_plot = user_params['structure_plot']
        output_folder = user_params['output_folder']

        try:
            os.mkdir(output_folder + 'structure')
        except FileExistsError:
            pass
        savedir = output_folder + 'structure/'
        filename = filepath.strip().split('/')[-1].split('.')[:-1]
        filename = '.'.join(filename)
        spf = savedir + filename

        filedata = read_csv(filepath, header=0, index_col=0)
        cluster, redim_array = Structure.redim_cluster(filedata, spf, user_params)
        filedata_index = filedata.index.values
        if structure_plot:
            Visualize.plot_structure(redim_array, spf, cluster, filedata_index, user_params)

    if args.engine:

        cv = user_params['cv']
        train = user_params['train']
        predict = user_params['predict']
        select_feature = user_params['select_feature']
        output_folder = user_params['output_folder']

        lgb_params_dict = Engine.get_params(user_params)

        try:
            os.mkdir(output_folder + 'engine')
        except FileExistsError:
            pass
        savedir = output_folder + 'engine/'

        if train:
            traingeno_data, trainphe_data = Engine.lgb_train(lgb_params_dict, savedir, user_params)
            if select_feature:
                Feature.exfeature(traingeno_data, trainphe_data, savedir, lgb_params_dict, user_params)

        if cv:
            Engine.lgb_cv(lgb_params_dict, user_params)

        if predict:
            Engine.lgb_predict(user_params, savedir)


if __name__ == '__main__':
    main()

