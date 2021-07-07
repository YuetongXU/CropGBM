# -*- coding: utf-8 -*-
import seaborn
import numpy as np
from pandas import read_csv
from scipy import optimize
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from cropgbm import Parameters as Param


colorlist = ['#B22222', '#F08080', '#FF0000', '#006400', '#3CB371', '#2E8B57', '#00FF7F', '#00FF00', '#7FFF00',
             '#FFFF00', '#EEB422', '#FF6666', '#FF9900', '#996633', '#836FFF', '#4876FF', '#0000FF', '#00008B',
             '#1C86EE', '#63B8FF', '#00BFFF', '#87CEFF', '#00688B', '#996666', '#FF00FF', '#8B008B', '#BF3EFF',
             '#912CEE', '#FF83FA', '#551A8B']
colorarray = np.array(colorlist)


def labcolor_dict(labels_set):
    """Randomly generate corresponding colors according to the label.
    Since there are a total of 30 colors in the built-in color library, if the number of types of labels exceeds 30,
    there may be cases where the colors of the various labels are the same.

    Parameter:
        labels_set: array or list,
            Contains the type and number of label

    Return: dict
        labelscolor_dict
    """

    labels_num = len(labels_set)

    if labels_num > 30:
        # Put back sampling
        color_index = np.random.randint(0, 30, size=labels_num)
    else:
        # No return sampling
        color_index = np.random.choice(range(0, 30), labels_num, replace=False)

    labelscolor_array = colorarray[color_index]
    labelscolor_dict = dict(zip(labels_set, labelscolor_array))

    return labelscolor_dict


def plot_structure(redim_array, savepath_prefix, cluster, index, user_params):
    """Output 3 scatter plots showing the dimensionality reduction and clustering results. '_plotredim.pdf' shows
    the 2-dimensional distribution of the reduced-dimensional data; '_plotpredict.pdf' shows the predicted category
    information of the clustered samples based on the reduced-dimensional data. The same category sample points
    have the same color (if The user provides the real category information of the sample, and the actual category
    information of the sample will be displayed in the same way on '_plotredim.pdf'); '_plotreach.pdf' shows
    the reachable distance of the clustering result, the same category of points have the same color.

    Parameters:
        redim_array: array,
            array after dimension reduction to be drawn.
        savepath_prefix: str,
            Indicate the storage path of the scatter plot drawn by the pca_data.
        cluster: instance,
            Contains the label of each sample calculated by the plot_array data under the
            OPTICS/Keams clustering method.
        user_params: dict
            params_name as key, params_value as value

   Return: None
    """
    print('The program only supports scatter plots for drawing 2D data. '
          'If the data dimension exceeds 2 dimensions, only plot the first 2 dimensions of data is drawn')

    redim_mode = user_params['redim_mode']
    cluster_mode = user_params['cluster_mode']

    with PdfPages(savepath_prefix + '_redim.pdf') as pdf:

        groupfile_path = user_params['sgroupfile_path']

        if groupfile_path:

            groupid_name = Param.check_params(user_params, 'sgroupid_name')

            groupfile_sep = user_params['sgroupfile_sep']

            groupfile_data = read_csv(groupfile_path, header=0, index_col=0, sep=groupfile_sep)
            groupfile_data = groupfile_data.loc[index, :]
            groupfile_array = groupfile_data.loc[:, groupid_name].values
            groupfile_array = groupfile_array.flatten()
            truegroup_set = sorted(list(set(groupfile_array)))
            labelscolor_dict = labcolor_dict(truegroup_set)

            plt.figure(figsize=(12, 9))
            for label in truegroup_set:
                label_index = np.where(groupfile_array == label)[0]
                label_array = redim_array[label_index, :]
                plt.scatter(label_array[:, 0], label_array[:, 1], c=labelscolor_dict[label], label=label)
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2)

        else:
            plt.figure(figsize=(9, 9))
            plt.scatter(redim_array[:, 0], redim_array[:, 1])
        plt.title('Dimensionality reduction')

        if redim_mode == 'pca':
            plt.xlabel('PC1')
            plt.ylabel('PC2')

        elif redim_mode == 'tsne':
            plt.xlabel('X')
            plt.ylabel('Y')

        plt.tight_layout()
        pdf.savefig()
        plt.close()

    with PdfPages(savepath_prefix + '_cluster.pdf') as pdf:

        predictgroup_set = set(cluster.labels_)

        # remove samples that are not clustered
        try:
            predictgroup_set.remove(-1)
        except KeyError:
            pass
        labelscolor_dict = labcolor_dict(predictgroup_set)
        plt.figure(figsize=(9, 9))

        for ilabel in labelscolor_dict:
            idata = redim_array[cluster.labels_ == ilabel]
            plt.scatter(idata[:, 0], idata[:, 1], c=labelscolor_dict[ilabel])
        plt.scatter(redim_array[cluster.labels_ == -1, 0], redim_array[cluster.labels_ == -1, 1],
                    c='#CCCCCC', marker='v', alpha=0.6, label='no_cluster')
        plt.title('Cluster')
        plt.legend()

        if redim_mode == 'pca':
            plt.xlabel('PC1')
            plt.ylabel('PC2')
        elif redim_mode == 'tsne':
            plt.xlabel('X')
            plt.ylabel('Y')

        plt.tight_layout()
        pdf.savefig()
        plt.close()

    if cluster_mode == 'optics':
        with PdfPages(savepath_prefix + '_reachability.pdf') as pdf:
            plt.figure(figsize=(9, 9))
            space = np.arange(redim_array.shape[0])
            reachability = cluster.reachability_[cluster.ordering_]
            labels = cluster.labels_[cluster.ordering_]
            for ilabel in labelscolor_dict:
                ix = space[labels == ilabel]
                ir = reachability[labels == ilabel]
                plt.plot(ix, ir, labelscolor_dict[ilabel], alpha=1)
            plt.plot(space[labels == -1], reachability[labels == -1], 'k.', alpha=0.5)
            plt.plot(space, np.full_like(space, 2., dtype=float), 'k-', alpha=0.5)
            plt.plot(space, np.full_like(space, 0.5, dtype=float), 'k-.', alpha=0.5)
            plt.ylabel('Reachability (epsilon distance)')
            plt.title('Reachability')
            plt.tight_layout()
            pdf.savefig()
            plt.close()


def plot_heatmap(heatmap_data, save_path, vmax=None, vmin=None, cmp='YlOrRd'):
    """Draw a heatmap based on the gain of the feature in each tree of the lightgbm model.
    Heatmap consists of N × M squares, where N = the number of trees in lightgbm, M = number of features in lightgbm.
    The shade of the square color represents the value of the gain of a particular feature in a particular tree.
    Heatmap shows the importance of different features in different trees of lightgbm.

    Parameters:
        heatmap_data: pandas Dataframe,
            Drawing data for heatmap.
        save_path: string,
            Indicate the save path of the heatmap.
        vmax: floats, optional, default=None
            Values to anchor the colormap, otherwise they are inferred from the data and other keyword arguments.
        vmin: same as vmax
        cmp: matplotlib colormap name or object, or list of colors, optional, default='YlOrRd'
            The mapping from data values to color space.

    Return: None
    """
    # heatmap_data.index = heatmap_data['feature_index'].values
    heatmap_data = heatmap_data.T
    pdf = PdfPages(save_path)
    plt.figure(figsize=(30, 15))
    seaborn.heatmap(heatmap_data, vmax=vmax, vmin=vmin, cmap=cmp)
    plt.xticks(size=10)
    plt.yticks(size=10)
    plt.tight_layout()
    pdf.savefig()
    plt.close()
    pdf.close()


def plot_scatter_prelab(prevalue_array, truevalue_array, save_path):
    """Draw a scatter plot with the predicted value as the x-axis and the true value as the y-axis,
    and observe the correlation between the predicted value and the true value by fitting a line.

    Parameters:
        prevalue_array: array,
            Lightgbm model prediction results for samples.
        truevalue_array: array,
            The true value of the sample.
        save_path: string,
            Indicate the save path of the boxplot.

    Return: None
    """
    def func(x, a, b):
        return a * x + b

    pdf = PdfPages(save_path)
    plt.figure(figsize=(15, 15))
    plt.scatter(prevalue_array, truevalue_array, c='k')
    a1, b1 = optimize.curve_fit(func, prevalue_array, truevalue_array)[0]
    x1 = np.array([np.min(prevalue_array), np.max(prevalue_array)])
    y1 = a1 * x1 + b1
    plt.plot(x1, y1, 'r-', label='fit: a=%5.3f, b=%5.3f' % (a1, b1))
    plt.xlabel('predict_value')
    plt.ylabel('true_value')
    plt.legend()

    pdf.savefig(save_path)
    plt.close()
    pdf.close()


def plot_hist(file_data, save_path, title=None):
    pdf = PdfPages(save_path)
    plt.figure(figsize=(15, 15))
    file_data.plot(kind='hist', xlim=(0, 1), bins=30)
    if title is not None:
        plt.title(title)
    plt.legend().remove()
    pdf.savefig()
    plt.close()
    pdf.close()
