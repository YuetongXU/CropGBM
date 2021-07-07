# -*- coding: utf-8 -*-
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import OPTICS
from sklearn.cluster import KMeans
from cropgbm import Parameters as Param


def calc_ws(geno_data, ws=20):
    """This function divides the data according to the provided sliding window size (parameter: ws),
    and counts the snp genotype in each window. Complete initial dimensional reduction of data.
    Reduce the computational cost of subsequent pca/tsne dimensionality reduction methods.
    The number of function output columns is the n/ws, where n is the number of input file columns.

    Parameters:
        geno_data: pandas Dataframe, Genotype data to be processed.
        ws: int, default=20, Indicate the window size.

    Return: pandas Dataframe
        ws_data
    """
    ws_data = pd.DataFrame()
    for i, start_loc in enumerate(range(0, geno_data.shape[1], ws)):
        end_loc = start_loc + ws
        data_ws = geno_data.iloc[:, start_loc:end_loc]
        data_ws = data_ws.sum(axis=1)
        ws_data.insert(0, i, data_ws)
    return ws_data


def redim_pca(data, explained_var=0.95):
    """This function uses sklearn.decomposition.PCA to complete the dimensionality reduction of data
    through PCA method.

    Parameters:
        data: pandas Dataframe, data to be processed.
        explained_var: int, default=0.95
            If 0<n_components<1, select the number of components such that the amount of variance
            that needs to be explained is greater than the percentage specified by n_components

    Return: array
         redim_array
    """
    pca = PCA(n_components=explained_var)
    redim_array = pca.fit_transform(data)
    return redim_array


def redim_tsne(data, dim=2):
    """This function uses sklearn.manifold.TSNE to complete the dimensionality reduction of data
    through t-SNE method.

    Parameters:
        data: pandas Dataframe, data to be processed.
        dim: int, default=2
            Dimension of the embedded space.

    Return: array
         redim_array
    """
    tsne = TSNE(n_components=dim, learning_rate=100)
    redim_array = tsne.fit_transform(data)
    return redim_array


def cluster_optics(data, min_samples=0.025, xi=.05, min_cluster_size=0.03):
    """This function uses sklearn.cluster.OPTICS to cluster data
    through OPTICS method.

    Parameters:
        data: array, pandas Dataframe, data to be processed.
        min_samples: int > 1 or float between 0 and 1, default=0.025
            The number of samples in a neighborhood for a point to be considered as a core point.
            Also, up and down steep regions can’t have more then min_samples consecutive non-steep points.
            Expressed as an absolute number or a fraction of the number of samples.
        xi: float, between 0 and 1, optional, default=0.05
            Determines the minimum steepness on the reachability plot that constitutes a cluster boundary.
        min_cluster_size: int > 1 or float between 0 and 1, default=0.03
            Minimum number of samples in an OPTICS cluster,
            expressed as an absolute number or a fraction of the number of samples.

    Return: a instance of OPTICS
        optics
    """
    optics = OPTICS(min_samples=min_samples, xi=xi, min_cluster_size=min_cluster_size)
    optics.fit(data)
    return optics


def cluster_kmeans(data, n_clusters):
    kmeans = KMeans(n_clusters=n_clusters)
    kmeans.fit(data)
    return kmeans


def redim_cluster(geno_data, savepath_prefix, user_params):
    """This function is used to reduce and cluster genetic data.

    Parameters:
        geno_data: pandas Dataframe, Genotype data to be processed.
        savepath_prefix: str
            Indicate the processed phenotype file storage path.
        user_params: dict
            params_name as key, params_value as value

    Return: a instance of OPTICS
        optics
    """

    redim_mode = user_params['redim_mode']
    cluster_mode = user_params['cluster_mode']
    window_size = user_params['window_size']
    pca_explained_var = user_params['pca_explained_var']
    optics_min_samples = user_params['optics_min_samples']
    optics_xi = user_params['optics_xi']
    optics_min_cluster_size = user_params['optics_min_cluster_size']

    # dimensionality reduction
    if redim_mode == 'pca':
        redim_array = redim_pca(geno_data, explained_var=pca_explained_var)
    elif redim_mode == 'tsne':
        ws_data = calc_ws(geno_data, ws=window_size)
        redim_array = redim_tsne(ws_data, dim=2)
    else:
        raise ValueError("The list of optional parameters for redim_mode  is ['pca', 'tsne']")

    redim_data = pd.DataFrame(redim_array)
    redim_data.index = geno_data.index.values
    redim_data.to_csv(savepath_prefix + '.redim', header=False, index=True)

    # clustering
    if cluster_mode == 'kmeans':
        n_clusters = int(Param.check_params(user_params, 'n_clusters'))
        cluster = cluster_kmeans(redim_array, n_clusters)
        label_array = cluster.labels_
    elif cluster_mode == 'optics':
        cluster = cluster_optics(redim_array, optics_min_samples, optics_xi, optics_min_cluster_size)
        label_array = cluster.labels_
    else:
        raise ValueError("The list of optional parameters for reduce dimension  is ['kmeans', 'optics']")

    sampleid_array = geno_data.index.values
    label_data = pd.DataFrame({'sampleid': sampleid_array, 'group': label_array})
    label_data.to_csv(savepath_prefix + '.cluster', header=True, index=False)
    return cluster, redim_array
