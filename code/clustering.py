import pandas as pd
import numpy as np
from sklearn.cluster import KMeans, SpectralClustering, DBSCAN
import matplotlib.pyplot as plt
import scanpy as sc
from sklearn.preprocessing import normalize
import anndata
import seaborn as sb
from sklearn.neighbors import kneighbors_graph, DistanceMetric
from sklearn.metrics import silhouette_score

def add_neighbors_and_umapCoords(adata, n_neighs, metric_name, min_d):
    ad1 = sc.AnnData(X=adata.X, var=adata.var, obs=adata.obs)

    sc.pp.neighbors(ad1, n_neighbors=n_neighs, use_rep="X", metric=metric_name)
    sc.tl.umap(ad1, min_dist=min_d, n_components=2)

    return ad1, metric_name

def cluster(ad, algorithm, n_cls, metricName, n, d):
    labels = []
    if algorithm.lower() == 'kmeans':
        clustering = KMeans(n_clusters=n_cls).fit(ad.X)
        labels = clustering.labels_
    elif algorithm.lower() == 'spectral':
        ds = ad.obsp["distances"]
        ds = pd.DataFrame.sparse.from_spmatrix(ds)
        ds_s = (ds + ds.T)/2
        clustering = SpectralClustering(n_clusters=n_cls, affinity='precomputed').fit(ds_s)
        labels = clustering.labels_
    y = pd.Series(labels)
    alg = algorithm.lower()
    labeling(y, ad, metricName, alg, n, d)


def labeling(y, ad, metricName, alg, nn, md):
    fig, (ax1) = plt.subplots(1, figsize=(20,4))


    colors1 = []

    for index, value in y.items():
        if value == 0:
            colors1.append('cl1')
        elif value == 1:
            colors1.append('cl2')
        elif value == 2:
            colors1.append('cl3')
        elif value == 3:
            colors1.append('cl4')
        elif value == 4:
            colors1.append('cl5')
        elif value == 5:
            colors1.append('cl6')
        elif value == 6:
            colors1.append('cl7')
        elif value == 7:
            colors1.append('cl8')
        elif value == 8:
            colors1.append('cl9')
        elif value == 9:
            colors1.append('cl10')
        elif value == 10:
            colors1.append('cl11')
        elif value == 11:
            colors1.append('cl12')


    names = pd.Series(ad1.obs_names)
    names_connect = pd.concat((names, y),axis=1)
    names_connect = names_connect.rename({0:'Cells', 1:'Cluster'},axis=1)


    s = silhouette_score(ad.X, y)
    fig.suptitle(metricName + '-' + alg + '-' + str(nn) + '-' + str(md)+'-' + str(s))
    print(s)

    ad.obs['ClassFound'] = colors1
    found = sc.pl.umap(ad, color="ClassFound", ax=ax1, show=False)


adata = sc.read_csv('x_y_gene_expression.tsv', delimiter='\t', first_column_names=True)

nn = [25, 130, 250]
md = [0.1, 0.4, 0.85]
for n in nn:
    for d in md:
        ad1, metric = add_neighbors_and_umapCoords(adata, n, 'canberra', d)
        cluster(ad1, 'kmeans', 11, metric, n, d)