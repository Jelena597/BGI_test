import pandas as pd
import numpy as np
from matplotlib import pyplot
import umap
from sklearn import cluster
from sklearn.decomposition import PCA


DATA_LOCATION_1 = 'E14.5_E1S3_Dorsal_Midbrain_GEM_CellBin_merge.tsv'
DATA_LOCATION_2 = 'gene_expressions.tsv'

pd.set_option('display.max_rows', 1000)

def _read_data(data_location):
    return pd.read_csv(data_location, sep='\t')


data_1 = _read_data(DATA_LOCATION_1)
data_2 = pd.read_csv(DATA_LOCATION_2, sep='\t', index_col=0)
df = data_1.drop(['x', 'y'], axis=1)
df.sort_values(by=['cell', 'geneID'])


for index, column in df.iterrows():
    data_2.loc[column['cell'], column['geneID']] = column['MIDCounts']

data_2 = data_2.fillna(0)
print(data_2)
data_2.to_csv('gene_expressions_filled.tsv', sep='\t')