import pandas as pd
import numpy as np


DATA_LOCATION = 'E14.5_E1S3_Dorsal_Midbrain_GEM_CellBin_merge.tsv'

pd.set_option('display.max_rows', 1000)

def _read_data(data_location):
    return pd.read_csv(data_location, sep='\t')

data = _read_data(DATA_LOCATION)
df = data.drop(['x', 'y'], axis=1)
df.sort_values(by=['cell', 'geneID'])
genes = [x for x in df['geneID']]
cells = [x for x in df['cell']]
cols = sorted(set(genes), key=genes.index)
rows = sorted(set(cells), key=cells.index)
new_d = pd.DataFrame(columns=list(cols), index=list(rows))
new_d.to_csv('gene_expressions.tsv', sep='\t', header=cols)