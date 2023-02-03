import pandas as pd
import numpy as np


DATA_LOCATION = 'gene_expressions_filled.tsv'
XY = 'x_y_arithmetic.tsv'

def _read_data(data_location):
    return pd.read_csv(data_location, sep='\t')


d = _read_data(DATA_LOCATION)
xy = _read_data(XY)
d = pd.merge(d, xy, how='left', left_on='Unnamed: 0', right_on='cell')
d = d.drop(['cell', 'Unnamed: 0_y'], axis=1)
d.to_csv('x_y_gene_expression.tsv', sep='\t')
