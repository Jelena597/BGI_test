import pandas as pd


DATA_LOCATION = 'E14.5_E1S3_Dorsal_Midbrain_GEM_CellBin_merge.tsv'

def _read_data(data_location):
    return pd.read_csv(data_location, sep='\t')


data = _read_data(DATA_LOCATION)
data.sort_values(by=['cell'])
i = 0
n = len(data.index)
first = data.iloc[i]['cell']
x = data.iloc[i]['x']
y = data.iloc[i]['y']
i = 1
sumX = 0
sumY = 0
k = 0
xArr = []
yArr = []
cellArr = []
nizK = []
p = 0
for index, column in data.iterrows():
    if p == 0:
        first = column['cell']
        x = column['x']
        y = column['y']
        sumX += x
        sumY += y
        k += 1
        cellArr.append(first)
    else:
        current = column['cell']
        x = column['x']
        y = column['y']
        if current == first:
            sumX += x
            sumY += y
            k += 1
        else:
            arX = sumX / k
            arY = sumY / k
            nizK.append(k)
            k = 0
            xArr.append(arX)
            yArr.append(arY)
            first = current
            cellArr.append(first)
            sumX = x
            sumY = y
            k += 1
    p += 1

arX = sumX / k
arY = sumY / k
nizK.append(k)
k = 0
xArr.append(arX)
yArr.append(arY)


print(len(xArr))
print(len(cellArr))
print(cellArr[-1])
dataXY = {'cell': cellArr, 'x': xArr, 'y': yArr}
df = pd.DataFrame(dataXY)
print(df)
df.to_csv('x_y_arithmetic.tsv', sep='\t')
