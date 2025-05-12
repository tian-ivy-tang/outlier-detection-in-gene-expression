import pandas as pd
import numpy as np

file_name = 'counts/Geuvadis.tsv'

COUNT_INT = np.uint64
if COUNT_INT == np.uint32:
    HIGHEST_COUNT = np.uint32(2 ** 24 - 1)
elif COUNT_INT == np.uint64:
    HIGHEST_COUNT = np.uint64(2 ** 24 - 1)


def tsv_to_df(file_name, dtype=COUNT_INT):
    # read only tsv gene count dataset into data frame 
    file = pd.read_csv(file_name, sep=' ', header=0, )
    print(file.shape)

    file = pd.read_csv(file_name, sep='\t', header=0, )
    print(file.shape)

    file = pd.read_table(file_name, sep=' ', header=0,)
    print(file.shape)

    file = pd.read_table(file_name, sep='\t', header=0,)
    print(file.shape)


tsv_to_df(file_name)