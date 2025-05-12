import argparse
import os
from pprint import pprint

import pandas as pd
import scipy
# import matplotlib.pyplot as plt
import numpy as np

import helpers as h
import optht

LINE = '----------------------------------'


# standardize = h.standardize
def standardize(data, axis=None):
    return h.transform(data, h._standardize, axis=axis, print_=False, mp=False)


def svd(data, U, s, VT, s_dimension):
    S = scipy.linalg.diagsvd(np.pad(s[:s_dimension], (0, len(s) - s_dimension)), *data.shape)
    return U.dot(S.dot(VT))


def main():
    parser = argparse.ArgumentParser(description='Train model.')
    parser.add_argument('file_name', metavar='file_name', type=str, nargs=1, help='file with count data')
    args = parser.parse_args()
    # print(args.corrected)
    file_name = args.file_name[0]
    process(file_name)

# change input from csv to dfz directly
def process(dfz_fast, file_base_name, unique_id):
    # get dataset name, remove '/counts' prefix & '.tsv'
    # file_base_name = os.path.splitext(file_name)[0][7:]  

    data = h.clean_zs(dfz_fast.values)

    U, s, VT = scipy.linalg.svd(data)

    s_dimension = optht.optht(data, sv=s, sigma=None)
 
    data_new = svd(data, U, s, VT, s_dimension)

    stds = h.std(data_new, axis=1)[:, 0].reshape((data_new.shape[0], 1))
    zs_outrider__ = (data - data_new) / stds
    _data2 = standardize(zs_outrider__, axis=1)

    dfp, dfp_adj = h.dfz_to_dfp_adj(pd.DataFrame(_data2, index=dfz_fast.index, columns=dfz_fast.columns))
    h.save_dfp_adj_to_csv(dfp, dfp_adj, file_base_name, unique_id)

    # # export summary statistics: dataset name, target rank & number of columns to csv file
    # summary_stats = pd.DataFrame(columns = ['Dataset', 'Time', 'Q', 'Column_nb'])
    # summary_stats.loc[0] = [file_base_name, 0.0, s_dimension, data.shape[1]]  # time to be overwritten in pipeline.py
    # summary_stats_filename = "results/summary_stats_" + file_base_name + '_Outsingle_' + unique_id + '.csv'
    # h.save_df_to_csv(summary_stats, summary_stats_filename, index = False)

    return s_dimension, data.shape[1]


if __name__ == '__main__':
    main()
