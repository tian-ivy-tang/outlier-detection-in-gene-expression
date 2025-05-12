import os
from io import StringIO
import math


import numpy as np
import pandas as pd
import scipy
import mpmath as mp
import multiprocess
#import statsmodels.api as sm  # for p-value FDR adjustment
import statsmodels.stats.multitest as multitest

COUNT_INT = np.uint64
if COUNT_INT == np.uint32:
    # HIGHEST_COUNT = np.uint32(2 ** 32 - 1)
    HIGHEST_COUNT = np.uint32(2 ** 24 - 1)
elif COUNT_INT == np.uint64:
    # HIGHEST_COUNT = np.uint64(2 ** 64 - 1)
    # HIGHEST_COUNT = np.uint64(2 ** 32 - 1)
    HIGHEST_COUNT = np.uint64(2 ** 24 - 1)


def tsv_to_df(file_name, dtype=COUNT_INT):
    # read only tsv gene count dataset into data frame 
    file = pd.read_csv(file_name, sep=' ', header=0, )
    # file = file.set_index('geneID')  # activate if on server dataset
    file = file.apply(pd.to_numeric, errors='coerce')

    # filter out genes 0 counts in all samples
    file = file.loc[(file != 0).any(axis=1)]

    # filter fpkm tbc

    return file


def csv_to_df(file_name, dtype=np.uint64):
    # read intermediate csv outputs into data frame, csv file is comma-delimited by default
    file = pd.read_csv(file_name, sep='\t', header=0, index_col=0)  
    file = file.apply(pd.to_numeric, errors='coerce')

    return file


def get_size_factors(df):
    _data = df.values
    N = _data.shape[1]
    _rows = []
    for _row in _data:
        if np.all(_row != 0):
            _rows.append(_row)
    J = len(_rows)
    data = np.array(_rows)
    data.shape = (J, N)
    counts_normalized = np.zeros(data.shape, dtype=np.float64)
    for j in range(data.shape[0]):
        _row = np.array(data[j, :])
        # counts = np.array([max(1, count) for count in row])
        counts = np.array([count for count in _row if count != 0])
        # print(counts)

        # Geometric mean for one gene (all samples)
        if len(counts) > 0:
            # denominator = math.exp(np.mean(np.log(counts)))
            denominator = mp_gmean(counts)
            # print(counts)
            # denominator = reduce(operator.mul, counts, 1) ** (1 / len(counts))
            # denominator = reduce(lambda x, y: x*y, counts)**(1.0/len(counts))
            # print(denominator)
        else:
            denominator = 0

        if denominator == 0:
            counts_normalized[j] = np.zeros(data.shape[1], dtype=np.float64)
        else:
            counts_normalized[j] = mp_fdiv(_row, denominator)
    # print(counts_normalized)
    size_factors = np.zeros(data.shape[1], dtype=np.float64)
    for i in range(data.shape[1]):
        column = np.array([count for count in counts_normalized[:, i] if count != 0])
        size_factors[i] = np.median(column)
    return size_factors


def save_df_to_csv(data_df, file_name, index = False):
    # if not os.path.exists(file_name):
    data_df.to_csv(file_name, index = False)

    #     with open(file_name, 'r') as f:
    #         # Removing the initial separator, as it makes problems for pd.read_csv
    #         text = f.read()[1:]

    #     with open(file_name, 'w') as f:
    #         f.write(text)
    # else:
    #     print('The file', file_name, 'already exists, not saving...')


def csv_to_df(file_name, dtype=COUNT_INT):
    # read tsv dataset into data frame
    file = pd.read_table(file_name, sep=' ', header=0, )
    file = file.apply(pd.to_numeric, errors='coerce')

    # filter out genes with 0 counts in all samples
    file = file.loc[(file != 0).any(axis=1)]

    return file

# seperate FDR adjusted p value calculation & export
def dfz_to_dfp_adj(dfz):
    # save_df_to_csv(dfz, filename)

    dfp = pd.DataFrame(convert_zscores_to_pvalues(dfz.values), index=dfz.index, columns=dfz.columns)
    
    dfp_adj = pd.DataFrame(0.0, columns=dfp.columns, index=dfp.index, dtype=np.float64)

    # Apply FDR adjustment to each column of p value matrix, method BY
    for n in range(dfp_adj.shape[1]):
        pv = dfp.iloc[:,n].to_list()
        pv_adj = multitest.multipletests(pv, alpha=0.01, method='fdr_by')[1]  # return adjusted p values
        dfp_adj.iloc[:, n] = pv_adj.reshape(dfp_adj.shape[0],1)

    return dfp, dfp_adj


def save_dfp_adj_to_csv(dfp, dfp_adj, file_base_name, unique_id):    

    # get aberrant genes per sample
    df_aberrant_per_sample = aberrant_per_sampe(dfp_adj)
    filename_aberrant_per_sample = 'results/aberrant_genes_Outsingle_' + file_base_name  + '_' + unique_id + '.csv'
    
    save_df_to_csv(df_aberrant_per_sample, filename_aberrant_per_sample)
    
    # get gene results table if Geuvadis dataset
    # if file_base_name == 'GTExSkinSmall':
    #     gene_results = format_gene_results(dfp, dfp_adj)
    #     gene_results_filename = 'results/gene_results_' + file_base_name + '_Outsingle_' + unique_id + '.csv'

    # get gene results table
    gene_results = format_gene_results(dfp, dfp_adj)
    gene_results_filename = 'results/gene_results_Outsingle_' + file_base_name + '_' + unique_id + '.csv'

    save_df_to_csv(gene_results, gene_results_filename)


def aberrant_per_sampe(dfp_adj):
    # get aberrant genes per sample with FDR < 0.05
    aberrant_counts = dfp_adj.apply(lambda x: (x < 0.05).sum())  # by default apply to each column
    df_aberrant_per_sample = pd.DataFrame(aberrant_counts, columns=['aberrant_per_sample'])
    df_aberrant_per_sample.insert(0, 'sampleID', df_aberrant_per_sample.index)  # insert rownames column
    df_aberrant_per_sample = df_aberrant_per_sample.sort_values(by=['aberrant_per_sample'])

    return df_aberrant_per_sample
 

def format_gene_results(dfp, dfp_adj): 
    # results table consists of 5 columns: geneID, sampleID, pValue, padjust, aberrant (TRUE/FALSE)
    # result table is sorted by padjust in ascending order

    # get melted p values table
    dfp.insert(0, "geneID", dfp.index)
    dfp_melted = pd.melt(dfp, id_vars='geneID', var_name='sampleID', value_name='pValue')
    dfp_melted['gene_sample'] = dfp_melted['geneID'].astype(str) + dfp_melted['sampleID']

    # get melted adjusted p values table & add aberrant column (TRUE/FALSE)
    dfp_adj.insert(0, "geneID", dfp_adj.index)
    dfp_adj_melted = pd.melt(dfp_adj, id_vars='geneID', var_name='sampleID', value_name='padjust')
    dfp_adj_melted['gene_sample'] = dfp_adj_melted['geneID'].astype(str) + dfp_adj_melted['sampleID']

    dfp_adj_melted['aberrant'] = dfp_adj_melted['padjust'].apply(lambda x: 'TRUE' if x < 0.05 else 'FALSE')
    dfp_adj_melted = dfp_adj_melted.drop(['geneID', 'sampleID'], axis=1)

    # merge 2 tables by gene_sample
    # add adjust p values column directly to p values table is much faster
    # but merge by gene_sample is the safest
    gene_results = pd.merge(dfp_melted, dfp_adj_melted, on='gene_sample').drop(['gene_sample'], axis=1)

    # sort by adjusted p values & export
    gene_results = gene_results.sort_values('padjust')

    return gene_results


def convert_zscores_to_pvalues(zs__):
    return 2 * scipy.stats.norm.cdf(-np.abs(zs__))


def clean_zs(data):
    _tmp = np.copy(data)
    _tmp[np.isinf(_tmp)] = 0
    _tmp[np.isneginf(_tmp)] = 0
    data[np.isinf(data)] = max(7, np.abs(_tmp).max())
    data[np.isneginf(data)] = min(-7, -np.abs(_tmp).max())
    return data


mp_power = np.frompyfunc(mp.power, 2, 1)
mp_fdiv = np.frompyfunc(mp.fdiv, 2, 1)


def mp_gmean(array):
    return mp_power(mp_fprod((mp.mpf(str(e)) for e in array)), (1.0 / len(array)))


def mp_fprod2(a, b):
    return mp.fprod([a, b])


def mp_fprod(list_):
    f = np.frompyfunc(mp_fprod2, 2, 1)
    res = mp.mpf('1')
    for e in list_:
        res = f(res, e)
    return res


def transform(a, transform_f, axis=None, print_=False, mp=True):
    # a = a.astype(np.float64)
    # a += 0.001
    a_new = np.zeros_like(a, dtype=np.float64)
    if len(a.shape) == 1 or axis is None:
        shape = a.shape
        a = np.ravel(a)
        a_new = transform_f(a)
        a_new.shape = shape
        return a_new
    elif axis == 0:  # Column-wise'
        def f(c, column):
            if print_:
                print('Processing column', c)
            return c, transform_f(column)
        # No need??? when called in optht_svd_zs.py, mp=False
        if mp:
            # multiprocessing
            # n_parts = math.floor(multiprocess.cpu_count() * 3 / 4) or 1
            n_parts = math.floor(multiprocess.cpu_count() * 5 / 12) or 1
            with multiprocess.Pool(processes=n_parts) as pool:
                results = pool.starmap(f, zip(range(a.shape[1]), (a[:, c] for c in range(a.shape[1]))))
            results.sort()
            for c in range(a.shape[1]):
                a_new[:, c] = results[c][1]
        else:
            # Single-processing
            for c in range(a.shape[1]):
                if print_:
                    print('Processing column', c)
                a_new[:, c] = transform_f(a[:, c])
        return a_new
    elif axis == 1:  # Row-wise
        def f(r, row):
            if print_:
                print('Processing row', r)
            return r, transform_f(row)
        if mp:
            # multiprocessing
            # n_parts = math.floor(multiprocess.cpu_count() * 3 / 4) or 1
            n_parts = math.floor(multiprocess.cpu_count() * 5 / 12) or 1
            with multiprocess.Pool(processes=n_parts) as pool:
                results = pool.starmap(f, zip(range(a.shape[0]), (a[r, :] for r in range(a.shape[0]))))
            results.sort()
            for r in range(a.shape[0]):
                a_new[r, :] = results[r][1]
        else:
            # Single-processing
            for r in range(a.shape[0]):
                if print_:
                    print('Processing row', r)
                a_new[r, :] = transform_f(a[r, :])
        return a_new
    else:
        raise Exception('Cannot happen')


def _std(data):
    assert len(data.shape) == 1
    N = len(data)
    data = clean_zs(data)
    # mu = data.mean()
    mu = np.median(data)
    try:
        c4 = np.sqrt(2/(N - 1)) * math.gamma(N/2) / math.gamma((N - 1)/2)
    except OverflowError: # N too big
        c4 = 1 - 1/4/N - 7/32/(N**2) - 19/128/(N**3)
    # std = data.std(ddof=1) / c4
    std = 1.4826 * mad(data) / c4
    # std = np.sqrt(((data - mu) ** 2).sum() / (data.size - 1))
    if std == 0:
        std = 0.000000000000001

    return std


def std(data, axis=None):
    return transform(data, _std, axis=axis, mp=False)


def mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation
    """
    arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    med = np.median(arr)
    return np.median(np.abs(arr - med))


def _standardize(data):
    assert len(data.shape) == 1
    N = len(data)
    data = clean_zs(data)
    mu = data.mean()
    # mu = np.median(data)
    try:
        c4 = np.sqrt(2/(N - 1)) * math.gamma(N/2) / math.gamma((N - 1)/2)
    except OverflowError: # N too big
        c4 = 1 - 1/4/N - 7/32/(N**2) - 19/128/(N**3)
    std = data.std(ddof=1) / c4
    # std = 1.4826 * mad(data) / c4
    # std = np.sqrt(((data - mu) ** 2).sum() / (data.size - 1))
    if std == 0:
        std = 0.000000000000001

    return (data - mu) / std