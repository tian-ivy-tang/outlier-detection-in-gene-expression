import os
import pandas as pd

import helpers as h
import fast_zscore_estimation as fast
import optht_svd_zs as svd

from datetime import datetime
from datetime import timedelta


# pls provide dataset file_name including prefix of "/counts/"
# extraction of file_base_name (e.g. 'Geuvadis') is done within functions
file_name = 'counts/GTExSkinSmall.tsv'

# # Use the user's system login name as a unique identifier or a timestamp
unique_id = os.getlogin()
# os.environ.get('USERNAME')
# import getpass / getpass.getuser()

def run_pipeline(file_name):
    start_time = datetime.now()

    file_base_name = os.path.split(file_name)[1]
    file_base_name = os.path.splitext(file_base_name)[0]

    dfz_fast = fast.run(file_name)
    s_dimension, col_nb = svd.process(dfz_fast, file_base_name, unique_id)

    time_taken = round(timedelta.total_seconds(datetime.now() - start_time) / 60.0, 2)

    # export summary statistics: dataset name, target rank & number of columns to csv file
    # summary_stats = pd.DataFrame(columns = ['Dataset', 'Time', 'Q', 'Column_nb'])
    # summary_stats.loc[0] = [file_base_name, time_taken, s_dimension, col_nb]  
    # summary_stats_filename = "results/summary_stats_" + file_base_name + '_Outsingle_' + unique_id + '.csv'
    # h.save_df_to_csv(summary_stats, summary_stats_filename, index = False)

    # export summary_stats
    summary_stats = pd.read_csv('results/summary_stats.csv', header=0, index_col=False)
    summary_stats.loc[len(summary_stats)] = [file_base_name, 'OutSingle', time_taken, s_dimension, col_nb] 
    h.save_df_to_csv(summary_stats, 'results/summary_stats.csv', index = False)


run_pipeline(file_name)