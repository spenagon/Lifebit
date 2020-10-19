
"""
Spyder Editor

This is a temporary script file.
"""

from pyspark.sql.types import Row, DoubleType
from datetime import datetime
from pyspark.sql import SparkSession
from pyspark.sql.window import Window
import pyspark.sql.functions as func
from pyspark import SparkConf, SparkContext
import pandas as pd
import threading
import os
import numpy as np
from multiprocessing.pool import ThreadPool

# import custom functions from other scripts
import main_functions
import DQ


begin_time = datetime.now()

# set configuration settings for spark
conf=SparkConf().setMaster('local[*]').setAppName("GEL data pipeline")
conf.set("spark.scheduler.mode",'FAIR')
conf.set("spark.driver.memory", "100g")
sc=SparkContext(conf=conf)

# start the spark session
spark = SparkSession.builder.appName("GEL data pipeline").getOrCreate()

# check spark configuration settings
print(spark.sparkContext.getConf().getAll())

'''
file_list = ['participant.tsv', 'COVID_data_synthetic.tsv', 'ons.tsv' ,'aggregate_gvcf_sample_stats.tsv',
             'cancer_participant_disease.tsv', 'cancer_participant_disease.tsv', 'death_details.tsv',
             'genome_file_paths_and_types.tsv', 'hes_ae.tsv', 'hes_apc.tsv', 'hes_cc.tsv', 'hes_op.tsv',
             'rare_diseases_family.tsv', 'rare_diseases_participant_disease.tsv', 'rare_diseases_participant_phenotype.tsv',
             'rare_diseases_pedigree.tsv', 'rare_diseases_pedigree_member.tsv', 'sequencing_report.tsv']
'''
# list of files to read in.
file_list = ['participant', 'death_details', 'rare_diseases_family', 'rare_diseases_pedigree', 'rare_diseases_pedigree_member',
                'rare_diseases_participant_disease', 'rare_diseases_participant_phenotype', 'cancer_participant_disease']

# import the data dictionary
data_dict = pd.read_csv(os.getcwd() + '/MetadataLifebit.csv')


# read files and union into one data frame
# main_functions.name_convert uses the data dictionary to update the names of all the columns
for file in file_list:
    df_file = spark.read.load('gel_' + file + '_20200402.tsv', format="csv", sep="\t", inferSchema="true", header="true")
    df_file = main_functions.name_convert(df_file, data_dict)
    try:
        df_out = main_functions.customUnion(df_file, df_out)
    except:
        df_out = df_file


# This step is used when merging data to existing data in the system
#df_out = main_dunctions.customUnion(new_merge, merge_all)

# Only keep consenting participants. Think about how to do this when not dealing with participant info
#df_out = df_out.filter(func.col('21') == 'Consenting')

# Update gender terms
df_out = DQ.fix_gender(data_dict, df_out)

# Collect Ids of each type
type_dict = {}
for d_type in data_dict['Type'].unique().tolist():
    list_ = data_dict.loc[data_dict['Type'] == d_type]['FieldID'].astype(str).tolist()
    type_dict[d_type] = list_

# Collect Ids for each value type
value_type_dict = {}
for v_type in data_dict['ValueType'].unique().tolist():
    list_ = data_dict.loc[data_dict['ValueType'] == v_type]['FieldID'].astype(str).tolist()
    value_type_dict[v_type] = list_


# Take 'string' out of all columns that are histogram type.
new_df = DQ.hist_strip(df_out, value_type_dict)

# Use the data type column to make sure the types of all columns are correct. This also changes 'i' to string.
new_df = DQ.update_type(new_df, type_dict)

'''
        here we are turning the data from current long format to wide.
        Thread pooling has been implemented here because we need to do the process on all the columns. So we are running the co-currently.
'''
#initialise some things
arrays = []
key = 'i'
data = new_df

# define the pool - chose 8 because we have 8 cores
pool = ThreadPool(8)

dataout=data.select(key).distinct()

def run_explode_df(col_name, dataout = dataout):
    array_no, data2 = main_functions.explode_df(data.select(key,col_name), key, col_name)
    #dataout=dataout.join(data2,dataout[key]==data2[key]).drop(data2[key])
    arrays.append({col_name:array_no})
    print('col done: ' + col_name)
    return data2

col_list = data.columns
col_list.remove('i')
dataout2 = pool.map(run_explode_df, col_list)
print(dataout2)


for df in dataout2:
        dataout = dataout.join(df,[key], 'left')
#dataout = dataout.join(dataout2[1], [key], 'left')
dataout.write.format('json').mode('overwrite').save(os.getcwd() + '/participant_json_20201005')

pool.close()
pool.join()


print(begin_time - datetime.now())
'''


def run_phenotypevaluesRV(col_, dataout = dataout):
        main_functions.phenotypevaluesRV(dataout, col_)


col_list2 = dataout.columns
col_list2.remove('i')
pool2 = ThreadPool(32)

pool2.map(run_phenotypevaluesRV, col_list2)
pool2.close()
pool2.join()
'''
'''
# define the pool - chose 8 because we have 8 cores
# define the pool - chose 8 because we have 8 cores
pool3 = ThreadPool(8)

df_cols = [i for i in dataout.columns if i.startswith('f')]
distinct_cols = set([col_.split('f')[1].split('i')[0] for col_ in df_cols])

def run_phenotypefields(col_name, dataout = dataout, data_dict = data_dict, df_cols = df_cols):
    data3 = main_functions.phenotypefields(dataout, data_dict, df_cols, col_name)
    #dataout=dataout.join(data2,dataout[key]==data2[key]).drop(data2[key])
    print('col done: ' + col_name)
    return data3

dataout3 = pool3.map(run_phenotypefields, list(distinct_cols)[1:10])
print(dataout3)


for df in dataout3:
    try:
        df2 = df2.append(df, ignore_index = True)
    except:
        df2 = df

df_pf = df2.replace(np.nan, '')

df_pf.to_json(os.getcwd() + '/phenotypefields_20201005.json', orient="records")


pool3.close()
pool3.join()



'''
    

