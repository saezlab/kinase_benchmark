#!/usr/bin/env python
# coding: utf-8

# In[1]:


if 'snakemake' in locals():
   input_file = snakemake.input[0]
   meta_data = snakemake.input[1]
   output_file = snakemake.output[0]
else:
   input_file = '../../../../results/hernandez/benchmark_files/KSEA_z-ptmsigdb.csv'
   meta_data = '../../../../results/hernandez/benchmark_files/obs_KSEA_z-ptmsigdb.csv'
   output_file = '../../../../results/hernandez/benchmark_res/bench_KSEA_z-ptmsigdb.csv'


# In[2]:


import pandas as pd
import numpy as np
import os
import re

import decoupler as dc


# In[3]:


scores = pd.read_csv(input_file)
scores = pd.DataFrame(scores)


# In[4]:


filename = os.path.basename(input_file)

# Use regular expression to extract the part before "-R2"
match = re.search(r'(.+)-', filename)
mth = match.group(1)
mth


# In[5]:


mthds = np.array([mth])


# In[6]:


exps = scores.iloc[:, 0].to_numpy()
scores = scores.iloc[:, 1:]


# In[7]:


scores.replace(0, np.nan, inplace=True)


# In[8]:


srcs = np.array(scores.columns)


# In[9]:


scores.index = exps
scores.columns = srcs


# In[10]:


res = {mthds[0]: scores}


# In[11]:


obs = pd.read_csv(meta_data)


# In[12]:


obs['perturb'] = obs['perturb'].str.split(';')


# In[13]:


obs.set_index('sample', inplace=True)


# In[15]:


bench_res = dc.get_performances(res, obs, groupby=None, by='experiment', metrics=['auroc', 'auprc', 'mcauroc', 'mcauprc'], n_iter=1000, min_exp=1)


# In[16]:


pd.DataFrame.to_csv(bench_res, output_file)

