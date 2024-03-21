#!/usr/bin/env python
# coding: utf-8

# In[1]:


if 'snakemake' in locals():
   input_file = snakemake.input[0]
   meta_data = snakemake.input[1]
   output_file = snakemake.output[0]
else:
   input_file = 'results/03_benchmark/hernandez/01_input_bench_subset/johnson/KSEA_z-phosphositeplus_johnson15.csv'
   meta_data = 'results/03_benchmark/hernandez/01_input_bench_subset/johnson/obs_KSEA_z-phosphositeplus_johnson15.csv'
   output_file = 'results/03_benchmark/hernandez/02_benchmark_res_subset/johnson/phosphositeplus_johnson15/bench_KSEA_z-phosphositeplus_johnson15.csv'


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

