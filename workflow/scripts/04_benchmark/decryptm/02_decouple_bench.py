#!/usr/bin/env python
# coding: utf-8

# In[1]:


if 'snakemake' in locals():
   input_file = snakemake.input[0]
   meta_data = snakemake.input[1]
   output_file = snakemake.output[0]
else:
   input_file = '../../../results/decryptm/benchmark_scores/ulm-R2_pEC50-jhonson.csv'
   meta_data = '../../../results/decryptm/benchmark_scores/obs_ulm-R2_pEC50-jhonson.csv'
   output_file = '../../../results/decryptm/benchmark/bench_ulm-R2_pEC50-jhonson.csv'


# In[21]:


import pandas as pd
import numpy as np
import os
import re

import decoupler as dc


# In[3]:


scores = pd.read_csv(input_file)
scores = pd.DataFrame(scores)


# In[23]:


filename = os.path.basename(input_file)

# Use regular expression to extract the part before "-R2"
match = re.search(r'(.+)-R2', filename)
mth = match.group(1)
mth


# In[4]:


mthds = np.array([mth])


# In[5]:


exps = scores.iloc[:, 0].to_numpy()
scores = scores.iloc[:, 1:]


# In[6]:


scores.replace(0, np.nan, inplace=True)


# In[7]:


srcs = np.array(scores.columns)


# In[8]:


scores.index = exps
scores.columns = srcs


# In[9]:


res = {mthds[0]: scores}


# In[10]:


obs = pd.read_csv(meta_data)


# In[11]:


obs['perturb'] = obs['perturb'].str.split(';')


# In[12]:


obs.set_index('sample', inplace=True)


# In[13]:


bench_res = dc.get_performances(res, obs, groupby=None, by='experiment', metrics=['auroc', 'auprc', 'mcauroc', 'mcauprc'], n_iter=1000, min_exp=1)


# In[ ]:


pd.DataFrame.to_csv(bench_res, output_file)

