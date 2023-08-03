#!/usr/bin/env python
# coding: utf-8

# In[4]:


if 'snakemake' in locals():
   input_file = snakemake.input[0]
   meta_data = snakemake.input[1]
   output_file = snakemake.output[0]
else:
   input_file = '../../../results/decryptm/benchmark_scores/ulm_R2_pEC50_GPS.csv'
   meta_data = '../../../results/decryptm/benchmark_scores/ulm_R2_pEC50_GPS_obs.csv'
   output_file = '../../../results/decryptm/benchmark/bench_ulm_R2_pEC50_GPS.csv'


# In[5]:


import pandas as pd
import numpy as np

import decoupler as dc


# In[6]:


scores = pd.read_csv(input_file)
scores = pd.DataFrame(scores)


# In[7]:

mth = input_file.split("/")[3].split("-")[0]
print(mth)
mthds = np.array([mth])


# In[8]:


exps = scores.iloc[:, 0].to_numpy()
scores = scores.iloc[:, 1:]


# In[9]:


scores.replace(0, np.nan, inplace=True)


# In[10]:


srcs = np.array(scores.columns)


# In[11]:


scores.index = exps
scores.columns = srcs


# In[12]:


res = {mthds[0]: scores}


# In[13]:


obs = pd.read_csv(meta_data)


# In[14]:


obs['perturb'] = obs['perturb'].str.split(';')


# In[15]:


obs.set_index('sample', inplace=True)


# In[16]:


bench_res = dc.get_performances(res, obs, groupby=None, by='experiment', metrics=['auroc', 'auprc', 'mcauroc', 'mcauprc'], n_iter=1000)


# In[17]:


pd.DataFrame.to_csv(bench_res, output_file)

