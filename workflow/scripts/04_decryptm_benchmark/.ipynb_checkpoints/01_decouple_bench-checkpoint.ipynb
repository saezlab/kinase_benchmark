{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 4,
   "id": "7d9be570-8c9e-4653-b146-058ca7f89c7c",
   "metadata": {},
   "outputs": [],
   "source": [
    "if 'snakemake' in locals():\n",
    "   input_file = snakemake.input[0]\n",
    "   meta_data = snakemake.input[1]\n",
    "   output_file = snakemake.output[0]\n",
    "else:\n",
    "   input_file = '../../../results/decryptm/benchmark_scores/ulm_R2_pEC50_GPS.csv'\n",
    "   mth = 'ulm'\n",
    "   meta_data = '../../../results/decryptm/benchmark_scores/ulm_R2_pEC50_GPS_obs.csv'\n",
    "   output_file = '../../../results/decryptm/benchmark/ulm_R2_pEC50_GPS_bench.csv'"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "id": "1b4137b0-5ed0-498d-8597-6089e666deb0",
   "metadata": {},
   "outputs": [],
   "source": [
    "import pandas as pd\n",
    "import numpy as np\n",
    "\n",
    "import decoupler as dc"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "id": "71b37b0b-52e9-4fdb-bb73-786c8c63836d",
   "metadata": {},
   "outputs": [],
   "source": [
    "scores = pd.read_csv(input_file)\n",
    "scores = pd.DataFrame(scores)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 7,
   "id": "bc480df6-fb55-42b1-b0f3-30de34fb4712",
   "metadata": {},
   "outputs": [],
   "source": [
    "mthds = np.array([mth])"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 8,
   "id": "eb227c64-7a4c-45ae-8466-8179c305e174",
   "metadata": {},
   "outputs": [],
   "source": [
    "exps = scores.iloc[:, 0].to_numpy()\n",
    "scores = scores.iloc[:, 1:]"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 9,
   "id": "1fcc8772-416c-4181-a720-31c6d5f6eb46",
   "metadata": {},
   "outputs": [],
   "source": [
    "scores.replace(0, np.nan, inplace=True)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 10,
   "id": "762530f6-d435-4413-b12a-88f89d2a01eb",
   "metadata": {},
   "outputs": [],
   "source": [
    "srcs = np.array(scores.columns)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 11,
   "id": "c50443ba-1298-42e0-8d35-08e4c088673e",
   "metadata": {},
   "outputs": [],
   "source": [
    "scores.index = exps\n",
    "scores.columns = srcs"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 12,
   "id": "d917bdda-5358-412a-a474-4a58f99182b2",
   "metadata": {},
   "outputs": [],
   "source": [
    "res = {mthds[0]: scores}"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 13,
   "id": "618b7397-8173-4cb6-848d-05c9b594894e",
   "metadata": {},
   "outputs": [],
   "source": [
    "obs = pd.read_csv(meta_data)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 14,
   "id": "fa5f1ded-c761-4ba4-b6ca-97bb40c83c68",
   "metadata": {},
   "outputs": [],
   "source": [
    "obs['perturb'] = obs['perturb'].str.split(';')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 15,
   "id": "d7ba7cd1-b8dc-47dd-82ca-102051f423b6",
   "metadata": {},
   "outputs": [],
   "source": [
    "obs.set_index('sample', inplace=True)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 16,
   "id": "f1925dbe-fd52-4580-ae57-e4a65193328e",
   "metadata": {},
   "outputs": [],
   "source": [
    "bench_res = dc.get_performances(res, obs, groupby=None, by='experiment', metrics=['auroc', 'auprc', 'mcauroc', 'mcauprc'], n_iter=1000)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 17,
   "id": "3133ed71-377b-4735-b6c2-bf30bf6effda",
   "metadata": {},
   "outputs": [],
   "source": [
    "pd.DataFrame.to_csv(bench_res, output_file)"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "benchmark",
   "language": "python",
   "name": "benchmark"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.11.4"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
