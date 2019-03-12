import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as plt
from matplotlib import pyplot as mp
import glob, os
file_path = "C:/Users/nealp/Dropbox (Partners HealthCare)/PNOIT2_singlecell_data/Single_cell_TCR_for_Neal"
os.chdir(file_path)

files = glob.glob("dist*.csv")
top_NN_5 = dict()

for i in files:
    data = pd.read_csv(file_path + "/" + i, index_col = 0)
    
    # Make the data into a matrix
    data = data.to_numpy()
    
    # Look at nearest neighbors for each sequence
    avg_list = []
    for j in range(0, np.size(data, 0)):
        row = data[j, :]
        sorted = np.sort(row)
        # Get rid of the leading 0
        sorted = sorted[1:len(sorted)] 
        # Get top neighbors 
        sorted = sorted[0:round(len(sorted) / (len(row) / 5))]
        # Get weights
        #weights = np.random.uniform(5, 0, len(sorted))
        #weights = -np.sort(-weights)
        # Get average
        avg = np.average(sorted)
        # Append to list
        avg_list.append(avg)
        
        top_NN_5["topNN_" + str(i)] = avg_list
        
ax = mp.axes()
ax.set_title("top 5 neighbors")
mp.xlabel("distance")
sns.distplot(top_NN_5["topNN_dist_mtx_96.csv"], hist = False, label = "96")
sns.distplot(top_NN_5["topNN_dist_mtx_95.csv"], hist = False, label = "95")
sns.distplot(top_NN_5["topNN_dist_mtx_90.csv"], hist = False, label = "90")
sns.distplot(top_NN_5["topNN_dist_mtx_93.csv"], hist = False, label = "93")
sns.distplot(top_NN_5["topNN_dist_mtx_97.csv"], hist = False, label = "97")
mp.savefig(file_path + '/figures/dens.plot.top5.png')

top_NN_10 = dict()
for i in files:
    data = pd.read_csv(file_path + "/" + i, index_col = 0)
    
    # Make the data into a matrix
    data = data.to_numpy()
    
    # Look at nearest neighbors for each sequence
    avg_list = []
    for j in range(0, np.size(data, 0)):
        row = data[j, :]
        sorted = np.sort(row)
        # Get rid of the leading 0
        sorted = sorted[1:len(sorted)] 
        # Get top neighbors 
        sorted = sorted[0:round(len(sorted) / (len(row) / 10))]
        # Get weights
        #weights = np.random.uniform(5, 0, len(sorted))
        #weights = -np.sort(-weights)
        # Get average
        avg = np.average(sorted)
        # Append to list
        avg_list.append(avg)
        
        top_NN_10["topNN_" + str(i)] = avg_list

ax = mp.axes()
ax.set_title("top 10 neighbors")
mp.xlabel("distance")
sns.distplot(top_NN_10["topNN_dist_mtx_96.csv"], hist = False, label = "96")
sns.distplot(top_NN_10["topNN_dist_mtx_95.csv"], hist = False, label = "95")
sns.distplot(top_NN_10["topNN_dist_mtx_90.csv"], hist = False, label = "90")
sns.distplot(top_NN_10["topNN_dist_mtx_93.csv"], hist = False, label = "93")
sns.distplot(top_NN_10["topNN_dist_mtx_97.csv"], hist = False, label = "97")
mp.savefig(file_path + '/figures/dens.plot.top10.png')