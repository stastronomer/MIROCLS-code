try:
  import numpy as np
  import matplotlib.pyplot as plt
  import math
  import os
  import pandas as pd
  from astropy.io import fits
  import csv
  from astropy.table import Table
  from matplotlib import colors
  from glob import glob

  from sklearn.cluster import DBSCAN, KMeans, MiniBatchKMeans, OPTICS
  from sklearn import metrics
  from sklearn.datasets import make_blobs
  from sklearn.preprocessing import StandardScaler, MinMaxScaler
  import matplotlib
  from sklearn.neighbors import NearestNeighbors

except ModuleNotFoundError:
  print("Some modules not installed, use: pip install <module-name>, critical modules: \
          sklearn, astroquery, matplotlib, astropy")




def cluster_field(data, xcen, ycen, r):
    '''
    Returns all the sources which lies within radius, r, of the given central points: (xcen, ycen)
    '''
    r = r/60
    return data[((data['Glon']-xcen)**2+(data['Glat']-ycen)**2<r**2)].shape[0]

def size(c,data1):
    '''
    Returns the rough radius of a cluster by averaging two farthest point from the cluster
    '''
    dist=[]

    l=data1.shape[0]
    for i in range(l):
        a=[data1.iloc[i]['Glon'], data1.iloc[i]['Glat']]
        dist.append(np.linalg.norm(np.array(c)-a))
    return max(dist)



def search():
    cluster_db = {'Glon': [], 'Glat': [], '3_6mag': [], 'e_3_6mag': [], '4_5mag': [], 'e_4_5mag': [], 'size': [], 'sources_cluster': [], 'tile': [],'color_bin': [], 'sources_20_arcmin': [], 'sources_tile': [], 'significance': []}




    #READ FILES
    print('enter folder path to tiles: ')
    path = input()+'/' 
    print('enter lower color limit (mag): ')
    lim1 = float(input())
    print('enter upper color limit (mag): ')
    lim2 = float(input())


    path1=path+'*.csv'
    data = glob(path1)
    print('OPTICS clusering algorithm will be run with tuned paramters')
    print('optimized params are: max_eps=0.05, min_samples=12, xi=0.05')
    counter = 1

    #CREATE OUTPUT DIRECTORY
    parent_dir = os.getcwd()
    directory = path+str(lim1)+'_'+str(lim2)
    path_op = os.path.join(parent_dir, directory)
    if not os.path.exists(path_op):
        os.mkdir(path_op)

    #RUNNING THE ALGORITHM
    for file1 in data:
        print('reading...', file1)
        df = pd.read_csv(file1)
        if df.shape[0]!=0:
            s=[]
            n = len(df['3_6mag'])
            for i in range(n):
                if math.isnan(df['3_6mag'][i]) and math.isnan(df['4_5mag'][i]):
                        s.append(i)
            data_good = df.drop(s)
            data_good=data_good.fillna(999)
            data_good=data_good.reset_index()
            data_good=data_good.drop(['index'],axis=1)
            idx=[]
            for i in range(data_good.shape[0]):
                if data_good['e_3_6mag'][i]>0.2 or data_good['e_4_5mag'][i]>0.2:
                    idx.append(i)

            data_good = data_good.drop(idx)
            data_good = data_good.reset_index()
            data_sdm_cs = data_good[(data_good['3_6mag']-data_good['4_5mag']>=lim1) & (data_good['3_6mag']-data_good['4_5mag']<=lim2)]
            X = data_sdm_cs[['Glon', 'Glat']]

            if X.shape[0]>12:
                
                model = OPTICS(max_eps=0.05, min_samples=12, xi=0.05)

                yhat = model.fit_predict(X)
    # retrieve unique clusters
                clusters = np.unique(yhat)
    # create scatter plot for samples from each cluster
                labels = model.labels_
                unique_labels = set(labels)
                colors=[plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]
                plt.figure(figsize=(10,10))
                plt.scatter(data_sdm_cs['Glon'],data_sdm_cs['Glat'],c=labels,cmap=matplotlib.colors.ListedColormap(colors),s=5)
                plt.title('OPTICS',fontsize=25)
                plt.xlabel('Gal Lon',fontsize=14)
                plt.ylabel('Gal Lat',fontsize=14)
                #plt.savefig(str(data_sdm_cs['Glon'])) #uncomment to save plots, but slower processing
                plt.close()
                print(path_op)
                print(file1[0:-4])

                data_sdm_cs['labels']=model.labels_
            
                if len(data_sdm_cs['labels'])>1:
                    for clidx in set(labels):
                        if clidx!=-1:
                            idx_cluster = np.where(data_sdm_cs['labels']==clidx)
                            data_clust = data_sdm_cs[data_sdm_cs.labels == clidx]
                            s20min = cluster_field(data_sdm_cs,np.median(data_clust['Glon']),np.median(data_clust['Glat']),10)
                            s_cluster = data_sdm_cs['labels'].value_counts()[clidx]
                            significance = s_cluster/(df.shape[0]*np.pi*size([np.median(data_clust['Glon']),np.median(data_clust['Glat'])],data_clust)**2)**0.5
                            cluster_db['Glon'].append(np.median(data_clust['Glon']))
                            cluster_db['Glat'].append(np.median(data_clust['Glat']))
                            cluster_db['3_6mag'].append(np.average(data_clust['3_6mag']))
                            cluster_db['e_3_6mag'].append(np.average(data_clust['e_3_6mag']))
                            cluster_db['4_5mag'].append(np.average(data_clust['4_5mag']))
                            cluster_db['e_4_5mag'].append(np.average(data_clust['e_4_5mag']))
                            cluster_db['size'].append(size([np.median(data_clust['Glon']),np.median(data_clust['Glat'])],data_clust))
                            cluster_db['tile'].append(file1)
                            cluster_db['color_bin'].append(str(lim1)+'-'+str(lim2))
                            cluster_db['sources_20_arcmin'].append(s20min)
                            cluster_db['sources_tile'].append(df.shape[0])
                            cluster_db['sources_cluster'].append(s_cluster)
                            cluster_db['significance'].append(significance)
        #print('check:', len(cluster_db['Glon']))
        print(100*counter/len(data),'% done...')
        counter=counter+1
    #plt.savefig('clust_run1.pdf')
    cluster_db = pd.DataFrame(cluster_db)
    #SAVE THE NEW POTENTIAL CLUSTER INFORMATION
    print('Saving as: ',path+'_'+str(lim1)+'_'+str(lim2)+'.csv')
    cluster_db.to_csv(path+'_'+str(lim1)+'_'+str(lim2)+'.csv')


    #print(path)
