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

  from sklearn.cluster import DBSCAN, KMeans, MiniBatchKMeans, OPTICS
  from sklearn import metrics
  from sklearn.preprocessing import StandardScaler, MinMaxScaler
  import matplotlib
  from sklearn.neighbors import NearestNeighbors

except ModuleNotFoundError:
  print("Some modules not installed, use: pip install <module-name>, critical modules: \
          sklearn, astroquery, matplotlib, astropy")

'''
This code will split the input data (as fits file) into smaller tiles of mentioned size or into $1deg^2$ tiles
if size not mentioned. Please input a two-dimensional fits file with (x,y) in degrees (l,b or ra,dec).
'''


def split(size=1):
  #Load Data
  print('enter file name: ')
  file = input()
  hdulist = fits.open(file)
  data = hdulist[1].data

  #Create destination folder
  print('creating new folder in current directory...')
  parent_dir = os.getcwd()
  directory = 'split_'+str(round(min(data['Glon'])))+'_'+str(round(max(data['Glon'])))
  path = os.path.join(parent_dir, directory)
  if not os.path.exists(path):  
    os.mkdir(path)


  #performing split operation
  temp1 = np.arange(min(data['Glon']),max(data['Glon']),size)
  temp2 = np.arange(min(data['Glat']),max(data['Glat']),size)

  l1, l2 = len(temp1), len(temp2)

  for i in range(l1):
    print(100*i/l1,'% completed...')
    for j in range(l2):
      df = []
      idx = []
      
      idx = np.where((data['Glon']>min(data['Glon'])+i) & (data['Glon']<min(data['Glon'])+i+1) & (data['Glat']>min(data['Glat'])+j) & (data['Glat']<min(data['Glat'])+j+1))
      df = pd.DataFrame(data[idx[0]])
      print("saving file as ", os.path.join(path,str(i)+'_'+str(j)+'.csv'))
      df.to_csv(os.path.join(path,str(i)+'_'+str(j)+'.csv'))
  print('100% done')
