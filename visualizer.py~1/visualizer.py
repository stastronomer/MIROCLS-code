from astroquery.vizier import Vizier as Vizier
import astropy.units as u
import astropy.coordinates as coord
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
import warnings
warnings.filterwarnings("ignore")
from astropy.coordinates import SkyCoord
import io
import requests
from astropy.table import Table
from PIL import Image

from io import BytesIO


from matplotlib import gridspec
from matplotlib.patches import Circle
cluster = pd.read_csv('/home/guptaa/GLIMPSE/OPTICS Cluster.csv')
annulus_log = {'Glon': [], 'Glat': [], 'annulus_color_med': [], 'annulus_color_rms': [], 
              'cluster_color_med': [], 'cluster_color_rms': [], 'c2a_ratio': []}
for i in range(cluster.shape[0]):
    #if i==1:
    #   break
    ra, dec, deg = cluster.iloc[i]['Glon'], cluster.iloc[i]['Glat'], cluster.iloc[i]['size']
    
    v = Vizier(columns=['Glon','Glat', 'RAJ200', 'DEJ200', 'Jmag', 'Hmag', 'Kmag',
                       'e_Jmag', 'e_Hmag', 'e_Kmag', '3.6mag', 'e_3.6mag', '4.5mag', 'e_4.5mag',
                       '5.8mag', 'e_5.8mag', '8.0mag', 'e_8.0mag'])
    v.ROW_LIMIT = -1
    #, column_filters={"e_3.6mag":"<0.2", "e_4.5mag":"<0.2", 'e_Kmag':'<0.5'}
    
    
    result = v.query_region(coord.SkyCoord(ra, dec,
                                                unit=(u.deg, u.deg),
                                                frame='galactic'),
                            radius=10**0.5*deg*u.deg, inner_radius=3*deg*u.deg,
                            catalog=["II/293/glimpse"])
    
    result_2 = v.query_region(coord.SkyCoord(ra, dec,
                                                unit=(u.deg, u.deg),
                                                frame='galactic'),
                            radius=deg*u.deg,
                            catalog=["II/293/glimpse"])
    
    data=result['II/293/glimpse']
    data_2 = result_2['II/293/glimpse']
    #fig=plt.figure(figsize=(3,10))
    #ax1 = plt.subplot(211)
    #ax2 = plt.subplot(212)
    #fig.subplots_adjust(hspace=0)
    #ax1.scatter(data['_3.6mag']-data['_4.5mag'], data['_3.6mag'],s=5, label='3X Cluster Background',c='black')
    #ax1.scatter(data_2['_3.6mag']-data_2['_4.5mag'], data_2['_3.6mag'],s=10, label='cluster members',c='r')
    #plt.xlabel('[3.6]-[4.5]',fontsize=18)
    #ax1.ylabel('[3.6]',fontsize=18)
    #ax1.gca().invert_yaxis()
    #ax2.scatter(data['Kmag']-data['_3.6mag'], data['Kmag'],s=5, label='3X Cluster Background',c='black')
    #ax2.scatter(data_2['Kmag']-data_2['_3.6mag'], data_2['Kmag'],s=5, c='r')
    #ax1.get_shared_x_axes().join(ax1, ax2)
    #ax1.set_ylabel('[3.6]',fontsize=18)
    #ax2.set_ylabel('[K]',fontsize=18)

    #ax2.set_xlabel('[3.6]-[4.5] / [K]-[3.6]',fontsize=12)
    #ax1.invert_yaxis()
    #ax2.invert_yaxis()
    #ax1.legend()

    
    fig=plt.figure(constrained_layout=True)
    fig.set_figheight(16)
    fig.set_figwidth(28)
    
    spec = gridspec.GridSpec(ncols=2, nrows=2,
                         width_ratios=[1, 1], wspace=0.5,
                         hspace=0.01, height_ratios=[1, 1])
    ax1 = fig.add_subplot(spec[0])
    ax2 = fig.add_subplot(spec[2])
    ax3 = fig.add_subplot(spec[1])
    ax4 = fig.add_subplot(spec[3])
    fig.subplots_adjust(hspace=0)
    ax1.scatter(data['_3.6mag']-data['_4.5mag'], data['_3.6mag'],s=10, label='circular annulus',c='black')
    ax1.scatter(data_2['_3.6mag']-data_2['_4.5mag'], data_2['_3.6mag'],s=20, label='cluster members',c='r')
    #plt.xlabel('[3.6]-[4.5]',fontsize=18)
    #ax1.ylabel('[3.6]',fontsize=18)
    #ax1.gca().invert_yaxis()
    ax2.scatter(data['Kmag']-data['_3.6mag'], data['Kmag'],s=10, label='3X Cluster Background',c='black')
    ax2.scatter(data_2['Kmag']-data_2['_3.6mag'], data_2['Kmag'],s=20, c='r')
    ax1.get_shared_x_axes().join(ax1, ax2)
    ax1.set_ylabel('[3.6]',fontsize=18)
    ax2.set_ylabel('[K]',fontsize=18)

    ax2.set_xlabel('[3.6]-[4.5] / [K]-[3.6]',fontsize=12)
    ax1.invert_yaxis()
    ax2.invert_yaxis()
    ax1.legend(fontsize=25)
    
    def rms(x):
        return np.sqrt(x.dot(x)/x.size)
    
    annulus_log['Glon'].append(ra)
    annulus_log['Glat'].append(dec)
    annulus_log['annulus_color_med'].append(np.median(data['_3.6mag'])-np.median(data['_4.5mag']))
    annulus_log['annulus_color_rms'].append(rms(data['_3.6mag']-data['_4.5mag']))
    annulus_log['cluster_color_med'].append(np.median(data_2['_3.6mag'])-np.median(data_2['_4.5mag']))
    annulus_log['cluster_color_rms'].append(rms(data_2['_3.6mag']-data_2['_4.5mag']))
    annulus_log['c2a_ratio'].append(round(len(data_2['_3.6mag'])/len(data['_3.6mag']),3))
    
    #plt.show()
    
    
    ra_gal, dec_gal, size = cluster['Glon'][i], cluster['Glat'][i], cluster['size'][i]

    gc = SkyCoord(l=ra_gal*u.degree, b=dec_gal*u.degree, frame='galactic', unit='deg')
    ra, dec = gc.fk5.ra.value, gc.fk5.dec.value
    url = 'https://sky.esa.int/esasky-tap/skyimage?target='+str(ra)+'%20'+str(dec)+'&fov=0.1&aspectratio=1.9437751004016063&hips=Spitzer%20cold%20SEIP%20IRAC-1-3-4%20RGB%20bright&size=512&fmt=PNG&_=1652798602280'
    response = requests.get(url)
    #img = Image.open(BytesIO(response.content))
    img=Image.open(BytesIO(response.content))

    arr=np.array(img)

    url_2mass = 'https://sky.esa.int/esasky-tap/skyimage?target='+str(ra)+'%20'+str(dec)+'&fov=0.1&aspectratio=1.9437751004016063&hips=2MASS%20color%20JHK&size=512&fmt=PNG&_=1652798602280'
    response_2mass = requests.get(url_2mass)
    #img = Image.open(BytesIO(response.content))
    img_2mass=Image.open(BytesIO(response_2mass.content))

    arr_2mass=np.array(img_2mass)





    avg_2mass = (arr_2mass[:,:,0]+arr_2mass[:,:,1]+arr_2mass[:,:,2])/3

    #plt.imshow(arr, extent=(ra_gal-0.05/2,ra_gal+0.05/2,dec_gal-0.05/2,dec_gal+0.05/2))
    #plt.xlim(ra_gal-cluster['size'][0]*3, ra_gal+cluster['size'][0]*3)
    #plt.ylim(dec_gal-cluster['size'][0]*3, dec_gal+cluster['size'][0]*3)


    
    avg = (arr[:,:,0]+arr[:,:,1]+arr[:,:,2])/3
    


    bg=2

    im=ax3.imshow(avg, extent=(ra_gal-0.05/2,ra_gal+0.05/2,dec_gal-0.05/2,dec_gal+0.05/2), aspect=1/2)
    ax3.set_xlim(ra_gal-cluster['size'][0]*bg, ra_gal+cluster['size'][0]*bg)
    ax3.set_ylim(dec_gal-cluster['size'][0]*bg, dec_gal+cluster['size'][0]*bg)
    ax3.set_title('GLIMPSE', fontsize=18)
    drawObject = Circle((ra_gal, dec_gal), cluster['size'][i], fill=False, color='r')
    ax3.add_patch(drawObject)
    drawObject = Circle((ra_gal, dec_gal), 3*cluster['size'][i], fill=False, color='white')
    ax3.add_patch(drawObject)
    drawObject = Circle((ra_gal, dec_gal), 10**.5*cluster['size'][i], fill=False, color='white')
    ax3.add_patch(drawObject)
    
    #ax3.Circle((ra_gal, dec_gal), 3*cluster['size'][i], fill=False, c='black' )
    #ax3.Circle((ra_gal, dec_gal), 10**.5*cluster['size'][i], fill=False, c='black' )
    #plt.colorbar(im,ax=ax3)

    im_2mass=ax4.imshow(avg_2mass, extent=(ra_gal-0.05/2,ra_gal+0.05/2,dec_gal-0.05/2,dec_gal+0.05/2), aspect=1/2)
    ax4.set_xlim(ra_gal-cluster['size'][0]*bg, ra_gal+cluster['size'][0]*bg)
    ax4.set_ylim(dec_gal-cluster['size'][0]*bg, dec_gal+cluster['size'][0]*bg)
    ax4.set_title('2MASS', fontsize=18)
    drawObject = Circle((ra_gal, dec_gal), cluster['size'][i], fill=False, color='r')
    ax4.add_patch(drawObject)
    drawObject = Circle((ra_gal, dec_gal), 3*cluster['size'][i], fill=False, color='black')
    ax4.add_patch(drawObject)
    drawObject = Circle((ra_gal, dec_gal), 10**.5*cluster['size'][i], fill=False, color='black')
    ax4.add_patch(drawObject)
    
    #plt.colorbar(im_2mass,ax=ax4)
    #ax2.imshow(arr[:,:,2])
    
    
    
    fig.suptitle('Coords: '+str(round(ra_gal,5))+', '+str(round(dec_gal,5))+'\n #sources in cluster: '
                 +str(len(data_2['_3.6mag']))+'\n #sources in annlulus: '+str(len(data['_3.6mag']))+
                 '\n cluster size: '+str(round(cluster['size'][i]*3600,1))+"'' ",
                fontsize=25)
    
    
    
    
    
    plt.savefig('/home/guptaa/GLIMPSE/split_15.0_35.0/img'+str(round(len(data_2['_3.6mag'])/len(data['_3.6mag']),2))+str(round(ra_gal,5))+'.jpeg')
    
    
    plt.close()
    print(100*(i+1)/cluster.shape[0], '% done...')

print('100% done')
annulus_log = pd.DataFrame(annulus_log)
annulus_log.to_csv('logfile.csv')
