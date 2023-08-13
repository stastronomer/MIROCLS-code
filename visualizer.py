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
import requests
from io import BytesIO
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import gridspec
from matplotlib.patches import Circle
import os


print("Input the csv file")
file = input()
cluster = pd.read_csv(file)

#CREATE OUTPUT DIRECTORY
parent_dir = os.getcwd()
directory = 'plot_images'
path_op = os.path.join(parent_dir, directory)
if not os.path.exists(path_op):
    os.mkdir(path_op)

print("Output directory set to: ", path_op)

# LOG FILE TO SAVE THE IMPORTANT INFORMATION ON POTENITAL CLUSTERS FOR FURTHER ANALYSIS, NOT IMPORTANT FOR PLOTTING
annulus_log = {'Glon': [], 'Glat': [], 'annulus_color_med': [], 'annulus_color_rms': [], 
              'cluster_color_med': [], 'cluster_color_rms': [], 'c2a_ratio': []}
clust_stat = cluster.copy()
clust_stat['overdensity'] = ""
for i in range(cluster.shape[0]):
    #if i==2:  #TO RUN A CHECK ON JUST ONE SOURCE UNCOMMENT if AND break
    #   break
    #enter l,b and size in degrees
    ra_gal, dec_gal, deg = cluster.iloc[i]['Glon'], cluster.iloc[i]['Glat'], cluster.iloc[i]['size']
    gc = SkyCoord(l=ra_gal*u.degree, b=dec_gal*u.degree, frame='galactic', unit='deg')
    ra, dec = gc.fk5.ra.value, gc.fk5.dec.value
    
    #load full data
    v = Vizier(columns=['Glon','Glat','_Glon', '_Glat' 'RAJ200', 'DEJ200', 'Jmag', 'Hmag', 'Kmag',
                       'e_Jmag', 'e_Hmag', 'e_Kmag', '3.6mag', 'e_3.6mag', '4.5mag', 'e_4.5mag',
                       '5.8mag', 'e_5.8mag', '8.0mag', 'e_8.0mag'])
    v.ROW_LIMIT = -1
    #, column_filters={"e_3.6mag":"<0.2", "e_4.5mag":"<0.2", 'e_Kmag':'<0.5'}
    
    xlim = 1.5*2.69**0.5*deg*u.deg
    result_1 = v.query_region(coord.SkyCoord(ra_gal, dec_gal,
                                                unit=(u.deg, u.deg),
                                                frame='galactic'),
                            radius=2*xlim,
                            catalog=["II/293/glimpse", "II/246/out"], cache=False)
    
    result = v.query_region(coord.SkyCoord(ra_gal, dec_gal,
                                                unit=(u.deg, u.deg),
                                                frame='galactic'),
                            radius=2.69**0.5*deg*u.deg, inner_radius=1.3*deg*u.deg,
                            catalog=["II/293/glimpse"], cache=False)

    
    result_2 = v.query_region(coord.SkyCoord(ra_gal, dec_gal,
                                                unit=(u.deg, u.deg),
                                                frame='galactic'),
                            radius=deg*u.deg,
                            catalog=["II/293/glimpse", "II/246/out"])
    
    
    data = result['II/293/glimpse']  #annular region wtih area same as cluster area
    data_1 = result_1['II/293/glimpse'] #field/background
    data_2 = result_2['II/293/glimpse'] #cluster region
    
    
    
    
    
    
    #load filtered data
    #-------------------
    lim1, lim2 = 0.6, 4.0
    idx_cl = np.where((data['_3.6mag']-data['_4.5mag']>lim1) & (data['_3.6mag']-data['_4.5mag']<=lim2))[0]
    idx_cl_2 = np.where((data_2['_3.6mag']-data_2['_4.5mag']>lim1) & (data_2['_3.6mag']-data_2['_4.5mag']<=lim2))[0]
    
    idx_cl1 = np.where((data_1['_3.6mag']-data_1['_4.5mag']>lim1) & (data_1['_3.6mag']-data_1['_4.5mag']<=lim2))[0]
    
    
    idx_cl2 = np.where((data['_3.6mag']-data['_4.5mag']>lim1) & (data['_3.6mag']-data['_4.5mag']<=lim2))[0]
    idx_cl_22 = np.where((data_2['_3.6mag']-data_2['_4.5mag']>lim1) & (data_2['_3.6mag']-data_2['_4.5mag']<=lim2))[0]

    
    fig=plt.figure(constrained_layout=True)
    
    fig.set_figheight(42)
    fig.set_figwidth(42.3)
    s = 100
    spec = gridspec.GridSpec(ncols=2, nrows=2,
                         width_ratios=[3.3,3.3], wspace=0.3,
                         hspace=0.3, height_ratios=[3,3])
    ax1 = fig.add_subplot(spec[0])
    ax2 = fig.add_subplot(spec[2])
    #ax3 = fig.add_subplot(spec[2])
    #ax4 = fig.add_subplot(spec[5])
    ax5 = fig.add_subplot(spec[1])
    ax5.grid()
    ax6 = fig.add_subplot(spec[3])
    ax6.grid()
    
    fig.subplots_adjust(hspace=0)
    
    ax1.scatter(data['_3.6mag']-data['_4.5mag'], 
                data['_3.6mag'],s=s, label='background annulus',c='black')
    ax1.scatter(data_2['_3.6mag']-data_2['_4.5mag'], 
                data_2['_3.6mag'],s=s, label='potential cluster members',c='r')
    ax1.set_xlabel('[3.6]-[4.5], mag',fontsize=75)
    ax1.tick_params(axis='x', labelsize=50)
    ax1.tick_params(axis='y', labelsize=50)
    ax1.xaxis.set_tick_params(width=5, size=25)
    ax1.yaxis.set_tick_params(width=5, size=25)
    ax1.axvline(lim1, c='blue', linestyle='--')
    ax1.axvline(lim2, c='blue', linestyle='--')
    #ax1.set_xticklabels(np.round(data_2['_3.6mag']-data_2['_4.5mag'],1), fontsize=35)
    #ax1.ylabel('[3.6]',fontsize=18)
    #ax1.gca().invert_yaxis()
    ax2.scatter(data['Kmag']-data['_3.6mag'], data['Kmag'],s=s,c='black')
    ax2.scatter(data_2['Kmag']-data_2['_3.6mag'], data_2['Kmag'],s=s, c='r')
    #
    ax2.scatter(data_2['Kmag'][idx_cl_22]-data_2['_3.6mag'][idx_cl_22], data_2['Kmag'][idx_cl_22],s=6*s, c='r', marker='x', label='sources with GLIMPSE counterparts')
    #ax2.scatter(data_2['Kmag'][idx_cl_22]-data_2['_3.6mag'][idx_cl_22], data_2['Kmag'][idx_cl2],s=6*s, c='black', marker='x')
    #
    ax1.get_shared_x_axes().join(ax1, ax2)
    ax1.set_ylabel('[3.6], mag',fontsize=75)
    ax2.set_ylabel('[K], mag',fontsize=75)
    ax2.tick_params(axis='x', labelsize=50)
    ax2.tick_params(axis='y', labelsize=50)
    ax2.xaxis.set_tick_params(width=5, size=25)
    ax2.yaxis.set_tick_params(width=5, size=25)
    #ax2.axvline(lim1, c='blue', linestyle='--')
    #ax2.axvline(lim2, c='blue', linestyle='--')
    ax2.set_xlabel('[K]-[3.6], mag',fontsize=75)
    ax1.invert_yaxis()
    ax2.invert_yaxis()
    ax1.set_xlim(-1, 4.5)
    ax2.set_xlim(-1, 4.5)
    ax1.legend(fontsize=50)
    ax2.legend(fontsize=40)
    
    def rms(x):
        return np.sqrt(x.dot(x)/x.size)
    try:
        annulus_log['Glon'].append(ra)
        annulus_log['Glat'].append(dec)
        annulus_log['annulus_color_med'].append(np.median(data['_3.6mag'])-np.median(data['_4.5mag']))
        annulus_log['annulus_color_rms'].append(rms(data['_3.6mag']-data['_4.5mag']))
        annulus_log['cluster_color_med'].append(np.median(data_2['_3.6mag'])-np.median(data_2['_4.5mag']))
        annulus_log['cluster_color_rms'].append(rms(data_2['_3.6mag']-data_2['_4.5mag']))
        annulus_log['c2a_ratio'].append(round(len(data_2['_3.6mag'])/len(data['_3.6mag']),3))
    except np.ma.core.MaskError:
        print('Skipped log')
        pass
    
    #plt.show()
    
    
    #GETTING A SPITZER AND 2MASS three color image for the sources using pyESASky, experimental, credits to Marcos Lopez, ESA/ESO SciOps workshop, commented out due to projection and website unreponsive

    '''
    url = 'https://sky.esa.int/esasky-tap/skyimage?target='+str(ra)+'%20'+str(dec)+'&fov=0.1&aspectratio=2.0&hips=Spitzer%20cold%20SEIP%20IRAC-1-3-4%20RGB%20bright&size=512&fmt=PNG&_=1652798602280'
    response = requests.get(url)
    #img = Image.open(BytesIO(response.content))
    img=Image.open(BytesIO(response.content))

    arr=np.array(img)

    url_2mass = 'https://sky.esa.int/esasky-tap/skyimage?target='+str(ra)+'%20'+str(dec)+'&fov=0.1&aspectratio=2.0&hips=2MASS%20color%20JHK&size=512&fmt=PNG&_=1652798602280'
    response_2mass = requests.get(url_2mass)
    #img = Image.open(BytesIO(response.content))
    img_2mass=Image.open(BytesIO(response_2mass.content))

    arr_2mass=np.array(img_2mass)





    avg_2mass = (arr_2mass[:,:,0]+arr_2mass[:,:,1]+arr_2mass[:,:,2])/3

    #plt.imshow(arr, extent=(ra_gal-0.05/2,ra_gal+0.05/2,dec_gal-0.05/2,dec_gal+0.05/2))
    #plt.xlim(ra_gal-cluster['size'][0]*3, ra_gal+cluster['size'][0]*3)
    #plt.ylim(dec_gal-cluster['size'][0]*3, dec_gal+cluster['size'][0]*3)


    
    avg = (arr[:,:,0]+arr[:,:,1]+arr[:,:,2])/3
    s=[]
    for j in range(len(data_2['Glat'])):
        if data_2['_3.6mag'][j]-data_2['_4.5mag'][j]>0.01:
            s.append((data_2['_3.6mag'][j]-data_2['_4.5mag'][j])*100)
        else:
            s.append(1)
    s=10

    bg=2

    im=ax3.imshow(arr, extent=(ra-0.1/2,ra+0.1/2,dec-0.1/2,dec+0.1/2), aspect=1/2)
    ax3.set_xlim(ra-cluster['size'][0]*bg, ra+cluster['size'][0]*bg)
    ax3.set_ylim(dec-cluster['size'][0]*bg, dec+cluster['size'][0]*bg)
    ax3.set_title('GLIMPSE', fontsize=18)
    drawObject = Circle((ra, dec), cluster['size'][i], fill=False, color='r')
    ax3.add_patch(drawObject)
    drawObject = Circle((ra, dec), 1.3*cluster['size'][i], fill=False, color='white')
    ax3.add_patch(drawObject)
    drawObject = Circle((ra, dec), 2.69**.5*cluster['size'][i], fill=False, color='white')
    ax3.add_patch(drawObject)
    ax3.scatter(clust_ra_2, clust_dec_2,s=s,edgecolors='r', facecolors='none')
    ax3.scatter(clust_ra, clust_dec,edgecolors='black', facecolors='none')
    #ax3.Circle((ra_gal, dec_gal), 3*cluster['size'][i], fill=False, c='black' )
    #ax3.Circle((ra_gal, dec_gal), 10**.5*cluster['size'][i], fill=False, c='black' )
    #plt.colorbar(im,ax=ax3)

    im_2mass=ax4.imshow(arr_2mass, extent=(ra-0.05/2,ra+0.05/2,dec-0.05/2,dec+0.05/2), aspect=1/2)
    ax4.set_xlim(ra-cluster['size'][0]*bg, ra+cluster['size'][0]*bg)
    ax4.set_ylim(dec-cluster['size'][0]*bg, dec+cluster['size'][0]*bg)
    ax4.set_title('2MASS', fontsize=18)
    drawObject = Circle((ra, dec), cluster['size'][i], fill=False, color='r')
    ax4.add_patch(drawObject)
    drawObject = Circle((ra, dec), 1.3*cluster['size'][i], fill=False, color='white')
    ax4.add_patch(drawObject)
    drawObject = Circle((ra, dec), 2.69**.5*cluster['size'][i], fill=False, color='white')
    ax4.add_patch(drawObject)
    ax4.scatter(clust_ra_2, clust_dec_2,s=s,edgecolors='r', facecolors='none')
    ax4.scatter(clust_ra, clust_dec,edgecolors='black', facecolors='none')
    #plt.colorbar(im_2mass,ax=ax4)
    #ax2.imshow(arr[:,:,2])
    '''
    
    s = 150
    ax5.scatter(data_1['_Glon'][idx_cl1], data_1['_Glat'][idx_cl1],c=data_1['_3.6mag'][idx_cl1]-data_1['_4.5mag'][idx_cl1],
                s=s, alpha=0.5, cmap='rainbow')
    sc1=ax5.scatter(data_2['_Glon'][idx_cl_2], data_2['_Glat'][idx_cl_2],s=s, facecolors='none', cmap='rainbow',
                    c=data_2['_3.6mag'][idx_cl_2]-data_2['_4.5mag'][idx_cl_2])
    #print(data_2['_3.6mag'][idx_cl_2]-data_2['_4.5mag'][idx_cl_2])
    ax5.scatter(data['_Glon'][idx_cl], data['_Glat'][idx_cl], facecolors='none', cmap='rainbow',
                c=data['_3.6mag'][idx_cl]-data['_4.5mag'][idx_cl], s=s)
    
    divider = make_axes_locatable(ax5)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar = fig.colorbar(sc1, cax=cax)
    cbar.ax.tick_params(labelsize=50) 
    cbar.ax.set_title('[3.6]-[4.5], mag',fontsize=50)
    
    drawObject = Circle((ra_gal, dec_gal), deg, fill=False, color='r')
    ax5.add_patch(drawObject)
    drawObject = Circle((ra_gal, dec_gal), 1.3*deg, fill=False, color='black')
    ax5.add_patch(drawObject)
    drawObject = Circle((ra_gal, dec_gal), 2.69**.5*deg, fill=False, color='black')
    ax5.add_patch(drawObject)
    ax5.tick_params(axis='x', labelsize=50)
    ax5.tick_params(axis='y', labelsize=50)
    ax5.xaxis.set_tick_params(width=5, size=25)
    ax5.yaxis.set_tick_params(width=5, size=25)
    ax5.set_xlabel('$l$', fontsize=75)
    ax5.set_ylabel('$b$', fontsize=75)
    ax5.set_xlim(ra_gal-2.69**0.5*deg*1.2,ra_gal+2.69**0.5*deg*1.2)
    ax5.set_ylim(dec_gal-2.69**0.5*deg*1.2,dec_gal+2.69**0.5*deg*1.2)
    
    
    ax6.scatter(data_1['_Glon'][idx_cl1], data_1['_Glat'][idx_cl1], edgecolors='black',
                s=0.1*s, alpha=0.5, c='black')
    ax6.scatter(data_1['_Glon'][idx_cl1], data_1['_Glat'][idx_cl1], edgecolors='black', c=data_1['Kmag'][idx_cl1]-data_1['_3.6mag'][idx_cl1],
                s=2*s, alpha=1, cmap='rainbow')
    sc2=ax6.scatter(data_2['_Glon'][idx_cl_22], data_2['_Glat'][idx_cl_22],s=2*s, edgecolors='black', 
                    c=data_2['Kmag'][idx_cl_22]-data_2['_3.6mag'][idx_cl_22], cmap='rainbow')
    ax6.scatter(data['_Glon'][idx_cl2], data['_Glat'][idx_cl2], 
                edgecolors='black', c=data['Kmag'][idx_cl2]-data['_3.6mag'][idx_cl2], cmap='rainbow', s=2*s)
    
    divider = make_axes_locatable(ax6)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar = fig.colorbar(sc2, cax=cax)
    cbar.ax.tick_params(labelsize=50) 
    cbar.ax.set_title('[K]-[3.6], mag',fontsize=50)
    drawObject = Circle((ra_gal, dec_gal), cluster['size'][i], fill=False, color='r')
    ax6.add_patch(drawObject)
    drawObject = Circle((ra_gal, dec_gal), 1.3*cluster['size'][i], fill=False, color='black')
    ax6.add_patch(drawObject)
    drawObject = Circle((ra_gal, dec_gal), 2.69**.5*cluster['size'][i], fill=False, color='black')
    ax6.add_patch(drawObject)
    ax6.tick_params(axis='x', labelsize=50)
    ax6.tick_params(axis='y', labelsize=50)
    ax6.xaxis.set_tick_params(width=5, size=25)
    ax6.yaxis.set_tick_params(width=5, size=25)
    ax6.set_xlabel('$l$', fontsize=75)
    ax6.set_ylabel('$b$', fontsize=75)
    ax6.set_xlim(ra_gal-2.69**0.5*deg*1.2,ra_gal+2.69**0.5*deg*1.2)
    ax6.set_ylim(dec_gal-2.69**0.5*deg*1.2,dec_gal+2.69**0.5*deg*1.2)
    
    
    df = pd.DataFrame({'_3.6mag': data['_3.6mag'], '_4.5mag': data['_4.5mag']})
    df = df.dropna()
    df_2 = pd.DataFrame({'_3.6mag': data_2['_3.6mag'], '_4.5mag': data_2['_4.5mag']})
    df_2 = df_2.dropna()
    df = df[(df['_3.6mag']-df['_4.5mag']<lim2) & (df['_3.6mag']-df['_4.5mag']>lim1)]
    df_2 = df_2[(df_2['_3.6mag']-df_2['_4.5mag']<lim2) & (df_2['_3.6mag']-df_2['_4.5mag']>lim1)]
    clust_count_r = df_2.shape[0]
    ann_count_r = df.shape[0]
    clust_od = (clust_count_r-ann_count_r)/np.sqrt(ann_count_r)
    clust_od = round(clust_od, 2)
    clust_stat['overdensity'] = clust_od
    fig.suptitle('Gal Coords: '+str(round(ra_gal,4))+', '+str(round(dec_gal,4))+' #cluster members: '
                 +str(clust_count_r)+'# nclust_optics: '+str(cluster['sources_cluster'][i])+' #field sources: '+str(ann_count_r)+
                 '\n cluster radius: '+str(round(cluster['size'][i]*3600,1))+"'' Overdensity: "+str(clust_od)+'$\sigma$',
                fontsize=50)
    
    
    print('saving '+str(ra_gal)+'..... ')
    os.chdir(path_op)
    plt.savefig(str(ra_gal)+'.jpeg')
    print(round(100*(i+1)/cluster.shape[0],3), '% done...')
    plt.close()

print('100% done')