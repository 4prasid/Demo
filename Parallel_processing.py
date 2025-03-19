#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np


# In[2]:


from astropy.table import Table
import multiprocessing, psutil
from sklearn.neighbors import BallTree


# In[3]:


path = r"3_million_stars.fits"
tab = Table.read(path) # 3 million sources 


# In[4]:


def fetch_table_element(colname, table):
    '''
    avoid table['col'].data vs table['col'].data.data problems with masked arrays in astropy tables 
    '''
    if type(colname) == str:
        if type(table[colname].data.data) == memoryview:
            dat_ = table[colname].data
        else:
            dat_ = table[colname].data.data
    elif type(colname) == list:
        dat_ = []
        for col in colname:
            dat_.append(fetch_table_element(col, table))
    return dat_


# In[5]:


ra, dec, pmra, pmdec, parallax, parallax_error, pmra_error, pmdec_error, G = fetch_table_element(
    ['ra', 'dec', 'pmra', 'pmdec', 'parallax', 'parallax_error', 'pmra_error', 'pmdec_error', 
     'phot_g_mean_mag'], tab )


# In[6]:


size_max_pc = 5 # max projected separation out to which to search
dispersion_max_kms = 5 # max velocity difference in kms
s_max_cluster = 206265*size_max_pc
theta_max_radians = s_max_cluster/(1000/parallax)/3600 * np.pi/180
coords = np.vstack([dec*np.pi/180, ra*np.pi/180,]).T


# In[7]:


tree = BallTree(coords[G < 18], leaf_size = 10, metric = 'haversine') 
# build tree of all stars brighter than 18
# as we have used here the haversine metric so we no need to define the function get_distance_arcsec 
# haversine automatically does that job.
#def get_distance_arcsec(ra1, dec1, ra2, dec2):
#     '''
#     angular separations. coords are assumed to be in degrees
#     '''
#     ra_rad1, dec_rad1 = ra1*np.pi/180, dec1 * np.pi/180
#     ra_rad2, dec_rad2 = ra2*np.pi/180, dec2 * np.pi/180
#     d_ra, d_dec = ra_rad1 - ra_rad2, dec_rad1 - dec_rad2
    
#     d_theta = 2*np.arcsin(np.sqrt(np.sin(0.5*d_dec)**2 + np.cos(dec_rad1)*np.cos(dec_rad2)*np.sin(0.5*d_ra)**2))
#     d_theta_deg = 180/np.pi*d_theta
#     d_theta_arcsec = d_theta_deg * 3600
#     return d_theta_arcsec


# In[8]:


# data for stars brighter than G = 18
ra_b, dec_b, pmra_b, pmdec_b, parallax_b, parallax_error_b, pmra_error_b, pmdec_error_b, G_b = ra[G < 18], dec[G < 18], pmra[G < 18], pmdec[G < 18], parallax[G < 18], parallax_error[G < 18], pmra_error[G < 18], pmdec_error[G < 18], G[G < 18]


# In[9]:


Nblock = 20000 # how many stars to process at once per core
Nmax = len(coords)//Nblock 
sigma_cut = 2 # how many sigma tolerance 


# In[10]:


def get_delta_mu_and_sigma(pmra1, pmdec1, pmra2, pmdec2, pmra_error1, 
    pmdec_error1, pmra_error2, pmdec_error2):
    '''
    Uses standard uncertainty propagation 
    Equations 4-5 of the paper. 

    assume that "1" is a float and "2" is an array
    '''
    delt_alpha, delt_delta = (pmra1 - pmra2)**2, (pmdec1 - pmdec2)**2
    delta_mu2 = delt_alpha + delt_delta
    
    try:
        lenn = len(pmra2) # checks whether pmra2 is an array 
        m = delta_mu2 == 0
        sigma2_delta_mu = np.zeros(len(pmra2))
        if np.sum(m):
            sigma2_delta_mu[m] = (pmra_error1**2 + pmra_error2[m]**2) + (pmdec_error1**2 + pmdec_error2[m]**2)
        if np.sum(~m):
            sigma2_delta_mu[~m] = ((pmra_error1**2 + pmra_error2[~m]**2) * (delt_alpha[~m]) +                 (pmdec_error1**2 + pmdec_error2[~m]**2)*delt_delta[~m])/delta_mu2[~m]
    except: # pmra2 is a float 
        if delta_mu2 == 0:
            sigma2_delta_mu = pmra_error1**2 + pmra_error2**2 + pmdec_error1**2 + pmdec_error2**2
        else:
            sigma2_delta_mu = ((pmra_error1**2 + pmra_error2**2) * (delt_alpha) +                 (pmdec_error1**2 + pmdec_error2**2)*delt_delta)/delta_mu2

    return np.sqrt(delta_mu2), np.sqrt(sigma2_delta_mu)


# In[11]:


def query_this_j(j):
    '''
    function to pass to multiprocessing pool. deal with Nblock stars. 
    '''
    # see how far along we are and make sure we aren't running out of memory.
    #print(j, j*Nblock/len(coords),  psutil.virtual_memory().percent)

    # find the stars in this block
    msk = (np.arange(len(coords)) >= int(j*Nblock)) & (np.arange(len(coords)) < int((j+1)*Nblock))
    
    # find their companions and angular distances
    these_inds, these_dists = tree.query_radius(coords[msk], r = theta_max_radians[msk], return_distance = True)      
    
    # copy astrometry of stars in this block  
    parallax_, parallax_error_, pmra_, pmra_error_, pmdec_, pmdec_error_ = parallax[msk], parallax_error[msk],         pmra[msk], pmra_error[msk], pmdec[msk], pmdec_error[msk] 

    # for each star, see how many of the companions within 5 pc (projected) have consistent parallax and similar        proper motion 
    N_neighbors = np.zeros(len(parallax_))
    for i, idxs in enumerate(these_inds):
        thetas_arcsec = these_dists[i]*180/np.pi*3600
        d_par_over_sigma = np.abs(parallax_[i] - parallax_b[idxs])/np.sqrt(parallax_error_[i]**2 +                                            parallax_error_b[idxs]**2)
        delta_mu, sigma_delta_mu = get_delta_mu_and_sigma(pmra1 = pmra_[i], pmdec1 = pmdec_[i], 
            pmra2 = pmra_b[idxs], pmdec2 = pmdec_b[idxs], pmra_error1 = pmra_error_[i], 
            pmdec_error1 = pmdec_error_[i], pmra_error2 = pmra_error_b[idxs], 
            pmdec_error2 = pmdec_error_b[idxs])
            
        mu_max = 0.21095*dispersion_max_kms*parallax_[i] # i.e. 1.05 * parallax           
        neighbors = (delta_mu < mu_max + sigma_cut*sigma_delta_mu) & (d_par_over_sigma < sigma_cut) &                                  (thetas_arcsec > 1e-3) # theta > 1e-3 arcsec to make sure you don't count yourself as a                                                                                                                 neighbor
        N_neighbors[i] = np.sum(neighbors) 
    return N_neighbors


# In[ ]:


pool = multiprocessing.Pool(multiprocessing.cpu_count())
all_result = pool.map(query_this_j,  np.arange(Nmax))
pool.close()
N_neighbors = np.concatenate(all_result)


# In[ ]:


# save these for later
np.savez('neighbor_counts_edr3_all.npz', source_id = fetch_table_element('source_id', tab), N_neighbors = N_neighbors)

