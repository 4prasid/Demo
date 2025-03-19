#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import astropy
import os


# In[2]:


files = [r"0h_90.csv",
        r"5h_-15_0.csv",
        r"6h.csv",
        r"6h_30.csv",
        r"9h_-30_-45.csv",
        r"18h.csv",
        r"18h_0_30.csv",
        r"18h_0h_90.csv",
        r"-30.csv",
        r"30_0h.csv",
        r"-30_6h.csv",
        r"-30_12.5h.csv",
        r"30_12h.csv",
        r"-30_12h.csv",
        r"-30_18h.csv",
        r"-45_-75_12h.csv",
        r"-60.csv",
        r"60_6h.csv",
        r"-90_6h.csv",
        r"above_18h_-30.csv"]


# In[3]:


data = [pd.read_csv(i) for i in files ]


# In[4]:


parallax = []
parallax_error = []
ra1 = []
dec1 = []
g_mag = []
source_id = []
pmra = []
pmdec = []
pmra_error = []
pmdec_error  = []

for i in range(20):
    parallax.append(np.array(data[i]['parallax']))
    parallax_error.append(np.array(data[i]['parallax_error']))
    ra1.append( np.array(data[i]['ra']))
    dec1.append(np.array(data[i]['dec']))
    g_mag.append(np.array(data[i]['phot_g_mean_mag']))
    source_id.append(np.array(data[i]['source_id']))

    pmra.append(np.array(data[i]['pmra']))
    pmdec.append(np.array(data[i]['pmdec']))
    pmra_error.append(np.array(data[i]['pmra_error']))
    pmdec_error.append(np.array(data[i]['pmdec_error']))


# In[5]:


valid_pairs1 = []

# Iterate over all pairs of indices for coordinate
for k in range(len(files)):
    
    temp_ra1 = ra1[k]
    temp_dec1 = dec1[k]
    temp_g_mag = g_mag[k]
    temp_parallax = parallax[k]
    temp_source_id = source_id[k]
    temp_valid_pairs1 = []
    
    for i in range(len(temp_source_id)):
        for j in range(i + 1, len(temp_source_id)):  
            # Calculate the absolute difference i.e. angluar separtion
            abs_diff = (np.sqrt((temp_ra1[i]-temp_ra1[j])**2 + (temp_dec1[i]-temp_dec1[j])**2))*3600
        
            # Find the maximum of g_mag[i] and g_mag[j]
            max_g_mag = max(temp_g_mag[i], temp_g_mag[j])
        
            # Get the corresponding index in parallax for the maximum g_mean
            index = np.argmax(temp_g_mag == max_g_mag)  # This will return the first index of the maximum
        
            # Check the condition
            if abs_diff <= 206.265 * temp_parallax[index]:
                temp_valid_pairs1.append((temp_source_id[i], temp_source_id[j],  i, j, k))
    
    valid_pairs1.append(temp_valid_pairs1)


# In[ ]:


column_names = ['Source_1','Source_2','Index_of_Source_1','Index_of_Source_2','Culster']

# Create a directory to store CSV files (if it doesn't exist)
output_dir1 = "valid_pairs1"
os.makedirs(output_dir1, exist_ok=True)

# Iterate through each sublist and write it to a separate CSV file
for idx, sublist in enumerate(valid_pairs1, start=1):
    # Convert the sublist of tuples into a DataFrame
    df1 = pd.DataFrame(sublist, columns=column_names)
    
    # Define the filename for each CSV
    file_name = f"sublist_{idx}.csv"
    file_path1 = os.path.join(output_dir1, file_name)
    
    # Write the DataFrame to CSV
    df1.to_csv(file_path1, index=False)


# In[ ]:


with open('No_of_Pairs.txt', 'a') as file:
    file.write('\n'+'No. of valid pairs by condition 1:')
for k in range(len(files)):
    with open('No_of_Pairs.txt', 'a') as file:
        file.write(f'\nIn cluster {k}: {len(valid_pairs1[k])}')


# In[ ]:


# Extracting the index of valid pairs by conditon 1:

index11 = []
index12 = []
for k in range(len(files)):
    temp_valid_pairs1 = valid_pairs1[k]
    index11.append([temp_valid_pairs1[i][2] for i in range(len(temp_valid_pairs1))])
    index12.append([temp_valid_pairs1[i][3] for i in range(len(temp_valid_pairs1))])


# In[ ]:


valid_pairs12 = []

# Iterate over all pairs of indices for coordinate
# for k in index11:


for k in range(len(files)):
    
    temp_valid_pairs1 = valid_pairs1[k]
    temp_ra1 = ra1[k]
    temp_dec1 = dec1[k]
    temp_parallax = parallax[k]
    temp_parallax_error = parallax_error[k]
    temp_source_id = source_id[k]
    temp_valid_pairs12 = []
    
    
    for m in range(len(temp_valid_pairs1)):  
        i = index11[k][m]
        j = index12[k][m]
        # Calculate the absolute difference i.e. angluar separtion
        abs_diff = (np.sqrt((temp_ra1[i]-temp_ra1[j])**2 + (temp_dec1[i]-temp_dec1[j])**2))*3600
        
        # Check the conditions
        # for theta < 4, b = 6
        if abs_diff < 4:
            abs(temp_parallax[i]-temp_parallax[j]) < 6*(np.sqrt((temp_parallax_error[i])**2 + (temp_parallax_error[j])**2))
            temp_valid_pairs12.append((temp_source_id[i], temp_source_id[j],  i, j, k+1))
        # for theta > 4, b = 3
        if abs_diff > 4:
            abs(temp_parallax[i]-temp_parallax[j]) < 3*(np.sqrt((temp_parallax_error[i])**2 + (temp_parallax_error[j])**2))
            temp_valid_pairs12.append((temp_source_id[i], temp_source_id[j],  i, j, k+1))
            
    valid_pairs12.append(temp_valid_pairs12)


# In[ ]:


# Create a directory to store CSV files (if it doesn't exist)
output_dir12 = "valid_pairs12"
os.makedirs(output_dir12, exist_ok=True)

# Iterate through each sublist and write it to a separate CSV file
for idx, sublist in enumerate(valid_pairs12, start=1):
    # Convert the sublist of tuples into a DataFrame
    df12 = pd.DataFrame(sublist, columns=column_names)
    
    # Define the filename for each CSV
    file_name = f"sublist_{idx}.csv"
    file_path12 = os.path.join(output_dir12, file_name)
    
    # Write the DataFrame to CSV
    df12.to_csv(file_path12, index=False)


# In[ ]:


with open('No_of_Pairs.txt', 'a') as file:
    file.write('\n'+'\n'+'no. of valid pairs by applying condition 2 on valid pairs of condition 1:')
for k in range(len(files)):
    with open('No_of_Pairs.txt', 'a') as file:
        file.write(f'\nIn cluster {k}: {len(valid_pairs12[k])}')


# In[ ]:


# Extracting the index of valid pairs by conditon 2:

index21 = []
index22 = []
for k in range(len(files)):
    temp_valid_pairs12 = valid_pairs12[k]
    index21.append([temp_valid_pairs12[i][2] for i in range(len(temp_valid_pairs12))])
    index22.append([temp_valid_pairs12[i][3] for i in range(len(temp_valid_pairs12))])


# In[ ]:


valid_pairs23 = []
 
# Iterate over all pairs of indices for coordinate
# for k in index11:
for k in range(len(files)):
    
    temp_valid_pairs12 = valid_pairs12[k]
    temp_ra1 = ra1[k]
    temp_dec1 = dec1[k]
    temp_parallax = parallax[k]
    temp_parallax_error = parallax_error[k]
    temp_g_mag = g_mag[k]
    temp_pmra = pmra[k]
    temp_pmdec = pmdec[k]
    temp_pmra_error = pmra_error[k]
    temp_pmdec_error = pmdec_error[k]
    temp_source_id = source_id[k]
    temp_valid_pairs23 = []
    
    
    for m in range(len(temp_valid_pairs12)):  
        i = index21[k][m]
        j = index22[k][m]
        # Calculate the absolute difference i.e. angluar separtion
        abs_diff = (np.sqrt((temp_ra1[i]-temp_ra1[j])**2 + (temp_dec1[i]-temp_dec1[j])**2))*3600
        # Calculate delta_mu
        delta_mu = np.sqrt((temp_pmra[i]*np.cos(temp_dec1[i]) - temp_pmra[j]*np.cos(temp_dec1[j]))**2 + 
                           (temp_pmdec[i]-temp_pmdec[j])**2)
        # Find the maximum of g_mag[i] and g_mag[j]
        max_g_mag = max(temp_g_mag[i], temp_g_mag[j])
        
        # Get the corresponding index in parallax for the maximum g_mean
        index = np.argmax(temp_g_mag == max_g_mag)  # This will return the first index of the maximum
        # Calculate delta_mu orbit
        delta_mu_orbit = 0.44*((temp_parallax[index])**1.5 )* (abs_diff**(-0.5))
        #Calculate sigma_delta_mu
        sigma_delta_mu = (np.sqrt((np.cos(temp_dec1[i])**2 * temp_pmra_error[i]**2 + 
                                   np.sin(temp_dec1[i])**2 * temp_pmra[i]**2 * temp_pmdec_error[i]**2 + 
                                   np.cos(temp_dec1[j])**2 * temp_pmra_error[j]**2 + 
                                   np.sin(temp_dec1[j])**2 * temp_pmra[j]**2 * temp_pmdec_error[j]**2) *
                                   (temp_pmra[i]*np.cos(temp_dec1[i]) - temp_pmra[j]*np.cos(temp_dec1[j]))**2 + 
                                   (temp_pmdec_error[i]**2 + temp_pmdec_error[j]**2)*(temp_pmdec[i]-temp_pmdec[j])**2 )
                         )/ delta_mu
        # Applying the condition
        if delta_mu <= delta_mu_orbit + 2*sigma_delta_mu:
            temp_valid_pairs23.append((temp_source_id[i], temp_source_id[j],  i, j, k+1))    
    
    valid_pairs23.append(temp_valid_pairs23)


# In[ ]:


# Create a directory to store CSV files (if it doesn't exist)
output_dir23 = "valid_pairs23"
os.makedirs(output_dir23, exist_ok=True)

# Iterate through each sublist and write it to a separate CSV file
for idx, sublist in enumerate(valid_pairs23, start=1):
    # Convert the sublist of tuples into a DataFrame
    df23 = pd.DataFrame(sublist, columns=column_names)
    
    # Define the filename for each CSV
    file_name = f"sublist_{idx}.csv"
    file_path23 = os.path.join(output_dir23, file_name)
    
    # Write the DataFrame to CSV
    df23.to_csv(file_path23, index=False)


# In[ ]:


with open('No_of_Pairs.txt', 'a') as file:
    file.write('\n'+'\n'+'no. of valid pairs by applying condition 3 on valid pairs of condition 2:')
for k in range(len(files)):
    with open('No_of_Pairs.txt', 'a') as file:
        file.write(f'\nIn cluster {k}: {len(valid_pairs23[k])}')

