#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 14:08:53 2022

@author: Xiaoyu Tong
email: xit321@lehigh.edu
"""

import scipy.io
import contrastive
import numpy as np
from scipy.io import savemat
from sklearn.decomposition import PCA

EEGband = 'ALPHA' # frequency band
dataName = '/home/xit321/yuzi20_123121/xit321/' + 'compiledData_Russ31_CPCA_' + EEGband + '.mat' # EEG FC data matrix
content = scipy.io.loadmat(dataName)
Fe_HC = content['Fe_HC']
Fe_ASD = content['Fe_ASD'] # customize your retrieval of patients and TD group. 

nPC = 300 # edit as needed

alpha_values = [0.5,1] # edit as needed
nAlpha = len(alpha_values)

foreground_train, foreground_test = Fe_ASD, Fe_ASD # use all ASD patients for learning the transformation
# assign data from different sets of ASD patients to foreground_train & foreground_test if you are doing
# this in a cross-validation manner or have independent test data
background_data = Fe_HC

Fe_testArchive = np.zeros([foreground_test.shape[0],nPC,len(alpha_values)])

for i in range(len(alpha_values)): # derive contrastive connectivity features
    alpha = alpha_values[i]
    mdl_CPCA = contrastive.CPCA(n_components=nPC)
    mdl_CPCA.fit(foreground_train, background_data) # you may run into issues if your dimensionality > 1000. 
    # either refer to the original GitHub repository of cPCA or I may be able to offer help
    Fe_train = mdl_CPCA.transform(foreground_train,alpha_selection='manual',alpha_value = alpha)
    Fe_test = mdl_CPCA.transform(foreground_test,alpha_selection='manual',alpha_value = alpha)
    Fe_testArchive[:,:,i] = Fe_test
    

# save contrastive connectivity features (edit as needed)
# 
# filename = 'CPC' + str(nPC) + '_allPatients_'+EEGband 
# filename.replace('.', '')
# mdic = {'Fe_testSpecific':Fe_testArchive}
# savemat(filename+'.mat',mdic)

print('Successfully extracted contrastive connectivity features')