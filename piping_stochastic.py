# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 17:48:39 2018

@author: presentatie
"""

#from scipy import io
#import numpy as np
#import matplotlib.pyplot as plt
#import openturns_piping_data as p_data
import piping_data_deltares as data
#import piping_stochastic_example_functions as p
#import scipy.stats as sp

#piping_data = io.loadmat('D:/Users/Nguyen/Desktop/Piping/Matlab/pipingmatrix.mat')
#
##1st row of data
#data = piping_data["pipingMatrix"][1,:]
#
##initialize mu and sigma
#mu = np.zeros(len(data[9:19]))
#sigma = np.zeros(len(data[9:19]))
#
##computes mu for lognormal distribution
#mu = data[9:19]
#
##computes sigma for lognormal distribution
#sigma = data[19:29]
#
##computes mu for lognormal distribution
#mu[0] = data[9] #h_exit
#mu[1:len(mu)] = np.log(data[10:19]/np.sqrt(1+data[20:29]/data[10:19]**2))
#mu[4]=data[13] #gamma_sat
#mu[6]=data[15] #theta
#mu[7]=data[16] #eta
#
##computes sigma for lognormal distribution
#sigma[0] = data[19]
#sigma[1:len(sigma)] = np.sqrt(np.log(1+data[20:29]/(data[10:19]**2)))
#sigma[4]=data[23]
#sigma[6]=data[25]
#sigma[7]=data[26]
#
##aantal iteraties
#n = 100000
#
##niveau buitenwaterstand 
#h = data[29]
#
##perform Monte Carlo and calculates mean and standard deviation
#uhp_data = p.MC_uhp(h,mu,sigma,n)
#
#uplift_mean = (1/n)*np.sum(uhp_data[:,0]) 
#heave_mean = (1/n)*np.sum(uhp_data[:,1]) 
#piping_mean = (1/n)*np.sum(uhp_data[:,2])
# 
#uplift_std = np.sqrt((1/(n-1))*np.sum((uhp_data[:,0]-uplift_mean)**2)) 
#heave_std = np.sqrt((1/(n-1))*np.sum((uhp_data[:,1]-heave_mean)**2))
#piping_std = np.sqrt((1/(n-1))*np.sum((uhp_data[:,2]-piping_mean)**2))
#
#uplift_fail = (1/n)*np.sum(uhp_data[:,0]<0) 
#heave_fail = (1/n)*np.sum(uhp_data[:,1]<0) 
#piping_fail = (1/n)*np.sum(uhp_data[:,2]<0) 
#
#print('The approximated mean for uplift: %f \n' %uplift_mean)
#print('The approximated mean for heave: %f \n' %heave_mean)
#print('The approximated mean for piping: %f \n' %piping_mean)
#print('The approximated standard deviation for uplift: %f \n' %uplift_std)
#print('The approximated standard deviation for heave: %f \n' %heave_std)
#print('The approximated deviation for piping: %f \n' %piping_std)
#
#print('The failure probability for uplift: %f \n' %uplift_fail)
#print('The failure probability for heave: %f \n' %heave_fail)
#print('The failure probability mean for piping: %f \n' %piping_fail)
#
#piping_fail_req = (1/n)*np.sum(np.sum(uhp_data<0,axis=1)==3)
#piping_fail_req_inverse = sp.norm.ppf(1-piping_fail_req)
#
#print('The failure probability mean for piping with requirements: %f \n' %piping_fail_req)
#print('Required reliability index (a = 0.75 and b = 300): %f' % piping_fail_req_inverse)
#
#plt.figure(0)
#x=np.linspace(-5,5,1000)
#plt.plot(x,sp.norm.pdf(x)) 
#plt.axvline(x=piping_fail_req_inverse) 
#x_new = np.linspace(-5,piping_fail_req_inverse,1000)
#plt.fill_between(x_new,0,sp.norm.pdf(x_new),color='grey')
#plt.show()

#%%

#number of iterations
dijk = data.dijk(3003005,302)
n = 1000000
print('Pf = ', dijk.mc(n))