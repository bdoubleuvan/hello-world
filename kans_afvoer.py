# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 10:41:18 2018

@author: nguyen
"""

import matplotlib.pyplot as plt
import numpy as np
import openturns as ot
import pandas as pd
import hrd
from scipy.interpolate import interp1d

import os
import glob
import urllib

a = glob.glob(r'C:\MyPrograms\Hydra-NL\werkmap\WBI2017_Benedenrijn_15-1_v04\*\Berekeningen\*\*hfreq.txt')
b = glob.glob(r'C:\MyPrograms\Hydra-NL\werkmap\WBI2017_Benedenrijn_15-1_v04\*\Berekeningen\*\*uitvoer.html')

#Hydra-NL

plt.close('all')

class water_probability(ot.PythonDistribution):
    def __init__(self, q, exc_prob):
        self.q = np.linspace(q.min(), q.max(), 1000)
        #np.r_[q, np.linspace(q[-1],q[-1]+(q[-1]-q[0])/3,len(q)//3)[1:]]
        self.exc_prob = exc_prob#We obtain P(H>h), hence this transformation is necessary
        self.interp = np.exp(interp1d(q.squeeze(), np.log(self.exc_prob.squeeze()), fill_value = 'extrapolate')(self.q))
#        self.interp1 = interp1d(q.squeeze(), self.exc_prob.squeeze(), fill_value = 'extrapolate')(self.q)
        super().__init__()
    
#path = r'D:\Users\Nguyen\Desktop\piping_data\hfreq_15001003.txt'
#data = pd.read_csv(path, header = 0, skiprows = 0, sep = '\s+')     
#h = np.expand_dims(data.index.values, axis = 1)
#cdf = data.values
#cdf_dis = water_probability(h, cdf)

path1 = r'C:\MyPrograms\Hydra-NL\data\invoer\Afvoer\Borgharen\Ovkans_Borgharen_piekafvoer_2017_metOnzHeid.txt' 
data1 = pd.read_csv(path1, comment = '*', header = None, sep = '\s+')       
q = data1.values[:,0]
q_prob = data1.values[:,1]
q_cdf = water_probability(q, q_prob)

plt.figure(0)
plt.plot(q, q_prob, marker='.', label='origineel')
plt.plot(q_cdf.q, q_cdf.interp, label='Geinterpoleerd')
#plt.plot(q_cdf.q, q_cdf.interp1, '--')
plt.xlabel('Debiet')
plt.ylabel('Kans')
plt.title('Kansverdelingsfunctie van debiet')

#plt.yscale('log')
plt.legend()
#plt.figure(1)
#plt.plot(q_cdf.q, np.log(q_cdf.interp), marker='.')
##plt.plot(q_cdf.q, np.log(q_cdf.interp1), '--', marker='.')
#plt.xlabel('Debiet')
#plt.ylabel('log(Kans)')
#plt.title('Kansverdelingsfunctie van debiet')
#
#plt.figure(2)
#plt.plot(q_cdf.q, np.log(1/q_cdf.interp), marker='.')
##plt.plot(q_cdf.q, np.log(1/q_cdf.interp1), '--')
#plt.xlabel('Debiet')
#plt.ylabel('log(Kans)')
#plt.title('Kansverdelingsfunctie van debiet')
#
path2 = r'D:\Users\Nguyen\Downloads\2018-08-14\18_bovenmaas_selectie.sqlite'
table = hrd.export_tables([2433], path2)
result = table[1][['Discharge Borgharen','h']].sort_values(by = ['h']).values
discharge = result[:,0]
h = result[:,1]

table_test = hrd.export_tables((0,0), path2)

#h_interp = interp1d(h, discharge, kind = "linear", fill_value = 'extrapolate')
#water_level = 44
#print('Probability with water level {}m = '.format(water_level), q_cdf.probability(h_interp(water_level)))        