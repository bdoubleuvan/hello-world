# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 14:12:34 2018

@author: presentatie
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#import matplotlib.pyplot as plt
import numpy as np
import piping_data_deltares as data
import piping_data_student as data1
import pandas as pd

#%% Deltares data
piping = data.dijk(3003005, 302)
Z_u = piping.Z_u(piping.h_exit, piping.D_cover, piping.m_u)
Z_h = piping.Z_h(piping.h_exit, piping.D_cover, piping.i_ch)
Z_p = piping.Z_p(piping.h_exit, piping.D_cover, piping.L, piping.D, piping.m_p , piping.k, piping.d70)
Z = np.logical_and(np.logical_and(Z_u < 0, Z_h < 0), Z_p < 0)

#%% Master student data

d = {}

for i in range(1, 15):
    str = 'DP-M {}'.format(i)
    piping1 = data1.dijk(str)
    
    Z_u1 = piping1.Z_u(piping1.h_exit, piping1.D_cover, piping1.m_u)
    Z_h1 = piping1.Z_h(piping1.h_exit, piping1.D_cover, piping1.i_ch)
    Z_p1 = piping1.Z_p(piping1.h_exit, piping1.D_cover, piping1.L, piping1.D, piping1.m_p , piping1.k, piping1.d70)
    Z1 = np.logical_and(np.logical_and(Z_u < 0, Z_h < 0), Z_p < 0)
    d[str] = Z_h1, Z_p1, Z_u1, Z1

df = pd.DataFrame(data = d, index = ['Z_u','Z_h','Z_p','Failure'])

#%% Deltares data
dijk_vak = data.dijk_vak(3003005)

#%%

#plt.figure(0)
#plt.hold(False)
#plt.plot(piping.h,Z_u<0,'o')
#plt.ylim(0,1)
#plt.xlabel('niveau buitenwaterstand')
#plt.ylabel('grenstoestandsfunctie')
#plt.title('Grenstoestandsfunctie uplift t.o.v. niveau buitenwaterstand')
#plt.show
#
#plt.figure(1)
#plt.plot(piping.h,Z_h<0,'o')
#plt.ylim(0,1)
#plt.xlabel('niveau buitenwaterstand')
#plt.ylabel('grenstoestandsfunctie')
#plt.title('Grenstoestandsfunctie heave t.o.v. niveau buitenwaterstand')
#plt.show
#
#plt.figure(2)
#plt.plot(piping.h,Z_p<0,'o')
#plt.ylim(0,1)
#plt.xlabel('niveau buitenwaterstand')
#plt.ylabel('grenstoestandsfunctie')
#plt.title('Grenstoestandsfunctie piping t.o.v. niveau buitenwaterstand')
#plt.show
#
#plt.figure(3)
#plt.plot(piping.h,Z,'o')
#plt.ylim(0,1)
#plt.xlabel('niveau buitenwaterstand')
#plt.ylabel('grenstoestandsfunctie')
#plt.title('Grenstoestandsfunctie piping t.o.v. niveau buitenwaterstand')
#plt.show