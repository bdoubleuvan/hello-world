# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 16:46:56 2018

@author: nguyen
"""

from __future__ import print_function
import pandas as pd
import openturns as ot
import math as m
import tqdm as tqdm
import numpy as np

class dijk:
         
    #pad met data
    pad = r'D:\Users\Nguyen\Desktop\piping_data\Piping berekeningen traject 15-1.xlsm'
    
    #Ophalen gegevens van excelsheet: Algemene pipinggegevens
    data = pd.read_excel(pad, sheet_name = ['Algemene pipinggegevens'], header = None, skiprows = 1, index_col = 0)#, index_col = 0, skipfooter = 3)
    data = data['Algemene pipinggegevens'] #Van een OrderDict een element pakken
    data = data.pivot_table(index = 0) #removes NaN values with index = 0
    data.columns = ['Waarden'] #headers aanpassen
    
    #Ophalen gegegevens excelsheet: Pipingsheet
    data1 = pd.read_excel(pad, sheet_name = ['Pipingsheet'], header = [1], skiprows = [0,1,2,5,6], index_col=[0,1])
    data1 = data1['Pipingsheet']
    
    #Sellmeijer model deterministische waarden
    d70m = data['Waarden']['D70m [m]'] #70% quantile of grain size
    eta = data['Waarden']['white coefficient (eta)'] #White's drag coefficient
    g = data['Waarden']['g [m/s2]'] #gravitational constant
    gamma_sp = data['Waarden']['eff korrelgew (gamma_p) [kN/m3]'] #volumetric weight of submerged sand particle
    gamma_water = 10 #volumetric weight of water
    r_c = 0.3 #reduction factor
    r_exit = 1 #damping factor at exit
    theta = data['Waarden']['rolweerstandshoek zandkorrels (theta)'] #bedding angle of sand grains
    v_water = data['Waarden']['kin_visc'] #kinematic viscosity of water
    
    #Sellmeijer model stochastische waarden
    i_ch_std = 0.1 #critical heave gradient sigma
    m_u = 1 #mean model factor for uplift
    m_u_std = 0.1 #standard deviation model factor for uplift
    m_p = 1 #mean model factor for piping
    m_p_std = 0.12 #standard deviation model factor for piping
    
    #Constructor, ophalen van data over dijk uit excelsheet
    def __init__(self, dijkpaal):
        self.dijkpaal = dijkpaal
        self.data_dijkpaal = self.data1[dijkpaal]
        
        #Sellmeijer model deterministische waarden
        self.h = self.data_dijkpaal['Algemene gegevens', 'Waterstand bij norm [m NAP]']
        self.gamma_sat = self.data_dijkpaal['Berekening Opbarsten', 'Volumegewicht deklaag verzadigd [kN/m3]'] #volumetric weight of saturated sand
        
        #Sellmeijer model stochastische waarden
        self.i_ch = self.data_dijkpaal['Berekening Heave', 'Kritieke heave gradient'] #critical heave gradient mu        
        self.D = self.data_dijkpaal['Piping gegevens', 'Dikte zandlaag [m]'] #upper sand cover layer
        self.D_std = 0
        self.D_cover = self.data_dijkpaal['Piping gegevens', 'Dikte deklaag [m]'] #aquifer cover layer
        self.D_cover_std = 0
        self.d70 = self.data_dijkpaal['Piping gegevens', 'Korreldiameter [m]'] #grain size
        self.d70_std = 0
        self. h_exit = self.data_dijkpaal['Piping gegevens', 'Kwelslootpeil [m NAP]'] #seepage ditch level
        self.h_exit_std = 0
        self.k = self.data_dijkpaal['Piping gegevens', 'Doorlatendheid zandlaag [m/s]'] #permeability sand
        self.k_std = 0
        self.L = self.data_dijkpaal['Piping gegevens', 'Aanwezige kwelweglengte [m]'] #seepage length
        self.L_std = 0
        
        self.Z_function = ot.PythonFunction(9, 1, self.Z)

#%%
        
    def print_data(self):
        print(self.data)
        print(self.data1[self.dijkpaal])
        
    def transform(self, mu, sigma):
        mu_t = m.log(mu / m.sqrt(1 + sigma / mu**2))
        sigma_t = m.sqrt(m.log(1 + sigma / mu**2))
        return [mu_t, sigma_t]
    
    def delta_phi_cu(self, D_cover):
        result = D_cover * (self.gamma_sat - self.gamma_water) / self.gamma_water
        #print('delta_phi_cu = ', result)
        return result
    
    def phi_exit(self, h_exit):
        result = h_exit + self.r_exit * (self.h - h_exit)
        #print('phi_exit = ', result)
        return result
    
    def F_res(self):
        result = self.eta * (self.gamma_sp / self.gamma_water) * m.tan(self.theta * m.pi / 180)
        #print('F_res = ', result)
        return result
    
    def F_scale(self,k, d70, L):
        kappa = (self.v_water / self.g) * k
        result = self.d70m / (kappa * L)**(1/3) * (d70 / self.d70m)**0.4
        #print('F_scale = ', result)
        return result
    
    def F_geo(self, D, L):
        if (D != L):
            result = 0.91 * (D / L)**(0.28 / ((D / L)**2.8 - 1) + 0.04)
        else:
            result = 1
        #print('F_geo = ', result)
        return result
    
    def H_c(self, k, d70, L, D):
        result = L * self.F_res() * self.F_scale(k, d70, L) * self.F_geo(D, L)
        #print('Kritiek verval = ', result)
        return result
        
    def Z_u(self, h_exit, D_cover, m_u):
        Z = m_u * self.delta_phi_cu(D_cover) - (self.phi_exit(h_exit) - h_exit)
        return Z
    
    def Z_h(self, h_exit, D_cover, i_ch):
        i = (self.phi_exit(h_exit) - h_exit) / D_cover
        #print('i = ', i)
        Z = i_ch - i
        return Z
    
    def Z_p(self, h_exit, D_cover, L, D, m_p , k, d70):
        Z = m_p * self.H_c(k, d70, L, D) - (self.h - h_exit - self.r_c * D_cover) 
        return Z
    
    def Z(self, X):
        h_exit, D_cover, m_u, i_ch, L, D, m_p, k, d70 = X
        Z_u_value = self.Z_u(h_exit, D_cover, m_u)
        Z_h_value = self.Z_h(h_exit, D_cover, i_ch)
        Z_p_value = self.Z_p(h_exit, D_cover, L, D, m_p , k, d70)
        result = int((Z_u_value < 0) and (Z_h_value < 0) and (Z_p_value < 0))
        print(result)
        return [result]
    
#    def Z_function(self):
#        return ot.PythonFunction(9, 1, self.Z)
    
    def mc(self,n):
        lognormdist = np.random.lognormal    
        result = np.zeros(3)
        result_p = 0
        [D_cover, D_cover_s] = self.transform(self.D_cover, self.D_cover_std)
        [L, L_s] = self.transform(self.L, self.L_std)
        [D, D_s] = self.transform(self.D, self.D_std)
        [m_u, m_u_s] = self.transform(self.m_u, self.m_u_std)
        [m_p, m_p_s] = self.transform(self.m_p, self.m_p_std)
        [i_ch, i_ch_s] = self.transform(self.i_ch, self.i_ch_std)
        for i in tqdm.tqdm(range(n)):
            h_exit_t = np.random.normal(self.h_exit, self.h_exit_std)
            #phi_polder = ot.Normal(self.phi_polder,self.phi_polder_s)
            D_cover_t = lognormdist(D_cover, D_cover_s)
            L_t = lognormdist(L, L_s)
            D_t = lognormdist(D, D_s)
            m_u_t = lognormdist(m_u, m_u_s)
            m_p_t = lognormdist(m_p, m_p_s)
            k_t = lognormdist(self.k, self.k_std)
            d70_t = lognormdist(self.d70, self.d70_std)
            i_ch_t = lognormdist(i_ch, i_ch_s)
            result = np.array([self.Z_u(h_exit_t, D_cover_t, m_u_t), self.Z_h(h_exit_t, D_cover_t, i_ch_t), self.Z_p(h_exit_t, D_cover_t, L_t, D_t, m_p_t, k_t, d70_t)])
            result_p += int(np.sum(result < 0) == 3)
        result_p = result_p / n
        return result_p