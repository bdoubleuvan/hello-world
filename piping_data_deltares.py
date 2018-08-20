# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 14:07:18 2018

@author: nguyen
"""

from scipy import stats
import pandas as pd
import numpy as np
import math as m
import tqdm
import openturns as ot
import hrd

class dijk:
    """
    This class uses dike data from deltares, data from LIWO overstromingskans and data from Hydra
    The functions in this class are based on the Sellmeijer Model to calculate dike failure probabilities.
    """
    
    path = r'D:\Users\Nguyen\Desktop\piping_data\pipingmatrix2.xlsx'
    path_dijk_lengte = r'D:\Users\Nguyen\Desktop\piping_data\static_information_geodata.geo_overstromingskans_2015_2020_faalkans_per_vak.xlsx'
    
    data = pd.read_excel(path, index_col = [0,1])
    data_dijk_lengte = pd.read_excel(path_dijk_lengte, index_col = 0, header = 0)
    
    hernoem_dict = {'%Seepage ditch waterlevel (42)_Mean' : 'h_exit', 
                '%Thickness of aquitard layer (44)_Mean' : 'D_cover',
                '%Seepage Length (48)_Mean': 'L',
                '%Thickness of upper sand layer (49)_Mean' : 'D',
                '%Volumic weight of the sand particles (50)_Mean' : 'gamma_sat',
                '%Model factor for Sellmeijer (51)_Mean' : 'm_p',
                '%Friction angle (Theta) (52)_Mean' : 'theta',
                '%Constant of "White" (53)_Mean' : 'eta',
                '%Permeability of the upper sand layer (55)_Mean' : 'k',
                '%Particle diameter (d70) (56)_Mean' : 'd70',
                '%Seepage ditch waterlevel (42)_Std' : 'h_exit_s',
                '%Thickness of aquitard layer (44)_Std' : 'D_cover_s',
                '%Seepage Length (48)_Std' : 'L_s',
                '%Thickness of upper sand layer (49)_Std' : 'D_s', 
                '%Volumic weight of the sand particles (50)_Std' : 'gamma_sat_s',
                '%Model factor for Sellmeijer (51)_Std' : 'm_p_s',
                '%Friction angle (Theta) (52)_Std' : 'theta_s',
                '%Constant of "White" (53)_Std' : 'eta_s',
                '%Permeability of the upper sand layer (55)_Std' : 'k_s',
                '%Particle diameter (d70) (56)_Std' : 'd70_s'
                }

    # Selecteer relevante kolommen
    perc_cols = [column for column in data.columns if column.startswith('%')]
    
#%% Data
    
    #shared dike data
    
    #deterministic
    d70m = 2.08e-4 #70% quantile of grain size
    g = 9.81 #gravitational constant
    gamma_water = 10 #volumetric weight of water
    r_c = 0.3 #reduction factor
    r_exit = 1 #damping factor at exit
    v_water = 1.33e-6 #kinematic viscosity of water
    
    #stochastic
    i_ch = 0.5 #critical heave gradient mu
    i_ch_s = 0.1 #critical heave gradient sigma
    m_u = 1.0 #model factor for uplift
    m_u_s = 0.1
    
    #reliability index data
    a = 0.75
    b = 300
    f = 0.24
    beta_norm = {
            300 : 2.71,
            1000 : 3.09,
            3000 : 3.4,
            10000 : 3.72,
            30000 : 3.99,
            100000 : 4.26
            }
    
    def __init__(self, vak_id, scenario):
        """ 
        This constructor has as input:
            vak_id, which states the dike section.
            scenario, which states the scenario used on the dike section.
            
        Moreover, constructor binds relevant data from the collected data resources. In our case that is Deltares data.
        """
        
        #dike segment
        self.vak_id = vak_id
        self.scenario = scenario
        
        #probabilites of water levels
        self.x = self.data['x-coordinate'][self.vak_id, self.scenario]
        self.y = self.data['y-coordinate'][self.vak_id, self.scenario]        
        path_h_freq = r'C:\MyPrograms\Hydra-NL\werkmap\WBI2017_Waddenzee_West_4-1_4-2_v03\Copy_WBI2017_Waddenzee_West_4-1_4-2_v03.sqlite'
        dir_h_freq = hrd.export_tables((self.x, self.y), path_h_freq)
        path_h_freq = r'C:\MyPrograms\Hydra-NL\werkmap\WBI2017_Waddenzee_West_4-1_4-2_v03\{}\Berekeningen\test\hfreq.txt'.format(dir_h_freq[0].values[0,2])
        data = pd.read_csv(path_h_freq, header = None, skiprows = 1, sep = "\s+")
        data.columns = ['waterstand', 'kans']    
        self.h_hydra = data['waterstand']
        self.p_h = data['kans']
                
        self.data_vak = self.data.loc[(self.vak_id, self.scenario), self.perc_cols]
        self.data_vak.rename(index = self.hernoem_dict, inplace = True)
        
        #Outside water level
        self.h = self.data.loc[(self.vak_id, self.scenario),'Water level with T = norm']
        self.h_norm = self.data.loc[(self.vak_id, self.scenario),'Water level with T = norm']
        
        #additional deterministic values
        self.eta = self.data_vak['eta'] #White's drag coefficient
        self. gamma_sp = self.data_vak['gamma_sat'] - self.gamma_water #volumetric weight of submerged sand particle
        self.gamma_sat = self.data_vak['gamma_sat'] #volumetric weight of saturated sand
        self.L_segm = self.data_dijk_lengte.loc[self.vak_id, 'LENGTE_M'] #(self.f / (stats.norm.cdf(-self.data.loc[(self.vak_id, self.scenario),'Required reliability index (a = 0.75 and b = 300)']) * self.T) - 1) * self.b / self.a 
        self.theta = self.data_vak['theta'] #bedding angle of sand grains
        
        #stochastic value
        self.D = self.data_vak['D'] #cover thickness of upper sand layer
        self.D_s = self.data_vak['D_s']
        self.D_cover = self.data_vak['D_cover'] #cover thickness of aquifer layer
        self.D_cover_s = self.data_vak['D_cover_s']
        self.d70 = self.data_vak['d70'] #grain size 
        self.d70_s = self.data_vak['d70_s']
        self.h_exit = self.data_vak['h_exit'] #seepage ditch level
        self.h_exit_s = self.data_vak['h_exit_s']
        self.k = self.data_vak['k'] #permeability of grain
        self.k_s = self.data_vak['k_s']
        self.L = self.data_vak['L'] #seepage length
        self.L_s = self.data_vak['L_s']
        self.m_p = self.data_vak['m_p'] #model factor for piping
        self.m_p_s = self.data_vak['m_p_s']  
        self.phi_polder = self.data_vak['h_exit'] #polder level == seepage ditch level
        self.phi_polder_s = self.data_vak['h_exit_s']
        
        #reliability index data
        self.T = self.data.loc[(self.vak_id, self.scenario),'Safety norm']
        self.p_norm = 1 / self.T
        self.beta_norm = self.beta_norm[self.T]
        self.p_T = self.f * self.p_norm
        self.p_T_cross = self.p_T_cross_function()
        self.beta_T = -stats.norm.ppf(self.p_T)
        self.beta_T_cross = self.beta_T_cross_function()
        
        self.Z_u_function = ot.PythonFunction(3, 1, self.Z_u_ot)
        self.Z_h_function = ot.PythonFunction(3, 1, self.Z_h_ot)
        self.Z_p_function = ot.PythonFunction(7, 1, self.Z_p_ot)
        self.Z_function = ot.PythonFunction(9, 1, self.Z)
        self.Z_MC_function = ot.PythonFunction(9, 1, self.Z_MC)
    
#%% Functions    
    def transform(self, mu, sigma):
        mu_t = m.log(mu / m.sqrt(1 + sigma / mu**2))
        sigma_t = m.sqrt(m.log(1 + sigma / mu**2))
        return [mu_t, sigma_t]

    def delta_phi_cu(self, D_cover):
        return D_cover * (self.gamma_sat - self.gamma_water) / self.gamma_water
    
    def phi_exit(self, h_exit):
        return h_exit + self.r_exit * (self.h - h_exit)
    
    def F_res(self):
        return self.eta * (self.gamma_sp / self.gamma_water) * m.tan(self.theta * m.pi / 180)
    
    def F_scale(self, d70, k, L):
        kappa = (self.v_water / self.g) * k
        return self.d70m / (kappa * L)**(1 / 3) * (d70 / self.d70m)**0.4
    
    def F_geo(self, D, L):
        if (D != L):
            result = 0.91 * (D / L)**(0.28 / ((D / L)**2.8 - 1) + 0.04)
        else:
            result = 1
        return result
    
    def H_c(self, D, d70, k, L):
        return L * self.F_res() * self.F_scale(d70, k, L) * self.F_geo(D, L)
        
    def Z_u(self, D_cover, h_exit, m_u):
        return m_u * self.delta_phi_cu(D_cover) - (self.phi_exit(h_exit) - h_exit)
    
    def Z_h(self, D_cover, h_exit, i_ch):
        return i_ch - (self.phi_exit(h_exit) - h_exit) / D_cover
    
    def Z_p(self, D, D_cover, d70, h_exit, k, L, m_p):
        return m_p * self.H_c(D, d70, k, L) - (self.h - h_exit - self.r_c * D_cover) 
        
#%% openTURNS function
    
    def Z_u_ot(self, X):
        D_cover, h_exit, m_u = X
        return [self.Z_u(D_cover, h_exit, m_u)]
    
    def Z_h_ot(self, X):
        D_cover, h_exit, i_ch = X
        return [self.Z_h(D_cover, h_exit, i_ch)]
    
    def Z_p_ot(self, X):
        D, D_cover, d70, h_exit, k, L, m_p = X
        return [self.Z_p(D, D_cover, d70, h_exit, k, L, m_p)]

    def Z(self, X):
        D, D_cover, d70, h_exit, i_ch, k, L, m_p, m_u = X
        return [max(self.Z_u(D_cover, h_exit, m_u), self.Z_h(D_cover, h_exit, i_ch), self.Z_p(D, D_cover, d70, h_exit, k, L, m_p))]
    
    def Z_MC(self, X):
        D, D_cover, d70, h_exit, i_ch, k, L, m_p, m_u = X
        return [max(self.Z_u(D_cover, h_exit, m_u), self.Z_h(D_cover, h_exit, i_ch), self.Z_p(D, D_cover, d70, h_exit, k, L, m_p)) >= 0]
    
            
#%% Reliability index functions
    
    def gamma_up(self):
        return self.D_cover * self.gamma_sp / (self.gamma_water * (self.h_norm - self.h_exit) * self.r_exit)
        
    def gamma_he(self):
        return self.D_cover * self.i_ch / ((self.h_norm - self.h_exit) * self.r_c)
    
    def gamma_pip(self):
        return self.H_c(self.k, self.d70, self.L, self.D) / (self.h_norm - self.h_exit - self.r_c * self.D_cover) 
    
    def beta_up(self):
        if self.gamma_up() != 0:
            result = 1 / 0.46 * (m.log(self.gamma_up() / 0.48) + 0.27 * self.beta_norm)
        else:
            result = -m.inf
        return result
    
    def beta_he(self):
        if self.gamma_he() != 0:
            result = 1 / 0.48 * (m.log(self.gamma_he() / 0.37) + 0.3 * self.beta_norm)
        else:
            result = -m.inf
        return result
    
    def beta_pip(self):
        if self.gamma_pip() != 0:
            result = 1 / 0.37 * (m.log(self.gamma_pip() / 1.04) + 0.43 * self.beta_norm)
        else:
            result = -m.inf
        return result
    
    def beta_f(self):
        return max(self.beta_up(), self.beta_he(), self.beta_pip())
    
    def p_f(self):
        return stats.norm.cdf(-self.beta_f())
    
    def p_T_cross_function(self):
        return self.f / (self.T * (1 + self.a * self.L_segm / self.b))
    
    def beta_T_cross_function(self):
        return -stats.norm.ppf(self.p_T_cross)
    
#%% Monte Carlo function
    
    def mc(self, n):
        lognormdist = np.random.lognormal    
        result = np.zeros(3)
        result_p = 0
        [D_cover, D_cover_s] = self.transform(self.D_cover, self.D_cover_s)
        [L, L_s] = self.transform(self.L, self.L_s)
        [D, D_s] = self.transform(self.D, self.D_s)
        [m_u, m_u_s] = self.transform(self.m_u, self.m_u_s)
        [m_p, m_p_s] = self.transform(self.m_p, self.m_p_s)
        [i_ch, i_ch_s] = self.transform(self.i_ch, self.i_ch_s)
        for i in tqdm.tqdm(range(n)):
            h_exit_t = np.random.normal(self.h_exit, self.h_exit_s)
            #phi_polder = ot.Normal(self.phi_polder,self.phi_polder_s)
            D_cover_t = lognormdist(D_cover, D_cover_s)
            L_t = lognormdist(L, L_s)
            D_t = lognormdist(D, D_s)
            m_u_t = lognormdist(m_u, m_u_s)
            m_p_t = lognormdist(m_p, m_p_s)
            k_t = lognormdist(self.k, self.k_s)
            d70_t = lognormdist(self.d70, self.d70_s)
            i_ch_t = lognormdist(i_ch, i_ch_s)
            result = np.array([self.Z_u(h_exit_t, D_cover_t, m_u_t), self.Z_h(h_exit_t, D_cover_t, i_ch_t), self.Z_p(h_exit_t, D_cover_t, L_t, D_t, m_p_t, k_t, d70_t)])
            result_p += int(np.sum(result < 0) == 3)
        return result_p / n
    
#%%
    def run(self, ds, form, mc, mc_l, sorm):
        """
        This function provides the following algorithms used from openTURNS:
            - Directional Sampling (DS)
            - First Order Reliability Method (FORM)
            - Monte Carlo Method (MC)
            - Second Order Reliability Method (SORM)
            
        The function run(ds, form, mc, sorm) has the following arguments
            ds, if assigned value "DS" then Directional Sampling is used, otherwise it is not used.
            form, if assigned "FORM" then FORM is used, otherwise it is not used.
            mc, if assigned "MC" then Monte Carlo method is used, otherwise it is not used.
            sorm, if assigned "SORM" then SORM is used, otherwise it is not used.
        """ 
#        self.path_cond = r'D:\Users\Nguyen\Desktop\piping_data\hfreq_{}.txt'.format(self.vak_id)
#        self.cond = pd.read_csv(self.path_cond, header = None, skiprows = 1, sep = "\s+")
#        self.cond.columns = ['waterhoogte', 'kans']    
#            
#        self.h_hydra = self.cond['waterhoogte'] #waterhoogte berekend via hydra
#        self.p_h = self.cond['kans'] #kans op waterhoogte hydra
        
        #self.variables = {} 
        #self.distribution = {}
        #for key in arg:
        #   self.variables['{}'.format(key)] = self.data_vak['{}'.format(key), '{}'.format(key + '_s')]
        #   self.distribution['{}'.format(key)] = ot.LogNormalMuSigma(self.variables['{}'.format(key)][0], self.variables['{}'.format(key)][1])
            
        #create joint distribution of parameters
        lognormdist = ot.LogNormalMuSigma
        self.distribution_D = lognormdist(self.D,self.D_s,0).getDistribution()
        self.distribution_D_cover = lognormdist(self.D_cover,self.D_cover_s).getDistribution()
        self.distribution_d70 = lognormdist(self.d70,self.d70_s,0).getDistribution()
        self.distribution_h_exit = ot.Normal(self.h_exit,self.h_exit_s)
        self.distribution_i_ch = lognormdist(self.i_ch,self.i_ch_s,0).getDistribution()
        self.distribution_k = lognormdist(self.k,self.k_s,0).getDistribution()
        self.distribution_L = lognormdist(self.L,self.L_s,0).getDistribution()
        self.distribution_m_p = lognormdist(self.m_p,self.m_p_s,0).getDistribution()
        self.distribution_m_u = lognormdist(self.m_u,self.m_u_s,0).getDistribution()
        self.marginals = [self.distribution_D, self.distribution_D_cover,
                          self.distribution_d70, self.distribution_h_exit,
                          self.distribution_i_ch, self.distribution_k,
                          self.distribution_L, self.distribution_m_p, 
                          self.distribution_m_u]
                
        if ds == "DS" or form == "FORM" or sorm == "SORM":
            #create copula
            self.RS = ot.CorrelationMatrix(len(self.marginals))
#            self.RS[2, 4] = 0.8 #correlation between permeability and grain size
#            self.RS[4, 2] = 0.8
            self.R = ot.NormalCopula.GetCorrelationFromSpearmanCorrelation(self.RS)
            self.copula = ot.NormalCopula(self.R)
            
            #create joint probability distribution
            self.distribution = ot.ComposedDistribution(self.marginals, self.copula)
            self.distribution.setDescription(['h_exit', 'D_cover', 'm_u', 'i_ch', 'L', 'D', 'm_p', 'k', 'd70'])
            
            #create the event we want to estimate the probability
            self.vect = ot.RandomVector(self.distribution)
            self.G = ot.CompositeRandomVector(self.Z_function, self.vect)
            self.event = ot.Event(self.G, ot.Less(), 0.0)
            self.event.setName('overall failure')
            
            if ds == "DS":
                #initialize list for conditional probability and amount of iterations 
                self.p_cond_DS = []
                self.iteration_DS = []
                
                #asks settings input for MC algorithm
#                self.default_DS = input('Would you like to use default  settings for DS? (True/False): ') or True
#                if self.default_DS == False:
#                    self.n_DS = input('How many iterations for MC: ') or 100000 #default 100000
#                    self.n_DS = int(self.n_DS)
#                    self.blocksize_DS = input('Set block size for MC: ') or 1 #default 1
#                    self.blocksize_DS = int(self.blocksize_DS)
#                    self.CoV_DS = input('Set the maximum for the coeffiecient of variation of MC: ') or 0.05 #default 0.05
#                    self.CoV_DS = float(self.CoV_DS)
#                else:
                self.n_DS = 100000 #default 100000
                self.blocksize_DS = 1 #default 1
                self.CoV_DS = 0.05 #default 0.05
                
                #root finding algorithm
                self.solver = ot.Bisection()
                self.rootStrategy = ot.RiskyAndFast(self.solver)
                
                #directional sampling algorithm
                self.samplingStrategy = ot.RandomDirection()
                    
                for h in self.h_hydra:
                    self.h = h
                    self.initialNumberOfCall_DS = self.Z_function.getEvaluationCallsNumber()
                        
                    #Using Directional sampling
                    self.algo_DS = ot.DirectionalSampling(self.event, self.rootStrategy, self.samplingStrategy)
                    self.algo_DS.setMaximumOuterSampling(self.n_DS)
                    self.algo_DS.setBlockSize(self.blocksize_DS)
                    self.algo_DS.setMaximumCoefficientOfVariation(self.CoV_DS)
                    #For statistics about the algorithm
                    self.initialNumberOfCall_DS = self.Z_function.getEvaluationCallsNumber()
                    
                    #Perform the analysis
                    self.algo_DS.run()
                    
                    self.result_DS = self.algo_DS.getResult()
                    self.p_cond_DS.append(self.result_DS.getProbabilityEstimate())
                    self.iteration_DS.append(self.Z_function.getEvaluationCallsNumber() - self.initialNumberOfCall_DS)
                
                #probability failure
                self.p_f_DS = sum([a * b for a, b in zip(self.p_cond_DS, self.p_h)]) / len(self.p_h)

                    
            
            #general settings for FORM and SORM
#            if form == "FORM" or sorm == "SORM":                
#                #Asks amount of max iterations
#                self.n = input('Set max amount of iterations for FORM and SORM: ') or 10000
#                self.n = int(self.n)
            self.n = 10000
            
            if form == "FORM":
                #initialize list for conditional probability and amount of iterations               
                self.p_cond_FORM = []
                self.iteration_FORM = []
                self.Hasofer_FORM = []
                self.standard_design_FORM = []
                self.physical_design_FORM = []
            
                #define a solver
                self.optimAlgo_FORM = ot.AbdoRackwitz()
#                self.optimAlgo_FORM = ot.Cobyla()
#                self.optimAlgo_FORM.setMaximumEvaluationNumber(self.n)
#                self.optimAlgo_FORM.setMaximumAbsoluteError(1.0e-10)
#                self.optimAlgo_FORM.setMaximumRelativeError(1.0e-10)
#                self.optimAlgo_FORM.setMaximumResidualError(1.0e-10)
#                self.optimAlgo_FORM.setMaximumConstraintError(1.0e-10)
                self.algo_FORM = ot.FORM(self.optimAlgo_FORM, self.event, self.distribution.getMean())
                
                for h in self.h_hydra:
                    self.h = h
                    self.initialNumberOfCall_FORM = self.Z_function.getEvaluationCallsNumber()
                                    
                    #run FORM
                    self.algo_FORM.run()
                   
                    #get result
                    self.result_FORM = self.algo_FORM.getResult()
                    self.p_cond_FORM.append(self.result_FORM.getEventProbability())
                    self.Hasofer_FORM.append(self.result_FORM.getHasoferReliabilityIndex())
                    self.standard_design_FORM.append(self.result_FORM.getStandardSpaceDesignPoint())
                    self.physical_design_FORM.append(self.result_FORM.getPhysicalSpaceDesignPoint())
                    self.iteration_FORM.append(self.Z_function.getEvaluationCallsNumber() - self.initialNumberOfCall_FORM)
                
                #probability failure
                self.p_f_FORM = sum([a * b for a, b in zip(self.p_cond_FORM, self.p_h)]) / len(self.p_h)

            
            if sorm == "SORM":
                #initialize list for conditional probability and amount of iterations   
                self.p_cond_Breitung = []
                self.p_cond_Hohen_Bichler = []
                self.p_cond_Tvedt = []
                self.iteration_SORM = []
                self.Breitung = []
                self.Hohen_Bichler = []
                self.Tvedt = []
                
                #define a solver
                self.optimAlgo_SORM = ot.AbdoRackwitz()
                self.optimAlgo_SORM.setMaximumEvaluationNumber(self.n)
                self.algo_SORM = ot.SORM(self.optimAlgo_SORM, self.event, self.distribution.getMean())
                
                for h in self.h_hydra:
                    self.h = h
                    self.initialNumberOfCall_SORM = self.Z_function.getEvaluationCallsNumber()
                    
                    #run SORM
                    self.algo_SORM.run()
                    
                    #get result
                    self.result_SORM = self.algo_SORM.getResult()
                    self.p_cond_Breitung.append( self.result_SORM.getEventProbabilityBreitung())
                    self.p_cond_Hohen_Bichler.append(self.result_SORM.getEventProbabilityHohenBichler())
                    self.p_cond_Tvedt.append(self.result_SORM.getEventProbabilityTvedt())
                    self.Breitung.append(self.result_SORM.getGeneralisedReliabilityIndexBreitung())
                    self.Hohen_Bichler.append(self.result_SORM.getGeneralisedReliabilityIndexHohenBichler())
                    self.Tvedt.append(self.result_SORM.getGeneralisedReliabilityIndexTvedt())
                    self.iteration_SORM.append(self.Z_function.getEvaluationCallsNumber() - self.initialNumberOfCall_SORM)
                    
                #probability failure
                self.p_f_Breitung = sum([a * b for a, b in zip(self.p_cond_Breitung, self.p_h)]) / len(self.p_h)
                self.p_f_Hohen_Bichler = sum([a * b for a, b in zip(self.p_cond_Hohen_Bichler, self.p_h)]) / len(self.p_h)
                self.p_f_Tvedt = sum([a * b for a, b in zip(self.p_cond_Tvedt, self.p_h)]) / len(self.p_h)

                
        if mc == "MC":
            #initialize list for conditional probability and amount of iterations               
            self.p_cond_MC = []
            self.iteration_MC = []
            
            #create combined event we want to estimate the probability
            self.distribution_MC = ot.ComposedDistribution(self.marginals)
            self.vect_MC = ot.RandomVector(self.distribution_MC)
            self.G_MC = ot.CompositeRandomVector(self.Z_function, self.vect_MC)
            self.event_MC = ot.Event(self.G_MC, ot.Less(), 0) 
            
            #asks settings input for MC algorithm
#            self.default_MC = input('Would you like to use default  settings for MC? (y/n): ') or 0
#            if self.default_MC == "n":
#                self.n_MC = input('How many iterations for MC: ') or 100000 #default 100000
#                self.n_MC = int(self.n_MC)
#                self.blocksize_MC = input('Set block size for MC: ') or 50 #default 50
#                self.blocksize_MC = int(self.blocksize_MC)
#                self.CoV_MC = input('Set the maximum for the coeffiecient of variation of MC: ') or 0.05 #default 0.05
#                self.CoV_MC = float(self.CoV_MC)
#            else:
            self.n_MC = 100000 #default 100000
            self.blocksize_MC = 50 #default 50
            self.CoV_MC = 0.05 #default 0.05
            
            #create Monte Carlo algorithm for combined mechanism
            self.experiment_MC = ot.MonteCarloExperiment()
            self.algo_MC = ot.ProbabilitySimulationAlgorithm(self.event_MC, self.experiment_MC)
            self.algo_MC.setMaximumCoefficientOfVariation(self.CoV_MC)
            self.algo_MC.setMaximumOuterSampling(self.n_MC)
            self.algo_MC.setBlockSize(self.blocksize_MC)
            
            for h in self.h_hydra:
                self.h = h
                self.initialNumberOfCall_MC = self.Z_function.getEvaluationCallsNumber()
                
                self.algo_MC.run()
            
                self.result_MC = self.algo_MC.getResult()
                self.p_cond_MC.append(self.result_MC.getProbabilityEstimate())
                self.iteration_MC.append(self.Z_function.getEvaluationCallsNumber() - self.initialNumberOfCall_MC)
        
            #probability failure
            self.p_f_MC = sum([a * b for a, b in zip(self.p_cond_MC, self.p_h)]) / len(self.p_h)
        
        if mc_l == "MC_L":
            #initialize list for conditional probability and amount of iterations               
            self.p_cond_MC_L = []
            self.iteration_MC_L = []
            
            #create combined event we want to estimate the probability
            self.distribution_MC_L = ot.ComposedDistribution(self.marginals)
            self.vect_MC_L = ot.RandomVector(self.distribution_MC)
            self.G_MC_L = ot.CompositeRandomVector(self.Z_MC_function, self.vect_MC)
            self.event_MC_L = ot.Event(self.G_MC, ot.Equal(), 0) 
            
#            #asks settings input for MC algorithm
#            self.default_MC_L = input('Would you like to use default  settings for MC? (y/n): ') or 0
#            if self.default_MC_L == "n":
#                self.n_MC_L = input('How many iterations for MC: ') or 100000 #default 100000
#                self.n_MC_L = int(self.n_MC_L)
#                self.blocksize_MC_L = input('Set block size for MC: ') or 50 #default 50
#                self.blocksize_MC_L = int(self.blocksize_MC_L)
#                self.CoV_MC_L = input('Set the maximum for the coeffiecient of variation of MC: ') or 0.05 #default 0.05
#                self.CoV_MC_L = float(self.CoV_MC_L)
#            else:
            self.n_MC_L = 100000 #default 100000
            self.blocksize_MC_L = 50 #default 50
            self.CoV_MC_L = 0.05 #default 0.05
            
            #create Monte Carlo algorithm for combined mechanism
            self.experiment_MC_L = ot.MonteCarloExperiment()
            self.algo_MC_L = ot.ProbabilitySimulationAlgorithm(self.event_MC_L, self.experiment_MC_L)
            self.algo_MC_L.setMaximumCoefficientOfVariation(self.CoV_MC_L)
            self.algo_MC_L.setMaximumOuterSampling(self.n_MC_L)
            self.algo_MC_L.setBlockSize(self.blocksize_MC_L)
            
            for h in self.h_hydra:
                self.h = h
                self.initialNumberOfCall_MC_L = self.Z_MC_function.getEvaluationCallsNumber()
                
                self.algo_MC_L.run()
            
                self.result_MC_L = self.algo_MC_L.getResult()
                self.p_cond_MC_L.append(self.result_MC_L.getProbabilityEstimate())
                self.iteration_MC_L.append(self.Z_MC_L_function.getEvaluationCallsNumber() - self.initialNumberOfCall_MC_L)
        
            #probability failure
            self.p_f_MC_L = sum([a * b for a, b in zip(self.p_cond_MC_L, self.p_h)]) / len(self.p_h)
        
    
 #%%
    def run_wbi(self, ds, form, mc, mc_l, sorm):
        #create joint distribution of parameters
        lognormdist = ot.LogNormalMuSigma
        self.distribution_D = lognormdist(self.D,self.D_s,0).getDistribution()
        self.distribution_D_cover = lognormdist(self.D_cover,self.D_cover_s).getDistribution()
        self.distribution_d70 = lognormdist(self.d70,self.d70_s,0).getDistribution()
        self.distribution_h_exit = ot.Normal(self.h_exit,self.h_exit_s)
        self.distribution_i_ch = lognormdist(self.i_ch,self.i_ch_s,0).getDistribution()
        self.distribution_k = lognormdist(self.k,self.k_s,0).getDistribution()
        self.distribution_L = lognormdist(self.L,self.L_s,0).getDistribution()
        self.distribution_m_p = lognormdist(self.m_p,self.m_p_s,0).getDistribution()
        self.distribution_m_u = lognormdist(self.m_u,self.m_u_s,0).getDistribution()
        self.marginals_u = [self.distribution_D_cover, self.distribution_h_exit, self.distribution_m_u]
        self.marginals_h = [self.distribution_D_cover, self.distribution_h_exit, self.distribution_i_ch]
        self.marginals_p = [self.distribution_D, self.distribution_D_cover,
                            self.distribution_d70, self.distribution_h_exit,
                            self.distribution_k, self.distribution_L, 
                            self.distribution_m_p]
                
        if ds == "DS" or form == "FORM" or sorm == "SORM":
            #create copula uplift
            self.RS_u = ot.CorrelationMatrix(len(self.marginals_u))
            self.R_u = ot.NormalCopula.GetCorrelationFromSpearmanCorrelation(self.RS_u)
            self.copula_u = ot.NormalCopula(self.R_u)
        
            #create copula heave
            self.RS_h = ot.CorrelationMatrix(len(self.marginals_h))
            self.R_h = ot.NormalCopula.GetCorrelationFromSpearmanCorrelation(self.RS_h)
            self.copula_h = ot.NormalCopula(self.R_h)
        
            #create copula piping
            self.RS_p = ot.CorrelationMatrix(len(self.marginals_p))
            self.R_p = ot.NormalCopula.GetCorrelationFromSpearmanCorrelation(self.RS_p)
            self.copula_p = ot.NormalCopula(self.R_p)
            
          
            #marginals uplift
            self.distribution_u = ot.ComposedDistribution(self.marginals_u)
            self.distribution_u.setDescription(['h_exit', 'D_cover', 'm_u'])
            #marginals heave
            self.distribution_h = ot.ComposedDistribution(self.marginals_h)
            self.distribution_h.setDescription(['h_exit', 'D_cover', 'i_ch'])
            #marginals piping
            self.distribution_p = ot.ComposedDistribution(self.marginals_p)
            self.distribution_p.setDescription(['h_exit', 'D_cover', 'L', 'D', 'm_p', 'k', 'd70'])
            
            #create the event we want to estimate the probability uplift
            self.vect_u = ot.RandomVector(self.distribution_u)
            self.G_u = ot.CompositeRandomVector(self.Z_u_function, self.vect_u)
            self.event_u = ot.Event(self.G_u, ot.Less(), 0.0)
            self.event_u.setName('uplift failure')
            
            #create the event we want to estimate the probability heave
            self.vect_h = ot.RandomVector(self.distribution_h)
            self.G_h = ot.CompositeRandomVector(self.Z_h_function, self.vect_h)
            self.event_h = ot.Event(self.G_h, ot.Less(), 0.0)
            self.event_h.setName('heave failure')
        
            #create the event we want to estimate the probability piping
            self.vect_p = ot.RandomVector(self.distribution_p)
            self.G_p = ot.CompositeRandomVector(self.Z_p_function, self.vect_p)
            self.event_p = ot.Event(self.G_p, ot.Less(), 0.0)
            self.event_p.setName('piping failure')
            
        #            self.uplift = self.event(self.vak_id, self.scenario, ['h_exit', 'D_cover', 'm_u'], 'uplift')
        #            self.heave = self.event(self.vak_id, self.scenario, ['h_exit', 'D_cover', 'i_ch'], 'heave')
        #            self.piping = self.event(self.vak_id, self.scenario, ['h_exit', 'D_cover', 'L', 'D', 'm_p', 'k', 'd70'], 'piping')
            
            #general settings for FORM and SORM
#            if form == "FORM" or sorm == "SORM":                
#                #Asks amount of max iterations
#                self.n = input('Set max amount of iterations for FORM and SORM: ') or 10000
#                self.n = int(self.n)
            self.n = 10000
            
            if form == "FORM":
                #initialize list for conditional probability and amount of iterations               
                self.p_cond_u_FORM = []
                self.p_cond_h_FORM = []
                self.p_cond_p_FORM = []
                self.iteration_u_FORM = []
                self.iteration_h_FORM = []
                self.iteration_p_FORM = []
                self.Hasofer_u_FORM = []
                self.Hasofer_h_FORM = []
                self.Hasofer_p_FORM = []
                self.standard_design_u_FORM = []
                self.standard_design_h_FORM = []
                self.standard_design_p_FORM = []
                self.physical_design_u_FORM = []
                self.physical_design_h_FORM = []
                self.physical_design_p_FORM = []
            
                #define a solver
                self.optimAlgo_FORM = ot.AbdoRackwitz()
#                self.optimAlgo_FORM = ot.Cobyla()
#                self.optimAlgo_FORM.setMaximumEvaluationNumber(self.n)
#                self.optimAlgo_FORM.setMaximumAbsoluteError(1.0e-10)
#                self.optimAlgo_FORM.setMaximumRelativeError(1.0e-10)
#                self.optimAlgo_FORM.setMaximumResidualError(1.0e-10)
#                self.optimAlgo_FORM.setMaximumConstraintError(1.0e-10)
                self.algo_u_FORM = ot.FORM(self.optimAlgo_FORM, self.event_u, self.distribution_u.getMean())
                self.algo_h_FORM = ot.FORM(self.optimAlgo_FORM, self.event_h, self.distribution_h.getMean())
                self.algo_p_FORM = ot.FORM(self.optimAlgo_FORM, self.event_p, self.distribution_p.getMean())
                
                for h in self.h_hydra:
                    self.h = h
                    self.initialNumberOfCall_u = self.Z_u_function.getEvaluationCallsNumber()
                    self.initialNumberOfCall_h = self.Z_h_function.getEvaluationCallsNumber()
                    self.initialNumberOfCall_p = self.Z_p_function.getEvaluationCallsNumber()
                    
                    #run FORM uplift
                    self.algo_u_FORM.run()
                    #get result
                    self.result_u_FORM = self.algo_u_FORM.getResult()
                    self.p_cond_u_FORM.append(self.result_u_FORM.getEventProbability())
                    self.Hasofer_u_FORM.append(self.result_u_FORM.getHasoferReliabilityIndex())
                    self.standard_design_u_FORM.append(self.result_u_FORM.getStandardSpaceDesignPoint())
                    self.physical_design_u_FORM.append(self.result_u_FORM.getPhysicalSpaceDesignPoint())
                    self.iteration_u_FORM.append(self.Z_u_function.getEvaluationCallsNumber() - self.initialNumberOfCall_u)
                
                    #run FORM heave
                    self.algo_h_FORM.run()
                    #get result
                    self.result_h_FORM = self.algo_h_FORM.getResult()
                    self.p_cond_h_FORM.append(self.result_h_FORM.getEventProbability())
                    self.Hasofer_h_FORM.append(self.result_h_FORM.getHasoferReliabilityIndex())
                    self.standard_design_h_FORM.append(self.result_h_FORM.getStandardSpaceDesignPoint())
                    self.physical_design_h_FORM.append(self.result_h_FORM.getPhysicalSpaceDesignPoint())
                    self.iteration_h_FORM.append(self.Z_h_function.getEvaluationCallsNumber() - self.initialNumberOfCall_h)
                    
                    #run FORM piping
                    self.algo_p_FORM.run()
                    #get result
                    self.result_p_FORM = self.algo_p_FORM.getResult()
                    self.p_cond_p_FORM.append(self.result_p_FORM.getEventProbability())
                    self.Hasofer_p_FORM.append(self.result_p_FORM.getHasoferReliabilityIndex())
                    self.standard_design_p_FORM.append(self.result_p_FORM.getStandardSpaceDesignPoint())
                    self.physical_design_p_FORM.append(self.result_p_FORM.getPhysicalSpaceDesignPoint())
                    self.iteration_p_FORM.append(self.Z_p_function.getEvaluationCallsNumber() - self.initialNumberOfCall_p)
                    
                
                #probability failure  
                self.p_cond_FORM = [min(a, b, c) for a, b, c in zip(self.p_cond_u_FORM, self.p_cond_h_FORM, self.p_cond_p_FORM)]
                self.p_f_FORM = sum([a * b for a, b in zip(self.p_cond_FORM, self.p_h)]) / len(self.p_h)
#                self.Hasofer_FORM = []
#                self.standard_design_FORM = []
#                self.physical_design_FORM = []
                self.iteration_FORM = [a + b + c for a, b, c in zip(self.iteration_u_FORM, self.iteration_h_FORM, self.iteration_p_FORM)]    
        
        if mc == "MC":
            #initialize list for conditional probability and amount of iterations               
            self.p_cond_u_MC = []
            self.p_cond_h_MC = []
            self.p_cond_p_MC = []
            self.iteration_u_MC = []
            self.iteration_h_MC = []
            self.iteration_p_MC = []
            
            #create combined event we want to estimate the probability
            #uplift
            self.distribution_u_MC = ot.ComposedDistribution(self.marginals_u)
            self.vect_u_MC = ot.RandomVector(self.distribution_u_MC)
            self.G_u_MC = ot.CompositeRandomVector(self.Z_u_function, self.vect_u_MC)
            self.event_u_MC = ot.Event(self.G_u_MC, ot.Less(), 0) 
            #heave
            self.distribution_h_MC = ot.ComposedDistribution(self.marginals_h)
            self.vect_h_MC = ot.RandomVector(self.distribution_h_MC)
            self.G_h_MC = ot.CompositeRandomVector(self.Z_h_function, self.vect_h_MC)
            self.event_h_MC = ot.Event(self.G_h_MC, ot.Less(), 0) 
            #piping
            self.distribution_p_MC = ot.ComposedDistribution(self.marginals_p)
            self.vect_p_MC = ot.RandomVector(self.distribution_p_MC)
            self.G_p_MC = ot.CompositeRandomVector(self.Z_p_function, self.vect_p_MC)
            self.event_p_MC = ot.Event(self.G_p_MC, ot.Less(), 0) 
            
            #asks settings input for MC algorithm
#            self.default_MC = input('Would you like to use default  settings for MC? (y/n): ') or 0
#            if self.default_MC == "n":
#                self.n_MC = input('How many iterations for MC: ') or 100000 #default 100000
#                self.n_MC = int(self.n_MC)
#                self.blocksize_MC = input('Set block size for MC: ') or 50 #default 50
#                self.blocksize_MC = int(self.blocksize_MC)
#                self.CoV_MC = input('Set the maximum for the coeffiecient of variation of MC: ') or 0.05 #default 0.05
#                self.CoV_MC = float(self.CoV_MC)
#            else:
            self.n_MC = 100000 #default 100000
            self.blocksize_MC = 50 #default 50
            self.CoV_MC = 0.05 #default 0.05
            
            self.experiment_MC = ot.MonteCarloExperiment()
            
            #create Monte Carlo algorithm for combined mechanism
            self.algo_u_MC = ot.ProbabilitySimulationAlgorithm(self.event_u_MC, self.experiment_MC)
            self.algo_u_MC.setMaximumCoefficientOfVariation(self.CoV_MC)
            self.algo_u_MC.setMaximumOuterSampling(self.n_MC)
            self.algo_u_MC.setBlockSize(self.blocksize_MC)
            
            #create Monte Carlo algorithm for combined mechanism
            self.algo_h_MC = ot.ProbabilitySimulationAlgorithm(self.event_h_MC, self.experiment_MC)
            self.algo_h_MC.setMaximumCoefficientOfVariation(self.CoV_MC)
            self.algo_h_MC.setMaximumOuterSampling(self.n_MC)
            self.algo_h_MC.setBlockSize(self.blocksize_MC)
            
            #create Monte Carlo algorithm for combined mechanism
            self.algo_p_MC = ot.ProbabilitySimulationAlgorithm(self.event_p_MC, self.experiment_MC)
            self.algo_p_MC.setMaximumCoefficientOfVariation(self.CoV_MC)
            self.algo_p_MC.setMaximumOuterSampling(self.n_MC)
            self.algo_p_MC.setBlockSize(self.blocksize_MC)
            
            for h in self.h_hydra:
                self.h = h
                
                self.initialNumberOfCall_u_MC = self.Z_u_function.getEvaluationCallsNumber()
                self.algo_u_MC.run()
                self.result_u_MC = self.algo_u_MC.getResult()
                self.p_cond_u_MC.append(self.result_u_MC.getProbabilityEstimate())
                self.iteration_u_MC.append(self.Z_u_function.getEvaluationCallsNumber() - self.initialNumberOfCall_u_MC)
                
                self.initialNumberOfCall_h_MC = self.Z_h_function.getEvaluationCallsNumber()
                self.algo_h_MC.run()
                self.result_h_MC = self.algo_h_MC.getResult()
                self.p_cond_h_MC.append(self.result_h_MC.getProbabilityEstimate())
                self.iteration_h_MC.append(self.Z_h_function.getEvaluationCallsNumber() - self.initialNumberOfCall_h_MC)
                
                self.initialNumberOfCall_p_MC = self.Z_p_function.getEvaluationCallsNumber()
                self.algo_p_MC.run()
                self.result_p_MC = self.algo_p_MC.getResult()
                self.p_cond_p_MC.append(self.result_p_MC.getProbabilityEstimate())
                self.iteration_p_MC.append(self.Z_p_function.getEvaluationCallsNumber() - self.initialNumberOfCall_p_MC)
        
            #probability failure
            self.p_cond_MC = [min(a, b, c) for a, b, c in zip(self.p_cond_u_MC, self.p_cond_h_MC, self.p_cond_p_MC)]
            self.p_f_MC = sum([a * b for a, b in zip(self.p_cond_MC, self.p_h)]) / len(self.p_h)
        

#%% Derived class for all scenarios   
    
class dijk_scenario(dijk):
    def __init__(self, vak_id):
        self.vak_id = vak_id
        self.scenario = {}
        self.index = [index for index in self.data.loc[self.vak_id].index]
        for index in self.index: 
            self.scenario[index] = dijk(self.vak_id, index)
        self.beta_cross = -stats.norm.ppf(self.p_cross())
            
    def p_cross(self):
        result = 0
        for i in self.index:
            result += self.scenario[i].p_f() 
        return result / len(self.index)
    
#%%
#        
#class dijk_hkv(dijk):
#    def __init__(self, vak_id, scenario):
#        super().__init__(vak_id, scenario)
#        
#    def run(self, ds, form, mc, mc_l, sorm):
#        """
#        This function provides the following algorithms used from openTURNS:
#            - Directional Sampling (DS)
#            - First Order Reliability Method (FORM)
#            - Monte Carlo Method (MC)
#            - Second Order Reliability Method (SORM)
#            
#        The function run(ds, form, mc, sorm) has the following arguments
#            ds, if assigned value "DS" then Directional Sampling is used, otherwise it is not used.
#            form, if assigned "FORM" then FORM is used, otherwise it is not used.
#            mc, if assigned "MC" then Monte Carlo method is used, otherwise it is not used.
#            sorm, if assigned "SORM" then SORM is used, otherwise it is not used.
#        """ 
#        self.path_cond = r'D:\Users\Nguyen\Desktop\piping_data\hfreq_{}.txt'.format(self.vak_id)
#        self.cond = pd.read_csv(self.path_cond, header = None, skiprows = 1, sep = "\s+")
#        self.cond.columns = ['waterhoogte', 'kans']    
#            
#        self.h_hydra = self.cond['waterhoogte'] #waterhoogte berekend via hydra
#        self.p_h = self.cond['kans'] #kans op waterhoogte hydra
#        
#        #self.variables = {} 
#        #self.distribution = {}
#        #for key in arg:
#        #   self.variables['{}'.format(key)] = self.data_vak['{}'.format(key), '{}'.format(key + '_s')]
#        #   self.distribution['{}'.format(key)] = ot.LogNormalMuSigma(self.variables['{}'.format(key)][0], self.variables['{}'.format(key)][1])
#            
#        #create joint distribution of parameters
#        lognormdist = ot.LogNormalMuSigma
#        self.distribution_D = lognormdist(self.D,self.D_s,0).getDistribution()
#        self.distribution_D_cover = lognormdist(self.D_cover,self.D_cover_s).getDistribution()
#        self.distribution_d70 = lognormdist(self.d70,self.d70_s,0).getDistribution()
#        self.distribution_h_exit = ot.Normal(self.h_exit,self.h_exit_s)
#        self.distribution_i_ch = lognormdist(self.i_ch,self.i_ch_s,0).getDistribution()
#        self.distribution_k = lognormdist(self.k,self.k_s,0).getDistribution()
#        self.distribution_L = lognormdist(self.L,self.L_s,0).getDistribution()
#        self.distribution_m_p = lognormdist(self.m_p,self.m_p_s,0).getDistribution()
#        self.distribution_m_u = lognormdist(self.m_u,self.m_u_s,0).getDistribution()
#        self.marginals = [self.distribution_D, self.distribution_D_cover,
#                          self.distribution_d70, self.distribution_h_exit,
#                          self.distribution_i_ch, self.distribution_k,
#                          self.distribution_L, self.distribution_m_p, 
#                          self.distribution_m_u]
#                
#        if ds == "DS" or form == "FORM" or sorm == "SORM":
#            #create copula
#            self.RS = ot.CorrelationMatrix(len(self.marginals))
##            self.RS[2, 4] = 0.8 #correlation between permeability and grain size
##            self.RS[4, 2] = 0.8
#            self.R = ot.NormalCopula.GetCorrelationFromSpearmanCorrelation(self.RS)
#            self.copula = ot.NormalCopula(self.R)
#            
#            #create joint probability distribution
#            self.distribution = ot.ComposedDistribution(self.marginals, self.copula)
#            self.distribution.setDescription(['h_exit', 'D_cover', 'm_u', 'i_ch', 'L', 'D', 'm_p', 'k', 'd70'])
#            
#            #create the event we want to estimate the probability
#            self.vect = ot.RandomVector(self.distribution)
#            self.G = ot.CompositeRandomVector(self.Z_function, self.vect)
#            self.event = ot.Event(self.G, ot.Less(), 0.0)
#            self.event.setName('overall failure')
#            
#            if ds == "DS":
#                #initialize list for conditional probability and amount of iterations 
#                self.p_cond_DS = []
#                self.iteration_DS = []
#                
#                #asks settings input for MC algorithm
#                self.default_DS = input('Would you like to use default  settings for DS? (True/False): ') or True
#                if self.default_DS == False:
#                    self.n_DS = input('How many iterations for MC: ') or 100000 #default 100000
#                    self.n_DS = int(self.n_DS)
#                    self.blocksize_DS = input('Set block size for MC: ') or 1 #default 1
#                    self.blocksize_DS = int(self.blocksize_DS)
#                    self.CoV_DS = input('Set the maximum for the coeffiecient of variation of MC: ') or 0.05 #default 0.05
#                    self.CoV_DS = float(self.CoV_DS)
#                else:
#                    self.n_DS = 100000 #default 100000
#                    self.blocksize_DS = 1 #default 1
#                    self.CoV_DS = 0.05 #default 0.05
#                
#                #root finding algorithm
#                self.solver = ot.Bisection()
#                self.rootStrategy = ot.RiskyAndFast(self.solver)
#                
#                #directional sampling algorithm
#                self.samplingStrategy = ot.RandomDirection()
#                    
#                for h in self.h_hydra:
#                    self.h = h
#                    self.initialNumberOfCall_DS = self.Z_function.getEvaluationCallsNumber()
#                        
#                    #Using Directional sampling
#                    self.algo_DS = ot.DirectionalSampling(self.event, self.rootStrategy, self.samplingStrategy)
#                    self.algo_DS.setMaximumOuterSampling(self.n_DS)
#                    self.algo_DS.setBlockSize(self.blocksize_DS)
#                    self.algo_DS.setMaximumCoefficientOfVariation(self.CoV_DS)
#                    #For statistics about the algorithm
#                    self.initialNumberOfCall_DS = self.Z_function.getEvaluationCallsNumber()
#                    
#                    #Perform the analysis
#                    self.algo_DS.run()
#                    
#                    self.result_DS = self.algo_DS.getResult()
#                    self.p_cond_DS.append(self.result_DS.getProbabilityEstimate())
#                    self.iteration_DS.append(self.Z_function.getEvaluationCallsNumber() - self.initialNumberOfCall_DS)
#                
#                #probability failure
#                self.p_f_DS = sum([a * b for a, b in zip(self.p_cond_DS, self.p_h)]) / len(self.p_h)
#
#                    
#            
#            #general settings for FORM and SORM
#            if form == "FORM" or sorm == "SORM":                
#                #Asks amount of max iterations
#                self.n = input('Set max amount of iterations for FORM and SORM: ') or 10000
#                self.n = int(self.n)
#            
#            if form == "FORM":
#                #initialize list for conditional probability and amount of iterations               
#                self.p_cond_FORM = []
#                self.iteration_FORM = []
#                self.Hasofer_FORM = []
#                self.standard_design_FORM = []
#                self.physical_design_FORM = []
#            
#                #define a solver
#                self.optimAlgo_FORM = ot.Cobyla()
#                self.optimAlgo_FORM.setMaximumEvaluationNumber(self.n)
#                self.optimAlgo_FORM.setMaximumAbsoluteError(1.0e-10)
#                self.optimAlgo_FORM.setMaximumRelativeError(1.0e-10)
#                self.optimAlgo_FORM.setMaximumResidualError(1.0e-10)
#                self.optimAlgo_FORM.setMaximumConstraintError(1.0e-10)
#                self.algo_FORM = ot.FORM(self.optimAlgo_FORM, self.event, self.distribution.getMean())
#                
#                for h in self.h_hydra:
#                    self.h = h
#                    self.initialNumberOfCall_FORM = self.Z_function.getEvaluationCallsNumber()
#                                    
#                    #run FORM
#                    self.algo_FORM.run()
#                   
#                    #get result
#                    self.result_FORM = self.algo_FORM.getResult()
#                    self.p_cond_FORM.append(self.result_FORM.getEventProbability())
#                    self.Hasofer_FORM.append(self.result_FORM.getHasoferReliabilityIndex())
#                    self.standard_design_FORM.append(self.result_FORM.getStandardSpaceDesignPoint())
#                    self.physical_design_FORM.append(self.result_FORM.getPhysicalSpaceDesignPoint())
#                    self.iteration_FORM.append(self.Z_function.getEvaluationCallsNumber() - self.initialNumberOfCall_FORM)
#                
#                #probability failure
#                self.p_f_FORM = sum([a * b for a, b in zip(self.p_cond_FORM, self.p_h)]) / len(self.p_h)
#
#            
#            if sorm == "SORM":
#                #initialize list for conditional probability and amount of iterations   
#                self.p_cond_Breitung = []
#                self.p_cond_Hohen_Bichler = []
#                self.p_cond_Tvedt = []
#                self.iteration_SORM = []
#                self.Breitung = []
#                self.Hohen_Bichler = []
#                self.Tvedt = []
#                
#                #define a solver
#                self.optimAlgo_SORM = ot.AbdoRackwitz()
#                self.optimAlgo_SORM.setMaximumEvaluationNumber(self.n)
#                self.algo_SORM = ot.SORM(self.optimAlgo_SORM, self.event, self.distribution.getMean())
#                
#                for h in self.h_hydra:
#                    self.h = h
#                    self.initialNumberOfCall_SORM = self.Z_function.getEvaluationCallsNumber()
#                    
#                    #run SORM
#                    self.algo_SORM.run()
#                    
#                    #get result
#                    self.result_SORM = self.algo_SORM.getResult()
#                    self.p_cond_Breitung.append( self.result_SORM.getEventProbabilityBreitung())
#                    self.p_cond_Hohen_Bichler.append(self.result_SORM.getEventProbabilityHohenBichler())
#                    self.p_cond_Tvedt.append(self.result_SORM.getEventProbabilityTvedt())
#                    self.Breitung.append(self.result_SORM.getGeneralisedReliabilityIndexBreitung())
#                    self.Hohen_Bichler.append(self.result_SORM.getGeneralisedReliabilityIndexHohenBichler())
#                    self.Tvedt.append(self.result_SORM.getGeneralisedReliabilityIndexTvedt())
#                    self.iteration_SORM.append(self.Z_function.getEvaluationCallsNumber() - self.initialNumberOfCall_SORM)
#                    
#                #probability failure
#                self.p_f_Breitung = sum([a * b for a, b in zip(self.p_cond_Breitung, self.p_h)]) / len(self.p_h)
#                self.p_f_Hohen_Bichler = sum([a * b for a, b in zip(self.p_cond_Hohen_Bichler, self.p_h)]) / len(self.p_h)
#                self.p_f_Tvedt = sum([a * b for a, b in zip(self.p_cond_Tvedt, self.p_h)]) / len(self.p_h)
#
#                
#        if mc == "MC":
#            #initialize list for conditional probability and amount of iterations               
#            self.p_cond_MC = []
#            self.iteration_MC = []
#            
#            #create combined event we want to estimate the probability
#            self.distribution_MC = ot.ComposedDistribution(self.marginals)
#            self.vect_MC = ot.RandomVector(self.distribution_MC)
#            self.G_MC = ot.CompositeRandomVector(self.Z_function, self.vect_MC)
#            self.event_MC = ot.Event(self.G_MC, ot.Less(), 0) 
#            
#            #asks settings input for MC algorithm
#            self.default_MC = input('Would you like to use default  settings for MC? (y/n): ') or 0
#            if self.default_MC == "n":
#                self.n_MC = input('How many iterations for MC: ') or 100000 #default 100000
#                self.n_MC = int(self.n_MC)
#                self.blocksize_MC = input('Set block size for MC: ') or 50 #default 50
#                self.blocksize_MC = int(self.blocksize_MC)
#                self.CoV_MC = input('Set the maximum for the coeffiecient of variation of MC: ') or 0.05 #default 0.05
#                self.CoV_MC = float(self.CoV_MC)
#            else:
#                self.n_MC = 100000 #default 100000
#                self.blocksize_MC = 50 #default 50
#                self.CoV_MC = 0.05 #default 0.05
#            
#            #create Monte Carlo algorithm for combined mechanism
#            self.experiment_MC = ot.MonteCarloExperiment()
#            self.algo_MC = ot.ProbabilitySimulationAlgorithm(self.event_MC, self.experiment_MC)
#            self.algo_MC.setMaximumCoefficientOfVariation(self.CoV_MC)
#            self.algo_MC.setMaximumOuterSampling(self.n_MC)
#            self.algo_MC.setBlockSize(self.blocksize_MC)
#            
#            for h in self.h_hydra:
#                self.h = h
#                self.initialNumberOfCall_MC = self.Z_function.getEvaluationCallsNumber()
#                
#                self.algo_MC.run()
#            
#                self.result_MC = self.algo_MC.getResult()
#                self.p_cond_MC.append(self.result_MC.getProbabilityEstimate())
#                self.iteration_MC.append(self.Z_function.getEvaluationCallsNumber() - self.initialNumberOfCall_MC)
#        
#            #probability failure
#            self.p_f_MC = sum([a * b for a, b in zip(self.p_cond_MC, self.p_h)]) / len(self.p_h)
#        
#        if mc_l == "MC_L":
#            #initialize list for conditional probability and amount of iterations               
#            self.p_cond_MC_L = []
#            self.iteration_MC_L = []
#            
#            #create combined event we want to estimate the probability
#            self.distribution_MC_L = ot.ComposedDistribution(self.marginals)
#            self.vect_MC_L = ot.RandomVector(self.distribution_MC)
#            self.G_MC_L = ot.CompositeRandomVector(self.Z_MC_function, self.vect_MC)
#            self.event_MC_L = ot.Event(self.G_MC, ot.Equal(), 0) 
#            
#            #asks settings input for MC algorithm
#            self.default_MC_L = input('Would you like to use default  settings for MC? (y/n): ') or 0
#            if self.default_MC_L == "n":
#                self.n_MC_L = input('How many iterations for MC: ') or 100000 #default 100000
#                self.n_MC_L = int(self.n_MC_L)
#                self.blocksize_MC_L = input('Set block size for MC: ') or 50 #default 50
#                self.blocksize_MC_L = int(self.blocksize_MC_L)
#                self.CoV_MC_L = input('Set the maximum for the coeffiecient of variation of MC: ') or 0.05 #default 0.05
#                self.CoV_MC_L = float(self.CoV_MC_L)
#            else:
#                self.n_MC_L = 100000 #default 100000
#                self.blocksize_MC_L = 50 #default 50
#                self.CoV_MC_L = 0.05 #default 0.05
#            
#            #create Monte Carlo algorithm for combined mechanism
#            self.experiment_MC_L = ot.MonteCarloExperiment()
#            self.algo_MC_L = ot.ProbabilitySimulationAlgorithm(self.event_MC_L, self.experiment_MC_L)
#            self.algo_MC_L.setMaximumCoefficientOfVariation(self.CoV_MC_L)
#            self.algo_MC_L.setMaximumOuterSampling(self.n_MC_L)
#            self.algo_MC_L.setBlockSize(self.blocksize_MC_L)
#            
#            for h in self.h_hydra:
#                self.h = h
#                self.initialNumberOfCall_MC_L = self.Z_MC_function.getEvaluationCallsNumber()
#                
#                self.algo_MC_L.run()
#            
#                self.result_MC_L = self.algo_MC_L.getResult()
#                self.p_cond_MC_L.append(self.result_MC_L.getProbabilityEstimate())
#                self.iteration_MC_L.append(self.Z_MC_L_function.getEvaluationCallsNumber() - self.initialNumberOfCall_MC_L)
#        
#            #probability failure
#            self.p_f_MC_L = sum([a * b for a, b in zip(self.p_cond_MC_L, self.p_h)]) / len(self.p_h)

#class set_event(dijk):
#    def __init__(self, vak_id, scenario, key, name):
#        """
#        This function creates an event for the openTURNS reliability methods (e.g. FORM, SORM, Monte Carlo)
#        
#        ---Parameters---
#        key, should be a list with the name of the variables used in the Sellmeijer model (D, D_cover, d70, h_exit, i_ch, k, L, m_p, m_u)
#        name, should be a string with the corresponding failure mechanism ('uplift', 'heave', 'piping', 'all')
#        
#        """
#        super().__init__(vak_id, scenario)
#        
#        #create joint distribution of parameters
#        lognormdist = ot.LogNormalMuSigma
#        self.distribution_D = lognormdist(self.D,self.D_s,0).getDistribution()
#        self.distribution_D_cover = lognormdist(self.D_cover,self.D_cover_s).getDistribution()
#        self.distribution_d70 = lognormdist(self.d70,self.d70_s,0).getDistribution()
#        self.distribution_h_exit = ot.Normal(self.h_exit,self.h_exit_s)
#        self.distribution_i_ch = lognormdist(self.i_ch,self.i_ch_s,0).getDistribution()
#        self.distribution_k = lognormdist(self.k,self.k_s,0).getDistribution()
#        self.distribution_L = lognormdist(self.L,self.L_s,0).getDistribution()
#        self.distribution_m_p = lognormdist(self.m_p,self.m_p_s,0).getDistribution()
#        self.distribution_m_u = lognormdist(self.m_u,self.m_u_s,0).getDistribution()
#        self.marginals_u = [self.distribution_D_cover, self.distribution_h_exit, self.distribution_m_u]
#        self.marginals_h = [self.distribution_D_cover, self.distribution_h_exit, self.distribution_i_ch]
#        self.marginals_p = [self.distribution_D, self.distribution_D_cover,
#                            self.distribution_d70, self.distribution_h_exit,
#                            self.distribution_k, self.distribution_L, 
#                            self.distribution_m_p]
#        
#        #dictionary
#        self.distribution_dict = {'D' : self.distribution_D,
#                             'D_cover' : self.distribution_D_cover,
#                             'd70' : self.distribution_d70,
#                             'h_exit' : self.distribution_h_exit,
#                             'i_ch' : self.distribution_i_ch,
#                             'k' : self.distribution_k,
#                             'L' : self.distribution_L,
#                             'm_p' : self.distribution_m_p,
#                             'm_u' : self.distribution_m_u
#                            }
#        
#        self.function_dict = {'uplift', self.Z_u_function,
#                              'heave', self.Z_p_function,
#                              'piping', self.Z_h_function,
#                              'all', self.Z_function 
#                             }
#                
#        self.marginals = []
#        for i in key:
#            self.marginals.append(self.distribution_dict[i])
#            
#        self.RS = ot.CorrelationMatrix(len(self.marginals))
#        self.R = ot.NormalCopula.GetCorrelationFromSpearmanCorrelation(self.RS)
#        self.copula = ot.NormalCopula(self.R)    
#            
#        self.distribution = ot.ComposedDistribution(self.marginals)
#        self.distribution.setDescription(key)
#        
#        self.vect = ot.RandomVector(self.distribution)
#        self.G = ot.CompositeRandomVector(self.function_dict[name], self.vect)
#        self.event = ot.Event(self.G, ot.Less(), 0.0)
#        self.event.setName('{} failure'.format(name))
        
##%%
#        
#class dijk_wbi(set_event):
#    def __init__(self, vak_id, scenario, key, name):
#        super().__init__(vak_id, scenario, key, name)
#    
#    def run_wbi(self, ds, form, mc, mc_l, sorm):
#        """
#        This function provides the following algorithms used from openTURNS:
#            - Directional Sampling (DS)
#            - First Order Reliability Method (FORM)
#            - Monte Carlo Method (MC)
#            - Second Order Reliability Method (SORM)
#            
#        The function run(ds, form, mc, sorm) has the following arguments
#            ds, if assigned value "DS" then Directional Sampling is used, otherwise it is not used.
#            form, if assigned "FORM" then FORM is used, otherwise it is not used.
#            mc, if assigned "MC" then Monte Carlo method is used, otherwise it is not used.
#            sorm, if assigned "SORM" then SORM is used, otherwise it is not used.
#        """ 
#        
#        #create joint distribution of parameters
#        lognormdist = ot.LogNormalMuSigma
#        self.distribution_D = lognormdist(self.D,self.D_s,0).getDistribution()
#        self.distribution_D_cover = lognormdist(self.D_cover,self.D_cover_s).getDistribution()
#        self.distribution_d70 = lognormdist(self.d70,self.d70_s,0).getDistribution()
#        self.distribution_h_exit = ot.Normal(self.h_exit,self.h_exit_s)
#        self.distribution_i_ch = lognormdist(self.i_ch,self.i_ch_s,0).getDistribution()
#        self.distribution_k = lognormdist(self.k,self.k_s,0).getDistribution()
#        self.distribution_L = lognormdist(self.L,self.L_s,0).getDistribution()
#        self.distribution_m_p = lognormdist(self.m_p,self.m_p_s,0).getDistribution()
#        self.distribution_m_u = lognormdist(self.m_u,self.m_u_s,0).getDistribution()
#        self.marginals_u = [self.distribution_D_cover, self.distribution_h_exit, self.distribution_m_u]
#        self.marginals_h = [self.distribution_D_cover, self.distribution_h_exit, self.distribution_i_ch]
#        self.marginals_p = [self.distribution_D, self.distribution_D_cover,
#                            self.distribution_d70, self.distribution_h_exit,
#                            self.distribution_k, self.distribution_L, 
#                            self.distribution_m_p]
#                
#        if ds == "DS" or form == "FORM" or sorm == "SORM":
#            #create copula uplift
#            self.RS_u = ot.CorrelationMatrix(len(self.marginals_u))
#            self.R_u = ot.NormalCopula.GetCorrelationFromSpearmanCorrelation(self.RS_u)
#            self.copula_u = ot.NormalCopula(self.R_u)
#
#            #create copula heave
#            self.RS_h = ot.CorrelationMatrix(len(self.marginals_h))
#            self.R_h = ot.NormalCopula.GetCorrelationFromSpearmanCorrelation(self.RS_h)
#            self.copula_h = ot.NormalCopula(self.R_h)
#
#            #create copula piping
#            self.RS_p = ot.CorrelationMatrix(len(self.marginals_p))
#            self.R_p = ot.NormalCopula.GetCorrelationFromSpearmanCorrelation(self.RS_p)
#            self.copula_p = ot.NormalCopula(self.R_p)
#            
#          
#            #marginals uplift
#            self.distribution_u = ot.ComposedDistribution(self.marginals_u)
#            self.distribution_u.setDescription(['h_exit', 'D_cover', 'm_u'])
#            #marginals heave
#            self.distribution_h = ot.ComposedDistribution(self.marginals_h)
#            self.distribution_h.setDescription(['h_exit', 'D_cover', 'i_ch'])
#            #marginals piping
#            self.distribution_p = ot.ComposedDistribution(self.marginals_p)
#            self.distribution_p.setDescription(['h_exit', 'D_cover', 'L', 'D', 'm_p', 'k', 'd70'])
#            
#            #create the event we want to estimate the probability uplift
#            self.vect_u = ot.RandomVector(self.distribution_u)
#            self.G_u = ot.CompositeRandomVector(self.Z_u_function, self.vect_u)
#            self.event_u = ot.Event(self.G_u, ot.Less(), 0.0)
#            self.event_u.setName('uplift failure')
#            
#            #create the event we want to estimate the probability heave
#            self.vect_h = ot.RandomVector(self.distribution_h)
#            self.G_h = ot.CompositeRandomVector(self.Z_h_function, self.vect_h)
#            self.event_h = ot.Event(self.G_h, ot.Less(), 0.0)
#            self.event_h.setName('heave failure')
#
#            #create the event we want to estimate the probability piping
#            self.vect_p = ot.RandomVector(self.distribution_p)
#            self.G_p = ot.CompositeRandomVector(self.Z_p_function, self.vect_p)
#            self.event_p = ot.Event(self.G_p, ot.Less(), 0.0)
#            self.event_p.setName('piping failure')
#            
##            self.uplift = self.event(self.vak_id, self.scenario, ['h_exit', 'D_cover', 'm_u'], 'uplift')
##            self.heave = self.event(self.vak_id, self.scenario, ['h_exit', 'D_cover', 'i_ch'], 'heave')
##            self.piping = self.event(self.vak_id, self.scenario, ['h_exit', 'D_cover', 'L', 'D', 'm_p', 'k', 'd70'], 'piping')
#            
#            #general settings for FORM and SORM
#            if form == "FORM" or sorm == "SORM":                
#                #Asks amount of max iterations
#                self.n = input('Set max amount of iterations for FORM and SORM: ') or 10000
#                self.n = int(self.n)
#            
#            if form == "FORM":
#                #initialize list for conditional probability and amount of iterations               
#                self.p_cond_u_FORM = []
#                self.p_cond_h_FORM = []
#                self.p_cond_p_FORM = []
#                self.iteration_u_FORM = []
#                self.iteration_h_FORM = []
#                self.iteration_p_FORM = []
#                self.Hasofer_u_FORM = []
#                self.Hasofer_h_FORM = []
#                self.Hasofer_p_FORM = []
#                self.standard_design__u_FORM = []
#                self.standard_design__h_FORM = []
#                self.standard_design_p_FORM = []
#                self.physical_design__u_FORM = []
#                self.physical_design_h_FORM = []
#                self.physical_design_p_FORM = []
#            
#                #define a solver
#                self.optimAlgo_FORM = ot.Cobyla()
#                self.optimAlgo_FORM.setMaximumEvaluationNumber(self.n)
#                self.optimAlgo_FORM.setMaximumAbsoluteError(1.0e-10)
#                self.optimAlgo_FORM.setMaximumRelativeError(1.0e-10)
#                self.optimAlgo_FORM.setMaximumResidualError(1.0e-10)
#                self.optimAlgo_FORM.setMaximumConstraintError(1.0e-10)
#                self.algo_FORM = ot.FORM(self.optimAlgo_FORM, self.event, self.distribution.getMean())
#                
#                for h in self.h_hydra:
#                    self.h = h
#                    self.initialNumberOfCall_u = self.Z_u_function.getEvaluationCallsNumber()
#                    self.initialNumberOfCall_h = self.Z_h_function.getEvaluationCallsNumber()
#                    self.initialNumberOfCall_p = self.Z_p_function.getEvaluationCallsNumber()
#                    
#                    #run FORM uplift
#                    self.algo_u.run()
#                    #get result
#                    self.result_u = self.algo_u.getResult()
#                    self.p_cond_u.append(self.result_u.getEventProbability())
#                    self.Hasofer_FORM.append(self.result_FORM.getHasoferReliabilityIndex())
#                    self.standard_design_FORM.append(self.result_FORM.getStandardSpaceDesignPoint())
#                    self.physical_design_FORM.append(self.result_FORM.getPhysicalSpaceDesignPoint())
#                    self.iteration_u.append(self.Z_u_function.getEvaluationCallsNumber() - self.initialNumberOfCall_u)
#                
#                    #run FORM heave
#                    self.algo_h.run()
#                    #get result
#                    self.result_h = self.algo_h.getResult()
#                    self.p_cond_h.append(self.result_h.getEventProbability())
#                    self.Hasofer_FORM.append(self.result_FORM.getHasoferReliabilityIndex())
#                    self.standard_design_FORM.append(self.result_FORM.getStandardSpaceDesignPoint())
#                    self.physical_design_FORM.append(self.result_FORM.getPhysicalSpaceDesignPoint())
#                    self.iteration_h.append(self.Z_h_function.getEvaluationCallsNumber() - self.initialNumberOfCall_h)
#                    
#                    #run FORM piping
#                    self.algo_p.run()
#                    #get result
#                    self.result_p = self.algo_p.getResult()
#                    self.p_cond_p.append(self.result_p.getEventProbability())
#                    self.Hasofer_FORM.append(self.result_FORM.getHasoferReliabilityIndex())
#                    self.standard_design_FORM.append(self.result_FORM.getStandardSpaceDesignPoint())
#                    self.physical_design_FORM.append(self.result_FORM.getPhysicalSpaceDesignPoint())
#                    self.iteration_p.append(self.Z_p_function.getEvaluationCallsNumber() - self.initialNumberOfCall_p)
#                    
#                    #run FORM
#                    self.algo_FORM.run()
#                   
#                    #get result
#                    self.result_FORM = self.algo_FORM.getResult()
#                    self.p_cond_FORM.append(self.result_FORM.getEventProbability())
#                    self.Hasofer_FORM.append(self.result_FORM.getHasoferReliabilityIndex())
#                    self.standard_design_FORM.append(self.result_FORM.getStandardSpaceDesignPoint())
#                    self.physical_design_FORM.append(self.result_FORM.getPhysicalSpaceDesignPoint())
#                    self.iteration_FORM.append(self.Z_function.getEvaluationCallsNumber() - self.initialNumberOfCall_FORM)
#                
#                #probability failure  
#                self.p_cond_FORM = [min(a, b, c) for a, b, c in zip(self.p_cond_u_FORM, self.p_cond_h_FORM, self.p_cond_p_FORM)]
#                self.p_f_FORM = sum([a * b for a, b in zip(self.p_cond_FORM, self.p_h)]) / len(self.p_h)
#                self.Hasofer_FORM = []
#                self.standard_design_FORM = []
#                self.physical_design_FORM = []
#                self.iteration_FORM = [a + b + c for a, b, c in zip(self.iteration_u_FORM, self.iteration_h_FORM, self.iteration_p_FORM)] 
#    