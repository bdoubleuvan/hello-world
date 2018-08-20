# -*- coding: utf-8 -*-
"""
Created on Fri Aug 10 10:24:00 2018

@author: nguyen
"""

import openturns as ot
import piping_data_deltares as data
import pandas as pd
import matplotlib.pyplot as plt

#create dike data
vak_id = 3003005
scenario = 302
dijk = data.dijk(vak_id, scenario)

#temporary implementation of water level and probability
path = r'D:\Users\Nguyen\Desktop\piping_data\hfreq.txt'
data = pd.read_csv(path, header = None, skiprows = 1, sep = "\s+")
data.columns = ['waterstand', 'kans']    

data_h = data['waterstand']
p_h = data['kans']


#initialize empty list for values in for-loop
p_cond = []
p_cond_u = []
p_cond_h = []
p_cond_p = []
iteration = []
iteration_u = []
iteration_h = []
iteration_p = []
Hasofer = []
standard_design = []
physical_design = []


#Resolution options:
cv = 0.05
n = int(1e5)

#create joint distribution of parameters
lognormdist = ot.LogNormalMuSigma
distribution_h_exit = ot.Normal(dijk.h_exit,dijk.h_exit_s)
#distribution_phi_polder = ot.Normal(dijk.phi_polder,dijk.phi_polder_s)
distribution_D_cover = lognormdist(dijk.D_cover,dijk.D_cover_s).getDistribution()
distribution_L = lognormdist(dijk.L,dijk.L_s,0).getDistribution()
distribution_D = lognormdist(dijk.D,dijk.D_s,0).getDistribution()
distribution_m_u = lognormdist(dijk.m_u,dijk.m_u_s,0).getDistribution()
distribution_m_p = lognormdist(dijk.m_p,dijk.m_p_s,0).getDistribution()
distribution_k = lognormdist(dijk.k,dijk.k_s,0).getDistribution()
distribution_d70 = lognormdist(dijk.d70,dijk.d70_s,0).getDistribution()
distribution_i_ch = lognormdist(dijk.i_ch,dijk.i_ch_s,0).getDistribution()

#create marginal distribution
marginals = [distribution_h_exit,
             distribution_D_cover, distribution_m_u,
             distribution_i_ch, distribution_L, 
             distribution_D, distribution_m_p, 
             distribution_k, distribution_d70]

#create copula
RS = ot.CorrelationMatrix(len(marginals))
R = ot.NormalCopula.GetCorrelationFromSpearmanCorrelation(RS)
copula = ot.NormalCopula(R)

#create joint probability distribution
distribution = ot.ComposedDistribution(marginals, copula)
distribution.setDescription(['h_exit','D_cover','m_u','i_ch','L','D','m_p','k','d70'])

#create the event we want to estimate the probability
vect = ot.RandomVector(distribution)
G = ot.CompositeRandomVector(dijk.Z_function, vect)
event = ot.Event(G, ot.Less(), 0.0)
event.setName('piping failure')

#root finding algorithm
solver = ot.Bisection()
rootStrategy = ot.SafeAndSlow(solver)

#directional sampling algorithm
samplingStrategy = ot.RandomDirection()
    
for h in data_h:
    dijk.h = h
    initialNumberOfCall = dijk.Z_function.getEvaluationCallsNumber()
        
    #Using Directional sampling
    algoDS = ot.DirectionalSampling(event, rootStrategy, samplingStrategy)
    algoDS.setMaximumOuterSampling(n)
    algoDS.setBlockSize(1)
    algoDS.setMaximumCoefficientOfVariation(cv)
    #For statistics about the algorithm
    initialNumberOfCall = dijk.Z_function.getEvaluationCallsNumber()
    
    #Perform the analysis
    algoDS.run()
    
    result = algoDS.getResult()
    p_cond.append(result.getProbabilityEstimate())
    iteration.append(dijk.Z_function.getEvaluationCallsNumber() - initialNumberOfCall)

#probability failure
print('P_f = ', sum([a * b for a, b in zip(p_cond, p_h)])/len(data_h))

#plot fragility curve
plt.figure(0)
plt.plot(data_h,p_cond)
plt.xlabel('waterhoogte in [m]')
plt.ylabel('kans')
plt.title('Fragility curve van dijk {}'.format(vak_id))  

##Results:
#result = algoDS.getResult()
#probability = result.getProbabilityEstimate()
#print('Number of executed iterations = ', result.getOuterSampling())
#print('Number of calls to the limit state = ', dijk.Z_function.getEvaluationCallsNumber())
#print('Pf = ', probability)
#print('CV = ', result.getCoefficientOfVariation())
#print('\n')
#algoDS.drawProbabilityConvergence()