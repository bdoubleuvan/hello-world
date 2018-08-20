# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 11:46:56 2018

@author: nguyen
"""

from __future__ import print_function
import openturns as ot
import piping_data_deltares as data
import matplotlib.pyplot as plt
import pandas as pd

#create dike data
vak_id = 4002003#15001015#3003005#15001003
scenario = 302 #303
dijk = data.dijk(vak_id, scenario)

#temporary implementation of water level and probability
path = r'D:\Users\Nguyen\Desktop\piping_data\hfreq_{}.txt'.format(vak_id)
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

#iterations
n = 10000

#create joint distribution of parameters
lognormdist = ot.LogNormalMuSigma
distribution_D = lognormdist(dijk.D,dijk.D_s,0).getDistribution()
distribution_D_cover = lognormdist(dijk.D_cover,dijk.D_cover_s).getDistribution()
distribution_d70 = lognormdist(dijk.d70,dijk.d70_s,0).getDistribution()
distribution_h_exit = ot.Normal(dijk.h_exit,dijk.h_exit_s)
distribution_i_ch = lognormdist(dijk.i_ch,dijk.i_ch_s,0).getDistribution()
distribution_k = lognormdist(dijk.k,dijk.k_s,0).getDistribution()
distribution_L = lognormdist(dijk.L,dijk.L_s,0).getDistribution()
distribution_m_p = lognormdist(dijk.m_p,dijk.m_p_s,0).getDistribution()
distribution_m_u = lognormdist(dijk.m_u,dijk.m_u_s,0).getDistribution()
#distribution_phi_polder = ot.Normal(dijk.phi_polder,dijk.phi_polder_s)

#marginals uplift
marginals_u = [distribution_D_cover, distribution_h_exit,
               distribution_m_u]
distribution_u = ot.ComposedDistribution(marginals_u)

#marginals heave
marginals_h = [distribution_D_cover, distribution_h_exit,
               distribution_i_ch]
distribution_h = ot.ComposedDistribution(marginals_h)

#marginals piping
marginals_p =[distribution_D, distribution_D_cover,
              distribution_d70, distribution_h_exit,
              distribution_k,
              distribution_L, distribution_m_p]
distribution_p = ot.ComposedDistribution(marginals_p)

#create marginal distribution
marginals = [distribution_D, distribution_D_cover,
             distribution_d70, distribution_h_exit,
             distribution_i_ch, distribution_k,
             distribution_L, distribution_m_p, 
             distribution_m_u]

#create copula uplift
RS_u = ot.CorrelationMatrix(len(marginals_u))
R_u = ot.NormalCopula.GetCorrelationFromSpearmanCorrelation(RS_u)
copula_u = ot.NormalCopula(R_u)

#create copula heave
RS_h = ot.CorrelationMatrix(len(marginals_h))
R_h = ot.NormalCopula.GetCorrelationFromSpearmanCorrelation(RS_h)
copula_h = ot.NormalCopula(R_h)

#create copula piping
RS_p = ot.CorrelationMatrix(len(marginals_p))
R_p = ot.NormalCopula.GetCorrelationFromSpearmanCorrelation(RS_p)
copula_p = ot.NormalCopula(R_p)

#create copula
RS = ot.CorrelationMatrix(len(marginals))
R = ot.NormalCopula.GetCorrelationFromSpearmanCorrelation(RS)
copula = ot.NormalCopula(R)

#create joint probability distribution uplift
distribution_u = ot.ComposedDistribution(marginals_u, copula_u)
distribution_u.setDescription(['h_exit', 'D_cover', 'm_u'])

#create joint probability distribution heave
distribution_h = ot.ComposedDistribution(marginals_h, copula_h)
distribution_h.setDescription(['h_exit', 'D_cover', 'i_ch'])

#create joint probability distribution
distribution_p = ot.ComposedDistribution(marginals_p, copula_p)
distribution_p.setDescription(['h_exit', 'D_cover', 'L', 'D', 'm_p', 'k', 'd70'])

#create joint probability distribution
distribution = ot.ComposedDistribution(marginals, copula)
distribution.setDescription(['h_exit', 'D_cover', 'm_u', 'i_ch', 'L', 'D', 'm_p', 'k', 'd70'])

#create the event we want to estimate the probability uplift
vect_u = ot.RandomVector(distribution_u)
G_u = ot.CompositeRandomVector(dijk.Z_u_function, vect_u)
event_u = ot.Event(G_u, ot.Less(), 0.0)
event_u.setName('uplift failure')

#create the event we want to estimate the probability heave
vect_h = ot.RandomVector(distribution_h)
G_h = ot.CompositeRandomVector(dijk.Z_h_function, vect_h)
event_h = ot.Event(G_h, ot.Less(), 0.0)
event_h.setName('heave failure')

#create the event we want to estimate the probability piping
vect_p = ot.RandomVector(distribution_p)
G_p = ot.CompositeRandomVector(dijk.Z_p_function, vect_p)
event_p = ot.Event(G_p, ot.Less(), 0.0)
event_p.setName('piping failure')

#create the event we want to estimate the probability
vect = ot.RandomVector(distribution)
G = ot.CompositeRandomVector(dijk.Z_function, vect)
event = ot.Event(G, ot.Less(), 0.0)
event.setName('overall failure')

#define a solver
optimAlgo = ot.Cobyla()
optimAlgo.setMaximumEvaluationNumber(n)
optimAlgo.setMaximumAbsoluteError(1.0e-10)
optimAlgo.setMaximumRelativeError(1.0e-10)
optimAlgo.setMaximumResidualError(1.0e-10)
optimAlgo.setMaximumConstraintError(1.0e-10)

#set-up algorithm
algo_u = ot.FORM(optimAlgo, event_u, distribution_u.getMean())
algo_h = ot.FORM(optimAlgo, event_h, distribution_h.getMean())
algo_p = ot.FORM(optimAlgo, event_p, distribution_p.getMean())
algo = ot.FORM(optimAlgo, event, distribution.getMean())

for h in data_h:
    dijk.h = h
            
    initialNumberOfCall_u = dijk.Z_u_function.getEvaluationCallsNumber()
    initialNumberOfCall_h = dijk.Z_h_function.getEvaluationCallsNumber()
    initialNumberOfCall_p = dijk.Z_p_function.getEvaluationCallsNumber()
    initialNumberOfCall = dijk.Z_function.getEvaluationCallsNumber()
    
    #%% FORM
    
    #run FORM uplift
    algo_u.run()
    #get result
    result_u = algo_u.getResult()
    p_cond_u.append(result_u.getEventProbability())
    iteration_u.append(dijk.Z_u_function.getEvaluationCallsNumber() - initialNumberOfCall_u)

    #run FORM heave
    algo_h.run()
    #get result
    result_h = algo_h.getResult()
    p_cond_h.append(result_h.getEventProbability())
    iteration_h.append(dijk.Z_h_function.getEvaluationCallsNumber() - initialNumberOfCall_h)
    
    #run FORM piping
    algo_p.run()
    #get result
    result_p = algo_p.getResult()
    p_cond_p.append(result_p.getEventProbability())
    iteration_p.append(dijk.Z_p_function.getEvaluationCallsNumber() - initialNumberOfCall_p)
    
    #run FORM
    algo.run()
    result = algo.getResult()
    #get result
    p_cond.append(result.getEventProbability())
    iteration.append(dijk.Z_function.getEvaluationCallsNumber() - initialNumberOfCall)
    Hasofer.append(result.getHasoferReliabilityIndex())
    standard_design.append(result.getStandardSpaceDesignPoint())
    physical_design.append(result.getPhysicalSpaceDesignPoint())

#probability failure  
p_cond_min = [min(a, b, c) for a, b, c in zip(p_cond_u, p_cond_h, p_cond_p)]
print('WBI, P_f = ', sum([a * b for a, b in zip(p_cond_min, p_h)])/len(data_h))

#plot fragility curve
plt.figure(0)
plt.plot(data_h, p_cond_min)
plt.xlabel('waterhoogte in [m]')
plt.ylabel('kans')
plt.title('Fragility curve van dijk {}'.format(vak_id)) 

#probability failure
print('P_f = ', sum([a * b for a, b in zip(p_cond, p_h)])/len(data_h))

#plot fragility curve
plt.figure(1)
plt.plot(data_h, p_cond)
plt.xlabel('waterhoogte in [m]')
plt.ylabel('kans')
plt.title('Fragility curve van dijk {}'.format(vak_id))   

gradient = []
for i in physical_design:
    gradient.append(dijk.Z_function.gradient(i))
    
#for i in standard_design:
#    print('Gradient value at standard design point = ',dijk.Z_function.gradient(i))    
#    #result statistic
#    print('Number of iterations = ', dijk.Z_function.getEvaluationCallsNumber() - initialNumberOfCall)
#    print('P_f = ', result.getEventProbability())
#    print('Hasofer reliability index = ', result.getHasoferReliabilityIndex())
#    print('Standard space design point = ', result.getStandardSpaceDesignPoint())
#    print('Physical space design point = ', result.getPhysicalSpaceDesignPoint())
#    print('\n\n\n')
#    result.drawImportanceFactors()
    
    
#    #error history
#    optimResult = result.getOptimizationResult()
#    graphErrors = optimResult.drawErrorHistory()
#    graphErrors.setLegendPosition('bottom')
#    graphErrors.setYMargin(0.0)
#    graphErrors
    
    
    
#    #%% SORM
#    
#    # Get additional results with SORM
#    initialNumberOfCall = dijk.Z_function.getEvaluationCallsNumber()
#    algo = ot.SORM(optimAlgo, event, distribution.getMean())
#    algo.run()
#    sorm_result = algo.getResult()
#    
#    print('Number of iterations = ', dijk.Z_function.getEvaluationCallsNumber() - initialNumberOfCall)
#    print('Breitung P_f = ', sorm_result.getEventProbabilityBreitung())
#    print('Breitung reliability index = ', sorm_result.getGeneralisedReliabilityIndexBreitung())
#    print('Hohen Bichler P_f = ', sorm_result.getEventProbabilityHohenBichler())
#    print('Hohen Bichler reliability index = ', sorm_result.getGeneralisedReliabilityIndexHohenBichler())
#    print('Tvedt P_f = ', sorm_result.getEventProbabilityTvedt())
#    print('Tvedt reliability index = ', sorm_result.getGeneralisedReliabilityIndexTvedt())