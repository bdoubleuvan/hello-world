# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 11:13:47 2018

@author: presentatie
"""

from __future__ import print_function
import openturns as ot
import piping_data_deltares as data
import pandas as pd
import matplotlib.pyplot as plt

#create dike data
vak_id = 15001003
scenario = 303
dijk = data.dijk(vak_id, scenario)

#iterations
n = 100000

#temporary implementation of water level and probability
path = r'D:\Users\Nguyen\Desktop\piping_data\hfreq.txt'
data = pd.read_csv(path, header = None, skiprows = 1, sep = "\s+")
data.columns = ['waterstand', 'kans']    

data_h = data['waterstand']
p_h = data['kans']

#initialize empty list for values in for-loop
p_cond = []
iteration = []

for h in data_h:
    dijk.h = h
    
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
    
    ##draw PDF
    #distribution_h_exit.drawPDF()
    #distribution_phi_polder.drawPDF()
    #distribution_D_cover.drawPDF()
    #distribution_L.drawPDF()
    #distribution_D.drawPDF()
    #distribution_m_u.drawPDF()
    #distribution_m_p.drawPDF()
    #distribution_k.drawPDF()
    #distribution_d70.drawPDF()
    #distribution_i_ch.drawPDF()
    
#    #marginals uplift
#    marginals_u = [distribution_h_exit,
#                   distribution_D_cover, distribution_m_u]
#    distribution_u = ot.ComposedDistribution(marginals_u)
#    
#    #marginals heave
#    marginals_h = [distribution_h_exit, distribution_D_cover,
#                   distribution_i_ch]
#    distribution_h = ot.ComposedDistribution(marginals_h)
#    
#    #marginals piping
#    marginals_p = [distribution_h_exit, distribution_D_cover,
#                   distribution_L, distribution_D, 
#                   distribution_m_p, distribution_k, 
#                   distribution_d70]
#    distribution_p = ot.ComposedDistribution(marginals_p)
    
    #marginals uplift, heave and piping combined
    marginals = [distribution_h_exit,
                 distribution_D_cover, distribution_m_u,
                 distribution_i_ch, distribution_L, 
                 distribution_D, distribution_m_p, 
                 distribution_k, distribution_d70]
    distribution = ot.ComposedDistribution(marginals)
    
    ##Analytical model definition uplift
    #Z_u = ot.SymbolicFunction(['h_exit','phi_polder','D_cover','m_u'],
    #                          ['m_u*(D_cover*({}-{})/{})-(phi_polder+{}*({}-phi_polder)-h_exit)'.format(gamma_sat,gamma_water,gamma_water,r_exit,h)])
    #
    ##Analytical model definition heave
    #Z_h = ot.SymbolicFunction(['h_exit','D_cover','i_ch'],
    #                          ['i_ch-{}*({}-h_exit)/D_cover'.format(r_exit,h)])
    #
    ##Analytical model definition piping
    #Z_p = ot.SymbolicFunction(['h_exit','D_cover','L','D','m_p','k','d70'],
    #                          ['m_p*L*{}*{}/{}*tan({}*pi_/180)*{}*({}*k*L/{})^(1/3)*(d70/{})^0.4*0.91*(D/L)^(0.28/((D/L)^2.8-1)+0.04)-({}-h_exit-{}*D_cover)'.format(eta,gamma_sp,gamma_water,theta,d70m,v_water,g,d70m,h,r_c)]) 
    
#    #create the uplift event we want to estimate the probability
#    vect_u = ot.RandomVector(distribution_u)
#    G_u = ot.RandomVector(dijk.Z_u_function, vect_u)
#    #event_u = ot.Event(G_u,ot.Less(),0)
#    
#    #create the heave event we want to estimate the probability
#    vect_h = ot.RandomVector(distribution_h)
#    G_h = ot.RandomVector(dijk.Z_h_function, vect_h)
#    #event_h = ot.Event(G_h,ot.Less(),0) 
#    
#    #create the piping event we want to estimate the probability
#    vect_p = ot.RandomVector(distribution_p)
#    G_p = ot.RandomVector(dijk.Z_p_function, vect_p)
#    
    #create combined event we want to estimate the probability
    vect = ot.RandomVector(distribution)
    G = ot.CompositeRandomVector(dijk.Z_function, vect)
    
#    event_u = ot.Event(G_u, ot.Less(), 0) 
#    event_h = ot.Event(G_h, ot.Less(), 0)
#    event_p = ot.Event(G_u, ot.Less(), 0) 
    event = ot.Event(G, ot.Less(), 0) 
    #initialNumberOfCall_u = functions.Z_u_function.getEvaluationCallsNumber()
    #initialNumberOfCall_h = functions.Z_h_function.getEvaluationCallsNumber()
    #initialNumberOfCall_p = functions.Z_p_function.getEvaluationCallsNumber()
    initialNumberOfCall = dijk.Z_function.getEvaluationCallsNumber()
    
    ##create Monte Carlo algorithm for uplift
    #experiment_u = ot.MonteCarloExperiment()
    #algo_u = ot.ProbabilitySimulationAlgorithm(event_u, experiment_u)
    #algo_u.setMaximumCoefficientOfVariation(0.05)
    #algo_u.setMaximumOuterSampling(n)
    #algo_u.setBlockSize(1)
    #algo_u.run()
    #
    ##create Monte Carlo algorithm for heave
    #experiment_h = ot.MonteCarloExperiment()
    #algo_h = ot.ProbabilitySimulationAlgorithm(event_h, experiment_h)
    #algo_h.setMaximumCoefficientOfVariation(0.05)
    #algo_h.setMaximumOuterSampling(n)
    #algo_h.setBlockSize(1)
    #algo_h.run()
    #
    ###create Monte Carlo algorithm for piping
    #experiment_p = ot.MonteCarloExperiment()
    #algo_p = ot.ProbabilitySimulationAlgorithm(event_p, experiment_p)
    #algo_p.setMaximumCoefficientOfVariation(0.05)
    #algo_p.setMaximumOuterSampling(n)
    #algo_p.setBlockSize(1)
    #algo_p.run()
    
    ##retrieve results for uplift
    #result_u = algo_u.getResult()
    #probability_u = result_u.getProbabilityEstimate()
    #print('P_failure_u = ', probability_u)
    #print('Number of iterations = ', functions.Z_u_function.getEvaluationCallsNumber() - initialNumberOfCall_u)
    #algo_u.drawProbabilityConvergence()
    #
    ##retrieve results for heave
    #result_h = algo_h.getResult()
    #probability_h = result_h.getProbabilityEstimate()
    #print('P_failure_h = ', probability_h)
    #print('Number of iterations = ', functions.Z_h_function.getEvaluationCallsNumber() - initialNumberOfCall_h)
    #algo_h.drawProbabilityConvergence()
    #
    ###retrieve results for piping
    #result_p = algo_p.getResult()
    #probability_p = result_p.getProbabilityEstimate()
    #print('P_failure_ = ', probability_p)
    #print('Number of iterations = ', functions.Z_p_function.getEvaluationCallsNumber() - initialNumberOfCall_p)
    #algo_p.drawProbabilityConvergence()
    
    ##create Monte Carlo algorithm for combined mechanism
    experiment = ot.MonteCarloExperiment()
    algo = ot.ProbabilitySimulationAlgorithm(event, experiment)
    algo.setMaximumCoefficientOfVariation(0.05)
    algo.setMaximumOuterSampling(n)
    algo.setBlockSize(1)
    algo.run()
    
    result = algo.getResult()
    p_cond.append(result.getProbabilityEstimate())

###retrieve results for combined mechanism
#result = algo.getResult()
#probability = result.getProbabilityEstimate()
#print('P_failure_ = ', probability)
#print('Number of iterations = ', dijk.Z_function.getEvaluationCallsNumber() - initialNumberOfCall)
#algo.drawProbabilityConvergence()
    
#probability failure
print('P_f = ', sum([a * b for a, b in zip(p_cond, p_h)])/len(data_h))

#plot fragility curve
plt.figure(0)
plt.plot(data_h,p_cond)
plt.xlabel('waterhoogte in [m]')
plt.ylabel('kans')
plt.title('Fragility curve van dijk {}'.format(vak_id))  