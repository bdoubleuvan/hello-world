# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 16:22:27 2018

@author: presentatie
"""

#from __future__ import print_function
import openturns as ot

#create joint distribution of parameters
distribution_R = ot.LogNormalMuSigma(300,30,0).getDistribution()
distribution_F = ot.Normal(75e3,5e3)
marginals = [distribution_R, distribution_F]
distribution = ot.ComposedDistribution(marginals)

#create the model
model = ot.SymbolicFunction(['R', 'F'], ['R-F/(pi_*100)'])

#create the event we want to estimate the probability
vect = ot.RandomVector(distribution)
G = ot.RandomVector(model, vect)
event = ot.Event(G, ot.Less(), 0.0)

initialNumberOfCall = model.getEvaluationCallsNumber()

#create a Monte Carlo algorithm
experiment = ot.MonteCarloExperiment()
algo = ot.ProbabilitySimulationAlgorithm(event, experiment)
algo.setMaximumCoefficientOfVariation(0.05)
algo.setMaximumOuterSampling(99)
algo.setBlockSize(5)
algo.run()

#retrieve results
result = algo.getResult()
probability = result.getProbabilityEstimate()
print('Pf=', probability)
print('#iterations = ', model.getEvaluationCallsNumber() - initialNumberOfCall)