# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 17:30:02 2018

@author: presentatie
"""

from __future__ import print_function
import openturns as ot

#dimension
dim = 2

#Analytical model definition
limitState = ot.SymbolicFunction(['R', 'F'], ['R-F/(pi_*100)'])

#Test of the limit state function
x = [300, 75000]
print('x=', x)
print('G(x)=', limitState(x))

#Stochastic model definition

#Create a first marginal : LogNormal distribution 1D, parametrized by
#its mean and standard deviation
R_dist = ot.LogNormalMuSigma(300, 30, 0).getDistribution()
R_dist.setName('Yield strength')
R_dist.setDescription('R')

#Graphical output of the PDF
R_dist.drawPDF()

#Create a second marginal : Normal distribution 1D
F_dist = ot.Normal(75000, 5000)
F_dist.setName('Traction_load')
F_dist.setDescription('F')

#Graphical output of the PDF
F_dist.drawPDF()

#Create a copula : IndependentCopula (no correlation)
aCopula = ot.IndependentCopula(dim)
aCopula.setName('Independent copula')

#Instanciate one distribution object
myDistribution = ot.ComposedDistribution([R_dist, F_dist], aCopula)
myDistribution.setName('myDist')

#We create a 'usual' RandomVector from the Distribution
vect = ot.RandomVector(myDistribution)

#We create a composite random vector
G = ot.RandomVector(limitState, vect)

#We create an Event from this RandomVector
myEvent = ot.Event(G, ot.Less(), 0.0)

#Using Monte Carlo simulations
cv = 0.05
NbSim = 100000

experiment = ot.MonteCarloExperiment()
algoMC = ot.ProbabilitySimulationAlgorithm(myEvent, experiment)
algoMC.setMaximumOuterSampling(NbSim)
algoMC.setBlockSize(1)
algoMC.setMaximumCoefficientOfVariation(cv)

#For statistics about the algorithm
initialNumberOfCall = limitState.getEvaluationCallsNumber()

#Perform the analysis
algoMC.run()

#Results:
result = algoMC.getResult()
probability = result.getProbabilityEstimate()
print('Monte Carlo result = ', result)
print('Number of executed iterations = ', result.getOuterSampling())
print('Number of calls to the limit state = ', limitState.getEvaluationCallsNumber() - initialNumberOfCall)
print('Pf = ', probability)
print('CV = ', result.getCoefficientOfVariation())
print('\n')
algoMC.drawProbabilityConvergence()



#Using FORM analysis

#We create a NearestPoint algorithm
myCobyla = ot.Cobyla()
#Resolution options
eps = 1e-3
myCobyla.setMaximumEvaluationNumber(100)
myCobyla.setMaximumAbsoluteError(eps)
myCobyla.setMaximumRelativeError(eps)
myCobyla.setMaximumResidualError(eps)
myCobyla.setMaximumConstraintError(eps)

#For statistics about the algorithm
initialNumberOfCall = limitState.getEvaluationCallsNumber()

#We create a FORM algorithm
# The first parameter is a NearestPointAlgorithm
#The second parameter is an event
#The third parameter is a starting point for the design point research

algoFORM = ot.FORM(myCobyla, myEvent, myDistribution.getMean())
#Perform the analysis
algoFORM.run()

#Results:
result = algoFORM.getResult()
print('Number of calls to the limite state = ', limitState.getEvaluationCallsNumber() - initialNumberOfCall)
print('Pf = ', result.getEventProbability())
print('\n')

#Graphical result output
result.drawImportanceFactors()



#Using Directional sampling

#Resolution options:
cv = 0.05
NbSim = int(1e5)

algoDS = ot.DirectionalSampling(myEvent)
algoDS.setMaximumOuterSampling(NbSim)
algoDS.setBlockSize(1)
algoDS.setMaximumCoefficientOfVariation(cv)
#For statistics about the algorithm
initialNumberOfCall = limitState.getEvaluationCallsNumber()

#Perform the analysis
algoDS.run()

#Results:
result = algoDS.getResult()
probability = result.getProbabilityEstimate()
print('Number of executed iterations = ', result.getOuterSampling())
print('Number of calls to the limit state = ', limitState.getEvaluationCallsNumber())
print('Pf = ', probability)
print('CV = ', result.getCoefficientOfVariation())
print('\n')
algoDS.drawProbabilityConvergence()
