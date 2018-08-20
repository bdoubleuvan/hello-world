# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 16:41:01 2018

@author: nguyen
"""

from __future__ import print_function
from tqdm import tqdm
import piping_data_deltares as data
import matplotlib.pyplot as plt

#create dike data
vak_id = 4002003#15001015#3003005#15001003
#scenario = 302 #303
#dijk = data.dijk(vak_id, scenario)

#initialize empty list for values in for-loop
p_cond_FORM = {}
p_f_FORM = {}
iteration_FORM = {}
p_cond_FORM_wbi = {}
p_f_FORM_wbi = {}
iteration_FORM_wbi = {}

#initialize empty list for values in for-loop
p_cond_MC = {}
p_f_MC = {}
iteration_MC = {}
p_cond_MC_wbi = {}
p_f_MC_wbi = {}
iteration_MC_wbi = {}


dijk_scenario = data.dijk_scenario(vak_id)
for index in tqdm(dijk_scenario.index):
    #hkv
    dijk_scenario.scenario[index].run(0,"FORM",0,0,0)
    p_cond_FORM[index] = dijk_scenario.scenario[index].p_cond_FORM
    p_f_FORM[index] = dijk_scenario.scenario[index].p_f_FORM
    iteration_FORM[index] = dijk_scenario.scenario[index].iteration_FORM
    #wbi
    dijk_scenario.scenario[index].run_wbi(0,"FORM",0,0,0)
    p_cond_FORM_wbi[index] = dijk_scenario.scenario[index].p_cond_FORM
    p_f_FORM_wbi[index] = dijk_scenario.scenario[index].p_f_FORM
    iteration_FORM_wbi[index] = dijk_scenario.scenario[index].iteration_FORM

for index in tqdm(dijk_scenario.index):
    #hkv
    dijk_scenario.scenario[index].run(0,0,"MC",0,0)
    p_cond_MC[index] = dijk_scenario.scenario[index].p_cond_MC
    p_f_MC[index] = dijk_scenario.scenario[index].p_f_MC
    iteration_MC[index] = dijk_scenario.scenario[index].iteration_MC
    #wbi
    dijk_scenario.scenario[index].run_wbi(0,0,"MC",0,0)
    p_cond_MC_wbi[index] = dijk_scenario.scenario[index].p_cond_MC
    p_f_MC_wbi[index] = dijk_scenario.scenario[index].p_f_MC
    iteration_MC_wbi[index] = dijk_scenario.scenario[index].iteration_MC