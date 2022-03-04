# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 10:59:45 2021

@author: Martin Johnson
"""
import numpy as np
from matplotlib import pyplot as plt
import sys
import pygmo as pg
sys.path.append("..")
from eMach import mach_opt as mo
from eMach import mach_eval as me
from typing import List,Tuple,Any
from copy import deepcopy
import MotorDesign
import MotorEvaluation as steps


class DataHandler:
    """Parent class for all data handlers"""
    def __init__(self):
        self.archive=[1,]
    
    def save_to_archive(self, x, design, full_results, objs):
        self.archive.append([x,design,full_results,objs])

    def save_designer(self, designer):
        pass

class DesignSpace():
    """Parent class for a optimization DesignSpace classes"""


    def check_constraints(self, full_results) -> bool:
        return True

    @property
    def n_obj(self) -> int:
        return 2

    def get_objectives(self, valid_constraints, full_results) -> tuple:
        last_results=full_results[-1]
        last_state=last_results[-1]
        r_ro=last_state.design.machine.r_ro
        l_st=last_state.design.machine.l_st
        v_tip=last_state.conditions.v_tip_max[0]
        B_delta=last_state.conditions.B_delta
        print(B_delta,last_state.design.machine.d_m,v_tip,
              r_ro)
        A_hat=100
        Omega=v_tip/r_ro
        V_r=np.pi*r_ro**2*l_st
        Torque=V_r*B_delta*A_hat
        Power=Omega*Torque
        return (-Omega,-Power)
    
    @property
    def bounds(self) -> tuple:
        # r_ro=x[1]
        # d_m=x[2]
        # d_ag=x[3]
        # l_tooth=x[4]
        # d_yoke=x[5]
        # k_tooth=x[6]
        bounds=((.01,0,.001,.01,.01,.01),
                (1,1,.01,1,1,1))
        return bounds
#%%

if __name__ == '__main__':
    
    #Create Designer
    des=MotorDesign.SPM_Designer

    #Create evaluation steps
    evalSteps=[steps.StructuralStep,
               steps.RDStep,
               steps.MagStep]#TODO define steps
    #Create Evaluator
    evaluator=me.MachineEvaluator(evalSteps)
    design_space=DesignSpace()
    dh=DataHandler()
    
    #set evaluation bounds
    bounds=([0,0,0],[1,1,1])
    #set number of objectives
    n_obj=3
    
    #Create Machine Design Problem
    machDesProb=mo.DesignProblem(des,evaluator,design_space,dh)
    
    #Run Optimization
    opt=mo.DesignOptimizationMOEAD(machDesProb)
    pop=opt.initial_pop(500)
    pop=opt.run_optimization(pop,10)
    fits, vectors = pop.get_f(), pop.get_x()
    ndf, dl, dc, ndr = pg.fast_non_dominated_sorting(fits) 
    plt.scatter(-fits[ndf[0],0],-fits[ndf[0],1])
    ax=plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')


