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
import os
import sys
import pickle



class DataHandler(mo.DataHandler):
    def __init__(self, archive_filepath, designer_filepath):
        self.archive_filepath = archive_filepath
        self.designer_filepath = designer_filepath

    def save_to_archive(self, x, design, full_results, objs):
        # assign relevant data to OptiData class attributes
        opti_data = mo.OptiData(
            x=x, design=design, full_results=full_results, objs=objs
        )
        # write to pkl file. 'ab' indicates binary append
        with open(self.archive_filepath, "ab") as archive:
            pickle.dump(opti_data, archive, -1)

    def load_from_archive(self):
        with open(self.archive_filepath, "rb") as f:
            while 1:
                try:
                    yield pickle.load(f)  # use generator
                except EOFError:
                    break

    def save_designer(self, designer):
        with open(self.designer_filepath, "wb") as des:
            pickle.dump(designer, des, -1)
# class DataHandler:
#     """Parent class for all data handlers"""
#     def __init__(self):
#         self.archive=[None,]
    
#     def save_to_archive(self, x, design, full_results, objs):
#         self.archive.append([x,design,full_results,objs])
#         print(self.archive)
#     def save_designer(self, designer):
#         pass

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
        v_tip=last_state.conditions.v_tip_max
        B_delta=last_state.conditions.B_delta
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
        bounds=((.001,0,.002,.01,.01,.01),
                (.025,1,.01,1,1,1))
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
    path = os.path.abspath("")
    arch_file = path + r"\opti_arch.pkl"  # specify path where saved data will reside
    des_file = path + r"\opti_designer.pkl"
    pop_file = path + r"\latest_pop.csv"
    dh=DataHandler(arch_file,des_file)
    
    #set evaluation bounds
    bounds=([0,0,0],[1,1,1])
    #set number of objectives
    n_obj=3
    
    #Create Machine Design Problem
    machDesProb=mo.DesignProblem(des,evaluator,design_space,dh)
    
    #Run Optimization
    opt=mo.DesignOptimizationMOEAD(machDesProb)
    pop_size=300
    pop=opt.initial_pop(pop_size)
    pop=opt.run_optimization(pop,40)
    archive=list(dh.load_from_archive())[-pop_size:]
    design_list=[None,]*pop_size
    full_results_list=[None,]*pop_size
    objs_list=[None,]*pop_size
    x_list=[None,]*pop_size
    for ind,data in enumerate(archive):
        design_list[ind]=data.design
        full_results_list[ind]=data.full_results
        objs_list[ind]=data.objs
        x_list[ind]=data.x
    fits, vectors = pop.get_f(), pop.get_x()
    objs_list=np.array(objs_list)
    x_list=np.array(x_list)
    full_results_list=np.array(full_results_list)
    design_list=np.array(design_list)
    ndf, dl, dc, ndr = pg.fast_non_dominated_sorting(objs_list) 
    plt.scatter(-objs_list[ndf[0],0]*60/(2*np.pi),-objs_list[ndf[0],1])
    ax=plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel('Power [W]')
    ax.set_xlabel('Speed [RPM]')
    B_delta_list=[None,]*len(ndf[0])
    for ind,jnd in enumerate(ndf[0]):
        B_delta_list[ind]=full_results_list[jnd,2,1][0]
    fig,axs =plt.subplots(4,1)
    axs[0].scatter(-objs_list[ndf[0],0]*60/(2*np.pi),x_list[ndf[0],0])
    axs[0].set_xticks([])
    axs[1].scatter(-objs_list[ndf[0],0]*60/(2*np.pi),x_list[ndf[0],1])
    axs[1].set_xticks([])
    axs[2].scatter(-objs_list[ndf[0],0]*60/(2*np.pi),full_results_list[ndf[0],0,1])
    axs[2].set_xticks([])
    axs[3].scatter(-objs_list[ndf[0],0]*60/(2*np.pi),B_delta_list)


