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
    def __init__(self,bounds):
        self._bounds=bounds
    def check_constraints(self, full_results) -> bool:
        return True

    @property
    def n_obj(self) -> int:
        return 2

    def get_objectives(self, valid_constraints, full_results) -> tuple:
        last_results=full_results[-1]
        last_state=last_results[-1]
        r_ro=last_state.design.machine.r_ro
        r_so=last_state.design.machine.r_so
        l_st=last_state.design.machine.l_st
        v_tip=last_state.design.settings.Omega*r_ro
        B_delta=last_state.design.machine.B_delta
        print(B_delta)
        J=last_state.conditions.J
        A_hat=last_state.design.machine.A_hat(J)
        V_r=np.pi*r_ro**2*l_st
        V_s=np.pi*r_so**2*l_st
        Torque=V_r*B_delta*A_hat
        Power=Omega*Torque
        Power_density=Power/V_s
        #print(Power)
        return (-Power,-Power_density)
    
    @property
    def bounds(self) -> tuple:

        return self._bounds
#%%


#
if __name__ == '__main__':
    #Settings
    h=100
    Omega=1001
    #Create Designer
    arch=MotorDesign.arch
    settings_handler=MotorDesign.SPM_SettingsHandler(Omega,h)
    des=me.MachineDesigner(arch,settings_handler)

    #Create evaluation steps
    evalSteps=[steps.RDStep,
               steps.SleeveDesignStep,
               steps.CoreLossStep,
               steps.ThermalStep]#TODO define steps
    #Create Evaluator
    evaluator=me.MachineEvaluator(evalSteps)
    
    #Create DesignSpace
    # r_ro=x[0]
    # d_m_norm=x[1]
    # L/r=x[2]
    # B_sy_nrom=x[3]
    # B_th_norm=x[4]
    # l_tooth_norm=x[5]
    bounds=((.01,0, 0, .1,1,0),
            (.25,1,20,10,10,5))
    design_space=DesignSpace(bounds)
    path = os.path.abspath("")
    arch_file = path + rf"\opti_arch{Omega}.pkl"  # specify path where saved data will reside
    des_file = path + rf"\opti_designer{Omega}.pkl"
    pop_file = path + rf"\latest_pop{Omega}.csv"
    dh=DataHandler(arch_file,des_file)
    # dh=DataHandler()
    
    #Create Machine Design Problem
    machDesProb=mo.DesignProblem(des,evaluator,design_space,dh)
    
    #Run Optimization
    opt=mo.DesignOptimizationMOEAD(machDesProb)
    pop_size=500
    pop=opt.initial_pop(pop_size)
    pop=opt.run_optimization(pop,10)
    archive=list(dh.load_from_archive())[-pop_size:]
    design_list=[None,]*pop_size
    full_results_list=[None,]*pop_size
    final_state_list=[None,]*pop_size
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
    plt.scatter(-objs_list[ndf[0],1],-objs_list[ndf[0],0])
    ax=plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel('Power [W]')
    ax.set_xlabel('Power Density [W/m^3]')
    ax.set_title(f'Omega = {Omega} rad/s')
    
    B_delta_list=[None,]*len(ndf[0])
    V_tip_list=[None,]*len(ndf[0])
    l_st_list=[None,]*len(ndf[0])
    A_hat_list=[None,]*len(ndf[0])
    J_list=[None,]*len(ndf[0])
    Power_list=[None,]*len(ndf[0])
    Power_den_list=[None,]*len(ndf[0])
    A_slot_list=[None,]*len(ndf[0])
    for ind,jnd in enumerate(ndf[0]):
        B_delta_list[ind]=design_list[jnd].machine.B_delta
        V_tip_list[ind]=Omega*design_list[jnd].machine.r_ro
        l_st_list[ind]=design_list[jnd].machine.l_st
        J_list[ind]=full_results_list[jnd][-1][-1].conditions.J
        A_hat_list[ind]=design_list[jnd].machine.A_hat(J_list[ind])
        Power_list[ind]=-objs_list[jnd,0]
        Power_den_list[ind]=-objs_list[jnd,1]
        A_slot_list[ind]=design_list[jnd].machine.A_slot
    fig,axs =plt.subplots(4,1,figsize=(4,4))
    axs[0].scatter(Power_list,B_delta_list,c=Power_den_list)
    axs[0].set_xticks([])
    axs[0].set_xscale('log')
    axs[0].set_xticks([])
    axs[1].scatter(Power_list,V_tip_list,c=Power_den_list)
    axs[1].set_xscale('log')
    axs[1].set_xticks([])
    axs[2].scatter(Power_list,l_st_list,c=Power_den_list)
    axs[2].set_xscale('log')
    axs[2].set_xticks([])
    axs[3].scatter(Power_list,A_hat_list,c=Power_den_list)
    axs[3].set_xscale('log')
    P_max=max(Power_list)
    N=design_list[0].settings.speed
    N*np.sqrt(P_max)
    fig,axs =plt.subplots(6,1,figsize=(4,4))
    axs[0].scatter(Power_list,x_list[ndf[0],0],c=Power_den_list)
    axs[0].set_xticks([])
    axs[0].set_xscale('log')
    axs[0].set_xticks([])
    axs[1].scatter(Power_list,x_list[ndf[0],1],c=Power_den_list)
    axs[1].set_xscale('log')
    axs[1].set_xticks([])
    axs[2].scatter(Power_list,x_list[ndf[0],2],c=Power_den_list)
    axs[2].set_xscale('log')
    axs[2].set_xticks([])
    axs[3].scatter(Power_list,x_list[ndf[0],3],c=Power_den_list)
    axs[3].set_xscale('log')
    axs[3].set_xticks([])
    axs[4].scatter(Power_list,x_list[ndf[0],4],c=Power_den_list)
    axs[4].set_xscale('log')
    axs[4].set_xticks([])
    axs[5].scatter(Power_list,x_list[ndf[0],5],c=Power_den_list)
    axs[5].set_xscale('log')
    axs[5].set_xticks([])
    
    
    fig,axs=plt.subplots(1,1)
    axs.scatter(B_delta_list,A_hat_list,c=A_slot_list)
    fig,axs=plt.subplots(1,1)
    axs.scatter(B_delta_list,A_slot_list,c=A_hat_list)