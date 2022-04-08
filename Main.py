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
        Omega=last_state.design.settings.Omega
        #v_tip=last_state.design.settings.Omega*r_ro
        B_delta=last_state.design.machine.B_delta
        #J=last_state.design.machine.J_hat
        A_hat=last_state.design.machine.A_hat_eff
        V_r=np.pi*(r_ro**2)*l_st
        V_s=np.pi*(r_so**2)*l_st
        Torque=V_r*B_delta*A_hat
        Power=Omega*Torque
        Power_density=Power/V_s
        return (-np.log10(Power),-np.log10(Omega))
    
    @property
    def bounds(self) -> tuple:

        return self._bounds
    
    @bounds.setter 
    def bounds(self,bounds_new):
        self._bounds=bounds_new
#%%


#
if __name__ == '__main__':
    #Settings
    h=10
    Omega=10000*2*np.pi/60
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
    # r_sy_norm=x[3]
    # r_so_norm=x[4]
    # w_tooth_norm=x[5]
    # J_hat=x[6] Note this is in A/mm^2
    # d_ag =x[7]
    # Omega =x[8] 1000-100000
    bounds=((.003, 0   , 0  ,  0.01 , 1.05  , 0.05    , 0  , .002, 1000),
            (.5 , 1   , 6 , .9    , 10    , .85     , 10 , 0.5, 1000000))
    design_space=DesignSpace(bounds)
    path = os.path.abspath("")
    arch_file = path + rf"\opti_arch_h_{h}_baseline.pkl"  # specify path where saved data will reside
    des_file = path + rf"\opti_designer_h_{h}_baseline.pkl"
    pop_file = path + rf"\latest_pop_h_{h}_baseline.csv"
    dh=DataHandler(arch_file,des_file)
    # dh=DataHandler()
    
    #Create Machine Design Problem
    machDesProb=mo.DesignProblem(des,evaluator,design_space,dh)
    
    #Run Optimization
    opt=mo.DesignOptimizationMOEAD(machDesProb)
    #opt=mo.DesignOptimizationSGA(machDesProb)
    pop_size=500
    pop=opt.initial_pop(pop_size)
    

    pop=opt.run_optimization(pop,500)
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
    Power_front=objs_list[ndf[0],0]
    Omega_front=objs_list[ndf[0],1]
    temp=zip(Omega_front,Power_front)
    sorted_pairs = sorted(temp)
    tuples = zip(*sorted_pairs)
    Omega_sorted, Power_sorted = [ list(tuple) for tuple in  tuples]
    Omega_sorted=np.array(Omega_sorted)
    Power_sorted=np.array(Power_sorted)
    Power_full_list=-objs_list[:,0]
    Power_density_full_list=-objs_list[:,0]/objs_list[:,1]
    #plt.scatter(-objs_list[ndf[0],0]/objs_list[ndf[0],1],-objs_list[ndf[0],0])
    #plt.scatter(60*10**(-objs_list[:,1])/(2*np.pi),10**(-objs_list[:,0]))
    #plt.scatter(60*10**(-objs_list[ndf[0],1])/(2*np.pi),10**(-objs_list[ndf[0],0]),marker='o')
    plt.plot(60*10**(-Omega_sorted)/(2*np.pi),10**(-Power_sorted),marker='.',color='r')
    N_vect=np.logspace(4,6,100)
    plt.plot(N_vect,(8E5/N_vect)**2*1000,'k')
    #plt.scatter(objs_list[ndf[0],1],-objs_list[ndf[0],0])
    ax=plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel('Power [W]')
    ax.set_xlabel('Speed [RPM]')
    ax.set_xlim([10**4,10**6])
    #ax.set_title(f'Omega = {Omega} rad/s')
    # obj_list2=np.array([-Power_density_full_list,-Power_full_list]).T
    # ndf, dl, dc, ndr = pg.fast_non_dominated_sorting(obj_list2) 
    plt.scatter(Power_density_full_list,Power_full_list)
    #plt.scatter(Power_density_full_list[ndf[0]],Power_full_list[ndf[0]])
    ax=plt.gca()
    # ax.set_yscale('log')
    # ax.set_xscale('log')
    ax.set_ylabel('Power [W]')
    ax.set_xlabel('Power Density [W/m^3]')
    ax.set_title(f'Omega = {Omega} rad/s')
    
    B_delta_list=[None,]*len(ndf[0])
    V_tip_list=[None,]*len(ndf[0])
    l_st_list=[None,]*len(ndf[0])
    A_hat_list=[None,]*len(ndf[0])
    A_hat_eff_list=[None,]*len(ndf[0])
    J_list=[None,]*len(ndf[0])
    Power_list=[None,]*len(ndf[0])
    Power_den_list=[None,]*len(ndf[0])
    A_slot_list=[None,]*len(ndf[0])
    L_tooth_list=[None,]*len(ndf[0])
    r_ro_list=[None,]*len(ndf[0])
    d_m_list=[None,]*len(ndf[0])
    B_sy_list=[None,]*len(ndf[0])
    B_sy_eff_list=[None,]*len(ndf[0])
    B_th_list=[None,]*len(ndf[0])
    B_th_eff_list=[None,]*len(ndf[0])
    Omega_list=[None,]*len(ndf[0])
    for ind,jnd in enumerate(ndf[0]):
        #jnd=ind
        Omega=10**(-objs_list[jnd,1])
        B_delta_list[ind]=design_list[jnd].machine.B_delta
        Omega_list[ind]=60*10**(-objs_list[jnd,1])/(2*np.pi)
        V_tip_list[ind]=Omega*design_list[jnd].machine.r_ro
        l_st_list[ind]=design_list[jnd].machine.l_st
        J_list[ind]=design_list[jnd].machine.J_hat
        A_hat_list[ind]=design_list[jnd].machine.A_hat
        A_hat_eff_list[ind]=design_list[jnd].machine.A_hat_eff
        Power_list[ind]=10**(-objs_list[jnd,0])
        Power_den_list[ind]=-objs_list[jnd,1]
        A_slot_list[ind]=design_list[jnd].machine.A_slot
        L_tooth_list[ind]=design_list[jnd].machine.l_tooth
        r_ro_list[ind]=design_list[jnd].machine.r_ro
        d_m_list[ind]=design_list[jnd].machine.d_m
        B_sy_list[ind]=design_list[jnd].machine.B_sy
        B_sy_eff_list[ind]=design_list[jnd].machine.B_sy_eff
        B_th_list[ind]=design_list[jnd].machine.B_th
        B_th_eff_list[ind]=design_list[jnd].machine.B_th_eff
    fig,axs =plt.subplots(4,1,figsize=(4,7))
    c1=axs[0].scatter(Omega_list,B_delta_list,c=Power_list)
    axs[0].set_xticks([])
    axs[0].set_xscale('log')
    axs[0].set_xticks([])
    axs[0].set_ylabel('B_delta')
    axs[0].set_ylim([0,1.3])
    
    axs[1].scatter(Omega_list,V_tip_list,c=Power_list)
    axs[1].set_xscale('log')
    axs[1].set_xticks([])
    axs[1].set_ylabel('V_tip')
    axs[1].set_ylim([150,450])
    
    axs[2].scatter(Omega_list,np.array(l_st_list)/np.array(r_ro_list),c=Power_list)
    axs[2].set_xscale('log')
    axs[2].set_xticks([])
    axs[2].set_ylabel('L_st/r_ro')
    axs[2].set_ylim([0,10])
    
    axs[3].scatter(Omega_list,A_hat_eff_list,c=Power_list)
    axs[3].set_xscale('log')
    axs[3].set_ylabel('A_hat')
    axs[3].set_ylim([0,200E3])
    fig.colorbar(c1, ax=axs, label='Power [W]')
    P_max=max(Power_list)
    N=design_list[0].settings.speed
    N*np.sqrt(P_max/1000)
    fig,axs =plt.subplots(8,1,figsize=(4,10))
    axs[0].scatter(Omega_list,x_list[ndf[0],0],c=Power_den_list)
    axs[0].set_xticks([])
    axs[0].set_xscale('log')
    axs[0].set_xticks([])
    axs[0].set_ylim([bounds[0][0],bounds[1][0]])
    axs[0].set_ylabel('r_ro')
    
    axs[1].scatter(Omega_list,x_list[ndf[0],1],c=Power_den_list)
    axs[1].set_xscale('log')
    axs[1].set_xticks([])
    axs[1].set_ylim([bounds[0][1],bounds[1][1]])
    axs[1].set_ylabel('d_m norm')
    
    axs[2].scatter(Omega_list,x_list[ndf[0],2],c=Power_den_list)
    axs[2].set_xscale('log')
    axs[2].set_xticks([])
    axs[2].set_ylim([bounds[0][2],bounds[1][2]])
    axs[2].set_ylabel('L/r')
    
    axs[3].scatter(Omega_list,x_list[ndf[0],3],c=Power_den_list)
    axs[3].set_xscale('log')
    axs[3].set_xticks([])
    axs[3].set_ylim([bounds[0][3],bounds[1][3]])
    axs[3].set_ylabel('r_sy norm')
    
    axs[4].scatter(Omega_list,x_list[ndf[0],4],c=Power_den_list)
    axs[4].set_xscale('log')
    axs[4].set_xticks([])
    axs[4].set_ylim([bounds[0][4],bounds[1][4]])
    axs[4].set_ylabel('r_so norm')
    
    axs[5].scatter(Omega_list,x_list[ndf[0],5],c=Power_den_list)
    axs[5].set_xscale('log')
    axs[5].set_xticks([])
    axs[5].set_ylim([bounds[0][5],bounds[1][5]])
    axs[5].set_ylabel('w_st norm')
    
    axs[6].scatter(Omega_list,x_list[ndf[0],6],c=Power_den_list)
    axs[6].set_xscale('log')
    axs[6].set_xticks([])
    axs[6].set_ylim([bounds[0][6],bounds[1][6]])
    axs[6].set_ylabel('J_hat')
    
    axs[7].scatter(Omega_list,x_list[ndf[0],7],c=Power_den_list)
    axs[7].set_xscale('log')
    axs[5].set_xticks([])
    axs[7].set_ylim([bounds[0][7],.20])
    axs[7].set_ylabel('d_ag')
    axs[7].set_xlabel('Speed [rad/s]')
    
    fig,axs=plt.subplots(1,1)
    axs.scatter(B_delta_list,A_hat_list,c=A_slot_list)
    fig,axs=plt.subplots(1,1)
    axs.scatter(B_delta_list,A_slot_list,c=A_hat_list)
    
    fig,axs=plt.subplots(1,1)
    V_r_list=np.pi*(np.array(r_ro_list)**2)*np.array(l_st_list)
    axs.scatter(Power_list,np.array(B_delta_list)*np.array(A_hat_list)*V_r_list)