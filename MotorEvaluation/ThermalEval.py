# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 12:29:38 2022

@author: Martin Johnson
"""

import numpy as np
from matplotlib import pyplot as plt
import sys
sys.path.append("..")
from eMach import mach_opt as mo
from eMach import mach_eval as me
from typing import List,Tuple,Any
from copy import deepcopy

__all__=['ThermalStep',]

class ThermalProblemDefinition(me.ProblemDefinition):
    """Class converts input state into a problem"""
    
    def get_problem(self,state:'me.State')->'me.Problem':
        """Returns Problem from Input State"""
        #TODO define problem definition
        g_sy=state.conditions.g_sy
        g_th=state.conditions.g_tooth
        w_st=state.design.machine.w_tooth
        l_st=state.design.machine.l_st
        l_tooth=state.design.machine.l_tooth####
        alpha_q=state.design.machine.alpha_q
        r_so=state.design.machine.r_so#####
        r_sy=state.design.machine.r_sy
        
        k_ins=state.design.machine.ins_mat['k']
        w_ins=state.design.machine.w_ins
        k_fe=state.design.machine.core_mat['core_therm_conductivity']
        h=state.design.settings.h
        alpha_slot=state.design.machine.alpha_slot
        T_coil=state.design.machine.coil_mat['Max_temp']
        
        r_si=state.design.machine.r_si
        Q=state.design.machine.Q
        sigma=state.design.machine.coil_mat['sigma']
        k_fill=state.design.machine.coil_mat['k_fill']
        y=state.design.machine.y
        k_ov=state.design.machine.coil_mat['k_ov']
        
        problem=ThermalProblem(g_sy,g_th,w_st,l_st,l_tooth,alpha_q,r_so,r_sy,k_ins,
                 w_ins,k_fe,h,alpha_slot,T_coil,r_si,Q,sigma,k_fill,y,k_ov)
        return problem

class ThermalProblem():
    """problem class utilized by the Analyzer
    
    Attributes:
        TODO
    """
    def __init__(self,g_sy,g_th,w_st,l_st,l_tooth,alpha_q,r_so,r_sy,k_ins,
                 w_ins,k_fe,h,alpha_slot,T_coil,r_si,Q,sigma,k_fill,y,k_ov):
        """Creates problem class
        
        Args:
            TODO
            
        """
        self.g_sy=g_sy
        self.g_th=g_th
        self.w_st=w_st
        self.l_st=l_st
        self.l_tooth=l_tooth####
        self.alpha_q=alpha_q
        self.r_so=r_so#####
        self.r_sy=r_sy
        
        self.k_ins=k_ins
        self.w_ins=w_ins
        self.k_fe=k_fe
        self.h=h
        self.alpha_slot=alpha_slot
        self.T_coil=T_coil
        self.r_si=r_si
        self.Q=Q
        self.sigma=sigma
        self.k_fill=k_fill
        self.y=y
        self.k_ov=k_ov
    
class ThermalAnalyzer(me.Analyzer):
    """"Class Analyzes the CubiodProblem  for volume and Surface Areas"""
    
    def analyze(self,problem:'me.Problem'):
        """Performs Analysis on a problem

        Args:
            problem (me.Problem): Problem Object

        Returns:
            results (Any): 
                Results of Analysis

        """
        g_sy=problem.g_sy
        g_th=problem.g_th
        w_st=problem.w_st
        l_st=problem.l_st
        l_tooth=problem.l_tooth####
        alpha_q=problem.alpha_q
        r_so=problem.r_so#####
        r_sy=problem.r_sy
        
        k_ins=problem.k_ins
        w_ins=problem.w_ins
        k_fe=problem.k_fe
        h=problem.h
        alpha_slot=problem.alpha_slot
        T_coil=problem.T_coil
        r_si=problem.r_si
        Q=problem.Q
        sigma=problem.sigma
        k_fill=problem.k_fill
        y=problem.y
        k_ov=problem.k_ov
        
        l_turn=2*(l_st+y*alpha_q*(r_si+r_sy)*k_ov)
        
        
        # r_vect=np.linspace(r_sy,r_so,100)
        
        Q_tooth=g_th*w_st*l_st*l_tooth/2
        Q_sy=g_sy*alpha_q/2*(r_so**2-r_sy**2)*l_st
        
        zeta=np.sqrt(2*k_ins/(w_st*w_ins*k_fe))

        
        
        M_th=w_st*l_st/(2*zeta)*np.tanh(zeta*l_tooth)
        R_sy=1/(h*r_so*alpha_q*l_st)+np.log(r_so/r_sy)/(k_fe*alpha_q*l_st)
        V_sy=l_st*(alpha_q/2)*(r_so**2-r_sy**2)
        M_sy=r_sy**2/(2*k_fe)*np.log(r_sy/r_so)+(r_so**2-r_sy**2)/(4*k_fe)+V_sy/(h*r_so*alpha_q*l_st)
        
        R_coil_st=w_ins*zeta/(k_ins*l_st*np.tanh(zeta*l_tooth))
        R_coil_sy=w_ins/(k_ins*r_sy*alpha_slot*l_st)
        R_coil=(2/R_coil_st + 1/R_coil_sy)**-1
        
        Q_coil=(T_coil -g_sy*M_sy - 2*Q_tooth*R_sy+M_th*g_th*R_coil-Q_tooth*R_coil)/(R_coil+R_sy)
        if Q_coil<0:
            raise mo.InvalidDesign(message='Negative Q_coil')
        A_slot=np.pi*(r_sy**2-r_si**2)/Q -w_st*(r_sy-r_si)        
        J_hat=np.sqrt(Q_coil*sigma*np.sqrt(2)/(l_turn*k_fill*A_slot))
        return [J_hat,Q_coil,Q_sy,Q_tooth]


class ThermalPostAnalyzer(me.PostAnalyzer):
    """Converts input state into output state for TemplateAnalyzer"""
    def get_next_state(self,results:Any,stateIn:'me.State')->'me.State':
        stateOut=deepcopy(stateIn)
        stateOut.conditions.J=results[0]
        stateOut.conditions.Q_coil=results[1]
        stateOut.conditions.Q_yoke=results[2]
        stateOut.conditions.Q_tooth=results[3]
        return stateOut
    

problem_def=ThermalProblemDefinition()
analyzer=ThermalAnalyzer()
post_analyzer = ThermalPostAnalyzer()

ThermalStep=me.AnalysisStep(problem_def, analyzer, post_analyzer)