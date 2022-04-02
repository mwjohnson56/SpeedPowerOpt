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

__all__=['CoreLossStep',]

class CoreLossProblemDefinition(me.ProblemDefinition):
    """Class converts input state into a problem"""
    
    def get_problem(self,state:'me.State')->'me.Problem':
        """Returns Problem from Input State"""
        #TODO define problem definition
      
        p=state.design.machine.p
        Omega=state.design.settings.Omega
        K_h=state.design.machine.core_mat['core_ironloss_Kh']
        a=state.design.machine.core_mat['core_ironloss_a']
        b=state.design.machine.core_mat['core_ironloss_b']
        K_e=state.design.machine.core_mat['core_ironloss_Ke']
        rho=state.design.machine.core_mat['core_material_density']
        k_stack=state.design.machine.core_mat['core_stacking_factor']
        B_sy=state.design.machine.B_sy
        B_tooth=state.design.machine.B_th
        problem=CoreLossProblem(p,Omega,K_h,a,b,K_e,rho,k_stack,B_sy,B_tooth)
        return problem

class CoreLossProblem():
    """problem class utilized by the Analyzer
    
    Attributes:
        TODO
    """
    def __init__(self,p,Omega,K_h,a,b,K_e,rho,k_stack,B_sy,B_tooth):
        """Creates problem class
        
        Args:
            TODO
            
        """
        self.p=p
        self.Omega=Omega
        self.K_h=K_h
        self.a=a
        self.b=b
        self.K_e=K_e
        self.rho=rho
        self.k_stack=k_stack
        self.B_sy=B_sy
        self.B_tooth=B_tooth
    
class CoreLossAnalyzer(me.Analyzer):
    """"Class Analyzes the CubiodProblem  for volume and Surface Areas"""
    
    def analyze(self,problem:'me.Problem'):
        """Performs Analysis on a problem

        Args:
            problem (me.Problem): Problem Object

        Returns:
            results (Any): 
                Results of Analysis

        """
        p=problem.p
        Omega=problem.Omega
        K_h=problem.K_h
        b=problem.a
        a=problem.b
        K_e=problem.K_e
        rho=problem.rho
        k_stack=problem.k_stack
        B_sy=problem.B_sy
        B_tooth=problem.B_tooth
        f=p*Omega/(2*np.pi)
        g_sy=(K_h*f**a*B_sy**b + K_e*(f*B_sy)**2)*k_stack
        g_tooth=(K_h*f**a*B_tooth**b + K_e*(f*B_tooth)**2)*k_stack
        results=[g_sy,g_tooth]
        return results
    


class CoreLossPostAnalyzer(me.PostAnalyzer):
    """Converts input state into output state for TemplateAnalyzer"""
    def get_next_state(self,results:Any,stateIn:'me.State')->'me.State':
        stateOut=deepcopy(stateIn)
        stateOut.conditions.g_sy=results[0]
        stateOut.conditions.g_tooth=results[1]
        return stateOut
    

problem_def=CoreLossProblemDefinition()
analyzer=CoreLossAnalyzer()
post_analyzer = CoreLossPostAnalyzer()

CoreLossStep=me.AnalysisStep(problem_def, analyzer, post_analyzer)