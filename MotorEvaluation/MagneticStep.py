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

__all__=['MagStep',]

class MagProblemDefinition(me.ProblemDefinition):
    """Class converts input state into a problem"""
    
    def get_problem(self,state:'me.State')->'me.Problem':
        """Returns Problem from Input State"""
        #TODO define problem definition
        magnet_mat=state.design.machine.magnet_mat
        core_mat=state.design.machine.core_mat
        
        B_r=magnet_mat.B_r
        mu_m=magnet_mat.mu_m
        mu_core=core_mat.mu_core
        
        
        d_m=state.design.machine.d_m
        A_m=np.pi*state.design.machine.r_m
        d_ag=state.design.machine.d_ag
        w_st=state.design.machine.w_st
        l_st=state.design.machine.l_st
        l_tooth=state.design.machine.l_tooth
        l_yoke=np.pi*(state.design.machine.r_sy+state.design.machine.r_so)/2
        d_yoke=state.design.machine.d_yoke

        
        problem=MagProblem(B_r,d_m,mu_m,A_m,d_ag,w_st,
                          l_st,mu_core,l_tooth,l_yoke,d_yoke)
        return problem

class MagProblem():
    """problem class utilized by the Analyzer
    
    Attributes:
        TODO
    """
    def __init__(self,B_r,d_m,mu_m,A_m,d_ag,w_st,
                 l_st,mu_core,l_tooth,l_yoke,d_yoke):
        """Creates problem class
        
        Args:
            TODO
            
        """
        self.B_r=B_r
        self.d_m=d_m
        self.mu_m=mu_m
        self.A_m=A_m
        self.d_ag=d_ag

        self.w_st=w_st
        self.l_st=l_st
        self.mu_core=mu_core
        self.l_tooth=l_tooth
        self.l_yoke=l_yoke
        self.d_yoke=d_yoke
    
class MagAnalyzer(me.Analyzer):
    """"Class Analyzes the CubiodProblem  for volume and Surface Areas"""
    
    def analyze(self,problem:'me.Problem'):
        """Performs Analysis on a problem

        Args:
            problem (me.Problem): Problem Object

        Returns:
            results (Any): 
                Results of Analysis

        """
        B_r=problem.B_r
        d_m=problem.d_m
        mu_m=problem.mu_m
        A_m=problem.A_m
        d_ag=problem.d_ag
        mu_0=4*np.pi*10**-7
        w_st=problem.w_st
        l_st=problem.l_st
        mu_core=problem.mu_core
        l_tooth=problem.l_tooth
        l_yoke=problem.l_yoke
        d_yoke=problem.d_yoke
        F=B_r*d_m/(mu_m*mu_0)
        R_mag=d_m/(mu_0*mu_m*A_m)
        R_ag=d_ag/(mu_0*A_m)
        R_tooth=l_tooth/(mu_0*mu_core*w_st*l_st)
        R_yoke=(l_yoke/2)/(mu_0*mu_core*d_yoke*l_st)
        slots_per_pole=3
        R_total=R_mag+R_ag+R_tooth/slots_per_pole+R_yoke
        phi=F/R_total
        
        B_tooth=phi/(slots_per_pole*w_st*l_st)
        if B_tooth > 1.6:
            B_tooth=1.6
            phi= B_tooth*slots_per_pole*w_st*l_st
        B_yoke=phi/(d_yoke*l_st)
        if B_yoke > 1.6:
            B_yoke=1.6
            phi= B_yoke*w_st*l_st
            B_tooth=phi/(w_st*l_st)
        B_delta=phi/(A_m)
        #B_delta=d_m*B_r/(mu_m*d_ag+d_m)
        results=[B_delta,B_tooth,B_yoke]
        return results
    


class MagPostAnalyzer(me.PostAnalyzer):
    """Converts input state into output state for TemplateAnalyzer"""
    def get_next_state(self,results:Any,stateIn:'me.State')->'me.State':
        stateOut=deepcopy(stateIn)
        stateOut.conditions.B_delta=results[0]
        stateOut.conditions.B_tooth=results[1]
        stateOut.conditions.B_yoke=results[2]
        return stateOut
    

problem_def=MagProblemDefinition()
analyzer=MagAnalyzer()
post_analyzer = MagPostAnalyzer()

MagStep=me.AnalysisStep(problem_def, analyzer, post_analyzer)