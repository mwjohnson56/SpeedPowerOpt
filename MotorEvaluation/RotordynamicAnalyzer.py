# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 11:35:30 2022

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

__all__=['RDStep',]

class RDProblemDefinition(me.ProblemDefinition):
    """Class converts input state into a problem"""
    
    def get_problem(self,state:'me.State')->'me.Problem':
        """Returns Problem from Input State"""
        #TODO define problem definition
        shaft_mat=state.design.machine.shaft_mat
        rho=shaft_mat.rho
        E=shaft_mat.E
        C_sh=np.sqrt(E/rho)
        r_ro=state.design.machine.r_ro
        alpha_r=state.design.machine.r_sh/r_ro
        alpha_l=2
        v_tip=state.conditions.v_tip_max
        problem=RDProblem(r_ro,alpha_r,C_sh,alpha_l,v_tip)
        return problem

class RDProblem():
    """problem class utilized by the Analyzer
    
    Attributes:
        TODO
    """
    def __init__(self,r_ro,alpha_r,C_sh,alpha_l,v_tip):
        """Creates problem class
        
        Args:
            TODO
            
        """
        self.r_ro=r_ro
        self.alpha_r=alpha_r
        self.C_sh = C_sh
        self.alpha_l=alpha_l
        self.v_tip=v_tip
    
class RDAnalyzer(me.Analyzer):
    """"Class Analyzes the CubiodProblem  for volume and Surface Areas"""
    
    def analyze(self,problem:'me.Problem'):
        """Performs Analysis on a problem

        Args:
            problem (me.Problem): Problem Object

        Returns:
            results (Any): 
                Results of Analysis

        """
        r_ro=problem.r_ro
        alpha_r=problem.alpha_r
        C_sh = problem.C_sh
        alpha_l=problem.alpha_l
        v_tip=problem.v_tip
        k_w=1.1
        l_st=r_ro*np.sqrt((alpha_r*4.71**2*C_sh)/(2*alpha_l**2*k_w*v_tip))
        results = [l_st]
        return results
    


class RDPostAnalyzer(me.PostAnalyzer):
    """Converts input state into output state for TemplateAnalyzer"""
    def get_next_state(self,results:Any,stateIn:'me.State')->'me.State':
        stateOut=deepcopy(stateIn)
        stateOut.design.machine.l_st=results
        return stateOut
    

problem_def=RDProblemDefinition()
analyzer=RDAnalyzer()
post_analyzer = RDPostAnalyzer()

RDStep=me.AnalysisStep(problem_def, analyzer, post_analyzer)