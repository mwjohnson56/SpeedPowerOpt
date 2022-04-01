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
from eMach import mach_opt as mo
from eMach import mach_eval as me

__all__=['RDStep',]

class RDProblemDefinition(me.ProblemDefinition):
    """Class converts input state into a problem"""
    
    def get_problem(self,state:'me.State')->'me.Problem':
        """Returns Problem from Input State"""
        #TODO define problem definition
        shaft_mat=state.design.machine.shaft_mat
        rho=shaft_mat['shaft_material_density']
        E=shaft_mat['shaft_youngs_modulus']
        C_sh=np.sqrt(E/rho)
        r_sh=state.design.machine.r_sh
        l_st=state.design.machine.l_st
        Omega=state.design.settings.Omega
        alpha_l=2
        problem=RDProblem(r_sh,C_sh,l_st,Omega,alpha_l)
        return problem

class RDProblem():
    """problem class utilized by the Analyzer
    
    Attributes:
        TODO
    """
    def __init__(self,r_sh,C_sh,l_st,Omega,alpha_l):
        """Creates problem class
        
        Args:
            TODO
            
        """
        self.r_sh=r_sh
        self.C_sh = C_sh
        self.l_st=l_st
        self.Omega=Omega
        self.alpha_l=alpha_l
    
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
        r_sh=problem.r_sh
        C_sh = problem.C_sh
        alpha_l=problem.alpha_l
        l_st=problem.l_st

        omega_n=(4.7**2)*(r_sh/(alpha_l*l_st))*C_sh
        if 1.2*problem.Omega>omega_n:
            raise mo.InvalidDesign(message='Critical Speed')
        else:
            return True
    


class RDPostAnalyzer(me.PostAnalyzer):
    """Converts input state into output state for TemplateAnalyzer"""
    def get_next_state(self,results:Any,stateIn:'me.State')->'me.State':
        stateOut=deepcopy(stateIn)
        stateOut.design.machine.l_s_valid=results
        return stateOut
    

problem_def=RDProblemDefinition()
analyzer=RDAnalyzer()
post_analyzer = RDPostAnalyzer()

RDStep=me.AnalysisStep(problem_def, analyzer, post_analyzer)