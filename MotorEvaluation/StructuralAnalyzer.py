# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 11:06:23 2022

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

 __all__=[StructuralStep,]
 
class StructuralProblemDefinition(me.ProblemDefinition):
    """Class converts input state into a problem"""
    
    def get_problem(self,state:'me.State')->'me.Problem':
        """Returns Problem from Input State"""
        #TODO define problem definition
        magnet_mat=state.design.machine.magnet_mat
        rho=magnet_mat.rho
        E=magnet_mat.E
        v=magnet_mat.v
        alpha=magnet_mat.alpha
        sigma_t_max=magnet_mat.sigma_t_max
        r_m=state.design.machine.r_m
        DT=10
        problem=StructuralProblem(rho,E,v,alpha,sigma_t_max,r_m,DT)
        return problem

class StructuralProblem():
    """problem class utilized by the Analyzer
    
    Attributes:
        TODO
    """
    def __init__(self,rho,E,v,alpha,sigma_t_max,r_m,DT):
        """Creates problem class
        
        Args:
            TODO
            
        """
        self.rho=rho
        self.E=E
        self.v = v
        self.alpha= alpha
        self.sigma_t_max = sigma_t_max
        self.r_m=r_m
        self.DT= DT
        
    
class StructuralAnalyzer(me.Analyzer):
    """"Class Analyzes the CubiodProblem  for volume and Surface Areas"""
    
    def analyze(self,problem:'me.Problem'):
        """Performs Analysis on a problem

        Args:
            problem (me.Problem): Problem Object

        Returns:
            results (Any): 
                Results of Analysis

        """
        rho=problem.rho
        E=problem.E
        v=problem.v
        alpha=problem.alpha
        C1=E/(1-v**2)
        C2=E*v/(1-v**2)
        C3=C1
        P_sl=problem.P_sl
        D_1=C1+C2
        D_2=C2-C1
        D_3=3*C1+C2
        D_4=C2+C1
        D_5=C1-C2
        D_6=3*C2+C1
        PSI=-rho/(8*C1)
        zeta_r=-(C1+C2)*alpha
        zeta_t=-(C2-C3)*alpha
        sigma_t_max=problem.sigma_t_max
        sigma_r_max=0
        r_m=problem.r_m
        DT=problem.DT
        v_tip_max_1=np.sqrt((sigma_r_max-P_sl)/(D_3*PSI*(r_m**2-1)))
        v_tip_max_2=np.sqrt((sigma_t_max-(D_4/D_1)*P_sl-(zeta_t-(D_4/D_1)*zeta_r)*DT)/(PSI*(D_6*r_m**2-(D_3*D_4)/D_1)))
        v_tip_max=np.minimum(v_tip_max_1,v_tip_max_2)
        results = [v_tip_max]
        return results
    


class StructuralPostAnalyzer(me.PostAnalyzer):
    """Converts input state into output state for TemplateAnalyzer"""
    def get_next_state(self,results:Any,stateIn:'me.State')->'me.State':
        stateOut=deepcopy(stateIn)
        stateOut.conditions.v_tip_max=results
        return stateOut
    

problem_def=StructuralProblemDefinition()
analyzer=StructuralAnalyzer()
post_analyzer = StructuralPostAnalyzer()

StructuralStep=me.AnalysisStep(problem_def, analyzer, post_analyzer)
