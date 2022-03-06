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
        core_mat=state.design.machine.core_mat
        coil_mat=state.design.machine.coil_mat
        ins_mat=state.design.machine.ins_mat



        l_st=state.design.machine.l_st
        alpha_q=state.design.machine.alpha_q

        T_coil = coil_mat.T_coil
        h = state.design.settings.h
        k_fe = core_mat.k
        k_ins = ins_mat.k
        w_ins = state.design.machine.w_ins
        w_st=state.design.machine.w_st
        r_so = state.design.machine.r_so
        r_sy = state.design.machine.r_sy
        r_si = state.design.machine.r_si
        r_ro = state.design.machine.r_ro
        A_slot=state.design.machine.A_slot
        alpha_q=alpha_q
        omega=state.conditions.v_tip_max/r_ro
        f=omega/(2*np.pi)
        k_stien=core_mat.k_stien
        a=core_mat.a
        b=core_mat.b
        B_yoke=state.conditions.B_yoke
        B_tooth=state.conditions.B_tooth
        Q=state.design.machine.Q
        k_ov=coil_mat.k_ov
        y=state.design.machine.y
        sigma_coil=coil_mat.sigma_coil
        k_fill=coil_mat.k_fill
        A_o = alpha_q*r_so*l_st
        
        problem=ThermalProblem(l_st,A_o,T_coil,h,k_fe,k_ins,w_ins,w_st,r_so,r_sy,r_si,alpha_q,
                 f,k_stien,a,b,B_yoke,B_tooth,Q,k_ov,y,sigma_coil,k_fill,A_slot)
        return problem

class ThermalProblem():
    """problem class utilized by the Analyzer
    
    Attributes:
        TODO
    """
    def __init__(self,l_st,A_o,T_coil,h,k_fe,k_ins,w_ins,w_st,r_so,r_sy,r_si,alpha_q,
                 f,k_stien,a,b,B_yoke,B_tooth,Q,k_ov,y,sigma_coil,k_fill,A_slot):
        """Creates problem class
        
        Args:
            TODO
            
        """
        self.l_st = l_st
        self.A_o = A_o
        self.T_coil = T_coil
        self.h = h
        self.k_fe = k_fe
        self.k_ins =k_ins
        self.w_ins = w_ins
        self.w_st = w_st
        self.r_so = r_so
        self.r_sy = r_sy
        self.r_si = r_si
        self.alpha_q=alpha_q
        self.f=f
        self.k_stien=k_stien
        self.a=a
        self.b=b
        self.B_yoke=B_yoke
        self.B_tooth=B_tooth
        self.Q=Q
        self.k_ov=k_ov
        self.y=y
        self.sigma_coil=sigma_coil
        self.k_fill=k_fill
        self.A_slot=A_slot
    
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
        l_st = problem.l_st
        A_o = problem.A_o
        T_coil = problem.T_coil
        h = problem.h
        k_fe = problem.k_fe
        k_ins = problem.k_ins
        w_ins = problem.w_ins
        w_st = problem.w_st
        r_so = problem.r_so
        r_sy = problem.r_sy
        alpha_q=problem.alpha_q
        f=problem.f
        k_stien=problem.k_stien
        a=problem.a
        b=problem.b
        B_yoke=problem.B_yoke
        B_tooth=problem.B_tooth
        Q=problem.Q
        k_ov=problem.k_ov
        y=problem.y
        sigma_coil=problem.sigma_coil
        k_fill=problem.k_fill
        A_slot=problem.A_slot
        
        zeta = np.sqrt(k_ins/(k_fe*w_st*w_ins))
        r_si = problem.r_si
        F_1 = -(2/zeta)*(1/np.sinh(zeta*(r_sy-r_si))-np.cosh(zeta*(r_sy-r_si))/np.sinh(zeta*(r_sy-r_si)))+alpha_q*r_sy        
        F_2 = (1/zeta)*(1/np.sinh(zeta*(r_sy-r_si))-np.cosh(zeta*(r_sy-r_si))/np.sinh(zeta*(r_sy-r_si)))+(r_sy-r_si)        
        V_yoke=alpha_q/2*(r_so**2-r_sy**2)*l_st
        V_tooth=w_st*l_st*(r_so-r_sy)
        g_yoke = k_stien*f**a*B_yoke**b
        g_tooth = k_stien*f**a*B_tooth**b
        P_tooth = g_tooth*V_tooth
        P_yoke = g_yoke*V_yoke
        
        L_mean=2*(l_st+k_ov*np.pi*y*(r_si+r_sy)/Q)
        P_coil = l_st*(4*A_o*F_1*T_coil*h*k_fe*k_ins - 2*A_o*F_1*g_yoke*h*k_ins*r_so**2*np.log(r_so/r_sy) + A_o*F_1*g_yoke*h*k_ins*r_so**2 - A_o*F_1*g_yoke*h*k_ins*r_sy**2 - 4*A_o*F_2*g_tooth*h*k_fe*w_ins*w_st - 4*F_1*P_tooth*h*k_ins*r_so*np.log(r_so/r_sy) - 4*F_1*P_tooth*k_fe*k_ins - 4*F_1*P_yoke*h*k_ins*r_so*np.log(r_so/r_sy) - 4*F_1*P_yoke*k_fe*k_ins)/(4*(A_o*h*k_fe*w_ins + F_1*h*k_ins*l_st*r_so*np.log(r_so/r_sy) + F_1*k_fe*k_ins*l_st))
        if P_coil< 0:
            raise mo.InvalidDesign(message='Iron losses too high')
        J=np.sqrt(P_coil*sigma_coil/(L_mean*k_fill*A_slot))
        return [J,P_coil,P_yoke,P_tooth]


class ThermalPostAnalyzer(me.PostAnalyzer):
    """Converts input state into output state for TemplateAnalyzer"""
    def get_next_state(self,results:Any,stateIn:'me.State')->'me.State':
        stateOut=deepcopy(stateIn)
        stateOut.conditions.J=results[0]
        stateOut.conditions.P_coil=results[1]
        stateOut.conditions.P_yoke=results[2]
        stateOut.conditions.P_tooth=results[3]
        return stateOut
    

problem_def=ThermalProblemDefinition()
analyzer=ThermalAnalyzer()
post_analyzer = ThermalPostAnalyzer()

ThermalStep=me.AnalysisStep(problem_def, analyzer, post_analyzer)