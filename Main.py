# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 10:59:45 2021

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


class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class SPM_MotorArchitect(me.Architect):
    """Class converts input tuple x into a machine object"""   
    def __init__(self,shaft_mat,magnet_mat,core_mat,coil_mat,sleeve_mat,MotorClass):
        self.shaft_mat=shaft_mat
        self.magnet_mat=magnet_mat
        self.core_mat=core_mat
        self.coil_mat=coil_mat
        self.sleeve_mat=sleeve_mat
        self.MotorClass=MotorClass
    def create_new_design(self,x:tuple)->"me.Machine":
        """
        converts x tuple into a machine object.

        Args:
            x (tuple): Input free variables.
            
        Returns:
            machine (me.Machine): Machine object
        """
        
        r_sh=x[0]
        r_ro=x[1]
        d_m=x[2]
        d_ag=x[3]
        l_tooth=x[4]
        l_yoke=x[5]
        k_tooth=x[6]
        machine=self.MotorClass(r_sh,r_ro,d_m,d_ag,l_tooth,l_yoke,k_tooth,
                 self.shaft_mat,self.magnet_mat,self.core_mat,
                 self.coil_mat,self.sleeve_mat)
        return machine
    
class TemplateSettingsHanlder(me.SettingsHandler):
    """Settings hanlder for design creation"""
    def __init__(self,cooling):
        self.cooling=cooling
    def getSettings(self,x):
        settings = self.cooling
        return settings
    
class Q6p1y1_SMP_Motor(me.Machine):
    """Class defines a Machine object 
    
    Attributes:
        TODO
    """
    
    def __init__(self,r_sh,r_ro,d_m,d_ag,l_tooth,l_yoke,k_tooth,
                 shaft_mat,magnet_mat,core_mat,coil_mat,sleeve_mat):
        """Creates a machine object.

        Args:
            TODO

        """
        
        self._r_sh=r_sh
        self._r_ro=r_ro
        self._d_m=d_m
        self._d_ag=d_ag
        self._l_tooth=l_tooth
        self._l_yoke= l_yoke
        self._k_tooth=k_tooth
        self._shaft_mat=shaft_mat
        self._magent_mat=magnet_mat
        self._core_mat=core_mat
        self._coil_mat=coil_mat
        self._sleeve_mat=sleeve_mat
        
        self._d_sl=self.d_ag/2
        self._l_st=1
        self.verify_design()
        
    def verify_design(self):
        if self.r_sh> self.r_ro-self.d_m:
            raise mo.InvalidDesign(message='Invalid rotor geometry')
    
    #Inputed properties   
    @property
    def r_sh(self):
        return self._r_sh
    @property
    def r_ro(self):
        return self._r_ro
    @property
    def d_m(self):
        return self._d_m
    @property
    def d_ag(self):
        return self._d_ag
    @property
    def l_tooth(self):
        return self._l_tooth
    @property
    def l_yoke(self):
        return self._l_yoke
    @property
    def k_tooth(self):
        return self._k_tooth
    @property
    def shaft_mat(self):
        return self._shaft_mat
    @property
    def magnet_mat(self):
        return self._magnet_mat
    @property
    def core_mat(self):
        return self._core_mat
    @property
    def coil_mat(self):
        return self._coil_mat
    @property
    def sleeve_mat(self):
        return self._sleeve_mat

    
    #Updated properties
    @property
    def d_sl(self):
        return self._d_sl     
    @d_sl.setter
    def d_sl(self,d_sl):
        self._d_sl=d_sl
    @property
    def l_st(self):
        return self._l_st     
    @l_st.setter
    def l_st(self,l_st):
        self._l_st=l_st
    
    #Design Constants
    @property
    def Q(self):
        return 6
    @property
    def p(self):
        return 1
    @property
    def alpha_q(self):
        return 3.1415/self.Q
    
    #Calculated Properties
    @property
    def r_m(self):
        return self.r_ro-self.d_m
    @property
    def r_si(self):
        return self.r_ro+self.d_ag
    @property
    def r_sy(self):
        return self.r_si+self.l_tooth
    @property
    def r_so(self):
        return self.r_sy+self.l_yoke
    @property
    def w_st(self):
        return self.k_tooth*self.r_si*self.alpha_q
        
    def check_required_properties(self):
        """Checks for required input properties"""
        #TODO 
        pass
    
    
    def get_missing_properties(self):
        """Returns missing input properites"""
        #TODO
        pass

    
class TemplateProblemDefinition(me.ProblemDefinition):
    """Class converts input state into a problem"""
    
    def getProblem(self,state:'me.State')->'me.Problem':
        """Returns Problem from Input State"""
        #TODO define problem definition
        args=None
        problem=TemplateProblem(args)
        return problem

class TemplateProblem():
    """problem class utilized by the Analyzer
    
    Attributes:
        TODO
    """
    def __init__(self,args):
        """Creates problem class
        
        Args:
            TODO
            
        """
        #TODO define problem 
        
    
class TemplateAnalyzer(me.Analyzer):
    """"Class Analyzes the CubiodProblem  for volume and Surface Areas"""
    
    def analyze(self,problem:'me.Problem'):
        """Performs Analysis on a problem

        Args:
            problem (me.Problem): Problem Object

        Returns:
            results (Any): 
                Results of Analysis

        """
        #TODO Define Analyzer
        results = []
        return results
    


class TemplatePostAnalyzer(me.PostAnalyzer):
    """Converts input state into output state for TemplateAnalyzer"""
    def getNextState(self,results:Any,stateIn:'me.State')->'me.State':
        stateOut=deepcopy(stateIn)
        #TODO define Post-Analyzer
        return stateOut
    

class ConstraintError(Error):
    """Error for violating optimization constraint"""
    def __init__(self,value):
        #TODO define error
        self.value=value

    
class TemplateConstraintEvaluationStep():
    """Constraint evaluation step template"""
    def step(self,stateIn):
        """Checks input state to see if constraint is violated
        
        Raises ConstraintError if violated, otherwise appends values to 
        State conditions and moves forward"""
        
        value = None #TODO define constraint
        if value >=0:
            raise ConstraintError(value)
        else:
            stateOut=deepcopy(stateIn)
            stateOut.stateConditions.constraintValue=value
            return [value,stateOut]
    
    
class DataHandler:
    def save(self,design,fullResults,objs):
        """Unimplented data handler"""
        #TODO Define datahandler
        pass

#%%

if __name__ == '__main__':
    
    #Create Designer
    des=me.MachineDesigner(TemplateArchitect(),TemplateSettingsHanlder())

    #Create evaluation steps
    evalSteps=[TemplateConstraintEvaluationStep(),
               me.AnalysisStep(TemplateProblemDefinition(),
                               TemplateAnalyzer(),
                               TemplatePostAnalyzer())]#TODO define steps
    #Create Evaluator
    evaluator=me.MachineEvaluator(evalSteps)
    objectives=TemplateObjective()
    dh=DataHandler()
    
    #set evaluation bounds
    bounds=([0,0,0],[1,1,1])
    #set number of objectives
    n_obj=3
    
    #Create Machine Design Problem
    machDesProb=do.DesignProblem(des,evaluator,objectives,dh,
                                        bounds,n_obj)
    
    #Run Optimization
    opt=do.DesignOptimizationMOEAD(machDesProb)
    pop=opt.run_optimization(496,10)
    fits, vectors = pop.get_f(), pop.get_x()
    ndf, dl, dc, ndr = pg.fast_non_dominated_sorting(fits) 