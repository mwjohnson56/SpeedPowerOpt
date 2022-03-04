# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 10:51:03 2022

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

__all__ = ['SPM_Designer',]


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
        d_yoke=x[5]
        k_tooth=x[6]
        machine=self.MotorClass(r_sh,r_ro,d_m,d_ag,l_tooth,d_yoke,k_tooth,
                 self.shaft_mat,self.magnet_mat,self.core_mat,
                 self.coil_mat,self.sleeve_mat)
        return machine
    
class SPM_SettingsHandler():
    """Settings hanlder for design creation"""
    def __init__(self,cooling):
        self.cooling=cooling
    def get_settings(self,x):
        settings = self.cooling
        return settings
    
class Q6p1y1_SMP_Motor(me.Machine):
    """Class defines a Machine object 
    
    Attributes:
        TODO
    """
    
    def __init__(self,r_sh,r_ro,d_m,d_ag,l_tooth,d_yoke,k_tooth,
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
        self._d_yoke= d_yoke
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
    def d_yoke(self):
        return self._d_yoke
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
        return self.r_sy+self.d_yoke
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

class MagnetMaterial:
    rho=7450
    E=160E9
    v=.24
    alpha=5E-6
    B_R=1.7
    mu_r=1
class ShaftMaterial:
    rho=7450
    E=160E9
    v=.24
    alpha=5E-6   
class CoreMaterial:
    mu_core=4000
    
    
arch=SPM_MotorArchitect(ShaftMaterial, MagnetMaterial, CoreMaterial, None, None, Q6p1y1_SMP_Motor)
settings_handler=SPM_SettingsHandler(None)
SPM_Designer=me.MachineDesigner(arch,settings_handler)