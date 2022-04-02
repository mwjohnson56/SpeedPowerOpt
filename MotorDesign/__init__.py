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
from .Materials import *
__all__ = ['arch','SPM_SettingsHandler']


class SPM_MotorArchitect(me.Architect):
    """Class converts input tuple x into a machine object"""   
    def __init__(self,shaft_mat,magnet_mat,core_mat,
                 coil_mat,sleeve_mat,ins_mat,MotorClass):
        self.shaft_mat=shaft_mat
        self.magnet_mat=magnet_mat
        self.core_mat=core_mat
        self.coil_mat=coil_mat
        self.sleeve_mat=sleeve_mat
        self.ins_mat=ins_mat
        self.MotorClass=MotorClass
    def create_new_design(self,x:tuple)->"me.Machine":
        """
        converts x tuple into a machine object.

        Args:
            x (tuple): Input free variables.
            
        Returns:
            machine (me.Machine): Machine object
        """
        
        
        r_ro=x[0]
        d_m_norm=x[1]
        d_m=d_m_norm*r_ro
        r_sh=r_ro-d_m
        d_ag=.002
        l_st=x[2]*r_ro
        B_sy_norm=x[3]
        B_th_norm=x[4]
        l_tooth_norm=x[5]
        
        Q= self.MotorClass.Q_motor
        p= self.MotorClass.p_motor
        
        r_si=r_ro+d_ag
        alpha_q=2*np.pi/Q
        
        B_delta=d_m*self.magnet_mat['B_r']/(self.magnet_mat['mu_r']*d_ag+d_m)
        B_sy=B_sy_norm*B_delta
        B_th=B_th_norm*B_delta
        d_yoke=np.pi*B_delta*r_si/(2*p*B_sy)
        w_tooth=B_delta*r_si*alpha_q/(B_th)
        l_tooth=l_tooth_norm*r_ro
        

        machine=self.MotorClass(r_sh,r_ro,d_m,d_ag,l_st,d_yoke,w_tooth,l_tooth,
                 self.shaft_mat,self.magnet_mat,self.core_mat,
                 self.coil_mat,self.sleeve_mat,self.ins_mat)
        return machine
    
class SPM_SettingsHandler():
    """Settings hanlder for design creation"""
    def __init__(self,Omega,h):
        self.Omega=Omega
        self.h=h
    def get_settings(self,x):
        settings = Settings(self.h,self.Omega)
        return settings

class Settings:
    def __init__(self,h,Omega):
        self.Omega=Omega
        self.h=h
    @property
    def speed(self):
        return self.Omega*60/(np.pi*2)
    @property
    def rotor_temp_rise(self):
        return 0
        
class Q6p1y1_SMP_Motor(me.Machine):
    """Class defines a Machine object 
    
    Attributes:
        TODO
    """
    Q_motor=6
    p_motor=1
    y_motor=3
    
    def __init__(self,r_sh,r_ro,d_m,d_ag,l_st,d_yoke,w_tooth,l_tooth,
                 shaft_mat,magnet_mat,core_mat,coil_mat,sleeve_mat,ins_mat):
        """Creates a machine object.

        Args:
            TODO

        """
        
        self._r_sh=r_sh
        self._r_ro=r_ro
        self._d_m=d_m
        self._d_ag=d_ag
        self._l_st=l_st
        self._d_yoke= d_yoke
        self._w_tooth=w_tooth
        self._shaft_mat=shaft_mat
        self._magnet_mat=magnet_mat
        self._core_mat=core_mat
        self._coil_mat=coil_mat
        self._sleeve_mat=sleeve_mat
        self.ins_mat=ins_mat
        self._d_sl=self.d_ag/2
        self._l_tooth=l_tooth
        self.verify_design()
        
    def verify_design(self):
        if self.B_sy> self.core_mat['core_saturation_feild']:
            raise mo.InvalidDesign(message='Saturated yoke feild')
        if self.B_th> self.core_mat['core_saturation_feild']:
            raise mo.InvalidDesign(message='Saturated tooth feild')
        if self.d_m < .0001:
            raise mo.InvalidDesign(message='Minnimum magnet thickness')
    
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
    def l_st(self):
        return self._l_st
    @property
    def d_yoke(self):
        return self._d_yoke
    @property
    def w_tooth(self):
        return self._w_tooth
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
    def l_tooth(self):
        return self._l_tooth     
    @l_tooth.setter
    def l_tooth(self,l_tooth):
        self._l_tooth=l_tooth
    
    #Design Constants
    @property
    def Q(self):
        return self.Q_motor
    @property
    def y(self):
        return self.y_motor
    @property
    def p(self):
        return self.p_motor
    @property
    def alpha_q(self):
        return 2*np.pi/self.Q
    @property
    def w_ins(self):
        return .001
    
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
    def B_delta(self):
        return self.d_m*self.magnet_mat['B_r']/(self.magnet_mat['mu_r']*self.d_ag+self.d_m)
    @property
    def B_sy(self):
        return np.pi*self.B_delta*self.r_si/(2*self.p*(self.r_so-self.r_sy))
    @property
    def B_th(self):
        return self.B_delta*self.r_si*self.alpha_q/(self.w_tooth)
    @property 
    def alpha_slot(self):
        return self.alpha_q-2*np.arctan(self.w_tooth/(2*self.r_sy))
    @property
    def k_w(self):
        alpha=np.pi*((self.Q-2*self.y)/(self.Q*self.p))
        n=self.Q/(2*self.p)
        m=self.Q/(6*self.p)
        Beta=np.pi/n
        k_w=np.cos(alpha/2)*(np.sin(m*Beta/2))/(m*np.sin(Beta/2))
        self._k_w=k_w
        return self._k_w
    @property
    def A_slot(self):
        return np.pi*(self.r_sy**2-self.r_si**2)/self.Q - \
            self.w_tooth*(self.r_sy-self.r_si)
        # A_slot=((1-self.k_tooth)*self.alpha_q/2)*(self.r_sy**2-self.r_si**2)+\
        #     (self.r_sy**2-self.r_si**2)*np.tan(self.k_tooth*self.alpha_q/2)
        # return A_slot

    def A_hat(self,J):
        A_hat=3*self.A_slot*J*(self.Q/3)*self.k_w*self.coil_mat['k_fill']/(np.pi*self.r_si)
        self._A_hat=A_hat
        return A_hat
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
    B_r=1.285
    mu_m=1.062
    sigma_t_max=80E6
class ShaftMaterial:
    rho=7450
    E=160E9
    v=.24
    alpha=5E-6   
class CoreMaterial:
    mu_core=4000
    k=40
    k_stien=.001
    a=1.5
    b=2
class CoilMaterial:
    T_coil=150
    k_ov=2
    sigma_coil=5.80E7
    k_fill=.38
class InsMaterial:
    k=1
class AirCooled:
    h=100

arch=SPM_MotorArchitect(Steel, N40H, M19Gauge29,
                        CoilCopper, CarbonFiber,Ins,
                        Q6p1y1_SMP_Motor)
Omega=1000
h=100
settings_handler=SPM_SettingsHandler(h,Omega)
SPM_Designer=me.MachineDesigner(arch,settings_handler)