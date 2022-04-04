# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 16:29:53 2022

@author: Martin Johnson
"""


CarbonFiber = {
    'sleeve_material'            : 'CarbonFiber',
    'sleeve_material_density'    : 1800, # kg/m3
    'sleeve_youngs_th_direction' : 125E9,  #Pa
    'sleeve_youngs_p_direction'  : 8.8E9,  #Pa
    'sleeve_poission_ratio_p'    :.015,
    'sleeve_poission_ratio_tp'   :.28,
    'sleeve_safety_factor'       : 1.5, #
    'sleeve_max_tan_stress'      : 1950E6, # Pa
    'sleeve_max_rad_stress'      : -100E6, # Pa
    'sleeve_therm_conductivity'  : 0.71, # W/m-k
    }

Steel = {
    'shaft_material'             : 'Steel',
    'shaft_material_density'     : 7870, # kg/m3
    'shaft_youngs_modulus'       : 206E9, #Pa
    'shaft_poission_ratio'       : .3, #[]
    'shaft_therm_conductivity'   : 51.9, # W/m-k  
    }

Copper = {
    'copper_material_cost'       : 73228, # $/m3
    'copper_elec_conductivity'   : 5.7773*1e7, # S/m
    }

Air = {
    'air_therm_conductivity'     :.02624, #W/m-K
    'air_viscosity'              :1.562E-5, #m^2/s
    'air_cp'                     :1, #kJ/kg
    'air_temp'                   :25, #[C]
    }

Hub = {
       'rotor_hub_therm_conductivity':205.0, #W/m-K       
      }

N40H = {
    'magnet_material'            : "Arnold/Reversible/N40H",
    'magnet_material_density'    : 7450, # kg/m3
    'magnet_youngs_modulus'      : 160E9, # Pa
    'magnet_poission_ratio'      :.24,
    'magnet_material_cost'       : 712756, # $/m3
    'magnetization_direction'    : 'Parallel',
    'B_r'                        : 1.285, # Tesla, magnet residual flux density
    'mu_r'                       : 1.062, # magnet relative permeability
    'magnet_max_temperature'     : 80, # deg C
    'magnet_max_rad_stress'      : 0, # Mpa  
    'magnet_therm_conductivity'  : 8.95, # W/m-k
    }
Arnon5 = {
    'core_material'              : 'Arnon5', 
    'core_material_density'      : 7650, # kg/m3
    'core_youngs_modulus'        : 185E9, # Pa
    'core_poission_ratio'        : .3, 
    'core_material_cost'         : 17087, # $/m3
    'core_ironloss_a'            : 1.58, 
    'core_ironloss_b'            : 1.17, 
    'core_ironloss_Kh'           : 78.94, # W/m3
    'core_ironloss_Ke'           : 0.0372, # W/m3
    'core_therm_conductivity'    : 28, # W/m-k
    'core_stacking_factor'       : 96 # percentage
    }

M19Gauge29 = {
    'core_material'              : 'M19Gauge29',
    'core_material_density'      : 7650, # kg/m3
    'core_youngs_modulus'        : 185E9, # Pa
    'core_poission_ratio'        : .3,
    'core_material_cost'         : 17087, # $/m3
    #'core_ironloss_a'            : 1.58,
    #'core_ironloss_b'            : 1.17,
    #'core_ironloss_Kh'           : 78.94, # W/m3
    #'core_ironloss_Ke'           : 0.0372, # W/m3
    'core_ironloss_a'            : 1.8,# Field
    'core_ironloss_b'            : 0.99,# Freq
    'core_ironloss_Kh'           : 497.25, # W/m3
    'core_ironloss_Ke'           : 0.3465, # W/m3
    'core_therm_conductivity'    : 28, # W/m-k
    'core_stacking_factor'       : 96, # percentage
    'core_saturation_feild'      : 1.6 #T
    }
CoilCopper = {
    'Max_temp'                   : 150, # Rise C
    'k_ov'                       : 1.8,
    'sigma'                      : 5.80E7,
    'k_fill'                     : .38}
Ins = {
       'k':1}