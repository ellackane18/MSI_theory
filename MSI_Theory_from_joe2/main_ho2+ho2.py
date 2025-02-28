# -*- coding: utf-8 -*-
"""

Pipeline for PAPR-MESS MSI code.

@author: Lei Lei
"""

import postprocessor as MSI
import os
import numpy as np

# Initialize MSI.PAPR_MESS:
# -- the first arguement is the file path for the nominal PAPR-MESS input files
# -- the second arguement is the name for the unperturbed (nominal) input file
# -- the third arguement is a dictionary specifying the nominal condition (to allow optimization starting from nonzero perturbations),
#    the value specifies temperature list, pressure list, and nominal perturbation
# -- the fourth arguement is a dictionary with elements specifying the perturbations (sensitivity analysis are based on the nominal conditions defined above),
#    the keys specify names of perturbation runs and values specify temperature, pressure, and perturbation
# -  the fifth arguement is a list indicating the channels of interest

Temperature_list = np.arange(200,700,10).tolist() + np.arange(700,3100,100).tolist() # K
Pressure_list = [[0.1, 0.5, 1., 2., 4.], ['[atm]']] # Torr will be converted into atm internally

channels = ['P1->W3', 'P1->P2', 'P1->P3', 'P1->P4', 'P1->P5',
            'P2->W3', 'P2->P3', 'P2->P4', 'P2->P5',
            'P3->W3', 'P3->P4', 'P3->P5',
            'P4->W3', 'P4->P5',
            'P5->W3']# channel-specific rate constants of interest

channels = ['P1->W3','W3->P1', 'P1->P3', 'P1->P4', 'P1->P5',
            'P3->W3','W3->P3', 'P3->P4', 'P3->P5',
            'P4->W3','W3->P4', 'P4->P5',
            'P5->W3','W3->P5'] # channel-specific rate constants of interest

pertubation_percent = 0.1 # percentage of perturbation
nominal_MESS_input_path = os.getcwd() + '/ho2+ho2/'
nominal_MESS_input = 'singlet.inp'
#nominal_MESS_input = 'combined.inp'

#nominal_MESS_input = 'ho2+ho2_Lei.inp'


model = MSI.PAPR_MESS(nominal_MESS_input_path, nominal_MESS_input,
                                    {'W3_Energy_1': [Temperature_list, Pressure_list, 0.00]},
                                    {#'colrel_FactorOne': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'colrel_PowerOne': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'colrel_Epsilons': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'colrel_Sigmas': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'W2a_Energy_1': [Temperature_list, Pressure_list, pertubation_percent * 349.759],
                                    # 'W2a_SymmetryFactor_1': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'W2a_Frequencies_1': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'W2b_Energy_1': [Temperature_list, Pressure_list, pertubation_percent * 349.759],
                                    # 'W2b_SymmetryFactor_1': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'W2b_Frequencies_1': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'W3_Energy_1': [Temperature_list, Pressure_list, pertubation_percent * 349.759],
                                    # 'W3_SymmetryFactor_1': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'W3_Frequencies_1': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'B3_SymmetryFactor_1': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'B3_Frequencies_1': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'B3_Energy_1' : [Temperature_list, Pressure_list, pertubation_percent * 349.759],
                                    # 'B3_ImaginaryFrequency_1': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'B4_SymmetryFactor_1': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'B4_Frequencies_1': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'B4_Energy_1': [Temperature_list, Pressure_list, pertubation_percent * 349.759],
                                    # 'B5_SymmetryFactor_1': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'B5_Frequencies_1': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'B5_Energy_1': [Temperature_list, Pressure_list, pertubation_percent * 349.759],
                                    # 'B6_SymmetryFactor_1': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'B6_Frequencies_1': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'B6_Energy_1' : [Temperature_list, Pressure_list, pertubation_percent * 349.759],
                                    # 'B6_ImaginaryFrequency_1': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'B7_SymmetryFactor_1': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'B7_Frequencies_1': [Temperature_list, Pressure_list, pertubation_percent],
                                     'B7_Energy_1' : [Temperature_list, Pressure_list, pertubation_percent * 349.759],
                                    # 'B9_SymmetryFactor_1': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'B9_Frequencies_1': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'B9_Energy_1' : [Temperature_list, Pressure_list, pertubation_percent * 349.759],

                                    # 'P1_SymmetryFactor_1': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'P1_Frequencies_1': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'P1_Energy_1' : [Temperature_list, Pressure_list, pertubation_percent * 349.759],

                                    # 'P3_SymmetryFactor_1': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'P3_Frequencies_1': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'P3_Energy_1' : [Temperature_list, Pressure_list, pertubation_percent * 349.759],

                                    # 'P4_SymmetryFactor_1': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'P4_Frequencies_1': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'P4_Energy_1' : [Temperature_list, Pressure_list, pertubation_percent * 349.759],

                                    # 'P5_SymmetryFactor_1': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'P5_Frequencies_1': [Temperature_list, Pressure_list, pertubation_percent],
                                    # 'P5_Energy_1' : [Temperature_list, Pressure_list, pertubation_percent * 349.759],                                    

                                    },
                                    channels)

# fit rate constants
n_P = 3  # number of pressure polynominals
n_T = 15 # number of temperature polynominals
same_line_result = True
model.Run() # execute PAPR-MESS perturbations
model.fit_Cheb_rates(n_P, n_T, P_min=0.0001, P_max=10.0, T_min=200.0, T_max=3000.0, same_line_result=same_line_result) # fit rate constants into Chebyshev expressions
model.Cheb_sens_coeff(same_line_result=same_line_result) # calculate sensitivity coefficients


# for debugging, fitting rate constants for a given trail
# target_dir = '/home/leil/MSI/PAPR-MESS_calculation/h+ho2=hooh/calculation_5/'
# model.fit_Cheb_rates(n_P, n_T, target_dir=target_dir)
# model.Cheb_sens_coeff()
