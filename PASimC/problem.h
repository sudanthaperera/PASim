/*
 * =====================================================================================
 *
 *       Filename:  problem.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  10/20/2014 12:56:39 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Sudantha Perera, 
 *   Organization:  OU-RIL
 *
 * =====================================================================================
 */

#include "pasdef.h"

// Recording data
#define SIM_EN TRUE
#define GEO_REC_EN TRUE
#define E_REC_EN TRUE
#define H_REC_EN TRUE
#define V_REC_EN TRUE
#define I_REC_EN TRUE
#define VS_REC_EN TRUE
#define IS_REC_EN TRUE
#define S_REC_EN TRUE

//defining the problem space parameters
#define NUMBER_OF_TIME_STEPS 1000

#define COURANT_FACTOR 0.9

#define NUMBER_OF_CELLS_PER_WAVELENGTH 20

#define DELTA_X 0.5e-3
#define DELTA_Y 0.5e-3
#define DELTA_Z 0.5e-3

#define BOUNDARY_TYPE_XN CPML
#define BOUNDARY_TYPE_XP CPML
#define BOUNDARY_TYPE_YN CPML
#define BOUNDARY_TYPE_YP CPML
#define BOUNDARY_TYPE_ZN CPML
#define BOUNDARY_TYPE_ZP CPML

#define BOUNDARY_AIR_BUFFER_NUMBER_OF_CELLS_XN 10
#define BOUNDARY_AIR_BUFFER_NUMBER_OF_CELLS_XP 10
#define BOUNDARY_AIR_BUFFER_NUMBER_OF_CELLS_YN 10
#define BOUNDARY_AIR_BUFFER_NUMBER_OF_CELLS_YP 10
#define BOUNDARY_AIR_BUFFER_NUMBER_OF_CELLS_ZN 10
#define BOUNDARY_AIR_BUFFER_NUMBER_OF_CELLS_ZP 10

#define BOUNDARY_CPML_NUMBER_OF_CELLS_XN 8
#define BOUNDARY_CPML_NUMBER_OF_CELLS_XP 8
#define BOUNDARY_CPML_NUMBER_OF_CELLS_YN 8
#define BOUNDARY_CPML_NUMBER_OF_CELLS_YP 8
#define BOUNDARY_CPML_NUMBER_OF_CELLS_ZN 8
#define BOUNDARY_CPML_NUMBER_OF_CELLS_ZP 8

#define BOUNDARY_CPML_ORDER 3
#define BOUNDARY_CPML_SIGMA_FACTOR 1.3
#define BOUNDARY_CPML_KAPPA_MAX 7
#define BOUNDARY_CPML_ALPHA_MIN 0
#define BOUNDARY_CPML_ALPHA_MAX 0.05

#define MATERIAL_TYPE_INDEX_AIR 1
#define MATERIAL_TYPE_INDEX_PEC 2
#define MATERIAL_TYPE_INDEX_PMC 3
#define MATERIAL_TYPE_INDEX_RT5880 4

#define NUMBER_OF_BRICKS 4

#define AIR_EPS_R 1
#define AIR_MU_R 1
#define AIR_SIGMA_E 0
#define AIR_SIGMA_M 0

#define PEC_EPS_R 1
#define PEC_MU_R 1
#define PEC_SIGMA_E 1e+10
#define PEC_SIGMA_M 0

#define PMC_EPS_R 1
#define PMC_MU_R 1
#define PMC_SIGMA_E 0
#define PMC_SIGMA_M 1e+10

#define RT5880_EPS_R 2.2
#define RT5880_MU_R 1
#define RT5880_SIGMA_E 0
#define RT5880_SIGMA_M 0

//defining the problem geometry

#define AIR_X_MIN -40e-3
#define AIR_Y_MIN -40e-3
#define AIR_Z_MIN -40e-3
#define AIR_X_MAX 120e-3
#define AIR_Y_MAX 120e-3
#define AIR_Z_MAX 70e-3

#define RT5880_X_MIN 4e-3
#define RT5880_Y_MIN 4e-3
#define RT5880_Z_MIN 0
#define RT5880_X_MAX 76e-3
#define RT5880_Y_MAX 76e-3
#define RT5880_Z_MAX 1.5e-3

#define GND_X_MIN 4e-3
#define GND_Y_MIN 4e-3
#define GND_Z_MIN 0
#define GND_X_MAX 76e-3
#define GND_Y_MAX 76e-3
#define GND_Z_MAX 0

#define PATCH_X_MIN 25e-3
#define PATCH_Y_MIN 25e-3
#define PATCH_Z_MIN 1.5e-3
#define PATCH_X_MAX 55e-3
#define PATCH_Y_MAX 55e-3
#define PATCH_Z_MAX 1.5e-3

//defining sources and lumped element components

#define WAVEFORMS_GAUSSIAN_NUMBER_OF_CELLS_PER_WAVELENGTH 15;

#define VOLTAGE_SOURCES_X_MIN 40e-3
#define VOLTAGE_SOURCES_Y_MIN 44.5e-3
#define VOLTAGE_SOURCES_Z_MIN 0
#define VOLTAGE_SOURCES_X_MAX 40e-3
#define VOLTAGE_SOURCES_Y_MAX 45.5e-3
#define VOLTAGE_SOURCES_Z_MAX 1.5e-3
#define VOLTAGE_SOURCES_DIRECTION ZP
#define VOLTAGE_SOURCES_RESISTANCE 50
#define VOLTAGE_SOURCES_MAGNITUDE 1
#define VOLTAGE_SOURCES_WAVEFORM_TYPE GAUSSIAN
#define VOLTAGE_SOURCES_WAVEFORM_INDEX 1

#define RUN_SIMULATION TRUE
#define SAVE_MESH TRUE
#define SAVE_SPACE TRUE

#define FARFIELD_FREQUENCIES 3.2e+9
#define FARFIELS_NUMBER_OF_CELLS_FROM_OUTER_BOUNDARY 13

#define FREQUENCY_DOMAIN_START 1.0e+8
#define FREQUENCY_DOMAIN_END 5.0e+9
#define FREQUENCY_DOMAIN_STEP 10e+6

#define SAMPLED_VOLTAGE_X_MIN 40e-3
#define SAMPLED_VOLTAGE_Y_MIN 44.5e-3
#define SAMPLED_VOLTAGE_Z_MIN 0
#define SAMPLED_VOLTAGE_X_MAX 40e-3
#define SAMPLED_VOLTAGE_Y_MAX 45.5e-3
#define SAMPLED_VOLTAGE_Z_MAX 1.5e-3
#define SAMPLED_VOLTAGE_DIRECTION ZP

#define SAMPLED_CURRENT_X_MIN 40e-3
#define SAMPLED_CURRENT_Y_MIN 44.5e-3
#define SAMPLED_CURRENT_Z_MIN 0
#define SAMPLED_CURRENT_X_MAX 40e-3
#define SAMPLED_CURRENT_Y_MAX 45.5e-3
#define SAMPLED_CURRENT_Z_MAX 1.5e-3
#define SAMPLED_CURRENT_DIRECTION ZP

#define PORT_SAMPLED_VOLTAGE_INDEX 1
#define PORT_SAMPLED_CURRENT_INDEX 1
#define PORT_IMPEDANCE 50
#define PORT_SOURCE 1


