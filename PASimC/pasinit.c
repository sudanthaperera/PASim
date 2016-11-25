/*
 * =====================================================================================
 *
 *       Filename:  passinit.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  09/09/2014 12:49:31 PM
 *       Compiler:  gcc
 *
 *         Author:  Sudantha Perera 
 *   Organization:  OU-RIL
 *
 * =====================================================================================
 */
#include<stdio.h>
#include "pasinit.h"
#include "problem.h"

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  calculate_domain_size
 *  Description:  
 * =====================================================================================
 */

void calculate_domain_size(int *nx, int *ny, int *nz){
	printf("Calculating the number of cells in the problem space...\n");
	
	double fdtd_domain_x_min = AIR_X_MIN - DELTA_X*BOUNDARY_AIR_BUFFER_NUMBER_OF_CELLS_XN;
	double fdtd_domain_y_min = AIR_Y_MIN - DELTA_Y*BOUNDARY_AIR_BUFFER_NUMBER_OF_CELLS_YN;
	double fdtd_domain_z_min = AIR_Z_MIN - DELTA_Z*BOUNDARY_AIR_BUFFER_NUMBER_OF_CELLS_ZN;
	double fdtd_domain_x_max = AIR_X_MAX - DELTA_X*BOUNDARY_AIR_BUFFER_NUMBER_OF_CELLS_XP;
	double fdtd_domain_y_max = AIR_Y_MAX - DELTA_Y*BOUNDARY_AIR_BUFFER_NUMBER_OF_CELLS_YP;
	double fdtd_domain_z_max = AIR_Z_MAX - DELTA_Z*BOUNDARY_AIR_BUFFER_NUMBER_OF_CELLS_ZP;

	double fdtd_domain_size_x = fdtd_domain_x_max - fdtd_domain_x_min;
	double fdtd_domain_size_y = fdtd_domain_y_max - fdtd_domain_y_min;
	double fdtd_domain_size_z = fdtd_domain_z_max - fdtd_domain_z_min;

	(*nx) = fdtd_domain_size_x/DELTA_X;	
	(*ny) = fdtd_domain_size_y/DELTA_Y;
	(*nz) = fdtd_domain_size_z/DELTA_Z;

	fdtd_domain_size_x = (*nx)*DELTA_X;
	fdtd_domain_size_y = (*nx)*DELTA_Y;
	fdtd_domain_size_z = (*nx)*DELTA_Z;
	
	fdtd_domain_x_max = fdtd_domain_x_min + fdtd_domain_size_x;	
	fdtd_domain_y_max = fdtd_domain_y_min + fdtd_domain_size_y;
	fdtd_domain_z_max = fdtd_domain_z_min + fdtd_domain_size_z;

	double fdtd_domain_cell_center_x[(*nx)][(*ny)][(*nz)];
	double fdtd_domain_cell_center_x[(*nx)][(*ny)][(*nz)];
	double fdtd_domain_cell_center_x[(*nx)][(*ny)][(*nz)];

	int i, j, k;
	for (i=1;i<(*nx);i++){
		fdtd_domain_cell_center_x[i][][] = (i - 0.5)*DELTA_X + fdtd_domain_x_min; 
	}
	for (j=1;j<(*ny);j++){
		fdtd_domain_cell_center_x[][j][] = (j - 0.5)*DELTA_Y + fdtd_domain_y_min;
	}
	for(k=1;k<(*nz);k++){
		fdtd_domain_cell_center_x[][][k] = (k - 0.5)*DELTA_Z + fdtd_domain_z_min;
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  create_bricks
 *  Description:  
 * =====================================================================================
 */

void create_bricks(void);{

}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  init_eps_array
 *  Description:  
 * =====================================================================================
 */


void init_eps_array(double ***eps_r_x, double ***eps_r_y, double ***eps_r_z, int nx, int ny, int nz){

}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  calculate_domain_size
 *  Description:  
 * =====================================================================================
 */


void init_mu_array(double ***mu_r_x, double ***mu_r_y, double ***mu_r_z, int nx, int ny, int nz){

}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  calculate_domain_size
 *  Description:  
 * =====================================================================================
 */


void int_sigma_e_array(double ***sigma_e_x, double ***sigma_e_y, double ***sigma_e_z, int nx, int ny, int nz){

}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  calculate_domain_size
 *  Description:  
 * =====================================================================================
 */


void int_sigma_m_array(double ***sigma_m_x, double ***sigma_m_y, double ***sigma_m_z, int nx, int ny, int nz){

}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  calculate_domain_size
 *  Description:  
 * =====================================================================================
 */


void calculate_eps_array(double ***eps_r_x, double ***eps_r_y, double ***eps_r_z, int nx, int ny, int nz){

}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  calculate_domain_size
 *  Description:  
 * =====================================================================================
 */


void calculate_mu_array(double ***mu_r_x, double ***mu_r_y, double ***mu_r_z, int nx, int ny, int nz){

}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  calculate_domain_size
 *  Description:  
 * =====================================================================================
 */


void calculate_sigma_e_array(double ***sigma_e_x, double ***sigma_e_y, double ***sigma_e_z, int nx, int ny, int nz){

}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  calculate_domain_size
 *  Description:  
 * =====================================================================================
 */


void calculate_sigma_m_array(double ***sigma_m_x, double ***sigma_m_y, double ***sigma_m_z, int nx, int ny, int nz){

}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  calculate_domain_size
 *  Description:  
 * =====================================================================================
 */


void create_PEC_plates(void){

}

