/*
 * =====================================================================================
 *
 *       Filename:  pasinit.h
 *
 *    Description:  (1) initializing FDTD material grid
 *						(a) calculate_domain_size
 *						(b) 
 *					(2) initializing FDTD parameters and arrays
 *					(3) initializing sources and lumped element components
 *					(4) initializing general updating coefficients
 *					(5) initialize boundary parameters
 *					(6) initializing the output parameters
 *					(7) initialize farfield arrays
 *					(8) 
 *
 *        Version:  1.0
 *        Created:  09/09/2014 11:24:01 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */

void calculate_domain_size(int *nx, int *ny, int *nz);
void create_bricks(void);

void init_eps_array(double ***eps_r_x, double ***eps_r_y, double ***eps_r_z, int nx, int ny, int nz);
void init_mu_array(double ***mu_r_x, double ***mu_r_y, double ***mu_r_z, int nx, int ny, int nz);
void int_sigma_e_array(double ***sigma_e_x, double ***sigma_e_y, double ***sigma_e_z, int nx, int ny, int nz);
void int_sigma_m_array(double ***sigma_m_x, double ***sigma_m_y, double ***sigma_m_z, int nx, int ny, int nz);

void calculate_eps_array(double ***eps_r_x, double ***eps_r_y, double ***eps_r_z, int nx, int ny, int nz);
void calculate_mu_array(double ***mu_r_x, double ***mu_r_y, double ***mu_r_z, int nx, int ny, int nz);
void calculate_sigma_e_array(double ***sigma_e_x, double ***sigma_e_y, double ***sigma_e_z, int nx, int ny, int nz);
void calculate_sigma_m_array(double ***sigma_m_x, double ***sigma_m_y, double ***sigma_m_z, int nx, int ny, int nz);

void create_PEC_plates(void);

