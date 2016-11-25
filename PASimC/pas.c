/*
 * =====================================================================================
 *
 *       Filename:  pas.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  9/17/2014 4:04:06 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#include "pastypes.h"
#include "pasdef.h"

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

double* build_1D_double_array(int size, double init_val);
double complex* build_1D_double_complex_array(int size, double complex init_val);
double complex** build_2D_double_complex_array(int n_x, int n_y, double complex init_val);
double complex**** build_4D_double_complex_array(int n_x, int n_y, int n_z, int n_u, double complex init_val);
double**** build_4D_double_array(int n_x, int n_y, int n_z, int n_u, double init_val);
double*** build_3D_double_array(int n_x, int n_y, int n_z, double init_val);
int*** build_3D_int_array(int n_x, int n_y, int n_z, int init_val);

double complex exp_complex(double complex number);

voltage_prob_t create_voltage_prob(double min_x, double min_y, double min_z, double max_x, double max_y, double max_z, int direction);
voltage_source_t create_voltage_source(double min_x, double min_y, double min_z, double max_x, double max_y, double max_z, int direction, double resistance, double magnitude, int waveform_type, int waveform_index);
brick_t create_brick(double min_x, double min_y, double min_z, double max_x, double max_y, double max_z, int material_type);
material_t create_material(int index, double eps_r, double mu_r, double sigma_e, double sigma_m);
void define_boundary(boundary_t *b, int air_buffer, int cpml);
double complex* time_to_frequency_domain(double *x, double dt, double *frequency_array,double time_shift, int number_of_time_steps,int number_of_frequencies);
void capture_sampled_voltages(int number_of_sampled_voltages, int time_step ,voltage_prob_t *sampled_voltages, double ***E_field_x, double ***E_field_y, double ***E_field_z);
void capture_electric_fields(int number_of_sampled_electric_fields, int time_step ,sampled_fields_t *sampled_electric_fields, double ***E_field_x, double ***E_field_y, double ***E_field_z);
void capture_sampled_currents(int number_of_samples, int time_step, current_prob_t *samples, double ***H_field_x, double ***H_field_y, double ***H_field_z, cell_t cell);
void out_of_memory_status(void* array);

void save_3D_int_array(char *file, int ***array, int size_x, int size_y, int size_z);
void save_1D_double_array(char *file, double *array, int size);
void save_3D_double_array(char *file, double ***array, int size_x, int size_y, int size_z);
void save_1D_complex_double_array(char *file, complex double *array, int size);
void save_2D_complex_double_array(char *file, complex double **array, int size_x, int size_y);

void fill_four_cell_average_eps_x(double*** array, material_t *mat, int*** index, int x_start, int x_end, int y_start, int y_end, int z_start, int z_end);
void fill_four_cell_average_eps_y(double*** array, material_t *mat, int*** index, int x_start, int x_end, int y_start, int y_end, int z_start, int z_end);
void fill_four_cell_average_eps_z(double*** array, material_t *mat, int*** index, int x_start, int x_end, int y_start, int y_end, int z_start, int z_end);

void fill_four_cell_average_sigma_x(double*** array,material_t *mat,int*** index, int x_start, int x_end, int y_start, int y_end, int z_start, int z_end);
void fill_four_cell_average_sigma_y(double*** array,material_t *mat,int*** index, int x_start, int x_end, int y_start, int y_end, int z_start, int z_end);
void fill_four_cell_average_sigma_z(double*** array,material_t *mat,int*** index, int x_start, int x_end, int y_start, int y_end, int z_start, int z_end);

void fill_two_cell_average_mu_x(double*** array, material_t *mat, int*** index, int x_start, int x_end, int y_start, int y_end, int z_start, int z_end);
void fill_two_cell_average_mu_y(double*** array, material_t *mat, int*** index, int x_start, int x_end, int y_start, int y_end, int z_start, int z_end);
void fill_two_cell_average_mu_z(double*** array, material_t *mat, int*** index, int x_start, int x_end, int y_start, int y_end, int z_start, int z_end);

void fill_two_cell_average_sigma_x(double*** array, material_t *mat, int*** index, int x_start, int x_end, int y_start, int y_end, int z_start, int z_end);
void fill_two_cell_average_sigma_y(double*** array, material_t *mat, int*** index, int x_start, int x_end, int y_start, int y_end, int z_start, int z_end);
void fill_two_cell_average_sigma_z(double*** array, material_t *mat, int*** index, int x_start, int x_end, int y_start, int y_end, int z_start, int z_end);
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  
 * =====================================================================================
 */

int main(){
	register unsigned int i,j,k;

    printf("Defining the problem space parameters...\n");
    int number_of_time_steps = 1000;
    double courant_factor =0.9;

    double dx = 0.5e-3; 
    double dy = 0.5e-3;
    double dz = 0.5e-3;

	cell_t unit_cell;
	unit_cell.dx = dx;
	unit_cell.dy = dy;
	unit_cell.dz = dz;

    boundary_t boundary;

	define_boundary(&boundary,10,8);

    material_t material_types[5];
 
    // vac
	material_types[0] = create_material(0, 1.0, 1.0, 1.0e-40, 1.0e-40);
    
    // air
	material_types[1] = create_material(1, 1.0, 1.0, 1.0e-40, 1.0e-40);

    // PEC : perfect electric conductor
	material_types[2] = create_material(2, 1.0, 1.0, 1.0e10, 1.0e-40);

    // PMC : perfect magnetic conductor
	material_types[3] = create_material(3, 1.0, 1.0, 1.0e-40, 1.0e10);

    // a dielectric 
	material_types[4] = create_material(4, 2.2, 1.0, 1.0e-40, 1.0e-40);

    printf("Defining the problem geometry...\n");

    brick_t bricks[4];

	//air
	bricks[0] = create_brick(-1e-3, -1e-3, -1e-3, 81e-3, 81e-3, 4e-3, 1);

	//dielectric
	bricks[1] = create_brick( 4e-3, 4e-3, 0, 76e-3, 76e-3, 1.5e-3, 4); 

	//ground
	bricks[2] = create_brick(4e-3, 4e-3, 0, 76e-3, 76e-3, 0, 2);
	
	//patch
	bricks[3] = create_brick(25e-3, 25e-3, 1.5e-3, 55e-3, 55e-3, 1.5e-3, 2);
  	
    printf("Defining sources and lumped element components...\n");
	int number_of_voltage_sources = 1;
	int number_of_current_sources = 0;
	
    voltage_source_t *voltage_sources;
	current_source_t *current_sources;

	if (number_of_voltage_sources > 0){
		voltage_sources = (voltage_source_t*)malloc(number_of_voltage_sources*sizeof(voltage_source_t));
		out_of_memory_status(voltage_sources);
	}

	if (number_of_current_sources > 0){
		current_sources = (current_source_t*)malloc(number_of_current_sources*sizeof(current_source_t));
		out_of_memory_status(current_sources);
	}

    printf("Defining source waveform types and parameters...\n");
	int number_of_gaussian_waveforms = 1;
	gaussian_t waveforms_gaussian[number_of_gaussian_waveforms];

    waveforms_gaussian[0].number_of_cells_per_wavelength = 20; 
    waveforms_gaussian[1].number_of_cells_per_wavelength = 25;

	waveforms_gaussian[0].waveform = build_1D_double_array(number_of_time_steps,0.0);
	waveforms_gaussian[1].waveform = build_1D_double_array(number_of_time_steps,0.0);

    //voltage sources
	if (number_of_voltage_sources > 0){
		voltage_sources[0] = create_voltage_source(40e-3, 44.5e-3, 0, 40e-3, 45.5e-3, 1.5e-3, ZP, 50, 1, GAUSSIAN, 1);
	}

    printf("Defining output parameters...\n");

	int number_of_sampled_electric_fields  = 0;
	int number_of_sampled_magnetic_fields  = 0;
	int number_of_sampled_voltages  = 1;
	int number_of_sampled_currents  = 1;
	int number_of_ports = 1;

    voltage_prob_t sampled_voltages[number_of_sampled_voltages];
    current_prob_t sampled_currents[number_of_sampled_currents];
    port_t ports[number_of_ports];
   
	farfield_t farfield;
	int number_of_farfield_frequencies = 5;
	farfield.frequencies = build_1D_double_array(number_of_farfield_frequencies,0.0);

    //far field calculation parameters
    farfield.frequencies[0] = 3.0e9;
    farfield.frequencies[1] = 3.1e9;
    farfield.frequencies[2] = 3.2e9;
    farfield.frequencies[3] = 3.3e9;
	farfield.frequencies[4] = 3.4e9;
    farfield.number_of_cells_from_outer_boundary = 13;

    //frequency domain parameters
	frequency_domain_t frequency_domain; 
    frequency_domain.start = 100e6;
    frequency_domain.end   = 5e9;
    frequency_domain.step  = 10e6;

    //define sampled voltages
	if (number_of_sampled_voltages > 0){
		sampled_voltages[0] = create_voltage_prob( 40e-3, 44.5e-3, 0, 40e-3, 45.5e-3, 1.5e-3, ZP);
	}

    //define sampled currents
    sampled_currents[0].min_x = 40e-3;
    sampled_currents[0].min_y = 44.5e-3;
    sampled_currents[0].min_z = 1e-3;
    sampled_currents[0].max_x = 40e-3;
    sampled_currents[0].max_y = 45.5e-3;
    sampled_currents[0].max_z = 1e-3;
    sampled_currents[0].direction = ZP;

    //define ports
    ports[0].sampled_voltage_index = 0;
    ports[0].sampled_current_index = 0;
    ports[0].impedance = 50;
    ports[0].is_source_port = TRUE;

	printf("Initializing FDTD material grid...\n");
	

	printf("Calculating the number of cells in the problem space...\n");

	int number_of_bricks  = 4;
	double fdtd_domain_min_x=0.0;
	double fdtd_domain_min_y=0.0;
	double fdtd_domain_min_z=0.0;
	double fdtd_domain_max_x=0.0;
	double fdtd_domain_max_y=0.0;
	double fdtd_domain_max_z=0.0;

	for(i=0;i<number_of_bricks;i++)
	{
		if(fdtd_domain_min_x > bricks[i].min_x)
			fdtd_domain_min_x = bricks[i].min_x; 
		if(fdtd_domain_min_y > bricks[i].min_y)
			fdtd_domain_min_y = bricks[i].min_y; 
		if(fdtd_domain_min_z > bricks[i].min_z)
			fdtd_domain_min_z = bricks[i].min_z; 
		if(fdtd_domain_max_x < bricks[i].max_x)
			fdtd_domain_max_x = bricks[i].max_x; 
		if(fdtd_domain_max_y < bricks[i].max_y)
			fdtd_domain_max_y = bricks[i].max_y; 
		if(fdtd_domain_max_z < bricks[i].max_z)
			fdtd_domain_max_z = bricks[i].max_z;
	} 

	printf("The problem space without air buffer\n\tmin_x = %f m\n\tmin_y = %f m\n\tmin_z = %f m\n\tmax_x = %f m\n\tmax_y = %f m\n\tmax_z = %f m\n",fdtd_domain_min_x,fdtd_domain_min_y,fdtd_domain_min_z,fdtd_domain_max_x,fdtd_domain_max_y,fdtd_domain_max_z);
	//Determine the problem space boundaries including air buffers 
	fdtd_domain_min_x = fdtd_domain_min_x - dx * boundary.air_buffer_number_of_cells_xn;
	fdtd_domain_min_y = fdtd_domain_min_y - dy * boundary.air_buffer_number_of_cells_yn;
	fdtd_domain_min_z = fdtd_domain_min_z - dz * boundary.air_buffer_number_of_cells_zn;
	fdtd_domain_max_x = fdtd_domain_max_x + dx * boundary.air_buffer_number_of_cells_xp;
	fdtd_domain_max_y = fdtd_domain_max_y + dy * boundary.air_buffer_number_of_cells_yp;
	fdtd_domain_max_z = fdtd_domain_max_z + dz * boundary.air_buffer_number_of_cells_zp;

	printf("The problem space with air buffer\n\tmin_x = %fm\n\tmin_y = %fm\n\tmin_z = %fm\n\tmax_x = %fm\n\tmax_y = %fm\n\tmax_z = %fm\n",fdtd_domain_min_x,fdtd_domain_min_y,fdtd_domain_min_z,fdtd_domain_max_x,fdtd_domain_max_y,fdtd_domain_max_z);
	//Determining the problem space size
	double fdtd_domain_size_x = fdtd_domain_max_x - fdtd_domain_min_x;
	double fdtd_domain_size_y = fdtd_domain_max_y - fdtd_domain_min_y;
	double fdtd_domain_size_z = fdtd_domain_max_z - fdtd_domain_min_z;

	printf("Problem Space size in x direction %fm\n",fdtd_domain_size_x);
	printf("Problem Space size in y direction %fm\n",fdtd_domain_size_y);
	printf("Problem Space size in z direction %fm\n",fdtd_domain_size_z);

	printf("Calculating number of cells in x, y, and z directions...\n");
	int nx = (int)round(fdtd_domain_size_x/dx);  
	int ny = (int)round(fdtd_domain_size_y/dy);
	int nz = (int)round(fdtd_domain_size_z/dz);

	printf("Calculating adjust domain size by snapping to cells...\n");
	fdtd_domain_size_x = nx * dx;
	fdtd_domain_size_y = ny * dy;
	fdtd_domain_size_x = nz * dz;

	fdtd_domain_max_x = fdtd_domain_min_x + fdtd_domain_size_x;
	fdtd_domain_max_y = fdtd_domain_min_y + fdtd_domain_size_y;
	fdtd_domain_max_z = fdtd_domain_min_z + fdtd_domain_size_z;

	//some frequently used auxiliary parameters 
	unsigned int nxp1 = nx+1, nyp1 = ny+1, nzp1 = nz+1, nxm1 = nx-1, nxm2 = nx-2, nym1 = ny-1, nym2 = ny-2, nzm1 = nz-1, nzm2 = nz-2;
	
	printf("Number of cells in X direction = %d\n",nx);
	printf("Number of cells in Y direction = %d\n",ny);
	printf("Number of cells in Z direction = %d\n",nz);
	printf("Number of cells = %d \n", nx*ny*nz);
	printf("Creating arrays storing the center coordinates of the cells...\n");
	
	double ***fdtd_domain_cell_center_coordinates_x;
	fdtd_domain_cell_center_coordinates_x = build_3D_double_array(nx,ny,nz,0.0);
	double ***fdtd_domain_cell_center_coordinates_y;
	fdtd_domain_cell_center_coordinates_y = build_3D_double_array(nx,ny,nz,0.0);
	double ***fdtd_domain_cell_center_coordinates_z;
	fdtd_domain_cell_center_coordinates_z = build_3D_double_array(nx,ny,nz,0.0);
	int ***material_3d_space;
	material_3d_space = build_3D_int_array(nx,ny,nz,1);

	for(k = 0;k < nz;k++){
		for(j = 0;j < ny;j++){
			for(i = 0;i < nx;i++){
				fdtd_domain_cell_center_coordinates_x[i][j][k] = (i+1 - 0.5) * dx + fdtd_domain_min_x;
			}
		}
	}

	for(i = 0;i < nx;i++){
		for(k = 0;k < nz;k++){
			for(j = 0;j < ny;j++){
				fdtd_domain_cell_center_coordinates_y[i][j][k] = (j+1 - 0.5) * dy + fdtd_domain_min_y;
			}
		}
	}

	for(i = 0;i < nx;i++){
		for(j = 0;j < ny;j++){
			for(k = 0;k < nz;k++){
				fdtd_domain_cell_center_coordinates_z[i][j][k] = (k+1 - 0.5) * dz + fdtd_domain_min_z;
			}
		}
	}
	
	printf("Creating bricks...\n");
	
	int blx,bly,blz,bux,buy,buz;
	int ind;
	for(ind = 1;ind < number_of_bricks;ind++){
		//convert brick end coordinates to node indices 
		blx = (int)round((bricks[ind].min_x - fdtd_domain_min_x)/dx); 
		bly = (int)round((bricks[ind].min_y - fdtd_domain_min_y)/dy); 
		blz = (int)round((bricks[ind].min_z - fdtd_domain_min_z)/dz); 

		bux = (int)round((bricks[ind].max_x - fdtd_domain_min_x)/dx); 
		buy = (int)round((bricks[ind].max_y - fdtd_domain_min_y)/dy); 
		buz = (int)round((bricks[ind].max_z - fdtd_domain_min_z)/dz); 

		//assign material type of the brick to the cells
		int i,j,k;
		for(i=blx;i < bux;i++){
			for(j=bly;j < buy;j++){
				for(k=blz;k < buz;k++){
					material_3d_space[i][j][k] = bricks[ind].material_type;
				}
			}
		}
	}

	save_3D_int_array("geometry_file.csv", material_3d_space, nx, ny, nz);

	printf("Creating material component arrays for a problem space...\n");

	double ***eps_r_x;
	eps_r_x = build_3D_double_array(nx,nyp1,nzp1,1.0);

	double ***eps_r_y;
	eps_r_y = build_3D_double_array(nxp1,ny,nzp1,1.0);

	double ***eps_r_z;
	eps_r_z = build_3D_double_array(nxp1,nyp1,nz,1.0);

	double ***mu_r_x;
	mu_r_x = build_3D_double_array(nxp1,ny,nz,1.0);

	double ***mu_r_y;
	mu_r_y = build_3D_double_array(nx,nyp1,nz,1.0);

	double ***mu_r_z;
	mu_r_z = build_3D_double_array(nx,ny,nzp1,1.0);

	double ***sigma_e_x;
	sigma_e_x = build_3D_double_array(nx, nyp1, nzp1, 0.0);

	double ***sigma_e_y;
	sigma_e_y = build_3D_double_array(nxp1, ny, nzp1, 0.0);

	double ***sigma_e_z;
	sigma_e_z = build_3D_double_array(nxp1, nyp1, nz, 0.0);

	double ***sigma_m_x;
	sigma_m_x = build_3D_double_array(nxp1, ny, nz, 0.0);

	double ***sigma_m_y;
	sigma_m_y = build_3D_double_array(nx, nyp1, nz, 0.0);
	
	double ***sigma_m_z;
	sigma_m_z = build_3D_double_array(nx, ny, nzp1, 0.0);

	printf("Filling material components arrays...\n");

	printf("\tCalculating eps_r_x...\n");
	fill_four_cell_average_eps_x(eps_r_x, material_types, material_3d_space, 0, nx, 1, ny, 1, nz);

	printf("\tCalculating eps_r_y...\n");
	fill_four_cell_average_eps_y(eps_r_y, material_types, material_3d_space, 1, nx, 0, ny, 1, nz);
                    
	printf("\tCalculating eps_r_z...\n");
	fill_four_cell_average_eps_z(eps_r_z, material_types, material_3d_space, 1, nx, 1, ny, 0, nz);

	printf("\tCalculating sigma_e_x...\n");
	fill_four_cell_average_sigma_x(sigma_e_x,material_types,material_3d_space, 0, nx, 1, ny, 1, nz);

	printf("\tCalculating sigma_e_y...\n");
	fill_four_cell_average_sigma_y(sigma_e_y,material_types,material_3d_space, 1, nx, 0, ny, 1, nz);
                    
	printf("\tCalculating sigma_e_z...\n");
	fill_four_cell_average_sigma_z(sigma_e_z,material_types,material_3d_space, 1, nx, 1, ny, 0, nz);
                    
	printf("\tCalculating mu_r_x...\n");
	fill_two_cell_average_mu_x(mu_r_x,material_types,material_3d_space, 1, nx, 0, ny, 0, nz);
                    
	printf("\tCalculating mu_r_y...\n");
	fill_two_cell_average_mu_y(mu_r_y,material_types,material_3d_space, 0, nx, 1, ny, 0, nz);
                    
	printf("\tCalculating mu_r_z...\n");
	fill_two_cell_average_mu_z(mu_r_z,material_types,material_3d_space, 0, nx, 0, ny, 1, nz);
                    
	printf("\tCalculating sigma_m_x...\n");
	fill_two_cell_average_sigma_x(sigma_m_x,material_types,material_3d_space, 1, nx, 0, ny, 0, nz);
                    
	printf("\tCalculating sigma_m_y...\n");
	fill_two_cell_average_sigma_y(sigma_m_y,material_types,material_3d_space, 0, nx, 1, ny, 0, nz);
                    
	printf("\tCalculating sigma_m_z...\n");
	fill_two_cell_average_sigma_z(sigma_m_z,material_types,material_3d_space, 0, nx, 0, ny, 1, nz);

	printf("Creating PEC plates on the material grid...\n");

	double sigma_pec;
	for(ind = 0;ind < number_of_bricks;ind++){
		sigma_pec = material_types[bricks[ind].material_type].sigma_e;
	
		// convert coordinates to node indices on the FDTD grid
		blx = (int)round((bricks[ind].min_x - fdtd_domain_min_x)/dx); 
		bly = (int)round((bricks[ind].min_y - fdtd_domain_min_y)/dy); 
		blz = (int)round((bricks[ind].min_z - fdtd_domain_min_z)/dz); 

		bux = (int)round((bricks[ind].max_x - fdtd_domain_min_x)/dx); 
		buy = (int)round((bricks[ind].max_y - fdtd_domain_min_y)/dy); 
		buz = (int)round((bricks[ind].max_z - fdtd_domain_min_z)/dz); 

		//find the zero thickness bricks
		int index_bx, index_by, index_bz;
		if (blx == bux)
		{
			for(index_by = bly;index_by <= buy;index_by++){
				for(index_bz = blz;index_bz <= buz;index_bz++){
					if(index_by < buy){  
						sigma_e_y[blx][index_by][index_bz] = sigma_pec;
					}
					if(index_bz < buz){	
						sigma_e_z[blx][index_by][index_bz] = sigma_pec;
					}
				}
			}
			
		}
		if (bly == buy)
		{
			for(index_bx = blx;index_bx <= bux;index_bx++){
				for(index_bz = blz;index_bz <= buz;index_bz++){
					if(index_bz < buz){
						sigma_e_z[index_bx][bly][index_bz] = sigma_pec;
					}
					if(index_bx < bux){
						sigma_e_x[index_bx][bly][index_bz] = sigma_pec;
					}
				}
			}
		}
		if (blz == buz)
		{
			for(index_bx = blx;index_bx <= bux;index_bx++){
				for(index_by = bly;index_by <= buy;index_by++){
					if(index_bx < bux){
						sigma_e_x[index_bx][index_by][blz] = sigma_pec;
					}
					if(index_by < buy){	
						sigma_e_y[index_bx][index_by][blz] = sigma_pec;
					}
				}
			}
		}
	}

	printf("Initializing FDTD parameters and arrays...\n");

	//constant parameters
	double eps_0 = 8.854187817620e-12;               
	double mu_0  = 4*M_PI*1e-7;                   
	double c = 1/sqrt(mu_0*eps_0);
	double eta_0 = sqrt(mu_0/eps_0);

	printf("Speed of light in free space = %g\n",c);
	printf("Intrinsic impedance of free space = %g\n",eta_0);

	//Duration of a time step in seconds
	double dt = courant_factor/(c*sqrt((1/(dx*dx))+(1/(dy*dy))+(1/(dz*dz))));

	printf("Duration of a time step in seconds = %g\n",dt);
	//time array
	double *time;
	time = build_1D_double_array(number_of_time_steps,0.0);
	
	int time_ind;
	for(time_ind=0;time_ind<number_of_time_steps;time_ind++){
		time[time_ind] = (0.5 + time_ind)*dt;
	}

	//Create and initialize field and current arrays
	printf("Creating field arrays...\n");
	
	double*** Hx;
	Hx = build_3D_double_array(nxp1,ny,nz,0.0);
	double*** Hy;   
	Hy = build_3D_double_array(nx,nyp1,nz,0.0);
	double*** Hz;   
	Hz = build_3D_double_array(nx,ny,nzp1,0.0);
	double*** Ex;
	Ex = build_3D_double_array(nx,nyp1,nzp1,0.0);
	double*** Ey;
	Ey = build_3D_double_array(nxp1,ny,nzp1,0.0);
	double*** Ez;
	Ez = build_3D_double_array(nxp1,nyp1,nz,0.0);
	
	printf("Initializing sources and lumped element components...\n");
	

	for (ind = 0; ind < number_of_gaussian_waveforms ; ind++){
		if (waveforms_gaussian[ind].number_of_cells_per_wavelength == 0){
			waveforms_gaussian[ind].number_of_cells_per_wavelength = DEFAULT_NUMBER_OF_CELLS_PER_WAVELENGTH;
		}

		waveforms_gaussian[ind].maximum_frequency = c/(waveforms_gaussian[ind].number_of_cells_per_wavelength*MAX(dx,MAX(dy,dz)));
		waveforms_gaussian[ind].tau = (waveforms_gaussian[ind].number_of_cells_per_wavelength*MAX(dx,MAX(dy,dz)))/(2*c);
		waveforms_gaussian[ind].t_0 = 4.5*waveforms_gaussian[ind].tau;

		for (time_ind = 0;time_ind < number_of_time_steps;time_ind++){
			waveforms_gaussian[ind].waveform[time_ind] = exp(-((time[time_ind] - waveforms_gaussian[ind].t_0)/waveforms_gaussian[ind].tau)* ((time[time_ind] - waveforms_gaussian[ind].t_0)/waveforms_gaussian[ind].tau));
		}
	}

	printf("initializing source waveforms...\n");
	
	int voltage_index;
	int is,js,ks,ie,je,ke;
	for(voltage_index = 0;voltage_index<number_of_voltage_sources;voltage_index++){
		is = (int)round((voltage_sources[voltage_index].min_x - fdtd_domain_min_x)/dx);
		js = (int)round((voltage_sources[voltage_index].min_y - fdtd_domain_min_y)/dy);
		ks = (int)round((voltage_sources[voltage_index].min_z - fdtd_domain_min_z)/dz);
		ie = (int)round((voltage_sources[voltage_index].max_x - fdtd_domain_min_x)/dx);
		je = (int)round((voltage_sources[voltage_index].max_y - fdtd_domain_min_y)/dy);
		ke = (int)round((voltage_sources[voltage_index].max_z - fdtd_domain_min_z)/dz);
    
		voltage_sources[voltage_index].is = is;
		voltage_sources[voltage_index].js = js;
		voltage_sources[voltage_index].ks = ks;
		voltage_sources[voltage_index].ie = ie;
		voltage_sources[voltage_index].je = je;
		voltage_sources[voltage_index].ke = ke;
		
		int n_fields;
		double r_magnitude_factor;
		double v_magnitude_factor;
	
		if(voltage_sources[voltage_index].direction == XN){
			n_fields = ie - is;
			r_magnitude_factor = ((1 + je - js) * (1 + ke - ks))/(ie - is); 
			v_magnitude_factor = (-1*voltage_sources[voltage_index].magnitude)/n_fields;
		}
		else if(voltage_sources[voltage_index].direction == XP){
			n_fields = ie - is;
			r_magnitude_factor = ((1 + je - js) * (1 + ke - ks))/(ie - is); 
			v_magnitude_factor = (1*voltage_sources[voltage_index].magnitude)/n_fields;
		}
		else if(voltage_sources[voltage_index].direction == YN){
			n_fields = je - js;
			r_magnitude_factor = ((1 + ie - is) * (1 + ke - ks))/(je - js); 
			v_magnitude_factor = (-1*voltage_sources[voltage_index].magnitude)/n_fields;
		}
		else if(voltage_sources[voltage_index].direction == YP){
			n_fields = je - js;
			r_magnitude_factor = ((1 + ie - is) * (1 + ke - ks))/(je - js); 
			v_magnitude_factor = (1*voltage_sources[voltage_index].magnitude)/n_fields;
		}
		else if(voltage_sources[voltage_index].direction == ZN){
			n_fields = ke - ks;
			r_magnitude_factor = ((1 + ie - is) * (1 + je - js))/(ke - ks); 
			v_magnitude_factor = (-1*voltage_sources[voltage_index].magnitude)/n_fields;
		}
		else if(voltage_sources[voltage_index].direction == ZP){
			n_fields = ke - ks;
			r_magnitude_factor = ((1 + ie - is) * (1 + je - js))/(ke - ks); 
			v_magnitude_factor = (1*voltage_sources[voltage_index].magnitude)/n_fields;
		}

		voltage_sources[voltage_index].resistance_per_component = r_magnitude_factor * voltage_sources[voltage_index].resistance;
    
		//copy waveform of the waveform type to waveform of the source
		
		voltage_sources[voltage_index].voltage_per_e_field = build_1D_double_array(number_of_time_steps,0.0);
		voltage_sources[voltage_index].waveform = build_1D_double_array(number_of_time_steps,0.0);

		for (time_ind=0;time_ind<number_of_time_steps;time_ind++){
			voltage_sources[voltage_index].voltage_per_e_field[time_ind] = (v_magnitude_factor)*(waveforms_gaussian[0].waveform[time_ind]);
			voltage_sources[voltage_index].waveform[time_ind] = (v_magnitude_factor)*(waveforms_gaussian[0].waveform[time_ind])*(n_fields);
		}
	}

/////	
	printf("initializing general updating coefficients....\n");
	
	double ***Cexe;
	Cexe = build_3D_double_array(nx,nyp1,nzp1,0.0);
	double ***Cexhz;
	Cexhz = build_3D_double_array(nx,nyp1,nzp1,0.0);
	double ***Cexhy;
	Cexhy = build_3D_double_array(nx,nyp1,nzp1,0.0);

	double ***Ceye;
	Ceye = build_3D_double_array(nxp1,ny,nzp1,0.0);
	double ***Ceyhx;
	Ceyhx = build_3D_double_array(nxp1,ny,nzp1,0.0);
	double ***Ceyhz;
	Ceyhz = build_3D_double_array(nxp1,ny,nzp1,0.0);
	
	double ***Ceze;
	Ceze = build_3D_double_array(nxp1, nyp1, nz,0.0);
	double ***Cezhy;
	Cezhy = build_3D_double_array(nxp1, nyp1, nz,0.0);
	double ***Cezhx;
	Cezhx = build_3D_double_array(nxp1, nyp1, nz,0.0);

	double ***Chxh;
	Chxh = build_3D_double_array(nxp1, ny, nz,0.0);
	double ***Chxez;
	Chxez = build_3D_double_array(nxp1, ny, nz,0.0);
	double ***Chxey;
	Chxey = build_3D_double_array(nxp1, ny, nz,0.0);
	
	double ***Chyh;
	Chyh = build_3D_double_array(nx, nyp1, nz,0.0);
	double ***Chyex;
	Chyex = build_3D_double_array(nx, nyp1, nz,0.0);
	double ***Chyez;
	Chyez = build_3D_double_array(nx, nyp1, nz,0.0);
	
	double ***Chzh;
	Chzh = build_3D_double_array(nx, ny, nzp1,0.0);
	double ***Chzey;
	Chzey = build_3D_double_array(nx, ny, nzp1,0.0);
	double ***Chzex;
	Chzex = build_3D_double_array(nx, ny, nzp1,0.0);


	// General electric field updating coefficients
	for(i=0;i<nx;i++){
		for(j=0;j<nyp1;j++){
			for(k=0;k<nzp1;k++){
				// Coeffiecients updating Ex
				Cexe[i][j][k]  =  (2*eps_r_x[i][j][k]*eps_0 - dt*sigma_e_x[i][j][k])/(2*eps_r_x[i][j][k]*eps_0 + dt*sigma_e_x[i][j][k]);
				Cexhz[i][j][k]=  (2*dt/dy)/(2*eps_r_x[i][j][k]*eps_0 + dt*sigma_e_x[i][j][k]);
				Cexhy[i][j][k]= -(2*dt/dz)/(2*eps_r_x[i][j][k]*eps_0 + dt*sigma_e_x[i][j][k]);
			}
		}
	}

	for(i=0;i<nxp1;i++){
		for(j=0;j<ny;j++){
			for(k=0;k<nzp1;k++){
				//Coeffiecients updating Ey
				Ceye[i][j][k]  =  (2*eps_r_y[i][j][k]*eps_0 - dt*sigma_e_y[i][j][k])/(2*eps_r_y[i][j][k]*eps_0 + dt*sigma_e_y[i][j][k]);
				Ceyhx[i][j][k] =  (2*dt/dz)/(2*eps_r_y[i][j][k]*eps_0 + dt*sigma_e_y[i][j][k]);
				Ceyhz[i][j][k] = -(2*dt/dx)/(2*eps_r_y[i][j][k]*eps_0 + dt*sigma_e_y[i][j][k]);
			}
		}
	}

	for(i=0;i<nxp1;i++){
		for(j=0;j<nyp1;j++){
			for(k=0;k<nz;k++){
				//Coeffiecients updating Ez
				Ceze[i][j][k]  =  (2*eps_r_z[i][j][k]*eps_0 - dt*sigma_e_z[i][j][k])/(2*eps_r_z[i][j][k]*eps_0 + dt*sigma_e_z[i][j][k]);
				Cezhy[i][j][k] =  (2*dt/dx)/(2*eps_r_z[i][j][k]*eps_0 + dt*sigma_e_z[i][j][k]);
				Cezhx[i][j][k] = -(2*dt/dy)/(2*eps_r_z[i][j][k]*eps_0 + dt*sigma_e_z[i][j][k]);
			}
		}
	}

	//General magnetic field updating coefficients
	for(i=0;i<nxp1;i++){
		for(j=0;j<ny;j++){
			for(k=0;k<nz;k++){
				//Coeffiecients updating Hx
				Chxh[i][j][k]  =  (2*mu_r_x[i][j][k]*mu_0 - dt*sigma_m_x[i][j][k])/(2*mu_r_x[i][j][k]*mu_0 + dt*sigma_m_x[i][j][k]);
				Chxez[i][j][k] = -(2*dt/dy)/(2*mu_r_x[i][j][k]*mu_0 + dt*sigma_m_x[i][j][k]);
				Chxey[i][j][k] =  (2*dt/dz)/(2*mu_r_x[i][j][k]*mu_0 + dt*sigma_m_x[i][j][k]);
			}
		}
	}

	for(i=0;i<nx;i++){
		for(j=0;j<nyp1;j++){
			for(k=0;k<nz;k++){
				//Coeffiecients updating Hy
				Chyh[i][j][k]  =  (2*mu_r_y[i][j][k]*mu_0 - dt*sigma_m_y[i][j][k])/(2*mu_r_y[i][j][k]*mu_0 + dt*sigma_m_y[i][j][k]);
				Chyex[i][j][k] = -(2*dt/dz)/(2*mu_r_y[i][j][k]*mu_0 + dt*sigma_m_y[i][j][k]);
				Chyez[i][j][k] =  (2*dt/dx)/(2*mu_r_y[i][j][k]*mu_0 + dt*sigma_m_y[i][j][k]);
			}
		}
	}

	for(i=0;i<nx;i++){
		for(j=0;j<ny;j++){
			for(k=0;k<nzp1;k++){
				//Coeffiecients updating Hz
				Chzh[i][j][k]  =  (2*mu_r_z[i][j][k]*mu_0 - dt*sigma_m_z[i][j][k])/(2*mu_r_z[i][j][k]*mu_0 + dt*sigma_m_z[i][j][k]);
				Chzey[i][j][k] = -(2*dt/dx)/(2*mu_r_z[i][j][k]*mu_0 + dt*sigma_m_z[i][j][k]);
				Chzex[i][j][k] =  (2*dt/dy)/(2*mu_r_z[i][j][k]*mu_0 + dt*sigma_m_z[i][j][k]);
			}
		}
	}

	printf("initializing voltage source updating coefficients...\n");
	
	double a_term; 
	double R;

	for (ind = 0; ind<number_of_voltage_sources;ind++){
		is = voltage_sources[ind].is;
		js = voltage_sources[ind].js;
		ks = voltage_sources[ind].ks;
		ie = voltage_sources[ind].ie;
		je = voltage_sources[ind].je;
		ke = voltage_sources[ind].ke;
		R = voltage_sources[ind].resistance_per_component;

		if (R == 0) {
			R = 1e-20;
		}
    
		if (voltage_sources[ind].direction==XP || voltage_sources[ind].direction==XN){
			a_term = (dt*dx)/(R*dy*dz);
			voltage_sources[ind].Cexs = build_3D_double_array(ie-is, je+1-js, ke+1-ks, 0.0);
			for(i = is;i<=ie-1;i++){
				for(j = js;j<=je;j++){
					for(k = ks;k<=ke;k++){
						Cexe[i][j][k] = (2*eps_0*eps_r_x[i][j][k] - dt*sigma_e_x[i][j][k] - a_term)/(2*eps_0*eps_r_x[i][j][k] + dt*sigma_e_x[i][j][k] + a_term);
						Cexhz[i][j][k] = (2*dt/dy)/(2*eps_r_x[i][j][k]*eps_0 + dt*sigma_e_x[i][j][k] + a_term);
						Cexhy[i][j][k]= -(2*dt/dz)/(2*eps_r_x[i][j][k]*eps_0 + dt*sigma_e_x[i][j][k] + a_term);
						voltage_sources[ind].Cexs[i-is][j-js][k-ks] = -(2*dt/(R*dy*dz))/(2*eps_r_x[i][j][k]*eps_0 + dt*sigma_e_x[i][j][k] + a_term);
					}
				}
			}
		}

		if (voltage_sources[ind].direction==YP || voltage_sources[ind].direction==YN){
			a_term = (dt*dy)/(R*dz*dx);
			voltage_sources[ind].Ceys = build_3D_double_array(ie+1-is, je-js, ke+1-ks, 0.0);
			for(i = is;i<=ie;i++){
				for(j = js;j<=je-1;j++){
					for(k = ks;k<=ke;k++){
						Ceye[i][j][k] = (2*eps_0*eps_r_y[i][j][k] - dt*sigma_e_y[i][j][k] - a_term)/(2*eps_0*eps_r_y[i][j][k] + dt*sigma_e_y[i][j][k] + a_term);
						Ceyhx[i][j][k] = (2*dt/dz)/(2*eps_r_y[i][j][k]*eps_0 + dt*sigma_e_y[i][j][k] + a_term);
						Ceyhz[i][j][k] = -(2*dt/dx)/(2*eps_r_y[i][j][k]*eps_0 + dt*sigma_e_y[i][j][k] + a_term);
						voltage_sources[ind].Ceys[i-is][j-js][k-ks] = -(2*dt/(R*dz*dx))/(2*eps_r_y[i][j][k]*eps_0 + dt*sigma_e_y[i][j][k] + a_term);
					}
				}
			}
		}

		if (voltage_sources[ind].direction==ZP || voltage_sources[ind].direction==ZN){
			a_term = (dt*dz)/(R*dx*dy);
			voltage_sources[ind].Cezs = build_3D_double_array(ie+1-is, je+1-js, ke-ks, 0.0);
			for(i = is;i<=ie;i++){
				for(j = js;j<=je;j++){
					for(k = ks;k<=ke-1;k++){
						Ceze[i][j][k] = (2*eps_0*eps_r_z[i][j][k] - dt*sigma_e_z[i][j][k] - a_term)/(2*eps_0*eps_r_z[i][j][k] + dt*sigma_e_z[i][j][k] + a_term);
						Cezhy[i][j][k]= (2*dt/dx)/ (2*eps_r_z[i][j][k]*eps_0 + dt*sigma_e_z[i][j][k] + a_term);
						Cezhx[i][j][k]= -(2*dt/dy)/ (2*eps_r_z[i][j][k]*eps_0 + dt*sigma_e_z[i][j][k] + a_term);      
						voltage_sources[ind].Cezs[i-is][j-js][k-ks] = -(2*dt/(R*dx*dy))/(2*eps_r_z[i][j][k]*eps_0 + dt*sigma_e_z[i][j][k] + a_term);
					}
				}
			}
        }
	}


	printf("initialize boundary parameters...\n");

	int is_cpml_xn = FALSE;
	int is_cpml_xp = FALSE;
	int is_cpml_yn = FALSE;
	int is_cpml_yp = FALSE;
	int is_cpml_zn = FALSE;
	int is_cpml_zp = FALSE;
	int is_any_side_cpml = FALSE;

	int n_cpml_xn = 0, n_cpml_xp = 0, n_cpml_yn = 0, n_cpml_yp = 0, n_cpml_zn = 0, n_cpml_zp = 0;

	if(boundary.type_xn == CPML)
	{
		is_cpml_xn = TRUE;
		n_cpml_xn = abs(boundary.cpml_number_of_cells_xn);
	}
	if(boundary.type_xp == CPML)
	{
		is_cpml_xp = TRUE;
		n_cpml_xp = abs(boundary.cpml_number_of_cells_xp);
	}
	if(boundary.type_yn == CPML)
	{
		is_cpml_yn = TRUE;
		n_cpml_yn = abs(boundary.cpml_number_of_cells_yn);
	}
	if(boundary.type_yp == CPML)
	{
		is_cpml_yp = TRUE;
		n_cpml_yp = abs(boundary.cpml_number_of_cells_yp);
	}
	if(boundary.type_zn == CPML)
	{
		is_cpml_zn = TRUE;
		n_cpml_zn = abs(boundary.cpml_number_of_cells_zn);
	}
	if(boundary.type_zp == CPML)
	{
		is_cpml_zp = TRUE;
		n_cpml_zp = abs(boundary.cpml_number_of_cells_zp);
	}
	if (is_cpml_xn == TRUE || is_cpml_xp == TRUE || is_cpml_yn == TRUE || is_cpml_yp == TRUE || is_cpml_zn == TRUE || is_cpml_zp == TRUE )
	{
		is_any_side_cpml = TRUE;
		printf("Call CPML initialization routine if any side is CPML...\n");	
	}
//////////
	double *sigma_pex_xn;
	double *sigma_pmx_xn;
	double *kappa_ex_xn;
	double *kappa_mx_xn;
	double *alpha_ex_xn;
	double *alpha_mx_xn;
	double *cpml_b_ex_xn;
	double *cpml_a_ex_xn;
	double *cpml_b_mx_xn;
	double *cpml_a_mx_xn;

	double ***Psi_eyx_xn; 
	double ***Psi_ezx_xn; 
	double ***Psi_hyx_xn; 
	double ***Psi_hzx_xn; 

	double ***CPsi_eyx_xn;
	double ***CPsi_ezx_xn;
	double ***CPsi_hyx_xn;
	double ***CPsi_hzx_xn;

	double *sigma_pex_xp;
	double *sigma_pmx_xp;
	double *kappa_ex_xp;
	double *kappa_mx_xp;
	double *alpha_ex_xp;
	double *alpha_mx_xp;
	double *cpml_b_ex_xp;
	double *cpml_a_ex_xp;
	double *cpml_b_mx_xp;
	double *cpml_a_mx_xp;

	double ***Psi_eyx_xp; 
	double ***Psi_ezx_xp; 
	double ***Psi_hyx_xp; 
	double ***Psi_hzx_xp; 

	double ***CPsi_eyx_xp;
	double ***CPsi_ezx_xp;
	double ***CPsi_hyx_xp;
	double ***CPsi_hzx_xp;

	double *sigma_pey_yn;
	double *sigma_pmy_yn;
	double *kappa_ey_yn;
	double *kappa_my_yn;
	double *alpha_ey_yn;
	double *alpha_my_yn;
	double *cpml_b_ey_yn;
	double *cpml_a_ey_yn;
	double *cpml_b_my_yn;
	double *cpml_a_my_yn;

	double ***Psi_ezy_yn; 
	double ***Psi_exy_yn; 
	double ***Psi_hzy_yn; 
	double ***Psi_hxy_yn; 

	double ***CPsi_ezy_yn;
	double ***CPsi_exy_yn;
	double ***CPsi_hzy_yn;
	double ***CPsi_hxy_yn;

	double *sigma_pey_yp;
	double *sigma_pmy_yp;
	double *kappa_ey_yp;
	double *kappa_my_yp;
	double *alpha_ey_yp;
	double *alpha_my_yp;
	double *cpml_b_ey_yp;
	double *cpml_a_ey_yp;
	double *cpml_b_my_yp;
	double *cpml_a_my_yp;

	double ***Psi_exy_yp; 
	double ***Psi_ezy_yp; 
	double ***Psi_hxy_yp; 
	double ***Psi_hzy_yp; 

	double ***CPsi_ezy_yp;
	double ***CPsi_exy_yp;
	double ***CPsi_hzy_yp;
	double ***CPsi_hxy_yp;

	double *sigma_pez_zn;
	double *sigma_pmz_zn;
	double *kappa_ez_zn;
	double *kappa_mz_zn;
	double *alpha_ez_zn;
	double *alpha_mz_zn;
	double *cpml_b_ez_zn;
	double *cpml_a_ez_zn;
	double *cpml_b_mz_zn;
	double *cpml_a_mz_zn;

	double ***Psi_eyz_zn; 
	double ***Psi_exz_zn; 
	double ***Psi_hyz_zn; 
	double ***Psi_hxz_zn; 

	double ***CPsi_eyz_zn;
	double ***CPsi_exz_zn;
	double ***CPsi_hyz_zn;
	double ***CPsi_hxz_zn;

	double *sigma_pez_zp;
	double *sigma_pmz_zp;
	double *kappa_ez_zp;
	double *kappa_mz_zp;
	double *alpha_ez_zp;
	double *alpha_mz_zp;
	double *cpml_b_ez_zp;
	double *cpml_a_ez_zp;
	double *cpml_b_mz_zp;
	double *cpml_a_mz_zp;

	double ***Psi_exz_zp; 
	double ***Psi_eyz_zp; 
	double ***Psi_hxz_zp; 
	double ***Psi_hyz_zp; 

	double ***CPsi_eyz_zp;
	double ***CPsi_exz_zp;
	double ***CPsi_hyz_zp;
	double ***CPsi_hxz_zp;


//////////
	if(is_any_side_cpml==TRUE)
	{
		//Initialize CPML boundary condition    

		int p_order = boundary.cpml_order; //order of the polynomial distribution
		double sigma_ratio = boundary.cpml_sigma_factor;
		double kappa_max = boundary.cpml_kappa_max;
		double alpha_min = boundary.cpml_alpha_min;
		double alpha_max = boundary.cpml_alpha_max;
		double sigma_max;

		int ncells;
		double *rho_e_xn;
		double *rho_m_xn;

		double *rho_e_xp;
		double *rho_m_xp;

		double *rho_e_yn;
		double *rho_m_yn;

		double *rho_e_yp;
		double *rho_m_yp;

		double *rho_e_zn;
		double *rho_m_zn;

		double *rho_e_zp;
		double *rho_m_zp;

		int ncells_ind;
		//Initialize cpml for xn region
		if (is_cpml_xn == TRUE)
		{
			//define one-dimensional temporary cpml parameter arrays 
			sigma_max = sigma_ratio*(p_order+1)/(150*M_PI*dx);
			ncells = n_cpml_xn;
	
			rho_e_xn = build_1D_double_array(ncells,0.0);
			rho_m_xn = build_1D_double_array(ncells,0.0);
			sigma_pex_xn = build_1D_double_array(ncells,0.0);
			sigma_pmx_xn = build_1D_double_array(ncells,0.0);
			kappa_ex_xn = build_1D_double_array(ncells,0.0);
			kappa_mx_xn = build_1D_double_array(ncells,0.0);
			alpha_ex_xn = build_1D_double_array(ncells,0.0);
			alpha_mx_xn = build_1D_double_array(ncells,0.0);
			cpml_b_ex_xn = build_1D_double_array(ncells,0.0);
			cpml_a_ex_xn = build_1D_double_array(ncells,0.0);
			cpml_b_mx_xn = build_1D_double_array(ncells,0.0);
			cpml_a_mx_xn = build_1D_double_array(ncells,0.0);

			for(ncells_ind = 0; ncells_ind < ncells; ncells_ind++){

				//sigma_pex_xn[ncells_ind] = 0.0;
				sigma_pmx_xn[ncells_ind] = 0.0;

				kappa_ex_xn[ncells_ind] = 0.0;
				kappa_mx_xn[ncells_ind] = 0.0;

				alpha_ex_xn[ncells_ind] = 0.0;
				alpha_mx_xn[ncells_ind] = 0.0;

				cpml_b_ex_xn[ncells_ind] = 0.0;
				cpml_a_ex_xn[ncells_ind] = 0.0;
				cpml_b_mx_xn[ncells_ind] = 0.0;
				cpml_a_mx_xn[ncells_ind] = 0.0;

				rho_e_xn[ncells_ind] = (((double)ncells - (double)ncells_ind) - 0.75)/((double)ncells);
				rho_m_xn[ncells_ind] = (((double)ncells - (double)ncells_ind) - 0.25)/((double)ncells);

				sigma_pex_xn[ncells_ind] = 1;
				sigma_pmx_xn[ncells_ind] = 1;
				
				int p_order_ind;
				for(p_order_ind = 1; p_order_ind <= p_order; p_order_ind++){
					sigma_pex_xn[ncells_ind] = sigma_pex_xn[ncells_ind]*rho_e_xn[ncells_ind];
					sigma_pmx_xn[ncells_ind] = sigma_pmx_xn[ncells_ind]*rho_m_xn[ncells_ind];
				}

				kappa_ex_xn[ncells_ind] = sigma_pex_xn[ncells_ind];
				kappa_mx_xn[ncells_ind] = sigma_pmx_xn[ncells_ind];

				sigma_pex_xn[ncells_ind] = sigma_max*sigma_pex_xn[ncells_ind];
				sigma_pmx_xn[ncells_ind] = sigma_max*sigma_pmx_xn[ncells_ind];

				kappa_ex_xn[ncells_ind] = 1 + (kappa_max - 1)*kappa_ex_xn[ncells_ind];
				kappa_mx_xn[ncells_ind] = 1 + (kappa_max - 1)*kappa_mx_xn[ncells_ind];
			
				sigma_pmx_xn[ncells_ind] = (mu_0/eps_0)*sigma_pmx_xn[ncells_ind];
				alpha_ex_xn[ncells_ind] = alpha_min + (alpha_max - alpha_min)*(1 - rho_e_xn[ncells_ind]);
				alpha_mx_xn[ncells_ind] = alpha_min + (alpha_max - alpha_min)*(1 - rho_m_xn[ncells_ind]);
				alpha_mx_xn[ncells_ind] = (mu_0/eps_0)*alpha_mx_xn[ncells_ind];

				//define one-dimensional cpml parameter arrays 
				cpml_b_ex_xn[ncells_ind] = exp((-dt/eps_0)*((sigma_pex_xn[ncells_ind]/kappa_ex_xn[ncells_ind]) + alpha_ex_xn[ncells_ind])); 
				cpml_a_ex_xn[ncells_ind] = (1/dx)*(cpml_b_ex_xn[ncells_ind] - 1.0)* sigma_pex_xn[ncells_ind]/(kappa_ex_xn[ncells_ind]*(sigma_pex_xn[ncells_ind] + kappa_ex_xn[ncells_ind]*alpha_ex_xn[ncells_ind]));
				cpml_b_mx_xn[ncells_ind] = exp((-dt/mu_0)*((sigma_pmx_xn[ncells_ind]/kappa_mx_xn[ncells_ind]) + alpha_mx_xn[ncells_ind])); 
				cpml_a_mx_xn[ncells_ind] = (1/dx)*(cpml_b_mx_xn[ncells_ind] - 1.0)*sigma_pmx_xn[ncells_ind]/(kappa_mx_xn[ncells_ind]*(sigma_pmx_xn[ncells_ind] + kappa_mx_xn[ncells_ind]*alpha_mx_xn[ncells_ind]));
			}
			//Create and initialize 2D cpml convolution parameters 
			Psi_eyx_xn = build_3D_double_array(ncells,ny,nzp1,0.0); 
			Psi_ezx_xn = build_3D_double_array(ncells,nyp1,nz,0.0); 
			Psi_hyx_xn = build_3D_double_array(ncells,nyp1,nz,0.0); 
			Psi_hzx_xn = build_3D_double_array(ncells,ny,nzp1,0.0);

			//Create and initialize 2D cpml convolution coefficients

			CPsi_eyx_xn = build_3D_double_array(ncells,ny,nzp1,0.0); 
			CPsi_hzx_xn = build_3D_double_array(ncells,ny,nzp1,0.0); 
			CPsi_ezx_xn = build_3D_double_array(ncells,nyp1,nz,0.0); 
			CPsi_hyx_xn = build_3D_double_array(ncells,nyp1,nz,0.0);
 
			for(i=0;i<ncells;i++){
				for(j=0;j<nyp1;j++){
					for(k=0;k<nzp1;k++){ 
						if(j<ny){	
							CPsi_eyx_xn[i][j][k] = Ceyhz[i+1][j][k]*dx;
							Ceyhz[i+1][j][k] = Ceyhz[i+1][j][k]/kappa_ex_xn[i];
							CPsi_hzx_xn[i][j][k] = Chzey[i][j][k]*dx;
							Chzey[i][j][k] = Chzey[i][j][k]/kappa_mx_xn[i];
						}
						if(k<nz){
							CPsi_ezx_xn[i][j][k] = Cezhy[i+1][j][k]*dx;
							Cezhy[i+1][j][k] = Cezhy[i+1][j][k]/kappa_ex_xn[i];
							CPsi_hyx_xn[i][j][k] = Chyez[i][j][k]*dx;
							Chyez[i][j][k] = Chyez[i][j][k]/kappa_mx_xn[i];
						}
					}
				}
			}
		}

		//Initialize cpml for xp region
		if (is_cpml_xp == TRUE)
		{
			//define one-dimensional temporary cpml parameter arrays 
			sigma_max = sigma_ratio*(p_order + 1)/(150*M_PI*dx);
			ncells = n_cpml_xp;

			rho_e_xp = build_1D_double_array(ncells,0.0);
			rho_m_xp = build_1D_double_array(ncells,0.0);
			sigma_pex_xp = build_1D_double_array(ncells,0.0);
			sigma_pmx_xp = build_1D_double_array(ncells,0.0);
			kappa_ex_xp = build_1D_double_array(ncells,0.0);
			kappa_mx_xp = build_1D_double_array(ncells,0.0);
			alpha_ex_xp = build_1D_double_array(ncells,0.0);
			alpha_mx_xp = build_1D_double_array(ncells,0.0);
			cpml_b_ex_xp = build_1D_double_array(ncells,0.0);
			cpml_a_ex_xp = build_1D_double_array(ncells,0.0);
			cpml_b_mx_xp = build_1D_double_array(ncells,0.0);
			cpml_a_mx_xp = build_1D_double_array(ncells,0.0);
			
			for(ncells_ind = 0;ncells_ind<ncells;ncells_ind++){
				rho_e_xp[ncells_ind] = (((double)ncells_ind + 1) - 0.75)/((double)ncells);
				rho_m_xp[ncells_ind] = (((double)ncells_ind + 1) - 0.25)/((double)ncells);

				sigma_pex_xp[ncells_ind] = 1;
				sigma_pmx_xp[ncells_ind] = 1;
				
				int p_order_ind;
				for(p_order_ind = 1 ;p_order_ind <= p_order; p_order_ind++){
					sigma_pex_xp[ncells_ind] = sigma_pex_xp[ncells_ind]*rho_e_xp[ncells_ind];
					sigma_pmx_xp[ncells_ind] = sigma_pmx_xp[ncells_ind]*rho_m_xp[ncells_ind];
				}
				kappa_ex_xp[ncells_ind] = sigma_pex_xp[ncells_ind];
				kappa_mx_xp[ncells_ind] = sigma_pmx_xp[ncells_ind];

				sigma_pex_xp[ncells_ind] = sigma_max*sigma_pex_xp[ncells_ind];
				sigma_pmx_xp[ncells_ind] = sigma_max*sigma_pmx_xp[ncells_ind];

				kappa_ex_xp[ncells_ind] = 1 + (kappa_max - 1)*kappa_ex_xp[ncells_ind];
				kappa_mx_xp[ncells_ind] = 1 + (kappa_max - 1)*kappa_mx_xp[ncells_ind];
			
				sigma_pmx_xp[ncells_ind] = (mu_0/eps_0)*sigma_pmx_xp[ncells_ind];
				alpha_ex_xp[ncells_ind] = alpha_min + (alpha_max - alpha_min)*(1 - rho_e_xp[ncells_ind]);
				alpha_mx_xp[ncells_ind] = alpha_min + (alpha_max - alpha_min)*(1 - rho_m_xp[ncells_ind]);
				alpha_mx_xp[ncells_ind] = (mu_0/eps_0)*alpha_mx_xp[ncells_ind];

				//define one-dimensional cpml parameter arrays 
				cpml_b_ex_xp[ncells_ind] = exp((-dt/eps_0)*((sigma_pex_xp[ncells_ind]/kappa_ex_xp[ncells_ind]) + alpha_ex_xp[ncells_ind]));
				cpml_a_ex_xp[ncells_ind] = (1/dx)*(cpml_b_ex_xp[ncells_ind] - 1.0)* sigma_pex_xp[ncells_ind]/(kappa_ex_xp[ncells_ind]*(sigma_pex_xp[ncells_ind] + kappa_ex_xp[ncells_ind]*alpha_ex_xp[ncells_ind]));
				cpml_b_mx_xp[ncells_ind] = exp((-dt/mu_0)*((sigma_pmx_xp[ncells_ind]/kappa_mx_xp[ncells_ind]) + alpha_mx_xp[ncells_ind]));
				cpml_a_mx_xp[ncells_ind] = (1/dx)*(cpml_b_mx_xp[ncells_ind] - 1.0)*sigma_pmx_xp[ncells_ind]/(kappa_mx_xp[ncells_ind]*(sigma_pmx_xp[ncells_ind] + kappa_mx_xp[ncells_ind]*alpha_mx_xp[ncells_ind]));
			}

			//Create and initialize 2D cpml convolution parameters 
			Psi_eyx_xp = build_3D_double_array(ncells,ny,nzp1,0.0); 
			Psi_ezx_xp = build_3D_double_array(ncells,nyp1,nz,0.0); 
			Psi_hyx_xp = build_3D_double_array(ncells,nyp1,nz,0.0); 
			Psi_hzx_xp = build_3D_double_array(ncells,ny,nzp1,0.0); 

			//Create and initialize 2D cpml convolution coefficients

			CPsi_eyx_xp = build_3D_double_array(ncells,ny,nzp1,0.0); 
			CPsi_hzx_xp = build_3D_double_array(ncells,ny,nzp1,0.0); 
			CPsi_ezx_xp = build_3D_double_array(ncells,nyp1,nz,0.0); 
			CPsi_hyx_xp = build_3D_double_array(ncells,nyp1,nz,0.0); 

			for(i=0;i<ncells;i++){
				for(j=0;j<nyp1;j++){
					for(k=0;k<nzp1;k++){ 
						if(j<ny){	
							CPsi_eyx_xp[i][j][k] = Ceyhz[nx - ncells + i][j][k]*dx;
							Ceyhz[nx - ncells + i][j][k] = Ceyhz[nx - ncells + i][j][k]/kappa_ex_xp[i];
							CPsi_hzx_xp[i][j][k] = Chzey[nx - ncells + i][j][k]*dx;
							Chzey[nx - ncells + i][j][k] = Chzey[nx - ncells + i][j][k]/kappa_mx_xp[i];
						}
						if(k<nz){
							CPsi_ezx_xp[i][j][k] = Cezhy[nx - ncells + i][j][k]*dx;
							Cezhy[nx - ncells + i][j][k] = Cezhy[nx - ncells + i][j][k]/kappa_ex_xp[i];
							CPsi_hyx_xp[i][j][k] = Chyez[nx - ncells + i][j][k]*dx;
							Chyez[nx - ncells + i][j][k] = Chyez[nx - ncells + i][j][k]/kappa_mx_xp[i];
						}
					}
				}
			}
  		}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//Initialize cpml for yn region
		if (is_cpml_yn == TRUE)
		{
			//define one-dimensional temporary cpml parameter arrays 
			sigma_max = sigma_ratio*(p_order+1)/(150*M_PI*dy);
			ncells = n_cpml_yn;
	
			rho_e_yn = build_1D_double_array(ncells,0.0);
			rho_m_yn = build_1D_double_array(ncells,0.0);
			sigma_pey_yn = build_1D_double_array(ncells,0.0);
			sigma_pmy_yn = build_1D_double_array(ncells,0.0);
			kappa_ey_yn = build_1D_double_array(ncells,0.0);
			kappa_my_yn = build_1D_double_array(ncells,0.0);
			alpha_ey_yn = build_1D_double_array(ncells,0.0);
			alpha_my_yn = build_1D_double_array(ncells,0.0);
			cpml_b_ey_yn = build_1D_double_array(ncells,0.0);
			cpml_a_ey_yn = build_1D_double_array(ncells,0.0);
			cpml_b_my_yn = build_1D_double_array(ncells,0.0);
			cpml_a_my_yn = build_1D_double_array(ncells,0.0);

			for(ncells_ind = 0;ncells_ind<ncells;ncells_ind++){
				rho_e_yn[ncells_ind] = (((double)ncells - (double)ncells_ind) - 0.75)/((double)ncells);
				rho_m_yn[ncells_ind] = (((double)ncells - (double)ncells_ind) - 0.25)/((double)ncells);

				sigma_pey_yn[ncells_ind] = 1;
				sigma_pmy_yn[ncells_ind] = 1;
				
				int p_order_ind;
				for(p_order_ind = 1; p_order_ind<=p_order; p_order_ind++){
					sigma_pey_yn[ncells_ind] = sigma_pey_yn[ncells_ind]*rho_e_yn[ncells_ind];
					sigma_pmy_yn[ncells_ind] = sigma_pmy_yn[ncells_ind]*rho_m_yn[ncells_ind];
				}

				kappa_ey_yn[ncells_ind] = sigma_pey_yn[ncells_ind];
				kappa_my_yn[ncells_ind] = sigma_pmy_yn[ncells_ind];

				sigma_pey_yn[ncells_ind] = sigma_max*sigma_pey_yn[ncells_ind];
				sigma_pmy_yn[ncells_ind] = sigma_max*sigma_pmy_yn[ncells_ind];

				kappa_ey_yn[ncells_ind] = 1 + (kappa_max - 1)*kappa_ey_yn[ncells_ind];
				kappa_my_yn[ncells_ind] = 1 + (kappa_max - 1)*kappa_my_yn[ncells_ind];
			
				sigma_pmy_yn[ncells_ind] = (mu_0 / eps_0)*sigma_pmy_yn[ncells_ind];
				alpha_ey_yn[ncells_ind] = alpha_min + (alpha_max - alpha_min)*(1 - rho_e_yn[ncells_ind]);
				alpha_my_yn[ncells_ind] = alpha_min + (alpha_max - alpha_min)*(1 - rho_m_yn[ncells_ind]);
				alpha_my_yn[ncells_ind] = (mu_0/eps_0)*alpha_my_yn[ncells_ind];

				//define one-dimensional cpml parameter arrays 
				cpml_b_ey_yn[ncells_ind] = exp((-dt/eps_0)*((sigma_pey_yn[ncells_ind]/kappa_ey_yn[ncells_ind]) + alpha_ey_yn[ncells_ind])); 
				cpml_a_ey_yn[ncells_ind] = (1/dx)*(cpml_b_ey_yn[ncells_ind] - 1.0)* sigma_pey_yn[ncells_ind]/(kappa_ey_yn[ncells_ind]*(sigma_pey_yn[ncells_ind] + kappa_ey_yn[ncells_ind]*alpha_ey_yn[ncells_ind]));
				cpml_b_my_yn[ncells_ind] = exp((-dt/mu_0)*((sigma_pmy_yn[ncells_ind]/kappa_my_yn[ncells_ind]) + alpha_my_yn[ncells_ind])); 
				cpml_a_my_yn[ncells_ind] = (1/dx)*(cpml_b_my_yn[ncells_ind] - 1.0)*sigma_pmy_yn[ncells_ind]/(kappa_my_yn[ncells_ind]*(sigma_pmy_yn[ncells_ind] + kappa_my_yn[ncells_ind]*alpha_my_yn[ncells_ind]));
			}
			//Create and initialize 2D cpml convolution parameters 
			Psi_ezy_yn = build_3D_double_array(nxp1,ncells,nz,0.0); 
			Psi_exy_yn = build_3D_double_array(nx,ncells,nzp1,0.0); 
			Psi_hzy_yn = build_3D_double_array(nx,ncells,nzp1,0.0); 
			Psi_hxy_yn = build_3D_double_array(nxp1,ncells,nz,0.0); 

			//Create and initialize 2D cpml convolution coefficients

			CPsi_ezy_yn = build_3D_double_array(nxp1,ncells,nz,0.0); 
			CPsi_hxy_yn = build_3D_double_array(nxp1,ncells,nz,0.0); 
			CPsi_exy_yn = build_3D_double_array(nx,ncells,nzp1,0.0); 
			CPsi_hzy_yn = build_3D_double_array(nx,ncells,nzp1,0.0);
 
			for(i=0;i<nxp1;i++){
				for(j=0;j<ncells;j++){
					for(k=0;k<nzp1;k++){ 
						if(i<nx){
							CPsi_exy_yn[i][j][k] = Cexhz[i][j+1][k]*dx;
							Cexhz[i][j+1][k] = Cexhz[i][j+1][k]/kappa_ey_yn[j];
							CPsi_hzy_yn[i][j][k] = Chzex[i][j][k]*dx;
							Chzex[i][j][k] = Chzex[i][j][k]/kappa_my_yn[j];

						}
						if(k<nz){
							CPsi_hxy_yn[i][j][k] = Chxez[i][j][k]*dx;
							Chxez[i][j][k] = Chxez[i][j][k]/kappa_my_yn[j];
							CPsi_ezy_yn[i][j][k] = Cezhx[i][j+1][k]*dx;
							Cezhx[i][j+1][k] = Cezhx[i][j+1][k]/kappa_ey_yn[j];

						}
					}
				}
			}
		}

		//Initialize cpml for yp region
		if (is_cpml_yp == TRUE)
		{
			//define one-dimensional temporary cpml parameter arrays 
			sigma_max = sigma_ratio*(p_order + 1)/(150*M_PI*dy);
			ncells = n_cpml_yp;

			rho_e_yp = build_1D_double_array(ncells,0.0);
			rho_m_yp = build_1D_double_array(ncells,0.0);
			sigma_pey_yp = build_1D_double_array(ncells,0.0);
			sigma_pmy_yp = build_1D_double_array(ncells,0.0);
			kappa_ey_yp = build_1D_double_array(ncells,0.0);
			kappa_my_yp = build_1D_double_array(ncells,0.0);
			alpha_ey_yp = build_1D_double_array(ncells,0.0);
			alpha_my_yp = build_1D_double_array(ncells,0.0);
			cpml_b_ey_yp = build_1D_double_array(ncells,0.0);
			cpml_a_ey_yp = build_1D_double_array(ncells,0.0);
			cpml_b_my_yp = build_1D_double_array(ncells,0.0);
			cpml_a_my_yp = build_1D_double_array(ncells,0.0);
			
			for(ncells_ind = 0;ncells_ind<ncells;ncells_ind++){
				rho_e_yp[ncells_ind] = (((double)ncells_ind + 1) - 0.75)/((double)ncells);
				rho_m_yp[ncells_ind] = (((double)ncells_ind + 1) - 0.25)/((double)ncells);

				sigma_pey_yp[ncells_ind] = 1;
				sigma_pmy_yp[ncells_ind] = 1;
				
				int p_order_ind;
				for(p_order_ind = 1; p_order_ind <= p_order; p_order_ind++){
					sigma_pey_yp[ncells_ind] = sigma_pey_yp[ncells_ind]*rho_e_yp[ncells_ind];
					sigma_pmy_yp[ncells_ind] = sigma_pmy_yp[ncells_ind]*rho_m_yp[ncells_ind];
				}
				kappa_ey_yp[ncells_ind] = sigma_pey_yp[ncells_ind];
				kappa_my_yp[ncells_ind] = sigma_pmy_yp[ncells_ind];

				sigma_pey_yp[ncells_ind] = sigma_max*sigma_pey_yp[ncells_ind];
				sigma_pmy_yp[ncells_ind] = sigma_max*sigma_pmy_yp[ncells_ind];

				kappa_ey_yp[ncells_ind] = 1 + (kappa_max - 1)*kappa_ey_yp[ncells_ind];
				kappa_my_yp[ncells_ind] = 1 + (kappa_max - 1)*kappa_my_yp[ncells_ind];
			
				sigma_pmy_yp[ncells_ind] = (mu_0/eps_0)*sigma_pmy_yp[ncells_ind];
				alpha_ey_yp[ncells_ind] = alpha_min + (alpha_max - alpha_min)*(1 - rho_e_yp[ncells_ind]);
				alpha_my_yp[ncells_ind] = alpha_min + (alpha_max - alpha_min)*(1 - rho_m_yp[ncells_ind]);
				alpha_my_yp[ncells_ind] = (mu_0/eps_0)*alpha_my_yp[ncells_ind];

				//define one-dimensional cpml parameter arrays 
				cpml_b_ey_yp[ncells_ind] = exp((-dt/eps_0)*((sigma_pey_yp[ncells_ind]/kappa_ey_yp[ncells_ind]) + alpha_ey_yp[ncells_ind]));
				cpml_a_ey_yp[ncells_ind] = (1/dx)*(cpml_b_ey_yp[ncells_ind] - 1.0)* sigma_pey_yp[ncells_ind]/(kappa_ey_yp[ncells_ind]*(sigma_pey_yp[ncells_ind] + kappa_ey_yp[ncells_ind]*alpha_ey_yp[ncells_ind]));
				cpml_b_my_yp[ncells_ind] = exp((-dt/mu_0)*((sigma_pmy_yp[ncells_ind]/kappa_my_yp[ncells_ind]) + alpha_my_yp[ncells_ind]));
				cpml_a_my_yp[ncells_ind] = (1/dx)*(cpml_b_my_yp[ncells_ind] - 1.0)*sigma_pmy_yp[ncells_ind]/(kappa_my_yp[ncells_ind]*(sigma_pmy_yp[ncells_ind] + kappa_my_yp[ncells_ind]*alpha_my_yp[ncells_ind]));
			}

			//Create and initialize 2D cpml convolution parameters 
			Psi_ezy_yp = build_3D_double_array(nxp1,ncells,nz,0.0); 
			Psi_exy_yp = build_3D_double_array(nx,ncells,nzp1,0.0); 
			Psi_hzy_yp = build_3D_double_array(nx,ncells,nzp1,0.0); 
			Psi_hxy_yp = build_3D_double_array(nxp1,ncells,nz,0.0); 

			//Create and initialize 2D cpml convolution coefficients
			
			CPsi_exy_yp = build_3D_double_array(nx,ncells,nzp1,0.0); 
			CPsi_hzy_yp = build_3D_double_array(nx,ncells,nzp1,0.0); 
			CPsi_ezy_yp = build_3D_double_array(nxp1,ncells,nz,0.0); 
			CPsi_hxy_yp = build_3D_double_array(nxp1,ncells,nz,0.0); 

			for(i=0;i<nxp1;i++){
				for(j=0;j<ncells;j++){
					for(k=0;k<nzp1;k++){ 
						if(i<nx){	
							CPsi_exy_yp[i][j][k] = Cexhz[i][ny - ncells + j][k]*dy;
							Cexhz[i][ny - ncells + j][k] = Cexhz[i][ny - ncells + j][k]/kappa_ey_yp[j];
							CPsi_hzy_yp[i][j][k] = Chzex[i][ny - ncells + j][k]*dy;
							Chzex[i][ny - ncells + j][k] = Chzex[i][ny - ncells + j][k]/kappa_my_yp[j];
						}
						if(k<nz){
							CPsi_ezy_yp[i][j][k] = Cezhx[i][ny - ncells + j][k]*dy;
							Cezhx[i][ny - ncells + j][k] = Cezhx[i][ny - ncells + j][k]/kappa_ey_yp[j];
							CPsi_hxy_yp[i][j][k] = Chxez[i][ny - ncells + j][k]*dy;
							Chxez[i][ny - ncells + j][k] = Chxez[i][ny - ncells + j][k]/kappa_my_yp[j];
						}
					}
				}
			}
		}

		//Initialize cpml for zn region
		if (is_cpml_zn == TRUE)
		{
			//define one-dimensional temporary cpml parameter arrays 
			sigma_max = sigma_ratio*(p_order+1)/(150*M_PI*dz);
			ncells = n_cpml_zn;
	
			rho_e_zn = build_1D_double_array(ncells,0.0);
			rho_m_zn = build_1D_double_array(ncells,0.0);
			sigma_pez_zn = build_1D_double_array(ncells,0.0);
			sigma_pmz_zn = build_1D_double_array(ncells,0.0);
			kappa_ez_zn = build_1D_double_array(ncells,0.0);
			kappa_mz_zn = build_1D_double_array(ncells,0.0);
			alpha_ez_zn = build_1D_double_array(ncells,0.0);
			alpha_mz_zn = build_1D_double_array(ncells,0.0);
			cpml_b_ez_zn = build_1D_double_array(ncells,0.0);
			cpml_a_ez_zn = build_1D_double_array(ncells,0.0);
			cpml_b_mz_zn = build_1D_double_array(ncells,0.0);
			cpml_a_mz_zn = build_1D_double_array(ncells,0.0);

			for(ncells_ind = 0;ncells_ind<ncells;ncells_ind++){
				rho_e_zn[ncells_ind] = (((double)ncells - (double)ncells_ind) - 0.75)/((double)ncells);
				rho_m_zn[ncells_ind] = (((double)ncells - (double)ncells_ind) - 0.25)/((double)ncells);

				sigma_pez_zn[ncells_ind] = 1;
				sigma_pmz_zn[ncells_ind] = 1;
				
				int p_order_ind;
				for(p_order_ind = 1; p_order_ind <= p_order; p_order_ind++){
					sigma_pez_zn[ncells_ind] = sigma_pez_zn[ncells_ind]*rho_e_zn[ncells_ind];
					sigma_pmz_zn[ncells_ind] = sigma_pmz_zn[ncells_ind]*rho_m_zn[ncells_ind];
				}

				kappa_ez_zn[ncells_ind] = sigma_pez_zn[ncells_ind];
				kappa_mz_zn[ncells_ind] = sigma_pmz_zn[ncells_ind];

				sigma_pez_zn[ncells_ind] = sigma_max*sigma_pez_zn[ncells_ind];
				sigma_pmz_zn[ncells_ind] = sigma_max*sigma_pmz_zn[ncells_ind];

				kappa_ez_zn[ncells_ind] = 1 + (kappa_max - 1)*kappa_ez_zn[ncells_ind];
				kappa_mz_zn[ncells_ind] = 1 + (kappa_max - 1)*kappa_mz_zn[ncells_ind];
			
				sigma_pmz_zn[ncells_ind] = (mu_0 / eps_0)*sigma_pmz_zn[ncells_ind];
				alpha_ez_zn[ncells_ind] = alpha_min + (alpha_max - alpha_min)*(1 - rho_e_zn[ncells_ind]);
				alpha_mz_zn[ncells_ind] = alpha_min + (alpha_max - alpha_min)*(1 - rho_m_zn[ncells_ind]);
				alpha_mz_zn[ncells_ind] = (mu_0/eps_0)*alpha_mz_zn[ncells_ind];

				//define one-dimensional cpml parameter arrays 
				cpml_b_ez_zn[ncells_ind] = exp((-dt/eps_0)*((sigma_pez_zn[ncells_ind]/kappa_ez_zn[ncells_ind]) + alpha_ez_zn[ncells_ind])); 
				cpml_a_ez_zn[ncells_ind] = (1/dx)*(cpml_b_ez_zn[ncells_ind] - 1.0)* sigma_pez_zn[ncells_ind]/(kappa_ez_zn[ncells_ind]*(sigma_pez_zn[ncells_ind] + kappa_ez_zn[ncells_ind]*alpha_ez_zn[ncells_ind]));
				cpml_b_mz_zn[ncells_ind] = exp((-dt/mu_0)*((sigma_pmz_zn[ncells_ind]/kappa_mz_zn[ncells_ind]) + alpha_mz_zn[ncells_ind])); 
				cpml_a_mz_zn[ncells_ind] = (1/dx)*(cpml_b_mz_zn[ncells_ind] - 1.0)*sigma_pmz_zn[ncells_ind]/(kappa_mz_zn[ncells_ind]*(sigma_pmz_zn[ncells_ind] + kappa_mz_zn[ncells_ind]*alpha_mz_zn[ncells_ind]));
			}
			//Create and initialize 2D cpml convolution parameters 
			Psi_eyz_zn = build_3D_double_array(nxp1,ny,ncells,0.0); 
			Psi_exz_zn = build_3D_double_array(nx,nyp1,ncells,0.0); 
			Psi_hyz_zn = build_3D_double_array(nx,nyp1,ncells,0.0); 
			Psi_hxz_zn = build_3D_double_array(nxp1,ny,ncells,0.0); 

			//Create and initialize 2D cpml convolution coefficients
			
			CPsi_eyz_zn = build_3D_double_array(nxp1,ny,ncells,0.0); 
			CPsi_hxz_zn = build_3D_double_array(nxp1,ny,ncells,0.0); 
			CPsi_exz_zn = build_3D_double_array(nx,nyp1,ncells,0.0); 
			CPsi_hyz_zn = build_3D_double_array(nx,nyp1,ncells,0.0); 
			
			for(i=0;i<nxp1;i++){
				for(j=0;j<nyp1;j++){
					for(k=0;k<ncells;k++){ 
						if(j<ny){	
							CPsi_eyz_zn[i][j][k] = Ceyhx[i][j][k+1]*dz;
							Ceyhx[i][j][k+1] = Ceyhx[i][j][k+1]/kappa_ez_zn[k];
							CPsi_hxz_zn[i][j][k] = Chxey[i][j][k]*dz;
							Chxey[i][j][k] = Chxey[i][j][k]/kappa_mz_zn[k];
						}
						if(i<nx){
							CPsi_exz_zn[i][j][k] = Cexhy[i][j][k+1]*dz;
							Cexhy[i][j][k+1] = Cexhy[i][j][k+1]/kappa_ez_zn[k];
							CPsi_hyz_zn[i][j][k] = Chyex[i][j][k]*dz;
							Chyex[i][j][k] = Chyex[i][j][k]/kappa_mz_zn[k];
						}
					}
				}
			}
		}

		//Initialize cpml for zp region
		if (is_cpml_zp == TRUE)
		{
			//define one-dimensional temporary cpml parameter arrays 
			sigma_max = sigma_ratio*(p_order + 1)/(150*M_PI*dz);
			ncells = n_cpml_zp;


			rho_e_zp = build_1D_double_array(ncells,0.0);
			rho_m_zp = build_1D_double_array(ncells,0.0);
			sigma_pez_zp = build_1D_double_array(ncells,0.0);
			sigma_pmz_zp = build_1D_double_array(ncells,0.0);
			kappa_ez_zp = build_1D_double_array(ncells,0.0);
			kappa_mz_zp = build_1D_double_array(ncells,0.0);
			alpha_ez_zp = build_1D_double_array(ncells,0.0);
			alpha_mz_zp = build_1D_double_array(ncells,0.0);
			cpml_b_ez_zp = build_1D_double_array(ncells,0.0);
			cpml_a_ez_zp = build_1D_double_array(ncells,0.0);
			cpml_b_mz_zp = build_1D_double_array(ncells,0.0);
			cpml_a_mz_zp = build_1D_double_array(ncells,0.0);
			
			for(ncells_ind = 0;ncells_ind<ncells;ncells_ind++){
				rho_e_zp[ncells_ind] = (((double)ncells_ind + 1) - 0.75)/((double)ncells);
				rho_m_zp[ncells_ind] = (((double)ncells_ind + 1) - 0.25)/((double)ncells);

				sigma_pez_zp[ncells_ind] = 1;
				sigma_pmz_zp[ncells_ind] = 1;
				
				int p_order_ind;
				for(p_order_ind=1; p_order_ind <= p_order; p_order_ind++){
					sigma_pez_zp[ncells_ind] = sigma_pez_zp[ncells_ind]*rho_e_zp[ncells_ind];
					sigma_pmz_zp[ncells_ind] = sigma_pmz_zp[ncells_ind]*rho_m_zp[ncells_ind];
				}
				kappa_ez_zp[ncells_ind] = sigma_pez_zp[ncells_ind];
				kappa_mz_zp[ncells_ind] = sigma_pmz_zp[ncells_ind];

				sigma_pez_zp[ncells_ind] = sigma_max*sigma_pez_zp[ncells_ind];
				sigma_pmz_zp[ncells_ind] = sigma_max*sigma_pmz_zp[ncells_ind];

				kappa_ez_zp[ncells_ind] = 1 + (kappa_max - 1)*kappa_ez_zp[ncells_ind];
				kappa_mz_zp[ncells_ind] = 1 + (kappa_max - 1)*kappa_mz_zp[ncells_ind];
			
				sigma_pmz_zp[ncells_ind] = (mu_0 / eps_0)*sigma_pmz_zp[ncells_ind];
				alpha_ez_zp[ncells_ind] = alpha_min + (alpha_max - alpha_min)*(1 - rho_e_zp[ncells_ind]);
				alpha_mz_zp[ncells_ind] = alpha_min + (alpha_max - alpha_min)*(1 - rho_m_zp[ncells_ind]);
				alpha_mz_zp[ncells_ind] = (mu_0/eps_0)*alpha_mz_zp[ncells_ind];

				//define one-dimensional cpml parameter arrays 
				cpml_b_ez_zp[ncells_ind] = exp((-dt/eps_0)*((sigma_pez_zp[ncells_ind]/kappa_ez_zp[ncells_ind]) + alpha_ez_zp[ncells_ind]));
				cpml_a_ez_zp[ncells_ind] = (1/dz)*(cpml_b_ez_zp[ncells_ind] - 1.0)* sigma_pez_zp[ncells_ind]/(kappa_ez_zp[ncells_ind]*(sigma_pez_zp[ncells_ind] + kappa_ez_zp[ncells_ind]*alpha_ez_zp[ncells_ind]));
				cpml_b_mz_zp[ncells_ind] = exp((-dt/mu_0)*((sigma_pmz_zp[ncells_ind]/kappa_mz_zp[ncells_ind]) + alpha_mz_zp[ncells_ind]));
				cpml_a_mz_zp[ncells_ind] = (1/dz)*(cpml_b_mz_zp[ncells_ind] - 1.0)*sigma_pmz_zp[ncells_ind]/(kappa_mz_zp[ncells_ind]*(sigma_pmz_zp[ncells_ind] + kappa_mz_zp[ncells_ind]*alpha_mz_zp[ncells_ind]));
			}

			//Create and initialize 2D cpml convolution parameters 
			Psi_exz_zp = build_3D_double_array(nx,nyp1,ncells,0.0); 
			Psi_eyz_zp = build_3D_double_array(nxp1,ny,ncells,0.0); 
			Psi_hxz_zp = build_3D_double_array(nxp1,ny,ncells,0.0); 
			Psi_hyz_zp = build_3D_double_array(nx,nyp1,ncells,0.0); 

			//Create and initialize 2D cpml convolution coefficients
			
			CPsi_eyz_zp = build_3D_double_array(nxp1,ny,ncells,0.0); 
			CPsi_hxz_zp = build_3D_double_array(nxp1,ny,ncells,0.0); 
			CPsi_exz_zp = build_3D_double_array(nx,nyp1,ncells,0.0); 
			CPsi_hyz_zp = build_3D_double_array(nx,nyp1,ncells,0.0); 

			for(i=0;i<nxp1;i++){
				for(j=0;j<nyp1;j++){
					for(k=0;k<ncells;k++){ 
						if(j<ny){	
							CPsi_eyz_zp[i][j][k] = Ceyhx[i][j][nz - ncells + k]*dz;
							Ceyhx[i][j][nz - ncells + k] = Ceyhx[i][j][nz - ncells + k]/kappa_ez_zp[k];
							CPsi_hxz_zp[i][j][k] = Chxey[i][j][nz - ncells + k]*dz;
							Chxey[i][j][nz - ncells + k] = Chxey[i][j][nz - ncells + k]/kappa_mz_zp[k];
						}
						if(i<nx){
							CPsi_exz_zp[i][j][k] = Cexhy[i][j][nz - ncells + k]*dz;
							Cexhy[i][j][nz - ncells + k] = Cexhy[i][j][nz - ncells + k]/kappa_ez_zp[k];
							CPsi_hyz_zp[i][j][k] = Chyex[i][j][nz - ncells + k]*dz;
							Chyex[i][j][nz - ncells + k] = Chyex[i][j][nz - ncells + k]/kappa_mz_zp[k];
						}
					}
				}
			}
		}
	}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	printf("initializing the output parameters...\n");

	
	//intialize frequency domain parameters
	frequency_domain.number_of_frequencies = (int)((frequency_domain.end - frequency_domain.start)/(frequency_domain.step));
	frequency_domain.frequencies = build_1D_double_array(frequency_domain.number_of_frequencies,0.0);

	int sample_time_ind;
	//initialize sampled electric field terms
	sampled_fields_t *sampled_electric_fields;
	
	if (number_of_sampled_electric_fields > 0){
		sampled_electric_fields = (sampled_fields_t*)malloc(number_of_sampled_electric_fields*sizeof(sampled_fields_t));
		out_of_memory_status(sampled_electric_fields);
		
		for (ind = 0; ind < number_of_sampled_electric_fields; ind++){  
			sampled_electric_fields[ind].is = (int)((sampled_electric_fields[ind].x - fdtd_domain_min_x)/dx)+1;
			sampled_electric_fields[ind].js = (int)((sampled_electric_fields[ind].y - fdtd_domain_min_y)/dy)+1;
			sampled_electric_fields[ind].ks = (int)((sampled_electric_fields[ind].z - fdtd_domain_min_z)/dz)+1;
			sampled_electric_fields[ind].sampled_value = build_1D_double_array(number_of_time_steps,0.0);
			sampled_electric_fields[ind].time = build_1D_double_array(number_of_time_steps,0.0);

			for (sample_time_ind = 0; sample_time_ind < number_of_time_steps; sample_time_ind++){
				sampled_electric_fields[ind].time[sample_time_ind] = (sample_time_ind + 1)*dt;
			}
		}
	}

	//initialize sampled magnetic field terms

	sampled_fields_t *sampled_magnetic_fields;

	if (number_of_sampled_magnetic_fields > 0){
		sampled_magnetic_fields = (sampled_fields_t*)malloc(number_of_sampled_magnetic_fields*sizeof(sampled_fields_t));
		out_of_memory_status(sampled_magnetic_fields);

		for (ind = 0; ind < number_of_sampled_magnetic_fields; ind++){  
			sampled_magnetic_fields[ind].is = (int)((sampled_magnetic_fields[ind].x - fdtd_domain_min_x)/dx)+1;
			sampled_magnetic_fields[ind].js = (int)((sampled_magnetic_fields[ind].y - fdtd_domain_min_y)/dy)+1;
			sampled_magnetic_fields[ind].ks = (int)((sampled_magnetic_fields[ind].z - fdtd_domain_min_z)/dz)+1;
			sampled_magnetic_fields[ind].sampled_value = build_1D_double_array(number_of_time_steps,0.0);
			sampled_magnetic_fields[ind].time = build_1D_double_array(number_of_time_steps,0.0);

			for (sample_time_ind = 0; sample_time_ind < number_of_time_steps; sample_time_ind++){
				sampled_magnetic_fields[ind].time[sample_time_ind] = (sample_time_ind + 0.5)*dt;
			}
		}
	}

	//initialize sampled voltage terms
	//
	if (number_of_sampled_voltages > 0){
		for (ind=0; ind < number_of_sampled_voltages; ind++){
			sampled_voltages[ind].is = (int)((sampled_voltages[ind].min_x - fdtd_domain_min_x)/dx)+1;
			sampled_voltages[ind].js = (int)((sampled_voltages[ind].min_y - fdtd_domain_min_y)/dy)+1;
			sampled_voltages[ind].ks = (int)((sampled_voltages[ind].min_z - fdtd_domain_min_z)/dz)+1;
			sampled_voltages[ind].ie = (int)((sampled_voltages[ind].max_x - fdtd_domain_min_x)/dx)+1;
			sampled_voltages[ind].je = (int)((sampled_voltages[ind].max_y - fdtd_domain_min_y)/dy)+1;
			sampled_voltages[ind].ke = (int)((sampled_voltages[ind].max_z - fdtd_domain_min_z)/dz)+1;
			sampled_voltages[ind].sampled_value = build_1D_double_array(number_of_time_steps,0.0);
			sampled_voltages[ind].time = build_1D_double_array(number_of_time_steps,0.0);

			if (sampled_voltages[ind].direction == XN || sampled_voltages[ind].direction == XP){
				sampled_voltages[ind].Csvf = -dx/((sampled_voltages[ind].je - sampled_voltages[ind].js + 1)*(sampled_voltages[ind].ke - sampled_voltages[ind].ks + 1));
			}
			if (sampled_voltages[ind].direction == YN || sampled_voltages[ind].direction == YP){
				sampled_voltages[ind].Csvf = -dy/((sampled_voltages[ind].ke - sampled_voltages[ind].ks + 1)*(sampled_voltages[ind].ie - sampled_voltages[ind].is + 1));
			}
			if (sampled_voltages[ind].direction == ZN || sampled_voltages[ind].direction == ZP){
				sampled_voltages[ind].Csvf = -dz/((sampled_voltages[ind].ie - sampled_voltages[ind].is + 1)*(sampled_voltages[ind].je - sampled_voltages[ind].js + 1));
			}
			if (sampled_voltages[ind].direction == XN || sampled_voltages[ind].direction == YN || sampled_voltages[ind].direction == ZN ){
				sampled_voltages[ind].Csvf =  -1 * sampled_voltages[ind].Csvf;
			}

			for (sample_time_ind = 0; sample_time_ind < number_of_time_steps; sample_time_ind++){
				sampled_voltages[ind].time[sample_time_ind] = (sample_time_ind + 1)*dt;
			}
		}
	}

	//initialize sampled current terms
	if (number_of_sampled_currents > 0){
		for (ind=1; ind < number_of_sampled_currents; ind++){  
			sampled_currents[ind].is = (int)((sampled_currents[ind].min_x - fdtd_domain_min_x)/dx)+1;
			sampled_currents[ind].js = (int)((sampled_currents[ind].min_y - fdtd_domain_min_y)/dy)+1;
			sampled_currents[ind].ks = (int)((sampled_currents[ind].min_z - fdtd_domain_min_z)/dz)+1;
			sampled_currents[ind].ie = (int)((sampled_currents[ind].max_x - fdtd_domain_min_x)/dx)+1;
			sampled_currents[ind].je = (int)((sampled_currents[ind].max_y - fdtd_domain_min_y)/dy)+1;
			sampled_currents[ind].ke = (int)((sampled_currents[ind].max_z - fdtd_domain_min_z)/dz)+1;
			sampled_currents[ind].sampled_value = build_1D_double_array(number_of_time_steps,0.0);
			sampled_currents[ind].time = build_1D_double_array(number_of_time_steps,0.0);

			for (sample_time_ind = 0; sample_time_ind < number_of_time_steps; sample_time_ind++){
				sampled_currents[ind].time[sample_time_ind] = (sample_time_ind + 0.5)*dt;
			}
		}
	}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//initialize farfield arrays
	double ****tjxyp;
	double ****tjxzp;
	double ****tjyxp;
	double ****tjyzp;
	double ****tjzxp;
	double ****tjzyp;
	double ****tjxyn;
	double ****tjxzn;
	double ****tjyxn;
	double ****tjyzn;
	double ****tjzxn;
	double ****tjzyn;
	double ****tmxyp;
	double ****tmxzp;
	double ****tmyxp;
	double ****tmyzp;
	double ****tmzxp;
	double ****tmzyp;
	double ****tmxyn;
	double ****tmxzn;
	double ****tmyxn;
	double ****tmyzn;
	double ****tmzxn;
	double ****tmzyn;
	double complex ****cjxyp;
	double complex ****cjxzp;
	double complex ****cjyxp;
	double complex ****cjyzp;
	double complex ****cjzxp;
	double complex ****cjzyp;
	double complex ****cjxyn;
	double complex ****cjxzn;
	double complex ****cjyxn;
	double complex ****cjyzn;
	double complex ****cjzxn;
	double complex ****cjzyn;
	double complex ****cmxyp;
	double complex ****cmxzp;
	double complex ****cmyxp;
	double complex ****cmyzp;
	double complex ****cmzxp;
	double complex ****cmzyp;
	double complex ****cmxyn;
	double complex ****cmxzn;
	double complex ****cmyxn;
	double complex ****cmyzn;
	double complex ****cmzxn;
	double complex ****cmzyn;

	int nc_farbuffer = farfield.number_of_cells_from_outer_boundary; 
	int li = nc_farbuffer + 1;
	int lj = nc_farbuffer + 1;
	int lk = nc_farbuffer + 1;
	int ui = nx - nc_farbuffer+1;
	int uj = ny - nc_farbuffer+1;
	int uk = nz - nc_farbuffer+1;
		
	double *farfield_w;

	if (number_of_farfield_frequencies > 0){
		farfield_w = build_1D_double_array(number_of_farfield_frequencies,0.0);

		for (ind = 0; ind < number_of_farfield_frequencies; ind++){
			farfield_w[ind] = 2*M_PI*farfield.frequencies[ind];
		}
		
		tjxyp = build_4D_double_array(1,ui-li,1,uk-lk, 0.0);
		out_of_memory_status(tjxyp);
		tjxzp = build_4D_double_array(1,ui-li,uj-lj,1, 0.0);
		out_of_memory_status(tjxzp);
		tjyxp = build_4D_double_array(1,1,uj-lj,uk-lk, 0.0);
		out_of_memory_status(tjyxp);
		tjyzp = build_4D_double_array(1,ui-li,uj-lj,1, 0.0);
		out_of_memory_status(tjyzp);
		tjzxp = build_4D_double_array(1,1,uj-lj,uk-lk, 0.0);
		out_of_memory_status(tjzxp);
		tjzyp = build_4D_double_array(1,ui-li,1,uk-lk, 0.0);
		out_of_memory_status(tjzyp);
		tjxyn = build_4D_double_array(1,ui-li,1,uk-lk, 0.0);
		out_of_memory_status(tjxyn);
		tjxzn = build_4D_double_array(1,ui-li,uj-lj,1, 0.0);
		out_of_memory_status(tjxzn);
		tjyxn = build_4D_double_array(1,1,uj-lj,uk-lk, 0.0);
		out_of_memory_status(tjyxn);
		tjyzn = build_4D_double_array(1,ui-li,uj-lj,1, 0.0);
		out_of_memory_status(tjyzn);
		tjzxn = build_4D_double_array(1,1,uj-lj,uk-lk, 0.0);
		out_of_memory_status(tjzxn);
		tjzyn = build_4D_double_array(1,ui-li,1,uk-lk, 0.0);
		out_of_memory_status(tjzyn);
		tmxyp = build_4D_double_array(1,ui-li,1,uk-lk, 0.0);
		out_of_memory_status(tmxyp);
		tmxzp = build_4D_double_array(1,ui-li,uj-lj,1, 0.0);
		out_of_memory_status(tmxzp);
		tmyxp = build_4D_double_array(1,1,uj-lj,uk-lk, 0.0);
		out_of_memory_status(tmyxp);
		tmyzp = build_4D_double_array(1,ui-li,uj-lj,1, 0.0);
		out_of_memory_status(tmyzp);
		tmzxp = build_4D_double_array(1,1,uj-lj,uk-lk, 0.0);
		out_of_memory_status(tmzxp);
		tmzyp = build_4D_double_array(1,ui-li,1,uk-lk, 0.0);
		out_of_memory_status(tmzyp);
		tmxyn = build_4D_double_array(1,ui-li,1,uk-lk, 0.0);
		out_of_memory_status(tmxyn);
		tmxzn = build_4D_double_array(1,ui-li,uj-lj,1, 0.0);
		out_of_memory_status(tmxzn);
		tmyxn = build_4D_double_array(1,1,uj-lj,uk-lk, 0.0);
		out_of_memory_status(tmyxn);
		tmyzn = build_4D_double_array(1,ui-li,uj-lj,1, 0.0);
		out_of_memory_status(tmyzn);
		tmzxn = build_4D_double_array(1,1,uj-lj,uk-lk, 0.0);
		out_of_memory_status(tmzxn);
		tmzyn = build_4D_double_array(1,ui-li,1,uk-lk, 0.0);	
		out_of_memory_status(tmzyn);
		cjxyp = build_4D_double_complex_array(ui-li,1,uk-lk,number_of_farfield_frequencies, 0.0 + 0.0*I);
		out_of_memory_status(cjxyp);
		cjxzp = build_4D_double_complex_array(ui-li,uj-lj,1,number_of_farfield_frequencies, 0.0 + 0.0*I);
		out_of_memory_status(cjxzp);
		cjyxp = build_4D_double_complex_array(1,uj-lj,uk-lk,number_of_farfield_frequencies, 0.0 + 0.0*I);
		out_of_memory_status(cjyxp);
		cjyzp = build_4D_double_complex_array(ui-li,uj-lj,1,number_of_farfield_frequencies, 0.0 + 0.0*I);
		out_of_memory_status(cjyzp);
		cjzxp = build_4D_double_complex_array(1,uj-lj,uk-lk,number_of_farfield_frequencies, 0.0 + 0.0*I);
		out_of_memory_status(cjzxp);
		cjzyp = build_4D_double_complex_array(ui-li,1,uk-lk,number_of_farfield_frequencies, 0.0 + 0.0*I);
		out_of_memory_status(cjzyp);
		cjxyn = build_4D_double_complex_array(ui-li,1,uk-lk,number_of_farfield_frequencies, 0.0 + 0.0*I);
		out_of_memory_status(cjxyn);
		cjxzn = build_4D_double_complex_array(ui-li,uj-lj,1,number_of_farfield_frequencies, 0.0 + 0.0*I);
		out_of_memory_status(cjxzn);
		cjyxn = build_4D_double_complex_array(1,uj-lj,uk-lk,number_of_farfield_frequencies, 0.0 + 0.0*I);
		out_of_memory_status(cjyxn);
		cjyzn = build_4D_double_complex_array(ui-li,uj-lj,1,number_of_farfield_frequencies, 0.0 + 0.0*I);
		out_of_memory_status(cjyzn);
		cjzxn = build_4D_double_complex_array(1,uj-lj,uk-lk,number_of_farfield_frequencies, 0.0 + 0.0*I);
		out_of_memory_status(cjzxn);
		cjzyn = build_4D_double_complex_array(ui-li,1,uk-lk,number_of_farfield_frequencies, 0.0 + 0.0*I);
		out_of_memory_status(cjzyn);
		cmxyp = build_4D_double_complex_array(ui-li,1,uk-lk,number_of_farfield_frequencies, 0.0 + 0.0*I);
		out_of_memory_status(cmxyp);
		cmxzp = build_4D_double_complex_array(ui-li,uj-lj,1,number_of_farfield_frequencies, 0.0 + 0.0*I);
		out_of_memory_status(cmxzp);
		cmyxp = build_4D_double_complex_array(1,uj-lj,uk-lk,number_of_farfield_frequencies, 0.0 + 0.0*I);
		out_of_memory_status(cmyxp);
		cmyzp = build_4D_double_complex_array(ui-li,uj-lj,1,number_of_farfield_frequencies, 0.0 + 0.0*I);
		out_of_memory_status(cmyzp);
		cmzxp = build_4D_double_complex_array(1,uj-lj,uk-lk,number_of_farfield_frequencies, 0.0 + 0.0*I);
		out_of_memory_status(cmzxp);
		cmzyp = build_4D_double_complex_array(ui-li,1,uk-lk,number_of_farfield_frequencies, 0.0 + 0.0*I);
		out_of_memory_status(cmzyp);
		cmxyn = build_4D_double_complex_array(ui-li,1,uk-lk,number_of_farfield_frequencies, 0.0 + 0.0*I);
		out_of_memory_status(cmxyn);
		cmxzn = build_4D_double_complex_array(ui-li,uj-lj,1,number_of_farfield_frequencies, 0.0 + 0.0*I);
		out_of_memory_status(cmxzn);
		cmyxn = build_4D_double_complex_array(1,uj-lj,uk-lk,number_of_farfield_frequencies, 0.0 + 0.0*I);
		out_of_memory_status(cmyxn);
		cmyzn = build_4D_double_complex_array(ui-li,uj-lj,1,number_of_farfield_frequencies, 0.0 + 0.0*I);
		out_of_memory_status(cmyzn);
		cmzxn = build_4D_double_complex_array(1,uj-lj,uk-lk,number_of_farfield_frequencies, 0.0 + 0.0*I);
		out_of_memory_status(cmzxn);
		cmzyn = build_4D_double_complex_array(ui-li,1,uk-lk,number_of_farfield_frequencies, 0.0 + 0.0*I);
		out_of_memory_status(cmzyn);
	}

	//Delete temporary arrays. These arrays will not be used any more.

	free(sigma_pex_xn);
	free(sigma_pmx_xn);
	free(kappa_ex_xn);
	free(kappa_mx_xn);
	free(alpha_ex_xn);
	free(alpha_mx_xn);
	free(sigma_pex_xp);
	free(sigma_pmx_xp);
	free(kappa_ex_xp);
	free(kappa_mx_xp);
	free(alpha_ex_xp);
	free(alpha_mx_xp);
	free(sigma_pey_yn);
	free(sigma_pmy_yn);
	free(kappa_ey_yn);
	free(kappa_my_yn);
	free(alpha_ey_yn);
	free(alpha_my_yn);
	free(sigma_pey_yp);
	free(sigma_pmy_yp);
	free(kappa_ey_yp);
	free(kappa_my_yp);
	free(alpha_ey_yp);
	free(alpha_my_yp);
	free(sigma_pez_zn);
	free(sigma_pmz_zn);
	free(kappa_ez_zn);
	free(kappa_mz_zn);
	free(alpha_ez_zn);
	free(alpha_mz_zn);
	free(sigma_pez_zp);
	free(sigma_pmz_zp);
	free(kappa_ez_zp);
	free(kappa_mz_zp);
	free(alpha_ez_zp);
	free(alpha_mz_zp);

	free(eps_r_x);
	free(eps_r_y);
	free(eps_r_z);
	free(mu_r_x);
	free(mu_r_y);
	free(mu_r_z);
	free(sigma_e_x);
	free(sigma_e_y);
	free(sigma_e_z);
	free(sigma_m_x);
	free(sigma_m_y);
	free(sigma_m_z);

	printf("Starting the time marching loop...\n");

	double current_time = 0.0;

	//////////////////////////////
	double svx;
	double svy;
	double svz;
	double sampled_value;
	int n_st;
	register unsigned int mi;

	clock_t start, end;
	double cpu_time_used;


	double complex zi = 0.0 + 1.0*I;
	double complex *exp_h;
	exp_h = build_1D_double_complex_array(number_of_farfield_frequencies, 0.0+0.0*I);
	double complex *exp_e;
	exp_e = build_1D_double_complex_array(number_of_farfield_frequencies, 0.0+0.0*I);

	/////////////////////////////
	int time_step;
	for (time_step = 0; time_step < number_of_time_steps; time_step++){
		start = clock();
		//////////////////// update magnetic fields////////////////////////////
		current_time  = current_time + dt/2;
		for (i=0; i<nxp1; i++){
			for (j=0; j<ny; j++){
				for (k=0; k<nz; k++){	   
					Hx[i][j][k] = Chxh[i][j][k]*Hx[i][j][k] + Chxey[i][j][k]*(Ey[i][j][k+1] - Ey[i][j][k]) + Chxez[i][j][k]*(Ez[i][j+1][k] - Ez[i][j][k]);
				}
			}
		}
		for (i=0; i<nx; i++){
			for (j=0; j<nyp1; j++){
				for (k=0; k<nz; k++){
					Hy[i][j][k] = Chyh[i][j][k]*Hy[i][j][k] + Chyez[i][j][k]*(Ez[i+1][j][k] - Ez[i][j][k]) + Chyex[i][j][k]*(Ex[i][j][k+1] - Ex[i][j][k]);
				}
			}
		}
		for (i=0; i<nx; i++){
			for (j=0; j<ny; j++){
				for (k=0; k<nzp1; k++){
					Hz[i][j][k] = Chzh[i][j][k]*Hz[i][j][k] + Chzex[i][j][k]*(Ex[i][j+1][k] - Ex[i][j][k]) + Chzey[i][j][k]*(Ey[i+1][j][k] - Ey[i][j][k]);
				}
			}
		}

		///////////////apply CPML to magnetic field components////////////////
		//is_cpml_xn == TRUE && k<nz
		for (i = 0; i<n_cpml_xn; i++){
			for(j = 0; j<nyp1; j++){
				for(k = 0; k<nz; k++){
					Psi_hyx_xn[i][j][k] = cpml_b_mx_xn[i]*Psi_hyx_xn[i][j][k] + cpml_a_mx_xn[i]*(Ez[i+1][j][k] - Ez[i][j][k]); 
					Hy[i][j][k] = Hy[i][j][k] + CPsi_hyx_xn[i][j][k]*Psi_hyx_xn[i][j][k];
				}
			}
		}

		//is_cpml_xn == TRUE && j<ny
		for (i = 0; i<n_cpml_xn; i++){
			for(j = 0; j<ny; j++){
				for(k = 0; k<nzp1; k++){
					Psi_hzx_xn[i][j][k] = cpml_b_mx_xn[i]*Psi_hzx_xn[i][j][k] + cpml_a_mx_xn[i]*(Ey[i+1][j][k] - Ey[i][j][k]); 
					Hz[i][j][k] = Hz[i][j][k] + CPsi_hzx_xn[i][j][k]*Psi_hzx_xn[i][j][k];    
				}
			}
		}

		//is_cpml_xp == TRUE && k<nz
		for (i = 0; i<n_cpml_xp; i++){
			for(j = 0; j<nyp1; j++){
				for(k = 0; k<nz; k++){
					Psi_hyx_xp[i][j][k] = cpml_b_mx_xp[i]*Psi_hyx_xp[i][j][k] + cpml_a_mx_xp[i]*(Ez[i+(nx - n_cpml_xp)+1][j][k] - Ez[i+(nx - n_cpml_xp)][j][k]); 
					Hy[(nx - n_cpml_xp)+i][j][k] = Hy[(nx - n_cpml_xp)+i][j][k] + CPsi_hyx_xp[i][j][k]*Psi_hyx_xp[i][j][k];
				}
			}
		}

		//is_cpml_xp == TRUE && j<ny
		for (i = 0; i<n_cpml_xp; i++){
			for(j = 0; j<ny; j++){
				for(k = 0; k<nzp1; k++){
					Psi_hzx_xp[i][j][k] = cpml_b_mx_xp[i]*Psi_hzx_xp[i][j][k] + cpml_a_mx_xp[i]*(Ey[i+(nx - n_cpml_xp)+1][j][k] - Ey[i+(nx - n_cpml_xp)][j][k]); 
					Hz[(nx - n_cpml_xp)+i][j][k] = Hz[(nx - n_cpml_xp)+i][j][k] + CPsi_hzx_xp[i][j][k]*Psi_hzx_xp[i][j][k];    
				}
			}
		}

		//is_cpml_yn == TRUE && i<nx
		for(i=0; i<nx; i++){
			for (j = 0; j<n_cpml_yn;j++){
				for (k = 0; k<nzp1; k++){
					Psi_hzy_yn[i][j][k] = cpml_b_my_yn[j]*Psi_hzy_yn[i][j][k] + cpml_a_my_yn[j]*(Ex[i][j+1][k] - Ex[i][j][k]); 
					Hz[i][j][k] = Hz[i][j][k] + CPsi_hzy_yn[i][j][k]*Psi_hzy_yn[i][j][k];
				}
			}
		}

		//is_cpml_yn == TRUE && k<nz
		for(i=0; i<nxp1; i++){
			for (j = 0; j<n_cpml_yn;j++){
				for (k = 0; k<nz; k++){
					Psi_hxy_yn[i][j][k] = cpml_b_my_yn[j]*Psi_hxy_yn[i][j][k] + cpml_a_my_yn[j]*(Ez[i][j+1][k] - Ez[i][j][k]); 
					Hx[i][j][k] = Hx[i][j][k] + CPsi_hxy_yn[i][j][k]*Psi_hxy_yn[i][j][k];
				}
			}
		}

		//is_cpml_yp == TRUE && i<nx
		for (i=0;i<nx;i++){
			for (j = 0; j<n_cpml_yp;j++){
				for (k = 0;k<nzp1;k++){
					Psi_hzy_yp[i][j][k] = cpml_b_my_yp[j] * Psi_hzy_yp[i][j][k] + cpml_a_my_yp[j]*(Ex[i][j+(ny - n_cpml_yp)+1][k] - Ex[i][j+(ny - n_cpml_yp)][k]); 
					Hz[i][(ny - n_cpml_yp)+j][k] = Hz[i][(ny - n_cpml_yp)+j][k] + CPsi_hzy_yp[i][j][k]*Psi_hzy_yp[i][j][k];
				}
			}
		}

		//is_cpml_yp == TRUE && k<nz
		for (i=0;i<nxp1;i++){
			for (j = 0; j<n_cpml_yp;j++){
				for (k = 0;k<nz;k++){
					Psi_hxy_yp[i][j][k] = cpml_b_my_yp[i] * Psi_hxy_yp[i][j][k] + cpml_a_my_yp[i]*(Ez[i][j+(ny - n_cpml_yp)+1][k] - Ez[i][j+(ny - n_cpml_yp)][k]); 
					Hx[i][(ny - n_cpml_yp)+j][k] = Hx[i][(ny - n_cpml_yp)+j][k] + CPsi_hxy_yp[i][j][k]*Psi_hxy_yp[i][j][k];    
				}
			}
		}

		//is_cpml_zn == TRUE && j<ny
		for (i = 0; i<nxp1;i++){
			for (j = 0; j<ny;j++){
				for (k = 0; k<n_cpml_zn; k++){
					Psi_hxz_zn[i][j][k] = cpml_b_mz_zn[k] * Psi_hxz_zn[i][j][k] + cpml_a_mz_zn[k]*(Ey[i][j][k+1] - Ey[i][j][k]);
					Hx[i][j][k] = Hx[i][j][k] + CPsi_hxz_zn[i][j][k]*Psi_hxz_zn[i][j][k];
				}
			}
		}

		//is_cpml_zn == TRUE && i<nx
		for (i = 0; i<nx;i++){
			for (j = 0; j<nyp1;j++){
				for (k = 0; k<n_cpml_zn; k++){
					Psi_hyz_zn[i][j][k] = cpml_b_mz_zn[k] * Psi_hyz_zn[i][j][k] + cpml_a_mz_zn[k]*(Ex[i][j][k+1] - Ex[i][j][k]);
					Hy[i][j][k] = Hy[i][j][k] + CPsi_hyz_zn[i][j][k]*Psi_hyz_zn[i][j][k];
				}	
			}
		}

		//is_cpml_zp ==TRUE && j<ny
		for (i = 0; i<nxp1; i++){
			for (j = 0; j<ny; j++){
				for (k = 0; k<n_cpml_zp; k++){
					Psi_hxz_zp[i][j][k] = cpml_b_mz_zp[k] * Psi_hxz_zp[i][j][k] + cpml_a_mz_zp[k]*(Ey[i][j][k+(nz - n_cpml_zp)+1] - Ey[i][j][k+(nz - n_cpml_zp)]); 
					Hx[i][j][(nz - n_cpml_zp)+k] = Hx[i][j][(nz - n_cpml_zp)+k] + CPsi_hxz_zp[i][j][k]*Psi_hxz_zp[i][j][k];
				}
			}
		}

		//is_cpml_zp ==TRUE && i<nx
		for (i = 0; i<nx; i++){
			for (j = 0; j<nyp1; j++){
				for (k = 0; k<n_cpml_zp; k++){					
					Psi_hyz_zp[i][j][k] = cpml_b_mz_zp[k] * Psi_hyz_zp[i][j][k] + cpml_a_mz_zp[k]*(Ex[i][j][k+(nz - n_cpml_zp)+1] - Ex[i][j][k+(nz - n_cpml_zp)]); 
					Hy[i][j][(nz - n_cpml_zp)+k] = Hy[i][j][(nz - n_cpml_zp)+k] + CPsi_hyz_zp[i][j][k]*Psi_hyz_zp[i][j][k];
				}    
			}
		}

		///////////////////////Capturing magnetic fields//////////////////////////////////
		for (ind=0; ind<number_of_sampled_magnetic_fields; ind++){
			sampled_value = 0.0;
			if (sampled_magnetic_fields[ind].component == X){
				for(i=sampled_magnetic_fields[ind].is; i<=sampled_magnetic_fields[ind].is; i++){
					for(j=sampled_magnetic_fields[ind].js-1; j<=sampled_magnetic_fields[ind].js; j++){
						for(k=sampled_magnetic_fields[ind].ks-1; k<=sampled_magnetic_fields[ind].ks;k++){
							sampled_value = sampled_value + Hx[i][j][k];
						}
					}
				}
			}
			else if (sampled_magnetic_fields[ind].component == Y){
				for(i=sampled_magnetic_fields[ind].is-1; i<=sampled_magnetic_fields[ind].is; i++){
					for(j=sampled_magnetic_fields[ind].js; j<=sampled_magnetic_fields[ind].js; j++){
						for(k=sampled_magnetic_fields[ind].ks-1; k<=sampled_magnetic_fields[ind].ks;k++){
							sampled_value = sampled_value + Hy[i][j][k];
						}
					}
				}
			}
			else if (sampled_magnetic_fields[ind].component == Z){
				for(i=sampled_magnetic_fields[ind].is-1; i<=sampled_magnetic_fields[ind].is; i++){
					for(j=sampled_magnetic_fields[ind].js-1; j<=sampled_magnetic_fields[ind].js; j++){
						for(k=sampled_magnetic_fields[ind].ks; k<=sampled_magnetic_fields[ind].ks;k++){
							sampled_value = sampled_value + Hz[i][j][k];
						}
					}
				}
			}	
			else if (sampled_magnetic_fields[ind].component == M){
				svx=svy=svz=0.0;
				for(i=sampled_magnetic_fields[ind].is-1; i<=sampled_magnetic_fields[ind].is; i++){
					for(j=sampled_magnetic_fields[ind].js-1; j<=sampled_magnetic_fields[ind].js; j++){
						for(k=sampled_magnetic_fields[ind].ks-1; k<=sampled_magnetic_fields[ind].ks;k++){
							if (i == sampled_magnetic_fields[ind].is){
								svx = svx + Hx[i+1][j][k];
							}
							if (j == sampled_magnetic_fields[ind].js){
								svy = svy + Hy[i][j+1][k];
							}
							if (k == sampled_magnetic_fields[ind].ks){
								svz = svz + Hz[i][j][k+1];
							}
						}
					}
				}
				svx = svx*svx;
				svy = svy*svy;
				svz = svz*svz;
	
				sampled_value = sqrt(svx + svy + svz);
			}
	
			sampled_magnetic_fields[ind].sampled_value[time_step] = sampled_value;
		}
	
		/////////////////////////Capturing sampled currents///////////////////////////////
		if (number_of_sampled_currents > 0){
			capture_sampled_currents(number_of_sampled_currents, time_step, sampled_currents, Hx, Hy, Hz, unit_cell);
		}
		//////////////update electric fields except the tangential components on the boundaries////////////////// 
	
		current_time  = current_time + dt/2;
	
		for (i=0;i<nx;i++){
			for (j=1;j<ny;j++){
				for (k=1;k<nz;k++){
					Ex[i][j][k] = Cexe[i][j][k]*Ex[i][j][k] + Cexhz[i][j][k]*(Hz[i][j][k] - Hz[i][j-1][k]) + Cexhy[i][j][k]*(Hy[i][j][k] - Hy[i][j][k-1]);   
				}
			}
		}
		for (i=1;i<nx;i++){
			for (j=0;j<ny;j++){
				for (k=1;k<nz;k++){
					Ey[i][j][k] = Ceye[i][j][k]*Ey[i][j][k] + Ceyhx[i][j][k]*(Hx[i][j][k] - Hx[i][j][k-1]) + Ceyhz[i][j][k]*(Hz[i][j][k] - Hz[i-1][j][k]); 
				}
			}
		}
		for (i=1;i<nx;i++){
			for (j=1;j<ny;j++){
				for (k=0;k<nz;k++){
					Ez[i][j][k] = Ceze[i][j][k]*Ez[i][j][k] + Cezhy[i][j][k]*(Hy[i][j][k] - Hy[i-1][j][k]) + Cezhx[i][j][k]*(Hx[i][j][k] - Hx[i][j-1][k]);
				}
			}
		}
		
	
		///////////////////////////////apply CPML to electric field components///////////////////////////////////
		//is_cpml_xn == TRUE && j<ny
		for (i = 0; i<n_cpml_xn; i++){
			for (j = 0; j<ny; j++){
				for (k = 0; k<nzp1;k++){
					Psi_eyx_xn[i][j][k] = cpml_b_ex_xn[i]*Psi_eyx_xn[i][j][k] + cpml_a_ex_xn[i]*(Hz[i+1][j][k] - Hz[i][j][k]);
					Ey[i+1][j][k] = Ey[i+1][j][k] + CPsi_eyx_xn[i][j][k]*Psi_eyx_xn[i][j][k];
				}
			}
		}

		//is_cpml_xn == TRUE && k<nz
		for (i = 0; i<n_cpml_xn; i++){
			for (j = 0; j<nyp1; j++){
				for (k = 0; k<nz;k++){ 
					Psi_ezx_xn[i][j][k] = cpml_b_ex_xn[i]*Psi_ezx_xn[i][j][k] + cpml_a_ex_xn[i]*(Hy[i+1][j][k] - Hy[i][j][k]); 
					Ez[i+1][j][k] = Ez[i+1][j][k] + CPsi_ezx_xn[i][j][k]*Psi_ezx_xn[i][j][k];
				}
			}
		}
	
		//is_cpml_xp == TRUE && j<ny
		n_st = nx - n_cpml_xp;
		for (i = 0; i<n_cpml_xp; i++){
			for (j = 0; j<ny; j++){
				for (k =0; k<nzp1; k++){
					Psi_eyx_xp[i][j][k] = cpml_b_ex_xp[i]*Psi_eyx_xp[i][j][k] + cpml_a_ex_xp[i]*(Hz[i+n_st][j][k] - Hz[i+n_st-1][j][k]);
					Ey[n_st+i][j][k] = Ey[n_st+i][j][k] + CPsi_eyx_xp[i][j][k]*Psi_eyx_xp[i][j][k];
				}
			}
		}

		//is_cpml_xp == TRUE && k<nz
		for (i = 0; i<n_cpml_xp; i++){
			for (j = 0; j<nyp1; j++){
				for (k =0; k<nz; k++){
					Psi_ezx_xp[i][j][k] = cpml_b_ex_xp[i]*Psi_ezx_xp[i][j][k] + cpml_a_ex_xp[i]*(Hy[i+n_st][j][k] - Hy[i+n_st-1][j][k]);
					Ez[n_st+i][j][k] = Ez[n_st+i][j][k] + CPsi_ezx_xp[i][j][k]*Psi_ezx_xp[i][j][k];
				}
			}
		}
	
		//is_cpml_yn == TRUE && k<nz
		for (i = 0; i<nxp1;i++){
		    for (j = 0; j<n_cpml_yn; j++){
			for (k = 0; k<nz; k++){
			    Psi_ezy_yn[i][j][k] = cpml_b_ey_yn[j] * Psi_ezy_yn[i][j][k] + cpml_a_ey_yn[j]*(Hx[i][j+1][k] - Hx[i][j][k]); 
			    Ez[i][j+1][k] = Ez[i][j+1][k] + CPsi_ezy_yn[i][j][k]* Psi_ezy_yn[i][j][k];
			}
		    }
		}
	
		//is_cpml_yn == TRUE && i<nx
		for (i = 0; i<nx;i++){
			for (j = 0; j<n_cpml_yn; j++){
				for (k = 0; k<nzp1; k++){
					Psi_exy_yn[i][j][k] = cpml_b_ey_yn[j] * Psi_exy_yn[i][j][k] + cpml_a_ey_yn[j]*(Hz[i][j+1][k] - Hz[i][j][k]); 
					Ex[i][j+1][k] = Ex[i][j+1][k] + CPsi_exy_yn[i][j][k]* Psi_exy_yn[i][j][k];
				}    
			}
		}
	
		//is_cpml_yp == TRUE && k<nz
		n_st = ny - n_cpml_yp;
		for(i = 0; i<nxp1; i++){
			for (j = 0; j<n_cpml_yp; j++){
				for (k = 0; k<nz; k++){
					Psi_ezy_yp[i][j][k] = cpml_b_ey_yp[j]*Psi_ezy_yp[i][j][k] + cpml_a_ey_yp[j]*(Hx[i][j+n_st][k] - Hx[i][j+n_st-1][k]); 
					Ez[i][n_st+j][k] = Ez[i][n_st+j][k] + CPsi_ezy_yp[i][j][k]*Psi_ezy_yp[i][j][k];
				}
			}
		}

		//is_cpml_yp == TRUE && i<nx
		for(i = 0; i<nx; i++){
			for (j = 0; j<n_cpml_yp; j++){
				for (k = 0; k<nzp1; k++){
					Psi_exy_yp[i][j][k] = cpml_b_ey_yp[j]*Psi_exy_yp[i][j][k] + cpml_a_ey_yp[j]*(Hz[i][j+n_st][k] - Hz[i][j+n_st-1][k]);
					Ex[i][n_st+j][k] = Ex[i][n_st+j][k] + CPsi_exy_yp[i][j][k]*Psi_exy_yp[i][j][k];
				}
			}
		}
	
		//is_cpml_zn == TRUE && i<nx
		for (i =0; i<nx;i++){
			for (j=0;j<nyp1;j++){
				for (k = 0; k<n_cpml_zn; k++){
					Psi_exz_zn[i][j][k] = cpml_b_ez_zn[k] * Psi_exz_zn[i][j][k] + cpml_a_ez_zn[k]*(Hy[i][j][k+1] - Hy[i][j][k]); 
					Ex[i][j][k+1] = Ex[i][j][k+1] + CPsi_exz_zn[i][j][k]*Psi_exz_zn[i][j][k];
				}
			}
		}

		//is_cpml_zn == TRUE && j<ny
		for (i =0; i<nxp1;i++){
			for (j=0;j<ny;j++){
				for (k = 0; k<n_cpml_zn; k++){
					Psi_eyz_zn[i][j][k] = cpml_b_ez_zn[k] * Psi_eyz_zn[i][j][k] + cpml_a_ez_zn[k]*(Hx[i][j][k+1] - Hx[i][j][k]); 
					Ey[i][j][k+1] = Ey[i][j][k+1] + CPsi_eyz_zn[i][j][k]*Psi_eyz_zn[i][j][k];    
				}
			}
		}
	
		//is_cpml_zp = TRUE && i<nx
		n_st = nz - n_cpml_zp;
		for (i = 0; i<nx;i++){
			for (j = 0; j<nyp1; j++){
				for (k = 0; k<n_cpml_zp; k++){
					Psi_exz_zp[i][j][k] = cpml_b_ez_zp[k]*Psi_exz_zp[i][j][k] + cpml_a_ez_zp[k]*(Hy[i][j][k+n_st]-Hy[i][j][k+n_st-1]); 	
					Ex[i][j][n_st+k] = Ex[i][j][n_st+k] + CPsi_exz_zp[i][j][k]*Psi_exz_zp[i][j][k];
				}
			}
		}

		//is_cpml_zp = TRUE && j<ny
		for (i = 0; i<nxp1;i++){
			for (j = 0; j<ny; j++){
				for (k = 0; k<n_cpml_zp; k++){
					Psi_eyz_zp[i][j][k] = cpml_b_ez_zp[k]*Psi_eyz_zp[i][j][k] + cpml_a_ez_zp[k]*(Hx[i][j][k+n_st]-Hx[i][j][k+n_st-1]);
					Ey[i][j][n_st+k] = Ey[i][j][n_st+k] + CPsi_eyz_zp[i][j][k]*Psi_eyz_zp[i][j][k];
				}
			}
		}
	
	///////////updating electric field components associated with the voltage sources/////////////////
	
		for (ind = 0; ind < number_of_voltage_sources; ind++){
			is = voltage_sources[ind].is;
			js = voltage_sources[ind].js;
			ks = voltage_sources[ind].ks;
			ie = voltage_sources[ind].ie;
			je = voltage_sources[ind].je;
			ke = voltage_sources[ind].ke;
	
			if (voltage_sources[ind].direction == XN || voltage_sources[ind].direction == XP){
				for (i=is; i < ie; i++){
					for (j=js; j < je; j++){
						for (k=ks; k < ke; k++){
							Ex[i][j][k] = Ex[i][j][k] + voltage_sources[ind].Cexs[i][j][k]*voltage_sources[ind].voltage_per_e_field[time_step];
						}
					}
				}
			}
			if (voltage_sources[ind].direction == YN || voltage_sources[ind].direction == YP){
				for (i=is; i < ie; i++){
					for (j=js; j < je; j++){
						for (k=ks; k < ke; k++){
							Ey[i][j][k] = Ey[i][j][k] + voltage_sources[ind].Ceys[i][j][k]*voltage_sources[ind].voltage_per_e_field[time_step];
						}
					}
				}
			}
			if (voltage_sources[ind].direction == ZN || voltage_sources[ind].direction == ZP){
				for (i=is; i <= ie; i++){
					for (j=js; j <= je; j++){
						for (k=ks; k <= ke; k++){
							Ez[i][j][k] = Ez[i][j][k] + voltage_sources[ind].Cezs[i-is][j-js][k-ks]*voltage_sources[ind].voltage_per_e_field[time_step];
						}
					}
				}
			}
		}
	
		/////////////////////Capturing electric fields/////////////////////////
		
		capture_electric_fields(number_of_sampled_electric_fields, time_step, sampled_electric_fields, Ex, Ey, Ez);
	
		////////////////////Capturing sampled voltages////////////////////////
	
		capture_sampled_voltages(number_of_sampled_voltages, time_step , sampled_voltages, Ex, Ey, Ez);
	
	
		//////////Calculate J and M on the imaginary farfiled surface/////////
	
		 
		for (mi=0; mi<number_of_farfield_frequencies; mi++){
			exp_h[mi] = dt*exp_complex(-zi*farfield_w[mi]*(time_step - 0.5)*dt);
			exp_e[mi] = dt*exp_complex(-zi*farfield_w[mi]*time_step*dt);
		}

		if (number_of_farfield_frequencies > 0){	
			for(i = li; i < ui; i++){
				for (j = lj; j < uj; j++){
					for (k = lk; k < uk; k++){
						if (i == ui){
							tmyxp[0][0][j-lj][k-lk] =  0.5*(Ez[ui][j][k] + Ez[ui][j+1][k]);
							tmzxp[0][0][j-lj][k-lk] = -0.5*(Ey[ui][j][k] + Ey[ui][j][k+1]);
		
							tjyxp[0][0][j-lj][k-lk] = -0.25*(Hz[ui][j][k] + Hz[ui][j][k+1] + Hz[ui-1][j][k] + Hz[ui-1][j][k+1]);
							tjzxp[0][0][j-lj][k-lk] =  0.25*(Hy[ui][j][k] + Hy[ui][j+1][k] + Hy[ui-1][j][k] + Hy[ui-1][j+1][k]);
		
							tmyxn[0][0][j-lj][k-lk] = -0.5*(Ez[li][j][k] + Ez[li][j+1][k]);
							tmzxn[0][0][j-lj][k-lk] =  0.5*(Ey[li][j][k] + Ey[li][j][k+1]);
		
							tjyxn[0][0][j-lj][k-lk] =  0.25*(Hz[li][j][k] + Hz[li][j][k+1] + Hz[li-1][j][k] + Hz[li-1][j][k+1]);
							tjzxn[0][0][j-lj][k-lk] = -0.25*(Hy[li][j][k] + Hy[li][j+1][k] + Hy[li-1][j][k] + Hy[li-1][j+1][k]);
		
							//fourier transform
							for (mi=0; mi<number_of_farfield_frequencies; mi++){
								cjyxp[0][j-lj][k-lk][mi] = cjyxp[0][j-lj][k-lk][mi] + exp_h[mi]*tjyxp[0][0][j-lj][k-lk];
								cjzxp[0][j-lj][k-lk][mi] = cjzxp[0][j-lj][k-lk][mi] + exp_h[mi]*tjzxp[0][0][j-lj][k-lk];
		
								cjyxn[0][j-lj][k-lk][mi] = cjyxn[0][j-lj][k-lk][mi] + exp_h[mi]*tjyxn[0][0][j-lj][k-lk];
								cjzxn[0][j-lj][k-lk][mi] = cjzxn[0][j-lj][k-lk][mi] + exp_h[mi]*tjzxn[0][0][j-lj][k-lk];
		
								cmyxp[0][j-lj][k-lk][mi] = cmyxp[0][j-lj][k-lk][mi] + exp_e[mi]*tmyxp[0][0][j-lj][k-lk];
								cmzxp[0][j-lj][k-lk][mi] = cmzxp[0][j-lj][k-lk][mi] + exp_e[mi]*tmzxp[0][0][j-lj][k-lk]; 
		
								cmyxn[0][j-lj][k-lk][mi] = cmyxn[0][j-lj][k-lk][mi] + exp_e[mi]*tmyxn[0][0][j-lj][k-lk];
								cmzxn[0][j-lj][k-lk][mi] = cmzxn[0][j-lj][k-lk][mi] + exp_e[mi]*tmzxn[0][0][j-lj][k-lk]; 
							}
						}
		
						if (j == uj){
							tmxyp[0][i-li][0][k-lk] = -0.5*(Ez[i][uj][k] + Ez[i+1][uj][k]);
							tmzyp[0][i-li][0][k-lk] =  0.5*(Ex[i][uj][k] + Ex[i][uj][k+1]);
		
							tjzyp[0][i-li][0][k-lk] = -0.25*(Hx[i][uj][k] + Hx[i+1][uj][k] + Hx[i][uj-1][k] + Hx[i+1][uj-1][k]);
							tjxyp[0][i-li][0][k-lk] =  0.25*(Hz[i][uj][k] + Hz[i][uj][k+1] + Hz[i][uj-1][k] + Hz[i][uj-1][k+1]);
		
							tmxyn[0][i-li][0][k-lk] =  0.5*(Ez[i][lj][k] + Ez[i+1][lj][k]);
							tmzyn[0][i-li][0][k-lk] = -0.5*(Ex[i][lj][k] + Ex[i][lj][k+1]);
		
							tjzyn[0][i-li][0][k-lk] =  0.25*(Hx[i][lj][k] + Hx[i+1][lj][k] + Hx[i][lj-1][k] + Hx[i+1][lj-1][k]);
							tjxyn[0][i-li][0][k-lk] = -0.25*(Hz[i][lj][k] + Hz[i][lj][k+1] + Hz[i][lj-1][k] + Hz[i][lj-1][k+1]);
		
							//fourier transform
							for (mi=0; mi<number_of_farfield_frequencies; mi++){
								cjxyp[i-li][0][k-lk][mi] = cjxyp[i-li][0][k-lk][mi] + exp_h[mi] * tjxyp[0][i-li][0][k-lk];
								cjzyp[i-li][0][k-lk][mi] = cjzyp[i-li][0][k-lk][mi] + exp_h[mi] * tjzyp[0][i-li][0][k-lk];
		
								cjxyn[i-li][0][k-lk][mi] = cjxyn[i-li][0][k-lk][mi] + exp_h[mi] * tjxyn[0][i-li][0][k-lk];
								cjzyn[i-li][0][k-lk][mi] = cjzyn[i-li][0][k-lk][mi] + exp_h[mi] * tjzyn[0][i-li][0][k-lk];
		
								cmxyp[i-li][0][k-lk][mi] = cmxyp[i-li][0][k-lk][mi] + exp_e[mi] * tmxyp[0][i-li][0][k-lk]; 
								cmzyp[i-li][0][k-lk][mi] = cmzyp[i-li][0][k-lk][mi] + exp_e[mi] * tmzyp[0][i-li][0][k-lk]; 
		
								cmxyn[i-li][0][k-lk][mi] = cmxyn[i-li][0][k-lk][mi] + exp_e[mi] * tmxyn[0][i-li][0][k-lk]; 
								cmzyn[i-li][0][k-lk][mi] = cmzyn[i-li][0][k-lk][mi] + exp_e[mi] * tmzyn[0][i-li][0][k-lk]; 
							}
						}
		
						if (k == uk){
							tmxzp[0][i-li][j-lj][0] =  0.5*(Ey[i][j][uk] + Ey[i+1][j][uk]);
							tmyzp[0][i-li][j-lj][0] = -0.5*(Ex[i][j][uk] + Ex[i][j+1][uk]);
		
							tjyzp[0][i-li][j-lj][0] =  0.25*(Hx[i][j][uk] + Hx[i+1][j][uk] + Hx[i][j][uk-1] + Hx[i+1][j][uk-1]);
							tjxzp[0][i-li][j-lj][0] = -0.25*(Hy[i][j][uk] + Hy[i][j+1][uk] + Hy[i][j][uk-1] + Hy[i][j+1][uk-1]);
		
							tmxzn[0][i-li][j-lj][0] = -0.5*(Ey[i][j][lk] + Ey[i+1][j][lk]);
							tmyzn[0][i-li][j-lj][0] =  0.5*(Ex[i][j][lk] + Ex[i][j+1][lk]);
		
							tjyzn[0][i-li][j-lj][0] = -0.25*(Hx[i][j][lk] + Hx[i+1][j][lk] + Hx[i][j][lk-1] + Hx[i+1][j][lk-1]);
							tjxzn[0][i-li][j-lj][0] =  0.25*(Hy[i][j][lk] + Hy[i][j+1][lk] + Hy[i][j][lk-1] + Hy[i][j+1][lk-1]);
			
							//fourier transform
							for (mi=0; mi<number_of_farfield_frequencies; mi++){
								cjxzp[i-li][j-lj][0][mi] = cjxzp[i-li][j-lj][0][mi] + exp_h[mi] * tjxzp[0][i-li][j-lj][0];
								cjyzp[i-li][j-lj][0][mi] = cjyzp[i-li][j-lj][0][mi] + exp_h[mi] * tjyzp[0][i-li][j-lj][0]; 
		
								cjxzn[i-li][j-lj][0][mi] = cjxzn[i-li][j-lj][0][mi] + exp_h[mi] * tjxzn[0][i-li][j-lj][0];
								cjyzn[i-li][j-lj][0][mi] = cjyzn[i-li][j-lj][0][mi] + exp_h[mi] * tjyzn[0][i-li][j-lj][0];
		
								cmxzp[i-li][j-lj][0][mi] = cmxzp[i-li][j-lj][0][mi] + exp_e[mi] * tmxzp[0][i-li][j-lj][0];
								cmyzp[i-li][j-lj][0][mi] = cmyzp[i-li][j-lj][0][mi] + exp_e[mi] * tmyzp[0][i-li][j-lj][0];
		
								cmxzn[i-li][j-lj][0][mi] = cmxzn[i-li][j-lj][0][mi] + exp_e[mi] * tmxzn[0][i-li][j-lj][0]; 
								cmyzn[i-li][j-lj][0][mi] = cmyzn[i-li][j-lj][0][mi] + exp_e[mi] * tmyzn[0][i-li][j-lj][0]; 
							}
						}
					}
				}
			}
		}
		if ((time_step+1)%100 == 0){
			end = clock();
			cpu_time_used = (number_of_time_steps - time_step -1)*((double) (end - start)) / CLOCKS_PER_SEC;
			printf("**** %d %% Completed, %f Seconds Remaining****\n",(int)((time_step+1)*100/number_of_time_steps), cpu_time_used);
		}	
	}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	printf("generating frequency domain outputs...\n");


	double *frequency_array;
	frequency_array = frequency_domain.frequencies;
	int number_of_frequencies = frequency_domain.number_of_frequencies;
	double time_shift;

	//sampled electric fields in frequency domain
	for (ind=0; ind < number_of_sampled_electric_fields; ind++){
		time_shift = 0.0;
		sampled_electric_fields[ind].frequency_domain_value = time_to_frequency_domain(sampled_electric_fields[ind].sampled_value, dt, frequency_array, time_shift, number_of_time_steps, number_of_frequencies);
		sampled_electric_fields[ind].frequencies = frequency_array;
	
		save_1D_complex_double_array("E_sample_F_domain.csv", sampled_electric_fields[ind].frequency_domain_value, number_of_frequencies);
		
	}

	//sampled magnetic fields in frequency domain
	for (ind=0; ind < number_of_sampled_magnetic_fields; ind++){
		time_shift = -dt/2;
		sampled_magnetic_fields[ind].frequency_domain_value = time_to_frequency_domain(sampled_magnetic_fields[ind].sampled_value, dt, frequency_array, time_shift, number_of_time_steps, number_of_frequencies);
		sampled_magnetic_fields[ind].frequencies = frequency_array;
			
		save_1D_complex_double_array("H_sample_F_domain.csv", sampled_magnetic_fields[ind].frequency_domain_value, number_of_frequencies);
	}

	//sampled voltages in frequency domain
	for (ind=0; ind < number_of_sampled_voltages; ind++){
		time_shift = 0.0;
		sampled_voltages[ind].frequency_domain_value = time_to_frequency_domain(sampled_voltages[ind].sampled_value, dt, frequency_array, time_shift, number_of_time_steps, number_of_frequencies);
		sampled_voltages[ind].frequencies = frequency_array;
			
		save_1D_complex_double_array("V_sample_F_domain.csv", sampled_voltages[ind].frequency_domain_value, number_of_frequencies);
		save_1D_double_array("V_sample_T_domain.csv", sampled_voltages[ind].sampled_value, number_of_time_steps);
	}

	//sampled currents in frequency domain
	for (ind=0; ind < number_of_sampled_currents; ind++){
	    time_shift = -dt/2;
	    sampled_currents[ind].frequency_domain_value = time_to_frequency_domain(sampled_currents[ind].sampled_value, dt, frequency_array, time_shift, number_of_time_steps, number_of_frequencies);
	    sampled_currents[ind].frequencies = frequency_array;
				
		save_1D_complex_double_array("I_sample_F_domain.csv", sampled_currents[ind].frequency_domain_value, number_of_frequencies);
	}

	//voltage sources in frequency domain
	for (ind=0; ind < number_of_voltage_sources; ind++){
		time_shift = 0.0;
		voltage_sources[ind].frequency_domain_value = time_to_frequency_domain(voltage_sources[ind].waveform, dt, frequency_array, time_shift, number_of_time_steps, number_of_frequencies);
		voltage_sources[ind].frequencies = frequency_array;
				
		save_1D_complex_double_array("VS_sample_F_domain.csv", voltage_sources[ind].frequency_domain_value, number_of_frequencies);
	}

	//current sources in frequency domain
	for (ind=0; ind < number_of_current_sources; ind++) {
		time_shift = 0.0;
		current_sources[ind].frequency_domain_value = time_to_frequency_domain(current_sources[ind].waveform, dt, frequency_array, time_shift, number_of_time_steps, number_of_frequencies);
		current_sources[ind].frequencies = frequency_array;
				
		save_1D_complex_double_array("IS_sample_F_domain.csv", current_sources[ind].frequency_domain_value, number_of_frequencies);
	}

	//calculation of S-parameters
	//calculate incident and reflected power waves
	int port_ind1,port_ind2,svi,sci;
	
	for (port_ind1=0; port_ind1 < number_of_ports; port_ind1++){
		svi = ports[port_ind1].sampled_voltage_index; 
		sci = ports[port_ind1].sampled_current_index;
		
		ports[port_ind1].a = build_1D_double_complex_array(number_of_frequencies, 0.0+0.0*I);
		ports[port_ind1].b = build_1D_double_complex_array(number_of_frequencies, 0.0+0.0*I);
 
		for (ind = 0; ind < number_of_frequencies; ind++){
			ports[port_ind1].a[ind] = 0.5*(sampled_voltages[svi].frequency_domain_value[ind] + ports[port_ind1].impedance*sampled_currents[sci].frequency_domain_value[ind])/sqrt(creal(ports[port_ind1].impedance));
			ports[port_ind1].b[ind] = 0.5*(sampled_voltages[svi].frequency_domain_value[ind] - conj(ports[port_ind1].impedance)*sampled_currents[sci].frequency_domain_value[ind])/sqrt(creal(ports[port_ind1].impedance));   
		}
		ports[port_ind1].frequencies = frequency_array;
	}

	//calculate the S-parameters
	for (port_ind1=0; port_ind1 < number_of_ports; port_ind1++){
		if (ports[port_ind1].is_source_port == TRUE) {
			ports[port_ind1].S = build_2D_double_complex_array(number_of_ports,number_of_frequencies, 0.0+0.0*I);
			for (port_ind2=0; port_ind2 < number_of_ports; port_ind2++){
				for (ind = 0; ind < number_of_frequencies; ind++){
					if (ports[port_ind1].a[ind] != 0.0 + 0.0*I){
						ports[port_ind1].S[port_ind2][ind] = ports[port_ind2].b[ind]/ports[port_ind1].a[ind];
					}
					else{
						ports[port_ind1].S[port_ind2][ind] = NAN + NAN*I;
					}
				}
				save_1D_complex_double_array("S_sample.csv", ports[port_ind1].S[port_ind2], number_of_frequencies);
			}
		}
	}
	
	//Calculate total radiated power
	
	double *radiated_power;
	radiated_power = build_1D_double_array(number_of_farfield_frequencies,0.0);

	double complex* sum_of_cxn;
	double complex* sum_of_cxp;
	double complex* sum_of_cyn;
	double complex* sum_of_cyp;
	double complex* sum_of_czn;
	double complex* sum_of_czp;

	if (number_of_farfield_frequencies > 0){
		sum_of_cxn = build_1D_double_complex_array(number_of_farfield_frequencies,0.0+0.0*I);
		sum_of_cxp = build_1D_double_complex_array(number_of_farfield_frequencies,0.0+0.0*I);
		sum_of_cyn = build_1D_double_complex_array(number_of_farfield_frequencies,0.0+0.0*I);
		sum_of_cyp = build_1D_double_complex_array(number_of_farfield_frequencies,0.0+0.0*I);
		sum_of_czn = build_1D_double_complex_array(number_of_farfield_frequencies,0.0+0.0*I);
		sum_of_czp = build_1D_double_complex_array(number_of_farfield_frequencies,0.0+0.0*I);
	}

	for (ind=1; ind < number_of_farfield_frequencies; ind++){

		for (i=0;i<ui-li;i++){
			for (j=0;j<uj-lj;j++){
				sum_of_czp[ind] = sum_of_czp[ind] + cmyzp[i][j][0][ind]*conj(cjxzp[i][j][0][ind]) - cmxzp[i][j][0][ind]*conj(cjyzp[i][j][0][ind]);
				sum_of_czn[ind] = sum_of_czn[ind] + cmyzn[i][j][0][ind]*conj(cjxzn[i][j][0][ind]) - cmxzn[i][j][0][ind]*conj(cjyzn[i][j][0][ind]);
			}
		}

		for (i=0;i<ui-li;i++){
			for (k=0;k<uk-lk;k++){
				sum_of_cyp[ind] = sum_of_cyp[ind] + cmxyp[i][0][k][ind]*conj(cjzyp[i][0][k][ind]) - cmzyp[i][0][k][ind]*conj(cjxyp[i][0][k][ind]);
				sum_of_cyn[ind] = sum_of_cyn[ind] + cmxyn[i][0][k][ind]*conj(cjzyn[i][0][k][ind]) - cmzyn[i][0][k][ind]*conj(cjxyn[i][0][k][ind]);
			}
		}

		for (j=0;j<uj-lj;j++){
			for (k=0;k<uk-lk;k++){
				sum_of_cxp[ind] = sum_of_cxp[ind] + cmzxp[0][j][k][ind]*conj(cjyxp[0][j][k][ind]) - cmyxp[0][j][k][ind]*conj(cjzxp[0][j][k][ind]);
				sum_of_cxn[ind] = sum_of_cxn[ind] + cmzxn[0][j][k][ind]*conj(cjyxn[0][j][k][ind]) - cmyxn[0][j][k][ind]*conj(cjzxn[0][j][k][ind]);
			}
		}

		radiated_power[ind] = 0.5*creal(dx*dy* sum_of_czp[ind] - dx*dy* sum_of_czn[ind] + dx*dz* sum_of_cyp[ind] - dx*dz* sum_of_cyn[ind] + dy*dz* sum_of_cxp[ind] - dy*dz* sum_of_cxn[ind]);

	}

	save_1D_double_array("Radiated_Powre.csv", radiated_power, number_of_frequencies);

	int number_of_angles = 360;
	int step_size = 10;			//increment between the rings in the polar grid
	int Nrings = 4;				//number of rings in the polar grid

	printf("Simulation is Completed!!\n");
    return 0;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  save_3D_int_array
 *  Description:  
 * =====================================================================================
 */

void save_3D_int_array(char *file, int ***array, int size_x, int size_y, int size_z){

	FILE *REC;
	register unsigned int i,j,k;

	REC = fopen(file,"w");
	if(REC == NULL){
		printf("***ERROR*** openning %s file\n", file);
	}
	else{
		for (i = 0; i < size_x; i++){
			for (j = 0; j < size_y; j++){
				for (k = 0; k < size_z; k++){
					fprintf(REC,"%d \n", array[i][j][k]);
				}
			}
		}

		fclose(REC);
		printf("Save %s file...\n", file);
	}
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  save_1D_double_array
 *  Description:  
 * =====================================================================================
 */

void save_1D_double_array(char *file, double *array, int size){

	FILE *REC;
	register unsigned int i;

	REC = fopen(file,"w");

	if(REC == NULL){
		printf("***ERROR*** openning %s file\n", file);
	}
	else{
		for (i = 0; i < size; i++){
			fprintf(REC,"%g \n", array[i]);
		}

		fclose(REC);
		printf("Save %s file...\n", file);
	}
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  save_3D_double_array
 *  Description:  
 * =====================================================================================
 */

void save_3D_double_array(char *file, double ***array, int size_x, int size_y, int size_z){

	FILE *REC;
	register unsigned int i,j,k;

	REC = fopen(file,"w");
	if(REC == NULL){
		printf("***ERROR*** openning %s file\n", file);
	}
	else{
		for (i = 0; i < size_x; i++){
			for (j = 0; j < size_y; j++){
				for (k = 0; k < size_z; k++){
					fprintf(REC,"%10.50f \n", array[i][j][k]);
				}
			}
		}

		fclose(REC);
		printf("Save %s file...\n", file);
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  save_2D_complex_double_array
 *  Description:  
 * =====================================================================================
 */

void save_2D_complex_double_array(char *file, complex double **array, int size_x, int size_y){

	FILE *REC;
	register unsigned int i,j;

	REC = fopen(file,"w");
	if(REC == NULL){
		printf("***ERROR*** openning %s file\n", file);
	}
	else{
		for (i = 0; i < size_x; i++){
			for (j = 0; j < size_y; j++){
				fprintf(REC,"%g , %g \n", creal(array[i][j]), cimag(array[i][j]));
			}
		}

		fclose(REC);
		printf("Save %s file...\n", file);
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  save_1D_complex_double_array
 *  Description:  
 * =====================================================================================
 */

void save_1D_complex_double_array(char *file, complex double *array, int size){

	FILE *REC;
	register unsigned int i;

	REC = fopen(file,"w");
	if(REC == NULL){
		printf("***ERROR*** openning %s file\n", file);
	}
	else{
		for (i = 0; i < size; i++){
			fprintf(REC,"%g , %g \n", creal(array[i]), cimag(array[i]));
		}

		fclose(REC);
		printf("Save %s file...\n", file);
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  create_prob
 *  Description:  
 * =====================================================================================
 */

voltage_prob_t create_voltage_prob(double min_x, double min_y, double min_z, double max_x, double max_y, double max_z, int direction){

	voltage_prob_t vp;

	vp.min_x = min_x;
	vp.min_y = min_y;
	vp.min_z = min_z;
	vp.max_x = max_x;
	vp.max_y = max_y;
	vp.max_z = max_z;
	vp.direction = direction;

	return vp;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  create_voltage_source
 *  Description:  
 * =====================================================================================
 */

voltage_source_t create_voltage_source(double min_x, double min_y, double min_z, double max_x, double max_y, double max_z, int direction, double resistance, double magnitude, int waveform_type, int waveform_index){

	voltage_source_t v;

	v.min_x = min_x;
	v.min_y = min_y;
	v.min_z = min_z;
	v.max_x = max_x;
	v.max_y = max_y;
	v.max_z = max_z; 
    v.direction = direction;
    v.resistance = resistance;
    v.magnitude = magnitude;
    v.waveform_type = waveform_type;
    v.waveform_index = waveform_index;

	return v;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  create_brick
 *  Description:  
 * =====================================================================================
 */

brick_t create_brick(double min_x, double min_y, double min_z, double max_x, double max_y, double max_z, int material_type){

	brick_t brick;

	brick.min_x = min_x;
	brick.min_y = min_y;
	brick.min_z = min_z;
	brick.max_x = max_x;
	brick.max_y = max_y;
	brick.max_z = max_z;
	brick.material_type = material_type;

	return brick;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  create_material
 *  Description:  
 * =====================================================================================
 */

material_t create_material(int index, double eps_r, double mu_r, double sigma_e, double sigma_m){

	material_t m;

	m.eps_r		= eps_r;
	m.mu_r		= mu_r;
	m.sigma_e	= sigma_e;
	m.sigma_m	= sigma_m;
	m.index		= index;

	return m;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  define_boundary
 *  Description:  
 * =====================================================================================
 */

void define_boundary(boundary_t *b, int air_buffer, int cpml){
    b->type_xn = CPML;
    b->air_buffer_number_of_cells_xn = air_buffer;
    b->cpml_number_of_cells_xn = cpml;

    b->type_xp = CPML;
    b->air_buffer_number_of_cells_xp = air_buffer;
    b->cpml_number_of_cells_xp = cpml;

    b->type_yn = CPML;
    b->air_buffer_number_of_cells_yn = air_buffer;
    b->cpml_number_of_cells_yn = cpml;

    b->type_yp = CPML;
    b->air_buffer_number_of_cells_yp = air_buffer;
    b->cpml_number_of_cells_yp = cpml;

    b->type_zn = CPML;
    b->air_buffer_number_of_cells_zn = air_buffer;
    b->cpml_number_of_cells_zn = cpml;

    b->type_zp = CPML;
    b->air_buffer_number_of_cells_zp = air_buffer;
    b->cpml_number_of_cells_zp = cpml;

    b->cpml_order = 3; 
    b->cpml_sigma_factor = 1.3;
    b->cpml_kappa_max = 7.0;
    b->cpml_alpha_min = 0.0;
    b->cpml_alpha_max = 0.05;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  time_to_frequency_domain
 *  Description:  
 * =====================================================================================
 */

double complex* time_to_frequency_domain(double *x, double dt, double *frequency_array,double time_shift, int number_of_time_steps,int number_of_frequencies){
	double complex *array;
	double complex j = 0.0 + 1.0*I;
	array = (double complex*)malloc(number_of_frequencies*sizeof(double complex));
	out_of_memory_status(array);

	int t_ind, f_ind;

	for (f_ind = 0; f_ind < number_of_frequencies; f_ind++){
		array[f_ind] = 0.0 + 0.0*I;
		for (t_ind = 0; t_ind < number_of_time_steps; t_ind++){
			array[f_ind] = array[f_ind] + x[t_ind]*exp(-j*(2*M_PI*frequency_array[f_ind])*(t_ind * dt + time_shift));
		}
		array[f_ind] = array[f_ind]*dt;
	}
	return array; 
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  exp_complex
 *  Description:  
 * =====================================================================================
 */

double complex exp_complex(double complex number){
	double r = creal(number);
	double i = cimag(number);
	return (double complex)((cos(i) + I*sin(i))*exp(r));
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  capture_sampled_currents
 *  Description:  
 * =====================================================================================
 */

void capture_sampled_currents(int number_of_samples, int time_step, current_prob_t *samples, double ***H_field_x, double ***H_field_y, double ***H_field_z, cell_t cell){
	double sampled_value = 0.0;
	double sampled_value1 = 0.0;
	double sampled_value2 = 0.0;
	double sampled_value3 = 0.0;
	double sampled_value4 = 0.0;

	double dx = cell.dx;
	double dy = cell.dy;
	double dz = cell.dz;

	register unsigned int i,j,k;
	int ind;	
	for (ind=1; ind<number_of_samples; ind++){
        int is = samples[ind].is;
        int js = samples[ind].js;
        int ks = samples[ind].ks;
        int ie = samples[ind].ie;
        int je = samples[ind].je;
        int ke = samples[ind].ke;

		if (samples[ind].direction == XN || samples[ind].direction == XP){
			for (j=js;j<=je;j++){
				sampled_value1 = sampled_value1 + dy*H_field_y[ie-1][j][ks-1];
				sampled_value3 = sampled_value3 + dy*H_field_y[ie-1][j][ke];
			}
			for (k=ks; k<=ke; k++){
				sampled_value2 = sampled_value2 + dz*H_field_z[ie-1][je][k];
				sampled_value4 = sampled_value4 + dz*H_field_z[ie-1][js-1][k];
			}
		}

		else if (samples[ind].direction == YN || samples[ind].direction == YP){
			for (i=is;i<=ie;i++){
				sampled_value1 = sampled_value1 + dx*H_field_x[i][je-1][ke];
				sampled_value3 = sampled_value3 + dx*H_field_x[i][je-1][ks-1];
			}
			for (k=ks; k<=ke; k++){
				sampled_value2 = sampled_value2 + dz*H_field_z[is-1][je-1][k];
				sampled_value4 = sampled_value4 + dz*H_field_z[ie][je-1][k];
			}
		}

		else if (samples[ind].direction == ZN || samples[ind].direction == ZP){
			for (i=is;i<=ie;i++){
				sampled_value1 = sampled_value1 + dx*H_field_x[i][js-1][ke-1];
				sampled_value3 = sampled_value3 + dx*H_field_x[i][je][ke-1];
			}
			for (k=ks; k<=ke; k++){
				sampled_value2 = sampled_value2 + dz*H_field_z[ie][j][ke-1];
				sampled_value4 = sampled_value4 + dz*H_field_z[is-1][j][ke-1];
			}
		}

		sampled_value = sampled_value + sampled_value1 + sampled_value2 - sampled_value3 - sampled_value4;
		
		if (samples[ind].direction == XN || samples[ind].direction == YN || samples[ind].direction == ZP){
			sampled_value = -1 * sampled_value;
		}
		samples[ind].sampled_value[time_step] = sampled_value;
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  capture_electric_fields
 *  Description:  
 * =====================================================================================
 */

void capture_electric_fields(int number_of_sampled_electric_fields, int time_step ,sampled_fields_t *sampled_electric_fields, double ***E_field_x, double ***E_field_y, double ***E_field_z){
	register unsigned int ind,i,j,k;
	for (ind=0; ind < number_of_sampled_electric_fields; ind++){
		int is = sampled_electric_fields[ind].is;
		int js = sampled_electric_fields[ind].js;
		int ks = sampled_electric_fields[ind].ks;

		double sampled_value = 0.0;
		double sampled_value_x = 0.0;
		double sampled_value_y = 0.0;
		double sampled_value_z = 0.0;
		
		if (sampled_electric_fields[ind].component == X){
			for (i=is-1;i<=is;i++){
				sampled_value = sampled_value + 0.5*E_field_x[i][js][ks];
			}
		}	
		if (sampled_electric_fields[ind].component == Y){
			for (j=js-1;j<=js;j++){
				sampled_value = sampled_value + 0.5*E_field_y[is][j][ks];
			}
		}	
		if (sampled_electric_fields[ind].component == Z){
			for (k=ks-1;k<=ks;k++){
				sampled_value = sampled_value + 0.5 * E_field_z[is][js][k];
			}
		}
		if (sampled_electric_fields[ind].component == M){

			for (i=is-1;i<=is;i++){
				sampled_value_x = sampled_value_x + 0.5*E_field_x[i][js][ks];
			}
			for (j=js-1;j<=js;j++){
				sampled_value_y = sampled_value_y + 0.5*E_field_y[is][j][ks];
			}
			for (k=ks-1;k<=ks;k++){
				sampled_value_z = sampled_value_z + 0.5 * E_field_z[is][js][k];
			}
			
            sampled_value = sqrt(sampled_value_x*sampled_value_x + sampled_value_y*sampled_value_y + sampled_value_z*sampled_value_z);
		}
		sampled_electric_fields[ind].sampled_value[time_step] = sampled_value;
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  capture_sampled_voltages
 *  Description:  
 * =====================================================================================
 */

void capture_sampled_voltages(int number_of_sampled_voltages, int time_step ,voltage_prob_t *sampled_voltages, double ***E_field_x, double ***E_field_y, double ***E_field_z){
	int ind;
	for (ind=0; ind < number_of_sampled_voltages;ind++){
		int is   = sampled_voltages[ind].is;
		int js   = sampled_voltages[ind].js;
		int ks   = sampled_voltages[ind].ks;
		int ie   = sampled_voltages[ind].ie;
		int je   = sampled_voltages[ind].je;
		int ke   = sampled_voltages[ind].ke;
		double Csvf = sampled_voltages[ind].Csvf;
		double sampled_value = 0.0;
		register unsigned int i,j,k;

		if (sampled_voltages[ind].direction == XN || sampled_voltages[ind].direction == XP){
			for (i=is;i<=ie;i++){
				for (j=js;j<=je;j++){
					for (k=ks;k<ke;k++){
						sampled_value = sampled_value + Csvf*E_field_x[i][j][k];
					}
				}
			}
		}
		if (sampled_voltages[ind].direction == YN || sampled_voltages[ind].direction == YP){
			for (i=is;i<=ie;i++){
				for (j=js;j<=je;j++){
					for (k=ks;k<ke;k++){
						sampled_value = sampled_value + Csvf*E_field_y[i][j][k];
					}
				}
			}
		}
		if (sampled_voltages[ind].direction == ZN || sampled_voltages[ind].direction == ZP){
			for (i=is;i<=ie;i++){
				for (j=js;j<=je;j++){
					for (k=ks;k<ke;k++){
						sampled_value = sampled_value + Csvf*E_field_z[i][j][k];
					}
				}
			}
		}

		sampled_voltages[ind].sampled_value[time_step] = sampled_value;
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  build_3D_int_array
 *  Description:  
 * =====================================================================================
 */

int*** build_3D_int_array(int n_x, int n_y, int n_z,int init_val)
{
	int*** array;
	array = (int***)malloc(n_x*sizeof(int**));

	out_of_memory_status(array);
	register unsigned int i_x,i_y,i_z;
	for(i_x = 0; i_x < n_x; i_x++)
	{
		array[i_x] = (int**)malloc(n_y * sizeof(int*));
		out_of_memory_status(array[i_x]);
		for(i_y = 0; i_y < n_y; i_y++)
		{
			array[i_x][i_y] = (int*)malloc(n_z * sizeof(int));
			out_of_memory_status(array[i_x][i_y]);
			for (i_z = 0; i_z < n_z; i_z++){
				array[i_x][i_y][i_z] = init_val;
			}
		}
	}

	return array;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  build_1D_double_array
 *  Description:  
 * =====================================================================================
 */

double* build_1D_double_array(int size, double init_val)
{
	double* array;
	array = (double*)malloc(size*sizeof(double));
	out_of_memory_status(array);

	return array;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  build_3D_double_array
 *  Description:  
 * =====================================================================================
 */

double*** build_3D_double_array(int n_x, int n_y, int n_z, double init_val)
{
	double*** array;
	array = (double***)malloc(n_x*sizeof(double**));

	out_of_memory_status(array);
	register unsigned int i_x,i_y,i_z;
	for(i_x = 0; i_x < n_x; i_x++)
	{
		array[i_x] = (double**)malloc(n_y * sizeof(double*));
		out_of_memory_status(array[i_x]);
		for(i_y = 0; i_y < n_y; i_y++)
		{
			array[i_x][i_y] = (double*)malloc(n_z * sizeof(double));
			out_of_memory_status(array[i_x][i_y]);
			for (i_z = 0; i_z < n_z; i_z++){
				array[i_x][i_y][i_z] = init_val;
			}

		}
	}

	return array;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  build_4D_double_array
 *  Description:  
 * =====================================================================================
 */

double**** build_4D_double_array(int n_x, int n_y, int n_z, int n_u, double init_val)
{
	double**** array;
	array = (double****)malloc(n_x*sizeof(double***));

	out_of_memory_status(array);
	register unsigned int i_x,i_y,i_z,i_u;
	for(i_x = 0; i_x < n_x; i_x++)
	{
		array[i_x] = (double***)malloc(n_y * sizeof(double**));
		out_of_memory_status(array[i_x]);
		for(i_y = 0; i_y < n_y; i_y++)
		{
			array[i_x][i_y] = (double**)malloc(n_z * sizeof(double*));
			out_of_memory_status(array[i_x][i_y]);
			for(i_z = 0; i_z < n_z; i_z++)
			{
				array[i_x][i_y][i_z] = (double*)malloc(n_u * sizeof(double));
				out_of_memory_status(array[i_x][i_y][i_z]);
				for (i_u = 0; i_u < n_u; i_u++){
					array[i_x][i_y][i_z][i_u] = init_val;
				}

			}
		}
	}

	return array;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  build_1D_double_complex_array
 *  Description:  
 * =====================================================================================
 */

double complex* build_1D_double_complex_array(int size, double complex init_val)
{
	double complex* array;
	array = (double complex*)malloc(size*sizeof(double complex));
	out_of_memory_status(array);

	return array;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  build_2D_double_complex_array
 *  Description:  
 * =====================================================================================
 */

double complex** build_2D_double_complex_array(int n_x, int n_y, double complex init_val)
{
	double complex** array;
	array = (double complex**)malloc(n_x*sizeof(double complex*));
	out_of_memory_status(array);

	register unsigned int i_x,i_y;
	for(i_x = 0; i_x < n_x; i_x++)
	{
		array[i_x] = (double complex*)malloc(n_y * sizeof(double complex));
		out_of_memory_status(array[i_x]);
		for(i_y = 0; i_y < n_y; i_y++)
		{
			array[i_x][i_y] = init_val;

		}
	}

	return array;
}
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  build_4D_double_complex_array
 *  Description:  
 * =====================================================================================
 */

double complex**** build_4D_double_complex_array(int n_x, int n_y, int n_z, int n_u, double complex init_val)
{
	double complex**** array;
	array = (double complex****)malloc(n_x*sizeof(double complex***));

	out_of_memory_status(array);
	register unsigned int i_x,i_y,i_z,i_u;
	for(i_x = 0; i_x < n_x; i_x++)
	{
		array[i_x] = (double complex***)malloc(n_y * sizeof(double complex**));
		out_of_memory_status(array[i_x]);
		for(i_y = 0; i_y < n_y; i_y++)
		{
			array[i_x][i_y] = (double complex**)malloc(n_z * sizeof(double complex*));
			out_of_memory_status(array[i_x][i_y]);
			for(i_z = 0; i_z < n_z; i_z++)
			{
				array[i_x][i_y][i_z] = (double complex*)malloc(n_u * sizeof(double complex));
				out_of_memory_status(array[i_x][i_y][i_z]);
				for (i_u = 0; i_u < n_u; i_u++){
					array[i_x][i_y][i_z][i_u] = init_val;
				}

			}
		}
	}

	return array;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  fill_two_cell_average_mu_x
 *  Description:  
 * =====================================================================================
 */

void fill_two_cell_average_mu_x(double*** array, material_t *mat, int*** index, int x_start, int x_end, int y_start, int y_end, int z_start, int z_end)
{
	register unsigned int i,j,k;
	for(i=x_start;i<x_end;i++){
		for(j=y_start;j<y_end;j++){
			for(k=z_start;k<z_end;k++){
				array[i][j][k] = (2.0*(mat[index[i][j][k]].mu_r)*(mat[index[i-1][j][k]].mu_r))/(mat[index[i][j][k]].mu_r + mat[index[i-1][j][k]].mu_r);
			}
		}
	}
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  fill_two_cell_average_mu_y
 *  Description:  
 * =====================================================================================
 */

void fill_two_cell_average_mu_y(double*** array, material_t *mat, int*** index, int x_start, int x_end, int y_start, int y_end, int z_start, int z_end)
{
	register unsigned int i,j,k;
	for(i=x_start;i<x_end;i++){
		for(j=y_start;j<y_end;j++){
			for(k=z_start;k<z_end;k++){
				array[i][j][k] = (2.0*(mat[index[i][j][k]].mu_r)*(mat[index[i][j-1][k]].mu_r))/(mat[index[i][j][k]].mu_r + mat[index[i][j-1][k]].mu_r);
			}
		}
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  fill_two_cell_average_mu_z
 *  Description:  
 * =====================================================================================
 */

void fill_two_cell_average_mu_z(double*** array, material_t *mat, int*** index, int x_start, int x_end, int y_start, int y_end, int z_start, int z_end)
{
	register unsigned int i,j,k;
	for(i=x_start;i<x_end;i++){
		for(j=y_start;j<y_end;j++){
			for(k=z_start;k<z_end;k++){
				array[i][j][k] = (2.0*(mat[index[i][j][k]].mu_r)*(mat[index[i][j][k-1]].mu_r))/(mat[index[i][j][k]].mu_r + mat[index[i][j][k-1]].mu_r);
			}
		}
	}
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  fill_two_cell_average_sigma_x
 *  Description:  
 * =====================================================================================
 */

void fill_two_cell_average_sigma_x(double*** array, material_t *mat, int*** index, int x_start, int x_end, int y_start, int y_end, int z_start, int z_end)
{
	register unsigned int i,j,k;
	for(i=x_start;i<x_end;i++){
		for(j=y_start;j<y_end;j++){
			for(k=z_start;k<z_end;k++){
				array[i][j][k] = (2.0*(mat[index[i][j][k]].sigma_m)*(mat[index[i-1][j][k]].sigma_m))/(mat[index[i][j][k]].sigma_m + mat[index[i-1][j][k]].sigma_m);
			}
		}
	}
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  fill_two_cell_average_sigma_y
 *  Description:  
 * =====================================================================================
 */

void fill_two_cell_average_sigma_y(double*** array, material_t *mat, int*** index, int x_start, int x_end, int y_start, int y_end, int z_start, int z_end)
{
	register unsigned int i,j,k;
	for(i=x_start;i<x_end;i++){
		for(j=y_start;j<y_end;j++){
			for(k=z_start;k<z_end;k++){
				array[i][j][k] = (2.0*(mat[index[i][j][k]].sigma_m)*(mat[index[i][j-1][k]].sigma_m))/(mat[index[i][j][k]].sigma_m + mat[index[i][j-1][k]].sigma_m);
			}
		}
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  fill_two_cell_average_sigma_z
 *  Description:  
 * =====================================================================================
 */

void fill_two_cell_average_sigma_z(double*** array, material_t *mat, int*** index, int x_start, int x_end, int y_start, int y_end, int z_start, int z_end)
{
	register unsigned int i,j,k;
	for(i=x_start;i<x_end;i++){
		for(j=y_start;j<y_end;j++){
			for(k=z_start;k<z_end;k++){
				array[i][j][k] = (2.0*(mat[index[i][j][k]].sigma_m)*(mat[index[i][j][k-1]].sigma_m))/(mat[index[i][j][k]].sigma_m + mat[index[i][j][k-1]].sigma_m);
			}
		}
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  fill_four_cell_average_eps_x
 *  Description:  
 * =====================================================================================
 */
 
void fill_four_cell_average_eps_x(double*** array,material_t *mat,int*** index, int x_start, int x_end, int y_start, int y_end, int z_start, int z_end)
{
	register unsigned int i,j,k;
	for(i=x_start;i<x_end;i++){
		for(j=y_start;j<y_end;j++){
			for(k=z_start;k<z_end;k++){
				array[i][j][k] = 0.25*(mat[index[i][j][k]].eps_r + mat[index[i][j-1][k]].eps_r + mat[index[i][j][k-1]].eps_r + mat[index[i][j-1][k-1]].eps_r);
			}
		}
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  fill_four_cell_average_eps_y
 *  Description:  
 * =====================================================================================
 */
 
void fill_four_cell_average_eps_y(double*** array,material_t *mat,int*** index, int x_start, int x_end, int y_start, int y_end, int z_start, int z_end)
{
	register unsigned int i,j,k;
	for(i=x_start;i<x_end;i++){
		for(j=y_start;j<y_end;j++){
			for(k=z_start;k<z_end;k++){
				array[i][j][k] = 0.25*(mat[index[i][j][k]].eps_r + mat[index[i-1][j][k]].eps_r + mat[index[i][j][k-1]].eps_r + mat[index[i-1][j][k-1]].eps_r);
			}
		}
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  fill_four_cell_average_eps_z
 *  Description:  
 * =====================================================================================
 */
 
void fill_four_cell_average_eps_z(double*** array,material_t *mat,int*** index, int x_start, int x_end, int y_start, int y_end, int z_start, int z_end)
{
	register unsigned int i,j,k;
	for(i=x_start;i<x_end;i++){
		for(j=y_start;j<y_end;j++){
			for(k=z_start;k<z_end;k++){
				array[i][j][k] = 0.25*(mat[index[i][j][k]].eps_r + mat[index[i-1][j][k]].eps_r + mat[index[i][j-1][k]].eps_r + mat[index[i-1][j-1][k]].eps_r);
			}
		}
	}
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  fill_four_cell_average_sigma_x
 *  Description:  
 * =====================================================================================
 */
 
void fill_four_cell_average_sigma_x(double*** array,material_t *mat,int*** index, int x_start, int x_end, int y_start, int y_end, int z_start, int z_end)
{
	register unsigned int i,j,k;
	for(i=x_start;i<x_end;i++){
		for(j=y_start;j<y_end;j++){
			for(k=z_start;k<z_end;k++){
				array[i][j][k] = 0.25*(mat[index[i][j][k]].sigma_e + mat[index[i][j-1][k]].sigma_e + mat[index[i][j][k-1]].sigma_e + mat[index[i][j-1][k-1]].sigma_e);
			}
		}
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  fill_four_cell_average_sigma_y
 *  Description:  
 * =====================================================================================
 */
 
void fill_four_cell_average_sigma_y(double*** array,material_t *mat,int*** index, int x_start, int x_end, int y_start, int y_end, int z_start, int z_end)
{
	register unsigned int i,j,k;
	for(i=x_start;i<x_end;i++){
		for(j=y_start;j<y_end;j++){
			for(k=z_start;k<z_end;k++){
				array[i][j][k] = 0.25*(mat[index[i][j][k]].sigma_e + mat[index[i-1][j][k]].sigma_e + mat[index[i][j][k-1]].sigma_e + mat[index[i-1][j][k-1]].sigma_e);
			}
		}
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  fill_four_cell_average_sigma_z
 *  Description:  
 * =====================================================================================
 */
 
void fill_four_cell_average_sigma_z(double*** array,material_t *mat,int*** index, int x_start, int x_end, int y_start, int y_end, int z_start, int z_end)
{
	register unsigned int i,j,k;
	for(i=x_start;i<x_end;i++){
		for(j=y_start;j<y_end;j++){
			for(k=z_start;k<z_end;k++){
				array[i][j][k] = 0.25*(mat[index[i][j][k]].sigma_e + mat[index[i-1][j][k]].sigma_e + mat[index[i][j-1][k]].sigma_e + mat[index[i-1][j-1][k]].sigma_e);
			}
		}
	}
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  out_of_memory_status
 *  Description:  
 * =====================================================================================
 */

void out_of_memory_status(void* array){
	if(array == NULL)
	{
		fprintf(stderr, "out of memory\n");
		exit(EXIT_FAILURE);
	}
}

