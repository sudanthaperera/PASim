/*
 * =====================================================================================
 *
 *       Filename:  pastypes.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  9/17/2014 2:56:02 PM
 *       Compiler:  gcc
 *
 *         Author:  Sudantha Perera, 
 *   Organization:  OU-RIL
 *
 * =====================================================================================
 */

#include <complex.h>

typedef struct{
	int type_xn;
	int air_buffer_number_of_cells_xn;
	int cpml_number_of_cells_xn;

	int type_xp;
	int air_buffer_number_of_cells_xp;
	int cpml_number_of_cells_xp;

	int type_yn;
	int air_buffer_number_of_cells_yn;
	int cpml_number_of_cells_yn;

	int type_yp;
	int air_buffer_number_of_cells_yp;
	int cpml_number_of_cells_yp;

	int type_zn;
	int air_buffer_number_of_cells_zn;
	int cpml_number_of_cells_zn;

	int type_zp;
	int air_buffer_number_of_cells_zp;
	int cpml_number_of_cells_zp;

	int cpml_order; 
	double cpml_sigma_factor;
	double cpml_kappa_max;
	double cpml_alpha_min;
	double cpml_alpha_max;
} boundary_t;

typedef struct{
    double dx;
    double dy;
    double dz;
} cell_t;

typedef struct{
	double min_x;
    double min_y;
    double min_z;
	double max_x;
	double max_y;
    double max_z;
    int direction;
    double resistance;
    double magnitude;
    int waveform_type;
    int waveform_index;
    int is;
    int js;
    int ks;
    int ie;
    int je;
	int ke;
    double resistance_per_component;
	double *voltage_per_e_field;
    double *waveform;
    double ***Cezs;
	double ***Ceys;
	double ***Cexs;
    double complex *frequency_domain_value;
    double *frequencies;
} voltage_source_t;

typedef struct{
    double min_x;
    double min_y;
    double min_z;
    double max_x;
    double max_y;
    double max_z;
	int direction;
    double complex *frequency_domain_value;
    double *frequencies;
    double *waveform;
} current_source_t;

typedef struct{
    double min_x;
    double min_y;
    double min_z;
    double max_x;
    double max_y;
    double max_z;
	int direction;
	int is;
	int js;
	int ks;
	int ie;
	int je;
	int ke;
	double *sampled_value;
	double Csvf;
	double *time;
	double complex *frequency_domain_value;
	double *frequencies;
} voltage_prob_t;

typedef struct{
    double min_x;
    double min_y;
    double min_z;
    double max_x;
    double max_y;
    double max_z;
	int direction;
	int is;
	int js;
	int ks;
	int ie;
	int je;
	int ke;
	double *sampled_value;
	double *time;
	double complex *frequency_domain_value;
	double *frequencies;
} current_prob_t;

typedef struct{
    int sampled_voltage_index;
    int sampled_current_index;
    double impedance;
    int is_source_port;
	double complex *a;
	double complex *b;
	double *frequencies;
	double complex **S; 
} port_t;

typedef struct{
    int index;
    double eps_r;
    double mu_r;
    double sigma_e;
    double sigma_m;
} material_t;

typedef struct{
    double min_x;
    double min_y;
    double min_z;
    double max_x;
    double max_y;
    double max_z;
    int material_type;
} brick_t;

typedef struct{
	int number_of_cells_per_wavelength;
    double maximum_frequency;
    double tau;
    double t_0;
    double *waveform;
} gaussian_t;

typedef struct{
	double *frequencies;
	double start;
	double step;
	double end;
	int number_of_frequencies;
} frequency_domain_t;

typedef struct{
	double x;
	double y;
	double z;
	int is;
	int js;
	int ks;
	double *sampled_value;
	int component;
	double *time;
	double complex *frequency_domain_value;
	double *frequencies;
} sampled_fields_t;

typedef struct{
	double *frequencies;
    int number_of_cells_from_outer_boundary;
} farfield_t;
