#ifndef UVH5_TOML_H
#define UVH5_TOML_H

#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include "toml.h"
#include "uvh5.h"
#include "uvh5_printer.h"

int UVH5toml_sexagesimal_in(toml_table_t* parent, const char* location, double* sexagesimal);
int UVH5toml_string_in(toml_table_t* parent, const char* location, char** string_out);
int UVH5toml_nstring_in(toml_table_t* parent, const char* location, char* string_out, size_t length);
int UVH5toml_double_in(toml_table_t* parent, const char* location, double* double_out);
int UVH5toml_float_in(toml_table_t* parent, const char* location, float* float_out);
int UVH5toml_int_in(toml_table_t* parent, const char* location, int* int_out);

int UVH5toml_string_at(toml_array_t* parent, const int location, char** string_out);
int UVH5toml_nstring_at(toml_array_t* parent, const int location, char* string_out, size_t length);
int UVH5toml_double_at(toml_array_t* parent, const int location, double* double_out);
int UVH5toml_float_at(toml_array_t* parent, const int location, float* float_out);
int UVH5toml_int_at(toml_array_t* parent, const int location, int* int_out);

int UVH5toml_antenna_table_in(
	toml_table_t* parent,
	int* ant_id,
	char** ant_name,
	double* ant_pos,
	float* ant_diameter
);

void UVH5toml_parse_telescope_info(char* file_path, UVH5_header_t* UVH5header);
void UVH5toml_parse_observation_info(char* file_path, UVH5_header_t* UVH5header);

#endif