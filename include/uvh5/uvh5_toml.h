#ifndef UVH5_TOML_H
#define UVH5_TOML_H

#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include "toml.h"
#include "uvh5_printer.h"

int uvh5_toml_sexagesimal_in(toml_table_t* parent, const char* location, double* sexagesimal);
int uvh5_toml_string_in(toml_table_t* parent, const char* location, char** string_out);
int uvh5_toml_nstring_in(toml_table_t* parent, const char* location, char* string_out, size_t length);
int uvh5_toml_double_in(toml_table_t* parent, const char* location, double* double_out);
int uvh5_toml_float_in(toml_table_t* parent, const char* location, float* float_out);
int uvh5_toml_int_in(toml_table_t* parent, const char* location, int* int_out);

int uvh5_toml_string_at(toml_array_t* parent, const int location, char** string_out);
int uvh5_toml_nstring_at(toml_array_t* parent, const int location, char* string_out, size_t length);
int uvh5_toml_double_at(toml_array_t* parent, const int location, double* double_out);
int uvh5_toml_float_at(toml_array_t* parent, const int location, float* float_out);
int uvh5_toml_int_at(toml_array_t* parent, const int location, int* int_out);

int uvh5_toml_antenna_table_in(
	toml_table_t* parent,
	int* ant_id,
	char** ant_name,
	double* ant_pos,
	float* ant_diameter
);

#endif