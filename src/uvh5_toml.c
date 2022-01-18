#include "uvh5_toml.h"

float sexagesimal_to_float(const char *sexagesiamal_string) {
		float f = 0.0;
		char *pend;
		char *pstart = (char *) sexagesiamal_string;
		if (*sexagesiamal_string == '-') {
			pstart += 1;
		}
		do {
			f += strtof(pstart, &pend);
			f *= 60;
			pstart = pend+1;
		} while(*(pend) == ':');
		f /= 60*60*60;
		if (*sexagesiamal_string == '-') {
			f *= -1.0;
		}
		return f;
}

int uvh5_toml_sexagesimal_in(toml_table_t* parent, const char* location, float* sexagesimal) {
	toml_datum_t toml_datum = toml_string_in(parent, location);
	if (!toml_datum.ok) {
		toml_datum = toml_double_in(parent, location);
		if (!toml_datum.ok) {
			uvh5_toml_error("cannot read location", location);
		}
		*sexagesimal = toml_datum.u.d;
		printf("%s: %f\n", location, *sexagesimal);
	}
	else {
		// process sexagesimal string to float
		*sexagesimal = sexagesimal_to_float(toml_datum.u.s);
		free(toml_datum.u.s);
		printf("%s: %f\n", location, *sexagesimal);
	}
	return 0;
}

int uvh5_toml_string_at(toml_array_t* parent, const int idx, char** string_out) {
	toml_datum_t toml_datum = toml_string_at(parent, idx);
	if (!toml_datum.ok) {
			uvh5_toml_error("cannot read index", NULL);
	}
	else {
		*string_out = malloc(strlen(toml_datum.u.s)+1);
		memcpy(*string_out, toml_datum.u.s, strlen(toml_datum.u.s)+1);
		free(toml_datum.u.s);
		printf("[%d]: %s\n", idx, *string_out);
	}
	return 0;
}

int uvh5_toml_nstring_at(toml_array_t* parent, const int location, char* string_out, size_t length) {
	toml_datum_t toml_datum = toml_string_at(parent, location);
	if (!toml_datum.ok) {
			uvh5_toml_error("cannot read index", NULL);
	}
	else {
		length = length <= strlen(toml_datum.u.s) ? length : strlen(toml_datum.u.s);
		memcpy(string_out, toml_datum.u.s, length);
		free(toml_datum.u.s);
		printf("[%d]: %s\n", location, string_out);
	}
	return 0;
}

int uvh5_toml_string_in(toml_table_t* parent, const char* location, char** string_out) {
	toml_datum_t toml_datum = toml_string_in(parent, location);
	if (!toml_datum.ok) {
			uvh5_toml_error("cannot read location", location);
	}
	else {
		*string_out = malloc(strlen(toml_datum.u.s)+1);
		memcpy(*string_out, toml_datum.u.s, strlen(toml_datum.u.s)+1);
		free(toml_datum.u.s);
		printf("%s: %s\n", location, *string_out);
	}
	return 0;
}

int uvh5_toml_nstring_in(toml_table_t* parent, const char* location, char* string_out, size_t length) {
	toml_datum_t toml_datum = toml_string_in(parent, location);
	if (!toml_datum.ok) {
			uvh5_toml_error("cannot read location", NULL);
	}
	else {
		memset(string_out, '\0', length);
		length = length <= strlen(toml_datum.u.s)+1 ? length : strlen(toml_datum.u.s)+1;
		memcpy(string_out, toml_datum.u.s, length);
		free(toml_datum.u.s);
		printf("%s: %s\n", location, string_out);
	}
	return 0;
}

int uvh5_toml_double_at(toml_array_t* parent, const int idx, double* double_out) {
	toml_datum_t toml_datum = toml_double_at(parent, idx);
	if (!toml_datum.ok) {
		uvh5_toml_error("cannot read index", NULL);
	}
	else {
		*double_out = toml_datum.u.d;
		printf("[%d]: %f\n", idx, *double_out);
	}
	return 0;
}
int uvh5_toml_float_at(toml_array_t* parent, const int idx, float* float_out) {
	double intermediate;
	if(uvh5_toml_double_at(parent, idx, &intermediate)){
		uvh5_toml_error("cannot parse double as float", NULL);
	}
	*float_out = (float) intermediate;
	return 0;
}

int uvh5_toml_double_in(toml_table_t* parent, const char* location, double* double_out) {
	toml_datum_t toml_datum = toml_double_in(parent, location);
	if (!toml_datum.ok) {
			uvh5_toml_error("cannot read location", location);
	}
	else {
		*double_out = toml_datum.u.d;
		printf("%s: %f\n", location, *double_out);
	}
	return 0;
}
int uvh5_toml_float_in(toml_table_t* parent, const char* location, float* float_out){
	double intermediate;
	if(uvh5_toml_double_in(parent, location, &intermediate)) {
		uvh5_toml_error("cannot parse double as float", NULL);
	}
	*float_out = (float) intermediate;
	return 0;
}

int uvh5_toml_int_at(toml_array_t* parent, const int idx, int* int_out) {
	toml_datum_t toml_datum = toml_int_at(parent, idx);
	if (!toml_datum.ok) {
		uvh5_toml_error("cannot read index", NULL);
	}
	else {
		*int_out = toml_datum.u.i;
		printf("[%d]: %d\n", idx, *int_out);
	}
	return 0;
}

int uvh5_toml_int_in(toml_table_t* parent, const char* location, int* int_out) {
	toml_datum_t toml_datum = toml_int_in(parent, location);
	if (!toml_datum.ok) {
		uvh5_toml_error("cannot read location", location);
	}
	else {
		*int_out = toml_datum.u.i;
		printf("%s: %d\n", location, *int_out);
	}
	return 0;
}

int uvh5_toml_antenna_table_in(
	toml_table_t* parent,
	int* ant_id,
	char** ant_name,
	float* ant_pos,
	float* ant_diameter
){
	if(ant_diameter != NULL){
		if(uvh5_toml_float_in(parent, "diameter", ant_diameter)) {
			uvh5_toml_error("cannot read diameter", NULL);
		}
	}
	if(uvh5_toml_string_in(parent, "name", ant_name)) {
		uvh5_toml_error("cannot read name", NULL);
	}
	if(uvh5_toml_int_in(parent, "number", ant_id)) {
		uvh5_toml_error("cannot read number", NULL);
	}
	toml_array_t* toml_ant_position = toml_array_in(parent, "position");
	if(toml_ant_position == NULL) {
		uvh5_toml_error("cannot read position", NULL);
	}

	int ant_pos_count = toml_array_nelem(toml_ant_position);
	for (size_t i = 0; i < ant_pos_count; i++)
	{
		if(uvh5_toml_float_at(toml_ant_position, i, ant_pos+i)) {
			uvh5_toml_error("cannot read float", NULL);
		}
	}
	return 0;
}