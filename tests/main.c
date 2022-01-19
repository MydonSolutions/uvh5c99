#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "uvh5.h"
#include "uvh5/uvh5_toml.h"

void uvh5_toml_parse_telescope_info(UVH5_header_t* uvh5_header, char* file_path) {
	FILE* fp;
	char errbuf[200];

	fp = fopen(file_path, "r");
	if (!fp) {
			uvh5_print_error(__FUNCTION__, "cannot open '%s' (%s)", file_path, strerror(errno));
			return;
	}

	toml_table_t* conf = toml_parse_file(fp, errbuf, sizeof(errbuf));
	fclose(fp);
	if (!conf) {
			uvh5_print_error(__FUNCTION__, "cannot parse %s", errbuf);
			return;
	}
	uvh5_toml_string_in(conf, "telescope_name", &uvh5_header->telescope_name);
	uvh5_toml_sexagesimal_in(conf, "latitude", &uvh5_header->latitude);
	uvh5_toml_sexagesimal_in(conf, "longitude", &uvh5_header->longitude);
	uvh5_toml_double_in(conf, "altitude", &uvh5_header->altitude);

	float universal_diameter = -1.0;
	uvh5_toml_float_in(conf, "antenna_diameter", &universal_diameter);
	
	toml_array_t* toml_antennas_array = toml_array_in(conf, "antennas");
	if (!toml_antennas_array) {
			uvh5_print_error(__FUNCTION__, "missing [[antennas]]");
	}
	else {
		uvh5_header->Nants_telescope = toml_array_nelem(toml_antennas_array);
		UVH5Halloc(uvh5_header); // Alloc Nants_telescope related
		uvh5_header->antenna_diameters = malloc(sizeof(float) * uvh5_header->Nants_telescope); // Optional

		for (size_t i = 0; i < uvh5_header->Nants_telescope; i++)
		{
			toml_table_t* toml_antenna_info = toml_table_at(toml_antennas_array, i);
			if(toml_antenna_info) {
				uvh5_print_verbose(__FUNCTION__, "Antenna: %ld", i);
				uvh5_toml_antenna_table_in(
					toml_antenna_info,
					uvh5_header->antenna_numbers + i,
					uvh5_header->antenna_names + i,
					uvh5_header->antenna_positions + 3*i,
					universal_diameter > 0.0f ? NULL : 
						uvh5_header->antenna_diameters != NULL ? uvh5_header->antenna_diameters + i : NULL
				);
				if(uvh5_header->antenna_diameters != NULL && universal_diameter > 0.0f) {
					uvh5_header->antenna_diameters[i] = universal_diameter;
				}
			}
			else {
				uvh5_print_error(__FUNCTION__, "cannot access antenna #%ld info", i);
			}
		}

		char *ant_pos_frame_str = NULL;
		char ant_pos_frame = FRAME_ECEF;
		if(uvh5_toml_string_in(conf, "antenna_position_frame", &ant_pos_frame_str)) {
			// Not specified
			double dist = uvh5_calc_hypotenuse(uvh5_header->antenna_positions, 3);
			if(dist < 6e6) {
				ant_pos_frame = FRAME_ENU;
			}
			else {
				ant_pos_frame = FRAME_ECEF;
			}
		}
		else {
			if(strcmp(ant_pos_frame_str, "ecef") == 0) {
				ant_pos_frame = FRAME_ECEF;
			}
			else if(strcmp(ant_pos_frame_str, "enu") == 0) {
				ant_pos_frame = FRAME_ENU;
			}
			else {
				uvh5_print_warn(__FUNCTION__, "Ignoring 'antenna_position_frame' specfication: '%s'.", ant_pos_frame_str);
			}
		}

		switch(ant_pos_frame) {
		 case FRAME_ECEF:
				uvh5_print_info(__FUNCTION__, "Translating from ECEF to XYZ!");
				uvh5_calc_position_to_xyz_frame_from_ecef(
					uvh5_header->antenna_positions,
					uvh5_header->Nants_telescope,
					uvh5_calc_deg2rad(uvh5_header->longitude),
					uvh5_calc_deg2rad(uvh5_header->latitude),
					uvh5_header->altitude
				);
				break;
			case FRAME_ENU:
				uvh5_print_info(__FUNCTION__, "Translating from ENU to XYZ!");
				uvh5_calc_position_to_xyz_frame_from_enu(
					uvh5_header->antenna_positions,
					uvh5_header->Nants_telescope,
					uvh5_calc_deg2rad(uvh5_header->longitude),
					uvh5_calc_deg2rad(uvh5_header->latitude),
					uvh5_header->altitude
				);
				break;
			case FRAME_XYZ:
				uvh5_print_info(__FUNCTION__, "Verbatim XYZ positions.");
				break;
		}
	}

	toml_free(conf);
}

// Populate baselines(ant_1/2_array)
void uvh5_toml_parse_obs_info(UVH5_header_t* uvh5_header, char* file_path) {
	FILE* fp;
	char errbuf[200];

	fp = fopen(file_path, "r");
	if (!fp) {
		uvh5_print_error(__FUNCTION__, "cannot open '%s' (%s)", file_path, strerror(errno));
		return;
	}

	toml_table_t* conf = toml_parse_file(fp, errbuf, sizeof(errbuf));
	fclose(fp);
	if (!conf) {
		uvh5_print_error(__FUNCTION__, "cannot parse %s", errbuf);
		return;
	}
	
	toml_array_t* toml_input_mapping = toml_array_in(conf, "input_map");
	if(toml_input_mapping) {
		int Npol_ant = 0; // Initial assumption
		uvh5_print_verbose(__FUNCTION__, "input_map length: %d", toml_array_nelem(toml_input_mapping));

		char *ant_name0 = NULL, *ant_name1 = NULL;
		char pols_ant[4] = {'\0'}; // at most 2 unique pols
		do {
			if (ant_name1 == NULL) {
				free(ant_name1);
			}
			toml_array_t* toml_input_ant_name_pol = toml_array_at(toml_input_mapping, Npol_ant);
			uvh5_toml_string_at(toml_input_ant_name_pol, 0, &ant_name1);
			uvh5_toml_nstring_at(toml_input_ant_name_pol, 1, pols_ant + Npol_ant, 1);
			if ( ant_name0 == NULL) {
				ant_name0 = malloc(strlen(ant_name1));
				memcpy(ant_name0, ant_name1, strlen(ant_name1)+1);
			}
		} while(strcmp(ant_name0, ant_name1) == 0 && ++Npol_ant <= 2);
		if (Npol_ant > 2) {
			uvh5_print_error(__FUNCTION__, "The first antenna is repeated more than twice indicating more than 2 polarisaitons.");
			return;
		}

		uvh5_print_info(__FUNCTION__, "Assuming every antenna has %d polarisations, as the first does.", Npol_ant);
		int Nants_data = toml_array_nelem(toml_input_mapping)/Npol_ant;
		uvh5_print_info(__FUNCTION__, "\tLeads to Nants_data: %d.", Nants_data);
		
		int* antenna_data_numbers = malloc(Nants_data*sizeof(int));
		antenna_data_numbers[0] = uvh5_header->antenna_numbers[find_antenna_index_by_name(uvh5_header, ant_name0)];
		antenna_data_numbers[1] = uvh5_header->antenna_numbers[find_antenna_index_by_name(uvh5_header, ant_name1)];

		free(ant_name1);
		for (size_t i = 2; i < Nants_data; i++) {
			toml_array_t* toml_input_ant_name_pol = toml_array_at(toml_input_mapping, i*Npol_ant);
			uvh5_toml_string_at(toml_input_ant_name_pol, 0, &ant_name1);
			antenna_data_numbers[i] = uvh5_header->antenna_numbers[find_antenna_index_by_name(uvh5_header, ant_name1)];
			free(ant_name1);
		}
		
		uvh5_header->Npols = Npol_ant*Npol_ant; // header->Npols is the pol-products
		uvh5_header->Nants_data = Nants_data;
		uvh5_header->Nbls = (uvh5_header->Nants_data*(uvh5_header->Nants_data+1))/2;
		UVH5Halloc(uvh5_header);

		char pol_product[3] = {'\0'};
		for (size_t i = 0; i < Npol_ant; i++) {
			pol_product[0] = pols_ant[i];
			for (size_t j = 0; j < Npol_ant; j++) {
				pol_product[1] = pols_ant[j];
				uvh5_header->polarization_array[i*2+j] = polarisation_string_key(pol_product, Npol_ant);
				uvh5_print_verbose(__FUNCTION__, "Pol-product '%s' with key %d.", pol_product, uvh5_header->polarization_array[i*2+j]);
			}
		}

		for (size_t i = 0; i < uvh5_header->Nants_data; i++)
		{
			uvh5_header->ant_1_array[i] = antenna_data_numbers[i];
			uvh5_header->ant_2_array[i] = antenna_data_numbers[i];
		}
		
		int ant_1_idx = 0;
		int ant_2_idx = 1;
		for (int bls_idx = uvh5_header->Nants_data; bls_idx < uvh5_header->Nbls; )
		{
			if(ant_1_idx != ant_2_idx) {
				uvh5_header->ant_1_array[bls_idx] = antenna_data_numbers[ant_1_idx];//ant_1_idx < ant_2_idx ? ant_1_idx : ant_2_idx];
				uvh5_header->ant_2_array[bls_idx] = antenna_data_numbers[ant_2_idx];//ant_1_idx > ant_2_idx ? ant_1_idx : ant_2_idx];
				bls_idx++;
			}
			ant_2_idx = (ant_2_idx + 1)%uvh5_header->Nants_data;
			if(ant_2_idx == 0) {
				ant_1_idx += 1;
			}
		}

		free(ant_name0);
		free(antenna_data_numbers);
	}
	else {
		uvh5_print_error(__FUNCTION__, "cannot read location 'input_map'");
	}

	toml_free(conf);
}

int main(int argc, const char * argv[]) {
	if (argc != 3) {
		fprintf(stderr, "Provide the telescope and observation info toml files.\n");
		return 1;
	}

	/* open a new file */
	UVH5_file_t uvh5 = {0};
	UVH5_header_t* uvh5_header = &uvh5.header;

	// set header scalar data
	uvh5_header->Ntimes = 0; // initially
	uvh5_header->Nfreqs = 192;
	uvh5_header->Nspws = 1;
	uvh5_header->Nblts = uvh5_header->Nbls * uvh5_header->Ntimes;

	uvh5_toml_parse_telescope_info(uvh5_header, (char*) argv[1]);
	UVH5Hadmin(uvh5_header);
	uvh5_toml_parse_obs_info(uvh5_header, (char*) argv[2]);

	uvh5_header->instrument = uvh5_header->telescope_name;
	uvh5_header->object_name = "zenith";
	uvh5_header->history = "None";
	uvh5_header->phase_type = "drift";

	if(strcmp(uvh5_header->phase_type, "drift") == 0) {
		memcpy(uvh5_header->_antenna_uvw_positions, uvh5_header->_antenna_enu_positions, sizeof(double)*uvh5_header->Nants_telescope*3);
		UVH5permutate_uvws(uvh5_header);
	}
	else if(strcmp(uvh5_header->phase_type, "phased") == 0) {
		memcpy(uvh5_header->_antenna_uvw_positions, uvh5_header->_antenna_enu_positions, sizeof(double)*uvh5_header->Nants_telescope*3);
		double hour_angle_rad = 0.0;
		double declination_rad = 0.0;
		uvh5_calc_position_to_uvw_frame_from_enu(
			uvh5_header->_antenna_uvw_positions,
			uvh5_header->Nants_data,
			hour_angle_rad,
			declination_rad,
			uvh5_calc_deg2rad(uvh5_header->latitude)
		);

		UVH5permutate_uvws(uvh5_header);
	}

	uvh5_header->flex_spw = UVH5_FALSE;

	float tau = 1.0;
	for (size_t i = 0; i < uvh5_header->Nbls; i++)
	{
		uvh5_header->time_array[i] = uvh5_calc_julian_date_from_guppi_param(16.0, 16*8192, 8192, 12371829, 0) + tau/(DAYSEC*2);
		uvh5_header->integration_time[i] = tau;
	}

	UVH5open("test_file.uvh5", &uvh5, UVH5TcreateCF32());

	while(uvh5_header->Ntimes < 3) {
		UVH5write_dynamic(&uvh5);
		for (size_t i = 0; i < uvh5_header->Nbls; i++)
		{
			uvh5_header->time_array[i] += tau/DAYSEC;
		}
	}

	free(uvh5_header->telescope_name);
	for (size_t i = 0; i < uvh5_header->Nants_telescope; i++)
	{
		free(uvh5_header->antenna_names[i]);
	}
	free(uvh5_header->antenna_names);
	uvh5_header->antenna_names = NULL;
	UVH5close(&uvh5);
	return 0;
}
