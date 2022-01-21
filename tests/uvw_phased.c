#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "uvh5.h"
#include "uvh5/uvh5_toml.h"
#include "uvw_phased.h"

const double _atol = 0.0001;

bool isclose(double a, double b) {return abs(a - b) < _atol;}

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

	UVH5toml_parse_telescope_info((char*) argv[1], uvh5_header);
	UVH5toml_parse_observation_info((char*) argv[2], uvh5_header);
	UVH5Hadmin(uvh5_header);

	uvh5_header->instrument = uvh5_header->telescope_name;
	uvh5_header->object_name = "test";
	uvh5_header->history = "None";
	uvh5_header->phase_type = "phased";

	bool failed = false;
	for (int i = 0; !failed && i < uvh5_header->Nants_data; i++)
	{
		for (int j = 0; !failed && j < 3; j++)
		{
			if(!isclose(uvh5_header->_antenna_enu_positions[i*3 + j], ENU[i*3 + j])) {
				failed = true;
			}
		}
		if(failed) {
			UVH5print_error(__FILE__, "enu[%d] discrepancy:\n\t{\n\t\t%f vs %f,\n\t\t%f vs %f,\n\t\t%f vs %f\n\t}",
				i,
				uvh5_header->_antenna_enu_positions[i*3 + 0], ENU[i*3 + 0],
				uvh5_header->_antenna_enu_positions[i*3 + 1], ENU[i*3 + 1],
				uvh5_header->_antenna_enu_positions[i*3 + 2], ENU[i*3 + 2]
			);
		}

		int ant_idx = uvh5_header->_antenna_num_idx_map[uvh5_header->ant_1_array[i]];
		UVH5print_info(__FILE__, "#%d Antenna[%d]: #%d, %s", i, ant_idx, uvh5_header->ant_1_array[i], uvh5_header->antenna_names[ant_idx]);
	}

	for (int i = 0; !failed && i < uvh5_header->Nbls; i++)
	{
		failed |= uvh5_header->ant_1_array[i] != ANTENNA_1_ARRAY[i];
		failed |= uvh5_header->ant_2_array[i] != ANTENNA_2_ARRAY[i];

		if(failed) {
			UVH5print_error(__FILE__, "bl[%d] discrepancy:\n\t%d -> %d (%s #%d -> %s #%d),\n\t\tvs\n\t%d -> %d (%s #%d -> %s #%d)\n",
				i,
				uvh5_header->_antenna_num_idx_map[uvh5_header->ant_1_array[i]],
				uvh5_header->_antenna_num_idx_map[uvh5_header->ant_2_array[i]],
				uvh5_header->antenna_names[uvh5_header->_antenna_num_idx_map[uvh5_header->ant_1_array[i]]],
				uvh5_header->ant_1_array[i],
				uvh5_header->antenna_names[uvh5_header->_antenna_num_idx_map[uvh5_header->ant_2_array[i]]],
				uvh5_header->ant_2_array[i],

				uvh5_header->_antenna_num_idx_map[ANTENNA_1_ARRAY[i]],
				uvh5_header->_antenna_num_idx_map[ANTENNA_2_ARRAY[i]],
				uvh5_header->antenna_names[uvh5_header->_antenna_num_idx_map[ANTENNA_1_ARRAY[i]]],
				ANTENNA_1_ARRAY[i],
				uvh5_header->antenna_names[uvh5_header->_antenna_num_idx_map[ANTENNA_2_ARRAY[i]]],
				ANTENNA_2_ARRAY[i]
			);
		}
	}

	if (!failed) {
		memcpy(uvh5_header->_antenna_uvw_positions, uvh5_header->_antenna_enu_positions, sizeof(double)*uvh5_header->Nants_telescope*3);
		double _ha_rad = 0.0;
		double _decli_rad = 0.0;
		UVH5calc_ha_dec_rad(
			ra_rad,
			dec_rad,
			UVH5calc_deg2rad(uvh5_header->longitude),
			UVH5calc_deg2rad(uvh5_header->latitude),
			uvh5_header->altitude,
			timemjd,
			dut1,
			&_ha_rad,
			&_decli_rad
		);

		if (_ha_rad != hour_angle_rad) {
			UVH5print_error(__FILE__, "_ha_rad %f != %f hour_angle_rad", _ha_rad, hour_angle_rad);
		}
		if (_decli_rad != declination_rad) {
			UVH5print_error(__FILE__, "_decli_rad %f != %f declination_rad", _decli_rad, declination_rad);
		}

		UVH5calc_position_to_uvw_frame_from_enu(
			uvh5_header->_antenna_uvw_positions,
			uvh5_header->Nants_data,
			_ha_rad,
			_decli_rad,
			UVH5calc_deg2rad(uvh5_header->latitude)
		);

		UVH5permutate_uvws(uvh5_header);

		for (int i = 0; !failed && i < uvh5_header->Nbls; i++)
		{
			for (int j = 0; !failed && j < 3; j++)
			{
				if(!isclose(uvh5_header->uvw_array[i*3 + j], UVW[i*3 + j])) {
					failed = true;
				}
			}
			if(failed) {
				UVH5print_error(__FILE__, "uvw_array[%d] discrepancy:\n\t{\n\t\t%f vs %f,\n\t\t%f vs %f,\n\t\t%f vs %f\n\t}",
					i,
					uvh5_header->uvw_array[i*3 + 0], UVW[i*3 + 0],
					uvh5_header->uvw_array[i*3 + 1], UVW[i*3 + 1],
					uvh5_header->uvw_array[i*3 + 2], UVW[i*3 + 2]
				);
			}
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
	if(failed) {
		exit(1);
	}
	return 0;
}
