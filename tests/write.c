#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "uvh5.h"
#include "uvh5/uvh5_toml.h"

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
	uvh5_header->object_name = "zenith";
	uvh5_header->history = "None";
	uvh5_header->phase_type = "drift";

	if(strcmp(uvh5_header->phase_type, "drift") == 0) {
		memcpy(uvh5_header->_antenna_uvw_positions, uvh5_header->_antenna_enu_positions, sizeof(double)*uvh5_header->Nants_telescope*3);
		UVH5permute_uvws(uvh5_header);
	}
	else if(strcmp(uvh5_header->phase_type, "phased") == 0) {
		uvh5_header->phase_center_ra = 0.0;
		uvh5_header->phase_center_dec = 0.0;
		uvh5_header->phase_center_epoch = 2000.0;
		uvh5_header->phase_center_frame = "icrs";

		memcpy(uvh5_header->_antenna_uvw_positions, uvh5_header->_antenna_enu_positions, sizeof(double)*uvh5_header->Nants_telescope*3);
		double hour_angle_rad = 0.0;
		double declination_rad = 0.0;
		calc_position_to_uvw_frame_from_enu(
			uvh5_header->_antenna_uvw_positions,
			uvh5_header->Nants_data,
			hour_angle_rad,
			declination_rad,
			calc_deg2rad(uvh5_header->latitude)
		);

		UVH5permute_uvws(uvh5_header);
	}

	uvh5_header->flex_spw = UVH5_FALSE;

	float tau = 1.0;
	for (size_t i = 0; i < uvh5_header->Nbls; i++)
	{
		uvh5_header->time_array[i] = calc_julian_date_from_guppi_param(16.0, 16*8192, 8192, 12371829, 0) + tau/(DAYSEC*2);
		uvh5_header->integration_time[i] = tau;
	}

	UVH5open("test_file.uvh5", &uvh5, UVH5TcreateCF32());

	while(uvh5_header->Ntimes < 2) {
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
