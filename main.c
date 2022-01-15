#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "uvh5.h"
#include "uvh5_toml.h"

void toml_load_telescope_data(UVH5_header_t* uvh5_header, char* file_path) {
	FILE* fp;
	char errbuf[200];

	fp = fopen(file_path, "r");
	if (!fp) {
			uvh5_toml_error("cannot open telescope info - ", strerror(errno));
			return;
	}

	toml_table_t* conf = toml_parse_file(fp, errbuf, sizeof(errbuf));
	fclose(fp);
	if (!conf) {
			uvh5_toml_error("cannot parse", errbuf);
			return;
	}
	uvh5_toml_string_in(conf, "telescope_name", &uvh5_header->telescope_name);
	uvh5_toml_sexagesimal_in(conf, "latitude", &uvh5_header->latitude);
	uvh5_toml_sexagesimal_in(conf, "longitude", &uvh5_header->longitude);
	uvh5_toml_float_in(conf, "altitude", &uvh5_header->altitude);

	float universal_diameter = -1.0;
	uvh5_toml_float_in(conf, "antenna_diameter", &universal_diameter);
	
	toml_array_t* toml_antennas_array = toml_array_in(conf, "antennas");
	if (!toml_antennas_array) {
			uvh5_toml_error("missing [[antennas]]", "");
	}
	else {
		uvh5_header->Nants_telescope = toml_array_nelem(toml_antennas_array);
		UVH5Halloc(uvh5_header); // Alloc Nants_telescope related
		uvh5_header->antenna_diameters = malloc(sizeof(float) * uvh5_header->Nants_telescope); // Optional

		int highest_antenna_number = -1;

		for (size_t i = 0; i < uvh5_header->Nants_telescope; i++)
		{
			toml_table_t* toml_antenna_info = toml_table_at(toml_antennas_array, i);
			if(toml_antenna_info) {
				fprintf(stderr, "Antenna: %ld\n", i);
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
				highest_antenna_number = highest_antenna_number < uvh5_header->antenna_numbers[i] ?
					uvh5_header->antenna_numbers[i] :
				highest_antenna_number;
			}
			else {
				uvh5_toml_error("cannot access antenna info", "");
			}
		}
		
		uvh5_header->_antenna_num_idx_map = malloc(sizeof(int) * (highest_antenna_number + 1)); // Administrative
		memset(uvh5_header->_antenna_num_idx_map, -1, sizeof(int) * (highest_antenna_number + 1));
		for (size_t i = 0; i < uvh5_header->Nants_telescope; i++) {
			uvh5_header->_antenna_num_idx_map[uvh5_header->antenna_numbers[i]] = i;
			fprintf(stderr, "Antenna number %d is at index %ld\n", uvh5_header->antenna_numbers[i], i);
		}
	}

	toml_free(conf);
}

// Populate baselines(ant_1/2_array)
void toml_load_obs_info(UVH5_header_t* uvh5_header, char* file_path) {
	FILE* fp;
	char errbuf[200];

	fp = fopen(file_path, "r");
	if (!fp) {
		uvh5_toml_error("cannot open observation info", strerror(errno));
		return;
	}
	fprintf(stderr, "Opened!\n");

	toml_table_t* conf = toml_parse_file(fp, errbuf, sizeof(errbuf));
	fclose(fp);
	if (!conf) {
		uvh5_toml_error("cannot parse", errbuf);
		return;
	}
	fprintf(stderr, "Parsed!\n");
	
	toml_array_t* toml_input_mapping = toml_array_in(conf, "input_map");
	if(toml_input_mapping) {
		int Npol = 0; // Initial assumption
		fprintf(stderr, "input_map length: %d\n", toml_array_nelem(toml_input_mapping));

		char *ant_name0 = NULL, *ant_name1 = NULL;
		toml_array_t* toml_input_ant0_name_pol = toml_array_at(toml_input_mapping, Npol++);
		fprintf(stderr, "AntInfo length: %d\n", toml_array_nelem(toml_input_ant0_name_pol));
		uvh5_toml_string_at(toml_input_ant0_name_pol, 0, &ant_name0);
		do {
			if(ant_name1 == NULL) {
				free(ant_name1);
			}
			toml_array_t* toml_input_ant_name_pol = toml_array_at(toml_input_mapping, Npol++);
			uvh5_toml_string_at(toml_input_ant_name_pol, 0, &ant_name1);
		} while(strcmp(ant_name0, ant_name1) == 0);

		int Nants_data = toml_array_nelem(toml_input_mapping)/(--Npol);
		fprintf(stderr, "Assuming every antenna has %d polarisations, as the first does.\n", Npol);
		fprintf(stderr, "\tLeads to Nants_data: %d.\n", Nants_data);
		
		int* ant_name_index_map = malloc(Nants_data*sizeof(int));
		ant_name_index_map[0] = find_antenna_index_by_name(uvh5_header, ant_name0);
		ant_name_index_map[1] = find_antenna_index_by_name(uvh5_header, ant_name1);

		free(ant_name1);
		for (size_t i = 2; i < Nants_data; i++) {
			toml_array_t* toml_input_ant_name_pol = toml_array_at(toml_input_mapping, i*Npol);
			uvh5_toml_string_at(toml_input_ant_name_pol, 0, &ant_name1);
			ant_name_index_map[i] = find_antenna_index_by_name(uvh5_header, ant_name1);
			free(ant_name1);
		}
		
		uvh5_header->Nants_data = Nants_data;
		uvh5_header->Npols = Npol;
		uvh5_header->Nbls = (uvh5_header->Nants_data*(uvh5_header->Nants_data+1))/2;
		UVH5Halloc(uvh5_header);

		for (size_t i = 0; i < uvh5_header->Nants_data; i++)
		{
			uvh5_header->ant_1_array[i] = ant_name_index_map[i];
			uvh5_header->ant_2_array[i] = ant_name_index_map[i];
		}
		
		int ant_1_idx = 0;
		int ant_2_idx = 1;
		for (int bls_idx = uvh5_header->Nants_data; bls_idx < uvh5_header->Nbls; )
		{
			if(ant_1_idx != ant_2_idx) {
				uvh5_header->ant_1_array[bls_idx] = ant_name_index_map[ant_1_idx];//ant_1_idx < ant_2_idx ? ant_1_idx : ant_2_idx];
				uvh5_header->ant_2_array[bls_idx] = ant_name_index_map[ant_2_idx];//ant_1_idx > ant_2_idx ? ant_1_idx : ant_2_idx];
				// fprintf(stderr, "#%d: Cross Correlation Baseline between %s and %s.\n", bls_idx, 
				// 	uvh5_header->antenna_names[uvh5_header->ant_1_array[bls_idx]],
				// 	uvh5_header->antenna_names[uvh5_header->ant_2_array[bls_idx]]
				// );
				bls_idx++;
			}
			ant_2_idx = (ant_2_idx + 1)%uvh5_header->Nants_data;
			if(ant_2_idx == 0) {
				ant_1_idx += 1;
			}
		}

		free(ant_name0);
		free(ant_name_index_map);
	}
	else {
		uvh5_toml_error("cannot read location", "input_map");
	}

	toml_free(conf);
}

int main() {
	/* open a new file */
	UVH5_file_t uvh5 = {0};
	UVH5_header_t* uvh5_header = &uvh5.header;

	// set header scalar data
	uvh5_header->Ntimes = 0; // initially
	uvh5_header->Nfreqs = 192;
	uvh5_header->Nspws = 1;
	uvh5_header->Nblts = uvh5_header->Nbls * uvh5_header->Ntimes;

	toml_load_telescope_data(uvh5_header, "./telinfo_ata.toml");
	toml_load_obs_info(uvh5_header, "./obsinfo.toml");

	uvh5_header->instrument = uvh5_header->telescope_name;
	uvh5_header->object_name = "zenith";
	uvh5_header->history = "None";
	uvh5_header->phase_type = "drift";

	uvh5_header->flex_spw = UVH5_FALSE;

	for (size_t i = 0; i < uvh5_header->Nbls; i++)
	{
		uvh5_header->uvw_array[i*3 + 0] = 0.0;
		uvh5_header->uvw_array[i*3 + 1] = 1.0;
		uvh5_header->uvw_array[i*3 + 2] = 2.0;
		uvh5_header->time_array[i] = 1.0;
		uvh5_header->integration_time[i] = 1.0;
	}

	uvh5.file_id = H5Fcreate("test_file.uvh5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	uvh5.DS_data_visdata.Tmem_id = UVH5TcreateCF32();
	uvh5.DS_data_visdata.Tsto_id = UVH5TcreateCF32();
	UVH5open(&uvh5);

	UVH5Hfree(uvh5_header);
	H5Fclose(uvh5.file_id);
	#if 0
	hid_t file_id = H5Fcreate("test_file.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	char string[] = "Don't wanna be an American idiot!!!";
	_H5DstringWrite(file_id, "string", 0, NULL, string);

	/* define array dimensions */
	int Nblts = 2;
	int Nfreqs = 3;
	int Npols = 2;
	hsize_t dims[3] = {Nblts, Nfreqs, Npols};

	/* initialize data array with FALSE values */
	UVH5_bool_t data[Nblts][Nfreqs][Npols];
	for (int i=0; i<Nblts; i++) {
		for (int j=0; j<Nfreqs; j++) {
			for (int k=0; k<Npols; k++) {
				data[i][j][k] = (i+k)%3 ? UVH5_TRUE : UVH5_FALSE;
			}
		}
	}

	/* make dataspace and write out data */
	_H5DboolWrite(file_id, "flags", 3, dims, (void*)data);
	
	float fdata[Nblts][Nfreqs][Npols];
	for (int i=0; i<Nblts; i++) {
		for (int j=0; j<Nfreqs; j++) {
			for (int k=0; k<Npols; k++) {
				fdata[i][j][k] = i*3.141 + j*1.618 + k*0.3333;
			}
		}
	}
	_H5DfloatWrite(file_id, "floats", 3, dims, (void*)fdata);

	/* close down */
	H5Fclose(file_id);
	#endif
	return 0;
}