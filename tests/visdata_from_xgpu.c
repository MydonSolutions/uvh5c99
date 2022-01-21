#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "uvh5.h"
#include "uvh5/uvh5_toml.h"
#include "matrix_h.h"
#include "visdata.h"

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
	uvh5_header->Nfreqs = VISDATA_DIMS[1];
	uvh5_header->Nspws = 1;
	uvh5_header->Nblts = uvh5_header->Nbls * uvh5_header->Ntimes;

	UVH5toml_parse_telescope_info((char*) argv[1], uvh5_header);
	UVH5toml_parse_observation_info((char*) argv[2], uvh5_header);
	UVH5Hadmin(uvh5_header);

	bool failed = false;
	// MATRIX_H [freq, xGPU-index, real_imag]
	// VISDATA  [Nbls, freq, Npols, real_imag]
	uvh5_header->Ntimes = 1;
	uvh5_header->Nblts = VISDATA_DIMS[2];
	if(VISDATA_DIMS[0] != uvh5_header->Nbls) {
		UVH5print_error("uvh5_header->Nbls", "...");
		failed = true;
	}

	if(VISDATA_DIMS[1] != MATRIX_H_DIMS[0]) {
		UVH5print_error("MATRIX_H_DIMS[0]", "...");
		failed = true;
	}

	if(VISDATA_DIMS[2] != uvh5_header->Npols) {
		UVH5print_error("uvh5_header->Npols", "...");
		failed = true;
	}

	if (!failed) {
		// UVH5visdata_from_xgpu_matrix_h comparison
		UVH5_CI32_t* xgpuOutput = (UVH5_CI32_t*) MATRIX_H;
		UVH5_CI32_t* visdata = (UVH5_CI32_t*) VISDATA;
		UVH5_CI32_t xgpuOutput_sample = {0};

		const int ant_pol_products = uvh5_header->Nbls*uvh5_header->Npols;
		const int visdata_bl_stride = uvh5_header->Nfreqs*uvh5_header->Npols;
		const int Npols = uvh5_header->Npols;
		const int Nfreqs = uvh5_header->Nfreqs;

		int visdata_offset; // = blidx*visdata_bl_stride + freq*Npols + pol_idx
		int* xgpu_idx = uvh5_header->_ant_pol_prod_xgpu_index;
		int* bl_idx = uvh5_header->_ant_pol_prod_bl_index;
		int* pol_idx = uvh5_header->_ant_pol_prod_pol_index;
		char* conjugate = uvh5_header->_ant_pol_prod_conj;
		char* is_auto = uvh5_header->_ant_pol_prod_auto;
		for (int approd_idx = 0; !failed && approd_idx < ant_pol_products; approd_idx++) {
			visdata_offset = (*bl_idx)*visdata_bl_stride + (*pol_idx);
			for (int freq = 0; !failed && freq < Nfreqs; freq++) {
				memcpy(&xgpuOutput_sample,
					xgpuOutput + (
						freq * MATRIX_H_DIMS[1] + (*xgpu_idx)),
					sizeof(UVH5_CF32_t)
				);
				if(*conjugate) {
					xgpuOutput_sample.i = -xgpuOutput_sample.i;
				}

				if(visdata[visdata_offset].r != xgpuOutput_sample.r ||
					 visdata[visdata_offset].i != xgpuOutput_sample.i 
				) {
					UVH5print_error(__FUNCTION__, "#(%d, %d): (xgpu = %d, blidx = %d, polidx = %d, isauto = %d, needsconj = %d)\n\t{%d, %di} vs {%d, %di}\n\t\t@ %d vs %d",
						approd_idx,
						freq,
						xgpu_idx[0],
						bl_idx[0],
						pol_idx[0],
						is_auto[0],
						conjugate[0],
						visdata[visdata_offset].r, visdata[visdata_offset].i,
					 	xgpuOutput_sample.r, xgpuOutput_sample.i,
						visdata_offset, freq * MATRIX_H_DIMS[1] + (*xgpu_idx)
					);
					failed = true;
				}
				if(*is_auto && *pol_idx == 1) {
					memcpy(&xgpuOutput_sample,
						xgpuOutput + (
							freq * MATRIX_H_DIMS[1] + (*xgpu_idx))-1,
						sizeof(UVH5_CF32_t)
					);
					if(*conjugate) {
						xgpuOutput_sample.i = -xgpuOutput_sample.i;
					}
					// If cross-pol autocorrelation:
					// "use some inside knowledge that
					//  should be publicized in XGPU documentation that the
					//  redundant cross-pol is at xgpuidx-1."
					if(visdata[visdata_offset+1].r != xgpuOutput_sample.r ||
					 visdata[visdata_offset+1].i != xgpuOutput_sample.i 
					) {
						UVH5print_error(__FUNCTION__, "%%(%d, %d): (xgpu = %d, blidx = %d, polidx = %d, isauto = %d, needsconj = %d)\n\t{%d, %di} vs {%d, %di} @ %d vs %d",
							approd_idx,
							freq,
							xgpu_idx[0],
							bl_idx[0],
							pol_idx[0],
							is_auto[0],
							conjugate[0],
							visdata[visdata_offset+1].r, visdata[visdata_offset+1].i,
							xgpuOutput_sample.r, xgpuOutput_sample.i,
							visdata_offset, freq * MATRIX_H_DIMS[1] + (*xgpu_idx)
						);
						failed = true;
					}
				}

				visdata_offset += Npols;
			}
			xgpu_idx++;
			bl_idx++;
			pol_idx++;
			conjugate++;
			is_auto++;
		}
	} // if compare

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
