#include "uvh5.h"

void UVH5Halloc(UVH5_header_t *header)
{
	if(header->ant_1_array == NULL && header->Nbls != 0) {
		header->ant_1_array = malloc(sizeof(int) * header->Nbls);
		memset(header->ant_1_array, 0, sizeof(int) * header->Nbls);
		UVH5print_verbose(__FUNCTION__, "'ant_1_array' allocated %ld bytes.", sizeof(int) * header->Nbls);
	}
	if(header->ant_2_array == NULL && header->Nbls != 0) {
		header->ant_2_array = malloc(sizeof(int) * header->Nbls);
		memset(header->ant_2_array, 0, sizeof(int) * header->Nbls);
		UVH5print_verbose(__FUNCTION__, "'ant_2_array' allocated %ld bytes.", sizeof(int) * header->Nbls);
	}
	if(header->antenna_numbers == NULL && header->Nants_telescope != 0) {
		header->antenna_numbers = malloc(sizeof(int) * header->Nants_telescope);
		memset(header->antenna_numbers, 0, sizeof(int) * header->Nants_telescope);
		UVH5print_verbose(__FUNCTION__, "'antenna_numbers' allocated %ld bytes.", sizeof(int) * header->Nants_telescope);
	}
	if(header->antenna_names == NULL && header->Nants_telescope != 0) {
		header->antenna_names = malloc(sizeof(char *) * header->Nants_telescope);
		memset(header->antenna_names, 0, sizeof(char *) * header->Nants_telescope);
		UVH5print_verbose(__FUNCTION__, "'antenna_names' allocated %ld bytes.", sizeof(char *) * header->Nants_telescope);
	}
	if(header->uvw_array == NULL && header->Nbls * 3 != 0) {
		header->uvw_array = malloc(sizeof(double) * header->Nbls * 3);
		memset(header->uvw_array, 0, sizeof(double) * header->Nbls * 3);
		UVH5print_verbose(__FUNCTION__, "'uvw_array' allocated %ld bytes.", sizeof(double) * header->Nbls * 3);
	}
	if(header->time_array == NULL && header->Nbls != 0) {
		header->time_array = malloc(sizeof(double) * header->Nbls);
		memset(header->time_array, 0, sizeof(double) * header->Nbls);
		UVH5print_verbose(__FUNCTION__, "'time_array' allocated %ld bytes.", sizeof(double) * header->Nbls);
	}
	if(header->integration_time == NULL && header->Nbls != 0) {
		header->integration_time = malloc(sizeof(double) * header->Nbls);
		memset(header->integration_time, 0, sizeof(double) * header->Nbls);
		UVH5print_verbose(__FUNCTION__, "'integration_time' allocated %ld bytes.", sizeof(double) * header->Nbls);
	}
	if(header->freq_array == NULL && header->Nfreqs != 0) {
		header->freq_array = malloc(sizeof(double) * header->Nfreqs);
		memset(header->freq_array, 0, sizeof(double) * header->Nfreqs);
		UVH5print_verbose(__FUNCTION__, "'freq_array' allocated %ld bytes.", sizeof(double) * header->Nfreqs);
	}
	if(header->channel_width == NULL && header->Nfreqs != 0) {
		header->channel_width = malloc(sizeof(double) * header->Nfreqs);
		memset(header->channel_width, 0, sizeof(double) * header->Nfreqs);
		UVH5print_verbose(__FUNCTION__, "'channel_width' allocated %ld bytes.", sizeof(double) * header->Nfreqs);
	}
	if(header->spw_array == NULL && header->Nspws != 0) {
		header->spw_array = malloc(sizeof(int) * header->Nspws);
		memset(header->spw_array, 0, sizeof(int) * header->Nspws);
		UVH5print_verbose(__FUNCTION__, "'spw_array' allocated %ld bytes.", sizeof(int) * header->Nspws);
	}
	if(header->polarization_array == NULL && header->Npols != 0) {
		header->polarization_array = malloc(sizeof(int) * header->Npols);
		memset(header->polarization_array, 0, sizeof(int) * header->Npols);
		UVH5print_verbose(__FUNCTION__, "'polarization_array' allocated %ld bytes.", sizeof(int) * header->Npols);
	}
	if(header->antenna_positions == NULL && header->Nants_telescope * 3 != 0) {
		header->antenna_positions = malloc(sizeof(double) * header->Nants_telescope * 3);
		memset(header->antenna_positions, 0, sizeof(double) * header->Nants_telescope * 3);
		UVH5print_verbose(__FUNCTION__, "'antenna_positions' allocated %ld bytes.", sizeof(double) * header->Nants_telescope * 3);
	}
	if(header->phase_center_id_array == NULL && header->Nbls != 0) {
		header->phase_center_id_array = malloc(sizeof(int) * header->Nbls);
		memset(header->phase_center_id_array, 0, sizeof(int) * header->Nbls);
		UVH5print_verbose(__FUNCTION__, "'phase_center_id_array' allocated %ld bytes.", sizeof(int) * header->Nbls);
	}
	if(header->phase_center_app_ra == NULL && header->Nbls != 0) {
		header->phase_center_app_ra = malloc(sizeof(double) * header->Nbls);
		memset(header->phase_center_app_ra, 0, sizeof(double) * header->Nbls);
		UVH5print_verbose(__FUNCTION__, "'phase_center_app_ra' allocated %ld bytes.", sizeof(double) * header->Nbls);
	}
	if(header->phase_center_app_dec == NULL && header->Nbls != 0) {
		header->phase_center_app_dec = malloc(sizeof(double) * header->Nbls);
		memset(header->phase_center_app_dec, 0, sizeof(double) * header->Nbls);
		UVH5print_verbose(__FUNCTION__, "'phase_center_app_dec' allocated %ld bytes.", sizeof(double) * header->Nbls);
	}
	if(header->phase_center_frame_pa == NULL && header->Nbls != 0) {
		header->phase_center_frame_pa = malloc(sizeof(double) * header->Nbls);
		memset(header->phase_center_frame_pa, 0, sizeof(double) * header->Nbls);
		UVH5print_verbose(__FUNCTION__, "'phase_center_frame_pa' allocated %ld bytes.", sizeof(double) * header->Nbls);
	}
	// Optional entries follow
	// if(header-> == NULL && header->Nbls != 0) {
	// 	header->lst_array = malloc(sizeof(double) * header->Nbls);
	// 	memset(header->lst_array, 0, sizeof(double) * header->Nbls);
	//	UVH5print_verbose(__FUNCTION__, "'lst_array' allocated %ld bytes.", sizeof(double) * header->Nbls);
	// }
	// Administrative entries follow
	if(header->_ant_pol_prod_xgpu_index == NULL && header->Nbls != 0 && header->Npols != 0) {
		header->_ant_pol_prod_xgpu_index = malloc(sizeof(int) * header->Nbls * header->Npols);
		memset(header->_ant_pol_prod_xgpu_index, 0, sizeof(int) * header->Nbls * header->Npols);
		UVH5print_verbose(__FUNCTION__, "'_ant_pol_prod_xgpu_index' allocated %ld bytes.", sizeof(int) * header->Nbls * header->Npols);
	}
	if(header->_ant_pol_prod_bl_index == NULL && header->Nbls != 0 && header->Npols != 0) {
		header->_ant_pol_prod_bl_index = malloc(sizeof(int) * header->Nbls * header->Npols);
		memset(header->_ant_pol_prod_bl_index, 0, sizeof(int) * header->Nbls * header->Npols);
		UVH5print_verbose(__FUNCTION__, "'_ant_pol_prod_bl_index' allocated %ld bytes.", sizeof(int) * header->Nbls * header->Npols);
	}
	if(header->_ant_pol_prod_pol_index == NULL && header->Nbls != 0 && header->Npols != 0) {
		header->_ant_pol_prod_pol_index = malloc(sizeof(int) * header->Nbls * header->Npols);
		memset(header->_ant_pol_prod_pol_index, 0, sizeof(int) * header->Nbls * header->Npols);
		UVH5print_verbose(__FUNCTION__, "'_ant_pol_prod_pol_index' allocated %ld bytes.", sizeof(int) * header->Nbls * header->Npols);
	}
	if(header->_ant_pol_prod_conj == NULL && header->Nbls != 0 && header->Npols != 0) {
		header->_ant_pol_prod_conj = malloc(sizeof(char) * header->Nbls * header->Npols);
		memset(header->_ant_pol_prod_conj, 0, sizeof(char) * header->Nbls * header->Npols);
		UVH5print_verbose(__FUNCTION__, "'_ant_pol_prod_conj' allocated %ld bytes.", sizeof(char) * header->Nbls * header->Npols);
	}
	if(header->_ant_pol_prod_auto == NULL && header->Nbls != 0 && header->Npols != 0) {
		header->_ant_pol_prod_auto = malloc(sizeof(char) * header->Nbls * header->Npols);
		memset(header->_ant_pol_prod_auto, 0, sizeof(char) * header->Nbls * header->Npols);
		UVH5print_verbose(__FUNCTION__, "'_ant_pol_prod_auto' allocated %ld bytes.", sizeof(char) * header->Nbls * header->Npols);
	}
	if(header->_antenna_numbers_data == NULL && header->Nants_data != 0) {
		header->_antenna_numbers_data = malloc(sizeof(int) * header->Nants_data);
		memset(header->_antenna_numbers_data, 0, sizeof(int) * header->Nants_data);
		UVH5print_verbose(__FUNCTION__, "'_antenna_numbers_data' allocated %ld bytes.", sizeof(int) * header->Nants_data);
	}
}

void UVH5Hmalloc_phase_center_catalog(UVH5_header_t *header, size_t catalog_length) {
	header->_phase_center_catalog_length = catalog_length;
	header->phase_center_catalog = malloc(header->_phase_center_catalog_length * sizeof(UVH5_phase_center_t));
	memset(header->phase_center_catalog, 0, header->_phase_center_catalog_length * sizeof(UVH5_phase_center_t));
}

/*
 * Generate `_antenna_num_idx_map` and `_antenna_enu_positions`.
 */
void UVH5Hadmin(UVH5_header_t *header) {
	int highest_antenna_number = 0;
	for (int i = 0; i < header->Nants_telescope; i++)
	{
		highest_antenna_number = header->antenna_numbers[i] > highest_antenna_number ?
			header->antenna_numbers[i] :
			highest_antenna_number
		;
	}
	
	//	Create antenna number -> index map
	header->_antenna_num_idx_map = malloc(sizeof(int) * (highest_antenna_number + 1));
	memset(header->_antenna_num_idx_map, -1, sizeof(int) * (highest_antenna_number + 1));
	for (size_t i = 0; i < header->Nants_telescope; i++) {
		header->_antenna_num_idx_map[header->antenna_numbers[i]] = i;
		UVH5print_verbose(__FUNCTION__, "Antenna number %d is at index %ld", header->antenna_numbers[i], i);
	}
	
	//	Create ENU antenna_positions from XYZ
	header->_antenna_enu_positions = malloc(sizeof(double) * header->Nants_telescope * 3);
	memcpy(header->_antenna_enu_positions, header->antenna_positions, sizeof(double) * header->Nants_telescope * 3);
	calc_position_to_enu_frame_from_xyz(
		header->_antenna_enu_positions,
		header->Nants_telescope,
		calc_rad_from_degree(header->longitude),
		calc_rad_from_degree(header->latitude),
		header->altitude
	);
	header->_antenna_uvw_positions = malloc(sizeof(double) * header->Nants_telescope * 3);
	memset(header->_antenna_uvw_positions, 0, sizeof(double) * header->Nants_telescope * 3);
}

void UVH5Hfree(UVH5_header_t *header)
{
	if (header->ant_1_array != NULL) {
		free(header->ant_1_array);
	}
	if (header->ant_2_array != NULL) {
		free(header->ant_2_array);
	}
	if (header->antenna_numbers != NULL) {
		free(header->antenna_numbers);
	}
	if (header->antenna_names != NULL) {
		free(header->antenna_names);
	}
	if (header->uvw_array != NULL) {
		free(header->uvw_array);
	}
	if (header->time_array != NULL) {
		free(header->time_array);
	}
	if (header->integration_time != NULL) {
		free(header->integration_time);
	}
	if (header->freq_array != NULL) {
		free(header->freq_array);
	}
	if (header->channel_width != NULL) {
		free(header->channel_width);
	}
	if (header->spw_array != NULL) {
		free(header->spw_array);
	}
	if (header->polarization_array != NULL) {
		free(header->polarization_array);
	}
	if (header->antenna_positions != NULL) {
		free(header->antenna_positions);
	}
	if (header->phase_center_catalog != NULL) {
		for(int catalog_index = 0; catalog_index < header->_phase_center_catalog_length; catalog_index++) {
			if(header->phase_center_catalog[catalog_index].type == UVH5_PHASE_CENTER_EPHEMERIS) {
				free(header->phase_center_catalog[catalog_index].ephem_lon);
				free(header->phase_center_catalog[catalog_index].ephem_lat);
				free(header->phase_center_catalog[catalog_index].times);
				free(header->phase_center_catalog[catalog_index].dist);
				free(header->phase_center_catalog[catalog_index].vrad);
			}
		}
		free(header->phase_center_catalog);
	}
	// Optional arrays follow
	if (header->antenna_diameters != NULL) {
		free(header->antenna_diameters);
	}
	// Administrative arrays follow
	if(header->_antenna_num_idx_map != NULL){
		free(header->_antenna_num_idx_map);
	}
	if(header->_antenna_enu_positions != NULL){
		free(header->_antenna_enu_positions);
	}
	if(header->_ant_pol_prod_xgpu_index != NULL){
		free(header->_ant_pol_prod_xgpu_index);
	}
	if(header->_ant_pol_prod_bl_index != NULL){
		free(header->_ant_pol_prod_bl_index);
	}
	if(header->_ant_pol_prod_pol_index != NULL){
		free(header->_ant_pol_prod_pol_index);
	}
	if(header->_ant_pol_prod_conj != NULL){
		free(header->_ant_pol_prod_conj);
	}
	if(header->_ant_pol_prod_auto != NULL){
		free(header->_ant_pol_prod_auto);
	}
	if(header->_antenna_numbers_data != NULL){
		free(header->_antenna_numbers_data);
	}
}

hid_t UVH5TcreateCI32()
{
	hid_t UVH5_CI32 = H5Tcreate(H5T_COMPOUND, sizeof(UVH5_CI32_t));
	H5Tinsert(UVH5_CI32, "r", HOFFSET(UVH5_CI32_t, r), H5T_NATIVE_INT32);
	H5Tinsert(UVH5_CI32, "i", HOFFSET(UVH5_CI32_t, i), H5T_NATIVE_INT32);
	return UVH5_CI32;
}

hid_t UVH5TcreateCF32()
{
	hid_t UVH5_CF32 = H5Tcreate(H5T_COMPOUND, sizeof(UVH5_CF32_t));
	H5Tinsert(UVH5_CF32, "r", HOFFSET(UVH5_CF32_t, r), H5T_NATIVE_FLOAT);
	H5Tinsert(UVH5_CF32, "i", HOFFSET(UVH5_CF32_t, i), H5T_NATIVE_FLOAT);
	return UVH5_CF32;
}

hid_t UVH5TcreateCF64()
{
	hid_t UVH5_CF64 = H5Tcreate(H5T_COMPOUND, sizeof(UVH5_CF64_t));
	H5Tinsert(UVH5_CF64, "r", HOFFSET(UVH5_CF64_t, r), H5T_NATIVE_DOUBLE);
	H5Tinsert(UVH5_CF64, "i", HOFFSET(UVH5_CF64_t, i), H5T_NATIVE_DOUBLE);
	return UVH5_CF64;
}

herr_t _H5DSphaseCenterCatalogWrite(
	hid_t file_id,
	UVH5_phase_center_t* catalog,
	int catalog_length
) {	
	herr_t status;
	
	hid_t parent_catalog_id = H5Gcreate(file_id, "/Header/phase_center_catalog", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (parent_catalog_id < 0)	{ UVH5print_error(__FUNCTION__, "failure opening catalog"); return -1; }

	// 'Header/phase_center_catalog' group
	char catalog_group_path[] = "/Header/phase_center_catalog/*********";
	for(int catalog_index = 0; catalog_index < catalog_length; catalog_index++) {
		sprintf(catalog_group_path, "/Header/phase_center_catalog/%d", catalog_index%1000000000);
		hid_t catalog_id = H5Gcreate(file_id, catalog_group_path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (catalog_id < 0)	{ UVH5print_error(__FUNCTION__, "failure opening catalog ID #%d", catalog_index); return -1; }

		status = H5DSstringWrite(catalog_id, "cat_name", 0, NULL, catalog[catalog_index].name);
		if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'cat_name' at catalog ID #%d", catalog_index); return status; }
		status = H5DSstringWrite(catalog_id, "cat_frame", 0, NULL, catalog[catalog_index].frame);
		if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'cat_frame' at catalog ID #%d", catalog_index); return status; }
		status = H5DSdoubleWrite(catalog_id, "cat_epoch", 0, NULL, &catalog[catalog_index].epoch);
		if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'cat_epoch' at catalog ID #%d", catalog_index); return status; }
		if(catalog[catalog_index].info_source) {
			status = H5DSstringWrite(catalog_id, "info_source", 0, NULL, catalog[catalog_index].info_source);
			if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'info_source' at catalog ID #%d", catalog_index); return status; }
		}

		switch(catalog[catalog_index].type) {
			case UVH5_PHASE_CENTER_SIDEREAL:
				status = H5DSstringWrite(catalog_id, "cat_type", 0, NULL, "sidereal");
				if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'cat_type' at catalog ID #%d", catalog_index); return status; }
				status = H5DSdoubleWrite(catalog_id, "cat_lon", 0, NULL, &catalog[catalog_index].lon);
				if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'cat_lon' at catalog ID #%d", catalog_index); return status; }
				status = H5DSdoubleWrite(catalog_id, "cat_lat", 0, NULL, &catalog[catalog_index].lat);
				if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'cat_lat' at catalog ID #%d", catalog_index); return status; }
				status = H5DSdoubleWrite(catalog_id, "pm_ra", 0, NULL, &catalog[catalog_index].pm_ra);
				if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'pm_ra' at catalog ID #%d", catalog_index); return status; }
				status = H5DSdoubleWrite(catalog_id, "pm_dec", 0, NULL, &catalog[catalog_index].pm_dec);
				if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'pm_dec' at catalog ID #%d", catalog_index); return status; }
				
				break;
			case UVH5_PHASE_CENTER_EPHEMERIS:
				status = H5DSstringWrite(catalog_id, "cat_type", 0, NULL, "ephem");
				if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'cat_type' at catalog ID #%d", catalog_index); return status; }

				const hsize_t dims1_Npts[] = {catalog[catalog_index]._ephem_points};
				status = H5DSdoubleWrite(catalog_id, "ephem_lon", 1, dims1_Npts, catalog[catalog_index].ephem_lon);
				if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'ephem_lon' at catalog ID #%d", catalog_index); return status; }
				status = H5DSdoubleWrite(catalog_id, "ephem_lat", 1, dims1_Npts, catalog[catalog_index].ephem_lat);
				if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'ephem_lat' at catalog ID #%d", catalog_index); return status; }
				status = H5DSdoubleWrite(catalog_id, "times", 1, dims1_Npts, catalog[catalog_index].times);
				if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'times' at catalog ID #%d", catalog_index); return status; }
				status = H5DSdoubleWrite(catalog_id, "dist", 1, dims1_Npts, catalog[catalog_index].dist);
				if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'dist' at catalog ID #%d", catalog_index); return status; }
				status = H5DSdoubleWrite(catalog_id, "vrad", 1, dims1_Npts, catalog[catalog_index].vrad);
				if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'vrad' at catalog ID #%d", catalog_index); return status; }
				break;
			case UVH5_PHASE_CENTER_DRIFTSCAN:
				status = H5DSstringWrite(catalog_id, "cat_type", 0, NULL, "driftscan");
				if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'cat_type' at catalog ID #%d", catalog_index); return status; }
				status = H5DSdoubleWrite(catalog_id, "cat_lon", 0, NULL, &catalog[catalog_index].lon);
				if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'cat_lon' at catalog ID #%d", catalog_index); return status; }
				status = H5DSdoubleWrite(catalog_id, "cat_lat", 0, NULL, &catalog[catalog_index].lat);
				if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'cat_lat' at catalog ID #%d", catalog_index); return status; }
				break;
			case UVH5_PHASE_CENTER_UNPROJECTED:
				status = H5DSstringWrite(catalog_id, "cat_type", 0, NULL, "unprojected");
				if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'cat_type' at catalog ID #%d", catalog_index); return status; }
				status = H5DSdoubleWrite(catalog_id, "cat_lon", 0, NULL, &catalog[catalog_index].lon);
				if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'cat_lon' at catalog ID #%d", catalog_index); return status; }
				status = H5DSdoubleWrite(catalog_id, "cat_lat", 0, NULL, &catalog[catalog_index].lat);
				if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'cat_lat' at catalog ID #%d", catalog_index); return status; }
				break;
			default:
				UVH5print_error(__FUNCTION__, "Unrecognised catalog type '%d' for catalog ID #%d", catalog[catalog_index].type, catalog_index);
				return -1;
		}

		status = H5Gclose(catalog_id);
		if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure closing catalog ID #%d", catalog_index); return status; }
	}
	status = H5Gclose(parent_catalog_id);
	if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure closing catalog"); return status; }

	return 0;	
}


void _UVH5_Hwrite_static(UVH5_file_t *UVH5file)
{
	// TODO return error code
	herr_t status;

	UVH5_header_t header = UVH5file->header;
	const hsize_t dims1_Nants_telescope[] = {header.Nants_telescope};
	const hsize_t dims1_Nspws[] = {header.Nspws};
	const hsize_t dims1_Nfreqs[] = {header.Nfreqs};
	const hsize_t dims1_Npols[] = {header.Npols};
	const hsize_t dims2_Nants_telescope_3[] = {header.Nants_telescope, 3};

	// TODO assert UVH5file->header_id

	status = H5DSdoubleWrite(UVH5file->header_id, "latitude", 0, NULL, &header.latitude);
	if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'latitude'"); return; }

	status = H5DSdoubleWrite(UVH5file->header_id, "longitude", 0, NULL, &header.longitude);
	if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'longitude'"); return; }

	status = H5DSdoubleWrite(UVH5file->header_id, "altitude", 0, NULL, &header.altitude);
	if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'altitude'"); return; }

	status = H5DSstringWrite(UVH5file->header_id, "telescope_name", 0, NULL, header.telescope_name);
	if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'telescope_name'"); return; }

	status = H5DSstringWrite(UVH5file->header_id, "instrument", 0, NULL, header.instrument);
	if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'instrument'"); return; }

	status = H5DSstringWrite(UVH5file->header_id, "history", 0, NULL, header.history);
	if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'history'"); return; }

	status = H5DSintWrite(UVH5file->header_id, "Nants_data", 0, NULL, &header.Nants_data);
	if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'Nants_data'"); return; }

	status = H5DSintWrite(UVH5file->header_id, "Nants_telescope", 0, NULL, &header.Nants_telescope);
	if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'Nants_telescope'"); return; }

	status = H5DSstringWrite(UVH5file->header_id, "antenna_names", 1, dims1_Nants_telescope, (char *)header.antenna_names);
	if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'antenna_names'"); return; }

	status = H5DSintWrite(UVH5file->header_id, "antenna_numbers", 1, dims1_Nants_telescope, header.antenna_numbers);
	if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'antenna_numbers'"); return; }

	status = H5DSintWrite(UVH5file->header_id, "Nbls", 0, NULL, &header.Nbls);
	if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'Nbls'"); return; }

	status = H5DSintWrite(UVH5file->header_id, "Nspws", 0, NULL, &header.Nspws);
	if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'Nspws'"); return; }

	status = H5DSintWrite(UVH5file->header_id, "Nfreqs", 0, NULL, &header.Nfreqs);
	if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'Nfreqs'"); return; }

	status = H5DSintWrite(UVH5file->header_id, "Npols", 0, NULL, &header.Npols);
	if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'Npols'"); return; }

	status = H5DSintWrite(UVH5file->header_id, "Nphase", 0, NULL, &header._phase_center_catalog_length);
	if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'Nphase'"); return; }

	status = H5DSdoubleWrite(UVH5file->header_id, "freq_array", 1, dims1_Nfreqs, header.freq_array);
	if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'freq_array'"); return; }

	status = H5DSdoubleWrite(UVH5file->header_id, "channel_width", 1, dims1_Nfreqs, header.channel_width);
	if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'channel_width'"); return; }

	status = H5DSintWrite(UVH5file->header_id, "spw_array", 1, dims1_Nspws, header.spw_array);
	if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'spw_array'"); return; }

	status = H5DSboolWrite(UVH5file->header_id, "flex_spw", 0, NULL, &header.flex_spw);
	if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'flex_spw'"); return; }

	status = H5DSintWrite(UVH5file->header_id, "polarization_array", 1, dims1_Npols, header.polarization_array);
	if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'polarization_array'"); return; }

	status = H5DSdoubleWrite(UVH5file->header_id, "antenna_positions", 2, dims2_Nants_telescope_3, header.antenna_positions);
	if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'antenna_positions'"); return; }

	status = _H5DSphaseCenterCatalogWrite(UVH5file->file_id, header.phase_center_catalog, header._phase_center_catalog_length);
	if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'phase_center_catalog'"); return; }

	status = H5DSstringWrite(UVH5file->header_id, "version", 0, NULL, "1.1");
	if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on 'version'"); return; }

	// Optional entries follow
	if(header.flex_spw_id_array){
		status = H5DSintWrite(UVH5file->header_id, "flex_spw_id_array", 1, dims1_Nfreqs, header.flex_spw_id_array);
		if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on optional 'flex_spw_id_array'"); }
	}

	if(header.dut1){
		status = H5DSdoubleWrite(UVH5file->header_id, "dut1", 0, NULL, &header.dut1);
		if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on optional 'dut1'"); }
	}

	if(header.earth_omega){
		status = H5DSdoubleWrite(UVH5file->header_id, "earth_omega", 0, NULL, &header.earth_omega);
		if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on optional 'earth_omega'"); }
	}

	if(header.gst0){
		status = H5DSdoubleWrite(UVH5file->header_id, "gst0", 0, NULL, &header.gst0);
		if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on optional 'gst0'"); }
	}

	if(header.rdate){
		status = H5DSstringWrite(UVH5file->header_id, "rdate", 0, NULL, header.rdate);
		if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on optional 'rdate'"); }
	}

	if(header.timesys){
		status = H5DSstringWrite(UVH5file->header_id, "timesys", 0, NULL, header.timesys);
		if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on optional 'timesys'"); }
	}

	if(header.x_orientation){
		status = H5DSstringWrite(UVH5file->header_id, "x_orientation", 0, NULL, header.x_orientation);
		if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on optional 'x_orientation'"); }
	}

	if(header.antenna_diameters){
		status = H5DSdoubleWrite(UVH5file->header_id, "antenna_diameters", 1, dims1_Nants_telescope, header.antenna_diameters);
		if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on optional 'antenna_diameters'"); }
	}

	if(header.uvplane_reference_time){
		status = H5DSintWrite(UVH5file->header_id, "uvplane_reference_time", 0, NULL, &header.uvplane_reference_time);
		if (status < 0)	{ UVH5print_error(__FUNCTION__, "failure on optional 'uvplane_reference_time'"); }
	}

}

void UVH5open(char* filepath, UVH5_file_t *UVH5file, hid_t Tvisdata)
{

	UVH5file->file_id = H5Fcreate(filepath, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	const hsize_t dim1_unlim[] = {H5S_UNLIMITED};
	const hsize_t dim1_nbls[] = {UVH5file->header.Nbls};
	const hsize_t dim2_unlim_3[] = {H5S_UNLIMITED, 3};
	const hsize_t dim2_nbls_3[] = {UVH5file->header.Nbls, 3};

	const hsize_t dim3_data_lim[] = {H5S_UNLIMITED, UVH5file->header.Nfreqs, UVH5file->header.Npols};
	const hsize_t dim3_data_chunk[] = {UVH5file->header.Nbls, UVH5file->header.Nfreqs, UVH5file->header.Npols};

	// 'Header' group
	UVH5file->header_id = H5Gcreate(UVH5file->file_id, "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	UVH5file->DS_header_Ntimes.name = "Ntimes";
	H5DSset(0, NULL, NULL, &UVH5file->DS_header_Ntimes);
	H5DSopenInt(UVH5file->header_id, &UVH5file->DS_header_Ntimes);

	UVH5file->DS_header_Nblts.name = "Nblts";
	H5DSset(0, NULL, NULL, &UVH5file->DS_header_Nblts);
	H5DSopenInt(UVH5file->header_id, &UVH5file->DS_header_Nblts);

	UVH5file->DS_header_ant_1_array.name = "ant_1_array";
	H5DSset(1, dim1_unlim, dim1_nbls, &UVH5file->DS_header_ant_1_array);
	H5DSopenInt(UVH5file->header_id, &UVH5file->DS_header_ant_1_array);

	UVH5file->DS_header_ant_2_array.name = "ant_2_array";
	H5DSset(1, dim1_unlim, dim1_nbls, &UVH5file->DS_header_ant_2_array);
	H5DSopenInt(UVH5file->header_id, &UVH5file->DS_header_ant_2_array);

	UVH5file->DS_header_uvw_array.name = "uvw_array";
	H5DSset(2, dim2_unlim_3, dim2_nbls_3, &UVH5file->DS_header_uvw_array);
	H5DSopenDouble(UVH5file->header_id, &UVH5file->DS_header_uvw_array);

	UVH5file->DS_header_time_array.name = "time_array";
	H5DSset(1, dim1_unlim, dim1_nbls, &UVH5file->DS_header_time_array);
	H5DSopenDouble(UVH5file->header_id, &UVH5file->DS_header_time_array);

	UVH5file->DS_header_integration_time.name = "integration_time";
	H5DSset(1, dim1_unlim, dim1_nbls, &UVH5file->DS_header_integration_time);
	H5DSopenDouble(UVH5file->header_id, &UVH5file->DS_header_integration_time);

	UVH5file->DS_phase_center_id_array.name = "phase_center_id_array";
	H5DSset(1, dim1_unlim, dim1_nbls, &UVH5file->DS_phase_center_id_array);
	H5DSopenInt(UVH5file->header_id, &UVH5file->DS_phase_center_id_array);

	UVH5file->DS_phase_center_app_ra.name = "phase_center_app_ra";
	H5DSset(1, dim1_unlim, dim1_nbls, &UVH5file->DS_phase_center_app_ra);
	H5DSopenDouble(UVH5file->header_id, &UVH5file->DS_phase_center_app_ra);

	UVH5file->DS_phase_center_app_dec.name = "phase_center_app_dec";
	H5DSset(1, dim1_unlim, dim1_nbls, &UVH5file->DS_phase_center_app_dec);
	H5DSopenDouble(UVH5file->header_id, &UVH5file->DS_phase_center_app_dec);

	UVH5file->DS_phase_center_frame_pa.name = "phase_center_frame_pa";
	H5DSset(1, dim1_unlim, dim1_nbls, &UVH5file->DS_phase_center_frame_pa);
	H5DSopenDouble(UVH5file->header_id, &UVH5file->DS_phase_center_frame_pa);

	if(UVH5file->header.lst_array){
		UVH5file->DS_header_lst_array.name = "lst_array";
		H5DSset(1, dim1_unlim, dim1_nbls, &UVH5file->DS_header_lst_array);
		H5DSopenDouble(UVH5file->header_id, &UVH5file->DS_header_lst_array);
	}

	_UVH5_Hwrite_static(UVH5file);
	int zero = 0;
	H5DSwrite(&UVH5file->DS_header_Ntimes, &zero);
	H5DSwrite(&UVH5file->DS_header_Nblts, &zero);

	// 'Data' group
	UVH5file->data_id = H5Gcreate(UVH5file->file_id, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	UVH5file->DS_data_visdata.name = "visdata";
	H5DSset(3, dim3_data_lim, dim3_data_chunk, &UVH5file->DS_data_visdata);
	H5DSopen(UVH5file->data_id, Tvisdata, H5Tcopy(Tvisdata), &UVH5file->DS_data_visdata);
	UVH5file->visdata = H5DSmalloc(&UVH5file->DS_data_visdata);

	UVH5file->DS_data_flags.name = "flags";
	H5DSset(3, dim3_data_lim, dim3_data_chunk, &UVH5file->DS_data_flags);
	UVH5file->DS_data_flags.filter_flag = H5_FILTER_FLAG_DEFLATE_3;
	H5DSopenBool(UVH5file->data_id, &UVH5file->DS_data_flags);
	UVH5file->flags = H5DSmalloc(&UVH5file->DS_data_flags);

	UVH5file->DS_data_nsamples.name = "nsamples";
	H5DSset(3, dim3_data_lim, dim3_data_chunk, &UVH5file->DS_data_nsamples);
	UVH5file->DS_data_flags.filter_flag = H5_FILTER_FLAG_DEFLATE_3;
	H5DSopenDouble(UVH5file->data_id, &UVH5file->DS_data_nsamples);
	UVH5file->nsamples = H5DSmalloc(&UVH5file->DS_data_nsamples);

	// 'extra_keywords' group
	UVH5file->keywords_id = H5Gcreate(UVH5file->file_id, "/Header/extra_keywords", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}

void UVH5close(UVH5_file_t *UVH5file)
{
	// 'Header' group
	if(UVH5file->header_id) {
		H5DSclose(&UVH5file->DS_header_Ntimes);
		H5DSclose(&UVH5file->DS_header_Nblts);
		H5DSclose(&UVH5file->DS_header_ant_1_array);
		H5DSclose(&UVH5file->DS_header_ant_2_array);
		H5DSclose(&UVH5file->DS_header_uvw_array);
		H5DSclose(&UVH5file->DS_header_time_array);
		H5DSclose(&UVH5file->DS_header_integration_time);
		H5DSclose(&UVH5file->DS_phase_center_id_array);
		H5DSclose(&UVH5file->DS_phase_center_app_ra);
		H5DSclose(&UVH5file->DS_phase_center_app_dec);
		H5DSclose(&UVH5file->DS_phase_center_frame_pa);

		if(UVH5file->header.lst_array){
			H5DSclose(&UVH5file->DS_header_lst_array);
		}
		H5Gclose(UVH5file->header_id);
	}
	UVH5Hfree(&UVH5file->header);

	// 'Data' group
	if(UVH5file->data_id) {
		H5DSclose(&UVH5file->DS_data_visdata);
		H5DSclose(&UVH5file->DS_data_flags);
		H5DSclose(&UVH5file->DS_data_nsamples);

		H5Gclose(UVH5file->data_id);
		free(UVH5file->visdata);
		free(UVH5file->flags);
		free(UVH5file->nsamples);
	}

	// 'extra_keywords' group
	if(UVH5file->keywords_id) {
		H5Gclose(UVH5file->keywords_id);
	}

	if(UVH5file->file_id) {
		H5Fclose(UVH5file->file_id);
		UVH5file->file_id = 0;
	}
}

int UVH5write_keyword_bool(UVH5_file_t* UVH5file, char* key, bool value) {
	// H5_bool_t value_ = value == 0 ? H5_FALSE : H5_TRUE;
	if (H5DSboolWrite(UVH5file->keywords_id, key, 0, NULL, &value) < 0) {
		UVH5print_error(__FUNCTION__, "failure on `extra_keywords/%s`", key);
		return -1;
	}
	return 0;
}

int UVH5write_keyword_double(UVH5_file_t* UVH5file, char* key, double value) {
	if (H5DSdoubleWrite(UVH5file->keywords_id, key, 0, NULL, &value) < 0) {
		UVH5print_error(__FUNCTION__, "failure on `extra_keywords/%s`", key);
		return -1;
	}
	return 0;
}

int UVH5write_keyword_int(UVH5_file_t* UVH5file, char* key, int value) {
	if (H5DSintWrite(UVH5file->keywords_id, key, 0, NULL, &value) < 0) {
		UVH5print_error(__FUNCTION__, "failure on `extra_keywords/%s`", key);
		return -1;
	}
	return 0;
}

int UVH5write_keyword_string(UVH5_file_t* UVH5file, char* key, char* value) {
	if (H5DSstringWrite(UVH5file->keywords_id, key, 0, NULL, value) < 0) {
		UVH5print_error(__FUNCTION__, "failure on `extra_keywords/%s`", key);
		return -1;
	}
	return 0;
}

int UVH5write_dynamic(UVH5_file_t* UVH5file) {
	herr_t status;
	// Header
	status = H5DSextend(&UVH5file->DS_header_ant_1_array);
	if (status < 0) { UVH5print_error(__FUNCTION__, "H5DSextend failure on 'header_ant_1_array'"); return -1;}
	status = H5DSwrite(&UVH5file->DS_header_ant_1_array, UVH5file->header.ant_1_array);
	if (status < 0) { UVH5print_error(__FUNCTION__, "H5DSwrite failure on 'header_ant_1_array'"); return -1;}

	status = H5DSextend(&UVH5file->DS_header_ant_2_array);
	if (status < 0) { UVH5print_error(__FUNCTION__, "H5DSextend failure on 'header_ant_2_array'"); return -1;}
	status = H5DSwrite(&UVH5file->DS_header_ant_2_array, UVH5file->header.ant_2_array);
	if (status < 0) { UVH5print_error(__FUNCTION__, "H5DSwrite failure on 'header_ant_2_array'"); return -1;}

	status = H5DSextend(&UVH5file->DS_header_uvw_array);
	if (status < 0) { UVH5print_error(__FUNCTION__, "H5DSextend failure on 'header_uvw_array'"); return -1;}
	status = H5DSwrite(&UVH5file->DS_header_uvw_array, UVH5file->header.uvw_array);
	if (status < 0) { UVH5print_error(__FUNCTION__, "H5DSwrite failure on 'header_uvw_array'"); return -1;}

	status = H5DSextend(&UVH5file->DS_header_time_array);
	if (status < 0) { UVH5print_error(__FUNCTION__, "H5DSextend failure on 'header_time_array'"); return -1;}
	status = H5DSwrite(&UVH5file->DS_header_time_array, UVH5file->header.time_array);
	if (status < 0) { UVH5print_error(__FUNCTION__, "H5DSwrite failure on 'header_time_array'"); return -1;}

	status = H5DSextend(&UVH5file->DS_header_integration_time);
	if (status < 0) { UVH5print_error(__FUNCTION__, "H5DSextend failure on 'header_integration_time'"); return -1;}
	status = H5DSwrite(&UVH5file->DS_header_integration_time, UVH5file->header.integration_time);
	if (status < 0) { UVH5print_error(__FUNCTION__, "H5DSwrite failure on 'header_integration_time'"); return -1;}

	status = H5DSextend(&UVH5file->DS_phase_center_id_array);
	if (status < 0) { UVH5print_error(__FUNCTION__, "H5DSextend failure on 'header_phase_center_id_array'"); return -1;}
	status = H5DSwrite(&UVH5file->DS_phase_center_id_array, UVH5file->header.phase_center_id_array);
	if (status < 0) { UVH5print_error(__FUNCTION__, "H5DSwrite failure on 'header_phase_center_id_array'"); return -1;}

	status = H5DSextend(&UVH5file->DS_phase_center_app_ra);
	if (status < 0) { UVH5print_error(__FUNCTION__, "H5DSextend failure on 'header_phase_center_app_ra'"); return -1;}
	status = H5DSwrite(&UVH5file->DS_phase_center_app_ra, UVH5file->header.phase_center_app_ra);
	if (status < 0) { UVH5print_error(__FUNCTION__, "H5DSwrite failure on 'header_phase_center_app_ra'"); return -1;}

	status = H5DSextend(&UVH5file->DS_phase_center_app_dec);
	if (status < 0) { UVH5print_error(__FUNCTION__, "H5DSextend failure on 'header_phase_center_app_dec'"); return -1;}
	status = H5DSwrite(&UVH5file->DS_phase_center_app_dec, UVH5file->header.phase_center_app_dec);
	if (status < 0) { UVH5print_error(__FUNCTION__, "H5DSwrite failure on 'header_phase_center_app_dec'"); return -1;}

	status = H5DSextend(&UVH5file->DS_phase_center_frame_pa);
	if (status < 0) { UVH5print_error(__FUNCTION__, "H5DSextend failure on 'header_phase_center_frame_pa'"); return -1;}
	status = H5DSwrite(&UVH5file->DS_phase_center_frame_pa, UVH5file->header.phase_center_frame_pa);
	if (status < 0) { UVH5print_error(__FUNCTION__, "H5DSwrite failure on 'header_phase_center_frame_pa'"); return -1;}

	// Data
	status = H5DSextend(&UVH5file->DS_data_visdata);
	if (status < 0) { UVH5print_error(__FUNCTION__, "H5DSextend failure on 'data_visdata'"); return -1;}
	status = H5DSwrite(&UVH5file->DS_data_visdata, UVH5file->visdata);
	if (status < 0) { UVH5print_error(__FUNCTION__, "H5DSwrite failure on 'data_visdata'"); return -1;}

	status = H5DSextend(&UVH5file->DS_data_flags);
	if (status < 0) { UVH5print_error(__FUNCTION__, "H5DSextend failure on 'data_flags'"); return -1;}
	status = H5DSwrite(&UVH5file->DS_data_flags, UVH5file->flags);
	if (status < 0) { UVH5print_error(__FUNCTION__, "H5DSwrite failure on 'data_flags'"); return -1;}

	status = H5DSextend(&UVH5file->DS_data_nsamples);
	if (status < 0) { UVH5print_error(__FUNCTION__, "H5DSextend failure on 'data_nsamples'"); return -1;}
	status = H5DSwrite(&UVH5file->DS_data_nsamples, UVH5file->nsamples);
	if (status < 0) { UVH5print_error(__FUNCTION__, "H5DSwrite failure on 'data_nsamples'"); return -1;}

	// Update static
	UVH5file->header.Ntimes += 1;
	UVH5file->header.Nblts += UVH5file->header.Nbls;
	status = H5DSwrite(&UVH5file->DS_header_Ntimes, &UVH5file->header.Ntimes);
	if (status < 0)	{ UVH5print_error(__FUNCTION__, "H5DSwrite failure on 'Ntimes'"); return -1; }

	status = H5DSwrite(&UVH5file->DS_header_Nblts, &UVH5file->header.Nblts);
	if (status < 0)	{ UVH5print_error(__FUNCTION__, "H5DSwrite failure on 'Nblts'"); return -1; }

	return 0;
}

int UVH5find_antenna_index_by_name(UVH5_header_t* header, char* name) {
	if(header->antenna_names == NULL) {
			UVH5print_error(__FUNCTION__, "UVH5_header_t->antenna_names is NULL");
			return -2;
	}
	for(size_t i = 0; i < header->Nants_telescope; i++) {
		if(strcmp(name, header->antenna_names[i]) == 0) { // Case sensitive
			UVH5print_verbose(__FUNCTION__, "Antenna '%s' is at index %ld.", name, i);
			return i;
		}
	}
	UVH5print_warn(__FUNCTION__, "UVH5_header_t->antenna_names does not contain '%s'.", name);
	return -1;
}

/*
 * https://github.com/david-macmahon/XGPU.jl/blob/1f08c9eb24268b0bd94659f4da78982ddb829ff0/src/XGPU.jl#L718-L732
 */
int xgpuInputPairOutputIndex(
	int ant_1_i, int ant_1_p,
	int ant_2_i, int ant_2_p,
	int npol
) {
	if (ant_1_i*npol + ant_1_p > ant_2_i*npol + ant_2_p) {
		int tmp = ant_1_i;
		ant_1_i = ant_2_i;
		ant_2_i = tmp;
		tmp = ant_1_p;
		ant_1_p = ant_2_p;
		ant_2_p = tmp;
	}

  int station_pair_index = (((ant_2_i+1)*ant_2_i) / 2) + ant_1_i;
  int pol_pair_index = ant_2_p * npol + ant_1_p;
  return station_pair_index * npol * npol + pol_pair_index;
}

/*
* This expects that the header has been alloced with:
*		header->Npols
*		header->Nants_data
*		header->Nbls
*
* The UVh5_inputpair_t array is used to determine actual information of the observation's
* baselines, expressed and captured as pairs of ant_numbers.
* This enables the population of `header->ant_1/2_array`.
* Furthermore, an index is populated for each baseline in the internal administrative
* arrays, which are critical to the `UVH5visdata_from_xgpu_int_output` function:
*  `header->_ant_pol_prod_xgpu_index`
*  `header->_ant_pol_prod_bl_index`
*  `header->_ant_pol_prod_pol_index`
*  `header->_ant_pol_prod_auto`
*  `header->_ant_pol_prod_conj`
*/
void UVH5parse_input_map(
	UVH5_header_t* header,
	UVH5_inputpair_t* inputs
) {
	int npols_in = header->Npols == 4 ? 2 : 1; // header->Npols is the pol-products
	int ninpairs = npols_in*header->Nants_data;

	//positive offset by most negative -8
	int pol_product_index_map_offset8[13] = {-1};
	for(int pol_idx = 0; pol_idx < header->Npols; pol_idx++) {
		pol_product_index_map_offset8[8+header->polarization_array[pol_idx]] = pol_idx;
	}
	char pol_product[3] = {'\0'};
	
	int* redundant_pol_index = malloc(sizeof(int)*header->Nants_data);
	memset(redundant_pol_index, -1, sizeof(int)*header->Nants_data);

	// Set the unique antenna_numbers for those in data
	int _ant_numbers_data_index = 0;

	int auto_bl_idx = -1;
	bool is_auto;
	int cross_bl_idx = header->Nants_data - 1;
	int cross_bl_idx_rerun_from;
	int idx = 0;
	int ant_1_idx, ant_1_numidx, ant_1_num, ant_2_idx, ant_2_numidx, ant_2_num;
	for(int inpair_1 = 0; inpair_1 < ninpairs; inpair_1++) {
		ant_1_idx = inpair_1/npols_in;
		pol_product[0] = inputs[inpair_1].polarization;
		
		if(inpair_1%npols_in == 0) {
			ant_1_numidx = UVH5find_antenna_index_by_name(header, inputs[inpair_1].antenna);
			if(ant_1_numidx >= 0) {
				ant_1_num = header->antenna_numbers[ant_1_numidx];
			}
			else {
				ant_1_num = -1;
			}
			header->_antenna_numbers_data[_ant_numbers_data_index] = ant_1_num;
			_ant_numbers_data_index++;

			cross_bl_idx_rerun_from = cross_bl_idx;
		}
		else {
			// rewind the cross_bl_idx
			cross_bl_idx = cross_bl_idx_rerun_from;
		}

		for(int inpair_2 = ant_1_idx*npols_in; inpair_2 < ninpairs; inpair_2++) {
			ant_2_idx = inpair_2/npols_in;
			pol_product[1] = inputs[inpair_2].polarization;
			if(inpair_2%npols_in == 0) {
				ant_2_numidx = UVH5find_antenna_index_by_name(header, inputs[inpair_2].antenna);
				if(ant_2_numidx >= 0) {
					ant_2_num = header->antenna_numbers[ant_2_numidx];
				}
				else {
					ant_2_num = -1;
				}
			}

			is_auto = ant_1_num == ant_2_num;

			if(inpair_1%npols_in == 0 && inpair_2%npols_in == 0) {
				if(is_auto) {
					auto_bl_idx++;
					header->ant_1_array[auto_bl_idx] = ant_1_num;
					header->ant_2_array[auto_bl_idx] = ant_2_num;
					UVH5print_verbose(__FUNCTION__, "Auto baseline #%d %s(%d)->(%d)%s", auto_bl_idx, inputs[inpair_1].antenna, ant_1_num, ant_2_num, inputs[inpair_2].antenna);
				}
				else {
					cross_bl_idx++;
					header->ant_1_array[cross_bl_idx] = ant_1_num;
					header->ant_2_array[cross_bl_idx] = ant_2_num;
					UVH5print_verbose(__FUNCTION__, "Cross baseline #%d %s(%d)->(%d)%s", cross_bl_idx, inputs[inpair_1].antenna, ant_1_num, ant_2_num, inputs[inpair_2].antenna);
				}
			}
			else if (!is_auto && inpair_1%npols_in != 0 && inpair_2%npols_in == 0) {
				// Re-run the cross baselines
				cross_bl_idx++;
			}

			header->_ant_pol_prod_auto[idx] = is_auto;
			header->_ant_pol_prod_xgpu_index[idx] = xgpuInputPairOutputIndex(
				ant_1_idx, inpair_1%npols_in,
				ant_2_idx, inpair_2%npols_in,
				npols_in
			);

			header->_ant_pol_prod_bl_index[idx] = is_auto ? auto_bl_idx : cross_bl_idx;
			header->_ant_pol_prod_pol_index[idx] = pol_product_index_map_offset8[
				8+UVH5polarisation_string_key(pol_product, npols_in)
			];
			header->_ant_pol_prod_conj[idx] = ! (is_auto || (
				(ant_1_idx < header->Nants_data-1) &&
				(ant_2_idx > ant_1_idx)
			));

			// If cross-pol autocorrelation:
			// "use some inside knowledge that
			//  should be publicized in XGPU documentation that the
			//  redundant cross-pol is at xgpuidx-1."
			if(is_auto) {
				if(header->_ant_pol_prod_pol_index[idx] == 1) {
					// Store Ant-pol-product index of 'XY'
					redundant_pol_index[ant_1_idx] = idx;
				}
				else if (header->_ant_pol_prod_pol_index[idx] == 2) {
					// Derive indices for 'YX' from those for 'XY'
					header->_ant_pol_prod_xgpu_index[idx] = header->_ant_pol_prod_xgpu_index[redundant_pol_index[ant_1_idx]] - 1;
				}
			}

			UVH5print_verbose(__FUNCTION__, "#%d: (xgpu = %d, blidx = %d, polidx = %d, isauto = %d, needsconj = %d) [%d, %d]",
				idx,
				header->_ant_pol_prod_xgpu_index[idx],
				header->_ant_pol_prod_bl_index[idx],
				header->_ant_pol_prod_pol_index[idx],
				header->_ant_pol_prod_auto[idx],
				header->_ant_pol_prod_conj[idx],
				inpair_1, inpair_2
			);
			idx++;
		}
	}
	free(redundant_pol_index);

}

void UVH5permute_uvws(UVH5_header_t* header) {
	int ant_1_idx, ant_2_idx;
	for (int bls_idx = 0; bls_idx < header->Nbls; bls_idx++) {
		ant_1_idx = header->_antenna_num_idx_map[
			header->ant_1_array[bls_idx]
		];
		ant_2_idx = header->_antenna_num_idx_map[
			header->ant_2_array[bls_idx]
		];

		UVH5print_verbose(__FUNCTION__, "Baseline #%d\n(#%d @ %d) -> (#%d @ %d)",
			bls_idx,
			header->ant_1_array[bls_idx],
			ant_1_idx,
			header->ant_2_array[bls_idx],
			ant_2_idx
		);
		UVH5print_verbose(__FUNCTION__, "\t%s -> %s",
			header->antenna_names[ant_1_idx],
			header->antenna_names[ant_2_idx]
		);
		for (size_t i = 0; i < 3; i++)
		{
			header->uvw_array[bls_idx*3 + i] = // ant_1 -> ant_2
				header->_antenna_uvw_positions[ant_2_idx*3 + i] -
				header->_antenna_uvw_positions[ant_1_idx*3 + i]
			;
			UVH5print_verbose(__FUNCTION__, "\t\t %f - %f = %f",
				header->_antenna_uvw_positions[ant_2_idx*3 + i],
				header->_antenna_uvw_positions[ant_1_idx*3 + i],
				header->uvw_array[bls_idx*3 + i]
			);
		}
		UVH5print_verbose(__FUNCTION__, "\t\t\tDist: %f",
			calc_hypotenuse(header->uvw_array + bls_idx*3, 3)
		);
	}
}

int UVH5polarisation_string_key(char* pol_string, int npols) {
	if(npols == 1) {
		switch (pol_string[0]) {
		case 'i':
		case 'I': return 1;
		case 'q':
		case 'Q': return 2;
		case 'u':
		case 'U': return 3;
		case 'v':
		case 'V': return 4;
		default:
			break;
		}
	}
	else if(npols == 2) {
		switch (pol_string[0]) {
			case 'r':
			case 'R':
				switch(pol_string[1]) { 	
					case 'r':
					case 'R': return -1; // 'RR'
					case 'l':
					case 'L': return -3; // 'RL'
					default:
						break;
				}
				break;
			case 'l':
			case 'L':
				switch(pol_string[1]) { 	
					case 'l':
					case 'L': return -2; // 'LL'
					case 'r':
					case 'R': return -4; // 'LR'
					default:
						break;
				}
				break;
			case 'x':
			case 'X':
				switch(pol_string[1]) { 	
					case 'x':
					case 'X': return -5; // 'XX'
					case 'y':
					case 'Y': return -7; // 'XY'
					default:
						break;
				}
				break;
			case 'y':
			case 'Y':
				switch(pol_string[1]) { 	
					case 'y':
					case 'Y': return -6; // 'YY'
					case 'x':
					case 'X': return -8; // 'YX'
					default:
						break;
				}
				break;
			default:
				break;
		}
	}
	return 0;
}

/*
 * Emulates Rawx.jl `for (xgpuidx, blidx, polidx, isauto, needsconj) in inpair_maps`
 *
 * Critically reliant on `UVH5_header_t->_ant_pol_prod_*` administrative arrays, 
 * which are populated currently in `UVH5toml_parse_input_map`;
 */
void UVH5visdata_from_xgpu_double_output(
	UVH5_CF64_t* xgpuOutput, // [freq, xgpu_antpol_prod]
	UVH5_CF64_t* visdata, // [bl, freq, antpol_prod]
	size_t xgpuOutputElements,
	UVH5_header_t* header
) {
	const int ant_pol_products = header->Nbls*header->Npols;
	const int Nfreqs = header->Nfreqs;
	const int Npols = header->Npols;
	const int visdata_bl_stride = Nfreqs*header->Npols;
	const int xgpu_freq_stride = xgpuOutputElements/Nfreqs;

	int visdata_offset; // = blidx*visdata_bl_stride + freq*Npols + pol_idx
	int* xgpu_idx = header->_ant_pol_prod_xgpu_index;
	int* bl_idx = header->_ant_pol_prod_bl_index;
	int* pol_idx = header->_ant_pol_prod_pol_index;
	char* conjugate = header->_ant_pol_prod_conj;
	char* is_auto = header->_ant_pol_prod_auto;
	for (int approd_idx = 0; approd_idx < ant_pol_products; approd_idx++) {
		visdata_offset = (*bl_idx)*visdata_bl_stride + (*pol_idx);

		for (int freq = 0; freq < Nfreqs; freq++) {
			visdata[visdata_offset] = xgpuOutput[freq * xgpu_freq_stride + (*xgpu_idx)];
			if(*conjugate) {
				visdata[visdata_offset].i = -visdata[visdata_offset].i;
			}

			visdata_offset += Npols;
		}
		xgpu_idx++;
		bl_idx++;
		pol_idx++;
		conjugate++;
		is_auto++;
	}
}

/*
 * Emulates Rawx.jl `for (xgpuidx, blidx, polidx, isauto, needsconj) in inpair_maps`
 *
 * Critically reliant on `UVH5_header_t->_ant_pol_prod_*` administrative arrays, 
 * which are populated currently in `UVH5toml_parse_input_map`;
 */
void UVH5visdata_from_xgpu_float_output(
	UVH5_CF32_t* xgpuOutput, // [freq, xgpu_antpol_prod]
	UVH5_CF32_t* visdata, // [bl, freq, antpol_prod]
	size_t xgpuOutputElements,
	UVH5_header_t* header
) {
	const int ant_pol_products = header->Nbls*header->Npols;
	const int Nfreqs = header->Nfreqs;
	const int Npols = header->Npols;
	const int visdata_bl_stride = Nfreqs*header->Npols;
	const int xgpu_freq_stride = xgpuOutputElements/Nfreqs;

	int visdata_offset; // = blidx*visdata_bl_stride + freq*Npols + pol_idx
	int* xgpu_idx = header->_ant_pol_prod_xgpu_index;
	int* bl_idx = header->_ant_pol_prod_bl_index;
	int* pol_idx = header->_ant_pol_prod_pol_index;
	char* conjugate = header->_ant_pol_prod_conj;
	char* is_auto = header->_ant_pol_prod_auto;
	for (int approd_idx = 0; approd_idx < ant_pol_products; approd_idx++) {
		visdata_offset = (*bl_idx)*visdata_bl_stride + (*pol_idx);

		for (int freq = 0; freq < Nfreqs; freq++) {
			visdata[visdata_offset] = xgpuOutput[freq * xgpu_freq_stride + (*xgpu_idx)];
			if(*conjugate) {
				visdata[visdata_offset].i = -visdata[visdata_offset].i;
			}

			visdata_offset += Npols;
		}
		xgpu_idx++;
		bl_idx++;
		pol_idx++;
		conjugate++;
		is_auto++;
	}
}

/*
 * Emulates Rawx.jl `for (xgpuidx, blidx, polidx, isauto, needsconj) in inpair_maps`
 *
 * Critically reliant on `UVH5_header_t->_ant_pol_prod_*` administrative arrays, 
 * which are populated currently in `UVH5toml_parse_input_map`;
 */
void UVH5visdata_from_xgpu_int_output(
	UVH5_CI32_t* xgpuOutput, // [freq, xgpu_antpol_prod]
	UVH5_CI32_t* visdata, // [bl, freq, antpol_prod]
	size_t xgpuOutputElements,
	UVH5_header_t* header
) {
	const int ant_pol_products = header->Nbls*header->Npols;
	const int Nfreqs = header->Nfreqs;
	const int Npols = header->Npols;
	const int visdata_bl_stride = Nfreqs*header->Npols;
	const int xgpu_freq_stride = xgpuOutputElements/Nfreqs;

	int visdata_offset; // = blidx*visdata_bl_stride + freq*Npols + pol_idx
	int* xgpu_idx = header->_ant_pol_prod_xgpu_index;
	int* bl_idx = header->_ant_pol_prod_bl_index;
	int* pol_idx = header->_ant_pol_prod_pol_index;
	char* conjugate = header->_ant_pol_prod_conj;
	char* is_auto = header->_ant_pol_prod_auto;
	for (int approd_idx = 0; approd_idx < ant_pol_products; approd_idx++) {
		visdata_offset = (*bl_idx)*visdata_bl_stride + (*pol_idx);

		for (int freq = 0; freq < Nfreqs; freq++) {
			visdata[visdata_offset] = xgpuOutput[freq * xgpu_freq_stride + (*xgpu_idx)];
			if(*conjugate) {
				visdata[visdata_offset].i = -visdata[visdata_offset].i;
			}
			visdata_offset += Npols;
		}
		xgpu_idx++;
		bl_idx++;
		pol_idx++;
		conjugate++;
		is_auto++;
	}
}