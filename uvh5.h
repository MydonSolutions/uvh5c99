#ifndef UVH5_H
#define UVH5_H

#include "hdf5/serial/hdf5.h"
#include <string.h>
#include <stdlib.h>

#include "uvh5_bool_t.h"
#include "h5_dataspace.h"

typedef struct
{
	float latitude;						 /* The latitude of the telescope site, in degrees. */
	float longitude;					 /* The longitude of the telescope site, in degrees. */
	float altitude;						 /* The altitude of the telescope site, in meters. */
	char *telescope_name;			 /* The name of the telescope used to take the data. The value is used to check that metadata
																is self-consistent for known telescopes in pyuvdata. */
	char *instrument;					 /* The name of the instrument, typically the telescope name. */
	char *object_name;				 /* The name of the object tracked by the telescope. For a driftscan antenna, this is
																typically "zenith". */
	char *history;						 /* The history of the data file. */
	char *phase_type;					 /* The phase type of the observation. Should be "phased" or "drift". Note that "drift" in
																this context more accurately means "unphased", in that baselines are computing using ENU
																coordinates, without any w-projection. Any other value is treated as an unrecognized
																type. */
	int Nants_data;						 /* The number of antennas that data in the file corresponds to. May be smaller than the
																number of antennas in the array. */
	int Nants_telescope;			 /* The number of antennas in the array. May be larger than the number of antennas with
																data corresponding to them. */
	int *ant_1_array;					 /* An array of the first antenna numbers corresponding to baselines present in the data.
																All entries in this array must exist in the antenna numbers array.
																This is a one-dimensional array of size [Nbls] (reiteratively written Ntimes, so that
																stored it has size [Nblts]). */
	int *ant_2_array;					 /* An array of the second antenna numbers corresponding to baselines
																present in the data. All entries in this array must exist in the antenna numbers array.
																This is a one-dimensional array of size [Nbls] (reiteratively written Ntimes, so that
																stored it has size [Nblts]). */
	int *antenna_numbers;			 /* An array of the numbers of the antennas present in the
																radio telescope (note that these are not indices, they do not need to start at zero or
																be continuous). This is a one-dimensional array of size Nants_telescope. Note there
																must be one entry for every unique antenna in ant 1 array and ant 2 array, but there
																may be additional entries. */
	char **antenna_names;			 /* An array of the names of antennas present in the radio telescope. This is a
																one-dimensional array of size Nants_telescope. Note there must be one entry for every
																unique antenna in ant 1 array and ant 2 array, but there may be additional entries. */
	int Nbls;									 /* The number of baselines present in the data. For full cross-correlation data (including
																auto-correlations), this should be Nants_data×(Nants_data+1)/2. */
	int Nblts;								 /* The number of baseline-times (i.e., the number of spectra) present in the data. Note
																that this value need not be equal to Nbls × Ntimes. */
	int Nspws;								 /* The number of spectral windows present in the data. */
	int Nfreqs;								 /* The total number of frequency channels in the data across all spectral windows. */
	int Npols;								 /* The number of polarization products in the data. */
	int Ntimes;								 /* The number of time samples present in the data. */
	float *uvw_array;					 /* An array of the uvw-coordinates corresponding to each observation in the data. This is a
																two-dimensional array of size [Nbls, 3] (reiteratively written Ntimes, so that stored it
																has size [Nblts, 3]). Units are in meters. */
	float *time_array;				 /* An array of the Julian Date corresponding to the temporal midpoint of the corresponding 
																baseline’s integration. This is a one-dimensional array of size [Nbls] (reiteratively
																written Ntimes, so that stored it has size [Nblts]). */
	float *integration_time;	 /* An array of the duration in seconds of an integration. This is a one-dimensional array
																of size [Nbls] (reiteratively written Ntimes, so that stored it has size [Nblts]). */
	float *freq_array;				 /* An array of all the frequencies (for all spectral windows) stored in the file in Hertz.
																This is a one-dimensional array of size [Nfreqs]. */
	float *channel_width;			 /* The width of frequency channels in the file in Hertz. This is a one-dimensional array of
																size [Nfreqs]. */
	int *spw_array;						 /* An array of the spectral windows in the file. This is a one-dimensional array of size
																[Nspws]. */
	UVH5_bool_t flex_spw;			 /* Whether the data are saved using flexible spectral windows. If more than one spectral
																window is present in the data, this must be True. */
	int *polarization_array;	 /* An array of the polarizations contained in the file. This is a one-dimensional array of
																size [Npols]. Note that the polarizations should be stored as an integer, and use the
																convention defined in AIPS Memo 117. */
	float *antenna_positions;	 /* An array of the antenna coordinates relative to the reference position of the radio
																telescope array, which is implicitly defined by the latitude, longitude, and altitude
																(LLA) parameters. More explicitly, these are the ECEF coordinates of individual antennas
																minus the ECEF coordinates of the reference telescope position, such that the telescope
																position plus the values stored in antenna positions equals the position of individual
																elements in ECEF. The conversion between LLA and ECEF is given by WGS84. This is a
																two-dimensional array of size [Nants_telescope, 3]. */
	// Optional entries follow
	int* flex_spw_id_array;		 /* The mapping of individual channels along the frequency axis
																to individual spectral windows, as listed in the spw array. This is a one-dimensional
																array of size [Nfreqs]. Note this is required if the file uses flexible spectral windows
																(see Sec. 5.2.2).*/
	float dut1;								 /* DUT1 (google it), in seconds. AIPS 117 calls it UT1UTC.*/
	float earth_omega;				 /* Earth's rotation rate in degrees per day. Note the dif`ference in
																units, which is inherited from the way this quantity is handled in UVFITS datasets
																(AIPS 117 calls it DEGPDY).*/
	float gst0;								 /* Greenwich sidereal time at midnight on reference date, in degrees. AIPS
																117 calls it GSTIA0. */
	char *rdate;							 /* Date for which GST0 (or whichever time saved in that field) applies.
																Note this is different from how UVFITS handles this quantity, which is saved as
																a float rather than a string. The user is encouraged to ensure it is being handled
																self-consistently for their desired application. */
	char* timesys;						 /* Time system. pyuvdata currently only supports UTC. */
	char* x_orientation;			 /* The orientation of the x-arm of a dipole antenna. It is assumed
																to be the same for all antennas in the dataset. For instance, "East" or "North" may
																be used. */
	float* antenna_diameters;	 /* An array of the diameters of the antennas in meters.
																This is a one-dimensional array of size [Nants_telescope]. */
	int uvplane_reference_time;/* The time at which the phase center is normal to the
																chosen UV plane for phasing. Used for interoperability with the FHD package. */
	float phase_center_ra;		 /* The right ascension of the phase center of the observation
																in radians. Required if phase type is "phased". */
	float phase_center_dec;		 /* The declination of the phase center of the observation in
																radians. Required if phase type is "phased". */
	float phase_center_epoch;	 /* The epoch year of the phase applied to the data (e.g.,
																2000.). Required if phase type is "phased". (phase center epoch)
																string phase center frame The frame the data and uvw array are phased to.
																Options are "gcrs" and "icrs", with default "icrs". These frames are defined as
																coordinate systems in astropy. */
	float* lst_array;					 /* An array corresponding to the local sidereal time of the center of
																each observation in the data in units of radians. If it is not specified, it is calculated
																from the latitude/longitude and the time array. Saving it in the file can be useful
																for files with many values in the time array, which would expensive to recompute. */
	
	// Administrative entries follow
	int* _antenna_num_idx_map; /* An array whose elements translate the corresponding index from antenna-number to index
																in the other antenna_* arrays. */

} UVH5_header_t;

static inline void UVH5Halloc(UVH5_header_t *header)
{
	if(header->ant_1_array == NULL && header->Nbls != 0) {
		header->ant_1_array = malloc(sizeof(int) * header->Nbls);
		memset(header->ant_1_array, 0, sizeof(int) * header->Nbls);
		fprintf(stderr, "UVH5: 'ant_1_array' allocated %ld bytes.\n", sizeof(int) * header->Nbls);
	}
	if(header->ant_2_array == NULL && header->Nbls != 0) {
		header->ant_2_array = malloc(sizeof(int) * header->Nbls);
		memset(header->ant_2_array, 0, sizeof(int) * header->Nbls);
		fprintf(stderr, "UVH5: 'ant_2_array' allocated %ld bytes.\n", sizeof(int) * header->Nbls);
	}
	if(header->antenna_numbers == NULL && header->Nants_telescope != 0) {
		header->antenna_numbers = malloc(sizeof(int) * header->Nants_telescope);
		memset(header->antenna_numbers, 0, sizeof(int) * header->Nants_telescope);
		fprintf(stderr, "UVH5: 'antenna_numbers' allocated %ld bytes.\n", sizeof(int) * header->Nants_telescope);
	}
	if(header->antenna_names == NULL && header->Nants_telescope != 0) {
		header->antenna_names = malloc(sizeof(char *) * header->Nants_telescope);
		memset(header->antenna_names, 0, sizeof(char *) * header->Nants_telescope);
		fprintf(stderr, "UVH5: 'antenna_names' allocated %ld bytes.\n", sizeof(char *) * header->Nants_telescope);
	}
	if(header->uvw_array == NULL && header->Nbls * 3 != 0) {
		header->uvw_array = malloc(sizeof(float) * header->Nbls * 3);
		memset(header->uvw_array, 0, sizeof(float) * header->Nbls * 3);
		fprintf(stderr, "UVH5: 'uvw_array' allocated %ld bytes.\n", sizeof(float) * header->Nbls * 3);
	}
	if(header->time_array == NULL && header->Nbls != 0) {
		header->time_array = malloc(sizeof(float) * header->Nbls);
		memset(header->time_array, 0, sizeof(float) * header->Nbls);
		fprintf(stderr, "UVH5: 'time_array' allocated %ld bytes.\n", sizeof(float) * header->Nbls);
	}
	if(header->integration_time == NULL && header->Nbls != 0) {
		header->integration_time = malloc(sizeof(float) * header->Nbls);
		memset(header->integration_time, 0, sizeof(float) * header->Nbls);
		fprintf(stderr, "UVH5: 'integration_time' allocated %ld bytes.\n", sizeof(float) * header->Nbls);
	}
	if(header->freq_array == NULL && header->Nfreqs != 0) {
		header->freq_array = malloc(sizeof(float) * header->Nfreqs);
		memset(header->freq_array, 0, sizeof(float) * header->Nfreqs);
		fprintf(stderr, "UVH5: 'freq_array' allocated %ld bytes.\n", sizeof(float) * header->Nfreqs);
	}
	if(header->channel_width == NULL && header->Nfreqs != 0) {
		header->channel_width = malloc(sizeof(float) * header->Nfreqs);
		memset(header->channel_width, 0, sizeof(float) * header->Nfreqs);
		fprintf(stderr, "UVH5: 'channel_width' allocated %ld bytes.\n", sizeof(float) * header->Nfreqs);
	}
	if(header->spw_array == NULL && header->Nspws != 0) {
		header->spw_array = malloc(sizeof(float) * header->Nspws);
		memset(header->spw_array, 0, sizeof(float) * header->Nspws);
		fprintf(stderr, "UVH5: 'spw_array' allocated %ld bytes.\n", sizeof(float) * header->Nspws);
	}
	if(header->polarization_array == NULL && header->Npols != 0) {
		header->polarization_array = malloc(sizeof(int) * header->Npols);
		memset(header->polarization_array, 0, sizeof(int) * header->Npols);
		fprintf(stderr, "UVH5: 'polarization_array' allocated %ld bytes.\n", sizeof(int) * header->Npols);
	}
	if(header->antenna_positions == NULL && header->Nants_telescope * 3 != 0) {
		header->antenna_positions = malloc(sizeof(float) * header->Nants_telescope * 3);
		memset(header->antenna_positions, 0, sizeof(float) * header->Nants_telescope * 3);
		fprintf(stderr, "UVH5: 'antenna_positions' allocated %ld bytes.\n", sizeof(float) * header->Nants_telescope * 3);
	}
	// Optional entries follow
	// if(header-> == NULL && header->Nbls != 0) {
	// 	header->lst_array = malloc(sizeof(float) * header->Nbls);
	// 	memset(header->lst_array, 0, sizeof(float) * header->Nbls);
	// 	fprintf(stderr, "UVH5: 'lst_array' allocated %ld bytes.\n", sizeof(float) * header->Nbls);
	// }
}

static inline void UVH5Hfree(UVH5_header_t *header)
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
	// Optional arrays follow
	if (header->antenna_diameters != NULL) {
		free(header->antenna_diameters);
	}
	// Administrative arrays follow
	free(header->_antenna_num_idx_map);
}

typedef struct
{
	hid_t file_id;
	hid_t header_id;
	UVH5_header_t header;
	H5_open_dataspace_t DS_header_Ntimes;
	H5_open_dataspace_t DS_header_Nblts;
	H5_open_dataspace_t DS_header_ant_1_array;
	H5_open_dataspace_t DS_header_ant_2_array;
	H5_open_dataspace_t DS_header_uvw_array;
	H5_open_dataspace_t DS_header_time_array;
	H5_open_dataspace_t DS_header_integration_time;
	H5_open_dataspace_t DS_header_lst_array;
	hid_t data_id;
	void *visdata;
	UVH5_bool_t *flags;
	float *nsamples;
	H5_open_dataspace_t DS_data_visdata;
	H5_open_dataspace_t DS_data_flags;
	H5_open_dataspace_t DS_data_nsamples;
} UVH5_file_t;

typedef struct __attribute__((__packed__))
{
	int32_t r;
	int32_t i;
} UVH5_CI32_t;
static hid_t UVH5TcreateCI32()
{
	hid_t UVH5_CI32 = H5Tcreate(H5T_COMPOUND, sizeof(UVH5_CI32_t));
	H5Tinsert(UVH5_CI32, "r", HOFFSET(UVH5_CI32_t, r), H5T_NATIVE_INT32);
	H5Tinsert(UVH5_CI32, "i", HOFFSET(UVH5_CI32_t, i), H5T_NATIVE_INT32);
	return UVH5_CI32;
}

typedef struct __attribute__((__packed__))
{
	float r;
	float i;
} UVH5_CF32_t;
static hid_t UVH5TcreateCF32()
{
	hid_t UVH5_CF32 = H5Tcreate(H5T_COMPOUND, sizeof(UVH5_CF32_t));
	H5Tinsert(UVH5_CF32, "r", HOFFSET(UVH5_CF32_t, r), H5T_NATIVE_FLOAT);
	H5Tinsert(UVH5_CF32, "i", HOFFSET(UVH5_CF32_t, i), H5T_NATIVE_FLOAT);
	return UVH5_CF32;
}

typedef struct __attribute__((__packed__))
{
	double r;
	double i;
} UVH5_CF64_t;
static hid_t UVH5TcreateCF64()
{
	hid_t UVH5_CF64 = H5Tcreate(H5T_COMPOUND, sizeof(UVH5_CF64_t));
	H5Tinsert(UVH5_CF64, "r", HOFFSET(UVH5_CF64_t, r), H5T_NATIVE_DOUBLE);
	H5Tinsert(UVH5_CF64, "i", HOFFSET(UVH5_CF64_t, i), H5T_NATIVE_DOUBLE);
	return UVH5_CF64;
}

static inline herr_t _H5DboolWrite(
		hid_t dest_id, char *d_name, const int rank, const hsize_t *dims, const UVH5_bool_t *data)
{
	H5_open_dataspace_t dataspace = {0};
	dataspace.name = d_name;
	H5DSset(rank, dims, NULL, &dataspace);
	if (rank > 0)
	{
		dataspace.filter_flag = 3;
	}
	H5DSopenBool(dest_id, &dataspace);

	herr_t status = H5DSwrite(&dataspace, data);
	status += H5DSclose(&dataspace);
	return status;
}

static inline herr_t _H5DdoubleWrite(
		hid_t dest_id, char *d_name, const int rank, const hsize_t *dims, const double *data)
{
	H5_open_dataspace_t dataspace = {0};
	dataspace.name = d_name;
	H5DSset(rank, dims, NULL, &dataspace);
	if (rank > 0)
	{
		dataspace.filter_flag = 3;
	}
	H5DSopenDouble(dest_id, &dataspace);

	herr_t status = H5DSwrite(&dataspace, data);
	status += H5DSclose(&dataspace);
	return status;
}

static inline herr_t _H5DfloatWrite(
		hid_t dest_id, char *d_name, const int rank, const hsize_t *dims, const float *data)
{
	H5_open_dataspace_t dataspace = {0};
	dataspace.name = d_name;
	H5DSset(rank, dims, NULL, &dataspace);
	if (rank > 0)
	{
		dataspace.filter_flag = 3;
	}
	H5DSopenFloat(dest_id, &dataspace);

	herr_t status = H5DSwrite(&dataspace, data);
	status += H5DSclose(&dataspace);
	return status;
}

static inline herr_t _H5DintWrite(
		hid_t dest_id, char *d_name, const int rank, const hsize_t *dims, const int *data)
{
	H5_open_dataspace_t dataspace = {0};
	dataspace.name = d_name;
	H5DSset(rank, dims, NULL, &dataspace);
	if (rank > 0)
	{
		dataspace.filter_flag = 3;
	}
	H5DSopenInt(dest_id, &dataspace);

	herr_t status = H5DSwrite(&dataspace, data);
	status += H5DSclose(&dataspace);
	return status;
}

static inline herr_t _H5DstringWrite(
		hid_t dest_id, char *d_name, const int rank, const hsize_t *dims, const char *data)
{
	H5_open_dataspace_t dataspace = {0};
	dataspace.name = d_name;
	H5DSset(rank, dims, NULL, &dataspace);

	hid_t T_id = H5Tcopy(H5T_C_S1);
	H5Tset_size(T_id, rank > 0 ? H5S_UNLIMITED : strlen(data));
	H5DSopen(dest_id, T_id, H5Tcopy(T_id), &dataspace);

	herr_t status = H5DSwrite(&dataspace, data);
	status += H5DSclose(&dataspace);
	return status;
}

static inline void _UVH5_Hwrite_static(UVH5_file_t *uvh5_file)
{
	// TODO return error code
	herr_t status;

	UVH5_header_t header = uvh5_file->header;
	const hsize_t dims1_Nants_telescope[] = {header.Nants_telescope};
	const hsize_t dims1_Nspws[] = {header.Nspws};
	const hsize_t dims1_Nfreqs[] = {header.Nfreqs};
	const hsize_t dims1_Npols[] = {header.Npols};
	const hsize_t dims2_Nants_telescope_3[] = {header.Nants_telescope, 3};

	// TODO assert uvh5_file->header_id

	status = _H5DfloatWrite(uvh5_file->header_id, "latitude", 0, NULL, &header.latitude);
	if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on 'latitude'\n"); return; }

	status = _H5DfloatWrite(uvh5_file->header_id, "longitude", 0, NULL, &header.longitude);
	if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on 'longitude'\n"); return; }

	status = _H5DfloatWrite(uvh5_file->header_id, "altitude", 0, NULL, &header.altitude);
	if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on 'altitude'\n"); return; }

	status = _H5DstringWrite(uvh5_file->header_id, "telescope_name", 0, NULL, header.telescope_name);
	if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on 'telescope_name'\n"); return; }

	status = _H5DstringWrite(uvh5_file->header_id, "instrument", 0, NULL, header.instrument);
	if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on 'instrument'\n"); return; }

	status = _H5DstringWrite(uvh5_file->header_id, "object_name", 0, NULL, header.object_name);
	if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on 'object_name'\n"); return; }

	status = _H5DstringWrite(uvh5_file->header_id, "history", 0, NULL, header.history);
	if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on 'history'\n"); return; }

	status = _H5DstringWrite(uvh5_file->header_id, "phase_type", 0, NULL, header.phase_type);
	if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on 'phase_type'\n"); return; }

	status = _H5DintWrite(uvh5_file->header_id, "Nants_telescope", 0, NULL, &header.Nants_telescope);
	if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on 'Nants_telescope'\n"); return; }

	status = _H5DstringWrite(uvh5_file->header_id, "antenna_names", 1, dims1_Nants_telescope, (char *)header.antenna_names);
	if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on 'antenna_names'\n"); return; }

	status = _H5DintWrite(uvh5_file->header_id, "antenna_numbers", 1, dims1_Nants_telescope, header.antenna_numbers);
	if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on 'antenna_numbers'\n"); return; }

	status = _H5DintWrite(uvh5_file->header_id, "Nbls", 0, NULL, &header.Nbls);
	if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on 'Nbls'\n"); return; }

	status = _H5DintWrite(uvh5_file->header_id, "Nspws", 0, NULL, &header.Nspws);
	if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on 'Nspws'\n"); return; }

	status = _H5DintWrite(uvh5_file->header_id, "Nfreqs", 0, NULL, &header.Nfreqs);
	if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on 'Nfreqs'\n"); return; }

	status = _H5DintWrite(uvh5_file->header_id, "Npols", 0, NULL, &header.Npols);
	if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on 'Npols'\n"); return; }

	status = _H5DfloatWrite(uvh5_file->header_id, "freq_array", 1, dims1_Nfreqs, header.freq_array);
	if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on 'freq_array'\n"); return; }

	status = _H5DfloatWrite(uvh5_file->header_id, "channel_width", 1, dims1_Nfreqs, header.channel_width);
	if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on 'channel_width'\n"); return; }

	status = _H5DintWrite(uvh5_file->header_id, "spw_array", 1, dims1_Nspws, header.spw_array);
	if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on 'spw_array'\n"); return; }

	status = _H5DboolWrite(uvh5_file->header_id, "flex_spw", 0, NULL, &header.flex_spw);
	if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on 'flex_spw'\n"); return; }

	status = _H5DintWrite(uvh5_file->header_id, "polarization_array", 1, dims1_Npols, header.polarization_array);
	if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on 'polarization_array'\n"); return; }

	status = _H5DfloatWrite(uvh5_file->header_id, "antenna_positions", 2, dims2_Nants_telescope_3, header.antenna_positions);
	if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on 'antenna_positions'\n"); return; }

	status = _H5DstringWrite(uvh5_file->header_id, "version", 0, NULL, "1.0");
	if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on 'version'\n"); return; }

	// Optional entries follow
	if(header.flex_spw_id_array){
		status = _H5DintWrite(uvh5_file->header_id, "flex_spw_id_array", 1, dims1_Nfreqs, header.flex_spw_id_array);
		if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on optional 'flex_spw_id_array'\n"); }
	}

	if(header.dut1){
		status = _H5DfloatWrite(uvh5_file->header_id, "dut1", 0, NULL, &header.dut1);
		if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on optional 'dut1'\n"); }
	}

	if(header.earth_omega){
		status = _H5DfloatWrite(uvh5_file->header_id, "earth_omega", 0, NULL, &header.earth_omega);
		if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on optional 'earth_omega'\n"); }
	}

	if(header.gst0){
		status = _H5DfloatWrite(uvh5_file->header_id, "gst0", 0, NULL, &header.gst0);
		if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on optional 'gst0'\n"); }
	}

	if(header.rdate){
		status = _H5DstringWrite(uvh5_file->header_id, "rdate", 0, NULL, header.rdate);
		if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on optional 'rdate'\n"); }
	}

	if(header.timesys){
		status = _H5DstringWrite(uvh5_file->header_id, "timesys", 0, NULL, header.timesys);
		if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on optional 'timesys'\n"); }
	}

	if(header.x_orientation){
		status = _H5DstringWrite(uvh5_file->header_id, "x_orientation", 0, NULL, header.x_orientation);
		if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on optional 'x_orientation'\n"); }
	}

	if(header.antenna_diameters){
		status = _H5DfloatWrite(uvh5_file->header_id, "antenna_diameters", 1, dims1_Nants_telescope, header.antenna_diameters);
		if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on optional 'antenna_diameters'\n"); }
	}

	if(header.uvplane_reference_time){
		status = _H5DintWrite(uvh5_file->header_id, "uvplane_reference_time", 0, NULL, &header.uvplane_reference_time);
		if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on optional 'uvplane_reference_time'\n"); }
	}

	if(header.phase_center_ra){
		status = _H5DfloatWrite(uvh5_file->header_id, "phase_center_ra", 0, NULL, &header.phase_center_ra);
		if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on optional 'phase_center_ra'\n"); }
	}

	if(header.phase_center_dec){
		status = _H5DfloatWrite(uvh5_file->header_id, "phase_center_dec", 0, NULL, &header.phase_center_dec);
		if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on optional 'phase_center_dec'\n"); }
	}

	if(header.phase_center_epoch){
		status = _H5DfloatWrite(uvh5_file->header_id, "phase_center_epoch", 0, NULL, &header.phase_center_epoch);
		if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on optional 'phase_center_epoch'\n"); }
	}

}

static inline void UVH5open(char* filepath, UVH5_file_t *uvh5_file, hid_t Tvisdata)
{
	
	uvh5_file->file_id = H5Fcreate(filepath, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	uvh5_file->DS_data_visdata.Tmem_id = Tvisdata;
	uvh5_file->DS_data_visdata.Tsto_id = H5Tcopy(Tvisdata);
	 
	const hsize_t dim1_unlim[] = {H5S_UNLIMITED};
	const hsize_t dim1_nbls[] = {uvh5_file->header.Nbls};
	const hsize_t dim2_unlim_3[] = {H5S_UNLIMITED, 3};
	const hsize_t dim2_nbls_3[] = {uvh5_file->header.Nbls, 3};

	const hsize_t dim3_data_lim[] = {H5S_UNLIMITED, uvh5_file->header.Nfreqs, uvh5_file->header.Npols};
	const hsize_t dim3_data_chunk[] = {uvh5_file->header.Nbls, uvh5_file->header.Nfreqs, uvh5_file->header.Npols};

	// 'Header' group
	uvh5_file->header_id = H5Gcreate(uvh5_file->file_id, "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	uvh5_file->DS_header_Ntimes.name = "Ntimes";
	H5DSset(0, NULL, NULL, &uvh5_file->DS_header_Ntimes);
	H5DSopenInt(uvh5_file->header_id, &uvh5_file->DS_header_Ntimes);

	uvh5_file->DS_header_Nblts.name = "Nblts";
	H5DSset(0, NULL, NULL, &uvh5_file->DS_header_Nblts);
	H5DSopenInt(uvh5_file->header_id, &uvh5_file->DS_header_Nblts);

	uvh5_file->DS_header_ant_1_array.name = "ant_1_array";
	H5DSset(1, dim1_unlim, dim1_nbls, &uvh5_file->DS_header_ant_1_array);
	H5DSopenInt(uvh5_file->header_id, &uvh5_file->DS_header_ant_1_array);

	uvh5_file->DS_header_ant_2_array.name = "ant_2_array";
	H5DSset(1, dim1_unlim, dim1_nbls, &uvh5_file->DS_header_ant_2_array);
	H5DSopenInt(uvh5_file->header_id, &uvh5_file->DS_header_ant_2_array);

	uvh5_file->DS_header_uvw_array.name = "uvw_array";
	H5DSset(2, dim2_unlim_3, dim2_nbls_3, &uvh5_file->DS_header_uvw_array);
	H5DSopenFloat(uvh5_file->header_id, &uvh5_file->DS_header_uvw_array);

	uvh5_file->DS_header_time_array.name = "time_array";
	H5DSset(1, dim1_unlim, dim1_nbls, &uvh5_file->DS_header_time_array);
	H5DSopenFloat(uvh5_file->header_id, &uvh5_file->DS_header_time_array);

	uvh5_file->DS_header_integration_time.name = "integration_time";
	H5DSset(1, dim1_unlim, dim1_nbls, &uvh5_file->DS_header_integration_time);
	H5DSopenFloat(uvh5_file->header_id, &uvh5_file->DS_header_integration_time);

	if(uvh5_file->header.lst_array){
		uvh5_file->DS_header_lst_array.name = "lst_array";
		H5DSset(1, dim1_unlim, dim1_nbls, &uvh5_file->DS_header_lst_array);
		H5DSopenFloat(uvh5_file->header_id, &uvh5_file->DS_header_lst_array);
	}

	_UVH5_Hwrite_static(uvh5_file);
	int zero = 0;
	H5DSwrite(&uvh5_file->DS_header_Ntimes, &zero);
	H5DSwrite(&uvh5_file->DS_header_Nblts, &zero);

	// 'Data' group
	uvh5_file->data_id = H5Gcreate(uvh5_file->file_id, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	uvh5_file->DS_data_visdata.name = "visdata";
	H5DSset(3, dim3_data_lim, dim3_data_chunk, &uvh5_file->DS_data_visdata);
	// zeros for Tmem/sto_id respect DS_data_visdata internal values
	H5DSopen(uvh5_file->data_id, 0, 0, &uvh5_file->DS_data_visdata);
	size_t data_visdata_nbytes = H5DSnelem_chunks(&uvh5_file->DS_data_visdata)*H5Tget_size(uvh5_file->DS_data_visdata.Tmem_id);
	uvh5_file->visdata = malloc(data_visdata_nbytes);
	memset(uvh5_file->visdata, 0, data_visdata_nbytes);
	fprintf(stderr, "UVH5: 'data_visdata' allocated %ld bytes.\n", data_visdata_nbytes);

	uvh5_file->DS_data_flags.name = "flags";
	H5DSset(3, dim3_data_lim, dim3_data_chunk, &uvh5_file->DS_data_flags);
	uvh5_file->DS_data_flags.filter_flag = H5_FILTER_FLAG_DEFLATE_3;
	H5DSopenBool(uvh5_file->data_id, &uvh5_file->DS_data_flags);
	size_t data_flags_nbytes = H5DSnelem_chunks(&uvh5_file->DS_data_flags)*H5Tget_size(uvh5_file->DS_data_flags.Tmem_id);
	uvh5_file->flags = malloc(data_flags_nbytes);
	memset(uvh5_file->flags, 0, data_flags_nbytes);
	fprintf(stderr, "UVH5: 'data_flags' allocated %ld bytes.\n", data_flags_nbytes);

	uvh5_file->DS_data_nsamples.name = "nsamples";
	H5DSset(3, dim3_data_lim, dim3_data_chunk, &uvh5_file->DS_data_nsamples);
	uvh5_file->DS_data_flags.filter_flag = H5_FILTER_FLAG_DEFLATE_3;
	H5DSopenFloat(uvh5_file->data_id, &uvh5_file->DS_data_nsamples);
	size_t data_nsamples_nbytes = H5DSnelem_chunks(&uvh5_file->DS_data_nsamples)*H5Tget_size(uvh5_file->DS_data_nsamples.Tmem_id);
	uvh5_file->nsamples = malloc(data_nsamples_nbytes);
	memset(uvh5_file->nsamples, 0, data_nsamples_nbytes);
	fprintf(stderr, "UVH5: 'data_nsamples' allocated %ld bytes.\n", data_nsamples_nbytes);
}

static inline void UVH5close(UVH5_file_t *uvh5_file)
{
	// 'Header' group
	H5DSclose(&uvh5_file->DS_header_Ntimes);
	H5DSclose(&uvh5_file->DS_header_Nblts);
	H5DSclose(&uvh5_file->DS_header_ant_1_array);
	H5DSclose(&uvh5_file->DS_header_ant_2_array);
	H5DSclose(&uvh5_file->DS_header_uvw_array);
	H5DSclose(&uvh5_file->DS_header_time_array);
	H5DSclose(&uvh5_file->DS_header_integration_time);

	if(uvh5_file->header.lst_array){
		H5DSclose(&uvh5_file->DS_header_lst_array);
	}
	H5Gclose(uvh5_file->header_id);
	UVH5Hfree(&uvh5_file->header);

	// 'Data' group
	H5DSclose(&uvh5_file->DS_data_visdata);
	H5DSclose(&uvh5_file->DS_data_flags);
	H5DSclose(&uvh5_file->DS_data_nsamples);

	H5Gclose(uvh5_file->data_id);
	free(uvh5_file->visdata);
	free(uvh5_file->flags);
	free(uvh5_file->nsamples);
	
	H5Fclose(uvh5_file->file_id);
}

int UVH5write_dynamic(UVH5_file_t* uvh5_file) {
	herr_t status;
	// Header
	status = H5DSextend(&uvh5_file->DS_header_ant_1_array);
	if (status < 0) { fprintf(stderr, "UVH5: H5DSextend failure on 'header_ant_1_array'\n"); return -1;}
	status = H5DSwrite(&uvh5_file->DS_header_ant_1_array, uvh5_file->header.ant_1_array);

	if (status < 0) { fprintf(stderr, "UVH5: H5DSwrite failure on 'header_ant_1_array'\n"); return -1;}
	status = H5DSextend(&uvh5_file->DS_header_ant_2_array);
	if (status < 0) { fprintf(stderr, "UVH5: H5DSextend failure on 'header_ant_2_array'\n"); return -1;}
	status = H5DSwrite(&uvh5_file->DS_header_ant_2_array, uvh5_file->header.ant_2_array);

	if (status < 0) { fprintf(stderr, "UVH5: H5DSwrite failure on 'header_ant_2_array'\n"); return -1;}
	status = H5DSextend(&uvh5_file->DS_header_uvw_array);
	if (status < 0) { fprintf(stderr, "UVH5: H5DSextend failure on 'header_uvw_array'\n"); return -1;}
	status = H5DSwrite(&uvh5_file->DS_header_uvw_array, uvh5_file->header.uvw_array);

	if (status < 0) { fprintf(stderr, "UVH5: H5DSwrite failure on 'header_uvw_array'\n"); return -1;}
	status = H5DSextend(&uvh5_file->DS_header_time_array);
	if (status < 0) { fprintf(stderr, "UVH5: H5DSextend failure on 'header_time_array'\n"); return -1;}
	status = H5DSwrite(&uvh5_file->DS_header_time_array, uvh5_file->header.time_array);

	if (status < 0) { fprintf(stderr, "UVH5: H5DSwrite failure on 'header_time_array'\n"); return -1;}
	status = H5DSextend(&uvh5_file->DS_header_integration_time);
	if (status < 0) { fprintf(stderr, "UVH5: H5DSextend failure on 'header_integration_time'\n"); return -1;}
	status = H5DSwrite(&uvh5_file->DS_header_integration_time, uvh5_file->header.integration_time);
	if (status < 0) { fprintf(stderr, "UVH5: H5DSwrite failure on 'header_integration_time'\n"); return -1;}

	// Data
	status = H5DSextend(&uvh5_file->DS_data_visdata);
	if (status < 0) { fprintf(stderr, "UVH5: H5DSextend failure on 'data_visdata'\n"); return -1;}
	status = H5DSwrite(&uvh5_file->DS_data_visdata, uvh5_file->visdata);

	if (status < 0) { fprintf(stderr, "UVH5: H5DSwrite failure on 'data_visdata'\n"); return -1;}
	status = H5DSextend(&uvh5_file->DS_data_flags);
	if (status < 0) { fprintf(stderr, "UVH5: H5DSextend failure on 'data_flags'\n"); return -1;}
	status = H5DSwrite(&uvh5_file->DS_data_flags, uvh5_file->flags);

	if (status < 0) { fprintf(stderr, "UVH5: H5DSwrite failure on 'data_flags'\n"); return -1;}
	status = H5DSextend(&uvh5_file->DS_data_nsamples);
	if (status < 0) { fprintf(stderr, "UVH5: H5DSextend failure on 'data_nsamples'\n"); return -1;}
	status = H5DSwrite(&uvh5_file->DS_data_nsamples, uvh5_file->nsamples);
	if (status < 0) { fprintf(stderr, "UVH5: H5DSwrite failure on 'data_nsamples'\n"); return -1;}

	// Update static
	uvh5_file->header.Ntimes += 1;
	uvh5_file->header.Nblts += uvh5_file->header.Nbls;
	status = H5DSwrite(&uvh5_file->DS_header_Ntimes, &uvh5_file->header.Ntimes);
	if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on 'Ntimes'\n"); return -1; }

	status = H5DSwrite(&uvh5_file->DS_header_Nblts, &uvh5_file->header.Nblts);
	if (status < 0)	{ fprintf(stderr, "UVH5 Header: failure on 'Nblts'\n"); return -1; }

	return 0;
}

static int find_antenna_index_by_name(UVH5_header_t* header, char* name) {
	if(header->antenna_names == NULL) {
			return -2;
	}
	for(size_t i = 0; i < header->Nants_telescope; i++) {
		if(strcmp(name, header->antenna_names[i]) == 0) { // Case sensitive
			fprintf(stderr, "Antenna '%s' is at index %ld.\n", name, i);
			return i;
		}
	}
	return -1;
}

#endif