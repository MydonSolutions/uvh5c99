#ifndef UVH5_H
#define UVH5_H

#include "hdf5.h"
#include <string.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "h5dsc99/h5_bool.h"
#include "h5dsc99/h5_dataspace.h"
#include "radiointerferometryc99.h"
#include "uvh5/uvh5_printer.h"
#include "uvh5/uvh5_phase_center.h"

typedef struct
{
    double latitude;							/* The latitude of the telescope site, in degrees. */
    double longitude;							/* The longitude of the telescope site, in degrees. */
    double altitude;							/* The altitude of the telescope site, in meters. */
    char *telescope_name;				 		/* The name of the telescope used to take the data. The value is used to check that metadata
                                                    is self-consistent for known telescopes in pyuvdata. */
    char *instrument;						 	/* The name of the instrument, typically the telescope name. */
    char *history;							 	/* The history of the data file. */
    int Nants_data;							 	/* The number of antennas that data in the file corresponds to. May be smaller than the
                                                    number of antennas in the array. */
    int Nants_telescope;				 		/* The number of antennas in the array. May be larger than the number of antennas with
                                                    data corresponding to them. */
    int *ant_1_array;						 	/* An array of the first antenna numbers corresponding to baselines present in the data.
                                                    All entries in this array must exist in the antenna numbers array.
                                                    This is a one-dimensional array of size [Nbls] (reiteratively written Ntimes, so that
                                                    stored it has size [Nblts]). */
    int *ant_2_array;						 	/* An array of the second antenna numbers corresponding to baselines
                                                    present in the data. All entries in this array must exist in the antenna numbers array.
                                                    This is a one-dimensional array of size [Nbls] (reiteratively written Ntimes, so that
                                                    stored it has size [Nblts]). */
    int *antenna_numbers;				 		/* An array of the numbers of the antennas present in the
                                                    radio telescope (note that these are not indices, they do not need to start at zero or
                                                    be continuous). This is a one-dimensional array of size Nants_telescope. Note there
                                                    must be one entry for every unique antenna in ant 1 array and ant 2 array, but there
                                                    may be additional entries. */
    char **antenna_names;				 		/* An array of the names of antennas present in the radio telescope. This is a
                                                    one-dimensional array of size Nants_telescope. Note there must be one entry for every
                                                    unique antenna in ant 1 array and ant 2 array, but there may be additional entries. */
    int Nbls;									/* The number of baselines present in the data. For full cross-correlation data (including
                                                    auto-correlations), this should be Nants_data×(Nants_data+1)/2. */
    int Nblts;									/* The number of baseline-times (i.e., the number of spectra) present in the data. Note
                                                    that this value need not be equal to Nbls × Ntimes. */
    int Nspws;									/* The number of spectral windows present in the data. */
    int Nfreqs;									/* The total number of frequency channels in the data across all spectral windows. */
    int Npols;									/* The number of polarization products in the data. */
    int Ntimes;									/* The number of time samples present in the data. */
    double *uvw_array;					 		/* An array of the uvw-coordinates corresponding to each observation in the data. This is a
                                                    two-dimensional array of size [Nbls, 3] (reiteratively written Ntimes, so that stored it
                                                    has size [Nblts, 3]). Units are in meters. */
    double *time_array;					 		/* An array of the Julian Date corresponding to the temporal midpoint of the corresponding
                                                    baseline's integration. This is a one-dimensional array of size [Nbls] (reiteratively
                                                    written Ntimes, so that stored it has size [Nblts]). */
    double *integration_time;		 			/* An array of the duration in seconds of an integration. This is a one-dimensional array
                                                    of size [Nbls] (reiteratively written Ntimes, so that stored it has size [Nblts]). */
    double *freq_array;					 		/* An array of all the frequencies (for all spectral windows) stored in the file in Hertz.
                                                    This is a one-dimensional array of size [Nfreqs]. */
    double *channel_width;			 			/* The width of frequency channels in the file in Hertz. This is a one-dimensional array of
                                                    size [Nfreqs]. */
    int *spw_array;							 	/* An array of the spectral windows in the file. This is a one-dimensional array of size
                                                    [Nspws]. */
    H5_bool_t flex_spw;				 	 		/* Whether the data are saved using flexible spectral windows. If more than one spectral
                                                    window is present in the data, this must be True. */
    int *polarization_array;					/* An array of the polarizations contained in the file. This is a one-dimensional array of
                                                    size [Npols]. Note that the polarizations should be stored as an integer, and use the
                                                    convention defined in AIPS Memo 117. */
    double *antenna_positions;	 				/* An array of the antenna coordinates relative to the reference position of the radio
                                                    telescope array, which is implicitly defined by the latitude, longitude, and altitude
                                                    (LLA) parameters. More explicitly, these are the ECEF coordinates of individual antennas
                                                    minus the ECEF coordinates of the reference telescope position, such that the telescope
                                                    position plus the values stored in antenna positions equals the position of individual
                                                    elements in ECEF. The conversion between LLA and ECEF is given by WGS84. This is a
                                                    two-dimensional array of size [Nants_telescope, 3]. */

    UVH5_phase_center_t* phase_center_catalog;	/* A list of phase-center details (UVH5_phase_center_t). The integer-indices
                                                    are the phase center catalog IDs, which are used to identify which
                                                    baseline-times are phased to which phase center via the
                                                    phase center id array.*/
    int* phase_center_id_array;					/* The cat id from the phase_center_catalog that each baseline-time is phased to.
                                                    A one dimensional array of length [Nbls] (reiteratively written Ntimes, so that
                                                    stored it has size [Nblts]). */
    double* phase_center_app_ra;  	 			/* Apparent right ascension of the phase center in the topocentric frame of the observatory,
                                                    in radians. This is a one-dimensional array of size [Nbls] (reiteratively written Ntimes,
                                                    so that stored it has size [Nblts]).
                                                    In the event that there are multiple phase centers, the phase center id array
                                                    can be used to identify which phase center is used for this calculation. For unprojected
                                                    phase types, this is just the apparent LST (LAST). */
    double* phase_center_app_dec;		 		/* Apparent declination of the phase center in the topocentric frame of the observatory, in
                                                    radians. This is a one-dimensional array of size [Nbls] (reiteratively written Ntimes, so
                                                    that stored it has size [Nblts]).
                                                    In the event that there are multiple phase centers, the phase center id array can be
                                                    used to identify which phase center is used for this calculation. For unprojected phase
                                                    types, this is just the telescope latitude. */
    double* phase_center_frame_pa;	 			/* Position angle between the hour circle (which is
                                                    a great circle that goes through the target postion and both poles) in the
                                                    apparent/topocentric frame, and the frame given in the phase center catalog under the cat
                                                    frame dataset. This is a one dimensional array of length [Nbls] (reiteratively written
                                                    Ntimes, so that stored it has size [Nblts]).
                                                    In the event that there are multiple phase centers with different frames, the
                                                    phase center id array can be used to identify which frame is used for each baseline-time
                                                    in this calculation.
                                                    This is set to zero for unprojected phase types. */

    // Optional entries follow
    int* flex_spw_id_array;						/* The mapping of individual channels along the frequency axis
                                                    to individual spectral windows, as listed in the spw array. This is a one-dimensional
                                                    array of size [Nfreqs]. Note this is required if the file uses flexible spectral windows
                                                    (see Sec. 5.2.2).*/
    double dut1;								/* DUT1 (google it), in seconds. AIPS 117 calls it UT1UTC.*/
    double earth_omega;					 		/* Earth's rotation rate in degrees per day. Note the dif`ference in
                                                    units, which is inherited from the way this quantity is handled in UVFITS datasets
                                                    (AIPS 117 calls it DEGPDY).*/
    double gst0;								/* Greenwich sidereal time at midnight on reference date, in degrees. AIPS
                                                    117 calls it GSTIA0. */
    char *rdate;								/* Date for which GST0 (or whichever time saved in that field) applies.
                                                    Note this is different from how UVFITS handles this quantity, which is saved as
                                                    a float rather than a string. The user is encouraged to ensure it is being handled
                                                    self-consistently for their desired application. */
    char* timesys;							 	/* Time system. pyuvdata currently only supports UTC. */
    char* x_orientation;				 		/* The orientation of the x-arm of a dipole antenna. It is assumed
                                                    to be the same for all antennas in the dataset. For instance, "East" or "North" may
                                                    be used. */
    double* antenna_diameters;	 				/* An array of the diameters of the antennas in meters.
                                                    This is a one-dimensional array of size [Nants_telescope]. */
    int uvplane_reference_time;	 				/* The time at which the phase center is normal to the
                                                    chosen UV plane for phasing. Used for interoperability with the FHD package. */
    double* lst_array;					 		/* An array corresponding to the local sidereal time of the center of
                                                    each observation in the data in units of radians. If it is not specified, it is calculated
                                                    from the latitude/longitude and the time array. Saving it in the file can be useful
                                                    for files with many values in the time array, which would expensive to recompute. */
    
    // Inconsequential extra-keyword entries follow
    
    char *reference_antenna_name;				 /* The name of the antenna to be used as a reference. Convenience entry read from telinfo
                                                    toml file if present under the same field name. */

    // Administrative entries follow
    int* _antenna_num_idx_map;	 				/* An array whose elements translate the corresponding index from antenna-number to index
                                                    in the other antenna_* arrays. */
    double* _antenna_enu_positions; 		 	/* The ENU-framed antenna_positions. */
    double* _antenna_uvw_positions; 		 	/* The UVW-framed antenna_positions. */
    int* _ant_pol_prod_xgpu_index; 			 	/* An array whose elements are the xGPU indices. Size of [Nbls*Npols]*/
    int* _ant_pol_prod_bl_index;	 			/* An array whose elements are the ant_1/2_array index. Size of [Nbls*Npols]*/
    int* _ant_pol_prod_pol_index; 			 	/* An array whose elements are the polarization index. Size of [Nbls*Npols]*/
    char* _ant_pol_prod_conj; 		 			/* An array whose elements indicate if the antenna-polarization-product needs
                                                    conjugation. Size of [Nbls*Npols]*/
    char* _ant_pol_prod_auto; 		 			/* An array whose elements indicate if the antenna-polarization-product is of auto-
                                                    correlation. Size of [Nbls*Npols]*/
    int* _antenna_numbers_data;  	 			/* An array of the numbers of the antennas present in the
                                                    observation (note that these are not indices, they do not need to start at zero or
                                                    be continuous). This is a one-dimensional array of size Nants_data. Note there
                                                    must be one entry for every unique antenna in ant 1 array and ant 2 array, and no
                                                    additional entries. */
    size_t _phase_center_catalog_length; 		/* The length of the phase_center_catalog.*/

} UVH5_header_t;

void UVH5Halloc(UVH5_header_t *header);

void UVH5Hmalloc_phase_center_catalog(UVH5_header_t *header, size_t catalog_length);

void UVH5Hadmin(UVH5_header_t *header);

void UVH5Hfree(UVH5_header_t *header);

typedef struct __attribute__((__packed__))
{
    int32_t r;
    int32_t i;
} UVH5_CI32_t;
hid_t UVH5TcreateCI32();

typedef struct __attribute__((__packed__))
{
    float r;
    float i;
} UVH5_CF32_t;
hid_t UVH5TcreateCF32();

typedef struct __attribute__((__packed__))
{
    double r;
    double i;
} UVH5_CF64_t;
hid_t UVH5TcreateCF64();

typedef struct
{
    hid_t file_id;
    hid_t header_id;
    hid_t keywords_id;
    UVH5_header_t header;
    H5_open_dataspace_t DS_header_Ntimes;
    H5_open_dataspace_t DS_header_Nblts;
    H5_open_dataspace_t DS_header_ant_1_array;
    H5_open_dataspace_t DS_header_ant_2_array;
    H5_open_dataspace_t DS_header_uvw_array;
    H5_open_dataspace_t DS_header_time_array;
    H5_open_dataspace_t DS_header_integration_time;
    H5_open_dataspace_t DS_header_lst_array;
    H5_open_dataspace_t DS_phase_center_id_array;
    H5_open_dataspace_t DS_phase_center_app_ra;
    H5_open_dataspace_t DS_phase_center_app_dec;
    H5_open_dataspace_t DS_phase_center_frame_pa;
    hid_t data_id;
    void *visdata;
    H5_bool_t *flags;
    double *nsamples;
    H5_open_dataspace_t DS_data_visdata;
    H5_open_dataspace_t DS_data_flags;
    H5_open_dataspace_t DS_data_nsamples;
} UVH5_file_t;

typedef struct
{
    char* antenna;
    char polarization;
} UVH5_inputpair_t;

void UVH5open(const char* filepath, UVH5_file_t *UVH5file, hid_t Tvisdata);
void UVH5open_with_fileaccess(const char* filepath, UVH5_file_t *UVH5file, hid_t Tvisdata, hid_t Pfapl);

void UVH5close(UVH5_file_t *UVH5file);

int UVH5write_dynamic(UVH5_file_t* UVH5file);

int UVH5write_keyword_double(UVH5_file_t* UVH5file, char* key, double value);
int UVH5write_keyword_bool(UVH5_file_t* UVH5file, char* key, bool value);
int UVH5write_keyword_int(UVH5_file_t* UVH5file, char* key, int value);
int UVH5write_keyword_string(UVH5_file_t* UVH5file, char* key, char* value);

int UVH5find_antenna_index_by_name(UVH5_header_t* header, char* name);

void UVH5parse_input_map(
    UVH5_header_t* header,
    UVH5_inputpair_t* inputs
);

void UVH5permute_uvws(UVH5_header_t* header);

int UVH5polarisation_string_key(char* pol_string, int npols);

void UVH5visdata_from_xgpu_double_output(
    UVH5_CF64_t* xgpuOutput, // [freq, xgpu_antpol_prod]
    UVH5_CF64_t* visdata, // [bl, freq, antpol_prod]
    size_t xgpuElements,
    UVH5_header_t* header
);

void UVH5visdata_from_xgpu_float_output(
    UVH5_CF32_t* xgpuOutput, // [freq, xgpu_antpol_prod]
    UVH5_CF32_t* visdata, // [bl, freq, antpol_prod]
    size_t xgpuElements,
    UVH5_header_t* header
);

void UVH5visdata_from_xgpu_int_output(
    UVH5_CI32_t* xgpuOutput, // [freq, xgpu_antpol_prod]
    UVH5_CI32_t* visdata, // [bl, freq, antpol_prod]
    size_t xgpuElements,
    UVH5_header_t* header
);

#ifdef __cplusplus
}
#endif

#endif