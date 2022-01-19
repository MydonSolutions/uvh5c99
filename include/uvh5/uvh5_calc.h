#ifndef _UVH5_CALC_H
#define _UVH5_CALC_H

#include <stddef.h>
#include <math.h>
#include "_geodesy.h"
#include "erfa.h"

#define DAYSEC 86400.0

enum position_frames {
	FRAME_ENU,
	FRAME_XYZ,
	FRAME_ECEF,
	FRAME_UVW
};

float julian_date_from_unix(float unix_sec);

float julian_date_from_guppi_param(
	float tbin,
	size_t sampleperblk,
	size_t piperblk,
	size_t synctime,
	size_t pktidx
);


float hypotenuse_f(float* position, int dims);
double hypotenuse(double* position, int dims);

// positions is of length position_count*3
void frames_translate(double* positions, int position_count, double translation[3]);

/*
 * https://github.com/JuliaGeo/Geodesy.jl/blob/dc2b3bd4d73a5fb4ed6f2f9c5462763ac54e5196/src/transformations.jl#L175-L188
 */
void ecef_from_lla(
	double ecef[3],
	const double longitude_rad,
	const double latitude_rad,
	const double altitude,
	const geodesy_t* geo
);

/*
 * Clockwise (right-hand curl) rotations
 * y' = cos*vec.y + sin*vec.z
 * z' = -sin*vec.y + cos*vec.z
 */
void rotate_around_x(
	double vec[3],
	double radians
);

/*
 * Clockwise (right-hand curl) rotations
 * x' = cos*vec.x - sin*vec.z
 * z' = sin*vec.x + cos*vec.z
 */
void rotate_around_y(
	double vec[3],
	double radians
);

/*
 * Clockwise (right-hand curl) rotations
 * x' = cos*vec.x - sin*vec.y
 * y' = sin*vec.x + cos*vec.y
 */
void rotate_around_z(
	double vec[3],
	double radians
);

/*
 * Subtracts ECEF(LLA, WGS84) from positions.
 */
void position_to_xyz_frame_from_ecef(
	double* positions,
	int position_count,
	double longitude_rad,
	double latitude_rad,
	double altitude
);

/*
 * Adds ECEF(LLA, WGS84) to positions.
 */
void position_to_ecef_frame_from_xyz(
	double* positions,
	int position_count,
	double longitude_rad,
	double latitude_rad,
	double altitude
);

/*
 * https://github.com/david-macmahon/RadioInterferometry.jl/blob/3a084d47919ddb422be109c41faade9d365c2a35/src/RadioInterferometry.jl#L659-L662
 * Rotates the ENU frame anticlockwise about the East (i.e. first)
 * axis by `lat_rad`, producing a (East,Z,X') frame, then rotates that frame
 * anticlockwise about the Z (i.e. second) axis by `-lon_rad`, producing a
 * (Y,Z,X) frame which is then permuted to (X,Y,Z).
 */
void position_to_xyz_frame_from_enu(
	double* positions,
	int position_count,
	double longitude_rad,
	double latitude_rad,
	double altitude // Not used
);

/*
 * https://github.com/david-macmahon/RadioInterferometry.jl/blob/3a084d47919ddb422be109c41faade9d365c2a35/src/RadioInterferometry.jl#L487-490
 * Rotates the XYZ frame anticlockwise about the Z (i.e. third)
 * axis by `lon_rad`, producing a (X',East,Z) frame, then rotates that frame
 * anticlockwise about the E (i.e. second) axis by `-lat_rad`, producing a
 * (U,E,N) frame which is then permuted to (E,N,U).
 */
void position_to_enu_frame_from_xyz(
	double* positions,
	int position_count,
	double longitude_rad,
	double latitude_rad,
	double altitude // Not used
);

/*
 * Effects `ecef -> xyz -> enu`.
 */
void position_to_enu_frame_from_ecef(
	double* positions,
	int position_count,
	double longitude_rad,
	double latitude_rad,
	double altitude
);

/*
 * Effects `enu -> xyz -> ecef`.
 */
void position_to_ecef_frame_from_enu(
	double* positions,
	int position_count,
	double longitude_rad,
	double latitude_rad,
	double altitude
);

/*
 * https://github.com/david-macmahon/RadioInterferometry.jl/blob/3a084d47919ddb422be109c41faade9d365c2a35/src/RadioInterferometry.jl#L584-L589
 * Rotates the ENU frame anticlockwise about the East (i.e. first)
 * axis by `lat_rad`, producing a (East,Z,X') frame, then rotates that frame
 * anticlockwise about the Z (i.e. second) axis by `-ha_rad`, producing a
 * (U,Z,X") frame, then rotates anticlockwise about the U (i.e. first) axis by
 * `-dec_rad`, producing the (U,V,W) frame where U is east, V is north, and W is
 * in the direction of projection.
 */
void position_to_uvw_frame_from_enu(
	double* positions,
	int position_count,
	double hour_angle_rad,
	double declination_rad,
	double latitude_rad
);

/*
 * https://github.com/david-macmahon/RadioInterferometry.jl/blob/3a084d47919ddb422be109c41faade9d365c2a35/src/RadioInterferometry.jl#L422-427
 * Rotates the XYZ frame anticlockwise about the Z (i.e. third)
 * axis by `lon_rad-ha_rad`, producing a (X',U,Z) frame, then rotates that frame
 * anticlockwise about the U (i.e. second) axis by `-dec_rad`, producing an
 * (W,U,V) frame which is then permuted to (U,V,W) where U is east, V is north,
 * and W is in the direction of the given hour angle and declination as seen from
 * the given longitude.
 */
void position_to_uvw_frame_from_xyz(
	double* positions,
	int position_count,
	double hour_angle_rad,
	double declination_rad,
	double longitude_rad
);

#endif