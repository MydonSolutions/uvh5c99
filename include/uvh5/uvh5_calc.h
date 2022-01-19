#ifndef _UVH5_CALC_H
#define _UVH5_CALC_H

#include <stddef.h>
#include <math.h>
#include "_geodesy.h"

#define DAYSEC 86400.0

static inline double deg2rad(double deg) {return (deg/180)*M_PI;};

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

void frame_translate(double* positions, int position_count, double translation[3]);

void ecef_from_lla(
	double ecef[3],
	const double longitude_rad,
	const double latitude_rad,
	const double altitude,
	const geodesy_t* geo
);

void rotate_around_x(double vec[3], double radians);
void rotate_around_y(double vec[3], double radians);
void rotate_around_z(double vec[3], double radians);

void position_to_xyz_frame_from_ecef(
	double* positions,
	int position_count,
	double longitude_rad,
	double latitude_rad,
	double altitude
);

void position_to_ecef_frame_from_xyz(
	double* positions,
	int position_count,
	double longitude_rad,
	double latitude_rad,
	double altitude
);

void position_to_xyz_frame_from_enu(
	double* positions,
	int position_count,
	double longitude_rad,
	double latitude_rad,
	double altitude // Not used
);

void position_to_enu_frame_from_xyz(
	double* positions,
	int position_count,
	double longitude_rad,
	double latitude_rad,
	double altitude // Not used
);

void position_to_enu_frame_from_ecef(
	double* positions,
	int position_count,
	double longitude_rad,
	double latitude_rad,
	double altitude
);

void position_to_ecef_frame_from_enu(
	double* positions,
	int position_count,
	double longitude_rad,
	double latitude_rad,
	double altitude
);

void position_to_uvw_frame_from_enu(
	double* positions,
	int position_count,
	double hour_angle_rad,
	double declination_rad,
	double latitude_rad
);

void position_to_uvw_frame_from_xyz(
	double* positions,
	int position_count,
	double hour_angle_rad,
	double declination_rad,
	double longitude_rad
);

#endif