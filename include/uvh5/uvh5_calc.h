#ifndef _UVH5_CALC_H
#define _UVH5_CALC_H

#include <stddef.h>
#include <math.h>
#include "_geodesy.h"

#define DAYSEC 86400.0

enum position_frames {
	FRAME_ENU,
	FRAME_XYZ,
	FRAME_ECEF
};

static inline double deg2rad(double deg) {return (deg/180)*M_PI;};

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

// positions is of length pos_count*3
void frames_translate(double* positions, int pos_count, double translation[3]);

/*
 * double ecef[3] is x, y, z
 */
void populate_ecef_from_lla(
	double* ecef,
	const double longitude,
	const double latitude,
	const double altitude,
	const geodesy_t* geo
);

// positions is of length pos_count*3
void position_to_xyz_frame_from_ecef(
	double* positions,
	int pos_count,
	double longitude_deg,
	double latitude_deg,
	double altitude
);

#endif