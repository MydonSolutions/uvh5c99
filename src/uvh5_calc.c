#include "uvh5/uvh5_calc.h"

float julian_date_from_unix(float unix_sec) {
	return (unix_sec / DAYSEC) / 2440587.5;
}

float julian_date_from_guppi_param(
	float tbin,
	size_t sampleperblk,
	size_t piperblk,
	size_t synctime,
	size_t pktidx
) {
	 float tperblk = sampleperblk * tbin;
	 float tperpktidx = tperblk / piperblk;
	 return julian_date_from_unix(synctime + tperpktidx*pktidx);
}

float hypotenuse(float* position, int dims) {
	double sum = 0.0;
	while(--dims > 0) {
		sum += position[dims]*position[dims];
	}
	return (float) sqrt(sum);
}

void frames_translate(double* positions, int pos_count, double translation[3]) {
	translation[0], translation[1], translation[2], pos_count
	);
	for (int i = 0; i < pos_count; i++)
	{
		positions[i*3+0] += translation[0];
		positions[i*3+1] += translation[1];
		positions[i*3+2] += translation[2];
	}
}

/*
 * double ecef[3] is x, y, z
 */
void populate_ecef_from_lla(
	double* ecef,
	const double longitude,
	const double latitude,
	const double altitude,
	const geodesy_t* geo
) {
    double sin_phi = sin(deg2rad(latitude));
		double cos_phi = cos(deg2rad(latitude));
    double sin_lambda = sin(deg2rad(longitude));
		double cos_lambda = cos(deg2rad(longitude));

    double N = geo->a / sqrt(1 - geo->e2 * sin_phi * sin_phi);  // Radius of curvature (meters)

    ecef[0] = (N + altitude) * cos_phi * cos_lambda;
    ecef[1] = (N + altitude) * cos_phi * sin_lambda;
    ecef[2] = (N * (1 - geo->e2) + altitude) * sin_phi;
}

void position_to_xyz_frame_from_ecef(
	double* positions,
	int pos_count,
	double longitude_deg,
	double latitude_deg,
	double altitude
) {
	double ecef[3];
	geodesy_t wgs84 = {0};
	geodesy_from_af(&wgs84, WGS84_A_METERS, WGS84_F);

	populate_ecef_from_lla(
		ecef,
		longitude_deg,
		latitude_deg,
		altitude,
		&wgs84
	);
	ecef[0] *= -1.0;
	ecef[1] *= -1.0;
	ecef[2] *= -1.0;
	frames_translate(positions, pos_count, ecef);
}
