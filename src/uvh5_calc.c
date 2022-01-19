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

float hypotenuse_f(float* position, int dims) {
	double sum = 0.0;
	while(--dims > 0) {
		sum += position[dims]*position[dims];
	}
	return (float) sqrt(sum);
}

double hypotenuse(double* position, int dims) {
	double sum = 0.0;
	while(--dims > 0) {
		sum += position[dims]*position[dims];
	}
	return sqrt(sum);
}

void frames_translate(double* positions, int position_count, double translation[3]) {
	for (int i = 0; i < position_count; i++)
	{
		positions[i*3+0] += translation[0];
		positions[i*3+1] += translation[1];
		positions[i*3+2] += translation[2];
	}
}

void ecef_from_lla(
	double ecef[3],
	const double longitude_rad,
	const double latitude_rad,
	const double altitude,
	const geodesy_t* geo
) {
    double sin_phi = sin(latitude_rad);
		double cos_phi = cos(latitude_rad);
    double sin_lambda = sin(longitude_rad);
		double cos_lambda = cos(longitude_rad);

    double N = geo->a / sqrt(1 - geo->e2 * sin_phi * sin_phi);  // Radius of curvature (meters)

    ecef[0] = (N + altitude) * cos_phi * cos_lambda;
    ecef[1] = (N + altitude) * cos_phi * sin_lambda;
    ecef[2] = (N * (1 - geo->e2) + altitude) * sin_phi;
}

static inline void _rotate_around_x_cached_trig(
	double vec[3],
	double sin_val,
	double cos_val
) {
	double y, z;
	y = vec[1];
	z = vec[2];
	vec[1] = cos_val*y - sin_val*z;
	vec[2] = sin_val*y + cos_val*z;
}
void rotate_around_x(
	double vec[3],
	double radians
) {
	_rotate_around_x_cached_trig(
		vec, sin(radians), cos(radians)
	);
}

static inline void _rotate_around_y_cached_trig(
	double vec[3],
	double sin_val,
	double cos_val
) {
	double x, z;
	x = vec[0];
	z = vec[2];
	vec[0] = cos_val*x + sin_val*z;
	vec[2] = -sin_val*x + cos_val*z;
}
void rotate_around_y(
	double vec[3],
	double radians
) {
	_rotate_around_y_cached_trig(
		vec, sin(radians), cos(radians)
	);
}

static inline void _rotate_around_z_cached_trig(
	double vec[3],
	double sin_val,
	double cos_val
) {
	double x, y;
	x = vec[0];
	y = vec[1];
	vec[0] = cos_val*x - sin_val*y;
	vec[1] = sin_val*x + cos_val*y;
}
void rotate_around_z(
	double vec[3],
	double radians
) {
	_rotate_around_z_cached_trig(
		vec, sin(radians), cos(radians)
	);
}

void rotate_enu_by_hda(
	double* enu,
	double hour_angle_rad,
	double declination_rad,
	double latitude_rad
) {
	//  rx(-dec_rad) * ry(-ha_rad) * rx(lat_rad)
	rotate_around_x(enu, -declination_rad); // anti-clockwise
	rotate_around_y(enu, -hour_angle_rad); // anti-clockwise
	rotate_around_x(enu, latitude_rad); // clockwise
}

void _rotate_enu_by_hda_cached_trig(
	double* enu,
	double sin_hour_angle,
	double cos_hour_angle,
	double sin_declination,
	double cos_declination,
	double sin_latitude,
	double cos_latitude
) {
	//  rx(-dec_rad) * ry(-ha_rad) * rx(lat_rad)
	_rotate_around_x_cached_trig(enu, -sin_declination, cos_declination); // anti-clockwise
	_rotate_around_y_cached_trig(enu, -sin_hour_angle, cos_hour_angle); // anti-clockwise
	_rotate_around_x_cached_trig(enu, sin_latitude, cos_latitude); // clockwise
}

void position_to_xyz_frame_from_ecef(
	double* positions,
	int position_count,
	double longitude_rad,
	double latitude_rad,
	double altitude
) {
	double ecef[3];
	geodesy_t wgs84 = {0};
	geodesy_from_af(&wgs84, WGS84_A_METERS, WGS84_F);

	ecef_from_lla(
		ecef,
		longitude_rad,
		latitude_rad,
		altitude,
		&wgs84
	);
	ecef[0] *= -1.0;
	ecef[1] *= -1.0;
	ecef[2] *= -1.0;
	frames_translate(positions, position_count, ecef);
}

void position_to_ecef_frame_from_xyz(
	double* positions,
	int position_count,
	double longitude_rad,
	double latitude_rad,
	double altitude
) {
	double ecef[3];
	geodesy_t wgs84 = {0};
	geodesy_from_af(&wgs84, WGS84_A_METERS, WGS84_F);

	ecef_from_lla(
		ecef,
		longitude_rad,
		latitude_rad,
		altitude,
		&wgs84
	);
	frames_translate(positions, position_count, ecef);
}

void position_to_xyz_frame_from_enu(
	double* positions,
	int position_count,
	double longitude_rad,
	double latitude_rad,
	double altitude // Not used
) {
	double sin_longitude = sin(longitude_rad);
	double cos_longitude = cos(longitude_rad);
	double sin_latitude = sin(latitude_rad);
	double cos_latitude = cos(latitude_rad);
	double tmp;
	for (int i = 0; i < position_count; i++)
	{
		// RotX(longitude) anti-clockwise
		_rotate_around_x_cached_trig(
			positions + i*3,
			-sin_latitude,
			cos_latitude
		);
		// RotY(longitude) clockwise
		_rotate_around_y_cached_trig(
			positions + i*3,
			sin_longitude,
			cos_longitude
		);
		// Permute (YZX) to (XYZ)
		tmp = positions[i*3 + 2]; // save X
		positions[i*3 + 2] = positions[i*3 + 1]; // move Z
		positions[i*3 + 1] = positions[i*3 + 0]; // move Y
		positions[i*3 + 0] = tmp; // move X
	}
}

void position_to_enu_frame_from_xyz(
	double* positions,
	int position_count,
	double longitude_rad,
	double latitude_rad,
	double altitude // Not used
) {
	double sin_longitude = sin(longitude_rad);
	double cos_longitude = cos(longitude_rad);
	double sin_latitude = sin(latitude_rad);
	double cos_latitude = cos(latitude_rad);
	double tmp;
	for (int i = 0; i < position_count; i++)
	{
		// RotZ(longitude) anti-clockwise
		_rotate_around_z_cached_trig(
			positions + i*3,
			-sin_longitude,
			cos_longitude
		);
		// RotY(longitude) clockwise
		_rotate_around_y_cached_trig(
			positions + i*3,
			sin_latitude,
			cos_latitude
		);
		// Permute (UEN) to (ENU)
		tmp = positions[i*3 + 0];
		positions[i*3 + 0] = positions[i*3 + 1];
		positions[i*3 + 1] = positions[i*3 + 2];
		positions[i*3 + 2] = tmp;
	}
}

void position_to_enu_frame_from_ecef(
	double* positions,
	int position_count,
	double longitude_rad,
	double latitude_rad,
	double altitude
) {
  position_to_xyz_frame_from_ecef(
		positions,
		position_count,
		longitude_rad,
		latitude_rad,
		altitude
	);
  position_to_enu_frame_from_xyz(
		positions,
		position_count,
		longitude_rad,
		latitude_rad,
		altitude
	);
}

void position_to_ecef_frame_from_enu(
	double* positions,
	int position_count,
	double longitude_rad,
	double latitude_rad,
	double altitude
) {
  position_to_xyz_frame_from_enu(
		positions,
		position_count,
		longitude_rad,
		latitude_rad,
		altitude
	);
  position_to_ecef_frame_from_xyz(
		positions,
		position_count,
		longitude_rad,
		latitude_rad,
		altitude
	);
}

void uvws_from_enu_radec_timemjd_lla(
	double* enu2uvws,
	int position_count,
	double ra_rad, double dec_rad,
	double time_mjd,
	double time_dut1,
	double longitude_rad,
	double latitude_rad,
	double altitude
) {
	double aob, zob, hob, dob, rob, eo;

	eraAtco13(
		ra_rad, dec_rad,
		0, 0, 0, 0,
		time_mjd, 0,
		time_dut1,
		longitude_rad, latitude_rad, altitude,
		0, 0,
		0, 0, 0, 0,
		&aob, &zob, &hob, &dob, &rob, &eo
	);
	
	double sin_hour_angle = sin(hob);
	double cos_hour_angle = cos(hob);
	double sin_declination = sin(dob);
	double cos_declination = cos(dob);
	double sin_latitude = sin(latitude_rad);
	double cos_latitude = cos(latitude_rad);

	while(--position_count >= 0) {
		_rotate_enu_by_hda_cached_trig(
			enu2uvws + position_count*3,
			sin_hour_angle,
			cos_hour_angle,
			sin_declination,
			cos_declination,
			sin_latitude,
			cos_latitude
		);
	}
}
