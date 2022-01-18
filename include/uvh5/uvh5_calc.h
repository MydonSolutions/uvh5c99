#ifndef _UVH5_CALC_H
#define _UVH5_CALC_H

#include <stddef.h>

#define DAYSEC 86400.0

float julian_date_from_unix(float unix_sec);

float julian_date_from_guppi_param(
	float tbin,
	size_t sampleperblk,
	size_t piperblk,
	size_t synctime,
	size_t pktidx
);

#endif