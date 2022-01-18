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
