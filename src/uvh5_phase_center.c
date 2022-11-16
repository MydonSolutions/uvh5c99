#include "uvh5/uvh5_phase_center.h"

void UVH5phase_center_malloc_ephemeris_points(UVH5_phase_center_t* phase_center, size_t points) {
  phase_center->_ephem_points = points;
  phase_center->ephem_lon = malloc(points * sizeof(float));
  phase_center->ephem_lat = malloc(points * sizeof(float));
  phase_center->times = malloc(points * sizeof(float));
  phase_center->dist = malloc(points * sizeof(float));
  phase_center->vrad = malloc(points * sizeof(float));
}