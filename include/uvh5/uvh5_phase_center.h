#ifndef UVH5_PHASE_CENTER_H
#define UVH5_PHASE_CENTER_H

enum UVH5_phase_center_type {
   UVH5_PHASE_CENTER_SIDEREAL=1,
   UVH5_PHASE_CENTER_EPHEMERIS,
   UVH5_PHASE_CENTER_DRIFTSCAN,
   UVH5_PHASE_CENTER_UNPROJECTED,
};

typedef struct
{
char* name;         /* The phase center catalog name. This does not have to be
                        unique, non-unique values can be used to indicate sets of phase centers that
                        make up a mosaic observation. */
enum UVH5_phase_center_type
      type;         /* One of four allowed values: (1) sidereal, (2) ephem, (3)
                        driftscan, (4) unprojected. Sidereal means a phase center that is fixed in RA
                        and Dec in a given celestial frame. Ephem means a phase center that has an
                        RA and Dec that moves with time. Driftscan means a phase center with a fixed
                        azimuth and elevation (note that this includes w-projection, even at zenith).
                        Unprojected means no phasing, including w-projection, has been applied. */
char* frame; 	 	/* The coordinate frame that the phase center coordinates are
                        defined in. It must be an astropy supported frame (e.g. fk4, fk5, icrs, gcrs, cirs,
                        galactic). */
double lon; 	 	/* The longitudinal coordinate of the phase center. This is commonly RA, but can also be
                        galactic longitude. It is azimuth for driftscan phase centers. */
double lat;         /* The latitudinal coordinate of the phase center. This is commonly Dec, but can also be
                        galactic latitude. It is elevation (altitude) for driftscan phase centers. */
// Optional
//  general
double epoch;       /* The epoch in years for the phase center coordinate. For most
                        frames this is the Julian epoch (e.g. 2000.0 for j2000) but for the FK4 frame
                        this will be treated as the Bessel-Newcomb epoch (e.g. 1950.0 for B1950). This
                        parameter is not used for frames without an epoch (e.g. ICRS) unless the there
                        is proper motion (specified in the cat pm ra and cat pm dec keys). */
char* info_source;  /* Information about provenance of the source details. Typically this is set either
                        to "file" if it originates from a file read operation, and "user" if it was added
                        because of a call to the phase() method in pyuvdata. But it can also be set to
                        contain more detailed information. */
//  sidereal
double pm_ra;  	 	/* (sidereal only) Proper motion in RA in milliarcseconds per
                        year for the source. */
double pm_dec; 	 	/* (sidereal only) Proper motion in Dec in milliarcseconds per
                        year for the source. */
//  ephemeris
double* ephem_lon;  /* The longitudinal coordinate of the phase center, a one dimensional array of length Npts (the number of ephemeris data
                        points) for ephem type phase centers. This is commonly RA, but can also be
                        galactic longitude. It is azimuth for driftscan phase centers. */
double* ephem_lat;  /* The latitudinal coordinate of the phase center, a one dimensional array of length Npts (the number of ephemeris data
                        points) for ephem type phase centers. This is commonly Dec, but can also be
                        galactic latitude. It is elevation (altitude) for driftscan phase centers. */
double* times;      /* Time in Julian Date for ephemeris points, a one dimensional
                        array of length Npts (the number of ephemeris data points). Only used for
                        ephem type phase centers. */
double* dist;  	 	/* Distance to the source in parsec (useful if parallax is important),
                        either a single value or a one dimensional array of length Npts (the number of
                        ephemeris data points) for ephem type phase centers. */
double* vrad;  	 	/* Radial velocity of the source in km/sec, either a single value or
                        a one dimensional array of length Npts (the number of ephemeris data points)
                        for ephem type phase centers. */

// Administrative
int _ephem_points;  /* (ephem only) The number of ephemeris data points (Npts). */
} UVH5_phase_center_t;

void UVH5phase_center_malloc_ephemeris_points(UVH5_phase_center_t* phase_center, int points);

#endif // UVH5_PHASE_CENTER_H