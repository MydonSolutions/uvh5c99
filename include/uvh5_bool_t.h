#ifndef UVH5_BOOL_T_H
#define UVH5_BOOL_T_H

#include "hdf5/serial/hdf5.h"

typedef enum {
	UVH5_FALSE,
	UVH5_TRUE
} UVH5_bool_t;

#define CPTR(VAR,CONST) ((VAR)=(CONST),&(VAR))
static hid_t UVH5TcreateBool() {
	hid_t UVH5_bool = H5Tcreate(H5T_ENUM, sizeof(UVH5_bool_t));
	UVH5_bool_t val;
	H5Tenum_insert(UVH5_bool, "FALSE", CPTR(val, UVH5_FALSE));
	H5Tenum_insert(UVH5_bool, "TRUE" , CPTR(val, UVH5_TRUE));
	return UVH5_bool;
}

#endif