#ifndef UVH5_BOOL_T_H
#define UVH5_BOOL_T_H

#include "hdf5/serial/hdf5.h"

typedef enum {
	UVH5_FALSE,
	UVH5_TRUE
} UVH5_bool_t;

#define CPTR(VAR,CONST) ((VAR)=(CONST),&(VAR))
hid_t UVH5TcreateBool();

#endif