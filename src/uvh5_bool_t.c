#include "uvh5/uvh5_bool_t.h"

hid_t UVH5TcreateBool() {
	hid_t UVH5_bool = H5Tcreate(H5T_ENUM, sizeof(UVH5_bool_t));
	UVH5_bool_t val;
	H5Tenum_insert(UVH5_bool, "FALSE", CPTR(val, UVH5_FALSE));
	H5Tenum_insert(UVH5_bool, "TRUE" , CPTR(val, UVH5_TRUE));
	return UVH5_bool;
}