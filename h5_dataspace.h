#ifndef H5_DATASPACE_H
#define H5_DATASPACE_H

#include "hdf5/serial/hdf5.h"
#include <string.h>
#include <stdlib.h>

#include "uvh5_bool_t.h"

// #define H5_DATASPACE_DEBUG

enum H5_FILTER_FLAG {
	H5_FILTER_FLAG_NONE = -1,
	H5_FILTER_FLAG_DEFLATE_0 = 0,
	H5_FILTER_FLAG_DEFLATE_1 = 1,
	H5_FILTER_FLAG_DEFLATE_2 = 2,
	H5_FILTER_FLAG_DEFLATE_3 = 3,
	H5_FILTER_FLAG_DEFLATE_4 = 4,
	H5_FILTER_FLAG_DEFLATE_5 = 5,
	H5_FILTER_FLAG_DEFLATE_6 = 6,
	H5_FILTER_FLAG_DEFLATE_7 = 7,
	H5_FILTER_FLAG_DEFLATE_8 = 8,
	H5_FILTER_FLAG_DEFLATE_9 = 9,
	H5_FILTER_FLAG_SHUFFLE = 10
};

typedef struct {
	int rank;
	hsize_t* dims;
	hsize_t* dimlims;
	hsize_t* dimchunks;
	char* name;
	hid_t S_id;
	hid_t D_id;
	hid_t P_id;
	hid_t Tmem_id;
	hid_t Tsto_id;
	char filter_flag; 
} H5_open_dataspace_t;

static inline void H5DSopen(
	hid_t dest_id, hid_t Tmem_id, hid_t Tsto_id, H5_open_dataspace_t* dataspace
) {
	#ifdef H5_DATASPACE_DEBUG
		fprintf(stderr, "H5DSopen: %s\n", dataspace->name);
	#endif

	if(Tmem_id > 0)
		dataspace->Tmem_id = Tmem_id;
	if(Tsto_id > 0)
		dataspace->Tsto_id = Tsto_id;

	if(!dataspace->P_id) {
		dataspace->P_id = H5Pcreate(H5P_DATASET_CREATE);
	}

	if(dataspace->rank > 0) {
		dataspace->S_id = H5Screate_simple(dataspace->rank, dataspace->dims, dataspace->dimlims);
	}
	else {
		dataspace->S_id = H5Screate(H5S_SCALAR);
		H5Pset_layout(dataspace->P_id, H5D_COMPACT);
	}

	if(dataspace->rank > 0 && dataspace->dimchunks != NULL) {
		H5Pset_chunk(dataspace->P_id, dataspace->rank, dataspace->dimchunks);
		// filters can only be used with chunked data
		// if(dataspace->filter_flag >= H5_FILTER_FLAG_NONE){
		// 	if(dataspace->filter_flag <= H5_FILTER_FLAG_DEFLATE_9){
		// 		H5Pset_deflate(dataspace->P_id, dataspace->filter_flag);
		// 	}
		// 	if(dataspace->filter_flag == H5_FILTER_FLAG_SHUFFLE){
		// 		H5Pset_shuffle(dataspace->P_id);
		// 	}
		// }
	}
	dataspace->D_id = H5Dcreate(dest_id, dataspace->name, dataspace->Tmem_id, dataspace->S_id, H5P_DEFAULT, dataspace->P_id, H5P_DEFAULT);
}

static inline herr_t H5DSclose(H5_open_dataspace_t* dataspace) {
	herr_t status = 0;
	if (dataspace->P_id) {
		status += H5Pclose(dataspace->P_id);
	}
	if (dataspace->Tmem_id) {
		status += H5Tclose(dataspace->Tmem_id);
	}
	if (dataspace->Tsto_id) {
		status += H5Tclose(dataspace->Tsto_id);
	}
	if (dataspace->D_id) {
		status += H5Dclose(dataspace->D_id);
	}
	if (dataspace->S_id) {
		status += H5Sclose(dataspace->S_id);
	}
	return status;
}

static inline herr_t H5DSfree(H5_open_dataspace_t* dataspace) {
	if (dataspace->dims) {
		free(dataspace->dims);
	}
	if (dataspace->dimlims) {
		free(dataspace->dimlims);
	}
	if (dataspace->dimchunks) {
		free(dataspace->dimchunks);
	}
	return 0;
}

static inline void H5DSset(int rank, const hsize_t* dimlims, const hsize_t* chunks, H5_open_dataspace_t* dataspace) {
	#ifdef H5_DATASPACE_DEBUG
		fprintf(stderr, "H5DSset: %s\n", dataspace->name);
	#endif
	
	dataspace->rank = rank;
	dataspace->dims = malloc(rank * sizeof(hsize_t));
	dataspace->dimlims = malloc(rank * sizeof(hsize_t));
	if(chunks != NULL) {
		dataspace->dimchunks = malloc(rank * sizeof(hsize_t));
	}
	for (size_t i = 0; i < rank; i++)
	{
		dataspace->dimlims[i] = dimlims[i];
		dataspace->dims[i] = dimlims[i];
		if(dimlims[i] == H5S_UNLIMITED) {
			dataspace->dims[i] = 0;
		}
		#ifdef H5_DATASPACE_DEBUG
			fprintf(stderr, "\tdim %ld (%llu/%llu)\n", i, dataspace->dims[i], dataspace->dimlims[i]);
		#endif
		if(chunks != NULL) {
			// TODO assert chunk[i] < dimlims[i]
			dataspace->dimchunks[i] = chunks[i];
			#ifdef H5_DATASPACE_DEBUG
				fprintf(stderr, "\t\t chunked %llu\n", dataspace->dimchunks[i]);
			#endif
		}
	}
}

static inline herr_t H5DSwrite(H5_open_dataspace_t* dataspace, const void* data) {
	if(dataspace->dimchunks != NULL){
 		return H5Dwrite(dataspace->D_id, dataspace->Tsto_id, H5S_ALL, dataspace->S_id,  H5P_DEFAULT, data);
	}
	else {
 		return H5Dwrite(dataspace->D_id, dataspace->Tsto_id, H5S_ALL, H5S_ALL,  H5P_DEFAULT, data);
	}
}

static inline herr_t H5DSextend(H5_open_dataspace_t* dataspace) {
	// for selecting extension of space
	hsize_t *start, *end;
	start = malloc(dataspace->rank * sizeof(hsize_t));
	end = malloc(dataspace->rank * sizeof(hsize_t));

	for (size_t i = 0; i < dataspace->rank; i++)
	{
		start[i] = 0;
		end[i] = dataspace->dims[i];
		// TODO assert at least one unlimited dimension
		if(dataspace->dimlims[i] == H5S_UNLIMITED) {
			dataspace->dims[i] += dataspace->dimchunks[i];
		}
	}
	
  herr_t status = H5Dset_extent(dataspace->D_id, dataspace->dims);
	// TODO is this necessary?
	H5Sclose(dataspace->S_id);
	dataspace->S_id = H5Dget_space(dataspace->D_id);

	// select extension of space (not old space)
	H5Sselect_hyperslab (dataspace->S_id, H5S_SELECT_NOTB, start, NULL, end, NULL);
	free(start);
	free(end);
	return status;
}

static inline void H5DSopenBool(
	hid_t dest_id, H5_open_dataspace_t* dataspace
) {
	H5DSopen(dest_id, UVH5TcreateBool(), UVH5TcreateBool(), dataspace);
}

static inline void H5DSopenDouble(
	hid_t dest_id, H5_open_dataspace_t* dataspace
) {
	H5DSopen(dest_id, H5Tcopy(H5T_NATIVE_DOUBLE), H5Tcopy(H5T_NATIVE_DOUBLE), dataspace);
}

static inline void H5DSopenFloat(
	hid_t dest_id, H5_open_dataspace_t* dataspace
) {
	H5DSopen(dest_id, H5Tcopy(H5T_NATIVE_FLOAT), H5Tcopy(H5T_NATIVE_FLOAT), dataspace);
}

static inline void H5DSopenInt(
	hid_t dest_id, H5_open_dataspace_t* dataspace
) {
	H5DSopen(dest_id, H5Tcopy(H5T_NATIVE_INT), H5Tcopy(H5T_NATIVE_INT), dataspace);
}

#endif