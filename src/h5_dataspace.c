#include "uvh5/h5_dataspace.h"

void H5DSopen(
	hid_t dest_id, hid_t Tmem_id, hid_t Tsto_id, H5_open_dataspace_t* dataspace
) {
	UVH5print_verbose(__FUNCTION__, "%s", dataspace->name);

	dataspace->P_id = 0;
	dataspace->Tmem_id = 0;
	dataspace->Tsto_id = 0;
	dataspace->D_id = 0;
	dataspace->S_id = 0;
	dataspace->C_id = 0;

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
		dataspace->C_id = H5Screate_simple(dataspace->rank, dataspace->dimchunks, NULL);
		// filters can only be used with chunked data
		if(dataspace->filter_flag >= H5_FILTER_FLAG_NONE){
			if(dataspace->filter_flag <= H5_FILTER_FLAG_DEFLATE_9){
				H5Pset_deflate(dataspace->P_id, dataspace->filter_flag);
			}
			if(dataspace->filter_flag == H5_FILTER_FLAG_SHUFFLE){
				H5Pset_shuffle(dataspace->P_id);
			}
		}
	}
	dataspace->D_id = H5Dcreate(dest_id, dataspace->name, dataspace->Tmem_id, dataspace->S_id, H5P_DEFAULT, dataspace->P_id, H5P_DEFAULT);
}

herr_t H5DSclose(H5_open_dataspace_t* dataspace) {
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
	if (dataspace->C_id) {
		status += H5Sclose(dataspace->C_id);
	}
	return status + H5DSfree(dataspace);
}

herr_t H5DSfree(H5_open_dataspace_t* dataspace) {
	if (dataspace->dims != NULL) {
		free(dataspace->dims);
	}
	if (dataspace->dimlims != NULL) {
		free(dataspace->dimlims);
	}
	if (dataspace->dimchunks != NULL) {
		free(dataspace->dimchunks);
	}
	return 0;
}

void* H5DSmalloc(H5_open_dataspace_t* dataspace) {
	size_t element_byte_size = H5Tget_size(dataspace->Tmem_id);
	size_t nelem = H5DSnelem_chunks(dataspace);
	if (nelem == 0) {
		nelem = H5DSnelem_lims(dataspace);
		nelem = nelem > 0 ? nelem : H5DSnelem(dataspace);
	}
	
	UVH5print_verbose(__FUNCTION__, "'%s' allocated %ld bytes.", dataspace->name, nelem*element_byte_size);
	return malloc(nelem*element_byte_size);
}

void H5DSset(
	int rank,
	const hsize_t* dimlims,
	const hsize_t* chunks,
	H5_open_dataspace_t* dataspace
) {
	UVH5print_verbose(__FUNCTION__, "%s", dataspace->name);
	
	dataspace->rank = rank;
	dataspace->dims = NULL;
	dataspace->dimlims = NULL;
	dataspace->dimchunks = NULL;
	if(rank > 0) {
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
			UVH5print_verbose(__FUNCTION__, "\tdim %ld (%llu/%llu)", i, dataspace->dims[i], dataspace->dimlims[i]);
			if(chunks != NULL) {
				// TODO assert chunk[i] < dimlims[i]
				dataspace->dimchunks[i] = chunks[i];
				UVH5print_verbose(__FUNCTION__, "\t\tchunked %llu", dataspace->dimchunks[i]);
			}
		}
	}
}

herr_t H5DSwrite(H5_open_dataspace_t* dataspace, const void* data) {
	if(dataspace->dimchunks != NULL){
 		return H5Dwrite(dataspace->D_id, dataspace->Tsto_id, dataspace->C_id, dataspace->S_id, H5P_DEFAULT, data);
	}
	else {
 		return H5Dwrite(dataspace->D_id, dataspace->Tsto_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	}
}

herr_t H5DSextend(H5_open_dataspace_t* dataspace) {
	// for selecting extension of space
	hsize_t *start = malloc(dataspace->rank * sizeof(hsize_t));
	memset(start, 0, dataspace->rank * sizeof(hsize_t));
	hsize_t *end = malloc(dataspace->rank * sizeof(hsize_t));
	memcpy(end, dataspace->dims, dataspace->rank * sizeof(hsize_t));

	for (size_t i = 0; i < dataspace->rank; i++)
	{
		// TODO assert at least one unlimited dimension
		if(dataspace->dimlims[i] == H5S_UNLIMITED) {
			dataspace->dims[i] += dataspace->dimchunks[i];
		}
	}
	
  herr_t status = H5Dset_extent(dataspace->D_id, dataspace->dims);
	H5Sclose(dataspace->S_id); // this
	dataspace->S_id = H5Dget_space(dataspace->D_id); // this

	// select extension of space
	H5Sselect_all(dataspace->S_id);
	H5Sselect_hyperslab(dataspace->S_id, H5S_SELECT_NOTB, start, NULL, end, NULL);

	free(start);
	free(end);
	return status;
}

void H5DSopenBool(
	hid_t dest_id, H5_open_dataspace_t* dataspace
) {
	H5DSopen(dest_id, UVH5TcreateBool(), UVH5TcreateBool(), dataspace);
}

void H5DSopenDouble(
	hid_t dest_id, H5_open_dataspace_t* dataspace
) {
	H5DSopen(dest_id, H5Tcopy(H5T_NATIVE_DOUBLE), H5Tcopy(H5T_NATIVE_DOUBLE), dataspace);
}

void H5DSopenFloat(
	hid_t dest_id, H5_open_dataspace_t* dataspace
) {
	H5DSopen(dest_id, H5Tcopy(H5T_NATIVE_FLOAT), H5Tcopy(H5T_NATIVE_FLOAT), dataspace);
}

void H5DSopenInt(
	hid_t dest_id, H5_open_dataspace_t* dataspace
) {
	H5DSopen(dest_id, H5Tcopy(H5T_NATIVE_INT), H5Tcopy(H5T_NATIVE_INT), dataspace);
}

size_t H5DSnelem(H5_open_dataspace_t* dataspace) {
	size_t nelem = 1;
	for (int i = 0; i < dataspace->rank; i++)
	{
		nelem *= dataspace->dims[i];
	}
	return nelem;
}

size_t H5DSnelem_chunks(H5_open_dataspace_t* dataspace) {
	if (dataspace->dimchunks == NULL) {
		return 0;
	}
	size_t nelem = 1;
	for (int i = 0; i < dataspace->rank; i++)
	{
		nelem *= dataspace->dimchunks[i];
	}
	return nelem;
}

size_t H5DSnelem_lims(H5_open_dataspace_t* dataspace) {
	size_t nelem = 1;
	for (int i = 0; i < dataspace->rank; i++)
	{
		if(dataspace->dimlims[i] == H5S_UNLIMITED && nelem > 0) {
			nelem *= -1;
		}
		else {
			nelem *= dataspace->dimlims[i];
		}
	}
	return nelem;
}
