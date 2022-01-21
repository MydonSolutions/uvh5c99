using HDF5

function scrape_matrix_h_and_visdata(
	uvh5_filepath::AbstractString
)
	h5open(uvh5_filepath, "r") do h5fio
		uvh5debug = h5fio["/Debug"]
		open("matrix_h.h", "w") do matfio
			matrix_h = read(uvh5debug, "matrix_h")[:,:,1]
			write(matfio, "#ifndef _MATRIX_H_H\n")
			write(matfio, "#define _MATRIX_H_H\n\n")
			write(matfio, "#define MATRIX_H_RANK ", string(length(size(matrix_h))+1), "\n")
			write(matfio, "const int MATRIX_H_DIMS[MATRIX_H_RANK] = {", join(reverse(size(matrix_h)), ","), ",2};\n\n")
			
			write(matfio, "const int MATRIX_H[] = {\n\t")
			map(i->write(matfio, replace(replace(string(i), r"\s\+?([-]?)\s*"=>s", \1"), "im"=>",\n\t")), matrix_h)
			write(matfio, "};\n\n#endif // _MATRIX_H_H\n")
		end
		uvh5data = h5fio["/Data"]
		open("visdata.h", "w") do vdfio
			visdata = read(uvh5data, "visdata")
			visdata = visdata[:,:,1:div(size(visdata, 3), 2)]
			write(vdfio, "#ifndef _VISDATA_H\n")
			write(vdfio, "#define _VISDATA_H\n\n")
			write(vdfio, "#define VISDATA_RANK ", string(length(size(visdata))+1), "\n")
			write(vdfio, "const int VISDATA_DIMS[VISDATA_RANK] = {", join(reverse(size(visdata)), ","), ",2};\n\n")
			
			write(vdfio, "const int VISDATA[] = {\n\t")
			map(i->write(vdfio, replace(replace(string(i), r"\s\+?([-]?)\s*"=>s", \1"), "im"=>",\n\t")), visdata)
			write(vdfio, "};\n\n#endif // _VISDATA_H\n")
		end
	end
end

scrape_matrix_h_and_visdata("/home/sonata/dev/hdf5_tests/guppi_59596_67211_481811_3c380_0001.debug.uvh5")