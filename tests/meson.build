test('write', executable(
  'write', ['write.c'],
	dependencies: lib_uvh5_dep,
	),
	args : files(['telinfo_ata.toml', 'obsinfo.toml']),
	is_parallel: false, timeout: 0
)
test('visdata_from_xgpu', executable(
  'visdata_from_xgpu', ['visdata_from_xgpu.c', 'matrix_h.h', 'visdata.h'],
	dependencies: lib_uvh5_dep,
	),
	args : files(['telinfo_ata.toml', 'obsinfo.toml']),
	is_parallel: false, timeout: 0
)
